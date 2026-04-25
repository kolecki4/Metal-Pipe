#include <iostream>
#include <cmath>
#include <vector>
#include "mpplib/line.h"
#include "mpplib/spectrumData.h"
#include "mpplib/atmInterp4D.h"
#include "mpplib/minimizer.h"
#include "mpplib/atmosphere.h"
#include "mpplib/linelistDefaults.h"
#include <thread>
#include <string>
#include <sstream>

#include <filesystem>
namespace fs = std::filesystem;


const int MAX_THREADS = 10;

// abundances from Asplund+2009 and Asplund2020
const double solarMOOG[96] = {0,
       12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
        6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
        3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
        3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
        1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
        1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
       -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
        0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
        0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
       -5.00,-0.54,-5.00,-5.00,-5.00
};
const double solar2020[96] = {0,
       12.00,10.91, 0.96, 1.38, 2.70, 8.46, 7.83, 8.69, 4.40, 8.06,
        6.22, 7.55, 6.30, 7.51, 5.41, 7.12, 5.31, 6.38, 5.07, 6.30,
        3.14, 4.97, 3.90, 5.62, 5.42, 7.46, 4.94, 6.20, 4.18, 4.56,
        3.02, 3.62, 2.30, 3.34, 2.54, 3.12, 2.32, 2.83, 2.21, 2.59,
        1.47, 1.88,-5.00, 1.75, 0.78, 1.57, 0.96, 1.71, 0.80, 2.02,
        1.01, 2.18, 1.55, 2.22, 1.08, 2.27, 1.11, 1.58, 0.75, 1.42,
       -5.00, 0.95, 0.52, 1.08, 0.31, 1.10, 0.48, 0.93, 0.11, 0.85,
        0.10, 0.85,-0.15, 0.79, 0.26, 1.35, 1.32, 1.61, 0.91, 1.17,
        0.92, 1.95, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.03,
       -5.00,-0.54,-5.00,-5.00,-5.00
};




int abundanceRunOnFile(std::string paramFile, std::vector<double> &abundances, std::vector<double> &vbroads, std::vector<double> &chi2s){


    bool moogatomformat = false;             // If false, lineListFile must be just a plain text list of wavelengths, one per line

    double nPixInstBroad = 3;                // Width of the PSF in pixels
    double maxChi2 = 99;                     // Max allowed chi2
    double maxVBroad = 25;                   // Max allowed v_broad

    std::vector<double> abVector;
    std::vector<double> vbVector;
    std::vector<double> x2Vector;


    //std::cout << argv[1];


    atmosphere atmosphereInfo;
    atmosphereInfo.readFromFile(paramFile);




    double minWave = 0;                  // 0 is default
//    if(atmosphereInfo.Teff < 4500){
//        minWave = 7000;
//    }

    double maxWave = 0;


    // Read in the complete observed stellar spectrum
    spectrumData wholeSpec = spectrumData(atmosphereInfo.workDir + atmosphereInfo.specFile, "obs");

    int nAlphaElementsFit = 0;
    double alpha=0;
    double stdevAlpha = 0;


    // Begin for loop; zoom through every line list file in the array so I don't have to manually retype it every time
    for(size_t j = 0; j < atmosphereInfo.elementString.size(); j++){ 
        abVector = {};
        vbVector = {};
        x2Vector = {};

        std::cout << "Complete Element list\n";
        for(size_t k = 0; k < atmosphereInfo.elementString.size(); k++){
            if(k ==j){std::cout << "*";}

            std::cout << atmosphereInfo.elementString[k] + " ";
        }

        std::cout << "\n";


        // Read in all lines of a given species
        std::string lineListFile = getDefaultLinelistName(int(std::stod(atmosphereInfo.elementString[j])));

        std::vector<double> wavelengthsToTest;
        if(moogatomformat){wavelengthsToTest = readMOOGlistWavelengths(lineListFile, atmosphereInfo.elementString[j]);}
        else{wavelengthsToTest = readListOfWavelengths(lineListFile);}




        // Don't test lines which are outside the range of our spectrum
        std::vector<double> wavelengths = wholeSpec.getColumn("wavelength");
        if(minWave == 0){minWave = wavelengths[0];}
        if(maxWave == 0){maxWave = wavelengths[wavelengths.size()-1];}
        while(wavelengthsToTest[0] < minWave){
            wavelengthsToTest.erase(wavelengthsToTest.begin());
        }
        while(wavelengthsToTest[wavelengthsToTest.size()-1] > maxWave){
            wavelengthsToTest.pop_back();
        }
        if(wavelengthsToTest.size() > 0){



            //Set the model atmosphere MOOG will use
            interpAtmosphere(atmosphereInfo.Teff,atmosphereInfo.logg,atmosphereInfo.v_micro,atmosphereInfo.MonH, atmosphereInfo.AonM);


            // Divide the line list equally among threads
            std::vector<double> threadLineLists[MAX_THREADS];
            for(size_t i = 0; i < wavelengthsToTest.size(); i++){
                threadLineLists[i%MAX_THREADS].push_back(wavelengthsToTest[i]);
            }
            wavelengthsToTest.clear();

            // Begin setup of lines to fit
            line newLine[MAX_THREADS];
            std::thread lineFittingThreads[MAX_THREADS];

            // For every line we want to fit, fit it
            for(size_t k = 0; k < threadLineLists[0].size(); k++){
                for(int i = 0; i < MAX_THREADS; i++){

                    if(k < threadLineLists[i].size()){
                        // Declare a line to synthesize, tell the program what observed spectrum we're comparing to
                        newLine[i] = line(threadLineLists[i][k],std::stoi(atmosphereInfo.elementString[j]), atmosphereInfo.useMolecules, wholeSpec);
                        newLine[i].lineInfo.maxAllowedChi2 = maxChi2;
                        newLine[i].lineInfo.maxAllowedVBroad = maxVBroad;
                        newLine[i].lineInfo.instBroadWidthPixels = nPixInstBroad;
                        
                        // Declare abundance offsets that differ from purely scaled solar
                        newLine[i].lineInfo.abundanceOffsets = atmosphereInfo.customAbundances;
                        newLine[i].parFileName = "MOOGout/thread" + std::to_string(i) + ".par";
                        newLine[i].stdOutFile = "MOOGout/stdOutThread"  + std::to_string(i) + ".txt";
                    
                        newLine[i].sumOutFile = "MOOGout/sumOutThread"  + std::to_string(i) + ".txt";
                        newLine[i].smoothedOutFile = "MOOGout/smoothedOutThread"  + std::to_string(i) + ".txt";
                        newLine[i].lineMakeSuffix = "/thread" + std::to_string(i);
                        //Fit Line
                        lineFittingThreads[i] = std::thread(fitLine,std::ref(newLine[i]),std::ref(atmosphereInfo), std::ref(abVector),std::ref(vbVector),std::ref(x2Vector), atmosphereInfo.outDataDir[j]);
                    }

                    else{
                        lineFittingThreads[i] = std::thread([](int a){return a;}, 0);
                    }
                }

                for(int i = 0; i < MAX_THREADS; i++){
                    lineFittingThreads[i].join();
                }
            }
            
            
            
            // Calculate median and standard deviation of fitted parameters
            double medAbundance = 0;
            double stdev = 0;
            if(abVector.size()>15){
                medAbundance= median(abVector,x2Vector,3);
                stdev = stdevMedian(abVector,x2Vector,3);
            }
            
            else{
                medAbundance= median(abVector);
                stdev = stdevMedian(abVector);
            }
            std::cout << "[M/H] = " << atmosphereInfo.MonH << "\n"; 

            std::cout << "[X/Fe] = " << medAbundance - atmosphereInfo.MonH << " +/- " << stdev << "\n"; 
            if(atmosphereInfo.elementString[j] == "26"  && abs(atmosphereInfo.MonH - medAbundance) > std::max(stdev, 0.024)){
                std::cout << "[Fe/H] not converged yet\n";
                std::cout << "Input was " << atmosphereInfo.MonH << "; Output was " << medAbundance << " +/- " << stdev << "\n";   
                return 1;
            }


            if(atmosphereInfo.elementString[j] == "22" || atmosphereInfo.elementString[j] == "20" ){
                alpha += (medAbundance-atmosphereInfo.MonH)/2;
                stdevAlpha += pow(stdev, 2);

                nAlphaElementsFit++;
                //std::cout << "Alpha alement" << nAlphaElementsFit << "fit\n";
                //std::cout << "ALPHA ELEMENT ABUNDANCE: " << medAbundance-atmosphereInfo.MonH << " - " << atmosphereInfo.AonM << " > " << pow(stdevAlpha,0.5) << "\n";

                if(nAlphaElementsFit ==2){
                    stdevAlpha = pow(stdevAlpha,0.5);
                    if(abs(atmosphereInfo.AonM - alpha) > std::max(stdevAlpha, 0.024) ){
                        std::cout << "[alpha/Fe] not converged yet\n";
                        std::cout << "Input was " << atmosphereInfo.AonM << "; Output was " << alpha << " +/- " << stdevAlpha << "\n";                        
                        return 1;
                    }
                }
            }

        atmosphereInfo.customAbundances.updateElement(std::stoi(atmosphereInfo.elementString[j]), medAbundance - atmosphereInfo.MonH - (atmosphereInfo.AonM)*(std::stoi(atmosphereInfo.elementString[j]) % 2 == 0 && std::stoi(atmosphereInfo.elementString[j]) > 7 && std::stoi(atmosphereInfo.elementString[j]) < 23));
        }
    } // End of for loop

    return 0;
}


int main(int argc, char* argv[]){
    std::string paramFile = argv[1];
    std::vector<double> abunds;
    std::vector<double> vbroads;
    std::vector<double> chi2s;
    
    // Create temp output folders for helper programs
    fs::create_directory("outlines");
    fs::create_directory("outsort");
    fs::create_directory("outtemp");
    fs::create_directory("MOOGout");

    int returnCode = abundanceRunOnFile(paramFile, abunds, vbroads, chi2s);

    // Remove temp folders
    fs::remove_all("outlines");
    fs::remove_all("outsort");
    fs::remove_all("outtemp");
    fs::remove_all("MOOGout");


    return returnCode;
}
