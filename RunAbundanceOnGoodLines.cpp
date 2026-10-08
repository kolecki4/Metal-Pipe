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

// Max number of threads available to Metal Pipe
// (No performance benefit for MAX_THREADS > # of threads on your CPU)
const int MAX_THREADS = 10;











// abundances from Asplund+2009 and Asplund+2021
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




int abundanceRunOnFile(std::string paramFile, bool noIter){

    // If false, lineListFile must be just a plain text list of wavelengths, one per line
    // (These days this should always be false. Set to true in rare cases)
    bool moogatomformat = false;             

    // Width of the PSF in pixels
    // (Should this be moved to the input file somehow?)
    double nPixInstBroad = 3;
    
    // Limits for the minimization algorithm
    // (Not super important but needed)
    double maxChi2 = 99;
    double maxVBroad = 25;


    // Declare output vectors here
    // ([X/H], v_broad, chi^2 per spectral line)
    std::vector<double> abVector;
    std::vector<double> vbVector;
    std::vector<double> x2Vector;

    // Wavelength limits to synthesize from/to
    // (0s mean use the whole spectrum)
    double minWave = 0;
    double maxWave = 0;

    // Special [alpha/Fe] measurement variables 
    // (Because it's the average of multiple abundances)
    int nAlphaElementsFit = 0;
    double alpha=0;
    double stdevAlpha = 0;


    // Construct the model atmosphere according to the params.txt file
    atmosphere currentAtmosphereModel = atmosphere(paramFile);

    // Read in the complete observed stellar spectrum
    spectrumData wholeSpec = spectrumData(currentAtmosphereModel.workDir + currentAtmosphereModel.specFile, "obs");



    // For every element listed in params.txt, fit each spectral line in its linelists file
    for(size_t j = 0; j < currentAtmosphereModel.elementString.size(); j++){
        // Clear output vectors before moving on to the next element
        abVector = {};
        vbVector = {};
        x2Vector = {};

        // Tell the user what element is being fit next
        std::cout << "Complete Element list\n";
        for(size_t k = 0; k < currentAtmosphereModel.elementString.size(); k++){
            if(k == j){std::cout << "*";}

            std::cout << currentAtmosphereModel.elementString[k] + " ";
        }
        std::cout << "\n";


        // Read in all lines of a given species
        std::string lineListFile = getDefaultLinelistName(int(std::stod(currentAtmosphereModel.elementString[j])));
        std::vector<double> wavelengthsToTest;
        if(moogatomformat){
            wavelengthsToTest = readMOOGlistWavelengths(lineListFile, currentAtmosphereModel.elementString[j]);
        }
        else{
            wavelengthsToTest = readListOfWavelengths(lineListFile);
        }



        // Trim the line list to the limits of the input spectrum (or to minWave maxWave)
        std::vector<double> wavelengths = wholeSpec.getColumn("wavelength");
        if(minWave == 0){minWave = wavelengths[0];}
        if(maxWave == 0){maxWave = wavelengths[wavelengths.size()-1];}
        while(std::abs(wavelengthsToTest[0]) < minWave){
            wavelengthsToTest.erase(wavelengthsToTest.begin());
        }
        while(std::abs(wavelengthsToTest[wavelengthsToTest.size()-1]) > maxWave){
            wavelengthsToTest.pop_back();
        }

        // Fit lines only if there are lines to fit
        if(wavelengthsToTest.size() > 0){

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
                
                // Give each thread one line
                for(int i = 0; i < MAX_THREADS; i++){
                    
                    // If there's a line to give this thread, give it
                    if(k < threadLineLists[i].size()){
                        // Declare a line to synthesize, tell the program what observed spectrum we're comparing to
                        newLine[i] = line( std::abs(threadLineLists[i][k]),std::stoi(currentAtmosphereModel.elementString[j]), currentAtmosphereModel.useMolecules || (threadLineLists[i][k] < 0) , wholeSpec);
                        newLine[i].lineInfo.maxAllowedChi2 = maxChi2;
                        newLine[i].lineInfo.maxAllowedVBroad = maxVBroad;
                        newLine[i].lineInfo.instBroadWidthPixels = nPixInstBroad;
                        
                        // Declare abundance offsets that differ from purely scaled solar
                        newLine[i].lineInfo.abundanceOffsets = currentAtmosphereModel.customAbundances;

                        // Set the name of the output files for this line
                        newLine[i].parFileName = "MOOGout/thread" + std::to_string(i) + ".par";
                        newLine[i].stdOutFile = "MOOGout/stdOutThread"  + std::to_string(i) + ".txt";
                        newLine[i].sumOutFile = "MOOGout/sumOutThread"  + std::to_string(i) + ".txt";
                        newLine[i].smoothedOutFile = "MOOGout/smoothedOutThread"  + std::to_string(i) + ".txt";
                        newLine[i].lineMakeSuffix = "/thread" + std::to_string(i);

                        //Fit Line
                        lineFittingThreads[i] = std::thread(fitLine,std::ref(newLine[i]),std::ref(currentAtmosphereModel), std::ref(abVector),std::ref(vbVector),std::ref(x2Vector), currentAtmosphereModel.outDataDir[j]);
                    }
                    
                    // Otherwise just tell the thread to sit tight
                    else{
                        lineFittingThreads[i] = std::thread([](int a){return a;}, 0);
                    }
                }

                // Wait for all the threads to conclude
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
            
            // Write out the abundance of this element
            std::cout << "[M/H] = " << currentAtmosphereModel.MonH << "\n"; 
            std::cout << "[X/Fe] = " << medAbundance - currentAtmosphereModel.MonH << " +/- " << stdev << "\n"; 
            
            // Check if [Fe/H]_in == [Fe/H]_out. If not, exit 1         
            if(  (currentAtmosphereModel.elementString[j] == "26"  && abs(currentAtmosphereModel.MonH - medAbundance) > std::max(stdev, 0.024)) && !noIter  ){
                std::cout << "[Fe/H] not converged yet\n";
                std::cout << "Input was " << currentAtmosphereModel.MonH << "; Output was " << medAbundance << " +/- " << stdev << "\n";   
                return 1;
            }


            // Add the alpha elements to the alpha abundance variable
            if(currentAtmosphereModel.elementString[j] == "22" || currentAtmosphereModel.elementString[j] == "20" ){
                alpha += (medAbundance-currentAtmosphereModel.MonH)/2;
                stdevAlpha += pow(stdev, 2);
                nAlphaElementsFit++;

                // If the alpha elements have been fit, check if [alpha/Fe]_in == [alpha/Fe]_out. If not, exit 1  
                if(nAlphaElementsFit ==2){
                    stdevAlpha = pow(stdevAlpha,0.5);
                    if( ( abs(currentAtmosphereModel.AonM - alpha) > std::max(stdevAlpha, 0.024) ) && !noIter ){
                        std::cout << "[alpha/Fe] not converged yet\n";
                        std::cout << "Input was " << currentAtmosphereModel.AonM << "; Output was " << alpha << " +/- " << stdevAlpha << "\n";                        
                        return 1;
                    }
                }
            }

        // Update the abundance of the element in the model atmosphere
        currentAtmosphereModel.customAbundances.updateElement(std::stoi(currentAtmosphereModel.elementString[j]), medAbundance - currentAtmosphereModel.MonH - (currentAtmosphereModel.AonM)*(std::stoi(currentAtmosphereModel.elementString[j]) % 2 == 0 && std::stoi(currentAtmosphereModel.elementString[j]) > 7 && std::stoi(currentAtmosphereModel.elementString[j]) < 23));

        }
    } 
    // End of for loop. All lines have been fit successfully

    return 0;
}


int main(int argc, char* argv[]){
    std::string paramFile = argv[1];
    bool noIter = false;

    if(argc > 2){
        if (std::string(argv[2]) == "-ni"){
            noIter = true;
        }
        else{
            std::cout << "Invalid args\n";
            return 1;
        }
    }

    // Create temp output folders for helper programs
    fs::create_directory("outlines");
    fs::create_directory("outsort");
    fs::create_directory("outtemp");
    fs::create_directory("MOOGout");

    // Run abundance analysis given a parameter file
    int returnCode = abundanceRunOnFile(paramFile, noIter);

    // Remove temp folders
    fs::remove_all("outlines");
    fs::remove_all("outsort");
    fs::remove_all("outtemp");
    fs::remove_all("MOOGout");
    return returnCode;
}


