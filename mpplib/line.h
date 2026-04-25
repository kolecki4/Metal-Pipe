#pragma once
#include <iostream>
#include <iomanip>
#include "spectrumData.h"
#include "statistics.h"
#include "externalInterfaces.h"
#include <string>
#include <fstream>
#include <vector>


struct abListEntry{
    int atomicNumber;
    double offset;
    abListEntry(int n, double o){
        atomicNumber = n;
        offset = o;
    }
};


class abundanceList{
    public:
        abListEntry getElement(int element, std::string mode);
        void addNewElement(int element);
        void addNewElement(int element, double newOffset);
        void updateElement(int element, double newOffset);
        void removeElement(int element);
        int size();
    private:
        int length = 0;
        std::vector<abListEntry> abundances;
};

abListEntry abundanceList::getElement(int element, std::string mode){
    if(length > 0){

        if(mode == "num"){
            for(int i = 0; i < length; i++){
                if (abundances[i].atomicNumber == element){
                    return abundances[i];
                }
            }
            return abListEntry(element, 0.0);
        }
        else if(mode == "ind"){
            return abundances[element];
        }
        else{return abListEntry(0,0);}
    }
    return abListEntry(0,0);
}

void abundanceList::addNewElement(int element){
    bool alreadyPresent = false;
    if(length > 0){
        for(int i = 0; i < length; i++){
            if (abundances[i].atomicNumber == element){
                alreadyPresent = true;
                break;
            }
        }
    }
    if(!alreadyPresent){
        abundances.push_back(abListEntry(element,0.0));
        length++;
    }
}

void abundanceList::addNewElement(int element, double offset){
    
    bool alreadyPresent = false;
    if(length > 0){
        for(int i = 0; i < length; i++){
            if (abundances[i].atomicNumber == element){
                alreadyPresent = true;
                break;
            }
        }
    }
    if(!alreadyPresent){
        abundances.push_back(abListEntry(element,offset));
        length++;
    }
}

void abundanceList::updateElement(int element, double offset){
        if(length > 0){
        for(int i = 0; i < length; i++){
            if (abundances[i].atomicNumber == element){
                abundances[i].offset = offset;
                break;
            }
        }
    }
}

void abundanceList::removeElement(int element){
    if(length > 0){
        for(int i = 0; i < length; i++){
            if (abundances[i].atomicNumber == element){
                abundances.erase(abundances.begin()+i);
                length--;
            }
        }
    }
}

int abundanceList::size(){
    return length;
}



















struct synthInfo{

    int species;
    double v_broad;
    abundanceList abundanceOffsets;
    
    double EP;
    double loggf;

    double LDCoeff;

    double centralWavelength;
    double widthOfSynthesis;
    double synthesisResolution;
    double instBroadWidthPixels;
    double instBroadFWHM;

    double maxAllowedChi2;
    double maxAllowedVBroad;

};

class line{

    public:
        bool considerMolecules;
        synthInfo lineInfo;

        std::string parFileName = "batch.par";
        std::string stdOutFile = "MOOGout/out1";
        std::string sumOutFile = "MOOGout/out2";
        std::string smoothedOutFile = "MOOGout/out3";
        std::string lineMakeSuffix = "";
        spectrumData obsWaveGrid;
        spectrumData synthWaveGrid;
        std::vector<std::string> listOfNearbyLines;
        std::vector<double> getFitRanges();
        line();
        line(double lineWav, int atomicNum, bool molecules,const spectrumData& cutFromSpectrum);

        void readObservedLine(std::string fileName);
        void cutOutObservedLine(double width);
        void readSynthesizedLine(std::string fileName);
        void writeParFile();
        void setLineList();
        void setEPandLoggf();
        void interpLDCoeff(std::vector<double> wave, std::vector<double> coeff);
        void MOOGitUp();


        void interpGridsToEachOther();
        void calculateFitRegions();
        void renormalizeObs();
        void crossCorrelateObs();
        void setWeightsForChi2(double vBroad);

        double obsGridChi2();
        double synGridChi2();
        double synGridGradientChi2();

        void outputBestFitAbundance(std::string fileName, double atmosphereMetallicity);

    private:
        spectrumData originalSpectrum;
        double rangeToChi2Fit = 0;
        double rangeToContinuumFit = 0;

        bool isInterpolated = false;
        bool isCrossCorrelated = false;
        bool isRenormalized = false;
        bool rangesSet = false;
        bool lineListSet = false;
        bool isWeighted = false;
};

line::line(){

}

line::line(double lineWav, int atomicNum, bool useMolecules,const spectrumData& cutFromSpectrum){
    lineInfo.centralWavelength = lineWav;
    lineInfo.species = atomicNum;
    lineInfo.instBroadFWHM = 0;
    lineInfo.v_broad = 0;
    lineInfo.widthOfSynthesis = 0;

    if(lineWav > 10000){
        lineInfo.synthesisResolution = 0.02;    
    }
    else{
        lineInfo.synthesisResolution = 0.01;    
    }

    considerMolecules = useMolecules;
    originalSpectrum = cutFromSpectrum;

}

std::vector<double> line::getFitRanges(){
    return {rangeToContinuumFit, rangeToChi2Fit};
}

void line::readObservedLine(std::string fileName){
    
    obsWaveGrid = spectrumData(fileName, "obs");
    isInterpolated = false;
    isRenormalized = false;
    isCrossCorrelated = false;
    rangesSet = false;

}

void line::interpLDCoeff(std::vector<double> wave, std::vector<double> coeff){

    lineInfo.LDCoeff = interp1DWrapper({lineInfo.centralWavelength},wave,coeff)[0];

}

void line::cutOutObservedLine(double width){

    std::vector<double> masterWavelengths = originalSpectrum.getColumn("wavelength");
    obsWaveGrid = spectrumData();
    for(int i = 0; i < originalSpectrum.size(); i++){
        if(masterWavelengths[i] >= (lineInfo.centralWavelength - (width/2)) && masterWavelengths[i] <= (lineInfo.centralWavelength + (width/2))){            
            obsWaveGrid.appendDataRow(originalSpectrum.getRow(i));
        }
    }

    if(obsWaveGrid.getColumn("wavelength").size() < 3){
        throw std::runtime_error("No data for this line");
}

    isInterpolated = false;
    isRenormalized = false;
    isCrossCorrelated = false;
}

void line::readSynthesizedLine(std::string fileName){
    synthWaveGrid = spectrumData(fileName, "syn");

    isInterpolated = false;
}

void line::writeParFile(){
    std::vector<std::string> parLines;

    std::ifstream inFile;
    std::string fLine;
    char tempString[128];
    inFile.open("mpplib/default.par");
    if(inFile.fail()){
        throw std::runtime_error("sfhsfghsdfghsdfghdfgh");
    }
    while(fLine.find("abundances") == std::string::npos){
        std::getline(inFile,fLine, '\n');


        if(fLine.find("standard_out") == 0){
            sprintf(tempString, "standard_out    '%s'", stdOutFile.c_str());
            fLine = tempString;
        }
        if(fLine.find("summary_out") == 0){
            sprintf(tempString, "summary_out     '%s'", sumOutFile.c_str());
            fLine = tempString;
        }
        if(fLine.find("smoothed_out") == 0){
            sprintf(tempString, "smoothed_out    '%s'", smoothedOutFile.c_str());
            fLine = tempString;
        }

        if(fLine.find("lines_in") == 0){
            sprintf(tempString, "lines_in        'outsort%s'", lineMakeSuffix.c_str());
            fLine = tempString;
        }

        parLines.push_back(fLine);
        if(inFile.eof()){
            throw std::runtime_error("Check that default.par exists and is formatted correctly\n");
        };
    }
    inFile.close();
    parLines.pop_back();


    sprintf(tempString, "abundances %6d    1", lineInfo.abundanceOffsets.size()  );
    parLines.push_back(tempString);

    for(int i = 0; i < lineInfo.abundanceOffsets.size(); i++){
        sprintf(tempString, "%5d%12.2f", lineInfo.abundanceOffsets.getElement(i, "ind").atomicNumber , lineInfo.abundanceOffsets.getElement(i, "ind").offset);
        parLines.push_back(tempString);
    }

    sprintf(tempString, "synlimits\n%10.2f%10.2f%10.2f%10.2f", lineInfo.centralWavelength - lineInfo.widthOfSynthesis/2, lineInfo.centralWavelength + lineInfo.widthOfSynthesis/2, lineInfo.synthesisResolution, 20.0);
    parLines.push_back(tempString);

    sprintf(tempString, "plotpars       1\n%10.2f%10.2f%10.2f%10.2f\n 0 0 0 1.0", lineInfo.centralWavelength - lineInfo.widthOfSynthesis/2 , lineInfo.centralWavelength + lineInfo.widthOfSynthesis/2 , 0.0, 1.3);
    parLines.push_back(tempString);

    sprintf(tempString, " r %.3f %.3f %.3f 0.0 0.0", lineInfo.instBroadFWHM, lineInfo.v_broad,lineInfo.LDCoeff);
    parLines.push_back(tempString);

    std::ofstream outFile;
    outFile.open(parFileName);
    for (size_t i = 0; i < parLines.size(); i++){
        if(parLines[i].length() > 0){outFile << parLines[i] << '\n';}
    }
    outFile.close();

}

void line::MOOGitUp(){
    writeParFile();
    std::string command = "printf '' | timeout 10 moog/MOOG " + parFileName + "  > MOOGout/MOOGconsoleout.txt";
    int moogStatus = std::system(command.c_str());
    if(moogStatus== 124){
        throw std::runtime_error("Took too long");
    }
}


void line::interpGridsToEachOther(){

    obsWaveGrid.setColumn("synthedFlux",interp1DWrapper(obsWaveGrid.getColumn("wavelength"),synthWaveGrid.getColumn("wavelength"), synthWaveGrid.getColumn("synthedFlux")));
    
    synthWaveGrid.setColumn("observedFlux",interp1DWrapper(synthWaveGrid.getColumn("wavelength"),obsWaveGrid.getColumn("wavelength"), obsWaveGrid.getColumn("observedFlux")));
    synthWaveGrid.setColumn("observedUncertainty",interp1DWrapper(synthWaveGrid.getColumn("wavelength"), obsWaveGrid.getColumn("wavelength"), obsWaveGrid.getColumn("observedUncertainty")));
    
    isInterpolated = true;
}

void line::calculateFitRegions(){
    double waveLow = lineInfo.centralWavelength - 0.05;
    double waveHigh = lineInfo.centralWavelength +0.05;
    
    std::string inLinemake;
    if(considerMolecules){
        inLinemake = "printf '" + std::to_string(waveLow) + "\n" + std::to_string(waveHigh) + "\ny\ny\ny\ny\ny\ny\ny\ny\ny\ny\ny\ny\n' | linemake/linemake " + lineMakeSuffix + " > linemakeconsoleout.txt";
    }
    else{
        inLinemake = "printf '" + std::to_string(waveLow) + "\n" + std::to_string(waveHigh) + "\nn\ny\n' | linemake/linemake " + lineMakeSuffix + " > linemake/linemakeconsoleout.txt";
    }
    int sysReturn = std::system(inLinemake.c_str());

    std::ifstream inFile;
    std::string fLine;
    int i = 0;
    inFile.open("linemake/linemakeconsoleout.txt");
    if(inFile.fail()){
        throw std::runtime_error("linemake/linemakeconsoleout.txt not found!");
    }
    else{

            while(!inFile.eof() && i < 10){
            std::getline(inFile,fLine, '\n');
            if(fLine.find("SORRY") != std::string::npos){
                throw std::runtime_error("Linemake won't work over this window. It can't handle the thousands digit of the wavelength range changing.\n");
            };
        }
    }
    inFile.close();


    lineInfo.widthOfSynthesis = 20.0;
    cutOutObservedLine(lineInfo.widthOfSynthesis);


    //lineInfo.instBroadFWHM = 3*lineInfo.centralWavelength/ obsWaveGrid.avgResolution();
    //lineInfo.v_broad = 3;
    
    MOOGitUp();

    readSynthesizedLine(smoothedOutFile);
    std::vector<double> point;
    double minWave=9999999;
    double maxWave=0;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        point = synthWaveGrid.getRow(i);
        if(point[3]< 0.9991){
            if(point[0]<minWave){
                minWave = point[0];
            }
            if(point[0]>maxWave){
                maxWave = point[0];
            }
        }
    } 
    rangeToContinuumFit = maxWave - minWave + 8;

    minWave = 9999999;
    maxWave = 0;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        point = synthWaveGrid.getRow(i);
        if(point[3]< 0.99){
            if(point[0]<minWave){
                minWave = point[0];
            }
            if(point[0]>maxWave){
                maxWave = point[0];
            }
        }
    } 
    rangeToChi2Fit = (maxWave - minWave + 1)*1.1;
    rangesSet = true;
}

void line::renormalizeObs(){
    if(!(rangesSet && lineListSet)){
        throw std::logic_error("Stop it you forgot something\n");
    }

    interpGridsToEachOther();
    
    std::vector<double> point;
    std::vector<double> continuumPointsObs;
    std::vector<double> continuumPointsSyn;
    std::vector<double> continuumPointsErr;

    for(int i = 0; i < obsWaveGrid.size(); i++){
        point = obsWaveGrid.getRow(i);

        if( (point[3] > 0.97 && point[3] < 1.03) && (point[1] > 0.95 && point[1] < 1.1) ) {
            continuumPointsObs.push_back(point[1]);
            continuumPointsSyn.push_back(point[3]);
            continuumPointsErr.push_back(point[2]);
        }
    }


    // If we end up with a chi^2 of inf or negative, give up
    // DoF = 3 so we need a size of at least 4
    if(continuumPointsObs.size() < 4){
        throw std::runtime_error("Length of continuum too short");
    }

    double optimalI = 1;
    double minChiSq = 1e31;
    double chiSq = 0;
    std::vector<double> temp = continuumPointsObs;

    for(double i = 0.8; i <= 1.2; i += 0.001){
        temp = continuumPointsObs;

        for(size_t j = 0; j < temp.size(); j++){
            temp[j] *= i;
        }


        chiSq = reducedChiSq(temp,continuumPointsSyn,continuumPointsErr);
        //std::cout << chiSq << std::endl;


        if(chiSq < minChiSq){
            optimalI = i;
            minChiSq = chiSq;
        }
    }

    std::vector<double> scaledObs = obsWaveGrid.getColumn("observedFlux");
    std::vector<double> scaledObsErr = obsWaveGrid.getColumn("observedUncertainty");
    for(int i = 0; i < obsWaveGrid.size(); i++){
        scaledObs[i] *= optimalI;
        scaledObsErr[i] *= optimalI;
    }
    obsWaveGrid.setColumn("observedFlux", scaledObs);
    obsWaveGrid.setColumn("observedUncertainty", scaledObsErr);
    isRenormalized = true;
}

void line::setLineList(){

    double  waveLow = lineInfo.centralWavelength - rangeToContinuumFit/2;
    double waveHigh = lineInfo.centralWavelength + rangeToContinuumFit/2;

    std::string inLinemake;
    if(considerMolecules){
        inLinemake = "printf '" + std::to_string(waveLow) + "\n" + std::to_string(waveHigh) + "\ny\ny\ny\ny\ny\ny\ny\ny\ny\ny\ny\ny\n' | linemake/linemake " + lineMakeSuffix + " > linemakeconsoleout.txt";
    }
    else{
        inLinemake = "printf '" + std::to_string(waveLow) + "\n" + std::to_string(waveHigh) + "\nn\ny\n' | linemake/linemake " + lineMakeSuffix + " > linemakeconsoleout.txt";
    }
    int sysReturn = std::system(inLinemake.c_str());

    std::ifstream inFile;
    std::string fLine;
    int i = 0;
    inFile.open("linemakeconsoleout.txt");
    while(!inFile.eof() && i < 10){
        std::getline(inFile,fLine, '\n');
        if(fLine.find("SORRY") != std::string::npos){
            throw std::runtime_error("Linemake won't work over this window. It can't handle the thousands digit of the wavelength range changing.\n");
        };
    }
    inFile.close();
    lineListSet = true;

}

void line::crossCorrelateObs(){
    if(!(isRenormalized && isInterpolated && rangesSet)){
        throw std::logic_error("Stop it you forgot something\n");
    }


    std::vector<double> weights;
    for(int i = 0; i < obsWaveGrid.size(); i++){
        weights.push_back(pow(std::min(double(i),double(obsWaveGrid.size()-i)),0.5));
    }

    double delay = 0;
    double maxCC = -999999;
    double sum;
    std::vector<double> obsGridWavelengthsShifted = obsWaveGrid.getColumn("wavelength");
    std::vector<double> point;
    for(double i = -.9; i < .9; i+=0.0005){
        sum = 0;


        for(int j = 0; j < obsWaveGrid.size(); j++){
            obsGridWavelengthsShifted[j] += i;
        }


        obsWaveGrid.setColumn("wavelength", obsGridWavelengthsShifted);
        interpGridsToEachOther();

        for(int j = 0; j < obsWaveGrid.size(); j++){
            point = obsWaveGrid.getRow(j);
            sum += pow(((1-point[3])*(1-point[1])),2)*weights[j];
            
        }
        if(sum > maxCC){
            delay = i;
            maxCC = sum;
        }

        for(int j = 0; j < obsWaveGrid.size(); j++){
            obsGridWavelengthsShifted[j] -= i;
        }
    }
    //std::cout << delay << "\n";

    for(int j = 0; j < obsWaveGrid.size(); j++){
        obsGridWavelengthsShifted[j] += delay;
    }


    obsWaveGrid.setColumn("wavelength", obsGridWavelengthsShifted);
    interpGridsToEachOther();
    isCrossCorrelated = true;
}

void line::setWeightsForChi2(double vBroad){
    if(!(isRenormalized && isInterpolated && rangesSet && isCrossCorrelated)){
        throw std::logic_error("Stop it you forgot something\n");
    }
    setLineList();
    
    line abundanceUp = *this;
    line abundanceDown = *this;
    abundanceUp.lineInfo.abundanceOffsets.updateElement(lineInfo.species,0.1);
    abundanceDown.lineInfo.abundanceOffsets.updateElement(lineInfo.species,-0.1);

    abundanceUp.lineInfo.v_broad = vBroad;
    abundanceDown.lineInfo.v_broad = vBroad;


    abundanceUp.writeParFile();
    abundanceUp.MOOGitUp();
    abundanceUp.readSynthesizedLine(smoothedOutFile);
    
    abundanceDown.writeParFile();
    abundanceDown.MOOGitUp();
    abundanceDown.readSynthesizedLine(smoothedOutFile);

    double weight = 0;
    double weightSum = 0;
    std::vector<double> weightVec;

    std::vector<double> downFlux = abundanceDown.synthWaveGrid.getColumn("synthedFlux");
    std::vector<double> upFlux = abundanceUp.synthWaveGrid.getColumn("synthedFlux");
    std::vector<double> Flux = synthWaveGrid.getColumn("observedFlux");
    std::vector<double> diffFlux;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        diffFlux.push_back(std::abs(upFlux[i] - downFlux[i]));
    }


    bool bigEnoughChange = false;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        if(i < 6 || i > synthWaveGrid.size()-6){
            weight = 0;
            weightVec.push_back(0);
        }
        else{
            weight = pow(diffFlux[i],1);
            if(weight > 0.00){bigEnoughChange = true;}
            weightVec.push_back(weight);
            weightSum += weight;
        }
    }
    if(!bigEnoughChange){
        throw std::runtime_error("Line doesn't change enough to fit well");
    }
    for(int i = 0; i < synthWaveGrid.size(); i++){
        weightVec[i] /= weightSum;
        weightVec[i] *= synthWaveGrid.size();
    }

/*
    double weightSum2 = 0;
    std::vector<double> weightVec2;

    for(int i = 0; i < synthWaveGrid.size(); i++){
        if(i < 7 || i > synthWaveGrid.size()-7){
            weight = 0;
            weightVec2.push_back(weight);
        }

        else{
            weight = pow(std::abs((Flux[i+2] - Flux[i+1]) - (Flux[i+1] - Flux[i])) ,1);
            //weight = pow(std::abs((Flux[i+1]) - (Flux[i])) ,1);
            //weight /= 2;
            weightVec2.push_back(weight);
            weightSum2 += weight;
        }
    }

    /////
    for(int i = 0; i < synthWaveGrid.size(); i++){
        weightVec2[i] /= weightSum2;
        weightVec2[i] *= synthWaveGrid.size();
    }    


    for(int i = 0; i < synthWaveGrid.size(); i++){
        weightVec[i] = (weightVec[i] * weightVec2[i]);
    }
    for(int i = 2; i < synthWaveGrid.size()-2; i++){
        weightVec2[i] += 0.32*weightVec[i-2] + 0.68*weightVec[i-1] + weightVec[i] + 0.68*weightVec[i+1] + 0.32*weightVec[i+2];
    }
    /////

    /////
    weightVec = weightVec2;
    weightSum = 0;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        weightSum +=weightVec[i];
    }

    for(int i = 2; i < synthWaveGrid.size()-2; i++){
        weightVec2[i] += 0.32*weightVec[i-2] + 0.68*weightVec[i-1] + weightVec[i] + 0.68*weightVec[i+1] + 0.32*weightVec[i+2];
    }
    /////

    weightVec = weightVec2;
    weightSum = 0;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        weightSum +=weightVec[i];
    }

    for(int i = 0; i < synthWaveGrid.size(); i++){
        weightVec[i] /= weightSum;
        weightVec[i] *= synthWaveGrid.size();
    }
*/


    synthWaveGrid.setColumn("chi2Weights", weightVec);
    obsWaveGrid.setColumn("chi2Weights", interp1DWrapper(obsWaveGrid.getColumn("wavelength"),synthWaveGrid.getColumn("wavelength"), weightVec));
    isWeighted = true;
}

double line::obsGridChi2(){
    if(!(isRenormalized && isInterpolated && rangesSet && isCrossCorrelated)){
        throw std::logic_error("Stop it you forgot something\n");
    }
    double chi2 = 0;
    std::vector<double> point;
    for(int i = 0; i < obsWaveGrid.size(); i++){
        point = obsWaveGrid.getRow(i);
        chi2 += ( ( pow(point[1] - point[3],2)/pow(point[2],2) )*point[4] ) / (obsWaveGrid.size() - 3);
    }
    return chi2;
}

double line::synGridChi2(){
    if(!(isRenormalized && isInterpolated && rangesSet && isCrossCorrelated)){
        throw std::logic_error("Stop it you forgot something\n");
    }
    double chi2 = 0;
    double err = 0;
    double synerr = 0;
    std::vector<double> point;
    for(int i = 0; i < synthWaveGrid.size(); i++){
        point = synthWaveGrid.getRow(i);
        synerr = pow(-0.05*log(point[3]), 1);
        err = pow(pow(point[2],2) + pow(synerr,2),0.5);
        chi2 += ( ( pow(point[1] - point[3],2)/pow(err,2) )*point[4] ) / (synthWaveGrid.size() - 3);
    }

    return chi2;
}

double line::synGridGradientChi2(){
    if(!(isRenormalized && isInterpolated && rangesSet && isCrossCorrelated)){
        throw std::logic_error("Stop it you forgot something\n");
    }

    std::vector<double> synFlux = synthWaveGrid.getColumn("synthedFlux");
    std::vector<double> obsFlux = synthWaveGrid.getColumn("observedFlux");
    std::vector<double> obsErr = synthWaveGrid.getColumn("observedUncertainty");
    std::vector<double> weight = synthWaveGrid.getColumn("chi2Weights");
    std::vector<double> wavelength = synthWaveGrid.getColumn("wavelength");

    std::vector<double> obsGrad;
    std::vector<double> synGrad;
    std::vector<double> gradErr;
    for(size_t i = 0; i < synFlux.size()-1; i++){
        obsGrad.push_back((obsFlux[i+1]-obsFlux[i])/(wavelength[i+1] - wavelength[i]));
        synGrad.push_back((synFlux[i+1]-synFlux[i])/(wavelength[i+1] - wavelength[i]));
        gradErr.push_back(pow(pow(obsErr[i+1],2)+pow(obsErr[i],2),0.5)/(wavelength[i+1] - wavelength[i]));
    }
    double chi2 = 0;
    for(size_t i = 0; i < synGrad.size(); i++){
        chi2 += ((pow(obsGrad[i] - synGrad[i],2)/pow(gradErr[i],2))*weight[i]) / (synGrad.size()-3);
    }


    return chi2;
}






void line::setEPandLoggf(){
    std::ifstream inFile;
    std::string inLine;
}

void line::outputBestFitAbundance(std::string fileName, double atmosphereMetallicity){
    std::ofstream outFile;
    std::string outLine;
    double wavelength = static_cast<double>(std::round(lineInfo.centralWavelength*100)/100.0);
    double solar[96] = {0,
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


    outFile.open(fileName, std::ios_base::app);
    if(!outFile){
        outFile.close();
        outFile.open(fileName);
    }

    //wavelength = ;
    outFile << std::setprecision(10);
    outFile << wavelength << " ";
    outFile << lineInfo.abundanceOffsets.getElement(lineInfo.species, "num").offset + atmosphereMetallicity + solar[lineInfo.species] << " ";
    outFile << lineInfo.v_broad << " ";
    outFile << synGridChi2() << "\n";


    outFile.close();
}
