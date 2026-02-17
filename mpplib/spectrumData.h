#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>

class spectrumData{

    public:
    //Construct either a blank array or build from file
    spectrumData();
    spectrumData(std::string fileName, std::string fileType);

    // Check whether was read from obs or syn file
    std::string getGridType();

    // Get data
    std::vector<double> getColumn(std::string colName);
    std::vector<double> getRow(int rowNum);
    double avgResolution();
    int size();
    // Add data
    void setColumn(std::string colName, std::vector<double> colData);
    void appendDataRow(double wl, double of, double oe, double sf, double wt);
    void appendDataRow(std::vector<double> dataRow);

    // Write to file
    void writeToTextFile(std::string fileName);

    private:
    int length = 0;

    std::string wavelengthGrid = "neither";
    std::vector<double> wavelength;
    std::vector<double> observedFlux;
    std::vector<double> observedUncertainty;
    std::vector<double> synthedFlux;
    std::vector<double> chi2Weights;
};


spectrumData::spectrumData(){}

spectrumData::spectrumData(std::string fileName, std::string fileType){

    std::ifstream inFile;
    std::string inLine;
    double value = 0;

    inFile.open(fileName);
    if(inFile.fail()){
        throw std::runtime_error("Spectrum file " + fileName + " not found");
    }


    if(fileType == "obs"){
        wavelengthGrid = fileType;
        int i = 0;
        while(!inFile.eof()){
            if (i % 3 == 0){
                std::getline(inFile, inLine, ' ');
                //std::cout << "inLine = " << inLine << "\n";
                if (inLine.length() > 0){
                    value = std::stod(inLine);
                    //std::cout << "value  = " << value << "\n";
                    wavelength.push_back(value);}
            }

            else if (i % 3 == 1){
                std::getline(inFile, inLine, ' ');
                if (inLine.length() > 0){
                    value = std::stod(inLine);
                    observedFlux.push_back(value);
                }
            }

            else{
                std::getline(inFile, inLine, '\n');
                if (inLine.length() > 0){
                    value = std::stod(inLine);
                    observedUncertainty.push_back(value);
                    synthedFlux.push_back(0);
                    chi2Weights.push_back(1);
                    length++;
                }
            }
            i++;
        }
    }


    if(fileType == "syn"){

        std::getline(inFile, inLine, '\n');
        std::getline(inFile, inLine, '\n');

        wavelengthGrid = fileType;
        int i = 0;
        while(!inFile.eof()){
            std::getline(inFile, inLine, ' ');
            if(inLine.length() > 0){

                if (i % 2 == 0){
                    value = std::stod(inLine);
                    wavelength.push_back(value);

                }
                else{
                    value = std::stod(inLine);                    
                    synthedFlux.push_back(value);
                    observedUncertainty.push_back(0);
                    observedFlux.push_back(0);
                    chi2Weights.push_back(1);
                    length++;
                }
                i++;
            }
            
        }
        inFile.close();

        // while(synthedFlux[0] == 1.0){
        //     wavelength.erase(wavelength.begin());
        //     synthedFlux.erase(synthedFlux.begin());
        //     length--;
        // }

        // while(synthedFlux[synthedFlux.size()-1] == 1.0){
        //     wavelength.pop_back();
        //     synthedFlux.pop_back();
        //     length--;
        // }
    }
}

std::string spectrumData::getGridType(){
    return wavelengthGrid;
}

std::vector<double> spectrumData::getColumn(std::string colName){

    if (colName == "wavelength"){
        return wavelength;
    }
    else if (colName == "observedFlux"){
        return observedFlux;
    }
    else if (colName == "observedUncertainty"){
        return observedUncertainty;
    }
    else if (colName == "synthedFlux"){
        return synthedFlux;
    }
    else if (colName == "chi2Weights"){
        return chi2Weights;
    }
    else{throw std::domain_error("Class spectrumData has no column of name \"" + colName + "\"\n");}
}

std::vector<double> spectrumData::getRow(int rowNum){
    if(rowNum>length-1){throw std::domain_error("spectrumData row index out of bounds\n");}
    else{return {wavelength[rowNum], observedFlux[rowNum], observedUncertainty[rowNum], synthedFlux[rowNum], chi2Weights[rowNum]};}
}

double spectrumData::avgResolution(){
    if(length <2){return 1000;}
    double avgDifference = 0;
    for(int i = 1; i < length; i++){
        avgDifference += (wavelength[i] - wavelength[i-1]);
    }
    avgDifference /= length-1;
    return wavelength[length/2]/avgDifference;
}

int spectrumData::size(){
    return length;
}

void spectrumData::setColumn(std::string colName, std::vector<double> colData){
    if(colData.size() != size_t(length)){
        throw(std::length_error("Input column does not match size of spectrum data\n"));
    }
    else if (colName == "wavelength"){
        wavelength = colData;
    }
    else if (colName == "observedFlux"){
        observedFlux = colData;
    }
    else if (colName == "observedUncertainty"){
        observedUncertainty = colData;
    }
    else if (colName == "synthedFlux"){
        synthedFlux = colData;
    }
    else if (colName == "chi2Weights"){
        chi2Weights = colData;
    }
    else{throw std::domain_error("Class spectrumData has no column of name \"" + colName + "\"\n");}
    
}

void spectrumData::appendDataRow(double wl, double of, double oe, double sf, double wt){
    wavelength.push_back(wl);
    observedFlux.push_back(of);
    observedUncertainty.push_back(oe);
    synthedFlux.push_back(sf);
    chi2Weights.push_back(wt);
    length++;
}

void spectrumData::appendDataRow(std::vector<double> dataRow){
    appendDataRow(dataRow[0], dataRow[1], dataRow[2], dataRow[3], dataRow[4]);
}

void spectrumData::writeToTextFile(std::string fileName){
    std::ofstream outFile;
    std::string outLine;
    std::vector<double> outValues;
    outFile.open(fileName); 
    outFile << std::setprecision(10);
    outFile << "wav, obs, sigma_obs, syn, weight\n";
    for(int i = 0; i < length; i++){
        outValues = getRow(i);

        for(size_t j = 0; j < outValues.size(); j++){
            
            outFile << outValues[j];
            if(size_t(j) == outValues.size()-1){
                outFile << "\n";
            }
            else{outFile << ", ";}
        }
    }
    outFile.close();
}
