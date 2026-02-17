#pragma once
#include <string>
#include <fstream>

std::vector<double> readMOOGlistWavelengths(std::string fileName, std::string elementStr){

    std::vector<double> wavelengths;
    double waveToAdd;
    std::ifstream inFile;
    std::string inLine;

    inFile.open(fileName);
    if(inFile.fail()){
        throw std::runtime_error("Spectrum file " + fileName + " not found");
    }


    while(!inFile.eof()){
        std::getline(inFile, inLine);
        if(inLine.find(" " + elementStr + " ") != std::string::npos){
            waveToAdd = std::stod(inLine.substr(0,12));
            
            if( (wavelengths.size() > 0) && (abs(waveToAdd - wavelengths[wavelengths.size()-1]) > 0.2)){
                wavelengths.push_back(waveToAdd);
            }

            else if(wavelengths.size() == 0){
                wavelengths.push_back(waveToAdd);
            }
        }
    }

    return wavelengths;
}

std::vector<double> readListOfWavelengths(std::string fileName){

    std::vector<double> wavelengths;
    std::ifstream inFile;
    std::string inLine;
    double wavelength;

    inFile.open(fileName);
    if(inFile.fail()){
        throw std::runtime_error("Spectrum file " + fileName + " not found");
    }


    while(!inFile.eof()){
        std::getline(inFile, inLine);
        if(inLine.length() > 0){
            wavelength = std::stod(inLine);
            //std::cout << inLine << " ";
            wavelengths.push_back(wavelength);}
    }

    return wavelengths;
}
