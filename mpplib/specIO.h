#pragma once
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "statistics.h"
/// @brief Structure to wrap up all the necessary items to send to the par file
struct parFileInfo {
    std::string parFile = "batch.par";
    int elToChange;
    double lineWav;
    double fitWindow;
    double contWindow;
    double windowSize;
    std::vector<int> elements;
    std::vector<double> offsets;
    double resolution;
    double FWHM = lineWav/resolution;
    double vsini;
    double delay;
    bool molecules;

    parFileInfo(){
        
    };

};





/// @brief Reads any spectrum (preferrably an observed stellar spectrum.) Input file must be of the form "wavelength flux" per line
/// @param fName spectrum file name 
/// @param wave empty vector to hold wavelength data
/// @param flux empty vector to hold flux data
void readObsSpec(std::string fName, std::vector<double> &wave, std::vector<double> &flux){

    wave.clear();
    flux.clear();

    std::ifstream inFile;
    std::string fLine;
    double value = 0;

    inFile.open(fName);
    if(!inFile.fail()){
        int i = 0;
        while(!inFile.eof()){
            if (i % 2 == 0){
                std::getline(inFile, fLine, ' ');
                if (fLine.length() > 0){
                    value = std::stod(fLine);
                    wave.push_back(value);
                }
            }
            else{
                std::getline(inFile, fLine, '\n');
                if (fLine.length() > 0){
                    value = std::stod(fLine);                    
                    flux.push_back(value);
                }
            }
            i++;
        }
        inFile.close();
    }
    else{
        throw std::runtime_error( "Observed spectrum file " + fName + " not found");
    }
}

/// @brief Reads any spectrum (preferrably an observed stellar spectrum.) Input file must be of the form "wavelength flux error" per line
/// @param fName spectrum file name 
/// @param wave empty vector to hold wavelength data
/// @param flux empty vector to hold flux data
/// @param error empty vector to hold uncertainty data
void readObsSpec(std::string fName, std::vector<double> &wave, std::vector<double> &flux, std::vector<double> &error){

    std::ifstream inFile;
    std::string fLine;
    double value = 0;

    wave.clear();
    flux.clear();
    error.clear();
    inFile.open(fName);
    if(!inFile.fail()){
        int i = 0;
        while(!inFile.eof()){
            if (i % 3 == 0){
                std::getline(inFile, fLine, ' ');
                if (fLine.length() > 0){
                    value = std::stod(fLine);
                    wave.push_back(value);
                }
            }
            else if (i % 3 == 1){
                std::getline(inFile, fLine, ' ');
                if (fLine.length() > 0){
                    value = std::stod(fLine);
                    flux.push_back(value);
                }
            }

            else{
                std::getline(inFile, fLine, '\n');
                if (fLine.length() > 0){
                    value = std::stod(fLine);
                    error.push_back(value);
                }
            }
            i++;
        }
        inFile.close();
    }
    else{
        throw std::runtime_error( "Observed spectrum file " + fName + " not found");
    }
}



/// @brief Reads a standard "smoothed_out" synthetic MOOG spectrum
/// @param fName spectrum file name
/// @param wave empty vector to hold wavelength data
/// @param flux empty vector to hold flux data
void readSynSpec(std::string fName, std::vector<double> &wave, std::vector<double> &flux){
    wave.clear();
    flux.clear();
    
    std::ifstream inFile;
    double value = 0;
    inFile.open(fName);

    std::string inLine;
    std::getline(inFile, inLine, '\n');
    std::getline(inFile, inLine, '\n');

    int i = 0;
    while(!inFile.eof()){
        std::getline(inFile, inLine, ' ');
        if(inLine.length() > 0){

            if (i % 2 == 0){
                try{
                    value = std::stod(inLine);
                    wave.push_back(value);
                }
                catch(...){break;}
                i++;

            }
            else{
                try{
                    value = std::stod(inLine);
                    flux.push_back(value);
                }
                catch(...){break;}
                i++;
            }
        }
        
    }
    inFile.close();

    while(flux[0] == 1.0){
        wave.erase(wave.begin());
        flux.erase(flux.begin());
    }

    while(flux[flux.size()-1] == 1.0){
        wave.pop_back();
        flux.pop_back();
    }
}

/// @brief Writes a spectrum (wave, flux only) to a file
/// @param fName 
/// @param wave 
/// @param flux 
/// @param err 
void writeSpec(std::string fName, std::vector<double> &wave, std::vector<double> &flux){
    std::ofstream outFile;
    outFile.open(fName);
    char tempString[128];

    for(int i = 0; i < wave.size(); i++){
        sprintf(tempString, "%10.4f%10.4f\n", wave[i], flux[i]);
        outFile << tempString;
    }
    outFile.close();
}

/// @brief Write a fitted line (Matched wavelemngths only!)
/// @param fName 
/// @param wave 
/// @param flux 
void writeSpecLineFit(std::string fName, std::vector<double> &wave, std::vector<double> &obsFlux, std::vector<double> &obsErr, std::vector<double> &synFlux){
    std::ofstream outFile;
    outFile.open(fName);
    char tempString[128];

    for(int i = 0; i < wave.size(); i++){
        sprintf(tempString, "%10.4f%10.4f%10.4f%10.4f\n", wave[i], obsFlux[i], obsErr[i], synFlux[i]);
        outFile << tempString;
    }
    outFile.close();
}

/// @brief Writes a spectrum with error bars to a file
/// @param fName 
/// @param wave 
/// @param flux 
/// @param err 
void writeSpec(std::string fName, std::vector<double> &wave, std::vector<double> &flux, std::vector<double> &err){
    std::ofstream outFile;
    outFile.open(fName);
    char tempString[128];

    for(int i = 0; i < wave.size(); i++){
        sprintf(tempString, "%10.4f%10.4f%10.4f\n", wave[i], flux[i], err[i]);
        outFile << tempString;
    }
    outFile.close();
}

/// @brief Cuts out a portion of a larger spectrum
/// @param wavLow lower bound of cut-out section
/// @param wavHi upper bound of cut-out section
/// @param inWave vector of wavelength data points to slice from
/// @param inFlux vector of flux data points to slice from
/// @param outWave empty vector to hold cropped wavelength data
/// @param outFlux empty vector to hold cropped flux data
void getObsSubset(double wavLow, double wavHi, std::vector<double> inWave, std::vector<double> inFlux, std::vector<double> &outWave, std::vector<double> &outFlux){
    outWave.clear();
    outFlux.clear();
    // Copies portion of larger spectra from wavLow to wavHi
    //to a separate array 

    for (int i = 0; i < inWave.size(); i++){
        if((inWave[i] >= wavLow) && (inWave[i] <= wavHi)){
            outWave.push_back(inWave[i]);
            outFlux.push_back(inFlux[i]);
        }
    }
}

/// @brief Cuts out a portion of a larger spectrum
/// @param wavLow lower bound of cut-out section
/// @param wavHi upper bound of cut-out section
/// @param inWave vector of wavelength data points to slice from
/// @param inFlux vector of flux data points to slice from
/// @param outWave empty vector to hold cropped wavelength data
/// @param outFlux empty vector to hold cropped flux data
void getObsSubset(double wavLow, double wavHi, std::vector<double> inWave, std::vector<double> inFlux, std::vector<double> inError,
                  std::vector<double> &outWave, std::vector<double> &outFlux, std::vector<double> &outError){

    // Copies portion of larger spectra from wavLow to wavHi
    //to a separate array 
    outWave.clear();
    outFlux.clear();
    outError.clear();
    for (int i = 0; i < inWave.size(); i++){
        if((inWave[i] >= wavLow) && (inWave[i] <= wavHi)){
            outWave.push_back(inWave[i]);
            outFlux.push_back(inFlux[i]);
            outError.push_back(inError[i]);

        }
    }
}


/// @brief Takes in observed wav,flx data and interpolates the synthetic spectrum to match the observed wavelength points
/// @param obsWav vector of observed (target) wavelengths
/// @param obsFlx vector of observed flux
/// @param synWav vector of wavelengths to be interpolated over
/// @param synFlx vector of flux points to be interpolated
/// @param synFlxOut empty vector to hold synthetic flux data, now matched to the wavelengths in 'obsWav'
void matchWavelengths(std::vector<double> &obsWav, std::vector<double> &synWav, std::vector<double> &synFlx, std::vector<double> &synFlxOut){

    // Takes in observed wav,flx data from a spectrograph and 
    // interpolates the synthetic flx data to the observed wavelength points
    // if(obsWav.size() == 0 || (synWav[0] < obsWav[0])){
    //     throw std::domain_error("ERROR: Interpolation will fail. Breaking to avoid GSL error handler");
    // }
    synFlxOut = interp1DWrapper(obsWav,synWav,synFlx);
}

/// @brief Updates MOOG par file with user-chosen values
/// @param MOOGParams The struct of parameters to be written
void updateParFile(parFileInfo updatedParams){
    //pUpdaterintf("updateParFile says offset is %.3f and vsini is %.3f\n", updatedParams.offsets[1], updatedParams.vsini);

    std::vector<std::string> parLines;

    std::ifstream inFile;
    std::string fLine;
    inFile.open(updatedParams.parFile);
    while(fLine.find("abundances") == std::string::npos){
        std::getline(inFile,fLine, '\n');
        parLines.push_back(fLine);
        if(inFile.eof()){
            throw std::runtime_error("Check that batch.par exists and is formatted correctly\n");
        };
    }
    inFile.close();
    parLines.pop_back();


    char tempString[128];
    sprintf(tempString, "abundances %6d    1", updatedParams.elements.size()  );
    parLines.push_back(tempString);

    for(int i = 0; i < updatedParams.elements.size(); i++){
        sprintf(tempString, "%5d%12.2f", updatedParams.elements[i], updatedParams.offsets[i]);
        parLines.push_back(tempString);
    }

    sprintf(tempString, "synlimits\n%10.2f%10.2f%10.2f%10.2f", updatedParams.lineWav - updatedParams.windowSize/2, updatedParams.lineWav + updatedParams.windowSize/2, 0.01, 20.0);
    parLines.push_back(tempString);

    sprintf(tempString, "plotpars       1\n%10.2f%10.2f%10.2f%10.2f\n 0 0 0 1.0", updatedParams.lineWav - updatedParams.windowSize/2 -0.01, updatedParams.lineWav + updatedParams.windowSize/2 + 0.01, 0.0, 1.3);
    parLines.push_back(tempString);

    sprintf(tempString, " r %.3f 0.0 0.0 %.3f 0.0", updatedParams.FWHM, updatedParams.vsini);
    parLines.push_back(tempString);

    std::ofstream outFile;
    outFile.open(updatedParams.parFile);
    for (int i = 0; i < parLines.size(); i++){
        if(parLines[i].length() > 0){outFile << parLines[i] << '\n';}
    }
    outFile.close();

}
