#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "line.h"



class atmosphere{
    public:
        std::string paramfileName;
        std::vector<std::string> elementString;
        std::string workDir = "";
        std::string specFile = "";
        std::vector<std::string> outDataDir;

        double Teff = INFINITY;
        double logg = INFINITY;
        double v_micro = INFINITY;
        double v_broad_init = INFINITY;
        double MonH = INFINITY;
        double AonM = INFINITY;

        std::vector<std::string> LD_Filters;
        std::vector<double> LD_Wavelengths;
        std::vector<double> LD_Coefficients;

        abundanceList customAbundances;
        bool useMolecules = false;
        
        atmosphere();
        void readFromFile(std::string filename);
        void updateFile();
        void updateStellarParameters(double newMetallicity);
};

atmosphere::atmosphere(){

}

void atmosphere::readFromFile(std::string fileName){

    paramfileName = fileName;

    std::ifstream inFile;
    std::string inLine;

    inFile.open(fileName);
    if(inFile.fail()){
        throw std::runtime_error("Atmosphere parameter file " + fileName + " not found");
    }

    while(!inFile.eof()){
        std::getline(inFile, inLine);
        if(inLine.substr(0,10).find("WorkDir") != std::string::npos){
            workDir = inLine.substr(11);
        }
        else if(inLine.substr(0,10).find("SpecFile") != std::string::npos){
            specFile = inLine.substr(11);
        }
        else if(inLine.substr(0,10).find("FitAtom") != std::string::npos){
            std::istringstream atomList(inLine.substr(10));
            std::string tempElement;
            int tempNum;
            std::string tempDirName;
            while(std::getline(atomList,tempElement, ',')){
                tempNum = std::stoi(tempElement);
                tempElement = std::to_string(tempNum);
                elementString.push_back(tempElement);
                customAbundances.addNewElement(std::stoi(tempElement));
                tempDirName = tempElement + "/";
                outDataDir.push_back(tempDirName);
                //std::cout << workDir + tempDirName << std::endl;
                if(std::filesystem::exists(workDir + tempDirName)){
                    std::filesystem::remove_all(workDir + tempDirName);
                }
                std::filesystem::create_directory(workDir + tempDirName);
            }
        }
        else if(inLine.substr(0,10).find("T_eff") != std::string::npos){
            Teff = std::stod(inLine.substr(10,10));
        }
        else if(inLine.substr(0,10).find("log(g)") != std::string::npos){
            logg = std::stod(inLine.substr(10,10));
        }
        else if(inLine.substr(0,10).find("v_micro") != std::string::npos){
            v_micro = std::stod(inLine.substr(10,10));
        }
        else if(inLine.substr(0,10).find("v_broad_0") != std::string::npos){
            v_broad_init = std::stod(inLine.substr(10,10));
        }
        else if(inLine.substr(0,10).find("[M/H]") != std::string::npos){
            MonH = std::stod(inLine.substr(10,10));
        }
        else if(inLine.substr(0,10).find("[alpha/M]") != std::string::npos){
            AonM = std::stod(inLine.substr(10,10));
        }
        else if(inLine.substr(0,10).find("LD_Coeff") != std::string::npos){
            LD_Wavelengths.push_back(std::stod(inLine.substr(20,10)));
            LD_Coefficients.push_back(std::stod(inLine.substr(30,10)));

        }

        else if(inLine.substr(0,10).find("SetAbund") != std::string::npos){
            customAbundances.addNewElement(std::stoi(inLine.substr(10,10)), std::stod(inLine.substr(20,10)));
        }
    }
    inFile.close();
    if(Teff < 4600){
        useMolecules = true;
    }

    if (outDataDir.size() == 0){
        outDataDir.push_back(elementString[0] + "/");
        std::filesystem::create_directory(workDir + outDataDir[0]);

    }
}
