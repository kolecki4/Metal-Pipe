#include <iostream>
#include <vector>

std::string getDefaultLinelistName(int element){
    std::ifstream inFile;
    std::string inLine;
    
    std::string fileName;
    inFile.open("linelists/defaults.txt");
    if(inFile.fail()){
        throw std::runtime_error("Default line lists file not found!");
    }
    while(!inFile.eof()){
        std::getline(inFile, inLine);
        if (std::stoi(inLine.substr(0,2)) == element){
            fileName = inLine.substr(4);
            inFile.close();
            return fileName;
        }
    }
    if (fileName.size() == 0){
        throw std::runtime_error("No default line list specified for element " + std::to_string(element));
    }
    inFile.close();
    return fileName;
}
