#pragma once
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <filesystem>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "statistics.h"
namespace fs = std::filesystem;

std::string buildAtmoFileName(double Teff, double logg, double M_H, double alpha){
    
    std::stringstream filename;
    filename << "PHOENIX/atm/";

    if(M_H > 0){filename << "+";}
    else if (M_H == 0){filename << "-";}
    filename  << std::fixed << std::setprecision(1)<< M_H;

    filename << "/alpha";
    if(alpha > 0){filename << "+";}
    else if (alpha == 0){filename << "-";}
    filename  << std::fixed << std::setprecision(1)<< alpha;


    filename << "/lte";
    filename << std::setfill('0') << std::setw(5) << (int(Teff));

    if(logg >= 0){filename << "-";}
    else{filename << "+";}
    filename << std::fixed << std::setprecision(2) << logg;

    if(M_H > 0){filename << "+";}
    else if (M_H == 0){filename << "-";}
    filename  << std::fixed << std::setprecision(1)<< M_H;

    if(alpha != 0){
        filename << ".Alpha=";
        if(alpha > 0){filename << "+";}
        filename  << std::fixed << std::setprecision(2)<< alpha;

    }

    filename << ".PHOENIX-ACES-AGSS-COND-2011.ATMOS.txt";

    return filename.str();
}

std::vector<std::vector<double>> readAtmoFile(double Teff, double logg, double M_H, double alpha){
    std::vector<std::vector<double> > atmosphere;
    std::ifstream inFile;
    std::string fLine;
    std::vector<double> buffer;
    inFile.open(buildAtmoFileName(Teff,logg,M_H, alpha));
    if(!inFile.fail()){
        while(!inFile.eof()){
            std::getline(inFile, fLine, '\n');
            while(fLine.size() > 0){
                buffer.push_back(std::stof(fLine.substr(0, 13)));
                if (fLine.size() > 12){fLine = fLine.substr(13, fLine.size()-13);}
                else{fLine = "";}
                
            }
            if(buffer.size() > 0){atmosphere.push_back(buffer);}
            buffer.clear();
        }
        inFile.close();
    }

    else{
        throw std::runtime_error( "Model atmosphere file " + buildAtmoFileName(Teff, logg, M_H, alpha) + " not found");
    }

    return atmosphere;
}

std::vector<double> getColOf2DVector(std::vector<std::vector<double>> &vector, int colNum){
    std::vector<double> column;
    for(size_t i = 0; i < vector.size(); i++){
        column.push_back(vector[i][colNum]);
    }
    return column;
}

double interpPointOn3DGrid(double x, double y, double z, std::vector<std::vector<std::vector<double> > > &grid){
    

    double c00 = grid[0][0][0]*(1-x) + grid[1][0][0]*x;
    double c01 = grid[0][0][1]*(1-x) + grid[1][0][1]*x;
    double c10 = grid[0][1][0]*(1-x) + grid[1][1][0]*x;
    double c11 = grid[0][1][1]*(1-x) + grid[1][1][1]*x;
    
    double c0 = c00 * (1-y) + c10*y;
    double c1 = c01 * (1-y) + c11*y;
    double c = c0 * (1-z) + c1*z;
    return c;
}

double interpPointOn4DGrid(double x, double y, double z, double a, std::vector<std::vector<std::vector<double> > > &grid0, std::vector<std::vector<std::vector<double> > > &grid1){
    
    // 3D Interpolation on grid where a = 0
    double c000 = grid0[0][0][0]*(1-x) + grid0[1][0][0]*x;
    double c010 = grid0[0][0][1]*(1-x) + grid0[1][0][1]*x;
    double c100 = grid0[0][1][0]*(1-x) + grid0[1][1][0]*x;
    double c110 = grid0[0][1][1]*(1-x) + grid0[1][1][1]*x;
    
    double c00 = c000 * (1-y) + c100*y;
    double c10 = c010 * (1-y) + c110*y;
    double c0 = c00 * (1-z) + c10*z;

    // 3D Interpolation on grid where a = 1
    double c001 = grid1[0][0][0]*(1-x) + grid1[1][0][0]*x;
    double c011 = grid1[0][0][1]*(1-x) + grid1[1][0][1]*x;
    double c101 = grid1[0][1][0]*(1-x) + grid1[1][1][0]*x;
    double c111 = grid1[0][1][1]*(1-x) + grid1[1][1][1]*x;
    
    double c01 = c001 * (1-y) + c101*y;
    double c11 = c011 * (1-y) + c111*y;
    double c1 = c01 * (1-z) + c11*z;
    

    // Return final dimension of interpolation
    return c0 * (1-a) + c1*a;
}

void interpAtmosphere(double Teff, double logg, double vmic, double M_H, double alpha){

    std::vector<std::vector<double> > model;

    bool noInterpT = false;
    bool noInterpG = false;
    bool noInterpM = false;
    bool noInterpA = false;
    // Set available values on the grid to interpolate
    // Check if we even need to interpolate or just
    // pluck a model from the grid
    std::vector<int> availTeff;
    for(int i = 0; i < 72; i++){
        availTeff.push_back(2600+(100*i));
        if(std::abs(2600+(100*i)-Teff) < 0.5){
            noInterpT = true;
        }
    }
    std::vector<float> availLogg;
    for(int i = 0; i < 14; i++){
        availLogg.push_back(-0.5+(0.5*i));
        if(std::abs(-0.5+(0.5*i)-logg) < 0.005){
            noInterpG = true;
        }
    }
    std::vector<float> availMetal = {-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0};
    for(size_t i = 0; i < availMetal.size(); i++){
        if(std::abs(availMetal[i]-M_H) < 0.005){
            noInterpM = true;
        }
    }


    std::vector<double> availAlpha;
    if(M_H <= 0){
        availAlpha = {-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2};
    }
    else{
        availAlpha = {-0.4,0.0,0.6,1.2};
    }
    for(size_t i = 0; i < availAlpha.size(); i++){
        if(std::abs(availAlpha[i]-alpha) < 0.005){
            noInterpA = true;
        }
    }


    //std::cout << noInterpT << noInterpG << noInterpM << "\n";
    //std::cout << (noInterpT && noInterpG && noInterpM) << "\n";

    // If no interpolation needed, get the file. Else interpolate
    if(noInterpT && noInterpG && noInterpM && noInterpA){
        model = readAtmoFile(Teff, logg, M_H, alpha);
    }
    else{
        std::vector<std::vector<std::vector<double> > > grid0 = {{{0,0},{0,0}},{{0,0},{0,0}}};
        std::vector<std::vector<double> > model10;
        std::vector<std::vector<double> > model20;
        std::vector<std::vector<double> > model30;
        std::vector<std::vector<double> > model40;
        std::vector<std::vector<double> > model50;
        std::vector<std::vector<double> > model60;
        std::vector<std::vector<double> > model70;
        std::vector<std::vector<double> > model80;
        std::vector<double> tauross;
        std::vector<double> tauross10;
        std::vector<double> tauross20;
        std::vector<double> tauross30;
        std::vector<double> tauross40;
        std::vector<double> tauross50;
        std::vector<double> tauross60;
        std::vector<double> tauross70;
        std::vector<double> tauross80;

        std::vector<std::vector<std::vector<double> > > grid1 = {{{0,0},{0,0}},{{0,0},{0,0}}};
        std::vector<std::vector<double> > model11;
        std::vector<std::vector<double> > model21;
        std::vector<std::vector<double> > model31;
        std::vector<std::vector<double> > model41;
        std::vector<std::vector<double> > model51;
        std::vector<std::vector<double> > model61;
        std::vector<std::vector<double> > model71;
        std::vector<std::vector<double> > model81;
        std::vector<double> tauross11;
        std::vector<double> tauross21;
        std::vector<double> tauross31;
        std::vector<double> tauross41;
        std::vector<double> tauross51;
        std::vector<double> tauross61;
        std::vector<double> tauross71;
        std::vector<double> tauross81;



        int tm1 = 0;
        int tp1 = 0;
        double gm1 = 0;
        double gp1 = 0;
        double mm1 = 0;
        double mp1 = 0;
        double am1 = 0;
        double ap1 = 0;

        for(size_t i = 0; i < availTeff.size() - 1; i++){
            if(Teff == availTeff[i]){
                tm1 = availTeff[i];
                tp1 = availTeff[i];
                break;
            }
            if((Teff - availTeff[i]) < (availTeff[i+1] - availTeff[i])){
                tm1 = availTeff[i];
                tp1 = availTeff[i+1];
                break;
            }
        }

        for(size_t i = 0; i < availLogg.size() - 1; i++){
            if(logg == availLogg[i]){
                gm1 = availLogg[i];
                gp1 = availLogg[i];
                break;
            }
            if((logg - availLogg[i]) < (availLogg[i+1] - availLogg[i])){
                gm1 = availLogg[i];
                gp1 = availLogg[i+1];
                break;
            }
        }
        for(size_t i = 0; i < availMetal.size() - 1; i++){
            if(M_H == availMetal[i]){
                mm1 = availMetal[i];
                mp1 = availMetal[i];
                break;
            }
            if((M_H - availMetal[i]) < (availMetal[i+1] - availMetal[i])){
                mm1 = availMetal[i];
                mp1 = availMetal[i+1];
                break;
            }
        }

        for(size_t i = 0; i < availAlpha.size() - 1; i++){
            if(alpha == availAlpha[i]){
                am1 = availAlpha[i];
                ap1 = availAlpha[i];
                break;
            }
            if((alpha - availAlpha[i]) < (availAlpha[i+1] - availAlpha[i])){
                am1 = availAlpha[i];
                ap1 = availAlpha[i+1];
                break;
            }
        }

        for(int i = 1; i < 9; i++){
            switch (i){
                case 1:
                    model = readAtmoFile(tm1,gm1,mm1,am1);
                    break;
                case 2:
                    model = readAtmoFile(tm1,gm1,mp1,am1);
                    break;
                case 3:
                    model = readAtmoFile(tm1,gp1,mm1,am1);
                    break;
                case 4:
                    model = readAtmoFile(tm1,gp1,mp1,am1);
                    break;
                case 5:
                    model = readAtmoFile(tp1,gm1,mm1,am1);
                    break;
                case 6:
                    model = readAtmoFile(tp1,gm1,mp1,am1);
                    break;
                case 7:
                    model = readAtmoFile(tp1,gp1,mm1,am1);
                    break;
                case 8:
                    model = readAtmoFile(tp1,gp1,mp1,am1);
                    break;
            }
            tauross = getColOf2DVector(model, 0);
                        switch (i){
                case 1:
                    model10 = model;
                    tauross10 = tauross;
                    break;
                case 2:
                    model20 = model;
                    tauross20 = tauross;
                    break;
                case 3:
                    model30 = model;
                    tauross30 = tauross;
                    break;
                case 4:
                    model40 = model;
                    tauross40 = tauross;
                    break;
                case 5:
                    model50 = model;
                    tauross50 = tauross;
                    break;
                case 6:
                    model60 = model;
                    tauross60 = tauross;
                    break;
                case 7:
                    model70 = model;
                    tauross70 = tauross;
                    break;
                case 8:
                    model80 = model;
                    tauross80 = tauross;
                    break;
            }

            //double bot_tauross = std::min_element(backTaus.begin(), backTaus.end());
        }

        for(int i = 1; i < 9; i++){
            switch (i){
                case 1:
                    model = readAtmoFile(tm1,gm1,mm1,ap1);
                    break;
                case 2:
                    model = readAtmoFile(tm1,gm1,mp1,ap1);
                    break;
                case 3:
                    model = readAtmoFile(tm1,gp1,mm1,ap1);
                    break;
                case 4:
                    model = readAtmoFile(tm1,gp1,mp1,ap1);
                    break;
                case 5:
                    model = readAtmoFile(tp1,gm1,mm1,ap1);
                    break;
                case 6:
                    model = readAtmoFile(tp1,gm1,mp1,ap1);
                    break;
                case 7:
                    model = readAtmoFile(tp1,gp1,mm1,ap1);
                    break;
                case 8:
                    model = readAtmoFile(tp1,gp1,mp1,ap1);
                    break;
            }
            tauross = getColOf2DVector(model, 0);
                        switch (i){
                case 1:
                    model11 = model;
                    tauross11 = tauross;
                    break;
                case 2:
                    model21 = model;
                    tauross21 = tauross;
                    break;
                case 3:
                    model31 = model;
                    tauross31 = tauross;
                    break;
                case 4:
                    model41 = model;
                    tauross41 = tauross;
                    break;
                case 5:
                    model51 = model;
                    tauross51 = tauross;
                    break;
                case 6:
                    model61 = model;
                    tauross61 = tauross;
                    break;
                case 7:
                    model71 = model;
                    tauross71 = tauross;
                    break;
                case 8:
                    model81 = model;
                    tauross81 = tauross;
                    break;
            }

            //double bot_tauross = std::min_element(backTaus.begin(), backTaus.end());
        }

        tauross = tauross10;
        std::vector<double> backTaus = {tauross10.back(),tauross20.back(),tauross30.back(),tauross40.back(),tauross50.back(),tauross60.back(),tauross70.back(),tauross80.back(), tauross11.back(),tauross21.back(),tauross31.back(),tauross41.back(),tauross51.back(),tauross61.back(),tauross71.back(),tauross81.back()};
        std::vector<double> frontTaus = {tauross10.front(),tauross20.front(),tauross30.front(),tauross40.front(),tauross50.front(),tauross60.front(),tauross70.front(),tauross80.front(), tauross11.front(),tauross21.front(),tauross31.front(),tauross41.front(),tauross51.front(),tauross61.front(),tauross71.front(),tauross81.front()};

        double botTauRoss = *std::min_element(backTaus.begin(), backTaus.end());
        double topTauRoss = *std::max_element(frontTaus.begin(), frontTaus.end());

        std::vector<double> tauRossNew;
        std::vector<double> layer;
        for(size_t i = 0; i < tauross.size(); i++){
            layer.push_back(i);
        }

        std::vector<double> interpX = layer;
        std::vector<double> interpY = tauross;
        while(interpY[0] < topTauRoss){
            interpY.erase(interpY.begin());
            interpX.erase(interpX.begin());
        }
        while(interpY[interpY.size()-1] > botTauRoss){
            interpY.pop_back();
            interpX.pop_back();
        }

        tauRossNew = interp1DWrapper(layer, interpX,interpY);


        double mapT = 0;
        double mapG = 0;
        double mapM = 0;
        double mapA = 0;
        if(tp1 != tm1){mapT = (Teff-tm1)/(tp1-tm1);}
        if(gp1 != gm1){mapG = (logg-gm1)/(gp1-gm1);}
        if(mp1 != mm1){mapM = (M_H - mm1)/(mp1-mm1);}
        if(ap1 != am1){mapA = (alpha - am1)/(ap1-am1);}









        for(size_t i = 0; i < layer.size(); i++){
            for(size_t j = 0; j < model[0].size(); j++){
                grid0[0][0][0] = interp1DWrapper({tauRossNew[i]}, tauross10,getColOf2DVector(model10,j))[0];
                grid0[0][0][1] = interp1DWrapper({tauRossNew[i]}, tauross20,getColOf2DVector(model20,j))[0];
                grid0[0][1][0] = interp1DWrapper({tauRossNew[i]}, tauross30,getColOf2DVector(model30,j))[0];
                grid0[0][1][1] = interp1DWrapper({tauRossNew[i]}, tauross40,getColOf2DVector(model40,j))[0];
                grid0[1][0][0] = interp1DWrapper({tauRossNew[i]}, tauross50,getColOf2DVector(model50,j))[0];
                grid0[1][0][1] = interp1DWrapper({tauRossNew[i]}, tauross60,getColOf2DVector(model60,j))[0];
                grid0[1][1][0] = interp1DWrapper({tauRossNew[i]}, tauross70,getColOf2DVector(model70,j))[0];
                grid0[1][1][1] = interp1DWrapper({tauRossNew[i]}, tauross80,getColOf2DVector(model80,j))[0];

                grid1[0][0][0] = interp1DWrapper({tauRossNew[i]}, tauross11,getColOf2DVector(model11,j))[0];
                grid1[0][0][1] = interp1DWrapper({tauRossNew[i]}, tauross21,getColOf2DVector(model21,j))[0];
                grid1[0][1][0] = interp1DWrapper({tauRossNew[i]}, tauross31,getColOf2DVector(model31,j))[0];
                grid1[0][1][1] = interp1DWrapper({tauRossNew[i]}, tauross41,getColOf2DVector(model41,j))[0];
                grid1[1][0][0] = interp1DWrapper({tauRossNew[i]}, tauross51,getColOf2DVector(model51,j))[0];
                grid1[1][0][1] = interp1DWrapper({tauRossNew[i]}, tauross61,getColOf2DVector(model61,j))[0];
                grid1[1][1][0] = interp1DWrapper({tauRossNew[i]}, tauross71,getColOf2DVector(model71,j))[0];
                grid1[1][1][1] = interp1DWrapper({tauRossNew[i]}, tauross81,getColOf2DVector(model81,j))[0];

                model[i][j] = interpPointOn4DGrid(mapT,mapG,mapM,mapA, grid0, grid1);
            }
        }

        //throw std::runtime_error( "Atmosphere interpolation not implemented yet. Choose values from grid.");

    }











    /* By here, interpolation is finished. Output atmosphere to a file*/


    // Abundances for editing based on alpha enhancement
    double solarMOOG[96] = {0,
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



    std::ofstream outFile;
    std::string fLine;
    std::vector<double> buffer;
    outFile.open("mpplib/modelatmosphere.txt");
    // Write the atmosphere to file    
    outFile << "GENERIC\n";
    outFile << "TEFF=" << Teff << "  GRAVITY=" << logg << " LTE \n";
    outFile << "NTAU          64\n";
    outFile << "5000\n";
    outFile << std::setprecision(4);

    for(size_t i = 0; i < model.size(); i++){
        //std::cout << "\nROW " << i <<":\nCOL: ";
        for(int j = 0; j < 4; j++){
            //std::cout << j;
            outFile << std::scientific << std::setw(12);
            outFile << model[i][j];
        }
        outFile << "\n";
    }
    outFile << std::fixed << std::setprecision(2) << std::setw(4);
    outFile << vmic << "\n";
    outFile << "NATOMS         8	" << M_H << "\n";
    for(int i = 8; i < 23; i+=2){
        outFile << i << " " << solarMOOG[i] + M_H + alpha << "    ";
    }
    outFile << "\n";
    outFile << "NMOL           20\n";
    outFile << "606.0 106.0 607.0 608.0 107.0 108.0 112.0 707.0 708.0 808.0 12.1 60808.0 10108.0 101.0 6.1 7.1 8.1 822.0 22.1 126.0\n";
    outFile.close();
}
