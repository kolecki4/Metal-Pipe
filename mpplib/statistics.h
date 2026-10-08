#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

std::vector<double> movingAverageSlope(std::vector<double> &data, int range){
    std::vector<double> slope;
    double m = 0;

    for(size_t i = 0; i < data.size()-range-1; i++){
        
        // Take average slope over next range points
        m = 0;
        for(int j = 1; j < range; j++){
            m += (data[i+j]-data[i])/j;
        }
        m /= range;

        slope.push_back(m);
    }

    while(slope.size() < data.size()){
        slope.push_back(0);
    }

    return slope;
}


std::vector<double> interp1DWrapper(std::vector<double> newXPoints, std::vector<double> xData, std::vector<double> yData){
    std::vector<double> newYPoints;
    double m = 0;

    if((xData.size() < 3) || (yData.size() < 3) || (newXPoints.size() == 0)){
        throw std::runtime_error("No data here to interpolate\n");
    }

    for(size_t i = 0; i < xData.size()-1; i++){
        if(xData[i+1]-xData[i] < 0){throw std::runtime_error("X points must be strictly increasing\n");}
    }

    gsl_interp_accel *accel_ptr = gsl_interp_accel_alloc();
    gsl_spline *spline_ptr = gsl_spline_alloc(gsl_interp_cspline, xData.size());
    gsl_spline_init(spline_ptr, xData.data(), yData.data(), xData.size());
    


    for (size_t i = 0; i < newXPoints.size(); i++){
        // If outside interp range, linearly extrapolate
        if(newXPoints[i] < xData[0]){
            m = (yData[1]-yData[0])/(xData[1]-xData[0]);
            newYPoints.push_back(1);
        }
        else if(newXPoints[i] > xData[xData.size()-1]){
            m = (yData[xData.size()-1]-yData[xData.size()-2])/(xData[xData.size()-1]-xData[xData.size()-2]);
            newYPoints.push_back(1);
        }
        
        else{newYPoints.push_back(gsl_spline_eval(spline_ptr, newXPoints[i], accel_ptr));}
    }

    gsl_spline_free(spline_ptr);
    gsl_interp_accel_free(accel_ptr);
    return newYPoints;
}


double crossCorlkjhrelateDelay(std::vector<double> obsWav,std::vector<double> obsFlx,std::vector<double> synWav,std::vector<double> synFlx){
    
    //TODO: Try interpolating observed data to synthetic wavelength grid
    std::vector<double> obsInterpd = interp1DWrapper(synWav,obsWav, obsFlx);
    std::vector<double> obsInterpdOffset = obsInterpd;

    std::vector<double> weights;
    for(size_t i = 0; i < obsInterpd.size(); i++){
        weights.push_back(std::min(double(i),double(obsInterpd.size()-i)));
    }

    double delay = 0;
    double maxCC = 0;
    double sum;
    for(double i = -.05; i < .05; i+=0.0005){
        sum = 0;

        for(size_t j = 0; j < obsWav.size(); j++){
            obsWav[j] += i;
        }
        obsInterpd = interp1DWrapper(synWav,obsWav, obsFlx);
        for(size_t j = 0; j < synWav.size(); j++){
            sum += ((1-synFlx[j])*(1-obsInterpd[j]));//*weights[j];
            
        }
        if(sum > maxCC){
            delay = i;
            maxCC = sum;
        }

        for(size_t j = 0; j < obsWav.size(); j++){
            obsWav[j] -= i;
        }
    }
    //std::cout << delay << "\n";
    return delay;
}

double reducedChiSq(std::vector<double> &obs, std::vector<double> &syn, std::vector<double> &err){
    double chiSq = 0;
    double DOF = obs.size() - 3.;

    for(size_t i = 0; i < obs.size(); i++){
        //printf("%.3f %.3f %.3f\n", obs[i], syn[i], err[i]);
        chiSq += ((pow(obs[i] - syn[i],2))/ pow(err[i],2));
    }

    return std::min(chiSq/ DOF, 999999.);
}

double reducedChiSq(std::vector<double> &obs, std::vector<double> &syn, std::vector<double> &err, std::vector<double> &weights){
    double chiSq = 0;
    double DOF = obs.size() - 3.;

    for(size_t i = 0; i < weights.size(); i++){
        if (weights[i] < 0){weights[i] = 0;}
    }

    for(size_t i = 0; i < obs.size(); i++){
        //printf("%.3f %.3f %.3f\n", obs[i], syn[i], err[i]);
        chiSq += ((pow(obs[i] - syn[i],2))/ pow(err[i],2))*weights[i];
    }

    return chiSq/ DOF;
}




double calculateResolutionAtLine(std::vector<double> &inWave, double lineCenter, double windowSize){
    std::vector<double> outWave;

    for (size_t i = 0; i < inWave.size(); i++){
        if((inWave[i] >= (lineCenter - windowSize/2.)) && (inWave[i] <= (lineCenter + windowSize/2.))){
            outWave.push_back(inWave[i]);
        }
    }

    double avgDifference = 0;
    for(size_t i = 1; i < outWave.size(); i++){
        avgDifference += (outWave[i] - outWave[i-1]);
    }
    avgDifference /= outWave.size()-1;

    return lineCenter/avgDifference;
}

double multNormFactor(std::vector<double> &obs, std::vector<double> &syn, std::vector<double> &err){
    std::vector<double> continuumPointsObs;
    std::vector<double> continuumPointsSyn;
    std::vector<double> continuumPointsErr;
    std::vector<double> ObsFluxSlope;
    for(size_t i = 0; i < obs.size() -3; i++){
        ObsFluxSlope.push_back(((obs[i+1] - obs[i]) + (obs[i+2] - obs[i+1]))/2);
    }
    ObsFluxSlope.push_back(0);
    ObsFluxSlope.push_back(0);
    // Grab continuum points from the spectra
    for(size_t i = 0; i < obs.size(); i++){
        if( abs(ObsFluxSlope[i]) < 555555.001 && (syn[i] > 0.95 && syn[i] < 1.01) && (obs[i] > 0.95 && obs[i] < 1.1)) {
            continuumPointsObs.push_back(obs[i]);
            continuumPointsSyn.push_back(syn[i]);
            continuumPointsErr.push_back(err[i]);
        }
    }
    //printf("Size of continuum region: %i", continuumPointsErr.size());
    // If we end up with a chi^2 of inf or negative, give up
    // DoF = 3 so we need a size of at least 4
    if(continuumPointsObs.size() < 4){
        //throw std::length_error("ERROR: Length of continuum too short");
        return 1.;
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
    return optimalI;
}

double median(std::vector<double> list){

    if(!(list.size() > 0)){
        return(0);
    }
    std::sort(list.begin(), list.end());

    if (list.size() % 2 ==0){
        return list[list.size()/2];
    }
    else{
        return 0.5*(list[list.size()/2] + list[list.size()/2 + 1]);
    }
}

double median(std::vector<double> list, std::vector<double> filter, double cutoff){
    if(!(list.size() > 0)){
        return(0);
    }
    
    if(list.size() != filter.size()){
        std::length_error("uh oh\n");
    }

    std::vector<double> filteredList;
    for(size_t i = 0; i < list.size(); i++){
        if(filter[i] < cutoff){

            filteredList.push_back(list[i]);
        }

    }
    

    std::sort(filteredList.begin(), filteredList.end());

    if (filteredList.size() % 2 ==1){
        return filteredList[filteredList.size()/2];
    }
    else{
        return 0.5*(filteredList[filteredList.size()/2 - 1] + filteredList[filteredList.size()/2]);
    }
}

double stdevMedian(std::vector<double> list){
    double med = median(list);
    double stdev = 0;
    for(size_t i = 0; i < list.size(); i++){
        stdev += pow((list[i] - med), 2)/list.size();
    }

    stdev = pow(stdev/list.size(), 0.5);
    return stdev; 
}

double stdevMedian(std::vector<double> list, std::vector<double> filter, double cutoff){


    std::vector<double> filteredList;
    for(size_t i = 0; i < list.size(); i++){
        if(filter[i] < cutoff){

            filteredList.push_back(list[i]);
        }

    }

    double med = median(filteredList);
    double stdev = 0;


    for(size_t i = 0; i < filteredList.size(); i++){
        stdev += pow((filteredList[i] - med), 2)/filteredList.size();
    }

    stdev = pow(stdev/list.size(), 0.5);
    return stdev;
}

double stdevMean(std::vector<double> list){
    double mean = std::accumulate(list.begin(), list.end(),0)/list.size();
    double stdev = 0;
    for(size_t i = 0; i < list.size(); i++){
        stdev += pow((list[i] - mean), 2)/list.size();
    }

    stdev = pow(stdev, 0.5);
    return stdev; 
}

double stdevMean(std::vector<double> list, std::vector<double> filter, double cutoff){


    std::vector<double> filteredList;
    for(size_t i = 0; i < list.size(); i++){
        if(filter[i] < cutoff){

            filteredList.push_back(list[i]);
        }

    }

    double mean = std::accumulate(filteredList.begin(), filteredList.end(),0)/filteredList.size();
    double stdev = 0;


    for(size_t i = 0; i < filteredList.size(); i++){
        stdev += pow((filteredList[i] - mean), 2)/filteredList.size();
    }

    stdev = pow(stdev, 0.5);
    return stdev;
}
