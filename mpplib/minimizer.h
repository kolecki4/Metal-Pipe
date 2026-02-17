#pragma once
#include "line.h"
#include "atmosphere.h"
#include <iostream>
#include <cmath>
#include <mutex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

double functionToMinimize(const gsl_vector *v, void *p){
    line *lineToFit = static_cast<line*>(p);

    double newOffset = gsl_vector_get(v,0);
    double newVBroad = gsl_vector_get(v,1);

    if(newVBroad < 0 || newVBroad > lineToFit->lineInfo.maxAllowedVBroad || newOffset > 1 || newOffset < -1){
        return 999999;
    }

    lineToFit->lineInfo.abundanceOffsets.updateElement(lineToFit->lineInfo.species, newOffset);
    lineToFit->lineInfo.v_broad = newVBroad;

    lineToFit->writeParFile();
    lineToFit->MOOGitUp();
    lineToFit->readSynthesizedLine(lineToFit->smoothedOutFile);
    lineToFit->interpGridsToEachOther();

    double chi2 = lineToFit->synGridChi2();
    
    
    if(std::isfinite(chi2)){return chi2;}
    else{return 999999;}
}

std::vector<double> chi2MinByAbundance(line &lineToFit){
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function F;
    F.n = 2;
    F.f = functionToMinimize;
    F.params = &lineToFit;

    // x=startingpoint
    x = gsl_vector_alloc(2);
    gsl_vector_set(x,0,0);
    gsl_vector_set(x,1,2);

    // ss = stepsize
    ss = gsl_vector_alloc(2);
    gsl_vector_set(ss,0,0.1);
    gsl_vector_set(ss,1,0.5);

    s = gsl_multimin_fminimizer_alloc(T,2);
    gsl_multimin_fminimizer_set(s, &F, x,ss);

    int iter = 0;
    int max_iter = 100;
    int status;
    double size = 1e8;
    double convergeSize = .05;
    do{


        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, convergeSize);
        //if (status == GSL_SUCCESS)
        //    printf ("Converged:\n");
        //printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",iter,gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1),s->fval, size);
        //printf("Iteration %i/%i for lambda = %.3f.\nSize = %.3f, need %.3f to converge\n", iter, max_iter, lineToFit.lineInfo.centralWavelength, size, convergeSize);
        //printf("gslmin says offset is %.3f and vsini is %.3f\n", gsl_vector_get(s->x,0), gsl_vector_get(s->x,1));

    }while(status == GSL_CONTINUE && iter < max_iter);
    if(size >convergeSize){
        throw std::runtime_error("Chi^2 minimization simplex did not reach size threshold for convergence");
    }
    double offset = gsl_vector_get(s->x,0);
    double  vsini = gsl_vector_get(s->x,1);
    double  chisq = s->fval;
    printf("%.3f: [X/Fe] = %6.3f; v_broad = %5.3f km/s\n", lineToFit.lineInfo.centralWavelength,offset, vsini);

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);

    return {offset,vsini,chisq};
}


int fitLine(line &lineToFit, atmosphere &atmosphereToFit, std::vector<double> &abundances, std::vector<double> &vbroads, std::vector<double> &chi2s, std::string outDataDir){
    std::mutex arrayWrite;
    bool badLine = false;

    // Find what the resolution of the spectrum is
    try{lineToFit.cutOutObservedLine(3);}
    catch(const std::runtime_error &e){
        badLine = true;
        std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
        return 1;
    }
    lineToFit.lineInfo.instBroadFWHM = lineToFit.lineInfo.instBroadWidthPixels*lineToFit.lineInfo.centralWavelength/ lineToFit.obsWaveGrid.avgResolution();

    // Find out how wide the synthesis region needs to be
    lineToFit.interpLDCoeff(atmosphereToFit.LD_Wavelengths,atmosphereToFit.LD_Coefficients);
    try{lineToFit.calculateFitRegions();}
    catch(const std::runtime_error &e){
        badLine = true;
        std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
        return 1;
    }
    // Set the width of synthesis
    lineToFit.lineInfo.widthOfSynthesis = lineToFit.getFitRanges()[1];

    // Set up variables for fitting
    if(lineToFit.lineInfo.widthOfSynthesis > 0){    
        // Set vsini
        lineToFit.lineInfo.v_broad = atmosphereToFit.v_broad_init;
        

        // Write info to par file, run linemake, then run MOOG
        lineToFit.writeParFile();
        try{lineToFit.setLineList();}
        catch(const std::runtime_error &e){
            badLine = true;
            std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
            return 1;
        }
        lineToFit.MOOGitUp();
        
        // Get the line data from the complete stellar spectrum
        try{lineToFit.cutOutObservedLine(lineToFit.lineInfo.widthOfSynthesis);}
        catch(const std::runtime_error &e){
            badLine = true;
            std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
            return 1;
        }

        // Get the synthesized line data from MOOG
        lineToFit.readSynthesizedLine(lineToFit.smoothedOutFile);


        // Set wavelength grids, renormalize spectra, calculate point fit weights 
        badLine = false;
        try{
            lineToFit.interpGridsToEachOther();
            lineToFit.renormalizeObs();
            lineToFit.crossCorrelateObs();
            lineToFit.setWeightsForChi2(atmosphereToFit.v_broad_init);
        }
        catch(const std::runtime_error &e){
            badLine = true;
            std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
            return 1;
        }

        /* THE LINE IS READY TO BE FIT */
        if(!badLine){

            // Fit the line
            try{
                std::vector<double> results = chi2MinByAbundance(lineToFit);
                lineToFit.lineInfo.abundanceOffsets.updateElement(lineToFit.lineInfo.species, results[0]);
                lineToFit.lineInfo.v_broad = results[1];
            }
        catch(const std::runtime_error &e){
            badLine = true;
            std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
            return 1;
        }

            /* THE LINE HAS BEEN FIT. OUTPUT RESULTING SPECTRUM TO FILE*/
            
            lineToFit.writeParFile();
            lineToFit.MOOGitUp();

            // Get the line data from the complete stellar spectrum
            try{lineToFit.cutOutObservedLine(lineToFit.lineInfo.widthOfSynthesis);}
            catch(const std::runtime_error &e){
                badLine = true;
                std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
                return 1;
            }
            // Get the synthesized line data from MOOG
            lineToFit.readSynthesizedLine(lineToFit.smoothedOutFile);


            // Make sure everything is accurate and then output the results
            try{
                lineToFit.interpGridsToEachOther();
                lineToFit.renormalizeObs();
                lineToFit.crossCorrelateObs();
                lineToFit.setWeightsForChi2(atmosphereToFit.v_broad_init);

                double minFlx = 99;
                double maxFlx = 0;
                std::vector<double> Flx = lineToFit.synthWaveGrid.getColumn("synthedFlux");
                for(size_t k = 0; k < Flx.size(); k++){
                    if(maxFlx < Flx[k]){maxFlx = Flx[k];}
                    if(minFlx > Flx[k]){minFlx = Flx[k];}
                }

                if(lineToFit.synGridChi2() < lineToFit.lineInfo.maxAllowedChi2 && minFlx > 0. && minFlx < 0.999 && lineToFit.lineInfo.v_broad < lineToFit.lineInfo.maxAllowedVBroad && !badLine){
                    arrayWrite.lock();
                    // If we're fitting an alpha element, we need to remember to add the alpha enhancement
                    if(lineToFit.lineInfo.species % 2 == 0 && lineToFit.lineInfo.species > 7 && lineToFit.lineInfo.species < 23){
                        lineToFit.outputBestFitAbundance(atmosphereToFit.workDir + outDataDir + "abundanceSummary.txt", atmosphereToFit.MonH + atmosphereToFit.AonM);
                        abundances.push_back(lineToFit.lineInfo.abundanceOffsets.getElement(lineToFit.lineInfo.species, "num").offset + atmosphereToFit.MonH + atmosphereToFit.AonM);
                    }
                    else{
                        lineToFit.outputBestFitAbundance(atmosphereToFit.workDir + outDataDir + "abundanceSummary.txt", atmosphereToFit.MonH);
                        abundances.push_back(lineToFit.lineInfo.abundanceOffsets.getElement(lineToFit.lineInfo.species, "num").offset + atmosphereToFit.MonH);
                    }
                    vbroads.push_back(lineToFit.lineInfo.v_broad);
                    chi2s.push_back(lineToFit.synGridChi2());
                    lineToFit.synthWaveGrid.writeToTextFile(atmosphereToFit.workDir + outDataDir + std::to_string(std::round(lineToFit.lineInfo.centralWavelength*100)/100));
                    arrayWrite.unlock();

                }
            }
            catch(const std::runtime_error &e){
                badLine = true;
                std::cout << "Bad line " << lineToFit.lineInfo.centralWavelength << " (" << e.what() <<")" << "\n" ;
                return 1;
            }
            // Write the resulting spectral data to files

        }
    }
    return 0;
}
