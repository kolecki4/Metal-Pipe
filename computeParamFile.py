import numpy as np
import sys
import os
from mpplib.newIsoLib import *
def setParams(starName, workDir, outFile, metal = 0, alpha = 0, runNum = 0):
    bands = "U,B,V,R,I,J,H,K,L,M".split(",")    
    filterWaves = "3650,4450,5510,6580,8060,12200,16300,21900,34500,47500".split(',')
    filterWaves = [float(i) for i in filterWaves]

    targ = starName
    outputFolder = workDir
    outputFile = outFile
    elements = np.array([26,20,22,12,14,6,7,8,16,11,13,19])
    #elements = np.array([26,6])

    specFile = ""
    files = os.listdir(outputFolder)
    for file in files:
        if "flattened" in file and ".txt" in file:
            specFile = file


    initOffsets = np.zeros_like(elements)
    initErrors = np.zeros_like(elements)

    MonHerror = 0
    AonMerror = 0

    lastMetallicity = 0
    lastAlpha = 0

    if runNum > 0:

        solarMOOG = [0,
            12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
             6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
             3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
             3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
             1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
             1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
            -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
             0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
             0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
            -5.00,-0.54,-5.00,-5.00,-5.00]



        lastMetallicity = readInputMetallicity(outputFolder+outputFile)
        lastAlpha= readInputAlphaEnhancement(outputFolder+outputFile)
        alpha = lastAlpha
        currentMetallicity = 0
        currentAlpha = 0

        elements, abundances, errors = readAbundances(outputFolder, outputFile, elements)
        #for i in range(len(initOffsets)):
        initOffsets = [abundances[i] - solarMOOG[int(elements[i])] for i in range(len(elements))]
        initErrors = [errors[i] for i in range(len(elements))]
        for i in range(len(elements)):
            if elements[i] == 26:
                metal = initOffsets[i]
                MonHerror = initErrors[i]
            elif elements[i] == 20 or elements[i] == 22:
                AonMerror += initErrors[i]**2
        AonMerror = AonMerror**0.5
        


        currentAlpha = np.average([  initOffsets[np.argmin(np.abs(elements - 20))]  ,  initOffsets[np.argmin(np.abs(elements - 22))]  ]) -(metal)
        if np.abs(currentAlpha-lastAlpha) > max(AonMerror,0.01):
            alpha = ((9*currentAlpha + lastAlpha)/10)
            for i in range(len(elements)):
                pass#initOffsets[i] = 0
        else:
            alpha = currentAlpha

        currentMetallicity = initOffsets[np.argmin(np.abs(elements - 26))]
        if np.abs(currentMetallicity-lastMetallicity) > max(MonHerror,0.01):
            metal = ((9*currentMetallicity + lastMetallicity)/10)
            for i in range(len(elements)):
                initOffsets[i] = 0
            alpha = lastAlpha
        else:
            metal = currentMetallicity
        



        print("[Fe/H]")
        print("this input: %2f" % lastMetallicity)
        print("this output: %2f" % currentMetallicity)
        print("next input : %2f" % metal)

        print("[A/Fe]")
        print("this input: %2f" % lastAlpha)
        print("this output: %2f" % currentAlpha)
        print("next input : %2f" % alpha)




    #saddgsdfgsdfw


    ##### GET STELLAR PARAMETERS FROM ISOCHRONES
    if "sun" in targ.lower():
        params = stellarParameters()
        if runNum == 0:
            params.metal = 0    
            params.alpha = 0 

        else:
            params.metal = metal    
            params.alpha = alpha   
        params.teff = 5777
        params.steff = 0
        params.logg = 4.44
        params.slogg = 0
        params.mass = 1
        params.smass = 0
        params.rad = 1
        params.srad = 0
        params.lum = 1
        params.slum = 0
   

    else:
        phot = queryPhotometry(targ)
        if metal in np.round(np.arange(-1.5,0.51,0.1),2):
            params =  MISTgetTandG(metal,"2", target = targ, SIMBADphot=phot)
            params.alpha = alpha 

        else:
            metLo = np.floor((metal)*10)/10
            metHi = np.ceil((metal)*10)/10
            paramsLo = MISTgetTandG(metLo,"2", target = targ, SIMBADphot=phot)
            paramsHi = MISTgetTandG(metHi,"2", target = targ, SIMBADphot=phot)
            
            params = stellarParameters()
            
            params.teff = ((paramsHi.teff - paramsLo.teff)/(metHi - metLo))*(metal - metLo) + paramsLo.teff
            params.steff = ((paramsHi.steff - paramsLo.steff)/(metHi - metLo))*(metal - metLo) + paramsLo.steff
            params.logg = ((paramsHi.logg - paramsLo.logg)/(metHi - metLo))*(metal - metLo) + paramsLo.logg
            params.slogg = ((paramsHi.slogg - paramsLo.slogg)/(metHi - metLo))*(metal - metLo) + paramsLo.slogg
            params.mass = ((paramsHi.mass - paramsLo.mass)/(metHi - metLo))*(metal - metLo) + paramsLo.mass
            params.smass = ((paramsHi.smass - paramsLo.smass)/(metHi - metLo))*(metal - metLo) + paramsLo.smass
            params.rad = ((paramsHi.rad - paramsLo.rad)/(metHi - metLo))*(metal - metLo) + paramsLo.rad
            params.srad = ((paramsHi.srad - paramsLo.srad)/(metHi - metLo))*(metal - metLo) + paramsLo.srad
            params.lum = ((paramsHi.lum - paramsLo.lum)/(metHi - metLo))*(metal - metLo) + paramsLo.lum
            params.slum = ((paramsHi.slum - paramsLo.slum)/(metHi - metLo))*(metal - metLo) + paramsLo.slum
            params.metal = metal    
            params.alpha = alpha    



    #### OUTPUT TO FILE
    print("Stellar Parameter fitting for target " + targ)
    print("%10s%10i%10i" % ("T_eff", params.teff, params.steff))
    print("%10s%10.2f%10.2f" % ("log(g)", params.logg, params.slogg))
    print("")
    print("%10s%10.2f%10.2f" % ("M/M_Sun", params.mass, params.smass))
    print("%10s%10.2f%10.2f" % ("R/R_Sun", params.rad, params.srad))
    print("%10s%10.2f%10.2f" % ("L/L_Sun", params.lum, params.slum))
    print("")
    print("%10s%10.2f%10.2f" % ("[M/H]",params.metal, MonHerror))
    print("%10s%10.2f%10.2f" % ("[alpha/M]", params.alpha, AonMerror))
    print("")

    print("This was iteration number %i" % runNum)

    with open(outputFolder + outputFile, "w") as f:
        print("%10s" % targ, file = f)
        print("*"*30, file = f)
        print("%10s %s" % ("WorkDir", outputFolder), file = f)
        print("%10s %s" % ("SpecFile", specFile), file = f)

        print("*"*30, file = f)
        print("%10s%10i" % ("FitAtom", elements[0]), file = f, end = "")
        for i in elements[1:]:   
            print(",%i" % i, file = f, end = "")

        print("\n" + "*"*30, file = f)
        print("%10s%10i%10i" % ("T_eff", params.teff, params.steff), file = f)
        print("%10s%10.2f%10.2f" % ("log(g)", params.logg, params.slogg), file = f)
        print("%10s%10.2f%10.2f"%("v_micro", 1.5,0.0), file = f)
        print("%10s%10.2f%10.2f"%("v_broad_0", 4.5,0.0), file = f)

        print("%10s%10.2f%10.2f" % ("[M/H]",params.metal, MonHerror), file = f)
        print("%10s%10.2f%10.2f" % ("[alpha/M]", params.alpha, AonMerror), file = f)
        print("", file = f)
        print("%10s%10.2f%10.2f" % ("M/M_Sun", params.mass, params.smass), file = f)
        print("%10s%10.2f%10.2f" % ("R/R_Sun", params.rad, params.srad), file = f)
        print("%10s%10.2f%10.2f" % ("L/L_Sun", params.lum, params.slum), file = f)
        print("", file = f)
        for i in range(len(bands)):
            print("%10s%10s%10i%10.3f" % ("LD_Coeff", bands[i], filterWaves[i], fitLDCurve(params.teff,params.logg,params.metal,bands[i])), file = f)
        print("", file = f)
        for i in range(len(elements)):
            print("%10s%10i%10.2f%10.2f" % ("SetAbund",elements[i],initOffsets[i]-(metal - lastMetallicity) - (alpha - lastAlpha)*(elements[i] % 2 == 0 and 7 < elements[i] < 23),initErrors[i]), file = f)
        os.rename("aaa.png", outputFolder+"IsochronalParameters.png")


if __name__ == "__main__":
    try:
        if len(sys.argv) == 4:
            setParams(sys.argv[1], sys.argv[2], sys.argv[3])
        elif len(sys.argv) == 7:
            setParams(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]))
        else:
            print("Usage: computeParamFile.py [Star Name] [Star Directory] [Output Param File] (Metallicity = 0) (alpha = 0) (runNum = 0)")
    except ValueError:
        exit(1)
