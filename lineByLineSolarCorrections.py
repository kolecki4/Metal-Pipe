import sys
import os

import pandas as pd
import numpy as np

#from mpplib.newIsoLib import *


def correctAbundances(workDir, inst):

    lblCorrections = pd.read_csv("SolarLBLData-trimmed.csv")
    with open(workDir + "params.txt") as f:
        paramFileLines = f.readlines()

    for line in paramFileLines:
        if "FitAtom" in line:
            elements = line.replace("FitAtom","").strip().split(",")
            break
    for element in elements:
        print(element)
        try:
            lineResults = pd.DataFrame(np.genfromtxt(workDir+element+"/abundanceSummary.txt", names = ["Wavelength","X_H","v_broad","chi2"], ndmin =1) ).sort_values("Wavelength").reset_index(drop = True)
            lblCorrectionsSubset = pd.merge(lblCorrections, lineResults,"outer", on="Wavelength").dropna(subset = [inst,"X_H"], how = "any").reset_index(drop = True)
            #print(lblCorrections)
            #print(lineResults)
            lblCorrectionsSubset["X_H"] = lblCorrectionsSubset["X_H"] - lblCorrectionsSubset[inst]  
            #print(lblCorrections)          
            lblCorrectionsSubset.to_csv(workDir+element+"/abundanceSummaryCorrected.csv", index = False, columns = ["Wavelength", "X_H", "v_broad", "chi2"])
            for i in range(len(paramFileLines)):
                if "%10s%10i" % ("SetAbund",int(element) ) in paramFileLines[i]:
                    print("ksfjldgjklgsdjlksf")
                    paramFileLines[i] = "%10s%10i%10.2f%10.2f\n" % ("SetAbund",int(element) ,np.nanmean(lblCorrectionsSubset["X_H"]),np.nanstd(lblCorrectionsSubset["X_H"])/len(lblCorrectionsSubset["X_H"])**0.5  )
                    break

        except FileNotFoundError:
            print("Skipping Element %s (no abundanceSummary.txt found)" % element)    
            for i in range(len(paramFileLines)):
                if "%10s%10i" % ("SetAbund",int(element) ) in paramFileLines[i]:
                    print("ksfjldgjklgsdjlksf")
                    paramFileLines[i] = "%10s%10i%10.2f%10.2f\n" % ("SetAbund",int(element) ,np.nan,np.nan)
                    break


    with open(workDir+"correctedParams.txt", "w") as f:
        f.writelines(paramFileLines)
if __name__ == "__main__":
    #try:
    if len(sys.argv) == 3:
        correctAbundances(sys.argv[1], sys.argv[2].upper())
    else:
        print("Usage: correctAbundances.py [workDir] [instrument]")
    #except Exception as e:
    #    print(e)        
    #    sys.exit(1)
