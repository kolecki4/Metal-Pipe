import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr, fit, t
import os
import sys

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
       -5.00,-0.54,-5.00,-5.00,-5.00
]

size = 24
plt.rcParams.update({'figure.dpi': 200})
plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'font.family': "FreeSans"})

plt.rcParams.update({'axes.linewidth': 4})
plt.rcParams.update({'xtick.major.width': 4})
plt.rcParams.update({'xtick.minor.width': 4})
plt.rcParams.update({'ytick.major.width': 4})
plt.rcParams.update({'ytick.minor.width': 4})

plt.rcParams.update({'xtick.major.size': 8})
plt.rcParams.update({'xtick.minor.size': 8})
plt.rcParams.update({'ytick.major.size': 8})
plt.rcParams.update({'ytick.minor.size': 8})

plt.rcParams.update({'ytick.labelcolor':  "#222"})
plt.rcParams.update({'xtick.labelcolor':  "#222"})
plt.rcParams.update({'ytick.color':  "#444"})
plt.rcParams.update({'xtick.color':  "#444"})

plt.rcParams.update({"text.color": "#222"})
plt.rcParams.update({"axes.edgecolor": "#444"})
plt.rcParams.update({"lines.markersize":10})

labelsize = 24


maxLines = 9
nRows = 3
plotWidth = 4

dataFolder = sys.argv[1]
folders = np.sort(os.listdir(dataFolder))
folders = folders[np.array([len(i) for i in folders]) < 3]
for folder in folders:
    files = np.sort(os.listdir(dataFolder  + folder))
    try:
        summary = np.genfromtxt(dataFolder  + folder + "/abundanceSummary.txt", names = ["wav", "XH", "v_broad", "chi2"], ndmin = 1)


        files = [i for i in files if "000" in i]

        newFiles = []
        for i in range(len(files)):
            if float(files[i]) in summary["wav"]:
                newFiles.append(files[i])
        files = np.array(newFiles)




        solar = [0,
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
            ]



        # Remove near-duplicate wavelengths
        files = files[np.unique(summary["wav"], return_index = True)[1]]
        summary = summary[np.unique(summary["wav"], return_index = True)[1]]
        size = len(summary["wav"])

        #print("files %i summary %i" % (len(files), len(summary)))
        #Remove obviously bad abundances
        if np.max(summary["XH"]) - np.min(summary["XH"]) > 1.5:
            files = files[(summary["XH"] > 0.3 + np.min(summary["XH"])) & (summary["XH"] < np.max(summary["XH"])-0.3) ]
            summary = summary[(summary["XH"] > 0.3 + np.min(summary["XH"])) & (summary["XH"] < np.max(summary["XH"])-0.3) ]

#        for i in range(len(files)):
#            print(files[i], np.argsort(summary["chi2"])[i])

        if len(files) > 9:
            Best = files[np.argsort(summary["chi2"])][0:9]
            BestAbs = summary["XH"][np.argsort(summary["chi2"])][0:9]

        else:
            Best = files
            BestAbs = summary["XH"]
        
        #print(Best)



        fig, ax = plt.subplots(nRows, 3, figsize = (plotWidth*nRows*2,plotWidth*nRows*1.5), layout = "constrained")
        fig.suptitle(folder, fontsize = 40)
        
        for i in range(len(Best)):
                 
            wav, obs, sigma_obs, syn, weight = np.genfromtxt(dataFolder + folder + "/" + Best[i], delimiter=",", skip_header = 1, unpack = True)
            sigma_syn = (-0.05*np.log(syn))**1
            err = np.sqrt(sigma_syn**2 + sigma_obs**2)
            swav = wav

            #plt.plot(sWavMatched,sFlxMatched, ls = '', marker = '.', color = "#dd0000")    
            unweightedChi2 = ((obs-syn)**2) /(err**2)/(len(wav) - 3)
            weightedChi2 = weight * ((obs-syn)**2) /(err**2) /(len(wav) - 3)

            ax[i%nRows][i//nRows].errorbar(wav, obs, yerr = np.abs(sigma_obs), ls = '-', marker = '', color = "#222222", label = "Observed Line", lw = 1)
            ax[i%nRows][i//nRows].fill_between(wav,obs - sigma_obs,obs + sigma_obs, ls = '-', color = "#222222", alpha = 0.4)

            ax[i%nRows][i//nRows].errorbar(wav,syn, yerr = np.abs(sigma_syn), ls = '-', marker = '.', color = "#dd0000", label = "Best-fit Synthesis", markersize = 2, lw = 1)
            ax[i%nRows][i//nRows].set_title( "[X/H] = %.2f, weighted chi2 = %.2f"% (BestAbs[i], np.sum(weightedChi2)))
            ax[i%nRows][i//nRows].fill_between(wav,syn - sigma_syn,syn + sigma_syn, ls = '-', color = "#dd0000", alpha = 0.2)


            ax[i%nRows][i//nRows].legend(fontsize = 18, framealpha = 0, loc = "lower left")

            if i%nRows == nRows - 1:
                ax[i%nRows][i//nRows].set_xlabel("Wavelength (\u212b)")

            if i//nRows == 0:
                ax[i%nRows][i//nRows].set_ylabel("Normalized Flux", fontsize = 32)
                #ax[1][i].set_ylabel(r"$\chi^2_{\nu}$" + " per point", usetex = True, fontsize = 32, fontfamily = "FreeSans")
                #ax[2][i].set_ylabel("Weight", fontsize = 32)
            
         
            #ax[i%nRows][i//nRows].set_xticks(np.arange(np.floor(np.min(wav)), np.ceil(np.max(wav)),0.2), labels = ["%.3g" % i for i in np.arange(np.floor(np.min(wav)), np.ceil(np.max(wav)),0.2)])
            
            
            
            
            ax[i%nRows][i//nRows].set_ylim(np.min(np.append(obs,syn))-0.1,1.1)
            ax[i%nRows][i//nRows].set_xlim(np.min(wav), np.max(wav))
        print(dataFolder + folder + "/diagnostic_plot.png")    
        plt.savefig(dataFolder + folder + "/diagnostic_plot.png")
        plt.close("all")




        plt.hist(summary["XH"], bins = np.arange(4.5,8.9, 0.05))
        plt.xlabel("[X/H]")
        print("---------")
        print(summary["XH"])
        print("---------")
        plt.xlim(np.percentile(summary["XH"],1)-0.1,np.percentile(summary["XH"],99)+ 0.1)
        plt.title("Median [X/H] = %.3f +/- %.3f" % (np.median(summary["XH"]),np.std(summary["XH"])/np.sqrt(len(summary))))

        plt.savefig(dataFolder + folder + "/Abundance Histogram.png")
        plt.close("all")

        plt.figure()
        plt.hist(summary["v_broad"], bins = np.arange(0,12.1, 0.1))
        plt.xlabel("v_broad")
        plt.title("Mean v_broad = %.3f +/- %.3f" % (np.average(summary["v_broad"]),np.std(summary["v_broad"])/np.sqrt(len(summary))))

        plt.savefig(dataFolder + folder + "/VBroad Histogram.png")
        plt.close("all")

        plt.figure()
        plt.hist(summary["chi2"], bins = np.arange(0.0,5.1, 0.1))
        plt.xlabel(r"Weighted $\chi^2_{\nu}$")
        plt.title("Median Chi^2 = %.3f +/- %.3f" % (np.median(summary["chi2"]),np.std(summary["chi2"])/np.sqrt(len(summary))))
            
        plt.savefig(dataFolder + folder + "/Chi2 Histogram.png")
        plt.close("all")
    except IOError:
        pass
