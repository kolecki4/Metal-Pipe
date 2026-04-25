import numpy as np
import numpy.lib.recfunctions as rfn
import sys
import os
from mpplib.newIsoLib import *
from scipy.interpolate import RegularGridInterpolator



#outputFolder = "/home/jared/Documents/Spectra/KPF/CAP4/34411/"
outputFolder = sys.argv[1]
outputFile = "params.txt"
element = sys.argv[2]
corrFile = element + ".txt"

elements = np.array([26,20,22,12,14,6,7,8,16,11,13,19])
elements, abundances, errors = readAbundances(outputFolder, outputFile, elements)




with open(outputFolder + outputFile) as f:
    paramFileLines = f.readlines()


teff = float(paramFileLines[7].split()[1])
logg = float(paramFileLines[8].split()[1])
xi = float(paramFileLines[9].split()[1])
mh = float(paramFileLines[11].split()[1])



nlteDir = "NLTE_Corrections/"
nlteCorrFiles = ["O.txt", "S.txt", "Ti.txt", "Ca.txt"]

if corrFile in ["C.txt", "O.txt", "S.txt"]:
    columns = ['Teff', 'logg', 'Fe', 'xi', 'leps', 'wav', 'adj']
    table = np.genfromtxt(nlteDir + corrFile, names = columns)
    table = np.sort(table, order='wav')
elif corrFile in ["K.txt"]:
    columns = ['Teff', 'logg', 'Fe', 'xi', 'leps', 'wav','x','y', 'adj']
    table = np.genfromtxt(nlteDir + corrFile, names = columns)
    table = np.sort(table, order='wav')

elif corrFile in ["Ti.txt"]:
    columns = ['info','wav1','wav2','wav3','wav4','wav5','wav6','wav7','wav8','wav9']
    table = np.genfromtxt(nlteDir + corrFile, names = columns, dtype = None)
    col1 = np.array([float(i[1:5]) for i in table["info"]])
    col2 = np.array([float(i[6:9]) for i in table["info"]])/100
    col3 = np.array([float(i[11:13]) for i in table["info"]])/10
    if table["info"][10] == "m": col3 = -col3
    col4 = np.ones_like(table["wav1"])*1.5
    col5 = np.ones_like(table["wav1"])*4.95 
    col6 = np.ones_like(table["wav1"])*5500
    col7 = np.array([     np.average([ table[ columns[i] ][j] for i in range(1,len(columns)) if table[ columns[i] ][j] < 10]) for j in range(len(col1)) ])
    
    table = rfn.append_fields(table,["Teff", "logg", "Fe", "xi", "leps", "wav", "adj"], [col1,col2,col3,col4,col5,col6,col7])

elif corrFile in ["Ca.txt"]:
    columns = ['info','wav1','wav2','wav3','wav4','wav5','wav6']
    table = np.genfromtxt(nlteDir + corrFile, names = columns, dtype = None)
    col1 = np.array([float(i[1:5]) for i in table["info"]])
    col2 = np.array([float(i[6:9]) for i in table["info"]])/100
    col3 = np.array([float(i[11:13]) for i in table["info"]])/10
    if table["info"][10] == "m": col3 = -col3
    col4 = np.ones_like(table["wav1"])*1.5
    col5 = np.ones_like(table["wav1"])*4.95 
    col6 = np.ones_like(table["wav1"])*5500
    col7 = np.array([     np.average([ table[ columns[i] ][j] for i in range(1,len(columns)) if table[ columns[i] ][j] < 10]) for j in range(len(col1)) ])
    
    table = rfn.append_fields(table,["Teff", "logg", "Fe", "xi", "leps", "wav", "adj"], [col1,col2,col3,col4,col5,col6,col7])





if corrFile == "O.txt":
    if 8 in elements:
        ab = abundances[np.argmax(elements == 8)]
    else:
        raise ArithmeticError
elif corrFile == "S.txt":
    if 16 in elements:
        ab = abundances[np.argmax(elements == 16)]
    else:
        raise ArithmeticError
elif corrFile == "K.txt":
    if 19 in elements:
        ab = abundances[np.argmax(elements == 19)]
    else:
        raise ArithmeticError
elif corrFile == "Ca.txt":
    if 20 in elements:
        ab = abundances[np.argmax(elements == 19)]
    else:
        raise ArithmeticError
elif corrFile == "Ti.txt":
    if 22 in elements:
        ab = abundances[np.argmax(elements == 19)]
    else:
        raise ArithmeticError



dfghd =0
for i in range(len(table['Teff'])):
    if 4900 < table['Teff'][i] < 5100:
        table['Teff'][i] = 5000
    if 5400 < table['Teff'][i] < 5600:
        table['Teff'][i] = 5500
    if 5900 < table['Teff'][i] < 6100:
        table['Teff'][i] = 6000
    if 6400 < table['Teff'][i] < 6600:
        table['Teff'][i] = 6500
    

teffs = np.unique(table['Teff'])
loggs = np.unique(table['logg'])
feabs = np.unique(table['Fe'])
vmics = np.unique(table['xi'])
lgeps = np.unique(table['leps'])
waves = np.unique(table['wav'])

teffs = teffs[(np.abs(teffs - teff) < 1500)]
loggs = loggs[(np.abs(loggs - logg) <= 2)]
feabs = feabs[(np.abs(feabs - mh) <= 2)]
vmics = vmics[(np.abs(vmics - xi) <= 0.5)]
lgeps = lgeps[(np.abs(lgeps - ab) < 1.5)]


def gridcorrection(teff,logg,fe,xi,leps,wav):
    T = teffs[np.argmin(np.abs(teffs-teff))]
    G = loggs[np.argmin(np.abs(loggs-logg))]
    M = feabs[np.argmin(np.abs(feabs-fe))]
    X = vmics[np.argmin(np.abs(vmics-xi))]
    E = lgeps[np.argmin(np.abs(lgeps-leps))]
    W = waves[np.argmin(np.abs(waves-wav))]
    #print(T,G,M,X,E,W)
    
    conds = [table['Teff'] == T,
                  table['logg'] == G,
                  table['Fe'] == M,
                  table['xi'] == X,
                  table['leps'] == E,
                  table['wav'] == W]

    gridcorr =  table['adj'][np.where(conds[0] & conds[1] & conds[2] & conds[3] & conds[4] & conds[5])]
    if len(gridcorr) == 0: gridcorr = np.append(gridcorr,np.nan)
    #print(gridcorr[0])
    return gridcorr[0]




NLTEgrid = np.zeros((len(teffs),len(loggs),len(feabs),len(vmics),len(lgeps),len(waves)))
for i in range(len(teffs)):
    for j in range(len(loggs)):
        for k in range(len(feabs)):
            for l in range(len(vmics)):
                for m in range(len(lgeps)):
                    for n in range(len(waves)):
                        NLTEgrid[i,j,k,l,m,n] = gridcorrection(teffs[i],loggs[j],feabs[k],vmics[l],lgeps[m],waves[n])
for i in range(len(teffs)):
    for j in range(len(loggs)):
        for k in range(len(feabs)):
            for l in range(len(vmics)):
                for n in range(len(waves)):
                    try: 
                        NLTEgrid[i,j,k,l,:,n] = np.interp(lgeps, lgeps[~np.isnan(NLTEgrid[i,j,k,l,:,n])], NLTEgrid[i,j,k,l,:,n][~np.isnan(NLTEgrid[i,j,k,l,:,n])])
                    except:
                        pass
                    
                    
def interpcorrection(T,G,X,M,O):                    
    try:
        Tm1 = max(teffs[np.where(teffs <= T)])
    except:
        Tm1 = min(teffs)
        print('Warning: Teff is', Tm1 - T, 'below the interpolation range')
    try:
        Tp1 = min(teffs[np.where(teffs >= T)])
    except:
        Tp1=max(teffs)
        print('Warning: Teff is', T - Tp1, 'above the interpolation range')
        
    try:
        Gm1 = max(loggs[np.where(loggs <= G)])
    except:
        Gm1=min(loggs)
        print('Warning: logg is', Gm1 - G, 'below the interpolation range')
    try:
        Gp1 = min(loggs[np.where(loggs >= G)])
    except:
        Gp1 = max(loggs)
        print('Warning: logg is', G - Gp1, 'above the interpolation range')
    try:
        Xm1 = max(vmics[np.where(vmics <= X)])
    except:
        Xm1 = min(vmics)
        print('Warning: xi is', Xm1 - X, 'below the interpolation range')
    try:
        Xp1 = min(vmics[np.where(vmics >= X)])
    except:
        Xp1 = max(vmics)
        print('Warning: xi is', X - Xp1, 'above the interpolation range')
    
    try:
        Mm1 = max(feabs[np.where(feabs <= M)])
    except:
        Mm1 = min(feabs)
        print('Warning: [Fe/H] is', Mm1 - M, 'below the interpolation range')
    try:
        Mp1 = min(feabs[np.where(feabs >= M)])
    except:
        Mp1 = max(feabs)
        print('Warning: [Fe/H] is', M - Mp1, 'above the interpolation range')
    
    try:
        Om1 = max(lgeps[np.where(lgeps <= O)])
    except:
        Om1 = min(lgeps)
        print('Warning: log(eps) is', Om1 - O, 'below the interpolation range')
    
    try:
        Op1 = min(lgeps[np.where(lgeps >= O)])
    except:
        Op1 = max(lgeps)
        print('Warning: log(eps) is', O - Op1, 'above the interpolation range')
    
    
    
    
    
    if Tp1 != Tm1:
        mapT = (T-Tm1)/(Tp1-Tm1)
    else:
        mapT = 0.5
    if Gp1 != Gm1:
        mapG = (G-Gm1)/(Gp1-Gm1)
    else:
        mapG = 0.5
    if Xp1 != Xm1:
        mapX = (X-Xm1)/(Xp1-Xm1)
    else:
        mapX = 0.5  
    if Mp1 != Mm1:
        mapM = (M-Mm1)/(Mp1-Mm1)
    else:
        mapM = 0.5  
    if Op1 != Om1:
        mapO = (O-Om1)/(Op1-Om1)
    else:
        mapO = 0.5  
    
    grid = np.zeros((2,2,2,2,2,len(waves)))
    


    corrections = []
    for i in range(len(waves)):
        grid[0,0,0,0,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,0,0,0,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,0,0,1,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,0,0,1,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,0,1,0,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,0,1,0,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,0,1,1,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,0,1,1,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,1,0,0,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,1,0,0,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,1,0,1,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,1,0,1,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,1,1,0,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,1,1,0,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[0,1,1,1,0,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[0,1,1,1,1,i] = NLTEgrid[np.argmin(np.abs(Tm1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,0,0,0,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,0,0,0,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,0,0,1,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,0,0,1,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,0,1,0,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,0,1,0,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,0,1,1,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,0,1,1,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gm1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,1,0,0,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,1,0,0,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,1,0,1,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,1,0,1,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mm1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,1,1,0,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,1,1,0,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xm1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        grid[1,1,1,1,0,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Om1-lgeps)),i]
        grid[1,1,1,1,1,i] = NLTEgrid[np.argmin(np.abs(Tp1-teffs)),np.argmin(np.abs(Gp1-loggs)),np.argmin(np.abs(Mp1-feabs)),np.argmin(np.abs(Xp1-vmics)),np.argmin(np.abs(Op1-lgeps)),i]
        myinterp = RegularGridInterpolator([[0,1],[0,1],[0,1],[0,1],[0,1]], grid[:,:,:,:,:,i])
        #if waves[i] == 777:
        #print(waves[i]*10,myinterp([mapT,mapG,mapM,mapX,mapO]))
        corrections.append(myinterp([mapT,mapG,mapM,mapX,mapO]))
    return waves, np.array(corrections)

print("Non-LTE corrections for " + element)
w, c = interpcorrection(teff,logg,xi,mh,ab)
if element == "O":


    with open(outputFolder + "nlte.txt", "w") as f:
        print("For the O I triplet at 7770A: delta = %.2f" % c[w == 777.0 ][0][0], file = f)
        print("    7772A: delta = %.3f" % c[w == 777.2 ][0][0], file = f)
        print("    7774A: delta = %.3f" % c[w == 777.4 ][0][0], file = f)
        print("    7775A: delta = %.3f" % c[w == 777.5 ][0][0], file = f)

    print("For the O I triplet at 7770A: delta = %.2f" % c[w == 777.0 ][0][0])
    print("    7772A: delta = %.3f" % c[w == 777.2 ][0][0])
    print("    7774A: delta = %.3f" % c[w == 777.4 ][0][0])
    print("    7775A: delta = %.3f" % c[w == 777.5 ][0][0])

if element == "K":
    with open(outputFolder + "nlte.txt", "a") as f:
        print("For the K I line at 7698A: delta = %.2f" % c[w == 7698.0 ][0][0], file = f)

    print("For the K I line at 7698A: delta = %.2f" % c[w == 7698.0 ][0][0])


if element == "Ca":
    with open(outputFolder + "nlte.txt", "a") as f:
        print("For Ca in the Optical: delta = %.2f" % c[w == 5500.0 ][0][0], file = f)

    print("For Ca in the Optical: delta = %.2f" % c[w == 5500.0 ][0][0])
if element == "Ti":
    with open(outputFolder + "nlte.txt", "a") as f:
        print("For Ti in the Optical: delta = %.2f" % c[w == 5500.0 ][0][0], file = f)

    print("For Ti in the Optical: delta = %.2f" % c[w == 5500.0 ][0][0])

