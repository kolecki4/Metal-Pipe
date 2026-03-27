import sys
import os

import numpy as np

from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

from scipy.ndimage import gaussian_filter
from scipy.stats import linregress
from scipy.integrate import simpson
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

from astropy.io import fits


#plt.style.use("/home/jared/pltstyles/jared.mplstyle")

class stellarParameters:    
    def __init__(self,teff = None,steff = None,logg = None,slogg = None,mass = None,smass = None,rad = None,srad = None,lum = None,slum = None,metal = None, alpha = None):
        self.teff = teff
        self.steff=steff
        self.logg = logg
        self.slogg = slogg
        self.mass = mass
        self.smass = smass
        self.rad = rad
        self.srad = srad
        self.lum = lum
        self.slum = slum
        self.metal = metal
        self.alpha = alpha

def queryPhotometry(star):
    """
    Uses astroquery to return an array containing 
    absolute magnitude photometric data for a target

    Parameters
    ----------
    star : string
        The SIMBAD name of a star.

    Returns
    -------
    starmags : ndarray
        An array which compiles UVBRI, 2MASS, Gaia, and WISE 
        photometric data for the star, with named columns
        according to the band name.

    """
    customSimbad=Simbad()
    customSimbad.add_votable_fields('parallax', 'plx_err')
    #customSimbad.add_votable_fields('flux_error(U)','flux_error(V)','flux_error(B)','flux_error(R)','flux_error(I)','flux_error(J)','flux_error(H)','flux_error(K)')
    customSimbad.add_votable_fields("flux")    
    columns = [('U',float),('B',float),('Bp',float),('V',float),
           ('R',float),('G',float),('Rp',float),('I',float),('J',float),('H',float),
           ('K',float),
           ('W1',float),('W2',float),('W3',float),('W4',float)]
    
    starmags = np.zeros(1, dtype = columns)
    magerrs = np.zeros(1, dtype = columns)
    
    #Query Object
    sim=customSimbad.query_object(star)
    ids=customSimbad.query_objectids(star, criteria="ident.id LIKE 'Gaia%'")
    #print(ids[0][0])
    try:
        gaia2 = None            
        
        for i in ids:
            if "Gaia DR2" in i[0]:
                gaia2 = i[0]
        gaia2 = gaia2.replace("Gaia DR2", "")
    except IndexError:
        pass

    try:
        gaia3 = ids[["Gaia DR3" in i for i in ids]][0]
        gaia3 = gaia3.replace("Gaia DR3", "") 
    except IndexError:
        gaia3 = None

    #print(sim["ID_GAIA_DR2"])
    if gaia2 is not None:
        gaiaid=gaia2
        gaia2Vizier = Vizier(columns = ["Source", "Plx", "e_Plx", "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag"])
        gaia2Vizier.ROW_LIMIT = -1
        gai=gaia2Vizier.query_object(star, catalog = 'I/345/gaia2', radius = 120*u.arcsec, coordinate_frame = "J2016")[0]
        gai = gai[np.argwhere(gai["Source"] == int(gaiaid))[0][0]]
        if type(gai['Plx']) == np.float64:
            plx=gai['Plx']
            eplx = gai['e_Plx']

            
        try:
            fuck
            cbj=Vizier.query_object(star, 'I/347/gaia2dis')
            d = cbj['I/347/gaia2dis']['rest'][np.where(cbj[0]['Source'] == int(gaiaid))[0][0]]
            #print('Bailer-Jones+2018 distance found:',d)
        except:
            #print('Using 1/parallax for distance value')
            d = 1/(plx/1000)
            ed = d*eplx/plx
        starmags['G']=gai['Gmag'] - 5*np.log10(d) + 5
        starmags['Bp']=gai['BPmag'] - 5*np.log10(d) + 5
        starmags['Rp']=gai['RPmag'] - 5*np.log10(d) + 5
        magerrs['G']=(gai['e_Gmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Bp']=(gai['e_BPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Rp']=(gai['e_RPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
    


    elif gaia3 is not None:
        gaiaid=sim['ID_GAIA_DR2'][0].replace('Gaia DR3 ', '')
        #print(gaiaid)
        gaia3Vizier = Vizier(columns = ["Source", "Plx", "e_Plx", "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag"])
        gaia3Vizier.ROW_LIMIT = -1
        gai=gaia3Vizier.query_object(star, catalog = 'I/355/gaiadr3', radius = 120*u.arcsec, coordinate_frame = "J2016")[0]
        gai = gai[np.argwhere(gai["Source"] == int(gaiaid))[0][0]]
        if type(gai['Plx']) == np.float64:
            plx=gai['Plx']
            eplx = gai['e_Plx']


        try:
            fuck
            cbj=Vizier.query_object(star, 'I/347/gaia2dis')
            d = cbj['I/347/gaia2dis']['rest'][np.where(cbj[0]['Source'] == int(gaiaid))[0][0]]
            #print('Bailer-Jones+2018 distance found:',d)
        except:
            #print('Using 1/parallax for distance value')
            d = 1/(plx/1000)
            ed = d*eplx/plx

        starmags['G']=gai['Gmag'] - 5*np.log10(d) + 5
        starmags['Bp']=gai['BPmag'] - 5*np.log10(d) + 5
        starmags['Rp']=gai['RPmag'] - 5*np.log10(d) + 5
        magerrs['G']=(0.002755**2 + gai['e_Gmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Bp']=(0.00279**2 + gai['e_BPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
        magerrs['Rp']=(0.003779**2 + gai['e_RPmag']**2 + (5*ed/(d*np.log(10)))**2)**0.5
            
    else:
        #print("No GAIA parallax found. Using parallax from SIMBAD sourced from " + sim['PLX_BIBCODE'][0])
        #print('Using 1/parallax for distance value')
        d = 1/(sim['plx_value'][0]/1000)
        ed = d*sim['plx_err'][0]/sim['plx_value'][0]
        
    
    badBibcodes = ["2012yCat.1322....0Z", "2005yCat.2263....0D", "2009yCat.1315....0Z"]
 
    for band in ["U","B","V","R","I","J","H","K"]:
        
        try: 
            if sim['flux.bibcode'][sim["flux.filter"] == band ][0] not in badBibcodes:
                starmags[band]=sim['flux'][sim["flux.filter"] == band ][0] - 5*np.log10(d) + 5
                magerrs[band]=(sim['flux_err'][sim["flux.filter"] == band ][0]**2 + (5*ed/(d*np.log(10)))**2)**0.5
        except IndexError: 
            pass    

    
    for band in ["U","B","V","R","I","J","H","K"]:
        if magerrs[band] == 0 and starmags[band] != 0:
            magerrs[band] = 0.1

    return starmags, magerrs

def MISTreadiso(age, metall):
    """
    Reads a Dartmouth theoretical isochrone into a numpy array
    
    Parameters
    ----------
    age : float
        The age in Gyr of the isochrone to read. (Rounded to the nearest available value)
    metall : float
        The metallicity ([Fe/H] relative to solar) of the isochrone to read.

    Returns
    -------
    combinediso : ndarray
        The requested isochrone with all synthetic photometry 
        and with all columns labelled.

    """
    isopath = 'isochrones/'
    
    if metall < -1.5:
        metall = -1.5
    if metall > 0.5:
        metall = 0.5
    
    if metall >= 0:                                               
        s1='p'                                                   
    else:                                                        
        s1='m'
    s2 = str(round(10*metall)/10).replace('.','').replace('-','')
    
    mstr = s1+s2
    
    outColumns = [('EEP',float),('M/Mo',float),('LogTeff',float),('LogG',float),
               ('LogL/Lo',float),('U',float),('B',float),('V',float),
               ('R',float),('I',float),('J',float),('H',float),
               ('K',float),('G',float),('Bp',float),('Rp',float),
               ('W1',float),('W2',float),('W3',float),('W4',float)]
    
    availAges = np.arange(5,10.11,0.05)
    
    ### GAIA/2MASS/UVBRI
    with open(isopath + mstr + 'G2U.iso.cmd', 'r') as f:
        logAge = np.log10(age)+9

        isoStartInds = []
        isoLen = -1
        nCols = -1
        
        targetIso = np.argmin(np.abs(availAges - logAge))
        
        j = 0
        while True:
            line = f.readline()
            j+=1
            if "number of isochrones =" in line:
                niso = int(line[-7:-1])
                break
            
        while True:
            line = f.readline()
            ### This denotes the start of an isochrone. Make note of where
            if 'number of EEPs, cols' in line:
                isoStartInds.append(j)
                isoLen = int(line[-9:-3])
                nCols = int(line[-3:-1])
            
            ### This implies we have reached the correct isochrone. Extract it
            if len(isoStartInds) == targetIso+1:
                
                ### Initialize the output array
                combinediso = np.zeros(isoLen, dtype = outColumns)
                ### Skip over unneeded text
                f.readline()
                
                ### Note the names of the columns in the MIST file
                columns = []
                for col in f.readline().split()[1:]:
                    columns.append(col)
                
                ### Line by line, put the data into the combined iso array
                for i in range(isoLen): 
                    line = f.readline().split()
                    combinediso[i] = (line[0],line[3],line[4],line[5],line[6],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[22],line[23],line[24],0,0,0,0)
                break
            j+=1
            
    ### WISE        
    with open(isopath + mstr + 'WISE.iso.cmd', 'r') as f:
    
        isoStartInds = []
        isoLen = -1
        nCols = -1
    
        j = 0
        while True:
            line = f.readline()
            j+=1
            if "number of isochrones =" in line:
                niso = int(line[-7:-1])
                break
            
        while True:
            line = f.readline()
            ### This denotes the start of an isochrone. Make note of where
            if 'number of EEPs, cols' in line:
                isoStartInds.append(j)
                isoLen = int(line[-9:-3])
                nCols = int(line[-3:-1])
            
            ### This implies we have reached the correct isochrone. Extract it
            if len(isoStartInds) == targetIso+1:
                ### Skip over unneeded text
                f.readline()
                
                ### Note the names of the columns in the MIST file
                columns = []
                for col in f.readline().split()[1:]:
                    columns.append(col)
                
                ### Line by line, put the data into the combined iso array
                for i in range(isoLen): 
                    line = f.readline().split()
                    combinediso[i][16],combinediso[i][17],combinediso[i][18],combinediso[i][19] = (line[9],line[10],line[11],line[12])
                break
            j+=1
        
    return combinediso

def gibbs_sampling(image, w_start, samples):

    image_pdf = image 
    height, width = image_pdf.shape
    result = []
    w_current = w_start

    for _ in range(samples):
        # sample height
        h_given_w = image_pdf[:, w_current] / image_pdf[:, w_current].sum()
        h_given_w = np.where(np.isnan(h_given_w), 0, h_given_w)
        
        if np.all(h_given_w == 0):
            h_given_w += 1
            h_given_w /= np.sum(h_given_w)

        h_current = np.random.choice(np.array(range(height)), size=1, p=h_given_w)[0]

        # sample width
        w_given_h = image_pdf[h_current, :] / image_pdf[h_current, :].sum()
        w_current = np.random.choice(np.array(range(width)), size=1, p=w_given_h)[0]

        result.append((h_current, w_current))

    return result

def rejection_sampling(image, approx_samples):

    image_pdf = image / image.sum()
    pdf_max = image_pdf.max()
    height, width = image_pdf.shape
    p_success = 1 / (height * width * pdf_max)
    actual_samples = min(int(approx_samples / p_success), int(1e8))

    samples_height = np.random.randint(0, high=height, size=actual_samples)
    samples_width = np.random.randint(0, high=width, size=actual_samples)
    samples_uniform = np.random.uniform(0, 1, size=actual_samples)

    result = [(h, w) for (h, w, u) in zip(samples_height, samples_width, samples_uniform) if
                    (image_pdf[h, w] >= pdf_max * u)]

    return result


def MISTgetTandG(metal, ageRange, target = None, SIMBADphot = None):
    """
    Derives an effective temperature and surface gravity
    for a target by comparing photometry to a MIST 
    theoretical isochrone of a given age and metallicity.
    Also required is an age range which is used to prevent
    inaccuracies
    
    Parameters
    ----------
    metall : float
        The metallicity ([Fe/H] relative to solar) of the isochrone to read.
    ageRange : int
        Flag for MIST/Dartmouth isochrone age range
            0: MIST standard age grid log(Age [yr]) = [5, 7.4]
            1: MIST standard age grid log(Age [yr]) = [7, 9.4]
            2: MIST standard age grid log(Age [yr]) = [9, 10.1]
    target : string
        The SIMBAD name of a star.
    SIMBADphot : ndarray
        The returned value of the querySIMBADphot function.
        A timesaver if you're doing multiple calls of this function
        on the same star.
    
    Returns
    -------
    T : float
        Effective Temperature.
    G : float
        Log(g).
    T_err : float
        The uncertainty on T.
    G_err : float
        The uncertainty on G.
    
    """
    
    if SIMBADphot == None and target != None: 
        SIMBADphot=queryPhotometry(target)
        starmags=SIMBADphot[0]
        magerrs=SIMBADphot[1]
    elif SIMBADphot != None:
        starmags=SIMBADphot[0]
        magerrs=SIMBADphot[1]
        
    
#    if ageRange == '0':
#        availages = 10**(np.arange(5,7.4,0.1)-9)
#        linages = np.arange(.0002,.025,.0002)
#    elif ageRange == '1':   
#        availages = 10**(np.arange(7,9.4,0.1)-9)
#        linages = np.arange(.01,2.5,.02)
#    elif ageRange == '2':
#        availages = 10**(np.arange(9,10.1,0.1)-9)
#        linages = np.arange(.1,13,.1)
#    else:
#        ageRange = ageRange.replace('[','').replace(']','').split(',')
#        ageRange = [float(i) for i in ageRange]
#        availages = 10**(np.arange(5,10.11,0.05)-9)
#        availages = availages[(availages > ageRange[0]) & (availages < ageRange[1])]
#        linages = np.linspace(ageRange[0],ageRange[1],125)
    # availages = np.append(np.arange(1,5,0.25),np.arange(5,13.1,0.5))
    availages = 10**(np.arange(8,10.1,0.05)-9)
    # availages = 10**(np.arange(7,9.4,0.1)-9)
    # availages = 10**(np.arange(5,7.4,0.1)-9)
    
    #availages = 10**(np.arange(6,9,0.05)-9)
    
    # Let the program know what data we have that's good to use
    bands = np.array(['U', 'B', 'Bp', 'V', 'G', 'R', 'Rp', 'I', 'J', 'H', 'K', 'W1', 'W2', 'W3', 'W4'])
    waves = np.array([365, 445,  532, 551, 640, 658,  797, 806,1250,1650,2150, 3400,4600, 12000,22000])
    
    goodPoints = np.zeros_like(waves)
    for i in bands:
        if starmags[i] != 0:
            print("%2s: %.2f %.2f" %(i, starmags[i], magerrs[i]))
            goodPoints[np.argwhere(bands == i)[0][0]] = 1

    bands = bands[goodPoints == 1]
    waves = waves[goodPoints == 1] 
    
    smags = np.array([])
    emags = np.array([])
    imags = np.array([])
    
    for i in bands:
        smags = np.append(smags, starmags[i])
        emags = np.append(emags, magerrs[i])
    
    teffs = np.zeros([len(availages),1500])
    loggs = np.zeros([len(availages),1500])
    mass = np.zeros([len(availages),1500])
    lum = np.zeros([len(availages),1500])
    rad = np.zeros([len(availages),1500])
    
    residuals = np.zeros(shape = [len(availages),1500])
    chi2 = np.zeros(shape = [len(availages),1500])+9e9
    likelihood = np.zeros(shape = [len(availages),1500])
    logl = np.zeros(shape = [len(availages),1500])
    imagslist = np.zeros(shape = [len(availages),1500], dtype = list)
    for k in range(len(availages)):
        isoc = MISTreadiso(availages[k],metal)
        
        for i in range(len(isoc['EEP'])):
            imags = np.array([])
            
            for j in range(len(bands)):
                imags = np.append(imags,isoc[bands[j]][i])
        
            resids = smags - imags
            chi2[k,i] = np.sum((resids**2)/(emags)**2)/(len(emags) - 2)

            logl[k,i] = -0.5 * np.sum( (resids/emags)**2  -np.log(2*np.pi*emags )  )
            
            residuals[k,i] = sum(np.abs(resids))
            teffs[k,i] = 10**isoc['LogTeff'][i]
            loggs[k,i] = isoc['LogG'][i]
            mass[k,i] = isoc['M/Mo'][i]
            lum[k,i] = 10**isoc['LogL/Lo'][i]
            rad[k,i] = np.sqrt((6.67408e-11 * isoc['M/Mo'][i]*1.988e30)/((10**isoc['LogG'][i] /100)))/6.957e8

            imagslist[k,i] = imags

    
    probabilities = np.exp(-0.5*(chi2/np.min(chi2) - 1 ))
    probabilities = np.where(np.isnan(probabilities),0, probabilities)
    probabilities = np.clip(probabilities,0,np.percentile(probabilities[probabilities > 0.001],99))
    probabilities = probabilities/np.nansum(probabilities)

#
#
#
#    plt.figure(figsize = (9,9), dpi = 200, layout = "constrained")
#    #plt.title("2D Isochrone Grid ([M/H]=0.3)")
    image = np.clip(np.log10(probabilities),-int(np.log10(np.count_nonzero(probabilities))), 6)
#    plt.imshow(image, origin = "lower", aspect = "auto", interpolation = "none")
#    plt.colorbar(label = "log10(P)")
#    #plt.contour(np.where(chi2 < twoSigmaChi2s[len(emags)-2]*np.min(chi2),1,0) + np.where(chi2 < oneSigmaChi2s[len(emags)-2]*np.min(chi2),1,0), origin = "lower", aspect = "auto", vmin = 0, vmax = 2, cmap = "Greens", interpolation = "none")
#    plt.contour(np.where(chi2 == 9e9,1,0), cmap = "Reds", interpolation = "none", levels = 1, vmin = -1e9,vmax = 1, origin = "lower")
#    plt.xlim(0,800)    
#    plt.xlabel("EEP")
#    plt.ylabel("Age (Gyr)")
#    plt.yticks(range(len(availages)), labels = ["%.1f" % i for i in availages])
    #plt.colorbar(label = "L")




    minPoint = np.where(chi2 == np.min(chi2))
    minPoint = (minPoint[0][0], minPoint[1][0])
    
    samples = rejection_sampling(probabilities, 2000)
    samplesx = [i[0] for i in samples]
    samplesy = [i[1] for i in samples]
    #st[samplesx,samplesy]
    #sloggs = loggs[samplesx,samplesy]


    isoc = MISTreadiso(availages[minPoint[0]],metal)
    
    imags = np.array([])
        
    for j in range(len(bands)):
        imags = np.append(imags,isoc[bands[j]][minPoint[1]])
    
#
#    plt.figure(figsize = (9,6), dpi = 200)
#    plt.plot(waves,imags, marker = 'o', label= "Best fit isochrone", lw =3, color = "#00beef")
#    for j in range(len(bands)):
#        plt.annotate(bands[j],(waves[j],imags[j]), (waves[j]+10,imags[j]+0.03))
#    plt.ylabel("Absolute Magnitude")
#    plt.xlabel("Filter Wavelength (nm)")
    
    #TODO : Check the labels on these fellas
#    plt.errorbar(waves,smags, yerr = emags, marker = "o", lw = 3, label= "Stellar Photometry", color = "#042069")
#    plt.legend()
#    
#    plt.figure()
#    plt.hist(availages[samplesx], bins = availages)
#    plt.xlim(0,13)
#    plt.xlabel("Age (Gyr)")
#    plt.show()

    fig = plt.figure(dpi = 200, figsize = (15,15), layout = "tight")
    gs = GridSpec(3,3,figure = fig)
    ax = []
    #plt.tight_layout(h_pad = -23.5, w_pad = -23.5)
    ax.append(fig.add_subplot(gs[1:3,0:2]))
    ax.append(fig.add_subplot(gs[0:1,0:2]))

    im=ax[0].imshow(image, origin = "lower", aspect = "auto", interpolation = "none")
    fig.colorbar(im, ax = ax[0], aspect = 30, label = "log10(P)")

    colors = ["#0000", "#888"]
    cmap = LinearSegmentedColormap.from_list("mylist", colors, N=2)
    ax[0].imshow(np.where(chi2 == 9e9,1,0), aspect = "auto", cmap = cmap, interpolation = "none", vmin = 0, vmax = 1, origin = "lower")
    ax[0].set_xlim(0,600)
    ax[0].set_xlabel("EEP")
    ax[0].set_ylabel("Age (Gyr)")
    ax[0].set_title("Probability Grid")
    ax[0].set_yticks(range(len(availages)), labels = ["%.1f" % i for i in availages])
    ax[1].plot(waves,imags, marker = 'o', label= "Best-fit Grid Point", lw =3, color = "#222")
    for j in range(len(bands)):
        ax[1].annotate(bands[j],(waves[j],imags[j]), (waves[j]-75,imags[j]-0.1))
    for j in range(100):
        ax[1].plot(waves, imagslist[samplesx[j],samplesy[j]], marker = 'o', lw = 3, color = "#222", alpha = 0.01)
    ax[1].set_ylabel("Absolute Magnitude")
    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    ax[1].set_xlabel("Filter Wavelength (nm)")
    ax[1].set_title("Photometric Fit")
    #TODO : Check the labels on these fellas
    ax[1].errorbar(waves,smags, yerr = emags, marker = "o", lw = 3, label= "Stellar Photometry", color = "#042069")
    ax[1].legend()

    ax.append(fig.add_subplot(gs[0:1,2:3]))
    ax.append(fig.add_subplot(gs[1:2,2:3]))
    ax.append(fig.add_subplot(gs[2:3,2:3]))
    ax[2].hist(availages[samplesx], bins = availages)
    ax[2].set_xlabel("Age (Gyr)")
    ax[2].set_xscale("log")
    ax[2].set_xticks([0.1,0.2,0.5,1,2,5,10], labels = ["%.3g" % i for i in [0.1,0.2,0.5,1,2,5,10]] )
    ax[2].set_xticks([],[], minor = True)
    ax[3].hist(teffs[samplesx,samplesy], bins = np.linspace(np.floor(np.min(teffs[samplesx,samplesy])/100)*100,np.ceil(np.max(teffs[samplesx,samplesy])/100)*100,20))
    ax[3].set_xlim(np.floor(np.min(teffs[samplesx,samplesy])/100)*100,np.ceil(np.max(teffs[samplesx,samplesy])/100)*100)
    ax[3].set_xlabel("Teff (K)")
    ax[4].hist(loggs[samplesx,samplesy], bins = np.linspace(np.floor(np.min(loggs[samplesx,samplesy])*10)/10,np.ceil(np.max(loggs[samplesx,samplesy])*10)/10,20))
    ax[4].set_xlim(np.floor(np.min(loggs[samplesx,samplesy])*10)/10,np.ceil(np.max(loggs[samplesx,samplesy])*10)/10)
    ax[4].set_xlabel("log(g) [cgs]")



    avTeff = np.average(teffs[samplesx,samplesy])
    avMass = np.average(mass[samplesx,samplesy])
    avLogg = np.average(loggs[samplesx,samplesy])
    avLum = np.average(lum[samplesx,samplesy])
    avRad = np.average(rad[samplesx,samplesy])

    stdTeff = np.std(teffs[samplesx,samplesy])
    stdLogg = np.std(loggs[samplesx,samplesy])
    stdMass = np.std(mass[samplesx,samplesy])
    stdLum = np.std(lum[samplesx,samplesy])
    stdRad = np.std(rad[samplesx,samplesy])

    ax[2].set_title("Age = %.1f +/- %.1f Gyr" % (np.average(availages[samplesx]), np.std(availages[samplesx])))
    ax[3].set_title("Teff = %i +/- %i K" % (avTeff, stdTeff))
    ax[4].set_title("logg = %.2f +/- %.2f" % (avLogg, stdLogg))
    fig.suptitle("Isochronal Parameters for %s\n ([M/H] = %.1f)"% (target,metal), fontsize = 36, fontweight = "bold")
    plt.savefig("aaa.png")
    plt.close("all")

    #print("Lowest Chi^2 Value: %.2f" % np.min(chi2))
    
    return stellarParameters(avTeff, stdTeff, avLogg, stdLogg, avMass, stdMass, avRad, stdRad, avLum, stdLum, metal)

def createLDFileName(teff,logg,metal,band):
    
    availmetals = np.array([1,0.5,-0.0,-0.5,-1,-1.5,-2,-3,-4])
    
    if teff > 4900 and teff < 5000:
        teff = 4900
    if teff >= 5000 and teff < 5100:
        teff = 5100


    filename = "PHOENIX/LD/"+band+"/Z"
    #if metal == 0: metal = -0.0
    filename += "%+.1f" % availmetals[np.argmin(np.abs(availmetals - metal))]
    filename += "/lte%05d" % round(teff, -2)
    filename += "%+4.2f" %(-round(logg*2,0)/2)
    filename += "%+3.1f" % availmetals[np.argmin(np.abs(availmetals - metal))]
    filename += ".PHOENIX-ACES-AGSS-COND-LIMBDARK-2011-"
    if band in ["U","B","V","R","I"]:
        filename += "Johnson_"
    filename += band
    filename += ".fits"
    return filename

def fitLDCurve(teff,logg,metal,band):
    ld = fits.open(createLDFileName(teff,logg,metal,band))
    fit = linregress(ld[1].data[ld[1].data > 0.1],ld[0].data[ld[1].data > 0.1]/np.max(ld[0].data))
    return np.nanmax([fit.slope,0])


def readAbundances(workDir, paramFile, elements):
    with open(workDir + paramFile) as f:
        lines = f.readlines()


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



    dirList = os.listdir(workDir)
    abundances = np.array([])
    errorbars = np.array([])

    with open(workDir + paramFile) as f:
        lines = f.readlines()

    for i in elements:
        try:
            summary = np.genfromtxt(workDir + "%i" % i + "/abundanceSummary.txt", names = ["wav", "[X/H]", "v_broad", "chi2"], ndmin = 2)

            if len(summary) > 15:
                summary = summary[summary["chi2"] < 3]
            
            abundances = np.append(abundances, np.median(summary["XH"]))
            errorbars = np.append(errorbars, np.std(summary["XH"])/np.sqrt(len(summary)))

        except:
            abundances = np.append(abundances, 0)
            errorbars = np.append(errorbars, 0)
    for i in range(len(elements)):
        if abundances[i] == 0:
            abundances[i] = solarMOOG[int(elements[i])]
    return elements, abundances, errorbars

def readInputMetallicity(paramFile):
    with open(paramFile) as f:
        lines = f.readlines()

    for line in lines:
        if "[M/H]" in line:
            return float(line[10:20]) 


def readInputAlphaEnhancement(paramFile):
    with open(paramFile) as f:
        lines = f.readlines()

    for line in lines:
        if "[alpha/M]" in line:
            return float(line[10:20]) 



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


