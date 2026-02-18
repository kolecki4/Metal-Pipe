# Metal Pipe
Metal Pipe is designed to be a highly accurate line-by-line stellar spectrum fitting code fit for use with as many wavelength ranges and stellar spectral types as possible.

Future work will focus on adding more line lists, especially for the near-infrared wavelength range, as support is extended to late-K and M type stars.

## Installation

### 0. Basic Prerequisites
The installation process requires the `gfortran` and `g++` compilers to be installed on your system. Your version of `g++` should be recent enough to support the C++17 standard.

If you are running Linux, you'll need to have `glibc` installed as well.

On Mac OS you'll need to install `coreutils` with homebrew and also make sure you have XCode installed.

### 1. GNU Scientific Library

This is the C++ library on which the core line-fitting algorithm is based. On my Ubuntu-based machine, I was simply able to run `sudo apt install libgsl-dev` and have everything done for me. If you have a Mac, you may run `brew install gsl`. 

On the other hand, if your favorite package manager does not have a GSL package ready, or you do not have write access to `/usr/include`, you will need to compile this package from source. Doing so in most cases should install GSL in the default (correct) location

**NICHE ISSUE:** Sometimes the install script will fail on compiling the C++ code complaining about multiple declarations of `isnan` and `isinf` functions. This is an issue as of the latest stable release of GSL (v2.8). To fix it, go to the ...`/include/gsl` directory, open `gsl_math.h`, and change line 22 from `#include <math.h>` to `#include <cmath>`.

### 2. Determine how many cores/threads your system has, and how many you want Metal Pipe to have access to
If you don't know how many threads your CPU has, you can look up the CPU model name online and find the information. Alternatively, try runnning `lscpu | grep "Model name\|CPU(s):\|Thread(s)\|Core(s)\|Socket(s)"` from the command line. This will give you a number of threads, listed as "CPU(s)" and tell you your CPU model name. 

Once you're aware of how many threads you'll be able to utilize, set the maximum number available to Metal Pipe by opening `RunAbundanceOnGoodLines.cpp` in a text editor and change the value of `MAX_THREADS` (on line 18) to be your desired value

### 3. Run INSTALL.sh
- If you installed GSL using apt, you can run INSTALL.sh as is.

- If you installed GSL using Homebrew, change line 9 of INSTALL.sh to the following:
   - `gppargs="-I/opt/homebrew/include -L/opt/homebrew/lib"`

- Otherwise, you're on your own determining where GSL was installed to. Change line 9 of INSTALL.sh to be identical to above. Then, change the -I flag to point to your include folder and change the -L flag to point to your lib folder 

- If you are running Mac OS with a conda environment activated (as in, you see "(base)" listed on your command line), you must `conda deactivate` before you run the install script. You can safely `conda activate` afterwards.



## Running Metal Pipe

### 1. Prepare the Stellar Spectrum
The observed stellar spectrum file **absolutely must:**
- Be normalized (continuum flux = 1 at all wavelengths)
- Be corrected for stellar RV and barycentric velocity 
- Be stored in plain-text format with one data point per line, in the format of `wavelength flux errorbar`. (If you don't have exact error bars, you can fill this column in with a constant value equal to the inverse of your signal-to-noise ratio. A non-zero error bar is required for the chi^2 fitting to work properly)
- Have wavelengths in air wavelengths (not vacuum wavelengths) and in units of Angstroms


As of now, if you name the spectrum file with the string "flattened" in it, the python script will automatically pick it up and add it to the parameters file.

**It is highly recommended** that you create a unique folder for each stellar spectrum, as Metal Pipe will dump output files to the same directory the spectrum is in.

### 2. Choose Line Lists
In its default configuration, Metal Pipe is built to fit lines of Fe, Ca, Ti, Mg, Si, C, O, Na, Al, and K, in that order, with a range of wavelengths chosen to coincide with KPF. These line lists are stored in the appropriately named `linelists` folder. This folder contains an additional file, `defaults.txt` which tells Metal Pipe what files should be used to fit which elemental abundances.

If this all sounds good to you (you're interested in fitting the stated elements using KPF spectra), great! You can stop reading the rest of this section and go to the next step.

...

Otherwise, you'll have to pick what lines you want to fit. I have some helper scripts for this which will be available on written request (they're ugly and not ready to be public just yet, but they work). You're also welcome to simply pick out your favorite lines, though the wavelength data must match the wavelength data provided in linemake's repository of lines, otherwise you risk Metal Pipe just not finding your lines. I recommend combing through linemake's `mooglists` directory and choosing the lines from there. This way, you can make sure that a) linemake actually has data for the lines you want, and b) the wavelengths you're telling Metal Pipe the lines are located at are actually correct.


### 3. Run `./runStar.sh`
First, open up a terminal in the Metal Pipe directory

The bash script requires a few command line arguments to kick everything off, but will carry the code from there to its completion. In order, they are as follows:

1. The name of your star, in a format resolvable by SIMBAD. This is used by the python script to get photometry and thus get stellar parameters
2. The folder which contains the spectrum file. This will be where all output files are dumped, and they are indeed dumped, without care for what else might have been in said folder, hence the recommendation to create a unique directory for each spectrum you wish to get abundances for. **NOTE: YOU MUST HAVE A SLASH AT THE END OF THIS FILEPATH, AS SHOWN BELOW.**
3. An initial guess at the metallicity of your star (\[M/H\]). The solar value of 0.0 is always a safe place to start, but you might wish to modify it, as a better guess can reduce the number of iterations required to fit the stellar parameters, and can significantly affect execution time

Here is an example of what those command line arguments look like:

`./runStar.sh "HD 10700" /home/jared/Documents/Spectra/KPF/CAP4/10700/ -0.41`


## Understanding the Output
Metal Pipe creates a directory for each element, named after its atomic number. In these folders, a text file is created for each line which was successfully fit. 

Generally, however, you will be most interested in the "params.txt" file. This file is used as both an input and an output file, and contains all the information you're likely to want out of Metal Pipe at a glance:

- Stellar parameters:
  - $T_{\textrm{eff}}$, $\log{g}$, $[M/H]$, $[\alpha/M]$, $M/M_\odot$, $R/R_\odot$, and $L/L_\odot$. Each is listed with their values and 1-sigma uncertainties.

- Limb darkening coefficients:
  - Reported for a linear limb darkening model, with the name of the photometric band, wavelength in angstroms, and the coefficient itself.

- Elemental Abundances:
  - Labeled according to the atomic number of their element, and output in the $[X/H]$ format (abundance relative to solar, where $X$ is any given element e.g. Fe, C, etc.). The rightmost column is again 1-sigma error bars.

#### Note:
If one of these error bars seems unusually large, you can try to go into the corresponding element's folder, open "abundanceSummary.txt" and remove lines from it that are significantly and obviously poor fits. This text file may be difficult to parse as of now, but contains best-fit data for each spectral line fit in the format 'wavelength abundance v_broad chi^2_nu' with one data point per line. Removing offensively bad lines from this file and then running 

`python computeParamFile.py "HD 10700" /home/jared/Documents/Spectra/KPF/CAP4/10700/ params.txt -0. 0. 4`

should improve things. 

**IMPORTANT: Before you do the above** As for runStar.sh, replace "HD 10700" with the name of your star, and the directory path with the path to your star's directory. Leave params.txt as is, but replace the two zeros in the above line with the \[M/H\] and \[alpha/M\] values respectively, from params.txt. 
##### I'm aware this is an extremely clunky way to perform manual sigma clipping of the abundances. I am certainly open to suggestions for improvement


# Acknowledgements
This code is directly based on multiple other codebases which have been written by other members of the astronomy community. Without them my work would not be possible.
- [PyMOOGi](https://github.com/madamow/pymoogi) (Monika Adamow) and of course [MOOG](https://www.as.utexas.edu/~chris/moog.html) (Chris Sneden)
- The line list data and the code for its manipulation: [Linemake](https://github.com/vmplacco/linemake) (Vinicius Placco)
- The atmosphere interpolation program which was based on PyKMOD, which itself was based on [KMOD](http://leda.as.utexas.edu/stools/) (Carlos Allende Prieto)
- The Python packages [Astropy](https://www.astropy.org/), [SciPy](https://scipy.org/), and [NumPy](https://numpy.org/)


