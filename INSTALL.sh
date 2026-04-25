#! /bin/bash



# IF YOU NEED TO SET G++ TO LINK TO GSL, DO SO HERE ###########################

# As an example, this is the correct setup for Mac OS 
# if GSL was installed with Homebrew
gppargs=""



# DO NOT CHANGE ANYTHING BELOW THIS POINT #####################################




# Unique requirements for Linux

autogppargs=""
if [[ $OSTYPE == *"linux"* ]];
then
    # Check that GCC is installed
    echo ""
    gcc --version | grep "gcc"
    if [ $? -ne 0 ];
    then
        echo "Check that you have gcc installed on your system"
        echo "Try running the following command:"
        echo ""
        echo "sudo apt install gcc"
        echo ""
        echo "and then re-run the install script"
        exit 1
    fi;

    # Check that GSL is installed
    echo ""
    echo "Package: libgsl-dev"
    dpkg -s libgsl-dev | grep "installed"
    if [ $? -ne 0 ];
    then
        echo "Searching for GSL..."
        find / -iname gsl 2>&1 | grep -v "Permission denied\|Invalid argument\|anaconda\|downloads\|Downloads"
        if [ $? -ne 0 ];
        then
            echo "Check that you have GSL installed on your system"
            echo "Try running the following command:"
            echo ""
            echo "sudo apt install libgsl-dev"
            echo ""
            echo "and then re-run the install script"
            exit 1
        else
            echo "GSL was found, but it was not installed using apt."
            echo "Make sure you set the right g++ compiler flags to point it to where GSL is installed"
        fi;
    else
        gppargs=" "
    fi;
    
    echo "" 

    





# Unique requirements for Mac OS
elif [[ $OSTYPE == *"darwin"* ]];
then
    # Check that coreutils is installed
    echo ""
    timeout 3 echo "coreutils installed"
    if [ $? -ne 0 ];
    then
        echo "Check that you have coreutils installed on your system"
        echo "Try running the following command:"
        echo ""
        echo "brew install coreutils"
        echo ""
        echo "and then re-run the install script"
        exit 1
    fi;  

    # Check that gcc is installed
    echo ""
    gcc --version | grep "gcc"
    if [ $? -eq 127 ];
    then
        echo "Check that you have gcc installed on your system"
        echo "Try running the following command:"
        echo ""
        echo "brew install gcc"
        echo ""
        echo "and then re-run the install script"
        exit 1
    fi;

    # Check that GSL is installed
    echo ""
    echo "Package: libgsl-dev"
    brew list gsl | grep "No such keg"
    if [ $? -eq 0 ];
    then
        echo "Searching for GSL..."
        find / -iname gsl 2>&1 | grep -v "Permission denied\|Invalid argument\|anaconda\|downloads\|Downloads"
        if [ $? -ne 0 ];
        then
            echo "Check that you have GSL installed on your system"
            echo "Try running the following command:"
            echo ""
            echo "brew install gsl"
            echo ""
            echo "and then re-run the install script"
            exit 1
        else
            echo "GSL was found, but it was not installed using apt."
            echo "Make sure you set the right g++ compiler flags to point it to where GSL is installed"
        fi;
    else
        autogppargs=" -I/opt/homebrew/include -L/opt/homebrew/lib"
    fi;


    # Check that there is no active conda environment (Apparently this breaks gfortran for some reason)
    conda env list | grep "*" -q
    if [ $? -eq 0 ];
    then
        echo "There is a conda environment active. This will mess up the installation process."
        echo "Please run 'conda deactivate' before proceeding. You are free to reactivate your conda environment once the installation script has been run"
        exit 1
    fi;
fi;

if [[ $gppargs == "" ]];
then
    gppargs=$autogppargs
fi;


# Compile RunAbundanceOnGoodLines.cpp 
echo "Compiling C++ code..."
echo ""
g++ RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -pthread -O2 -std=c++17 $gppargs
if [ $? -ne 0 ];
then
    g++ RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -pthread -O2 -std=c++2a $gppargs
    if [ $? -ne 0 ];
    then
        echo "Something went wrong when compiling the C++ code"
        echo ""
        echo "Check that your version of gcc is recent enough to support the c++20 standard"
        echo ""
        echo "Also check that you have g++ installed and that you've pointed the compiler to the correct directory where libgsl is located, using the compiler flag -L[/path/to/libgsl]"
        exit 1
    fi;
fi;





# Try installing python packages if they aren't already
pip show numpy | grep "Name\|Version"
if [ $? -eq 1 ];
then
    echo ""
    echo "Installing numpy..."
    pip install numpy
fi;
echo ""
echo "NumPy installed"


pip show astropy | grep "Name\|Version"
if [ $? -eq 1 ];
then
    echo ""
    echo "Installing astropy..."
    pip install astropy
fi;
echo ""
echo "AstroPy installed"

pip show scipy | grep "Name\|Version"
if [ $? -eq 1 ];
then
    echo ""
    echo "Installing scipy..."
    pip install scipy
fi;
echo ""
echo "SciPy installed"

pip show astroquery | grep "Name\|Version"
if [ $? -eq 1 ];
then
    echo ""
    echo "Installing astroquery..."
    pip install astroquery
fi;
echo ""
echo "AstroQuery installed"






# Compile linemake 
echo ""
echo "Compiling linemake..."
gfortran linemake/linemake.f -o linemake/linemake -O2 -ffixed-line-length-none
if [ $? -ne 0 ];
then
    echo "Something went wrong when compiling Linemake"
    echo "Check that you have gfortran installed"
    exit 1
fi;





# Compile MOOG 
cd moog
rm *.o
make | grep "gfortran"
cd ..





# Download and extract large data files
if [ ! -d "mooglists" ];
then
    echo ""
    echo "Downloading mooglists"
    curl "https://drive.usercontent.google.com/download?id=1zZEePXb1mW17dlddEtWwYrrRPksNKfek&confirm=xxx" -o mooglists.tar.gz
    echo "Extracting mooglists..."
    tar -xf mooglists.tar.gz
    rm mooglists.tar.gz
fi;

if [ ! -d "isochrones" ];
then
    echo ""
    echo "Downloading isochrones"
    curl "https://drive.usercontent.google.com/download?id=1FS8aIk9e-7hbHdNQTs0e2u-8vTNzM0T9&confirm=xxx" -o isochrones.tar.gz
    echo "Extracting isochrones..."
    tar -xf isochrones.tar.gz
    rm isochrones.tar.gz
fi;

if [ ! -d "PHOENIX" ];
then
    echo ""
    echo "Downloading PHOENIX data"
    curl "https://drive.usercontent.google.com/download?id=17Mpy5yd1zTaYQgtBdA_H2ae8P17CWLvM&confirm=xxx" -o PHOENIX.tar.gz
    echo "Extracting PHOENIX data..."
    tar -xf PHOENIX.tar.gz
    rm PHOENIX.tar.gz
fi;





# Exit successfully
echo ""
echo "DONE"
echo ""
echo "Try executing the following command:"
echo ""
echo "./runStar.sh \"tau Ceti\" example-tauCet/ -0.41"
echo ""
exit 0
