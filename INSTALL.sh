#! /bin/bash



# IF YOU NEED TO SET G++ TO LINK TO GSL, DO SO HERE ###########################


gslflags=""
  
# DO NOT CHANGE ANYTHING BELOW THIS POINT #####################################

echo """
    Welcome to the install script for Metal Pipe! 
    To successfully install Metal Pipe on your system, 
    a few things will need to be in place before we start.

    Metal Pipe Prerequisites:
    ╠── Compilers
    ║   ╠── gfortran version >4.0
    ║   ╠── g++ version >7.0
    ║   ╚── XCode (MacOS only)
    ╠── Python environment
    ║   ╠── Python version >3.11
    ║   ╚── Required python packages can be installed with: 
    ║       'pip install pandas numpy matplotlib scipy astropy astroquery'
    ╠── Additional libraries
    ║   ╠── GSL (GNU Scientific Library)
    ║   ║   ╠── On Ubuntu, try 'sudo apt install libgsl-dev'
    ║   ║   ╚── On MacOS, try 'brew install gsl'
    ║   ╠── GLIBC 
    ║   ║   ╠── On Ubuntu, this should be installed already
    ║   ║   ╚── On MacOS, try 'brew install glibc'
    ║   ╚── coreutils (MacOS only)
    ║       ╚── On MacOS, try 'brew install coreutils'
    ╚── Final checks
        ╠── You MUST disable your python environment before running  
        ║   the install script (e.g. by running 'conda deactivate'). 
        ║   You can safely 'conda activate' again once it's done.
        ╚── You MUST verify that 'gslflags' (line 8 of INSTALL.sh)
            points to the location of the gsl library on your system. 
            If you installed gsl with apt, you can leave this blank.
            If you installed gsl with brew on MacOS, change it to
            '-I/opt/homebrew/include -L/opt/homebrew/lib'. Elsewise,
            you'll have to figure this one out yourself.
"""

userConfirm=""
while [[ !(${userConfirm,,} =~ ((^yes$)|(^y$)|(^no$)|(^n$))) ]];
do
    echo "I confirm that these prerequisites are met. Compile Metal Pipe (yes/no?): "
    read userConfirm

    if [[ ${userConfirm,,} =~ ((^yes$)|(^y$)) ]];
    then
        continue
    elif [[ ${userConfirm,,} =~ ((^no$)|(^n$)) ]];
    then
        exit 1
    fi

done;


# Compile RunAbundanceOnGoodLines.cpp 
echo "Compiling C++ code..."
echo ""
g++ RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -pthread -O2 -std=c++17 $gslflags
if [ $? -ne 0 ];
then
    g++ RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -pthread -O2 -std=c++2a $gslflags
    if [ $? -ne 0 ];
    then
        echo "Something went wrong when compiling the C++ code"
        echo ""
        echo "Check that your version of gcc is recent enough to support the c++17 standard"
        echo ""
        echo "Also check that you have g++ installed and that you've pointed the compiler to the correct directory where libgsl is located, using the compiler flag -L[/path/to/libgsl]"
        exit 1
    fi;
fi;




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

if [ ! -d "photometryCache" ];
then
    mkdir photometryCache
fi;
if [ ! -f "photometryCache/HD - SIMBADxGaia.csv" ];
then
    echo ""
    echo "Downloading photometry cache for all HD stars"
    curl "https://drive.usercontent.google.com/download?id=1Vg65ggTl1zNstwRclvHFIuCTao-l46fE&confirm=xxx" -o "photometryCache/HD - SIMBADxGaia.csv.tar.gz"
    echo "Extracting photometry cache..."
    tar -xf "photometryCache/HD - SIMBADxGaia.csv.tar.gz" -C photometryCache
    rm "photometryCache/HD - SIMBADxGaia.csv.tar.gz"
fi;






# Exit successfully
echo ""
echo "DONE"
echo ""
echo "Try executing the following command:"
echo ""
echo "./runStar.sh \"tau Cet\" example-tauCet/ -0.41 KPF"
echo ""
exit 0
