#! /bin/bash
ARGC=$#

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                              ║"
echo "║                                               ,,,,,,,,;;;;''\                ║"
echo "║                                 ,,,,,,,;;;;'''''''           |               ║"
echo "║                 _,,,,,,,;;;;''''''''                        /                ║"
echo "║                / \'''                        ,,,,,,;;;;'''''                 ║"
echo "║               |   |            ,,,,,,,;;;;;'''''                             ║"
echo "║                \_/,,,,,;;;''''''''                                           ║"
echo "║                                                                              ║"
echo "╟──────────────────────────────────────────────────────────────────────────────╢"
echo "║ ░░░░░░░▒▒▒▒▒▒▒▓▓▓▓▓▓▓███████ Beta Version 2.0.0 ███████▓▓▓▓▓▓▓▒▒▒▒▒▒▒░░░░░░░ ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"


if [ $ARGC -ne 4 ];
then
    echo """

    This script runs a single stellar spectrum through Metal Pipe. 

USAGE: 
./runStar.sh [Star Name] [Working Directory] [Initial Metallicity] [Instrument]

    - [Star Name]: SIMBAD-resolvable star id. Used to query stellar photometry

    - [Working Directory]: A full path to your star's working directory.
      The final character of this path MUST BE A SLASH. This is the folder 
      in which your input data is located and where output files will be dumped.
   
    - [Initial Metallicity]: A \"starting guess\" for the code.
      If you have absolutely no idea, a value of 0 will work fine enough
 
    - [Instrument]: This dictates which solar abundance corrections are used.
      Valid choices are ESPADONS, ESPRESSO, FEROS, HIRES, KPF, NEID, PEPSI.
   
"""
else
    starName=$1
    workDir=$2
    initMetal=$3
    inst=$4
    
    # Check if user-input instrument choice is valid
    if [[ !(${inst^^} =~ ((^ESPADONS$)|(^ESPRESSO$)|(^FEROS$)|(^HIRES$)|(^KPF$)|(^PEPSI$)|(^NEID$))) ]];
    then
        echo No line-by-line solar abundances from ${inst}. Metal Pipe cannot run.
        exit 1;
    fi;

    # Check if workDir exists
    if [ ! -d $workDir ];
    then
        echo No such directory: $workDir. Metal Pipe cannot run.
        exit 1;
    fi;

    # Set 0th iteration
    i=0

    # Try to computeParamFile.py. If failure, note the folder in FailedLog
    python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
    statusComputeParams=$?   
    if [[ $statusComputeParams -ne 0 ]];
    then
        echo "There was a problem executing computeParamFile.py"
        echo $(date): ${starName} ${workDir} Could not create params.txt on iteration ${i} >> FailedLog
        exit 1;
    fi;

    # Until the output abundances match the input parameters, iterate
    statusRunAbundances=1
    ((i++))
    until [[ $i -gt 0 && statusRunAbundances -eq 0 ]];
    do

        # If we've been iterating on this star for a while, give up and note the folder in FailedLog
        if [[ $i -gt 15 ]];
        then
            echo "Too many iterations. I give up."
            echo $(date): ${starName} ${workDir} Reached iteration number limit >> FailedLog
            exit 1;
        fi;


        #valgrind [ARG] -s --leak-check=full --show-leak-kinds=all --verbose --track-origins=yes --log-file=valgrind-out.txt 

        # Run line-fitting procedure on the given spectrum
        # if statusRunAbundances returns 0, then the output is converged with the input
        # if statusRunAbundances returns 1, then we must keep iterating
        ./RunAbundanceOnGoodLines ${workDir}params.txt
        statusRunAbundances=$?

        # Try to computeParamFile.py. If failure, note the folder in FailedLog
        python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
        statusComputeParams=$?
        if [[ $statusComputeParams -ne 0 ]];
        then
            echo "There was a problem executing computeParamFile.py"
            echo $(date): ${starName} ${workDir} Could not create params.txt on iteration ${i} >> FailedLog
            exit 1;
        fi;

        ((i++))
    done

    # Run line by line solar abundance correction script
    python lineByLineSolarCorrections.py ${workDir} ${inst}

    # Make some diagnostic plots of the line fits we just ran
    python plotSomeLines.py ${workDir}
fi;
