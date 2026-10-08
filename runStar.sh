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
    echo "USAGE: ./runStar.sh [Star Name] [Working Directory] [Initial Metallicity] [Instrument]"
    echo ""    
    echo "  - Star Name: Must be resolvable by SIMBAD. 
    Used to query stellar photometry"
    echo ""
    echo "  - Working Directory: A full path to your star's working directory
    The final character of this path MUST BE A SLASH. This is 
    the folder in which your input data is located and where 
    output files will be dumped."
    echo ""    
    echo "  - Initial Metallicity: A \"starting guess\" for the code.
    If you have absolutely no idea, a value of 0 will work fine enough"
    echo ""    
    echo "  - Instrument: Valid choices are the following:
    ESPADONS, ESPRESSO, FEROS, HIRES, KPF, NEID, PEPSI. 
    This dictates which solar abundance corrections will be used."    

else
    starName=$1
    workDir=$2
    initMetal=$3
    if [ $ARGC -eq 4 ];
    then
        inst=$4
        if [[ !(${inst^^} =~ ((^ESPADONS$)|(^ESPRESSO$)|(^FEROS$)|(^HIRES$)|(^KPF$)|(^PEPSI$)|(^NEID$))) ]];
        then
            echo No line-by-line solar abundances from ${inst}. Metal Pipe cannot run.
            exit 1;
        fi;
    fi;
    i=0
    python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
    statusComputeParams=$?   

    if [[ $statusComputeParams -ne 0 ]];
    then
        echo "There was a problem executing computeParamFile.py"
        echo ${starName} ${workDir} Could not create params.txt on iteration ${i} >> FailedLog
        exit 1;
    fi;

    statusRunAbundances=1
    ((i++))
    until [[ $i -gt 0 && statusRunAbundances -eq 0 ]];
    do
        

        #valgrind [ARG] -s --leak-check=full --show-leak-kinds=all --verbose --track-origins=yes --log-file=valgrind-out.txt 
        ./RunAbundanceOnGoodLines ${workDir}params.txt
        statusRunAbundances=$?
        python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
        statusComputeParams=$?
        if [[ $statusComputeParams -ne 0 ]];
        then
            echo "There was a problem executing computeParamFile.py"
            echo ${starName} ${workDir} Could not create params.txt on iteration ${i} >> FailedLog
            exit 1;
        fi;
        if [[ $i -gt 15 ]];
        then
            echo "Too many iterations. I give up."
            echo ${starName} ${workDir} Reached iteration number limit >> FailedLog
            exit 1;
        fi;
        ((i++))
    done
    #python correctForNLTE.py ${workDir} O
    #python correctForNLTE.py ${workDir} S
    #/scr/jkolecki/miniconda3/bin/python3 correctForNLTE.py ${workDir} K
    #python correctForNLTE.py ${workDir} Ca
    #python correctForNLTE.py ${workDir} Ti
    python lineByLineSolarCorrections.py ${workDir} ${inst}
    python plotSomeLines.py ${workDir}
fi;
