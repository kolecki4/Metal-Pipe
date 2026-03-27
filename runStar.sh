#! /bin/bash
ARGC=$#
#echo " ╔═════════════════════════════════════════════════════════════════════════════╗"
#echo " ║     *__________* ___    __  * ______    ______  * ______  *  _      _     * ║"
#echo " ║     ___________ │   ╲  ╱  ╲  ╱  __  ╲ *╱  __  ╲  ╱  ____╲  _│ │_ *_│ │_     ║"
#echo " ║    ____________ │  \ ╲╱    ││  ╱  ╲  ││  ╱  ╲  ││  ╱  ___ │_\`  _││_\`  _│    ║"
#echo " ║ * _____________ │  │╲  ╱│  ││ │    │ ││ │    │ ││ │  │_  │  │_│    │_│      ║"
#echo " ║  ______________ │  │ ╲╱ │  ││  ╲__╱  ││  ╲__╱  ││  ╲__╱  │    *     *       ║"
#echo " ║ _______________ │__│ *  │__│ ╲______╱* ╲______╱  ╲______╱ *            *    ║"
#echo " ║   *        *                 *                  *               *           ║"
#echo " ╟─────────────────────────────────────────────────────────────────────────────╢"
#echo " ║ ░░░░░░░▒▒▒▒▒▒▒▓▓▓▓▓▓▓████████ Beta Version 1.1 ███████▓▓▓▓▓▓▓▒▒▒▒▒▒▒░░░░░░░ ║"
#echo " ╚═════════════════════════════════════════════════════════════════════════════╝"
#
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
echo "║ ░░░░░░░▒▒▒▒▒▒▒▓▓▓▓▓▓▓███████ Beta Version 1.5.1 ███████▓▓▓▓▓▓▓▒▒▒▒▒▒▒░░░░░░░ ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"

n=$(($RANDOM % 100))
#echo $n
if [ $n -lt 10 ];
then
    ffplay -autoexit -nodisp "Destruction_Metal_Pole_L_Wave_2_0_0.wav" &>/dev/null &

fi;

if [ $ARGC -ne 3 ];
then
    echo "USAGE: ./runStar.sh [Star Name] [Working Directory] [Initial Metallicity]"
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
    
else
    starName=$1
    workDir=$2
    initMetal=$3

    i=0
    python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
    statusComputeParams=$?   

    if [[ $statusComputeParams -ne 0 ]];
    then
        echo "There was a problem executing computeParamFile.py"
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
            exit 1;
        fi;
        ((i++))
    done
    python correctForNLTE.py ${workDir} O
    python correctForNLTE.py ${workDir} S
    #/scr/jkolecki/miniconda3/bin/python3 correctForNLTE.py ${workDir} K
    python correctForNLTE.py ${workDir} Ca
    python correctForNLTE.py ${workDir} Ti
    python plotSomeLines.py ${workDir}
fi;
