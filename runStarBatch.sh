ARGC=$#
if [ $ARGC -ne 1 ];
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

    yourfilenames=`ls $1`
    for folder in "$1"/*/
    do
        needsRun=1       
        if [[ $folder != *"."*  &&  $folder != *"plots"* ]];
        then 

            for file in "$folder"/*
            do
            if [[ $file == *"params.txt" ]];
            then
                needsRun=0
            fi;
            done




            starname=${folder/$1/""}
            starname=${starname//"/"/""}

            #echo "${starname}"

            if [[ $starname =~ ^[0-9]+$|^[0-9]+[A-D]$ ]];
            then            
                starname="HD_$starname"
           
            elif [[ $starname =~ ^T00[0-9]+$ ]];
            then            
                starname=${starname/"T00"/"TOI-"}

            elif [[ $starname =~ ^K00[0-9]+$ ]];
            then            
                starname=${starname/"K00"/"KOI-"}
            elif [[ ${starname^^} =~ ((^HD)|(^GJ)|(^GL))([0-9]+$|[0-9]+[A-D]$) ]];
            then
                starname=${starname^^}
                starname=${starname:0:2}"_"${starname:2}

            elif [[ ${starname^^} =~ (^HIP)|(^KIC)|(^TIC)([0-9]+$|[0-9]+[A-D]$)} ]]
            then
                starname=${starname^^}
                starname=${starname:0:3}"_"${starname:3}

            elif [[ ${starname^^} =~ ((^KEPLER-)|(^KEPLER)|(^KOI-)|(^TOI-)|(^K2-)|(^KELT-)|(^WASP-))([0-9]+$|[0-9]+[A-D]$) ]];
            then
                starname=$starname
            else
                #echo \"$starname\" is not a valid name for a star. I\'m ignoring this folder.
                echo Folder name \"$starname\" assumed SIMBAD-resolvable
        	#needsRun=1
            fi;
            
            if [[ $needsRun == 1 ]];
            then
               ./runStar.sh ${starname} ${folder} 0
            fi;
        fi;
    done
#    starName=$1
#    workDir=$2
#    initMetal=$3
#
#    i=0
#    python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
#    notConverged=1
#    problem=0
#    ((i++))
#    until [[ $i -gt 0 && notConverged -eq 0 ]];
#    do
#        
#
#        #valgrind -s --leak-check=full --show-leak-kinds=all --verbose --track-origins=yes --log-file=valgrind-out.txt 
#        ./RunAbundanceOnGoodLines ${workDir}params.txt 
#        notConverged=$?
#        python computeParamFile.py "${starName}" ${workDir} params.txt ${initMetal} 0.0 ${i}
#        problem=$?
#        if [[ $problem -ne 0 ]];
#        then
#            return 1;
#        fi;
#        ((i++))
#    done
#    python correctForNLTE.py ${workDir} O
#    python correctForNLTE.py ${workDir} S
#    python correctForNLTE.py ${workDir} K
#    python plotSomeLines.py ${workDir}
fi;
