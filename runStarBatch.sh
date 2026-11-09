ARGC=$#
if [ $ARGC -ne 2 ];
then
    echo """

    This script runs a batch of stellar spectrum folders through Metal Pipe.

    USAGE: ./runStar.sh [Working Directory] [Instrument]

    - [Working Directory]: A full path to your folder of star folders
      The final character of this path MUST BE A SLASH. Metal Pipe will be run
      on every folder in this directory, attempting to translate folder names 
      into valid star names.

    - [Instrument]: This dictates which solar abundance corrections are used.
      Valid choices: ESPADONS, ESPRESSO, FEROS, HIRES, KPF, NEID, PEPSI.
    """
    
else
    # Check if user-input instrument choice is valid
    inst=$2
    if [[ !(${inst^^} =~ ((^ESPADONS$)|(^ESPRESSO$)|(^FEROS$)|(^HIRES$)|(^KPF$)|(^PEPSI$)|(^NEID$))) ]];
    then
        echo No line-by-line solar abundances from ${inst}. Metal Pipe cannot run.
        exit 1;
    fi;

    # Check if user-input directory exists
    if [ ! -d $1 ];
    then
        echo No such directory: $1. Metal Pipe cannot run.
        exit 1;
    fi;


    # For folder in user-input directory, check if it contains a 
    # correctedParams.txt file (signifying an already-complete run
    # of Metal Pipe). If not, then attempt to resolve the
    # folder name into a star name, and run Metal Pipe on it
    yourfilenames=`ls $1`
    for folder in "$1"/*/
    do
        needsRun=1       
        if [[ $folder != *"."*  &&  $folder != *"plots"* ]];
        then 
            # Check if folder/correctedParams.txt exists. If it does,
            # the folder doesn't need to be re-run through Metal Pipe
            for file in "$folder"/*
            do
            if [[ $file == *"correctedParams.txt" ]];
            then
                needsRun=0
            fi;
            done

    
            if [[ $needsRun == 1 ]];
            then

                # Parse the folder name into a star name (adding spaces
                # and hyphens where necessary, reformatting common CLS
                # IDs into the actual ones, etc.)
                starname=${folder/$1/""}
                starname=${starname//"/"/""}

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
                    echo "Folder name \"$starname\" assumed SIMBAD-resolvable (or Solar spectrum)"
                fi;

                # Run runStar.sh on the folder
               ./runStar.sh ${starname} ${folder} 0 ${inst}
            fi;
        fi;
    done
fi;
