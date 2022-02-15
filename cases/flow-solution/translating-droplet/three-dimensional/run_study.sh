#!/bin/bash

print_usage() {
    echo "Usage: ./run_study.sh param_file study_pattern study_runner no_density_eqn" 
    echo "param_file : PyFoam .parameter file."
    echo "study_pattern: string pattern that distinguishes parameter studies."
    echo "study_runner: "
    echo "    slurm (cluster)"
    echo "    foamJob (background)"  
    echo "    commandLine (interactive)" 
    echo "    none (only case generation)."
    echo "no_rho_eqn: if non-empty, runs lentFoam -no-density-equation "
    exit 0
}

if [[ "$1" == "help" || "$1" == "-h" || -z $1 ]]; 
then 
    print_usage
fi

PARAM_FILE=$1
PATTERN=$2
RUNNER=$3
NO_RHO_EQN=$4

if [ ! ${PARAM_FILE##*.} == "parameter" ]; 
then 
    echo "Error: the parameter file does not end with .parameter" 
    print_usage
    exit 1
fi

if [[ ("$RUNNER" == "slurm") || \
      ("$RUNNER" == "foamJob") || \
      ("$RUNNER" == "commandLine") || \
      ("$RUNNER" == "none") ]]; 
then 
    # Source LENT's python modules for parameter variations
    source ../../../../cases/scripts/bashrc
    
    # Prepare and initialize the study directories
    if [ ! -z "$PATTERN" ]; 
    then 
	lent_prepare_study_variants.py "$PARAM_FILE" -p "$PATTERN" -m block  
    else 
	lent_prepare_study_variants.py "$PARAM_FILE" -m block  
    fi

    DIR_PATTERN="" 

    if [ ! -z "$PATTERN" ]; 
    then 
	DIR_PATTERN=$PATTERN-${PARAM_FILE}_00
    else
	DIR_PATTERN=${PARAM_FILE}_00
    fi

    if [ ! "$RUNNER" == "none" ]; 
    then 
	# Run each variant 
	for dir in "$DIR_PATTERN"*; do
	    echo $dir 
	    cd $dir
	    if [ "$RUNNER" == "slurm" ]; 
	    then 
	        echo "Submitting $(pwd)"
	        if [ -z "$NO_RHO_EQN" ]; 
	        then 
	            sbatch ../lentFoam.sbatch
	        else 	
	            sbatch ../lentFoamNoRhoEqn.sbatch
	        fi
	        # Do not overburden SLURM workload manager
	        sleep 1
	    elif [ "$RUNNER" == "foamJob" ]; 
	    then 
	        if [ -z "$NO_RHO_EQN" ]; 
	        then 
	            foamJob lentFoam 	
	        else 
	            foamJob lentFoam -no-density-equation
	        fi
	    elif [ "$RUNNER" == "commandLine" ]; 
	    then 
	        if [ -z "$NO_RHO_EQN" ]; 
	        then 
	            lentFoam 
	        else 
	            lentFoam -no-density-equation
	        fi
	    fi
	    cd .. 
	done
    fi
else 
    echo 'Use following options: "slurm" (cluster), "foamJob" (background), \
	  "commandLine" (interactive), "none" (only case generation).'
    exit 1
fi
