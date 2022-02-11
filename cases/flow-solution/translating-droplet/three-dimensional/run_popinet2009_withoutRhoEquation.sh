#!/bin/bash

RUNNER=$1

if [[ ("$RUNNER" == "slurm") || ("$RUNNER" == "foamJob") || ("$RUNNER" == "commandLine") ]]; 
then 
    # Source LENT's python modules for parameter variations
    source ../../../../cases/scripts/bashrc

    # Prepare and initialize the study directories
    lent_prepare_study_variants.py popinet2009.parameter -p without-rhoEquation -m block

    # Run each variant 
    for dir in without-rhoEquation-popinet2009*; do
        cd $dir
        if [ "$RUNNER" == "slurm" ]; 
        then 
            echo "Submitting $(pwd)"
            sbatch ../lentFoamNoRhoEqn.sbatch
            # Do not overburden SLURM workload manager
            sleep 1
        elif [ "$RUNNER" = "foamJob" ]; 
        then 
            foamJob lentFoam -no-density-equation 
        elif [ "$RUNNER" = "commandLine" ]; 
        then 
            lentFoam -no-density-equation 
        fi
        cd ..
    done
else
    echo 'Use following options: "slurm" (cluster), "foamJob" (background), or "commandLine".'
    exit 1
fi
