#!/bin/bash

if [ -z "$LENTSCRIPTSSOURCED" ]
then
    echo "Error: LENT scripts not sourced. Source LENT/cases/scripts/bashrc before using this script."
    exit 1
fi

# Prepare and initialize the study directories
lent_prepare_study_variants.py siliconoilAirDropletTranslation.parameter -m block

# Submit a job for each variant
for dir in siliconoilAirDropletTranslation.parameter_*; do
    cd $dir
    echo "Submitting $(pwd)"
    sbatch ../lentFoam.sbatch
    cd ..
    # Do not overburden SLURM workload manager
    sleep 1
done
