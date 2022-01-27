#!/bin/bash

# Source LENT's python modules for parameter variations
source ../../../../cases/scripts/bashrc

# Prepare and initialize the study directories
lent_prepare_study_variants.py mercuryAirDropletTranslation.parameter -m block

# Submit a job for each variant
for dir in mercuryAirDropletTranslation.parameter_*; do
    cd $dir
    echo "Submitting $(pwd)"
    sbatch ../lentFoam.sbatch
    cd ..
    # Do not overburden SLURM workload manager
    sleep 1
done
