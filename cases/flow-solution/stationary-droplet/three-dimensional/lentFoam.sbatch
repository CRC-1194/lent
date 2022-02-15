#!/bin/bash
#SBATCH -J lentFoam 
#SBATCH -A special00004
#SBATCH -n 1 
#SBATCH --mem-per-cpu=3000 
#SBATCH -t 23:00:00     

# -------------------------------

srun lentFoam 
# -no-density-equation
