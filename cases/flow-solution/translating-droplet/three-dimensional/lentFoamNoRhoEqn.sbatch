#!/bin/bash
#SBATCH	-A special00005
#SBATCH -J lentFoam
#SBATCH -n 1 
#SBATCH --mem-per-cpu=6000 
#SBATCH -t 15:00:00     

# -------------------------------

srun lentFoam -no-density-equation 
