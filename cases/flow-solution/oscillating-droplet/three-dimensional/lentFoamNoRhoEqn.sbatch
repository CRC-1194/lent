#!/bin/bash
#SBATCH -p test24
#SBATCH -J lentNoRho 
#SBATCH -A special00005
#SBATCH -n 1 
#SBATCH --mem-per-cpu=6000 
#SBATCH -t 23:00:00     

# -------------------------------

srun lentFoam -no-density-equation
