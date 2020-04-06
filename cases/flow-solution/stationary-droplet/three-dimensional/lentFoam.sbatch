#!/bin/bash
#SBATCH -J lentFoam
#SBATCH -n 1 
#SBATCH --mem-per-cpu=1000 
#SBATCH -t 10:00:00     

# -------------------------------

srun lentFoam 
