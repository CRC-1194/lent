#!/bin/bash
#SBATCH -J lentFoam

#SBATCH -n 1 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 13:00:00     

# -------------------------------

srun lentFoam 
