#!/bin/bash
#SBATCH -J lentFoam
#SBATCH -n 1 
#SBATCH --mem-per-cpu=6000 
#SBATCH -t 23:00:00     

# -------------------------------

srun lentFoam 
