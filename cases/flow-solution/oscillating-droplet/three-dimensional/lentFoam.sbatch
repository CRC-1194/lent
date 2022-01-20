#!/bin/bash
#SBATCH -p test7d
#SBATCH -J lentFoam 
#SBATCH -A special00005
#SBATCH -n 1 
#SBATCH --mem-per-cpu=6000 
#SBATCH -t 06-23:59:00     

# -------------------------------

srun lentFoam 
