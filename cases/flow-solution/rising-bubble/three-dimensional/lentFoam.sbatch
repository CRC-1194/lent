#!/bin/bash
#SBATCH -A special00004
#SBATCH -J lentFoam
#SBATCH -n 1 
#SBATCH --mem-per-cpu=6000 
#SBATCH -t 23:00:00     

# -------------------------------

srun lentFoam # -relative-frame
