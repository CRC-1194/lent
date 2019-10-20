#!/bin/bash
#SBATCH -J lentFoam
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tolle@mma.tu-darmstadt.de
#SBATCH -n 1 
#SBATCH --mem-per-cpu=8000   
#SBATCH -t 24:00:00     

# -------------------------------

srun lentFoam
