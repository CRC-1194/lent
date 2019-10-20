#!/bin/bash
#SBATCH -J createFront
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tolle@mma.tu-darmstadt.de
#SBATCH -n 1       
#SBATCH --mem-per-cpu=6000   
#SBATCH -t 00:10:00     

# -------------------------------

srun lentCreateFront
