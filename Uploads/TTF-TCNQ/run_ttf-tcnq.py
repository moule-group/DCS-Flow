#!/bin/bash
#SBATCH -J TCNQ_1000steps
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -C knl
#SBATCH --output=out.out
#SBATCH --error=err.out
#SBATCH --open-mode=append

#run the application:                            
export OMP_NUM_THREADS=4
export KMP_AFFINITY=disabled
export PYTHONUNBUFFERED=1

eval $'cnss -T train'