#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 1                     
#SBATCH -o job.out
#SBATCH -p west
srun west partdiff.sh

echo "fertig" >> job.out
