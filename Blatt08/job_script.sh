#!/bin/bash
#SBATCH -N 4                      # Anzahl der Knoten
#SBATCH -n 5                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o timempi.out
#SBATCH -p west
srun -N 4 -n 5 -p west partdiff.sh

echo "fertig" >> job.out
