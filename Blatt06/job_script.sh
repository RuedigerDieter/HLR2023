#!/bin/bash
#SBATCH -N 1                      # Anzahl der Knoten
#SBATCH -n 5                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o timempi.out
#SBATCH -p west
srun -N 1 -n 5 -p west timempi.sh

echo "fertig" >> job.out
