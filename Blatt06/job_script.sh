#!/bin/bash
#SBATCH -N 5                      # Anzahl der Knoten
#SBATCH -n 25                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o partdiff.out
#SBATCH -p west
srun -N 3 -n 9 -p west timempi.sh

echo "fertig" >> timescript.out
