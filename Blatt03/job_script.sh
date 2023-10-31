#!/bin/bash
#SBATCH -N 5                      # Anzahl der Knoten
#SBATCH -n 25                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)

srun -N 5 -n 25 -p west timescript.sh >> timescript.out

echo "fertig" >> timescript.out
