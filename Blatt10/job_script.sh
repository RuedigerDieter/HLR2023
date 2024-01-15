#!/bin/bash
#SBATCH -N 3                      # Anzahl der Knoten
#SBATCH -n 8                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o output.out
#SBATCH -p west
#SBATCH --export=ALL
mpirun ./run > para.txt

