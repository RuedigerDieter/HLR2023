#!/bin/bash
#SBATCH -N 3                      # Anzahl der Knoten
#SBATCH -n 8                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o output.out
#SBATCH -p west
#SBATCH --export=ALL
./run.sh > para.txt
./run1.sh > para1.txt
./run2.sh > para2.txt
