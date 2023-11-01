#!/bin/bash
#SBATCH -N 5                      
#SBATCH -n 25                     
srun -N 5 -n 25 -p west timescript.sh >> timescript.out

echo "fertig" >> timescript.out
