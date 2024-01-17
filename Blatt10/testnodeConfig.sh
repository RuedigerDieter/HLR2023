#!/bin/bash
#SBATCH -o output.out
#SBATCH -p west
#SBATCH --export=ALL

for i in 1 2 3 4
do
    for j in 1 2 3 4
    do
        echo "Nodes: $i, Tasks per Node: $j"
        srun -nodes $i --ntasks-per-node $j ./run.sh | grep "Norm"
    done
done