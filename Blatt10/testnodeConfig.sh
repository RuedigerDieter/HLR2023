#!/bin/bash
#SBATCH -o output.out
#SBATCH -p west
#SBATCH --export=ALL

for i in 1 2 3 4
do
    for j in 1 2 3 4
    do
        k = $((i*j))
        echo "Nodes: $i, Tasks per Node: $j, Tasks: $k"
        mpiexec -p west --nodes $i --ntasks-per-node $j --ntasks $k ./run.sh | grep "Norm"
    done
done