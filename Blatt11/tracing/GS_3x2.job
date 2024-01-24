#!/bin/bash
#SBATCH -N 2                      # Anzahl der Knoten
#SBATCH --ntasks=3                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o GS_3x2.out
#SBATCH -p west
#SBATCH --export=ALL

. /opt/spack/20220821/share/spack/setup-env.sh
spack load scorep

srun bash -c 'SCOREP_ENABLE_TRACING=true mpiexec -n 3 ./partdiff 1 1 64 2 2 20'
