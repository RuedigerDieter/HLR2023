#!/bin/bash
#SBATCH -N 1                      # Anzahl der Knoten
#SBATCH -n 1                     # Anzahl der Prozesse (5 Knoten * 5 Prozesse)
#SBATCH -o output.out
#SBATCH -p west
#SBATCH --export=ALL

. /opt/spack/20220821/share/spack/setup-env.sh
spack load scorep

srun bash -c 'SCOREP_ENABLE_TRACING=true ./cicle.sh'
