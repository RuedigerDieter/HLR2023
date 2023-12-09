#!/bin/sh
SCOREP_ENABLE_TRACING=true mpirun -n 4 ./circle 10
