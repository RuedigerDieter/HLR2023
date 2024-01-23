#!/bin/sh

sbatch --wait  GS_3x2.job
sbatch --wait GS_5x4.job
sbatch --wait JA_3x2.job
sbatch --wait JA_5x4.job
