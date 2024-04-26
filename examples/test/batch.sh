#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 72
#SBATCH --time 04:00:00

export OMP_NUM_THREADS=72
ulimit -s unlimited
../../build/shy sod_LF_3materials_mix.json
