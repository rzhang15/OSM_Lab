#!/bin/bash

#SBATCH --job-name=dot_prod-openmp
#SBATCH --output=dot_prod-openmp.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sandyb

export OMP_NUM_THREADS=600

./dot_prod.exec
