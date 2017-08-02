#!/bin/bash

#SBATCH --job-name=BS_openmp

#SBATCH --output=BS_openmp.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --partition=sandyb

export OMP_NUM_THREADS=1000

./BS.exec
