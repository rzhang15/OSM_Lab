#!/bin/bash

#SBATCH --job-name=normalize_vec_openmp

#SBATCH --output=normalize_vec_openmp.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --partition=sandyb

export OMP_NUM_THREADS=500

./normalize_vec.exec
