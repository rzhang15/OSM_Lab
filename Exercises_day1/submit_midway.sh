#!/bin/bash

#SBATCH --job-name=compute_pi

#SBATCH --output=compute_pi.out

#SBATCH --mail-type=BEGIN,END,DAIL

#SBATCH --ntasks=1

#SBATCH --partition=sandyb

./compute_pi.cpp.exec
