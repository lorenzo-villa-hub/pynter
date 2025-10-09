#!/bin/sh

#SBATCH --account=p0025190
#SBATCH --error=err.txt
#SBATCH --job-name=SO_bulk_s2
#SBATCH --mem-per-cpu=2500
#SBATCH --ntasks=96
#SBATCH --output=out.txt
#SBATCH --time=24:00:00

srun /home/la303207/vasp_std
