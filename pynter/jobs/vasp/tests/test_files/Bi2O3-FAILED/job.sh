#!/bin/sh

#SBATCH --account=p0024766
#SBATCH --error=err.txt
#SBATCH --job-name=Bi32O48_DFT_062
#SBATCH --mem-per-cpu=2400M
#SBATCH --ntasks=64
#SBATCH --output=out.txt
#SBATCH --time=24:00:00

srun /home/la303207/vasp_std