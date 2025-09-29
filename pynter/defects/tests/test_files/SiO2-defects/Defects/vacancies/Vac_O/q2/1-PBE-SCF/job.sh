#!/bin/sh

#SBATCH --account=p0025190
#SBATCH --error=err.txt
#SBATCH --job-name=SO_Vac_O_q2_PBE-rel_1
#SBATCH --mem-per-cpu=2500
#SBATCH --ntasks=96
#SBATCH --output=out.txt
#SBATCH --time=24:00:00

srun /home/la303207/vasp_std
source ~/.bashrc
pynter automation vasp --chgcar --contcar