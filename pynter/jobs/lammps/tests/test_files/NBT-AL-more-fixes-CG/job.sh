#!/bin/sh

#SBATCH --account=special00003
#SBATCH --error=err.txt
#SBATCH --job-name=NBT_Bi4Na4Ti8O24_L0001
#SBATCH --mem-per-cpu=2400M
#SBATCH --ntasks=64
#SBATCH --output=out.txt
#SBATCH --time=01:00:00

module purge
module load lammps_ml
srun lmp -in input.in