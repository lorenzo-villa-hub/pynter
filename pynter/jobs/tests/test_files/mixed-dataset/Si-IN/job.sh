#!/bin/sh
#SBATCH --account=special00003
#SBATCH --error=err.%j
#SBATCH --job-name=Si
#SBATCH --mem-per-cpu=3500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out.%j
#SBATCH --time=01:00:00

module purge
ml lammps_ml 
unset OMP_PROC_BIND
unset OMP_PLACES

srun lmp -in input.in
