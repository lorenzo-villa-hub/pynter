#!/bin/sh
#SBATCH -A special00003
#SBATCH --job-name=Si-test_Sub_P_on_Si-Sub_B_on_Si_q0_PBE-rel_1
#SBATCH --mail-user=villa@mm.tu-darmstadt.de
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --mem-per-cpu=3500

module purge
ml intel/2020.4 
ml intelmpi/2020.4 
ml fftw/3.3.10 

srun /home/groups/da_mm/codes/vasp/build/vasp.5.4.4/bin/vasp_std

automation_vasp.py --contcar --chgcar --wavecar --error-check --check-kpoints
