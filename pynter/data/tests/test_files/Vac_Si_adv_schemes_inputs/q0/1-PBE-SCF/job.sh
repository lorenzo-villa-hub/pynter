#!/bin/sh
#SBATCH --account=project0000
#SBATCH --error=err.%j
#SBATCH --job-name=Si_adv_schemes_Vac_Si_q0_PBE-rel_1
#SBATCH --mail-user=test@pynter.com
#SBATCH --mem-per-cpu=3500
#SBATCH --nodes=1
#SBATCH --ntasks=384
#SBATCH --output=out.%j
#SBATCH --partition=deflt
#SBATCH --time=00:30:00

module purge
ml intel/2019.2 
ml intel/2019.3 
ml intelmpi/2019.3 
ml fftw/3.3.8 

srun /home/vasp-5-3-3

automation_vasp.py --contcar --chgcar --wavecar --check-kpoints --error-check
