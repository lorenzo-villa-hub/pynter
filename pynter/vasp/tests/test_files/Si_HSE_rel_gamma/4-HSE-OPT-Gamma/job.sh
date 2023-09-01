#!/bin/sh
#SBATCH --account=special00003
#SBATCH --error=err.%j
#SBATCH --job-name=Si_schemes_HSE-rel-gamma-ext_4
#SBATCH --mail-user=villa@mm.tu-darmstadt.de
#SBATCH --mem-per-cpu=3500
#SBATCH --nodes=1
#SBATCH --ntasks=384
#SBATCH --output=out.%j
#SBATCH --partition=deflt
#SBATCH --time=24:00:00
#SBATCH --array=1-7%1

module purge
ml intel/2019.2 
ml intel/2019.3 
ml intelmpi/2019.3 
ml fftw/3.3.8 

if [ ! -f POSCAR_initial ] ; then
    cp POSCAR POSCAR_initial
fi
if [ -f CONTCAR ]
then
    cp CONTCAR POSCAR
fi

srun /home/vasp-5-3-3

pynter analysis vasprun --convergence > convergence.txt
if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then
    automation_vasp.py --contcar --chgcar --wavecar --check-kpoints --error-check
    scancel ${SLURM_ARRAY_JOB_ID}_*
fi
