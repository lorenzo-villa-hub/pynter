#!/bin/sh
#SBATCH -A project01596
#SBATCH --job-name=NN_Pbcm_NEB_O-vacancy_q2_NEB_4
#SBATCH --mail-user=villa@mm.tu-darmstadt.de
#SBATCH --mail-type=ALL
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --time=24:00:00
#SBATCH -p test24
#SBATCH --exclusive
#SBATCH --mem-per-cpu=3500

module purge
ml intel/2020.2 
ml intelmpi/2020.2 
ml fftw/3.3.8 
I_MPI_PMI_LIBRARY=/opt/slurm/current/lib/libpmi2.so
srun /home/lv51dypu/vasp.5.4.4_vtstcode/bin/vasp_std
