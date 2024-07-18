#!/bin/bash
#SBATCH -J FPMsimu_test
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH --output=FPMsimu_test.out
#SBATCH --error=FPMsimu_test.err

module purge
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0 mathlib/gsl/2.7/intel

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src/run_fastpm.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/input.ini
