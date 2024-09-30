#!/bin/bash
#SBATCH -J FPM_L1000_N1024_1000cosmo
#SBATCH -p kshcnormal
#SBATCH -N 32
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/FPM_L1000_N1024_1000cosmo.out
#SBATCH --error=job_outputs/FPM_L1000_N1024_1000cosmo.err

start=`date +%s`
module purge
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0 mathlib/gsl/2.7/intel

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src/run_fastpm.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/input.ini 1024
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec