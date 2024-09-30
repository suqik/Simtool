#!/bin/bash
#SBATCH -J AUX_TASK
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --output=AUX_TASK.out
#SBATCH --error=AUX_TASK.err

start=`date +%s`
module purge
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0 mathlib/gsl/2.7/intel

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/aux/show_pk.py
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec