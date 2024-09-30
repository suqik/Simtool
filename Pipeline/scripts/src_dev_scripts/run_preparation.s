#!/bin/bash
#SBATCH -J PRE_L1000_N1024_fixpair
#SBATCH -p kshcnormal
#SBATCH -n 1
#SBATCH --output=job_outputs/PRE_L1000_N1024_fixpair.out
#SBATCH --error=job_outputs/PRE_L1000_N1024_fixpair.err


start=`date +%s`
module purge
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0 mathlib/gsl/2.7/intel

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/preparation.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/fixpair_input.ini
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec