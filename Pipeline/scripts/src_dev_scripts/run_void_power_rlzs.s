#!/bin/bash
#SBATCH -J VOID_POWER_L1000_N1024
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --output=job_outputs/VOID_POWER_L1000_N1024.out
#SBATCH --error=job_outputs/VOID_POWER_L1000_N1024.err

start=`date +%s`
module purge
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export OMP_NUM_THREADS=32

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/run_power_rlzs.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/input.ini -s 1 -e 10
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
