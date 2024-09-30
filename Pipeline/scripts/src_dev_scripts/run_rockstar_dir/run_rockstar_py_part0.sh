#!/bin/bash
#SBATCH -J Rockstar_L1000_N1024_validation_part0
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --output=job_outputs/run_rockstar_dir/Rockstar_L1000_N1024_validation_part0.out
#SBATCH --error=job_outputs/run_rockstar_dir/Rockstar_L1000_N1024_validation_part0.err

start=`date +%s`
module purge
module load compiler/devtoolset/7.3.1 compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/run_rockstar.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/val_input.ini -s 10 -e 11
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
