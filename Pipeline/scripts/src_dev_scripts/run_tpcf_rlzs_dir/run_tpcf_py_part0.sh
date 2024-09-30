#!/bin/bash
#SBATCH -J run_tpcf_L1000_N1024_rlzs_part0
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/run_tpcf_dir/L1000_N1024_rlzs_part0.out
#SBATCH --error=job_outputs/run_tpcf_dir/L1000_N1024_rlzs_part0.err

start=`date +%s`
module purge
module load compiler/gnu/9.3.0 mathlib/gsl/2.7/intel compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export LD_LIBRARY_PATH=/public/home/suchen/applications/CGAL-5.6.1/lib:$LD_LIBRARY_PATH

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/run_tpcf.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/input.ini -i 0
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
