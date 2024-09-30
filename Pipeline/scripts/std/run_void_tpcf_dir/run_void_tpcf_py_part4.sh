#!/bin/bash
#SBATCH -J run_void_tpcf_L1000_N1024_1000cosmo_part4
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/run_void_tpcf_dir/L1000_N1024_1000cosmo_part4.out
#SBATCH --error=job_outputs/run_void_tpcf_dir/L1000_N1024_1000cosmo_part4.err

start=`date +%s`
module purge
module load compiler/gnu/9.3.0 mathlib/gsl/2.7/intel compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export LD_LIBRARY_PATH=/public/home/suchen/applications/CGAL-5.6.1/lib:$LD_LIBRARY_PATH

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src/run_tpcf.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/input.ini -s 200 -e 250
end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
