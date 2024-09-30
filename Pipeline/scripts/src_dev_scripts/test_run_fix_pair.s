#!/bin/bash
#SBATCH -J test_L1000_N1024_fix_pair
#SBATCH -p kshcnormal
#SBATCH -N 32
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/test_L1000_N1024_fix_pair.out
#SBATCH --error=job_outputs/test_L1000_N1024_fix_pair.err

start=`date +%s`
module purge
module load compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0 mathlib/gsl/2.7/intel
export OMP_NUM_THREADS=1                                                                                                              

mpirun -np 1024 /public/home/suchen/applications/fastpm_intel/src/fastpm /public/home/suchen/Programs/Simtool/Pipeline/cfgs/fix_pair/fiducial/fastpm/fpm_inv.lua

end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec