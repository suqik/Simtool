#!/bin/bash
#SBATCH -J GAL_VOID_L1000_N1024_validation
#SBATCH -p kshcnormal
#SBATCH -N 4
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/run_GAL_VOID_L1000_N1024_validational_void.out
#SBATCH --error=job_outputs/GAL_VOID_L1000_N1024_validation.err

module purge
module load compiler/gnu/9.3.0 mathlib/gsl/2.7/intel compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export LD_LIBRARY_PATH=/public/home/suchen/applications/CGAL-5.6.1/lib:$LD_LIBRARY_PATH

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

mpirun -np 100 $PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/run_gal_void.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/val_input.ini
