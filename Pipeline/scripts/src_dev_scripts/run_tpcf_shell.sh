#!/bin/bash
ITE_NUM=10
for ((i=0; i<$ITE_NUM; i++)); do

BASENAME=run_tpcf_val
if [[ ! -d ${BASENAME}_dir ]]; then
    mkdir ${BASENAME}_dir
fi
if [[ ! -d job_outputs/${BASENAME}_dir ]]; then
    mkdir -p job_outputs/${BASENAME}_dir
fi

COSMO_START=$(( $i*10 ))
COSMO_END=$(( ($i+1)*10 ))

FNAME=${BASENAME}_dir/${BASENAME}_py_part$i.s
cat <<EOF >${FNAME}
#!/bin/bash
#SBATCH -J ${BASENAME}_L1000_N1024_rlzs_part$i
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/${BASENAME}_dir/L1000_N1024_rlzs_part$i.out
#SBATCH --error=job_outputs/${BASENAME}_dir/L1000_N1024_rlzs_part$i.err

start=\`date +%s\`
module purge
module load compiler/gnu/9.3.0 mathlib/gsl/2.7/intel compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export LD_LIBRARY_PATH=/public/home/suchen/applications/CGAL-5.6.1/lib:\$LD_LIBRARY_PATH

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

\$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/run_tpcf.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/val_input.ini -s $COSMO_START -e $COSMO_END
end=\`date +%s\`
dif=\$[ end - start ]
echo running time: \$dif sec
EOF

sbatch <${FNAME}
done
