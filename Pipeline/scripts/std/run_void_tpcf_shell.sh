#!/bin/bash
IND_SIM_NUM=50
TOTAL_SIM=1000
ITE_NUM=$(( $TOTAL_SIM/$IND_SIM_NUM ))
for ((i=1; i<$ITE_NUM; i++)); do

COSMO_START=$(( $i*$IND_SIM_NUM ))
COSMO_END=$(( ($i+1)*$IND_SIM_NUM ))
BASENAME=run_void_tpcf
if [[ ! -d ${BASENAME}_dir ]]; then
    mkdir ${BASENAME}_dir
fi
if [[ ! -d job_outputs/${BASENAME}_dir ]]; then
    mkdir -p job_outputs/${BASENAME}_dir
fi

FNAME=${BASENAME}_dir/${BASENAME}_py_part$i.sh
cat <<EOF >${FNAME}
#!/bin/bash
#SBATCH -J ${BASENAME}_L1000_N1024_1000cosmo_part$i
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=job_outputs/${BASENAME}_dir/L1000_N1024_1000cosmo_part$i.out
#SBATCH --error=job_outputs/${BASENAME}_dir/L1000_N1024_1000cosmo_part$i.err

start=\`date +%s\`
module purge
module load compiler/gnu/9.3.0 mathlib/gsl/2.7/intel compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export LD_LIBRARY_PATH=/public/home/suchen/applications/CGAL-5.6.1/lib:\$LD_LIBRARY_PATH

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

\$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src/run_tpcf.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/input.ini -s $(( $i*$IND_SIM_NUM )) -e $(( ($i+1)*$IND_SIM_NUM ))
end=\`date +%s\`
dif=\$[ end - start ]
echo running time: \$dif sec
EOF

sbatch <${FNAME}
done