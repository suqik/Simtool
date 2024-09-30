#!/bin/bash

for ((i=9; i<10; i++)); do

COSMO_START=$i
COSMO_END=100
# COSMO_START=$(( $i*100 ))
# COSMO_END=$(( ($i+1)*100 ))
if [[ ! -d run_rockstar_dir ]]; then
    mkdir run_rockstar_dir
fi
if [[ ! -d job_outputs/run_rockstar_dir ]]; then
    mkdir job_outputs/run_rockstar_dir
fi

FNAME=run_rockstar_dir/run_rockstar_py_part$i.sh
cat <<EOF >${FNAME}
#!/bin/bash
#SBATCH -J Rockstar_L1000_N1024_validation_part$i
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --output=job_outputs/run_rockstar_dir/Rockstar_L1000_N1024_validation_part$i.out
#SBATCH --error=job_outputs/run_rockstar_dir/Rockstar_L1000_N1024_validation_part$i.err

start=\`date +%s\`
module purge
module load compiler/devtoolset/7.3.1 compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0

PYTHON=/public/home/suchen/miniconda3/envs/nbodykit-env/bin/python

\$PYTHON /public/home/suchen/Programs/Simtool/Pipeline/src_dev/run_rockstar.py /public/home/suchen/Programs/Simtool/Pipeline/cfgs/val_input.ini -s $COSMO_START -e $COSMO_END
end=\`date +%s\`
dif=\$[ end - start ]
echo running time: \$dif sec
EOF

sbatch <${FNAME}
done