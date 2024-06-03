#!/bin/bash
start=`date +%s`

### begin of the job
HOME=/home/suqikuai777
WDIR=${HOME}/Program/SBI_Void/simulator/fastpm
PYTHON=python
MPIRUN=/usr/local/bin/mpirun
MKCFG=${WDIR}/Simtool/make-cfgs.py
FCFC=${HOME}/applications/Measurements/FCFC/FCFC_2PT
DDIR=${HOME}/Data/BOSS_dr12

cfgpath=${WDIR}/cfgs
if [ ! -d $cfgpath ]; then
    mkdir $cfgpath
fi

### FCFC CONFIGURATION SETUP
export OMP_NUM_THREADS=1

fcfc_type=fcfc_sky

binning=2 # 0 for 2pcf, 1 for multipoles, 2 for proj_cf
smin=0
smax=75
size=3

if [ $binning -eq 1 ]; then
	mu="[0,2,4]"
	mubin=100
else
	mu=None
	mubin=None
fi

if [ $binning -eq 2 ]; then
	pj=T
	pimin=20
	pimax=200
else
	pj=F
	pimin=None
	pimax=None
fi

selection=None

for ((i=1; i<2; i++)) do
	# DM input file
    if [ $i -lt 10 ];
	then
		idx=000${i}
	elif [ $i -lt 100 ];
	then
		idx=00${i}
	elif [ $i -lt 1000 ];
	then 
		idx=0${i}
	elif [ $i -lt 10000 ];
	then
		idx=${i}
	fi
    dfile=${DDIR}/mocks/mock_galaxy_DR12_LOWZ_N_QPM_${idx}.rdzw 
    # Void input file
    rfile=${DDIR}/mocks/mock_random_DR12_LOWZ_N_50x1.rdzw
	if [ ! -d ${WDIR}/multipoles/ ]; then
		mkdir -p ${WDIR}/multipoles
	fi

	opfilebase=${WDIR}/multipoles/mock_LOWZ_N_QPM_${idx}
	cffile=${opfilebase}.clu
	if [ $binning -eq 1 ]; then
		mufile=${opfilebase}.multi
	else
		mufile=None
	fi
	if [ $binning -eq 2 ]; then
		pjfile=${opfilebase}.proj
	else
		pjfile=None
	fi
	cfgfile=$cfgpath/fcfc_${idx}.conf

	$PYTHON $MKCFG $fcfc_type \
	-CATALOG "[$dfile, $rfile]" -CATALOG_LABEL "[D, R]" \
	-ASCII_FORMATTER "[\"%f %f %f %f %f\", \"%f %f %f %f\"]" -POSITION "[\$1, \$2, \$3, \$1, \$2, \$3]" -WEIGHT "[\$4 * \$5, \$4]" -SELECTION $selection \
	-OMEGA_M 0.31 -OMEGA_LAMBDA 0.69 -DE_EOS_W m1 -CMVDST_ERR "1e-10" \
	-BINNING_SCHEME $binning -PAIR_COUNT "[DD, DR, RR]" -PAIR_COUNT_FILE "[pairs.dd, pairs.dr, pairs.rr]" -CF_ESTIMATOR "(DD - 2*DR + RR)/RR" -CF_OUTPUT_FILE "$cffile" \
	-SEP_BIN_MIN $smin -SEP_BIN_MAX $smax -SEP_BIN_SIZE $size \
	-MULTIPOLE $mu -MULTIPOLE_FILE $mufile -MU_BIN_NUM $mubin \
	-PROJECTED_CF $pj -PROJECTED_FILE $pjfile -PI_BIN_MIN $pimin -PI_BIN_MAX $pimax -PI_BIN_SIZE $size \
	-ofile $cfgfile

	$FCFC -c $cfgfile
	# rm pairs.dd pairs.dr 

done
### end of the job

end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
