#!/bin/bash
start=`date +%s`

### begin of the job
HOME=/home/suqikuai777
WDIR=${HOME}/Program/SBI_Void/simulator/fastpm
PYTHON=python
MPIRUN=/usr/local/bin/mpirun
MKCFG=${WDIR}/Simtool/make-cfgs.py
FCFC=${HOME}/applications/Measurements/FCFC_intel/FCFC_2PT
DDIR=${HOME}/Data/BOSS_dr12

cfgpath=${WDIR}/cfgs
if [ ! -d $cfgpath ]; then
    mkdir $cfgpath
fi

### FCFC CONFIGURATION SETUP

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
	pimin=0
	pimax=75
else
	pj=F
	pimin=None
	pimax=None
fi

basename=DR12v5_LOWZ_North

dfile=${DDIR}/galaxy_${basename}.fits
rfile=${DDIR}/random0_${basename}.fits

resultdir=${WDIR}/results
if [ ! -d ${resultdir}/ ]; then
	mkdir -p ${resultdir}
fi

opfilebase=${resultdir}/galaxy_${basename}
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
cfgfile=$cfgpath/fcfc_data.conf

$PYTHON $MKCFG $fcfc_type \
-CATALOG "[$dfile, $rfile]" -CATALOG_TYPE "[1, 1]" -CATALOG_LABEL "[D, R]" \
-POSITION "[\${RA}, \${DEC}, \${Z}, \${RA}, \${DEC}, \${Z}]" -WEIGHT "[\${WEIGHT_SYSTOT} * \${WEIGHT_NOZ} * \${WEIGHT_CP} * \${WEIGHT_FKP}, \${WEIGHT_FKP}]" -SELECTION "[\${Z}<0.43 && \${Z}>0.2, \${Z}<0.43 && \${Z}>0.2]" \
-OMEGA_M 0.31 -OMEGA_LAMBDA 0.69 -DE_EOS_W m1 -CMVDST_ERR "1e-10" \
-BINNING_SCHEME $binning -PAIR_COUNT "[DD, DR, RR]" -PAIR_COUNT_FILE "[${opfilebase}.dd, ${opfilebase}.dr, ${opfilebase}.rr]" -CF_ESTIMATOR "(DD - 2*DR + RR)/RR" -CF_OUTPUT_FILE "$cffile" \
-SEP_BIN_MIN $smin -SEP_BIN_MAX $smax -SEP_BIN_SIZE $size \
-MULTIPOLE $mu -MULTIPOLE_FILE $mufile -MU_BIN_NUM $mubin \
-PROJECTED_CF $pj -PROJECTED_FILE $pjfile -PI_BIN_MIN $pimin -PI_BIN_MAX $pimax -PI_BIN_SIZE $size \
-ofile $cfgfile

$FCFC -c $cfgfile
# rm pairs.dd pairs.dr 

### end of the job

end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
