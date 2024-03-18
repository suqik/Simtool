#!/bin/bash
start=`date +%s`

### Note: `dir` means dictionaries, `path` means files.
wdir=/home/suqikuai777/Program/SBI_Void/simulator/fastpm # working dictionary
catdir=${wdir}/catalog/calibration # catalog dictionary
cfgdir=${wdir}/cfgs/fastpm_cfgs/calibration

if [ ! -d $catdir ]; then
    mkdir $catdir
fi
if [ ! -d $cfgdir ]; then
    mkdir $cfgdir
fi

### cosmology
Hubble=0.6727
OmegaM=0.3156
OmegaL=0.6844
Omegab=0.0491
NS=1.0
SIGMA8=0.831

### Flags
PKLIN_FLAG=1
FPMPART_FLAG=1
VOID_FLAG=1
CONVERT_FLAG=0
LCPART_FLAG=0
RAYPART_FLAG=0

##### FastPM Part #####

### calculate linear power spectrum
PYTHON=python
MKPKLIN=/home/suqikuai777/Program/SBI_Void/simulator/fastpm/Simtool/make-pklin.py # execute file of making linear pk
powerpath=${cfgdir}/powerspec.txt

if [ $PKLIN_FLAG -eq 1 ]; then
    $PYTHON $MKPKLIN -hubble $Hubble -OmegaM $OmegaM -Omegab $Omegab -sigma8 $SIGMA8 -ofile $powerpath
    echo -e "Calculate linear power spectrum [\e[32mDONE\e[0m]"
fi

MKZSLC=${wdir}/Simtool/make-redshifts.py # execute file of calculating redshifts
MKCFG=${wdir}/Simtool/make-cfgs.py # execute file of generating configuration
FASTPM=/home/suqikuai777/applications/Simulations/fastpm/src/fastpm # execute file of FastPM simulation
DIVE=/home/suqikuai777/applications/special/DIVE/DIVE

### simulation general settings
nc=256 # simulation particle number
boxsize=400 # box length in Mpc/h

iseed=100 # initial random seed
seedinfo=${cfgdir}/jiajun/seed_${iseed}

if [ ! -d ${cfgdir}/jiajun/ ]; then
    mkdir ${cfgdir}/jiajun/
fi
if [ ! -f $seedinfo ]; then
    touch $seedinfo
fi

if [ $FPMPART_FLAG -eq 1 ]; then
    # calculate output redshifts
    ZSLC="{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}"

    # generate configure file
    seed=$(($iseed))

    cfgpath=${cfgdir}/jiajun.lua
    snapdir=${catdir}/jiajun/

    $PYTHON $MKCFG fastpm -nc $nc boxsize $boxsize -output_redshifts "$ZSLC" -Omega_m $OmegaM -hubble $Hubble -read_powerspectrum $powerpath -random_seed $seed -write_rfof ${snapdir}/a -ofile $cfgpath
    echo -e "Generate FastPM configuration file [\e[32mDONE\e[0m]"

    # run FastPM
    export OMP_NUM_THREADS=2
    mpirun -np 16 $FASTPM $cfgpath
    echo -e "Run FastPM [\e[32mDONE\e[0m]"

    # convert to gadget format
    if [ $VOID_FLAG -eq 1 ]; then
        I=0
        files=$(ls ${snapdir}/ | grep a_) 
        CONVERT=${wdir}/Simtool/convert-to-txt.py # execute file of converting to txt
        for file in $files
        do
            $PYTHON $CONVERT ${snapdir}/${file}/ halo ${snapdir}/${file}/halo_tmp.txt
            $DIVE -i ${snapdir}/${file}/halo_tmp.txt -o ${snapdir}/${file}/Void.txt -u $boxsize
            rm ${snapdir}/${file}/halo_tmp.txt
        done
        echo -e "Run void finder [\e[32mDONE\e[0m]"
    fi
fi

end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
