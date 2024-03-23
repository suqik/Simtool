#!/bin/bash
start=`date +%s`

### Note: `dir` means dictionaries, `path` means files.
wdir=/home/suqikuai777/Program/SBI_Void/simulator/fastpm # working dictionary
catdir=${wdir}/catalog/match_Nbody # catalog dictionary

### cosmology
Hubble=0.6727
OmegaM=0.3156
OmegaL=0.6844
Omegab=0.0491
NS=1.0
SIGMA8=0.831

### Flags
PKLIN_FLAG=0
FPMPART_FLAG=1
CONVERT_FLAG=1
LCPART_FLAG=0
RAYPART_FLAG=0

##### FastPM Part #####

### calculate linear power spectrum
PYTHON=python
MKPKLIN=/home/suqikuai777/Program/SBI_Void/simulator/fastpm/Simtool/make-pklin.py # execute file of making linear pk
powerpath=${wdir}/cfgs/fastpm_cfgs/match_Nbody/powerspec_nbody.txt

if [ $PKLIN_FLAG -eq 1 ]; then
    $PYTHON $MKPKLIN -hubble $Hubble -OmegaM $OmegaM -OmegaL $OmegaL -sigma8 $SIGMA8 -ofile $powerpath
    echo -e "Calculate linear power spectrum [\e[32mDONE\e[0m]"
fi

### slice distance arrays
nslice=23
nbox=5
rlow=(40 440 840 1240 1640)
rup=(361 761 1161 1561 1840)
dr=80


MKZSLC=${wdir}/Simtool/make-redshifts.py # execute file of calculating redshifts
MKCFG=${wdir}/Simtool/make-cfgs.py # execute file of generating configuration
FASTPM=/home/suqikuai777/applications/Simulations/fastpm/src/fastpm # execute file of FastPM simulation
DIVE=/home/suqikuai777/applications/special/DIVE/DIVE

### simulation general settings
nc=256 # simulation particle number
boxsize=400 # box length in Mpc/h
for (( irlz=0; irlz<10; irlz++ )) do
    ### FastPM params

    iseed=$(($irlz+1))00 # initial random seed

    cfgdir=${wdir}/cfgs/fastpm_cfgs/match_Nbody/rlz$irlz/
    seedinfo=${cfgdir}/seed_${iseed}_$(($iseed+$nbox-1))

    if [ ! -d $cfgdir ]; then
        mkdir $cfgdir
    fi
    if [ ! -f $seedinfo ]; then
        touch $seedinfo
    fi

    if [ $FPMPART_FLAG -eq 1 ]; then
        for ((ibox=0; ibox<$nbox; ibox++)) do
            # calculate redshifts of slices
            ZSLC=$($PYTHON $MKZSLC -a "(${rlow[$ibox]} ${rup[$ibox]} $dr)" -hubble $Hubble -Om $OmegaM -OL $OmegaL)

            # generate configure file
            seed=$(($iseed+$ibox))
            
            cfgpath=${cfgdir}/box$ibox.lua
            snapdir=${catdir}/rlz$irlz/box$ibox

            $PYTHON $MKCFG fastpm -nc $nc boxsize $boxsize -output_redshifts "$ZSLC" -Omega_m $OmegaM -hubble $Hubble -read_powerspectrum $powerpath -random_seed $seed -write_snapshot ${snapshot}/a -write_rfof ${snapdir}/a -ofile $cfgpath
            echo -e "Generate FastPM configuration file [\e[32mDONE\e[0m]"

            # run FastPM
            export OMP_NUM_THREADS=2

            mpirun -np 16 $FASTPM $cfgpath
            echo -e "Run FastPM [\e[32mDONE\e[0m]"

            # convert to gadget format
            if [ $CONVERT_FLAG -eq 1 ]; then
                I=0

                files=$(ls ${snapdir}/ | grep a_) 
                CONVERT=${wdir}/Simtool/convert-to-txt.py # execute file of converting to gadget format
                for file in $files
                do
                    if [ ! -d ${snapdir}/${file}/Void ]; then
                        mkdir ${snapdir}/${file}/Void
                    fi
                    $PYTHON $CONVERT ${snapdir}/${file}/ halo ${snapdir}/${file}/halo_tmp.txt
                    $DIVE -i ${snapdir}/${file}/halo_tmp.txt -o ${snapdir}/${file}/Void.txt -u $boxsize
                    rm ${snapdir}/${file}/halo_tmp.txt
                done

                CONVERT=/home/suqikuai777/applications/Simulations/fastpm/python/convert-to-gadget-1.py # execute file of converting to gadget format
                if [ ! -d ${snapdir}/gadget/ ]; then
                    mkdir ${snapdir}/gadget
                fi
                for file in $files
                do
                    if [ $I -lt 10 ]; 
                    then
                        $PYTHON $CONVERT ${snapdir}/${file} ${snapdir}/gadget/snapshot_00$I --nperfile $((${nc}*${nc}*${nc}))
                    elif [ $I -lt 100 ]; 
                    then
                        $PYTHON $CONVERT ${snapdir}/${file} ${snapdir}/gadget/snapshot_0$I --nperfile $((${nc}*${nc}*${nc}))
                    elif [ $I -lt 1000 ]; 
                    then
                        $PYTHON $CONVERT ${snapdir}/${file} ${snapdir}/gadget/snapshot_$I --nperfile $((${nc}*${nc}*${nc}))
                    fi
                    I=$(($I + 1))
                done
                echo -e "Convert to GADGET format [\e[32mDONE\e[0m]"
            fi
        done
    fi

    ##### Raytracing Part #####
    cfgdir=${wdir}/cfgs/raytrace_cfgs
    GENERATOR=/home/suqikuai777/applications/Simulations/gluons_parallel/preproc/generator
    GLUONS=/home/suqikuai777/applications/Simulations/gluons_parallel/gluons

    ### generate LC configuration file
    cfgpath=${cfgdir}/preproc.param
    lcdir=${catdir}/rlz$irlz/LC
    if [ ! -d $lcdir ]; then
        mkdir $lcdir
    fi

    if [ $LCPART_FLAG -eq 1 ]; then
        $PYTHON $MKCFG preproc -SIMPATH $catdir -SIMBASE rlz${irlz}/box -SNAPBASE gadget/snapshot -SIMNUM $nbox -NPART $nc -FILENUM 1 -SLICEPATH ${lcdir} -HALOPATH ${catdir}/HALO -LIST 1 -b $nbox -s $nslice -w $dr -ofile $cfgpath
        echo -e "Generate LC configuration file [\e[32mDONE\e[0m]"

        ### run lightcone
        $GENERATOR $cfgpath
        echo -e "Run construct LightCone [\e[32mDONE\e[0m]"
    fi

    ### generate Raytracing configuration file
    cfgpath=${cfgdir}/raytrace.param
    lensdir=${catdir}/rlz$irlz/lens
    if [ ! -d $lensdir ]; then
        mkdir $lensdir
    fi

    ### Raytracing params
    rng=36000.0
    nray=1024
    cra=0.5
    cdec=0.5
    zsrc=0.75
    fov=$(echo "scale=2;${rng}/3600" | bc)

    if [ $RAYPART_FLAG -eq 1 ]; then
        $PYTHON $MKCFG gluons -OMEGAX $OmegaL -OMEGAM $OmegaM -HUBBLE $Hubble -NUMSLC $nslice -PTHSLC ${lcdir}/slice -OUTPUT ${lensdir}/fov${fov} -RNGDEC $rng -RNGRA $rng -NUMDEC $nray -NUMRA $nray -CTRDEC $cdec -CTRRA $cra -ZSOURC $zsrc -ofile $cfgpath
        echo -e "Generate Raytracing configuration file [\e[32mDONE\e[0m]"

        ### Run Raytracing 
        $GLUONS $cfgpath
        echo -e "Run Raytracing [\e[32mDONE\e[0m]"
    fi

    ### delete gadget files to save disk space
    for ((j=0; j<5; j++)) do
        if [ -d ${catdir}/rlz$irlz/box$j/gadget/ ]; then
            rm -r ${catdir}/rlz$irlz/box$j/gadget/
        fi
    done
done

end=`date +%s`
dif=$[ end - start ]
echo running time: $dif sec
