import os
import numpy as np
from convert import Convert
from sham import load_halo_head_data, SHAM_sigma_model
##===================== Driver Functions ========================##

def fastpm_driver(
        fpm_cfgpath, 
        FastPM_exec = "/public/home/suchen/applications/fastpm_intel/src/fastpm", 
        nCPUs = 32
        ):
    
    cmd = "mpirun -np "+f"{nCPUs} "+FastPM_exec+" "+fpm_cfgpath
    print(cmd, flush=True)
    os.system(cmd)

def rockstar_driver(
        snappath,
        nfile, 
        precision,
        rstar_cfgpath,
        gadgetpath,
        halopath, 
        tmp_script_name,
        Rockstar_exec = "/public/home/suchen/applications/rockstar/rockstar",
        FindPAR_exec = None,
        boxsize = None,
):
    
    Convert(snappath, gadgetpath, nfile, precision)

    ### I don't know why but this way works ...
    f = open(f"{tmp_script_name}", "w+")
    f.write("#!/bin/bash\n")
    f.write(f"RSTAR={Rockstar_exec}\n")
    if FindPAR_exec is not None and boxsize is not None:
        f.write(f"FINDPAR={FindPAR_exec}\n")
        f.write(f"BOXSIZE={boxsize}\n")
    f.write(f"CFG={rstar_cfgpath}\n")
    f.write(f"IDIR={gadgetpath}\n")
    f.write(f"ODIR={halopath}\n")
    f.write("$RSTAR -c $CFG &\n")
    if FindPAR_exec is not None:
        f.write("$FINDPAR ${ODIR}/out_0.list ${BOXSIZE} >${ODIR}/out_0_wsub.list\n")
    f.write("export ODIR\n")
    f.write("perl -e \'sleep 1 while (!(-e \"$ENV{ODIR}/auto-rockstar.cfg\"))\'\n")
    f.write("$RSTAR -c ${ODIR}/auto-rockstar.cfg\n")
    f.write("rm -r ${ODIR}/halos* ${ODIR}/*.cfg ${ODIR}/profiling/\n")
    f.write("rm ${IDIR}*\n")

    f.close()
    os.system(f"bash {tmp_script_name}")
    os.system(f"rm {tmp_script_name}")

def gal_void_driver(
        sham_param_names, 
        sham_param_vals, 
        halopath,
        halofinder:str,
        feature, 
        ref_num_den, 
        seed,
        zspace,
        galpath,
        voidpath,
        DIVE_exec
):
    
    halo = load_halo_head_data(halopath, halofinder, feature, zspace)
    pos = np.c_[halo["x"], halo["y"], halo["z"]]
    boxsize = halo.meta['boxsize']

    sham_param_names = np.atleast_1d(sham_param_names)
    sham_param_vals = np.atleast_1d(sham_param_vals)

    if len(sham_param_names) == 1 and sham_param_names[0] == "sigma":
        gsamples = SHAM_sigma_model(sham_param_vals[0], pos, feature, boxsize, ref_num_den, seed)
        np.savetxt(galpath, gsamples[:,:3], fmt="%.3f %.3f %.3f")
    
    os.system(f"{DIVE_exec} -i {galpath} -o {voidpath} -u {boxsize}")
    