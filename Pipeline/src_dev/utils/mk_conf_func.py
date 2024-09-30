import os
import configparser
import numpy as np
from nbodykit.cosmology import LinearPower, Cosmology
from nbodykit.lab import LinearMesh
from .cfg_params import fastpm_default, rockstar_default

def mk_ini_Pk(hubble, Om0, Ob0, ns, sigma8, output):
    MYcosmo = Cosmology(m_ncdm=[],
                        Omega0_b=Ob0,
                        Omega0_cdm=Om0 - Ob0,
                        h=hubble, 
                        n_s=ns)\
                        .match(sigma8=sigma8)

    pklin = LinearPower(MYcosmo, redshift=0)
    k = np.logspace(-3,2,10000, endpoint=True)
    np.savetxt(output, np.c_[k, pklin(k)])

    return None

def mk_ini_lineark(hubble, Om0, Ob0, ns, sigma8, 
                   Nmesh, boxsize, output, 
                   unitary_amplitude=False, 
                   inverted_phase=False, seed=0):
    MYcosmo = Cosmology(m_ncdm=[],
                        Omega0_b=Ob0,
                        Omega0_cdm=Om0 - Ob0,
                        h=hubble, 
                        n_s=ns)\
                        .match(sigma8=sigma8)
    
    pklin = LinearPower(MYcosmo, redshift=0)
    mesh = LinearMesh(Plin=pklin, BoxSize=boxsize, Nmesh=Nmesh, seed=seed, unitary_amplitude=unitary_amplitude, inverted_phase=inverted_phase)
    mesh.save(output, dataset="LinearDensityK", mode="complex")

    return mesh

def mk_fastpm_conf(conf, icosmo, Om0, hubble, pkpath, output, relic="", lineark=False):
    ### FPM general params
    fpm_params = fastpm_default.copy()
    fpm_params["boxsize"] = conf["FastPM"].getfloat("boxsize")
    fpm_params["nc"] = conf["FastPM"].getint("npart")
    fpm_params["time_step"] = str(conf.get("FastPM", "time_step")).strip("\"")
    redshifts = list(map(str, conf.get("FastPM", "redshifts").split(", ")))
    fpm_params["output_redshifts"] = "{{{}}}".format(" ".join(redshifts))

    ### FPM varied params
    fpm_seedini = conf["FastPM"].getint("seedini")
    fpm_seed = fpm_seedini + icosmo

    snapdir = str(conf.get("FastPM", "snapdir")).strip("\"")
    snapbase = str(conf.get("FastPM", "snapbase")).strip("\"")
    if not os.path.isdir(snapdir):
        os.makedirs(snapdir)

    FOF = conf["FastPM"].getboolean("FOF")
    if FOF:
        if "fof_nmin" in conf.options["FastPM"]:
            fof_nmin = conf["FastPM"].getint("fof_nmin")
        else:
            fof_nmin = 20

    fpm_params["Omega_m"] = Om0
    fpm_params["hubble"] = hubble
    if lineark:
        fpm_params["read_lineark"] = repr(pkpath)
        fpm_params["read_powerspectrum"] = None
    else:
        fpm_params["read_powerspectrum"] = repr(pkpath)
        fpm_params["random_seed"] = fpm_seed
        if "IFFIX" in conf["FastPM"]:
            fpm_params["remove_cosmic_variance"]

    snappath = snapdir+snapbase+f"{relic}/a"
    fpm_params["write_snapshot"] = repr(snappath)
    if FOF:
        fpm_params["write_fof"] = repr(snappath)
        fpm_params["fof_nmin"] = fof_nmin
    
    f = open(output, "w+", encoding="utf-8")
    for key in fpm_params.keys():
        if fpm_params[key] is not None:
            if key == 'hubble':
                f.write(f"h = {fpm_params[key]}\n")
            else:
                f.write(f"{key} = {fpm_params[key]}\n")

    f.close()

    return None

def mk_rockstar_conf(conf, cvt_opbase, cvt_nfile, filename, output, relic="", scale_factor=None, redshift=None):
    rstar_params = rockstar_default.copy()
    if scale_factor is None and redshift is None:
        print("Must input scale factor OR redshift of snapshot!")
        exit()
    elif scale_factor is not None and redshift is not None:
        if np.abs((1./(1+redshift))-scale_factor) > 1e-5:
            print("Input redshift and scale factor do not match! Use redshift.")
            scale_factor = 1./(1+redshift)
    elif scale_factor is None:
        scale_factor = 1./(1+redshift)

    op_base = conf.get("ROCKSTAR", "outputbase").strip("\"")
    snapdir = str(conf.get("FastPM", "snapdir")).strip("\"")
    snapbase = str(conf.get("FastPM", "snapbase")).strip("\"")
    snappath = snapdir+snapbase+"{}/a_{:.4f}/".format(relic, scale_factor)
    rstar_params["FORCE_RES"]    = conf["ROCKSTAR"].getfloat("force_res")
    rstar_params["PARALLEL_IO"]  = conf["ROCKSTAR"].getint("parallel")
    rstar_params["FILENAME"]     = repr(filename)
    rstar_params["INBASE"]       = repr(snappath+cvt_opbase)
    rstar_params["OUTBASE"]      = repr(snappath+op_base)

    if rstar_params["PARALLEL_IO"]:
        rstar_params["NUM_BLOCKS"] = cvt_nfile
        rstar_params["NUM_READERS"] = min(cvt_nfile, conf["ROCKSTAR"].getint("num_readers"))
        rstar_params["NUM_WRITERS"] = conf["ROCKSTAR"].getint("num_writers")
        rstar_params["FORK_PROCESSORS_PER_MACHINE"]  = conf["ROCKSTAR"].getint("process_per_machine")

    f = open(output, "w+", encoding="utf-8")
    for key in rstar_params.keys():
        if rstar_params[key] is not None:
            f.write(f"{key} = {rstar_params[key]}\n")

    f.close()
