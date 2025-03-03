import os
import configparser
import numpy as np
from nbodykit.cosmology import LinearPower, Cosmology
from .cfg_params import fastpm_default, rockstar_default

def conf_get_list(conf, section, key, type, sep=", "):
    return list(map(type, conf.get(section, key).split(sep)))

def get_cosmo_params(conf):
    cosmo_param_dict = {}

    # check if `input` in `General` section
    if "input" in conf.options("General"):
        return NotImplementedError

    if "n_fix_params" in conf.options("General"):
        n_fix_params = conf["General"].getint("n_fix_params")
        
        if n_fix_params == 5:
            return NotImplementedError
        else:
            ### load fixed params
            fix_cosmo_names = conf_get_list(conf, "General", "fix_cosmo_names", str, ", ")
            fix_cosmo_vals = conf_get_list(conf, "General", "fix_cosmo_vals", float, ", ")

            cosmo_param_dict["fix"] = {}
            for i in range(n_fix_params):
                cosmo_param_dict["fix"][fix_cosmo_names[i]] = fix_cosmo_vals[i]
            
            ### load varied params
            n_vari_params = conf["General"].getint("n_vari_params")
            vari_cosmo_names = conf_get_list(conf, "General", "vari_cosmo_names", str, ", ")
            prior_low = conf_get_list(conf, "General", "prior_low", float, ", ")
            prior_up = conf_get_list(conf, "General", "prior_up", float, ", ")
            seed = conf["General"].getint("seed")
            ncosmo = conf["General"].getint("ncosmo")
            cosmo_param_rng = np.random.default_rng(seed=seed)
            vari_cosmo_vals = cosmo_param_rng.uniform(low=prior_low, high=prior_up, size=(ncosmo, n_vari_params))

            cosmo_param_dict["vari"] = {}
            for i in range(n_vari_params):
                cosmo_param_dict["vari"][vari_cosmo_names[i]] = vari_cosmo_vals[:,i]

            ### save parameters
            cfgbase = str(conf.get("General", "cfgbase")).strip("\"")
            output = str(conf.get("General", "output")).strip("\"")
            f = open(os.path.join(cfgbase,output), "w+", encoding="utf-8")
            f.write("# {}\n".format(" ".join(vari_cosmo_names)))
            np.savetxt(f, vari_cosmo_vals, fmt="%3f %3f")
            f.close()

            return cosmo_param_dict, ncosmo

def mk_ini_Pk(cosmo_dict_input:dict, output):
    if "sigma8" not in cosmo_dict_input.keys() and "S8" in cosmo_dict_input.keys():
        sigma8 = cosmo_dict_input["S8"]/np.sqrt(cosmo_dict_input["OmegaM"]/0.3)
    else:
        sigma8 = cosmo_dict_input["sigma8"]

    MYcosmo = Cosmology(m_ncdm=[],
                        Omega0_b=cosmo_dict_input["Omegab"],
                        Omega0_cdm=cosmo_dict_input["OmegaM"] - cosmo_dict_input["Omegab"],
                        h=cosmo_dict_input["hubble"], 
                        n_s=cosmo_dict_input["ns"])\
                        .match(sigma8=sigma8)

    pklin = LinearPower(MYcosmo, redshift=0)
    k = np.logspace(-3,2,10000, endpoint=True)
    np.savetxt(output, np.c_[k, pklin(k)])

    return None

def mk_fastpm_conf(conf, seed, cosmo_dict_input:dict, snappath, pkpath, output):
    ### FPM general params
    fpm_params = fastpm_default.copy()
    fpm_params["boxsize"] = conf["FastPM"].getfloat("boxsize")
    fpm_params["nc"] = conf["FastPM"].getint("npart")
    fpm_params["time_step"] = str(conf.get("FastPM", "time_step")).strip("\"")
    redshifts = conf_get_list(conf, "FastPM", "redshifts", str, sep=", ")
    fpm_params["output_redshifts"] = "{{{}}}".format(" ".join(redshifts))

    ### FPM varied params
    # fpm_seedini = conf["FastPM"].getint("seedini")
    fpm_seed = seed

    # snapdir = str(conf.get("FastPM", "snapdir")).strip("\"")
    # # snapbase = str(conf.get("FastPM", "snapbase")).strip("\"")
    # if not os.path.isdir(snapdir):
    #     os.makedirs(snapdir)

    FOF = conf["FastPM"].getboolean("FOF")
    if FOF:
        if "fof_nmin" in conf.options["FastPM"]:
            fof_nmin = conf["FastPM"].getint("fof_nmin")
        else:
            fof_nmin = 20

    fpm_params["Omega_m"] = cosmo_dict_input["OmegaM"]
    fpm_params["hubble"] = cosmo_dict_input["hubble"]
    fpm_params["read_powerspectrum"] = repr(pkpath)
    fpm_params["random_seed"] = fpm_seed
    # snappath = snapdir+snapbase+f"{icosmo}/a"
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

def mk_rockstar_conf(conf, snappath, cvt_opbase, cvt_nfile, filename, output, scale_factor=None, redshift=None):
    rstar_params = rockstar_default.copy()

    op_base = conf.get("ROCKSTAR", "outputbase").strip("\"")
    rstar_params["FORCE_RES"]    = conf["ROCKSTAR"].getfloat("force_res")
    rstar_params["PARALLEL_IO"]  = conf["ROCKSTAR"].getint("parallel")
    rstar_params["FILENAME"]     = repr(filename)
    rstar_params["INBASE"]       = repr(os.path.join(snappath,cvt_opbase))
    rstar_params["OUTBASE"]      = repr(os.path.join(snappath,op_base))

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
