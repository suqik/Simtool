import os
import configparser
import numpy as np
import json
from utils.io_func import conf_get_list

def write_cosmo_params(
        output, 
        fix_cosmo_names, 
        fix_cosmo_vals, 
        vari_cosmo_names,
        vari_cosmo_vals, 
):
    with open(output, "w+") as f:
        f.write("# {}\n".format(" ".join([f"{fix_cosmo_names[i]}={fix_cosmo_vals[i]:.4f}" for i in range(len(fix_cosmo_names))])))
        f.write("# {}\n".format(" ".join(vari_cosmo_names)))
        np.savetxt(f, vari_cosmo_vals, fmt="%3f %3f")

def write_sham_params(
        output, 
        sham_param_dict
):
    with open(output, "w+") as f:
        json.dump(sham_param_dict, f)
    
def gen_cosmo_params(
        output, 
        fix_cosmo_names, 
        fix_cosmo_vals, 
        vari_cosmo_names, 
        prior_low, 
        prior_up, 
        ncosmo, 
        seed, 
):
    cosmo_param_dict = {}
    n_fix_params = len(fix_cosmo_names)
    cosmo_param_dict["fix"] = {}
    for i in range(n_fix_params):
        cosmo_param_dict["fix"][fix_cosmo_names[i]] = fix_cosmo_vals[i]
    
    n_vari_params = len(vari_cosmo_names)
    cosmo_param_rng = np.random.default_rng(seed=seed)
    vari_cosmo_vals = cosmo_param_rng.uniform(low=prior_low, high=prior_up, size=(ncosmo, n_vari_params))

    cosmo_param_dict["vari"] = {}

    for i in range(n_vari_params):
        cosmo_param_dict["vari"][vari_cosmo_names[i]] = vari_cosmo_vals[:,i].tolist()
    
    cosmo_param_dict["ncosmo"] = ncosmo
    with open(output, "w+") as f:
        json.dump(cosmo_param_dict, f)

def gen_sham_params(
        output, 
        param_name, 
        prior_low, 
        prior_up, 
        ncosmo, 
        nrlz_per_cosmo,
        redshifts,
        nsham, 
        seedini
):
    param_name = np.atleast_1d(param_name)
    prior_low = np.atleast_1d(prior_low)
    prior_up  = np.atleast_1d(prior_up)
    redshifts = np.atleast_1d(redshifts)

    n_sham_param = len(param_name)

    sham_param_dict = {}
    seed_offset = 0
    for icosmo in range(ncosmo):
        sham_param_dict[f"cosmo{icosmo}"] = {}
        for irlz in range(nrlz_per_cosmo):
            sham_param_dict[f"cosmo{icosmo}"][f"rlz{irlz}"] = {}
            for zi in redshifts:
                rng = np.random.default_rng(seed=seedini+seed_offset)
                params = rng.uniform(prior_low, prior_up, size=(nsham, n_sham_param))
                sham_param_dict[f"cosmo{icosmo}"][f"rlz{irlz}"][f"z{zi:.2f}"] = params.tolist()
                # for i in range(n_sham_param):
                #     sham_param_dict[f"cosmo{icosmo}"][f"rlz{irlz}"][f"z{zi:.2f}"][param_name[i]] = params[:,i].tolist()
                seed_offset += 1
    
    sham_param_dict["sham_param_names"] = param_name.tolist()
    sham_param_dict["nsham"] = nsham
    with open(output, "w+") as f:
        json.dump(sham_param_dict, f)

if __name__ == "__main__":
    import sys
    conf = configparser.ConfigParser()
    conf.read(sys.argv[1])

    ncosmo = 1

    cfgbase = str(conf.get("General", "cfgbase")).strip("\"")
    if not os.path.isdir(cfgbase):
        os.makedirs(cfgbase)
    coutput = str(conf.get("General", "input")).strip("\"")
    
    gen_cosmo_params(
        os.path.join(cfgbase, coutput), 
        fix_cosmo_names=["Omegab", "hubble", "ns"], 
        fix_cosmo_vals=[0.0491, 0.676, 0.97], 
        vari_cosmo_names=["OmegaM", "S8"], 
        prior_low=[0.2, 0.6], 
        prior_up=[0.4, 0.9], 
        ncosmo=ncosmo, 
        seed=4321
    )

    nrlz_per_cosmo = conf["FastPM"].getint("nrlzs")
    redshifts = conf_get_list(conf, "FastPM", "redshifts", float)
    nsham = 1
    
    soutput = str(conf.get("SHAM", "input")).strip("\"")

    gen_sham_params(
        os.path.join(cfgbase, soutput), 
        param_name=["sigma"],
        prior_low=[0.0], 
        prior_up=[5.0], 
        ncosmo=ncosmo, 
        nrlz_per_cosmo=nrlz_per_cosmo, 
        redshifts=redshifts, 
        nsham=nsham, 
        seedini=1234
    )
    
    # seed_offset = 0
    # for icosmo in range(ncosmo):
    #     for irlz in range(nrlz_per_cosmo):
    #         shambase = os.path.join(cfgbase,cfgsubbase+f"{icosmo}/rlz{irlz}/sham")
    #         if not os.path.isdir(shambase):
    #             os.makedirs(shambase)
    #         for zi in redshifts:
    #             gen_sham_params(
    #                 os.path.join(shambase, f"{zi:.2f}_{soutput}"), 
    #                 param_name=["sigma"], 
    #                 prior_low=[0.0], 
    #                 prior_up=[5.0], 
    #                 nsham=nsham,
    #                 seed=1234+seed_offset
    #             )
    #             seed_offset += 1