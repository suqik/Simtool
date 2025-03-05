import os
import json

import configparser

def conf_get_list(conf, section, key, type, sep=", "):
    return list(map(type, conf.get(section, key).split(sep)))

# >===============   Get Parameters   =================<
def get_cosmo_params(conf:configparser.ConfigParser):
    if "input" in conf.options("General"):
        cfgbase = conf.get("General", "cfgbase").strip("\"")
        fnamebase = conf.get("General", "input").strip("\"")
        fname = os.path.join(cfgbase, fnamebase)

        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File {fname} does not exist!")
        
        with open(fname, "r") as f:
            cosmo_param_dict = json.load(f)
        
        return cosmo_param_dict

    if "n_fix_params" in conf.options("General"):
        return NotImplementedError
    
def get_sham_params(conf:configparser.ConfigParser):
    if "input" in conf.options("SHAM"):
        cfgbase = conf.get("General", "cfgbase").strip("\"")
        fnamebase = conf.get("SHAM", "input").strip("\"")
        fname = os.path.join(cfgbase, fnamebase)

        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File {fname} does not exist!")
        
        with open(fname, "r") as f:
            sham_param_dict = json.load(f)
        
        return sham_param_dict

# >====================================================<

def get_snappath(conf, icosmo, irlz, zi):
    snapdir = conf.get("FastPM", "snapdir").strip("\"")
    snapbase = conf.get("FastPM", "snapbase").strip("\"")
    snappath = os.path.join(snapdir, snapbase+f"{icosmo:d}/rlz{irlz:d}/a_{(1./(1+zi)):.4f}/")
    return snappath

def get_fpm_cfgpath(conf, icosmo, irlz):
    base = conf.get("General","cfgbase").strip("\"")
    subbase = conf.get("General","cfgsubbase").strip("\"")
    fpm_cfgpath = os.path.join(base, subbase+f"{icosmo}/rlz{irlz}/fastpm/fpm.lua")
    return fpm_cfgpath

def get_gadgetpath(conf, icosmo, irlz, idx, zi):
    snappath = get_snappath(conf, icosmo, irlz, zi)
    gadgetbase = conf.get("Convert", "outputbase").strip("\"")
    gadgetpath = snappath+gadgetbase+"{:03d}".format(idx)
    return gadgetpath

def get_rstar_cfgpath(conf, icosmo, irlz, zi):
    cfgbase = conf.get("General","cfgbase").strip("\"")
    cfgsubbase = conf.get("General","cfgsubbase").strip("\"")
    rstar_cfgpath = os.path.join(cfgbase, cfgsubbase+f"{icosmo:d}/rlz{irlz:d}/rockstar/z{zi:.2f}.cfg")
    return rstar_cfgpath

def get_halopath(conf, icosmo, irlz, zi):
    snappath = get_snappath(conf, icosmo, irlz, zi)
    halobase = conf.get("ROCKSTAR", "outputbase").strip("\"")
    halopath = os.path.join(snappath, halobase)
    return halopath

def get_galpath(conf, icosmo, icrlz, zi, isham, isrlz, feature):
    snappath = get_snappath(conf, icosmo, icrlz, zi)
    galbase = conf.get("SHAM", "outputbase").strip("\"")
    galpath = os.path.join(snappath, galbase+f"{isham}/Gal_{feature}_rlz{isrlz}.txt")
    return galpath

def get_voidpath(conf, icosmo, icrlz, zi, isham, isrlz, feature):
    snappath = get_snappath(conf, icosmo, icrlz, zi)
    galbase = conf.get("SHAM", "outputbase").strip("\"")
    voidpath = os.path.join(snappath, galbase+f"{isham}/Void_{feature}_rlz{isrlz}.txt")
    return voidpath

def get_start_end(conf, args, label):
    if label == "cosmo":
        cdict = get_cosmo_params(conf)
        ntot = cdict["ncosmo"]
        start = args.cosmo_start
        end = args.cosmo_end
    elif label == "crlz":
        ntot = conf["General"].getint("nrlzs")
        try:
            start = args.rlz_start
            end   = args.rlz_end
        except:
            start = args.crlz_start
            end   = args.crlz_end
    elif label == "sham":
        sdict = get_sham_params(conf)
        ntot = sdict["nsham"]
        start = args.sham_start
        end = args.sham_end
    elif label == "srlz":
        ntot = conf["SHAM"].getint("nrlzs")
        start = args.srlz_start
        end   = args.srlz_end

    if end < 0:
        end = ntot

    return start, end

