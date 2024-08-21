import os, sys
import configparser
import numpy as np
from utils.mk_conf_func import mk_ini_Pk, mk_fastpm_conf, mk_rockstar_conf

####################################################
### Initial: reading pipeline configuration file ###
####################################################
if len(sys.argv) != 2:
    print("Usage: python preparation.py cfg_file")
    exit()

conf_file = sys.argv[1]
if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
    print("Configure file does not exist!")
    exit()

conf = configparser.ConfigParser()
conf.read(conf_file)

### get configuration file path
cfgbase = conf.get("General", "cfgbase").strip("\"")
cfgsubbase = conf.get("General", "cfgsubbase").strip("\"")

########################################
### generate cosmological parameters ###
########################################

### setup
nparams = conf["General"].getint("nparams")
param_names = list(map(str, conf.get("General", "param_names").split(", ")))
prior_low = list(map(float, conf.get("General", "prior_low").split(", ")))
prior_up  = list(map(float, conf.get("General", "prior_up").split(", ")))

fix_cosmo_names = list(map(str, conf.get("General", "fix_cosmo_names").split(", ")))
fix_cosmo_vals  = list(map(float, conf.get("General", "fix_cosmo_vals").split(", ")))
fix_cosmo = {}
for i in range(len(fix_cosmo_names)):
    fix_cosmo[fix_cosmo_names[i]] = fix_cosmo_vals[i]

VAR_H = False
if "hubble" in fix_cosmo_names:
    hubble = fix_cosmo["hubble"]
else:
    VAR_H = True

VAR_OM = False
if "OmegaM" in fix_cosmo_names:
    Om0 = fix_cosmo["OmegaM"]
else:
    VAR_OM = True

VAR_OB = False
if "Omegab" in fix_cosmo_names:
    Ob0 = fix_cosmo["Omegab"]
else:
    VAR_OB = True

VAR_NS = False
if "ns" in fix_cosmo_names:
    ns = fix_cosmo["ns"]
else:
    VAR_NS = True

VAR_S8 = False
if "S8" in fix_cosmo_names:
    S8 = fix_cosmo["S8"]
else:
    VAR_S8 = True

ncosmo = conf["General"].getint("ncosmo")
cosmo_param_seed = conf["General"].getint("seed")
output = conf.get("General", "output").strip("\"")

### execute
cosmo_param_rng = np.random.default_rng(seed=cosmo_param_seed)
cosmo_params = cosmo_param_rng.uniform(low=prior_low, high=prior_up, size=(ncosmo, nparams))

param_dict = {}
for i in range(nparams):
    param_dict[param_names[i]] = cosmo_params[:,i]

f = open(cfgbase+output, "w+", encoding="utf-8")
f.write("# {}\n".format(" ".join(param_names)))
np.savetxt(f, cosmo_params, fmt="%3f %3f")
f.close()

###TODO: Support read cosmo params from input file.

####################################
### generate configuration files ###
####################################

### Load gadget file number & path (used by rockstar cfg).
ROCKSTAR = conf["FastPM"].getboolean("ROCKSTAR")
cvt_nfile  = conf["Convert"].getint("nfile")
cvt_opbase = conf.get("Convert", "outputbase").strip("\"")

for icosmo in range(ncosmo):
    ### Generate FastPM configuration file ###
    fpm_cfgpath = cfgbase+cfgsubbase+f"{icosmo}/fastpm/"
    if not os.path.isdir(fpm_cfgpath):
        os.makedirs(fpm_cfgpath)

    # calculate initial matter power spectrum
    if VAR_H:
        hubble = param_dict["hubble"][icosmo]
    if VAR_OM:
        Om0 = param_dict["OmegaM"][icosmo]
    if VAR_OB:
        Ob0 = param_dict["Omegab"][icosmo]
    if VAR_NS:
        ns  = param_dict["ns"][icosmo]
    if VAR_S8:
        S8  = param_dict["S8"][icosmo]
    
    sigma8 = S8/np.sqrt(Om0/0.3)

    # mk_ini_Pk(hubble, Om0, Ob0, ns, sigma8, fpm_cfgpath+"Pkini.txt")

    # write fastpm conf
    mk_fastpm_conf(conf, icosmo, Om0, hubble, fpm_cfgpath+"Pkini.txt", fpm_cfgpath+"fpm.lua")

    if ROCKSTAR:
        rstar_cfgpath = cfgbase+cfgsubbase+f"{icosmo}/rockstar/"
        if not os.path.isdir(rstar_cfgpath):
            os.makedirs(rstar_cfgpath)
        redshifts = list(map(float, conf.get("FastPM", "redshifts").split(", ")))
        for idx, zi in enumerate(redshifts):
            mk_rockstar_conf(conf, icosmo, cvt_opbase.split("/")[0], 
                            cvt_nfile, cvt_opbase.split("/")[1]+"{:03d}.<block>".format(idx), 
                            output=rstar_cfgpath+"z{:.2f}.cfg".format(zi), redshift=zi)