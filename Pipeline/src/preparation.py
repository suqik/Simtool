import os
import argparse
import configparser
from utils.mk_conf_func import *
from utils.io_func import get_cosmo_params

argpar = argparse.ArgumentParser()
argpar.add_argument("-c", "--conf", help="Pipeline config", type=str)

args = argpar.parse_args()

conf_file = args.conf
if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
    print("Configure file does not exist!")
    exit()

conf = configparser.ConfigParser()
conf.read(conf_file)

### get configuration file path
cfgbase = conf.get("General", "cfgbase").strip("\"")
if not os.path.isdir(cfgbase):
    os.makedirs(cfgbase)
cfgsubbase = conf.get("General", "cfgsubbase").strip("\"")

########################################
### generate cosmological parameters ###
########################################

cosmo_param_dict = get_cosmo_params(conf)
ncosmo = (cosmo_param_dict["ncosmo"])

####################################
### generate configuration files ###
####################################

### Load FastPM configuration file info
fpm_seedini = conf["FastPM"].getint("seedini")
nrlzs_per_cosmo = conf["FastPM"].getint("nrlzs")
snapdir = str(conf.get("FastPM", "snapdir")).strip("\"")

### Load gadget file number & path (used by rockstar cfg).
ROCKSTAR = conf["FastPM"].getboolean("ROCKSTAR")
cvt_nfile  = conf["Convert"].getint("nfile")
cvt_opbase = conf.get("Convert", "outputbase").strip("\"")

cosmo_dict_input = {}
for keys, vals in cosmo_param_dict["fix"].items():
    cosmo_dict_input[keys] = vals

for icosmo in range(ncosmo):
    for irlz in range(nrlzs_per_cosmo):
        ### Generate FastPM configuration file ###
        fpm_cfgpath = os.path.join(cfgbase, cfgsubbase+f"{icosmo}/rlz{irlz}/fastpm/")
        if not os.path.isdir(fpm_cfgpath):
            os.makedirs(fpm_cfgpath)

        # update cosmological parameters
        for key, val in cosmo_param_dict["vari"].items():
            cosmo_dict_input[key] = val[icosmo]

        # calculate initial matter power spectrum
        mk_ini_Pk(cosmo_dict_input, fpm_cfgpath+"Pkini.txt")

        # write fastpm conf
        fpm_seed = fpm_seedini + icosmo # FIXME: can be changed
        snappath = os.path.join(snapdir, cfgsubbase+f"{icosmo}/rlz{irlz}/a")
        mk_fastpm_conf(conf, fpm_seed, cosmo_dict_input, snappath, fpm_cfgpath+"Pkini.txt", fpm_cfgpath+"fpm.lua")

        ### Generate Rockstar configuration file ###
        if ROCKSTAR:
            rstar_cfgpath = os.path.join(cfgbase, cfgsubbase+f"{icosmo}/rlz{irlz}/rockstar/")
            if not os.path.isdir(rstar_cfgpath):
                os.makedirs(rstar_cfgpath)
            redshifts = conf_get_list(conf, "FastPM", "redshifts", float, sep=", ")
            for idx, zi in enumerate(redshifts):
                scale_factor = 1./(1+zi)
                mk_rockstar_conf(conf, snappath+f"_{scale_factor:.4f}", cvt_opbase.split("/")[0], 
                                cvt_nfile, cvt_opbase.split("/")[1]+f"{idx:03d}.<block>", 
                                output=os.path.join(rstar_cfgpath,f"z{zi:.2f}.cfg"))