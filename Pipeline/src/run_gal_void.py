import argparse
import configparser
from utils.io_func import *
from core.drivers import gal_void_driver

argpar = argparse.ArgumentParser()
argpar.add_argument("-c", "--conf", help="Pipeline config", type=str)
argpar.add_argument("-cs", "--cosmo_start", help="Staring label of cosmology", type=int, default=0)
argpar.add_argument("-ce", "--cosmo_end", help="Ending label of cosmology, minus means running all", type=int, default=-1)
argpar.add_argument("-crs", "--crlz_start", help="Staring label of cosmo realization", type=int, default=0)
argpar.add_argument("-cre", "--crlz_end", help="Ending label of cosmo realization, minus means running all", type=int, default=-1)
argpar.add_argument("-ss", "--sham_start", help="Staring label of cosmology", type=int, default=0)
argpar.add_argument("-se", "--sham_end", help="Ending label of cosmology, minus means running all", type=int, default=-1)
argpar.add_argument("-srs", "--srlz_start", help="Staring label of sham realization", type=int, default=0)
argpar.add_argument("-sre", "--srlz_end", help="Ending label of sham realization, minus means running all", type=int, default=-1)
argpar.add_argument("--DIVE_PATH", help="absolute path of DIVE executable file.", type=str)

args = argpar.parse_args()

conf = configparser.ConfigParser()
conf.read(args.conf)

sham_param_dict = get_sham_params(conf)
nsham = sham_param_dict["nsham"]
sham_param_names = sham_param_dict["sham_param_names"]

cosmo_param_dict = get_cosmo_params(conf)
ncosmo = cosmo_param_dict["ncosmo"]
nrlzs_per_cosmo = conf["FastPM"].getint("nrlzs")
redshifts = conf_get_list(conf, "FastPM", "redshifts", float, sep=", ")

nrlzs_per_sham = conf["SHAM"].getint("nrlzs")
seedini = conf["SHAM"].getint("seedini")
ref_num_den = conf["SHAM"].getfloat("ref_num_den")
feature = conf.get("SHAM", "feature").strip("\"")
zspace = conf["SHAM"].getboolean("z_space")
outputbase = conf.get("SHAM", "outputbase").strip("\"")

DIVE_exec = args.DIVE_PATH

cosmo_start, cosmo_end = get_start_end(conf, args, "cosmo")
crlz_start, crlz_end = get_start_end(conf, args, "crlz")
sham_start, sham_end = get_start_end(conf, args, "sham")
srlz_start, srlz_end = get_start_end(conf, args, "srlz")

for icosmo in range(cosmo_start, cosmo_end):
    for icrlz in range(crlz_start, crlz_end):
        for zi in redshifts:
            halopath = get_halopath(conf, icosmo, icrlz, zi)
            for isham in range(sham_start, sham_end):
                sham_param_vals = sham_param_dict[f"cosmo{icosmo}"][f"rlz{icrlz}"][f"z{zi:.2f}"][isham]
                # FIXME: Can be moved to very beginning
                seed_offset = 0
                for isrlz in range(srlz_start, srlz_end):
                    galpath = get_galpath(conf, icosmo, icrlz, zi, isham, isrlz, feature)
                    voidpath = get_voidpath(conf, icosmo, icrlz, zi, isham, isrlz, feature)
                    gal_void_dir = os.path.dirname(galpath)
                    if not os.path.isdir(gal_void_dir):
                        os.makedirs(gal_void_dir)
                    gal_void_driver(
                        sham_param_names, 
                        sham_param_vals, 
                        halopath, 
                        halofinder="rockstar", 
                        feature=feature, 
                        ref_num_den=ref_num_den, 
                        seed=seedini + seed_offset, 
                        zspace=zspace, 
                        galpath=galpath, 
                        voidpath=voidpath, 
                        DIVE_exec=DIVE_exec,
                    )