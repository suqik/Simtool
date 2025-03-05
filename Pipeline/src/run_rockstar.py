import os
import argparse
import configparser
from utils.mk_conf_func import conf_get_list
from utils.io_func import *
from core.convert import Convert
from core.drivers import rockstar_driver

argpar = argparse.ArgumentParser()
argpar.add_argument("-c", "--conf", help="Pipeline config", type=str)
argpar.add_argument("-cs", "--cosmo_start", help="Staring label of cosmology", type=int, default=0)
argpar.add_argument("-ce", "--cosmo_end", help="Ending label of cosmology, minus means running all", type=int, default=-1)
argpar.add_argument("-rs", "--rlz_start", help="Staring label of realization", type=int, default=0)
argpar.add_argument("-re", "--rlz_end", help="Ending label of realization, minus means running all", type=int, default=-1)

args = argpar.parse_args()

args = argpar.parse_args()
conf = configparser.ConfigParser()
conf.read(args.conf)

### Convert params
nfile = conf["Convert"].getint("nfile")
precision = conf.get("Convert", "precision").strip("\"")
# gadgetbase = conf.get("Convert", "outputbase").strip("\"")

### Snapshot path
# snapdir = conf.get("FastPM", "snapdir").strip("\"")
# snapbase = conf.get("FastPM", "snapbase").strip("\"")
ncosmo = conf["General"].getint("ncosmo")
nrlzs_per_cosmo = conf["FastPM"].getint("nrlzs")
redshifts = conf_get_list(conf, "FastPM", "redshifts", float, sep=", ")

### Rockstar part
# cfgbase = conf.get("General","cfgbase").strip("\"")
# cfgsubbase = conf.get("General","cfgsubbase").strip("\"")
# halobase = conf.get("ROCKSTAR", "outputbase").strip("\"")

Rockstar_exec = "/public/home/suchen/applications/rockstar/rockstar"
FindPAR_exec = "/public/home/suchen/applications/rockstar/util/find_parents"

os.environ["OMP_NUM_THREADS"] = "1"

cosmo_start, cosmo_end = get_start_end(conf, args, "cosmo")
rlz_start, rlz_end = get_start_end(conf, args, "crlz")

for icosmo in range(cosmo_start, cosmo_end):
    for irlz in range(rlz_start, rlz_end):
        ### convert bigfile catalog to gadget
        for idx, zi in enumerate(redshifts):
            # snappath = os.path.join(snapdir, snapbase+f"{icosmo:d}/rlz{irlz:d}/a_{(1./(1+zi)):.4f}/")
            snappath = get_snappath(conf, icosmo, irlz, zi)
            # gadgetpath = snappath+gadgetbase+"{:03d}".format(idx)
            gadgetpath = get_gadgetpath(conf, icosmo, irlz, idx, zi)

            ### execute convert
            # Convert(snappath, gadgetpath, nfile, precision)

            ### run rockstar
            # rstar_cfgpath = os.path.join(cfgbase, cfgsubbase+f"{icosmo:d}/rlz{irlz:d}/rockstar/z{zi:.2f}.cfg")
            rstar_cfgpath = get_rstar_cfgpath(conf, icosmo, irlz, zi)
            halopath = get_halopath(conf, icosmo, irlz, zi)
            if not os.path.isdir(halopath):
                os.mkdir(halopath)
            
            tmp_script_name = f"tmp_run_rstar_{icosmo}_{irlz}_{idx}.sh"
            rockstar_driver(rstar_cfgpath, gadgetpath, halopath, Rockstar_exec, FindPAR_exec)

