import argparse
import configparser
from src.utils.io_func import get_fpm_cfgpath, get_start_end
from src.core.drivers import fastpm_driver

argpar = argparse.ArgumentParser()
argpar.add_argument("-c", "--conf", help="Pipeline config", type=str)
argpar.add_argument("-cs", "--cosmo_start", help="Staring label of cosmology", type=int, default=0)
argpar.add_argument("-ce", "--cosmo_end", help="Ending label of cosmology, minus means running all", type=int, default=-1)
argpar.add_argument("-rs", "--rlz_start", help="Staring label of realization", type=int, default=0)
argpar.add_argument("-re", "--rlz_end", help="Ending label of realization, minus means running all", type=int, default=-1)

args = argpar.parse_args()
conf = configparser.ConfigParser()
conf.read(args.conf)

cosmo_start, cosmo_end = get_start_end(conf, args, "cosmo")
rlz_start, rlz_end = get_start_end(conf, args, "crlz")

for icosmo in range(cosmo_start, cosmo_end):
    for irlz in range(rlz_start, rlz_end):
        # base = conf.get("General","cfgbase").strip("\"")
        # subbase = conf.get("General","cfgsubbase").strip("\"")
        # fpm_cfgpath = os.path.join(base, subbase+f"{icosmo}/rlz{irlz}/fastpm/fpm.lua")
        fpm_cfgpath = get_fpm_cfgpath(conf, icosmo, irlz)
        fastpm_driver(conf, fpm_cfgpath)