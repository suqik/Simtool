import argparse
import configparser
from utils.io_func import get_fpm_cfgpath, get_start_end
from core.drivers import fastpm_driver

argpar = argparse.ArgumentParser()
argpar.add_argument("-c", "--conf", help="Pipeline config", type=str)
argpar.add_argument("-cs", "--cosmo_start", help="Staring label of cosmology", type=int, default=0)
argpar.add_argument("-ce", "--cosmo_end", help="Ending label of cosmology, minus means running all", type=int, default=-1)
argpar.add_argument("-crs", "--crlz_start", help="Staring label of realization", type=int, default=0)
argpar.add_argument("-cre", "--crlz_end", help="Ending label of realization, minus means running all", type=int, default=-1)

args = argpar.parse_args()
conf = configparser.ConfigParser()
conf.read(args.conf)

cosmo_start, cosmo_end = get_start_end(conf, args, "cosmo")
rlz_start, rlz_end = get_start_end(conf, args, "crlz")

FastPM_exec = "/home/suqikuai777/applications/Simulations/fastpm/src/fastpm"

for icosmo in range(cosmo_start, cosmo_end):
    for irlz in range(rlz_start, rlz_end):
        fpm_cfgpath = get_fpm_cfgpath(conf, icosmo, irlz)
        fastpm_driver(fpm_cfgpath, FastPM_exec=FastPM_exec, nCPUs=16)