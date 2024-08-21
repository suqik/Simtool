import os, sys
import argparse
import configparser
import numpy as np
from mpi4py import MPI
from utils.sham import *

# if len(sys.argv) != 2:
#     print("Usage: python run_rockstar.py pipeline_conf")
#     exit()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-s", "--start", help="start index of cosmology", type=int, default=0)
parser.add_argument("-e", "--end", help="end index of cosmology. Minus means the maximum of the index", type=int, default=-1)
args = parser.parse_args()

conf_file = args.conf
if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
    print("Configure file does not exist!")
    exit()

conf = configparser.ConfigParser()
conf.read(conf_file)

### load SHAM params
ncats = conf["SHAM"].getint("ncats")
nrlzs = conf["SHAM"].getint("nrlzs")
seedini = conf["SHAM"].getint("seedini")

nparams = conf["SHAM"].getint("nparams")
param_names = conf.get("SHAM", "param_names").strip("\"")
prior_low = list(map(float, conf.get("SHAM", "prior_low").split(", ")))
prior_up = list(map(float, conf.get("SHAM", "prior_up").split(", ")))

ref_num_den = conf["SHAM"].getfloat("ref_num_den")
feature = conf.get("SHAM", "feature").strip("\"")
z_space = conf["SHAM"].getboolean("z_space")

galbase = conf.get("SHAM", "outputbase").strip("\"")

### Void Finder path
DIVE_exec = conf.get("DIVE", "DIVE_EXEC").strip("\"")

### load other necessary params
ipdir = conf.get("FastPM", "snapdir").strip("\"")
ipbase = conf.get("FastPM", "snapbase").strip("\"")
redshifts = list(map(float, conf.get("FastPM", "redshifts").split(", ")))

boxsize = conf["FastPM"].getfloat("boxsize")
ncosmo = conf["General"].getint("ncosmo")
halobase = conf.get("ROCKSTAR", "outputbase").strip("\"")

cfgbase = conf.get("General", "cfgbase").strip("\"")
cfgsubbase = conf.get("General", "cfgsubbase").strip("\"")

if args.end < 0 or args.end > ncosmo:
    args.end = ncosmo

### execute ###
# for icosmo in range(args.start, args.end):
icosmo = args.start + rank
while (icosmo < args.end):
    SHAM_paramfile = cfgbase+cfgsubbase+f"{icosmo}/SHAM_list.txt"
    ### FIXME: Check if different snapshots should use the same seed ###
    ### setup random generator ###
    # ### Use spawn, but need python>3.9 and numpy>1.25
    # seed = seedini + icosmo
    # parent_rng = np.random.default_rng(seed=seed)
    # param_rng, rlz_rng = parent_rng.spawn(2)

    ### seed interval eqs to 2
    seed1 = seedini + icosmo*2
    seed2 = seedini + icosmo*2 + 1
    param_rng = np.random.default_rng(seed=seed1)
    rlz_rng = np.random.default_rng(seed=seed2)

    param_list = param_rng.uniform(low=prior_low, high=prior_up, size=(ncats, nparams))
    f = open(SHAM_paramfile, "w+", encoding="utf-8")
    f.write(f"# {param_names}\n")
    np.savetxt(f, param_list)
    f.close()
    
    for zi in (redshifts):
        snappath = ipdir+ipbase+"{:d}/a_{:.4f}/".format(icosmo, 1./(1+zi))
        halo_file = snappath+halobase+"out_0.list"
        halo = load_rockstar_halo(halo_file, feature=feature, z_space=z_space)
        for ipar, iparam in enumerate(param_list):
            gal_void_filebase = snappath+halobase+galbase+f"{ipar}/"
            if not os.path.isdir(gal_void_filebase):
                os.mkdir(gal_void_filebase)
            for irlz in range(nrlzs):
                ### populate galaxy
                igal = SHAM_model(iparam, halo, feature, z_space, ref_num_den, rng=rlz_rng)
                gal_file = gal_void_filebase+f"gal_{feature}_rlz{irlz}.txt"
                fo = open(gal_file, "w+", encoding="utf-8")
                ### TODO: DIVE does not support header line, maybe change io funs
                # if z_space:
                #     fo.write("# x y zrsd\n")
                # else:
                #     fo.write("# x y z\n")
                np.savetxt(fo, igal[:,:3], fmt="%.3f %.3f %.3f")
                fo.close()
                ### find voids
                void_file = gal_void_filebase+f"void_{feature}_rlz{irlz}.txt"
                os.system(f"{DIVE_exec} -i {gal_file} -o {void_file} -u {boxsize}")
    
    icosmo += size
