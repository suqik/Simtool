import os, sys
import argparse
import configparser
import numpy as np
from utils.power import get_powerspec_cat

# if len(sys.argv) != 2:
#     print("Usage: python run_rockstar.py pipeline_conf")
#     exit()

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-s", "--start", help="start index of cosmology", type=int, default=0)
parser.add_argument("-e", "--end", help="end index of cosmology. Minus means the maximum of the index", type=int, default=-1)
parser.add_argument("-n", "--nthreads", help="Threads used in measuring power spectrum", type=int, default=1)
args = parser.parse_args()

conf_file = args.conf
if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
    print("Configure file does not exist!")
    exit()

conf = configparser.ConfigParser()
conf.read(conf_file)

### power spectrum info
Nmesh = conf["POWER"].getint("Nmesh") # position to grid
MAS = conf.get("POWER", "MAS") # Mass Alignment Scheme

kmin = conf["POWER"].getfloat("min_k")
kmax = conf["POWER"].getfloat("max_k")
nk   = conf["POWER"].getint("nk")
multipoles = list(map(int, conf.get("POWER", "multipoles").split(", ")))

### output file
outputbase = conf.get("POWER", "outputbase").strip("\"")

### input file
ncosmo = conf["General"].getint("ncosmo")
nsham_per_cosmo = conf["SHAM"].getint("ncats")
nrlzs_per_sham = conf["SHAM"].getint("nrlzs")
boxsize = conf["FastPM"].getfloat("boxsize")

snapdir = conf.get("FastPM", "snapdir").strip("\"")
snapbase = conf.get("FastPM", "snapbase").strip("\"")
redshifts = list(map(float, conf.get("FastPM", "redshifts").split(", ")))
halobase = conf.get("ROCKSTAR", "outputbase").strip("\"")
shambase = conf.get("SHAM", "outputbase").strip("\"")
feature = conf.get("SHAM", "feature").strip("\"")

if args.end < 0 or args.end > ncosmo:
    args.end = ncosmo
nthreads = args.nthreads

for icosmo in range(args.start, args.end):
    for zi in redshifts:
        snappath = snapdir+snapbase+"{:d}/a_{:.4f}/".format(icosmo,1./(1.+zi))
        for isham in range(nsham_per_cosmo):
            outputpath = outputbase+snapbase+"{:d}/a_{:.4f}/SHAM{:d}/".format(icosmo,1./(1.+zi),isham)
            if not os.path.isdir(outputpath):
                os.makedirs(outputpath)

            Pkout_mean = []
            for irlz in range(nrlzs_per_sham):
                galpath = snappath+halobase+shambase+f"{isham}/gal_{feature}_rlz{irlz}.txt"
                gal = np.loadtxt(galpath)

                kout, Pkout = get_powerspec_cat(gal, Nmesh, boxsize, MAS, kmin, kmax, nk, multipoles, nthreads)
                Pkout_mean.append(Pkout)
            Pkout_mean = np.mean(np.asarray(Pkout_mean), axis=0)

            f = open(outputpath+f"power.txt", "w+")
            f.write("# k (h/Mpc) ")
            for ell in multipoles:
                f.write(f"Pk{ell} ")
            f.write("\n")
            np.savetxt(f, np.c_[kout, Pkout_mean.T])
            # np.savetxt(f, Pkout_mean)
            f.close()
                
