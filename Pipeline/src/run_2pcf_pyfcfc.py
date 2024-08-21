import os, sys
import configparser
import numpy as np
from nbodykit.source.catalog import BigFileCatalog
from utils.tpcf import get_stacked_vprof

if len(sys.argv) != 2:
    print("Usage: python run_rockstar.py pipeline_conf")
    exit()

conf_file = sys.argv[1]
if not os.path.isfile(os.path.join(os.getcwd(), conf_file)):
    print("Configure file does not exist!")
    exit()

conf = configparser.ConfigParser()
conf.read(conf_file)

pyfcfc_path = conf.get("FCFC", "pyFCFC_PATH").strip("\"")
sys.path.append(pyfcfc_path)

### void size bin info
Rmin = conf["FCFC"].getfloat("RVmin")
Rmax = conf["FCFC"].getfloat("RVmax")
dRV  = conf["FCFC"].getfloat("dRV")

Rmins = np.arange(Rmin, Rmax, dRV)
Rmaxs = np.append(Rmins[1:], Rmax)

### separation bins info
min_sep_inRv = conf["FCFC"].getfloat("min_sep_inRv")
max_sep_inRv = conf["FCFC"].getfloat("max_sep_inRv")
n_sep_bins  = conf["FCFC"].getint("n_sep_bins")

### output file
outputbase = conf.get("FCFC", "outputbase").strip("\"")

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

for icosmo in range(ncosmo):
    for zi in redshifts:
        snappath = snapdir+snapbase+"{:d}/a_{:.4f}/".format(icosmo,1./(1.+zi))
        tmp = BigFileCatalog(snappath, dataset="1/", header="Header")
        dm = tmp['Position'].compute()
        del tmp
        for isham in range(nsham_per_cosmo):
            outputpath = outputbase+snapbase+"{:d}/a_{:.4f}/SHAM{:d}/".format(icosmo,1./(1.+zi),isham)
            if not os.path.isdir(outputpath):
                os.makedirs(outputpath)
            xi_iso_mean = []
            xi_iso_stacked_mean = []
            for irlz in range(nrlzs_per_sham):
                voidpath = snappath+halobase+shambase+f"{isham}/void_{feature}_rlz{irlz}.txt"
                void = np.loadtxt(voidpath)
                results_list = get_stacked_vprof(dm, void, boxsize, Rmins, Rmaxs,
                                          min_sep_inRv=min_sep_inRv, 
                                          max_sep_inRv=max_sep_inRv, 
                                          nbins=n_sep_bins,
                                          ) #ind_file=outputpath, stack_file=outputpath

                xi_iso_mean.append(results_list["xi_iso_list"])
                xi_iso_stacked_mean.append(results_list["xi_iso_stacked"])
            
            xi_iso_mean = np.mean(np.asarray(xi_iso_mean), axis=0)
            xi_iso_stacked_mean = np.mean(np.asarray(xi_iso_stacked_mean), axis=0)

            for iR, ind_xi_iso in enumerate(xi_iso_mean):
                f = open(outputpath+f"rvbin{iR}.txt", "w+")
                f.write("# sep (Mpc/h) xi_iso\n")
                np.savetxt(f, ind_xi_iso)
                f.close()
            f = open(outputpath+"stacked.txt", "w+")
            f.write("# Rmin={:.2f} Rmax={:.2f} Nbins={:d}\n".format(Rmin, Rmax, len(Rmins)))
            f.write("# sep (Rv) xi_iso\n")
            np.savetxt(f, xi_iso_stacked_mean)
            f.close()
