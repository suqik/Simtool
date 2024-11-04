import os, sys
import argparse
import configparser
import numpy as np
import random
from nbodykit.source.catalog import BigFileCatalog
from utils.tpcf import get_stacked_vprof
from utils.power import mk_mesh_cat

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

if args.end < 0 or args.end > ncosmo:
    args.end = ncosmo

for icosmo in range(args.start, args.end):
    for zi in redshifts:
        snappath = snapdir+snapbase+"{:d}/a_{:.4f}/".format(icosmo,1./(1.+zi))
        tmp = BigFileCatalog(snappath, dataset="1/", header="Header")
        data_size = tmp.csize
        dm = tmp['Position'].compute()
        del tmp
        ######################################
        ### for test !!! 
        # dm_w = np.ones(len(dm))
        ### test downsample
        dsample_rate = 0.01
        random.seed(0)
        sample_idx = random.sample(list(np.arange(data_size)), int(dsample_rate*data_size))
        dm = dm[sample_idx]
        dm_w = np.ones(len(sample_idx))
        ### test mesh cat
        # dm_mesh, dm_w = mk_mesh_cat(dm, 512, 1000)
        # del dm
        # dm = dm_mesh
        for isham in np.arange(0,1):
        ######################################
            print(f"Measuring 2pcf of cosmo{icosmo}, redshift={zi}, SHAM{isham}", flush=True)
            outputpath = outputbase+snapbase+"{:d}/a_{:.4f}/SHAM{:d}/".format(icosmo,1./(1.+zi),isham)
            if not os.path.isdir(outputpath):
                os.makedirs(outputpath)
                
            xi_iso_mean = []
            xi_iso_stacked_mean = []
            for irlz in range(nrlzs_per_sham):
                voidpath = snappath+halobase+shambase+f"{isham}/void_{feature}_rlz{irlz}.txt"
                void = np.loadtxt(voidpath)
                results_list = get_stacked_vprof(dm, void, boxsize, Rmins, Rmaxs,
                                        wdm=dm_w,
                                        min_sep_inRv=min_sep_inRv, 
                                        max_sep_inRv=max_sep_inRv, 
                                        nbins=n_sep_bins,
                                        ) #ind_file=outputpath, stack_file=outputpath

                xi_iso_mean.append(results_list["xi_iso_list"])
                xi_iso_stacked_mean.append(results_list["xi_iso_stacked"])
            
            xi_iso_mean = np.mean(np.asarray(xi_iso_mean), axis=0)
            xi_iso_stacked_mean = np.mean(np.asarray(xi_iso_stacked_mean), axis=0)
            print(f"Measuring 2pcf of cosmo{icosmo}, redshift={zi}, SHAM{isham} Done.", flush=True)

            print(f"Saving 2pcf of cosmo{icosmo}, redshift={zi}, SHAM{isham}", flush=True)
            for iR, ind_xi_iso in enumerate(xi_iso_mean):
                f = open(outputpath+f"rvbin{iR}_1percent.txt", "w+")
                f.write("# sep (Mpc/h) xi_iso\n")
                Rv = 0.5*(Rmins[iR]+Rmaxs[iR])
                sep_list = np.logspace(
                    np.log10(min_sep_inRv*Rv),
                    np.log10(max_sep_inRv*Rv),
                    n_sep_bins
                )
                np.savetxt(f, np.c_[sep_list, ind_xi_iso])
                f.close()
            f = open(outputpath+"stacked_1percent.txt", "w+")
            f.write("# Rmin={:.2f} Rmax={:.2f} Nbins={:d}\n".format(Rmin, Rmax, len(Rmins)))
            f.write("# sep (Rv) xi_iso\n")
            sep_list = np.logspace(
                np.log10(min_sep_inRv),
                np.log10(max_sep_inRv),
                n_sep_bins
            )
            np.savetxt(f, np.c_[sep_list, xi_iso_stacked_mean])
            f.close()
            print(f"Saving 2pcf of cosmo{icosmo}, redshift={zi}, SHAM{isham} Done.", flush=True)
