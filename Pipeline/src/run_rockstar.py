import os, sys
import configparser
import time
from utils.convert import Convert

if len(sys.argv) != 2:
    print("Usage: python run_rockstar.py pipeline_conf")
    exit()

conf_name = sys.argv[1]
conf = configparser.ConfigParser()
conf.read(conf_name)

### Convert params
nfile = conf["Convert"].getint("nfile")
precision = conf.get("Convert", "precision").strip("\"")
gadgetbase = conf.get("Convert", "outputbase").strip("\"")

### Input & output
ipdir = conf.get("FastPM", "snapdir").strip("\"")
ipbase = conf.get("FastPM", "snapbase").strip("\"")
redshifts = list(map(float, conf.get("FastPM", "redshifts").split(", ")))

### Rockstar part
boxsize = conf["FastPM"].getfloat("boxsize")
base = conf.get("General","cfgbase").strip("\"")
subbase = conf.get("General","cfgsubbase").strip("\"")
ncosmo = conf["General"].getint("ncosmo")
rstaropbase = conf.get("ROCKSTAR", "outputbase").strip("\"")

Rockstar_exec = "/public/home/suchen/applications/rockstar/rockstar"
FindPAR_exec = "/public/home/suchen/applications/rockstar/util/find_parents"

os.environ["OMP_NUM_THREADS"] = "1"
for icosmo in range(1):
    ### convert bigfile catalog to gadget
    for idx, zi in enumerate(redshifts):
        snappath = ipdir+ipbase+"{:d}/a_{:.4f}/".format(icosmo, 1./(1+zi))
        gadgetpath = snappath+gadgetbase+"{:03d}".format(idx)

        ### execute convert
        Convert(snappath, gadgetpath, nfile, precision)

        ### run rockstar
        rstar_cfgpath = base+subbase+"{:d}/rockstar/z{:.2f}.cfg".format(icosmo, zi)
        if not os.path.isdir(snappath+rstaropbase):
            os.mkdir(snappath+rstaropbase)

        ### I don't know why but this way works ...
        ### TODO: add find parent and remove gadget snapshots/tmp halo file/tmp script
        
        f = open("tmp_run_rstar.sh", "w+")
        f.write("#!/bin/bash\n")
        f.write(f"RSTAR={Rockstar_exec}\n")
        f.write(f"CFG={rstar_cfgpath}\n")
        f.write(f"ODIR={snappath+rstaropbase}\n")
        f.write("$RSTAR -c $CFG &\n")
        f.write("export ODIR\n")
        f.write("perl -e \'sleep 1 while (!(-e \"$ENV{ODIR}/auto-rockstar.cfg\"))\'\n")
        f.write("$RSTAR -c ${ODIR}/auto-rockstar.cfg\n")
        f.close()
        os.system("bash tmp_run_rstar.sh")

