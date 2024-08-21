import os, sys
import argparse
import configparser
# from mpi4py import MPI
from utils.convert import Convert
# from nbodykit.source.catalog import BigFileCatalog
# from utils.my_convert import Convert

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

Rockstar_exec = conf.get("ROCKSTAR", "ROCKSTAR_EXEC").strip("\"")
# FindPAR_exec = "/public/home/suchen/applications/rockstar/util/find_parents"

if args.end < 0 or args.end > ncosmo:
    args.end = ncosmo

os.environ["OMP_NUM_THREADS"] = "1"
if __name__ == "__main__":
    for icosmo in range(args.start, args.end):
        ### convert bigfile catalog to gadget
        for irds, zi in enumerate(redshifts):
            snappath = ipdir+ipbase+"{:d}/a_{:.4f}/".format(icosmo, 1./(1+zi))
            gadgetpath = snappath+gadgetbase+"{:03d}".format(irds)

            ### execute convert
            if not os.path.isfile(gadgetpath+".0"):
                Convert(snappath, gadgetpath, nfile, precision)

            ### run rockstar
            rstar_cfgpath = base+subbase+"{:d}/rockstar/z{:.2f}.cfg".format(icosmo, zi)
            if not os.path.isdir(snappath+rstaropbase):
                os.mkdir(snappath+rstaropbase)

            ### I don't know why but this way works ... 
            tmp_fname = "tmp_run_rstar_cosmo{}.sh".format(icosmo)
            f = open(tmp_fname, "w+")
            
            f.write("#!/bin/bash\n")
            f.write(f"RSTAR={Rockstar_exec}\n")
            # f.write(f"FINDPAR={FindPAR_exec}\n")
            f.write(f"BOXSIZE={boxsize}\n")
            f.write(f"CFG={rstar_cfgpath}\n")
            f.write(f"IDIR={gadgetpath}\n")
            f.write(f"ODIR={snappath+rstaropbase}\n")
            f.write("$RSTAR -c $CFG &\n")
            f.write("export ODIR\n")
            f.write("perl -e \'sleep 1 while (!(-e \"$ENV{ODIR}/auto-rockstar.cfg\"))\'\n")
            f.write("$RSTAR -c ${ODIR}/auto-rockstar.cfg\n")
            # f.write("$FINDPAR ${ODIR}/out_0.list $BOXSIZE >${ODIR}/out_0_wsub.list\n")
            f.write("rm -r ${ODIR}/halos* ${ODIR}/*.cfg ${ODIR}/profiling/\n")
            f.write("rm ${IDIR}*\n")
        
            f.close()
            os.system(f"bash {tmp_fname}")
            os.system(f"rm {tmp_fname}")

