import os, sys
import argparse
import numpy as np
from utils.runner import POWER_Runner

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-s", "--start", help="start index of cosmology", type=int, default=0)
parser.add_argument("-e", "--end", help="end index of cosmology. Minus means the maximum of the index", type=int, default=-1)
parser.add_argument("-n", "--nthreads", help="Threads used in measuring power spectrum", type=int, default=1)
args = parser.parse_args()

conf_file = args.conf
runner = POWER_Runner()
runner.load_config_file(conf_file=conf_file)
runner.set_params(
    outputbase = "/public/home/suchen/Programs/Simtool/Pipeline/results/Gpower_rlzs/", 
    snapdir = "/public/share/ace66so15x/suchen/L1000_N1024_rlzs/", 
    snapbase = "rlz", 
    nsham_per_cosmo = 1,
    nrlzs_per_sham = 10,
)

runner.declare()

for irlz in range(int(args.start), int(args.end)):
    runner.run(relic=f"{irlz}")