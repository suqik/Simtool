import argparse
import configparser
import numpy as np
from mpi4py import MPI
from utils.runner import GAL_VOID_Runner

# if len(sys.argv) != 2:
#     print("Usage: python run_rockstar.py pipeline_conf")
#     exit()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration file of pipeline")
args = parser.parse_args()

conf = args.conf

runner = GAL_VOID_Runner()
runner.load_config_file(conf)

runner.set_params()

tmp = np.loadtxt("/public/home/suchen/Programs/Simtool/Pipeline/cfgs/validation/SHAM_tot_list.txt")
SHAM_param_list = np.array([tmp[rank]])
if rank == 50:
    print(f"cosmo50 has input SHAM parameter of {tmp[rank]}", flush=True)
runner.run(SHAM_param_list=SHAM_param_list, snapname_relic=f"{rank}", seed2=1234)