import argparse
import numpy as np
from mpi4py import MPI
from utils.runner import GAL_VOID_Runner

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration file of pipeline")
args = parser.parse_args()

conf = args.conf

runner = GAL_VOID_Runner()
runner.load_config_file(conf)

runner.set_params()

SHAM_param_list = np.array([1.5])
# if rank == 50:
#     print(f"cosmo50 has input SHAM parameter of {tmp[rank]}", flush=True)
runner.run(SHAM_param_list=SHAM_param_list, snapname_relic=f"{rank}", seed2=1234)