import os, sys
import argparse
from utils.runner import ROCKSTAR_Runner

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-s", "--start", help="start index of cosmology", type=int, default=0)
parser.add_argument("-e", "--end", help="end index of cosmology. Minus means the maximum of the index", type=int, default=-1)
args = parser.parse_args()

conf_file = args.conf

runner = ROCKSTAR_Runner()
runner.load_config_file(conf_file)
runner.set_params()

os.environ["OMP_NUM_THREADS"] = "1"
if __name__ == "__main__":
    for icosmo in range(args.start, args.end):
        runner.run(snapname_relic=f"{icosmo}")