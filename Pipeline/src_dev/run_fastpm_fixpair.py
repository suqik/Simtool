import argparse
from utils.runner import FASTPM_Runner

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-n", "--nCPUs", help="number of CPUs", type=int)
parser.add_argument("-i", "--idx", help="index of realization", type=int, default=-1)
args = parser.parse_args()

conf_file = args.conf

runner = FASTPM_Runner()
runner.load_config_file(conf_file)
runner.set_params()

runner.run(nCPUs=args.nCPUs, snapname_relic="", iteration=False)

# if args.idx < 0:
#     print("Will automatically run ncosmo simulations given by configuration file.")
#     runner.run(nCPUs=args.nCPUs, iteration=True)
# else:
#     runner.run(nCPUs=args.nCPUs, snapname_relic=f"{args.idx}", iteration=False)