import argparse
from utils.runner import FASTPM_Runner

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-n", "--nCPUs", help="number of CPUs", type=int)
parser.add_argument("-i", "--idx", help="index of realization", type=int, default=-1)
parser.add_argument("-s", "--start", help="start index of realization", type=int, default=-1)
parser.add_argument("-e", "--end", help="end index of realization", type=int, default=-1)
args = parser.parse_args()

conf_file = args.conf

runner = FASTPM_Runner()
runner.load_config_file(conf_file)
runner.set_params()

if args.idx >= 0:
    runner.run(nCPUs=args.nCPUs, snapname_relic=f"{args.idx}", iteration=False)
elif args.start >= 0:
    if args.end > args.start:
        for run_idx in range(args.start, args.end+1):
            runner.run(nCPUs=args.nCPUs, snapname_relic=f"{run_idx}", iteration=False)
    else:
        print("End idx is smaller than start idx. Will only run simulation with start idx.")
        runner.run(nCPUs=args.nCPUs, snapname_relic=f"{args.start}", iteration=False)
else:
    print("Will automatically run ncosmo simulations given by configuration file.")
    runner.run(nCPUs=args.nCPUs, iteration=True)