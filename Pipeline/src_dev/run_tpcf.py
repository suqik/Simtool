import argparse
from utils.runner import TPCF_Runner

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
parser.add_argument("-s", "--start", help="start index of realization", type=int, default=0)
parser.add_argument("-e", "--end", help="end index of realization", type=int, default=-1)
args = parser.parse_args()

conf_file = args.conf

runner = TPCF_Runner()
runner.load_config_file(conf_file)

runner.set_params()

runner.declare()
start = int(args.start)
end = int(args.end)
if end < 0:
    end = runner.ncosmo
    
for idx in range(start, end):
    runner.run(snapname_relic=f"{idx}")

# runner = TPCF_Runner()
# runner.load_config_file(conf_file)

# runner.set_params(
#     outputbase="/public/home/suchen/Programs/Simtool/Pipeline/results/JK_tests/",
#     nsham_per_cosmo=1,
#     nrlzs_per_sham=1,
#     snapdir="/public/share/ace66so15x/suchen/L1000_N1024_rlzs/",
#     snapbase="rlz",
# )

# runner.declare()
# idx = int(args.idx)
# runner.run(snapname_relic=f"{idx}")