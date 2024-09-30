import argparse
from utils.runner import PREPARE_Runner

parser = argparse.ArgumentParser()
parser.add_argument("conf", help="configuration of pipeline")
args = parser.parse_args()

conf_file = args.conf

runner = PREPARE_Runner()
runner.load_config_file(conf_file)
runner.set_params()
runner.run()