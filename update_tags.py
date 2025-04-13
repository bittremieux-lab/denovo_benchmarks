import argparse
from dataset_utils import collect_dataset_tags
from dataset_config import get_config


# Paths parsing
parser = argparse.ArgumentParser()
parser.add_argument("--config_path", type=str, help="path to the dataset config file")
args = parser.parse_args()

# Config setup
config = get_config(args.config_path)
dset_id = config.download.dset_id
dset_name = config.name

print(f"Updating tags for {dset_name} ({dset_id})")

# Add dataset tags to the dataset_tags file
collect_dataset_tags(config)