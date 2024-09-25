# read dataset tags from 
# DATASET_TAGS_PATH = os.path.join(ROOT, "denovo_benchmarks", "dataset_tags.tsv")
# get tags for a specific object (passed as argument)
# and pass them to the bash script as a KEY=VALUE pairs
# (in the bash script, $KEY=VALUE variables will be created for each pair
# for (optional) subsequent use within the make_predictions logic)

import os
import argparse
# TODO: should be installed in all the algorithm containers 
# (if this script is used in make_predictions.sh)
import pandas as pd

DATASET_TAGS_PATH = os.environ['DATASET_TAGS_PATH']

parser = argparse.ArgumentParser()
parser.add_argument(
    "--dataset",
    help="Name of the dataset (folder with .mgf files).",
)
args = parser.parse_args()

# Extract properties tags for the dataset
df = pd.read_csv(DATASET_TAGS_PATH, sep='\t')
df = df.set_index("dataset")
dset_tags = df.loc[args.dataset]
dset_tags = dset_tags.to_dict()

# Print the extracted values in a key=value format
# (Expected to be read by make_predictions.sh script)
for key, value in dset_tags.items():
    if key == "proteome":
        # print(f"{key}={value}")
        print("{}={}".format(key, value))
    else:
        # print(f"{key}={int(value)}")
        print("{}={}".format(key, int(value)))
