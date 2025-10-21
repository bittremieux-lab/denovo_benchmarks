"""Precompute parameters of de novo predicted spectra for faster evaluation."""

import argparse
import os
import numpy as np
import pandas as pd

from . import utils
# TODO: move to a separate file?
from .spectrum_prediction import (
    N_CALIBRATION_PSMS,
    FRAGMENT_MASS_TOL,
    supported_mods_I,
    supported_mods_rt,
    check_supported_by_model,
    predict_intensities,
    predict_RT,
    calculate_SA,
)


DATASET_TAGS_PATH = os.environ['DATASET_TAGS_PATH'] 


parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_dir",
    help="The path to the output directory containing algorithm outputs.",
)
parser.add_argument(
    "--data_dir", help="The path to the input data directory with spectra in mgf/ subdirectory."
)
parser.add_argument(
    "--algo_name",
    help="The name of the algorithm (used in the output file name). If not provided, "
         "all algorithms in output_dir will be processed.",
    default=None,
)
parser.add_argument(
    "--force",
    help="Force re-computation even if pred_RT and SA already exist in the output file.",
    action="store_true",
)
args = parser.parse_args()

# output_dir="$output_root_dir/$dset_name" always contains dataset_name as the last part
dataset_name = os.path.basename(os.path.normpath(args.output_dir))

# Get database_path from dataset tags (proteome column, by dataset_name)
tags_df = pd.read_csv(DATASET_TAGS_PATH, sep='\t')
tags_df = tags_df.set_index("dataset")
dataset_tags = tags_df.loc[dataset_name]
dataset_tags = tuple(dataset_tags.index[dataset_tags == 1])

output_file = f"{args.algo_name}_output.csv"
output_path = os.path.join(args.output_dir, output_file)
output_data = pd.read_csv(output_path)

# Skip if output_data already has pred_RT and SA columns and not args.force
if ("pred_RT" in output_data.columns) and ("SA" in output_data.columns) and not args.force:
    print(f"Output data {output_path} already has pred_RT and SA columns. Skipping...")
    exit(0)

use_cols = ["sequence", "score", "aa_scores", "spectrum_id"]
output_data = output_data[use_cols]

# Add experimental spectra data (we only need charge and precursor mass)
dataset_path = os.path.join(args.data_dir, "mgf")
spectra_params = utils.extract_spectra_params(dataset_path)
output_data = output_data.join(spectra_params.set_index("spectrum_id"), on="spectrum_id")
# Find predicted sequences supported by the intensity prediction model
supported_I_idx = output_data.apply(
    lambda row: check_supported_by_model(row["sequence"], row["charge"], supported_mods_I),
    axis=1,
)
# Find predicted sequences supported by the RT prediction model
supported_rt_idx = output_data.apply(
    lambda row: check_supported_by_model(row["sequence"], row["charge"], supported_mods_rt),
    axis=1,
)
print("DEBUG: de novo peptides supported by intensity prediction model")
print(supported_I_idx.value_counts(), "\n")
print("DEBUG: de novo peptides supported by RT prediction model")
print(supported_rt_idx.value_counts(), "\n")

output_data["SA"] = np.nan
if any(supported_I_idx):
    # Get intensity predictions for de novo peptides
    predictions_mz, predictions_I = predict_intensities(output_data[supported_I_idx], dataset_tags)
    # Calculate spectral angles with experimental spectra
    for filename in output_data["filename"].value_counts().index:
        print(filename)
        file_mask = (output_data["filename"] == filename) & supported_I_idx
        # spec_idx: index - psm idx in dataframe, value - 0-based spectrum idx in mgf file
        spec_idxs = output_data[file_mask]["idx"].astype(np.int64)

        # Load mgf file, iterate through experimental spectra, calculate SA
        mgf_path = os.path.join(dataset_path, filename + ".mgf")
        output_data.loc[file_mask, "SA"] = calculate_SA(spec_idxs, predictions_mz, predictions_I, mgf_path)
else:
    print("No de novo peptides supported by intensity prediction model. Skipping intensity and SA prediction.")

output_data["pred_RT"] = np.nan
if any(supported_rt_idx):
    # Get RT predictions for de novo peptides and store them in output_data
    output_data.loc[supported_rt_idx, "pred_RT"] = predict_RT(output_data[supported_rt_idx])
else:
    print("No de novo peptides supported by RT prediction model. Skipping RT prediction.")

# Save updated output data with pred_RT and SA
output_data.to_csv(output_path, index=False)