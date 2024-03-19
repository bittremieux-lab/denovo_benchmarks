import argparse
import os
import depthcharge
import pandas as pd
from pyteomics import mgf

from ground_truth_mapper import convert_GT_to_output_format
from metrics import aa_match_metrics, aa_match_batch


parser = argparse.ArgumentParser()
parser.add_argument("output_dir", help="")
parser.add_argument("data_dir", help="")
args = parser.parse_args()


input_paths = [
    os.path.join(args.data_dir, x) for x in os.listdir(args.data_dir) if ".mgf" in x
]
input_paths = sorted(input_paths)

# convert ground truth labels to the common format
sequences_true = {k: [] for k in ["seq", "scans"]}
for file_i, mgf_path in enumerate(input_paths):
    spectra = mgf.read(mgf_path)
    for spectrum in spectra:
        sequences_true["seq"].append(spectrum["params"]["seq"])
        sequences_true["scans"].append(f'F{file_i}:{spectrum["params"]["scans"]}')
sequences_true = pd.DataFrame(sequences_true)
sequences_true["seq"] = sequences_true["seq"].apply(convert_GT_to_output_format)

# load predictions data, match to GT by scan id if available
output_metrics = {}

for output_file in os.listdir(args.output_dir):
    algo_name = output_file.split("_")[0]
    output_path = os.path.join(args.output_dir, output_file)

    output_data = pd.read_csv(output_path)
    if "scans" in output_data.columns:
        output_data = pd.merge(
            sequences_true, output_data[["scans", "sequence"]], on="scans", how="left"
        )
        output_data = output_data.rename({"seq": "sequence_true"}, axis=1)
    else:
        output_data = output_data[["sequence"]]
        output_data["sequence_true"] = sequences_true["seq"].values

    # calculate metrics
    sequenced_idx = output_data["sequence"].notnull()

    aa_matches_batch, n_aa1, n_aa2 = aa_match_batch(
        output_data["sequence_true"][sequenced_idx],
        output_data["sequence"][sequenced_idx],
        depthcharge.masses.PeptideMass("massivekb").masses,
    )

    aa_precision, aa_recall, pep_precision = aa_match_metrics(
        aa_matches_batch, n_aa1, n_aa2
    )
    output_metrics[algo_name] = {
        "N sequences": sequenced_idx.size,
        "N predicted": sequenced_idx.sum(),
        "AA precision": aa_precision,
        "AA recall": aa_recall,
        "Pep precision": pep_precision,
    }

output_metrics = pd.DataFrame(output_metrics).T
output_metrics.to_csv("./metrics.csv")
