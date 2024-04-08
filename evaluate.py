"""Evaluating collected algorithms predictions with respect to the 
ground truth labels."""

import argparse
import os
import depthcharge
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyteomics import mgf
from sklearn.metrics import auc

from ground_truth_mapper import convert_GT_to_output_format
from metrics import aa_match_metrics, aa_match_batch


parser = argparse.ArgumentParser()
parser.add_argument(
    "output_dir",
    help="The path to the directory containing algorithm predictions stored in `algorithm_outputs.csv` files.",
)
parser.add_argument(
    "data_dir", help="The path to the input data with ground truth labels."
)
parser.add_argument(
    "--results_dir",
    default="results/",
    help="The path to save evaluation results (default: 'results/').",
)
args = parser.parse_args()


input_paths = [
    os.path.join(args.data_dir, x)
    for x in os.listdir(args.data_dir)
    if ".mgf" in x
]
input_paths = sorted(input_paths)

# convert ground truth labels to the common format
sequences_true = {k: [] for k in ["seq", "scans", "scan_indices"]}
for file_i, mgf_path in enumerate(input_paths):
    spectra = mgf.read(mgf_path)
    for spectrum_i, spectrum in enumerate(spectra):
        sequences_true["seq"].append(spectrum["params"]["seq"])
        sequences_true["scans"].append(
            f'F{file_i}:{spectrum["params"]["scans"]}'
        )
        sequences_true["scan_indices"].append(f"F{file_i}:{spectrum_i}")
sequences_true = pd.DataFrame(sequences_true)
sequences_true["seq"] = sequences_true["seq"].apply(
    convert_GT_to_output_format
)

# load predictions data, match to GT by scan id or scan index if available
width = 8  # TODO rename and move to constants?
plt.figure(figsize=(width, width))
output_metrics = {}

for output_file in os.listdir(args.output_dir):
    algo_name = output_file.split("_")[0]
    output_path = os.path.join(args.output_dir, output_file)

    output_data = pd.read_csv(output_path)
    if "scans" in output_data.columns:
        output_data = pd.merge(
            sequences_true,
            output_data[["scans", "sequence"]],
            on="scans",
            how="left",
        )
        output_data = output_data.rename({"seq": "sequence_true"}, axis=1)
    elif "scan_indices" in output_data.columns:
        output_data = pd.merge(
            sequences_true,
            output_data[["scan_indices", "sequence", "score"]],
            on="scan_indices",
            how="left",
        )
        output_data = output_data.rename({"seq": "sequence_true"}, axis=1)
    else:
        output_data = output_data[["sequence"]]
        output_data["sequence_true"] = sequences_true["seq"].values

    # calculate metrics
    output_data = output_data.sort_values("score", ascending=False)
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

    # plot the precision–coverage curve.
    peptide_matches = np.asarray(
        [aa_match[1] for aa_match in aa_matches_batch]
    )
    precision = np.cumsum(peptide_matches) / np.arange(
        1, len(peptide_matches) + 1
    )
    coverage = np.arange(1, len(peptide_matches) + 1) / len(peptide_matches)
    threshold = np.argmax(output_data["score"] < 0)
    plt.plot(
        coverage,
        precision,
        label=f"{algo_name} AUC = {auc(coverage, precision):.3f}",
    )

os.makedirs(args.results_dir, exist_ok=True)

plt.xlim(0, 1), plt.ylim(0, 1)
plt.xlabel("Coverage"), plt.ylabel("Precision")
plt.legend(loc="upper right")
plt.savefig(
    os.path.join(args.results_dir, "precision_coverage.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.close()

output_metrics = pd.DataFrame(output_metrics).T
output_metrics.to_csv(os.path.join(args.results_dir, "metrics.csv"))
