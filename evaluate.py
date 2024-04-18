"""Evaluating collected algorithms predictions with respect to the 
ground truth labels."""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyteomics import mgf
from sklearn.metrics import auc
from tqdm import tqdm

from ground_truth_mapper import convert_GT_to_output_format
from metrics import aa_match_metrics, aa_match_batch


# TODO: move to a separate file? 
AA_MASSES = {
    'G': 57.021463735,
    'A': 71.037113805,
    'S': 87.032028435,
    'P': 97.052763875,
    'V': 99.068413945,
    'T': 101.047678505,
    'C+57.021': 160.030644505,
    'L': 113.084064015,
    'I': 113.084064015,
    'N': 114.04292747,
    'D': 115.026943065,
    'Q': 128.05857754,
    'K': 128.09496305,
    'E': 129.042593135,
    'M': 131.040484645,
    'H': 137.058911875,
    'F': 147.068413945,
    'R': 156.10111105,
    'Y': 163.063328575,
    'W': 186.07931298,
    '+42.011': 42.010565,
    '+43.006': 43.005814,
    '-17.027': -17.026549,
    '+43.006-17.027': 25.980265,
    'M+15.995': 147.03539964499998,
    'N+0.984': 115.02694346999999,
    'Q+0.984': 129.04259353999998
 }


def parse_scores(aa_scores: str) -> list[float]:
    """
    TODO.
    * assumes that AA confidence scores always come
    as a string of float numbers separated by a comma.
    """

    aa_scores = aa_scores.split(",")
    aa_scores = list(map(float, aa_scores))
    return aa_scores


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


# Load ground truth labels and convert to the common output format
input_paths = [
    os.path.join(args.data_dir, x)
    for x in os.listdir(args.data_dir)
    if ".mgf" in x
]
input_paths = sorted(
    input_paths
)  # TODO remove? not needed if match by scan or scan_index

sequences_true = {k: [] for k in ["seq", "scans", "scan_indices"]}
for file_i, mgf_path in enumerate(input_paths):
    spectra = mgf.read(mgf_path)
    for spectrum_i, spectrum in tqdm(enumerate(spectra)):
        sequences_true["seq"].append(spectrum["params"]["seq"])
        sequences_true["scans"].append(
            f'F{file_i}:{spectrum["params"]["scans"]}'
        )
        sequences_true["scan_indices"].append(f"F{file_i}:{spectrum_i}")
sequences_true = pd.DataFrame(sequences_true)
sequences_true["seq"] = sequences_true["seq"].apply(
    convert_GT_to_output_format
)

# Load predictions data, match to GT by scan id or scan index if available
width = 8
fig, (pep_ax, aa_ax) = plt.subplots(1, 2, figsize=(2 * width, width))

output_metrics = {}
for output_file in os.listdir(args.output_dir):
    algo_name = output_file.split("_")[0]
    output_path = os.path.join(args.output_dir, output_file)

    output_data = pd.read_csv(output_path)
    use_cols = ["sequence", "score"]
    if "aa_scores" in output_data.columns:
        use_cols.append("aa_scores")

    if "scans" in output_data.columns:
        use_cols.append("scans")
        output_data = pd.merge(
            sequences_true,
            output_data[use_cols],
            on="scans",
            how="left",
        )
    elif "scan_indices" in output_data.columns:
        use_cols.append("scan_indices")
        output_data = pd.merge(
            sequences_true,
            output_data[use_cols],
            on="scan_indices",
            how="left",
        )
    else:  # TODO: keep or replace with exception+skip algorithm?
        output_data = output_data[use_cols]
        output_data["seq"] = sequences_true["seq"].values
    output_data = output_data.rename({"seq": "sequence_true"}, axis=1)

    # Calculate metrics
    output_data = output_data.sort_values("score", ascending=False)
    sequenced_idx = output_data[
        "sequence"
    ].notnull()  # TODO: indicate number of not sequenced peptides?
    aa_matches_batch, n_aa1, n_aa2 = aa_match_batch(
        output_data["sequence"][sequenced_idx],
        output_data["sequence_true"][sequenced_idx],
        AA_MASSES,
    )

    # Collect metrics
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

    # Plot the peptide precision–coverage curve
    pep_matches = np.array([aa_match[1] for aa_match in aa_matches_batch])
    precision = np.cumsum(pep_matches) / np.arange(1, len(pep_matches) + 1)
    coverage = np.arange(1, len(pep_matches) + 1) / len(pep_matches)
    pep_ax.plot(
        coverage,
        precision,
        label=f"{algo_name} AUC = {auc(coverage, precision):.3f}",
    )

    # Plot the amino acid precision–coverage curve (if aa_scores available)
    if "aa_scores" in output_data.columns:
        aa_scores = np.concatenate(
            list(
                map(
                    parse_scores,
                    output_data["aa_scores"][sequenced_idx].values.tolist(),
                )
            )
        )
        sort_idx = np.argsort(aa_scores)[::-1]
        aa_matches_pred = np.concatenate(
            [aa_match[2][0] for aa_match in aa_matches_batch]
        )
        precision = np.cumsum(aa_matches_pred[sort_idx]) / np.arange(
            1, len(aa_matches_pred) + 1
        )
        coverage = np.arange(1, len(aa_matches_pred) + 1) / len(
            aa_matches_pred
        )
        aa_ax.plot(
            coverage,
            precision,
            label=f"{algo_name} AUC = {auc(coverage, precision):.3f}",
        )

pep_ax.set_title("Peptide precision & coverage")
pep_ax.set_xlim(0, 1), pep_ax.set_ylim(0, 1)
pep_ax.set_xlabel("Coverage"), pep_ax.set_ylabel("Precision")
pep_ax.legend(loc="upper right")

aa_ax.set_title("Amino acid precision & coverage")
aa_ax.set_xlim(0, 1), aa_ax.set_ylim(0, 1)
aa_ax.set_xlabel("Coverage"), aa_ax.set_ylabel("Precision")
aa_ax.legend(loc="upper right")

# Save results
os.makedirs(args.results_dir, exist_ok=True)
plt.savefig(
    os.path.join(args.results_dir, "precision_coverage.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.close()
output_metrics = pd.DataFrame(output_metrics).T
output_metrics.to_csv(os.path.join(args.results_dir, "metrics.csv"))
