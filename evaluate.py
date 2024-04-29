"""Evaluating collected algorithms predictions with respect to the 
ground truth labels."""

import argparse
import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pyteomics import mgf
from sklearn.metrics import auc
from tqdm import tqdm

from ground_truth_mapper import format_sequence as format_sequence_GT
from metrics import aa_match_metrics, aa_match_batch
from token_masses import AA_MASSES


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
    help="""
    The path to the directory containing algorithm predictions 
    stored in `algorithm_outputs.csv` files.
    """,
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


# Define dataset name and path to store evaluation results
dataset_name = os.path.basename(os.path.normpath(args.output_dir))
print(f"Evaluating results for {dataset_name}.")

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
    format_sequence_GT
)

# Load predictions data, match to GT by scan id or scan index if available
PLOT_N_POINTS = 10000
PLOT_HEIGHT = 600
PLOT_WIDTH = int(600 * 1.2)

layout = go.Layout(
    height=PLOT_HEIGHT, width=PLOT_WIDTH,
    title_x=0.5,
    margin_t=50,
    xaxis_title="Coverage",
    yaxis_title="Precision",
    xaxis_range=[0, 1],
    yaxis_range=[0, 1],
    legend=dict(
        y=0.01,
        x=0.01,
    )
)
pep_fig = go.Figure(layout=layout)
pep_fig.update_layout(title_text="<b>Peptide precision & coverage</b>")
aa_fig = go.Figure(layout=layout)
aa_fig.update_layout(title_text="<b>AA precision & coverage</b>")
    
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
            sequences_true, output_data[use_cols], on="scans", how="left",
        )
    elif "scan_indices" in output_data.columns:
        use_cols.append("scan_indices")
        output_data = pd.merge(
            sequences_true, output_data[use_cols], on="scan_indices", how="left",
        )
    else: # TODO: keep or replace with exception+skip algorithm? 
        output_data = output_data[use_cols]
        output_data["seq"] = sequences_true["seq"].values
    output_data = output_data.rename({"seq": "sequence_true"}, axis=1)

    # Calculate metrics
    output_data = output_data.sort_values("score", ascending=False)
    sequenced_idx = output_data["sequence"].notnull() # TODO: indicate number of not sequenced peptides?
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
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(np.int64)
    pep_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs], y=precision[plot_idxs],
            mode="lines",
            name=f"{algo_name} AUC = {auc(coverage, precision):.3f}")
    )
    
    # Plot the amino acid precision–coverage curve (if aa_scores available)
    if "aa_scores" not in output_data.columns:
        # define proxy aa_scores from the peptide score
        # TODO: move to algo_name/output_mapper
        output_data.loc[sequenced_idx, "aa_scores"] = output_data[sequenced_idx].apply(
            lambda row: ",".join([str(row["score"]),] * len(re.split(r"(?<=.)(?=[A-Z])", row["sequence"]))), 
            axis=1,
        )
        
    aa_scores = np.concatenate(list(map(
        parse_scores, 
        output_data["aa_scores"][sequenced_idx].values.tolist()
    )))
    sort_idx = np.argsort(aa_scores)[::-1]

    aa_matches_pred = np.concatenate([aa_match[2][0] for aa_match in aa_matches_batch])
    precision = np.cumsum(aa_matches_pred[sort_idx]) / np.arange(1, len(aa_matches_pred) + 1)
    coverage = np.arange(1, len(aa_matches_pred) + 1) / len(aa_matches_pred)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(np.int64)
    aa_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs], y=precision[plot_idxs],
            mode="lines",
            name=f"{algo_name} AUC = {auc(coverage, precision):.3f}")
    )

# Save results
dataset_results_dir = os.path.join(args.results_dir, dataset_name)
os.makedirs(dataset_results_dir, exist_ok=True)

pep_fig.write_html(os.path.join(dataset_results_dir, "peptide_precision_coverage.html"))
aa_fig.write_html(os.path.join(dataset_results_dir, "AA_precision_coverage.html"))

output_metrics = pd.DataFrame(output_metrics).T
output_metrics.to_csv(os.path.join(dataset_results_dir, "metrics.csv"))
