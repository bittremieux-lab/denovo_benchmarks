"""Evaluating collected algorithms predictions with respect to the 
ground truth labels."""

import argparse
import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from functools import partial
from pyteomics import mgf, fasta
from sklearn.metrics import auc
from tqdm import tqdm

from ground_truth_mapper import format_sequence as format_sequence_GT
from metrics import aa_match_metrics, aa_match_batch
from token_masses import AA_MASSES

# TODO: remove? and pass full path to ref proteome.fasta instead? 
VSC_SCRATCH = "/scratch/antwerpen/209/vsc20960/"
ROOT = os.path.join(VSC_SCRATCH, "benchmarking")
PROTEOMES_DIR = os.path.join(ROOT, "proteomes")
DATASET_TAGS_PATH = os.path.join(ROOT, "denovo_benchmarks", "dataset_tags.tsv")


def parse_scores(aa_scores: str) -> list[float]:
    """
    TODO.
    * assumes that AA confidence scores always come
    as a string of float numbers separated by a comma.
    """

    aa_scores = aa_scores.split(",")
    aa_scores = list(map(float, aa_scores))
    return aa_scores


def remove_ptms(sequence, ptm_pattern="[^A-Z]"):
    return re.sub(ptm_pattern, "", sequence)


def isoleucine_to_leucine(sequence):
    return sequence.replace("I", "L")


def read_fasta(filename):
    # Initialize an empty list to store sequences
    sequences = {}
    # Read the sequences from the FASTA file
    with fasta.read(filename) as fasta_file:
        for record in fasta_file:
            # entry is a tuple (description, sequence)
            description, sequence = record
            
            sequence = isoleucine_to_leucine(sequence)
            sequences[description] = sequence
    return sequences


def find_match_in_proteome(peptide_seq, proteome):
    for ref_id, ref_seq in proteome.items():
        if peptide_seq in ref_seq:
            return True
    return False


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

labels_path = os.path.join(args.data_dir, "labels.csv")
sequences_true = pd.read_csv(labels_path)
sequences_true["seq"] = sequences_true["seq"].apply(format_sequence_GT)

# get database_path from dataset tags (proteome column, by dataset_name)
tags_df = pd.read_csv(DATASET_TAGS_PATH, sep='\t')
tags_df = tags_df.set_index("dataset")
database_path = tags_df.loc[dataset_name, "proteome"]
proteome = read_fasta(os.path.join(PROTEOMES_DIR, database_path))
print(f"Reference proteome length: {len(proteome)} proteins.")
proteome = "|".join(list(proteome.values()))

# Load predictions data, match to GT by scan id or scan index if available
PLOT_N_POINTS = 10000
PLOT_HEIGHT = 400
PLOT_WIDTH = int(PLOT_HEIGHT * 1.2)

layout = go.Layout(
    height=PLOT_HEIGHT,
    width=PLOT_WIDTH,
    title_x=0.5,
    margin_t=50,
    xaxis_title="Coverage",
    yaxis_title="Precision",
    xaxis_range=[0, 1],
    yaxis_range=[0, 1],
    legend=dict(
        y=0.01,
        x=0.01,
    ),
)
prot_match_fig = go.Figure(layout=layout)
prot_match_fig.update_layout(
    title_text="<b>Number of proteome matches vs. score</b>",
    xaxis_title="Score threshold",
    yaxis_title="Number of matches",
    yaxis_range=None,
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
    use_cols = ["sequence", "score", "aa_scores", "spectrum_id"]
    output_data = pd.merge(
        sequences_true,
        output_data[use_cols],
        on="spectrum_id",
        how="outer",
    )
    output_data = output_data.rename({"seq": "sequence_true"}, axis=1)

    # Calculate metrics
    output_data = output_data.sort_values("score", ascending=False)
    sequenced_idx = output_data["sequence"].notnull() # TODO: indicate number of not sequenced peptides?
    labeled_idx = output_data["sequence_true"].notnull()

    output_data["sequence_no_ptm"] = np.nan
    output_data.loc[sequenced_idx, "sequence_no_ptm"] = output_data.loc[sequenced_idx, "sequence"].apply(
        partial(remove_ptms, ptm_pattern='[^A-Z]')
    )
    matches = [
        (isinstance(seq, str)) and (isoleucine_to_leucine(seq) in proteome) # TODO: move to one place with remove_ptms? 
        for seq in tqdm(output_data["sequence_no_ptm"].tolist())
    ]
    output_data["proteome_match"] = matches

    aa_matches_batch, n_aa1, n_aa2 = aa_match_batch(
        output_data["sequence"][sequenced_idx * labeled_idx],
        output_data["sequence_true"][sequenced_idx * labeled_idx],
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
        "N proteome matches": output_data["proteome_match"][sequenced_idx].mean()
    }

    # Plot the peptide precision–coverage curve
    prot_matches = output_data["proteome_match"][sequenced_idx].values
    n_matches = np.cumsum(prot_matches) # ?
    pep_score = output_data["score"][sequenced_idx].values # ?
    plot_idxs = np.linspace(0, len(pep_score) - 1, PLOT_N_POINTS).astype(
        np.int64
    )
    prot_match_fig.add_trace(
        go.Scatter(
            x=pep_score[plot_idxs],
            y=n_matches[plot_idxs],
            mode="lines",
            name=f"{algo_name}",
        )
    )

    # Plot the peptide precision–coverage curve
    pep_matches = np.array([aa_match[1] for aa_match in aa_matches_batch])
    precision = np.cumsum(pep_matches) / np.arange(1, len(pep_matches) + 1)
    coverage = np.arange(1, len(pep_matches) + 1) / len(pep_matches)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(
        np.int64
    )
    pep_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs],
            y=precision[plot_idxs],
            mode="lines",
            name=f"{algo_name} AUC = {auc(coverage, precision):.3f}",
        )
    )

    # Plot the amino acid precision–coverage curve
    aa_scores = np.concatenate(
        list(
            map(
                parse_scores,
                output_data["aa_scores"][sequenced_idx * labeled_idx].values.tolist(),
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
    coverage = np.arange(1, len(aa_matches_pred) + 1) / len(aa_matches_pred)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(
        np.int64
    )
    aa_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs],
            y=precision[plot_idxs],
            mode="lines",
            name=f"{algo_name} AUC = {auc(coverage, precision):.3f}",
        )
    )

# Save results
dataset_results_dir = os.path.join(args.results_dir, dataset_name)
os.makedirs(dataset_results_dir, exist_ok=True)

prot_match_fig.write_html(
    os.path.join(dataset_results_dir, "number_of_proteome_matches.html")
)
pep_fig.write_html(
    os.path.join(dataset_results_dir, "peptide_precision_coverage.html")
)
aa_fig.write_html(
    os.path.join(dataset_results_dir, "AA_precision_coverage.html")
)

output_metrics = pd.DataFrame(output_metrics).T
output_metrics.to_csv(os.path.join(dataset_results_dir, "metrics.csv"))
