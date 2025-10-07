"""Evaluating collected algorithms predictions with respect to the 
ground truth labels."""

import argparse
import os
import re
import shutil
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import subprocess
from functools import partial
from pyteomics import mgf, proforma
from pyteomics.mass.unimod import Unimod
from sklearn.metrics import auc
from tqdm import tqdm

from . import utils
from . import mmseqs
from . import ground_truth_mapper
from .spectrum_prediction import (
    N_CALIBRATION_PSMS,
    FRAGMENT_MASS_TOL,
    WINDOW_SIZE,
    supported_mods_I,
    map_mods_delta_mass_to_unimod,
    check_supported_by_model,
    predict_intensities,
    predict_RT,
    get_calibration_model,
    calculate_SA,
)
from .metrics import aa_match_metrics, aa_match_batch
from token_masses import AA_MASSES


DATASET_TAGS_PATH = os.environ['DATASET_TAGS_PATH'] 
PROTEOMES_DIR = os.environ['PROTEOMES_DIR']
MMSEQS2_ARGS = [
    "--seed-sub-mat VTML40.out",
    "--comp-bias-corr 0 --mask 0",
    "--spaced-kmer-mode 0",
    "-k 5",
]


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

# Get database_path from dataset tags (proteome column, by dataset_name)
tags_df = pd.read_csv(DATASET_TAGS_PATH, sep='\t')
tags_df = tags_df.set_index("dataset")
database_path = tags_df.loc[dataset_name, "proteome"]
dataset_tags = tags_df.loc[dataset_name]
dataset_tags = tuple(dataset_tags.index[dataset_tags == 1])

# Create directories for MMseqs2 proteome matches search
(
    search_tmp_dir, tmp_files_dir, target_db_dir, query_db_dir, search_result_dir, search_result_path
) = mmseqs.setup_mmseqs_dirs(search_tmp_dir = "./mmseqs2_tmp")

# create a database from a reference proteome
reference_proteome_path = os.path.join(PROTEOMES_DIR, database_path)
contam_path = os.path.join(PROTEOMES_DIR, "crap.fasta")
target_fasta_path = os.path.join(search_tmp_dir, "proteome.fasta")
mmseqs.create_target_fasta(reference_proteome_path, contam_path, target_fasta_path)


# Load GT peptide labels
labels_path = os.path.join(args.data_dir, "labels.csv")
sequences_true = pd.read_csv(labels_path)
sequences_true["seq"] = sequences_true["seq"].apply(ground_truth_mapper.format_sequence)

# Get experimental spectra params
dataset_path = os.path.join(args.data_dir, "mgf")
spectra_params = utils.extract_spectra_params(dataset_path)

# Predict intensities and RT for GT peptides
true_psms = sequences_true.join(spectra_params[["spectrum_id", "charge", "precursor_mass", "true_RT"]].set_index("spectrum_id"), on="spectrum_id")
true_psms[["filename", "idx"]] = true_psms.spectrum_id.str.split(":", expand=True)
# true_psms["seq_unimod"] = true_psms["seq"].apply(map_mods_delta_mass_to_unimod) # only for Prosit models
true_psms["seq_unimod"] = true_psms["seq"].copy() # for other models, keep original format
true_psms_supported_idx = true_psms.apply(
    lambda row: check_supported_by_model(row["seq_unimod"], row["charge"], supported_mods_I),
    axis=1,
)
print("DEBUG: GT peptides supported by prediction model")
print(true_psms_supported_idx.value_counts(), "\n")

# if any(true_psms_supported_idx): TODO: handle this case

# Get intensity predictions for de novo peptides
gt_predictions_mz, gt_predictions_I = predict_intensities(
    true_psms[true_psms_supported_idx].rename({"seq_unimod": "sequence"}, axis=1),
    dataset_tags,
)
# Get RT predictions for de novo peptides and store them in the dataframe
true_psms.loc[true_psms_supported_idx, "pred_RT"] = predict_RT(
    true_psms[true_psms_supported_idx].rename({"seq_unimod": "sequence"}, axis=1)
)
# Calculate spectral angles and RT differences (on calibrated RT)
true_psms["SA"] = 0.
true_psms["true_RT_calib"] = 0.
spectra_params["true_RT_calib"] = 0.
for filename in true_psms["filename"].value_counts().index:
    print(filename)
    # Select calibration PSMs (from GT PSMs)
    true_psms_file_mask = (true_psms["filename"] == filename) & true_psms_supported_idx
    # spec_idx: index - psm idx in dataframe, value - 0-based spectrum idx in mgf file
    spec_idxs = true_psms[true_psms_file_mask]["idx"].astype(np.int64)

    calib_psms = true_psms[true_psms_file_mask]
    calib_psms = calib_psms.sample(n=min(N_CALIBRATION_PSMS, len(calib_psms)), replace=False, random_state=0)
    # Get predictions for calibration PSMs
    calib_psms["pred_RT"] = predict_RT(calib_psms[["seq_unimod"]].rename({"seq_unimod": "sequence"}, axis=1))
    # Train calibration model (for this particular file)
    rt_calib_reg = get_calibration_model(calib_psms)
    # Calibrate true_RT to iRT
    true_psms.loc[true_psms_file_mask, "true_RT_calib"] = rt_calib_reg.predict(
        true_psms.loc[true_psms_file_mask, "true_RT"].values[:, None]
    )[:, 0]
    # Calibrate true_RT for all spectra in the file
    all_spectra_file_mask = (spectra_params["filename"] == filename)
    spectra_params.loc[all_spectra_file_mask, "true_RT_calib"] = rt_calib_reg.predict(
        spectra_params.loc[all_spectra_file_mask, "true_RT"].values[:, None]
    )[:, 0]

    # Load mgf file, iterate through experimental spectra, calculate SA
    mgf_path = os.path.join(dataset_path, filename + ".mgf")
    true_psms.loc[true_psms_file_mask, "SA"] = calculate_SA(spec_idxs, gt_predictions_mz, gt_predictions_I, mgf_path)

# Calculate RT differences (on calibrated RT)
max_true_irt = true_psms["true_RT_calib"].max()
true_psms["RT_diff"] = (true_psms["pred_RT"] - true_psms["true_RT_calib"]).abs() / max_true_irt

# TODO: add database search baselines to RT_diff, SA plots
# TODO: add shuffled RT baseline to RT_diff plot

# Load predictions data, match to GT by scan id or scan index if available
PLOT_N_POINTS = 10000
PLOT_HEIGHT = 440
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
        bgcolor="rgba(255,255,255,0.6)",  # translucent legend background
        font=dict(size=10),
    ),
)
pep_fig = go.Figure(layout=layout)
pep_fig.update_layout(title_text="<b>Peptide precision & coverage</b>")
aa_fig = go.Figure(layout=layout)
aa_fig.update_layout(title_text="<b>AA precision & coverage</b>")
prot_match_fig = go.Figure(layout=layout)
prot_match_fig.update_layout(
    title_text="<b>Number of proteome matches\nvs. number of peptides</b>",
    xaxis_title="Number of predicted peptides",
    yaxis_title="Number of matches",
    xaxis_range=None,
    yaxis_range=None,
) # plot number of matches versus number of predictions? (above some score value?)
rt_diff_fig = go.Figure(layout=layout)
rt_diff_fig.update_layout(
    title_text="<b>Absolute difference between\npredicted peptide RT and experimental RT</b>",
    yaxis_title="RT difference",
)
sa_fig = go.Figure(layout=layout)
sa_fig.update_layout(
    title_text="<b>SA between predicted peptide spectrum\nand experimental spectrum</b>",
    yaxis_title="Spectral angle",
)


output_metrics = {}
for output_file in os.listdir(args.output_dir):
    # algo_name = output_file.split("_")[0]
    algo_name = "_".join(output_file.split("_")[:-1])
    print("EVALUATE", algo_name)

    # Load tool predictions, match with ground truth
    output_path = os.path.join(args.output_dir, output_file)
    output_data = utils.load_predictions(output_path, sequences_true)

    # Get idxs of GT labeled peptides & sequenced peptides (in correct output format)
    print(algo_name)
    print("NaN sequences:", output_data["score"].isnull().sum())
    output_data = output_data.sort_values("score", ascending=False)
    labeled_idx = output_data["sequence_true"].notnull()  
    sequenced_idx = utils.get_sequenced_idx(output_data)
    output_data.loc[~sequenced_idx, "sequence"] = ""
    output_data.loc[~sequenced_idx, "aa_scores"] = ""

    # Use precalculated pred_RT and SA
    # Add experimental spectra data
    output_data = output_data.join(spectra_params.set_index("spectrum_id"), on="spectrum_id")
    # Find predicted sequences supported by the model (Prosit or other)
    # MB no need to recalculate supported_idx, derive it from augment_predictions.py?
    supported_idx = output_data.apply(
        lambda row: check_supported_by_model(row["sequence"], row["charge"], supported_mods_I),
        axis=1,
    )

    # (true_RT_calib is already in spectrum_params (after step for the true_psms above))
    # Calculate RT differences (on calibrated RT) and normalize by max iRT
    output_data["RT_diff"] = (output_data["pred_RT"] - output_data["true_RT_calib"]).abs() / max_true_irt

    # TODO: this is incorrect, should be done on true_psms? (once for all files)
    ## Calculate RT differences for shuffled sequences (supported by the model)
    # true_RT_shuffled = np.random.permutation(output_data.loc[supported_idx, "true_RT_calib"])
    # rt_shuffled_diff = (output_data.loc[supported_idx, "true_RT_calib"] - true_RT_shuffled).abs() / max_true_irt
    # rt_shuffled_diff_wma = np.convolve(
    #     rt_shuffled_diff, 
    #     np.ones(WINDOW_SIZE) / WINDOW_SIZE, 
    #     mode='valid'
    # )

    # Calculate amino acid and peptide-level precision and recall
    # Prepare output sequences for metrics calculation
    output_data.loc[sequenced_idx, ["sequence", "aa_scores"]] = output_data.loc[sequenced_idx].apply(
        lambda row: utils.ptms_to_delta_mass(row["sequence"], row["aa_scores"]),
        axis=1,
        result_type="expand",
    ).values
    # Calculate metrics (aa precision, recall, peptide precision)
    aa_matches_batch, n_aa1, n_aa2 = aa_match_batch(
        output_data["sequence"][labeled_idx],
        output_data["sequence_true"][labeled_idx],
        AA_MASSES,
    )
    aa_precision, aa_recall, pep_precision = aa_match_metrics(aa_matches_batch, n_aa1, n_aa2)

    # Calculate number of proteome matches
    # Create "database" of de novo predicted peptides
    # output_data["sequence_no_ptm"] = np.nan
    query_fasta_path = os.path.join(search_tmp_dir, "denovo_predicted_peptides.fasta")
    output_data.loc[sequenced_idx, "sequence_no_ptm"] = output_data.loc[sequenced_idx, "sequence"].apply(
        partial(remove_ptms, ptm_pattern='[^A-Z]')
    )
    mmseqs.create_query_fasta(output_data.loc[sequenced_idx], query_fasta_path)
    # Run mmseqs search
    search_df = mmseqs.run_mmseqs(
        target_fasta_path,
        query_fasta_path,
        target_db_dir,
        query_db_dir,
        search_result_dir,
        search_result_path,
        tmp_files_dir,
        args=MMSEQS2_ARGS,
    )
    # Compute number of proteome matches
    output_data["proteome_match"] = False
    output_data["proteome_match"].loc[search_df.index] = True
    output_data = output_data.join(search_df) # ["qaln", "taln", "mismatch", "fident", "evalue"]
    n_proteome_matches = output_data["proteome_match"].sum() # TODO: use number or fraction?
    
    # [Debug] Check number of GT peptide matches
    pep_matches = np.array([aa_match[1] for aa_match in aa_matches_batch])
    output_data["pep_match"] = False
    output_data.loc[labeled_idx, "pep_match"] = pep_matches
    
    # Collect metrics
    output_metrics[algo_name] = {
        "N sequences": sequenced_idx.size,
        "N predicted": sequenced_idx.sum(),
        "AA precision": aa_precision,
        "AA recall": aa_recall,
        "Pep precision": pep_precision,
        "N proteome matches": n_proteome_matches, # TODO: use number or fraction of matches?
    }

    # PLOTTING
    # TODO: how can we store calculated values/metrics to not recalculate them?
    # (can we store them in output_data and save to csv?)

    # Plot the RT difference curve
    rt_diff = output_data[supported_idx].sort_values("score", ascending=False)["RT_diff"]
    rt_diff_wma = np.convolve(
        rt_diff, 
        np.ones(WINDOW_SIZE) / WINDOW_SIZE, 
        mode='valid'
    )
    coverage = np.arange(1, len(rt_diff_wma) + 1) / len(rt_diff_wma)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(np.int64)
    rt_diff_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs],
            y=rt_diff_wma[plot_idxs],
            mode="lines",
            name=f"{algo_name}",
        )
    )

    # Plot the SA curve
    SA = output_data[supported_idx].sort_values("score", ascending=False)["SA"]
    SA_wma = np.convolve(
        SA, 
        np.ones(WINDOW_SIZE) / WINDOW_SIZE, 
        mode='valid'
    )
    coverage = np.arange(1, len(SA_wma) + 1) / len(SA_wma)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(np.int64)
    sa_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs],
            y=SA_wma[plot_idxs],
            mode="lines",
            name=f"{algo_name}",
        )
    )

    # Plot the proteome matches vs number of predictions curve
    prot_matches = output_data["proteome_match"][sequenced_idx].values
    n_matches = np.cumsum(prot_matches)
    n_sequenced = np.arange(sequenced_idx.sum())
    plot_idxs = np.linspace(0, len(n_sequenced) - 1, PLOT_N_POINTS).astype(np.int64)
    prot_match_fig.add_trace(
        go.Scatter(
            x=n_sequenced[plot_idxs],
            y=n_matches[plot_idxs],
            mode="lines",
            name=f"{algo_name}",
        )
    )

    # Plot the peptide precision–coverage curve
    pep_matches = np.array([aa_match[1] for aa_match in aa_matches_batch])
    precision = np.cumsum(pep_matches) / np.arange(1, len(pep_matches) + 1)
    coverage = np.arange(1, len(pep_matches) + 1) / len(pep_matches)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(np.int64)
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
                output_data["aa_scores"][labeled_idx].values.tolist(),
            )
        )
    )
    sort_idx = np.argsort(aa_scores)[::-1]

    aa_matches_pred = np.concatenate([aa_match[2][0] for aa_match in aa_matches_batch])
    precision = np.cumsum(aa_matches_pred[sort_idx]) / np.arange(1, len(aa_matches_pred) + 1)
    coverage = np.arange(1, len(aa_matches_pred) + 1) / len(aa_matches_pred)
    plot_idxs = np.linspace(0, len(coverage) - 1, PLOT_N_POINTS).astype(np.int64)
    aa_fig.add_trace(
        go.Scatter(
            x=coverage[plot_idxs],
            y=precision[plot_idxs],
            mode="lines",
            name=f"{algo_name} AUC = {auc(coverage, precision):.3f}",
        )
    )
    
    # [Debug] display number of peptide matches & proteome matches
    print("DEBUG: N GT peptide matches:", output_data["pep_match"].sum())
    print("DEBUG: N proteome matches:", n_proteome_matches)
    idx = output_data["pep_match"] & ~output_data["proteome_match"]
    print("DEBUG: GT peptide matches w/o proteome match:", idx.sum())
    idx = ~output_data["pep_match"] & output_data["proteome_match"]
    print("DEBUG: Proteome matches w/o GT peptide match:", idx.sum())
    print("\n", "=" * 100, "\n")


# Save results
dataset_results_dir = os.path.join(args.results_dir, dataset_name)
os.makedirs(dataset_results_dir, exist_ok=True)

pep_fig.write_html(
    os.path.join(dataset_results_dir, "peptide_precision_coverage.html")
)
aa_fig.write_html(
    os.path.join(dataset_results_dir, "AA_precision_coverage.html")
)
prot_match_fig.write_html(
    os.path.join(dataset_results_dir, "number_of_proteome_matches.html")
)
rt_diff_fig.write_html(
    os.path.join(dataset_results_dir, "RT_difference.html")
)
sa_fig.write_html(
    os.path.join(dataset_results_dir, "SA.html")
)

output_metrics = pd.DataFrame(output_metrics).T
output_metrics.to_csv(os.path.join(dataset_results_dir, "metrics.csv"))


# Clean tmp folders
shutil.rmtree(search_tmp_dir)
