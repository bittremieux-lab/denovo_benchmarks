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

from ground_truth_mapper import format_sequence as format_sequence_GT
from metrics import aa_match_metrics, aa_match_batch
from token_masses import AA_MASSES

DATASET_TAGS_PATH = os.environ['DATASET_TAGS_PATH'] 
PROTEOMES_DIR = os.environ['PROTEOMES_DIR']

UNIMOD_DB = Unimod()
ptm_masses = {}

def _transform_match_ptm(match: re.Match) -> str:
    """
    TODO
    """
    ptm_idx = int(match.group(1))
    
    if ptm_idx not in ptm_masses:
        ptm_masses[ptm_idx] = UNIMOD_DB.get(ptm_idx).monoisotopic_mass
        print(ptm_masses)
    
    ptm_mass = str(ptm_masses[ptm_idx])
    if not ptm_mass.startswith("-"):
        ptm_mass = "+" + ptm_mass
    return f"[{ptm_mass}]"

def ptms_to_delta_mass(sequence):
    """Convert PTM representation from Unimod notation to delta mass notation."""
    
    PTM_PATTERN = r"\[UNIMOD:([0-9]+)\]" # find ptms
    sequence = re.sub(PTM_PATTERN, _transform_match_ptm, sequence)
    return sequence


def parse_scores(aa_scores: str) -> list[float]:
    """
    TODO.
    * assumes that AA confidence scores always come
    as a string of float numbers separated by a comma.
    """
    if not aa_scores:
        return []
    aa_scores = aa_scores.split(",")
    aa_scores = list(map(float, aa_scores))
    return aa_scores


def remove_ptms(sequence, ptm_pattern="[^A-Z]"):
    return re.sub(ptm_pattern, "", sequence)


def isoleucine_to_leucine(sequence):
    return sequence.replace("I", "L") if isinstance(sequence, str) else sequence


def get_n_tokens(sequence: str) -> int:
    """Calculate number of tokens in a sequence in ProForma notation."""
    seq = proforma.parse(sequence.replace("][", ", "))
    n_tokens = len(seq[0])
    if seq[1]["n_term"]:
        n_tokens += len(seq[1]["n_term"])
    if seq[1]["c_term"]:
        n_tokens += len(seq[1]["c_term"])
    return n_tokens


# Auxiliary methods to validate output format of predicted sequences and aa scores

def get_n_scores(scores: str) -> int:
    return len(scores.split(","))


def validate_spectrum_id(spectrum_id: str) -> bool:
    SPECTRUM_ID_PATTERN = r"[^:]+:\d+"
    return bool(re.fullmatch(SPECTRUM_ID_PATTERN, spectrum_id))


def validate_sequence(sequence: str) -> bool:
    try:
        # Merge subsequent modifications together, because pyteomics.proforma 
        # cannot parse sequence with multiple N-term modifications
        # TODO: maybe only merge _terminal_ modifications?
        seq = proforma.parse(sequence.replace("][", ", "))
    except:
        return False
    return True


def validate_token_scores(scores: str, sequence: str) -> bool:
    n_tokens = get_n_tokens(sequence)
    return get_n_scores(scores) == n_tokens


def load_predictions(output_path, sequences_true):
    """Load de novo predictions and combine them with ground truth sequences."""
    use_cols = ["sequence", "score", "aa_scores", "spectrum_id"]
    
    output_data = pd.read_csv(output_path)
    output_data = pd.merge(
        sequences_true,
        output_data[use_cols],
        on="spectrum_id",
        how="outer",
    )
    output_data = output_data.rename({"seq": "sequence_true"}, axis=1)
    return output_data


def get_sequenced_idx(output_data):
    """Validate predictions and collect indices of correctly formatted predicted sequences."""
    
    sequenced_idx = output_data["sequence"].notnull()
    failed_idx = []
    print("Failed sequences:\n")
    for row_idx, row in output_data[sequenced_idx].iterrows():

        if not validate_sequence(row["sequence"]):
            failed_idx.append(row_idx)
            print(f"Spectrum id: {row['spectrum_id']}")
            print(f"Predicted sequence: {row['sequence']}")
            print(f"FAILED: sequence is not in the ProForma format.\n")

        elif not validate_token_scores(row["aa_scores"], row["sequence"]):
            failed_idx.append(row_idx)
            print(f"{row['spectrum_id']}")
            print(f"Predicted sequence: {row['sequence']}, number of tokens: {get_n_tokens(row['sequence'])}")
            print(f"Predicted scores: {row['aa_scores']}, number of scores: {get_n_scores(row['aa_scores'])}")
            print(f"FAILED: number of per-token scores (','-separated scores in `aa_scores`) does not match number of individual tokens in the predicted sequence.\n")

    sequenced_idx[failed_idx] = False
    return sequenced_idx


def create_query_fasta(output_data, query_fasta_path):
    """
    Create query database of de novo predicted peptides 
    (with I replaced by L)
    """
    
    with open(query_fasta_path, "w") as f:
        f.write(
            "\n".join([
                f">{idx}\n{isoleucine_to_leucine(peptide)}" # I -> L
                for idx, peptide 
                in output_data.sequence_no_ptm.reset_index().values
            ])
        )
    print("queryDB (fasta):", query_fasta_path)


def run_mmseqs(
    target_fasta_path,
    query_fasta_path,
    target_db_dir,
    query_db_dir,
    search_result_dir,
    search_result_path,
    tmp_files_dir,
    args=[],
):
    SEARCH_COLS = "query,target,qaln,taln,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    QUERY_KEY = "target"
    
    print("\n CLEAN EXISTING RESULTS")
    cmd = ["rm -rf", search_result_dir + "/*"]
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    print("\n CREATE TARGET DB")
    target_db_path = os.path.join(target_db_dir, "targetDB")
    cmd = [
        "mmseqs",
        "createdb",
        target_fasta_path,
        target_db_path,
        "-v 1",
    ]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    print("\n CREATE QUERY DB")
    query_db_path = os.path.join(query_db_dir, "queryDB")
    cmd = [
        "mmseqs",
        "createdb",
        query_fasta_path,
        query_db_path,
        "-v 1",
    ]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    print("\n SEARCH/MAP SEQUENCES")
    # mmseqs map <i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir> [options]
    cmd = [
        "mmseqs",
        "map",
        target_db_path,
        query_db_path,
        os.path.join(search_result_dir, "search_result"),
        tmp_files_dir,
        "-a", # backtrace alignments 
        "-e inf",
        "--cov-mode 1",
#         "-v 1",
    ] + args
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)

    print("\n CONVERT RESULTS TO .m8")
    # mmseqs convertalis queryDB targetDB resultDB resultDB.m8
    cmd = [
        "mmseqs",
        "convertalis",
        target_db_path,
        query_db_path,
        os.path.join(search_result_dir, "search_result"),
        search_result_path,
        f"--format-output {SEARCH_COLS}",
        "-v 1"
    ]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    search_df = pd.read_csv(search_result_path, sep="\t", names=SEARCH_COLS.split(","))
    search_df = search_df.drop_duplicates(QUERY_KEY)
    search_df = search_df.set_index(search_df[QUERY_KEY].apply(int))
    search_df = search_df[search_df.mismatch <= 2] # filter s.t. n_mismatches <= 2
    return search_df


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

# Load GT peptide labels
labels_path = os.path.join(args.data_dir, "labels.csv")
sequences_true = pd.read_csv(labels_path)
sequences_true["seq"] = sequences_true["seq"].apply(format_sequence_GT)


# Create directories for MMseqs2 proteome matches search
# TODO: move to a function
# dir for MMseqs2 files
search_tmp_dir = "./mmseqs2_tmp"
os.makedirs(search_tmp_dir, exist_ok=True) 
# dir for MMseqs2 runs tmp files
tmp_files_dir = os.path.join(search_tmp_dir, "tmp")
os.makedirs(tmp_files_dir, exist_ok=True)
# dir target DB (multiple files)
target_db_dir = os.path.join(search_tmp_dir, "targetDB")
os.makedirs(target_db_dir, exist_ok=True)
# dir query DB (multiple files)
query_db_dir = os.path.join(search_tmp_dir, "queryDB")
os.makedirs(query_db_dir, exist_ok=True)
# dir raw search results (multiple files)
search_result_dir = os.path.join(search_tmp_dir, "results")
os.makedirs(search_result_dir, exist_ok=True)
# path to search results in .m8 format
search_result_path = os.path.join(search_tmp_dir, "search_result.m8")

# create a database from a reference proteome
# - join refernce proteome with common contaminants
# - replace I with L (indistinguishable by mass)
reference_proteome_path = os.path.join(PROTEOMES_DIR, database_path)
contam_path = os.path.join(PROTEOMES_DIR, "crap.fasta")
target_fasta_path = os.path.join(search_tmp_dir, "proteome.fasta")

cmd = [
    "cat",
    reference_proteome_path,
    contam_path,
    "| sed -e 's/I/L/g'",
    ">",
    target_fasta_path,
]
subprocess.run(" ".join(cmd), shell=True, check=True)
print("targetDB (fasta):", target_fasta_path)


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
prot_match_fig = go.Figure(layout=layout)
prot_match_fig.update_layout(
    title_text="<b>Number of proteome matches\nvs. number of peptides</b>",
    xaxis_title="Number of predicted peptides",
    yaxis_title="Number of matches",
    xaxis_range=None,
    yaxis_range=None,
) # plot number of matches versus number of predictions? (above some score value?)
pep_fig = go.Figure(layout=layout)
pep_fig.update_layout(title_text="<b>Peptide precision & coverage</b>")
aa_fig = go.Figure(layout=layout)
aa_fig.update_layout(title_text="<b>AA precision & coverage</b>")


MMSEQS2_ARGS = [
    "--seed-sub-mat VTML40.out",
    "--comp-bias-corr 0 --mask 0",
    "--spaced-kmer-mode 0",
    "-k 5",
]
output_metrics = {}
for output_file in os.listdir(args.output_dir):
    algo_name = output_file.split("_")[0]
    print("EVALUATE", algo_name)

    # Load tool predictions, match with ground truth
    output_path = os.path.join(args.output_dir, output_file)
    output_data = load_predictions(output_path, sequences_true)

    # Get idxs of GT labeled peptides & sequenced peptides (in correct output format)
    print(algo_name, output_data["score"].isnull().sum())
    output_data = output_data.sort_values("score", ascending=False)
    labeled_idx = output_data["sequence_true"].notnull()  
    sequenced_idx = get_sequenced_idx(output_data)

    # Prepare output sequences for metrics calculation
    output_data.loc[~sequenced_idx, "sequence"] = ""
    output_data.loc[~sequenced_idx, "aa_scores"] = ""
    output_data.loc[sequenced_idx, "sequence"] = output_data.loc[sequenced_idx, "sequence"].apply(
        ptms_to_delta_mass
    )
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
    create_query_fasta(output_data.loc[sequenced_idx], query_fasta_path)
    
    # Run mmseqs search
    search_df = run_mmseqs(
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
    
    # [Debug] display number of peptide matches & proteome matches
    print("N GT peptide matches:", output_data["pep_match"].sum())
    print("N proteome matches:", n_proteome_matches)
    idx = output_data["pep_match"] & ~output_data["proteome_match"]
    print("GT peptide matches w/o proteome match:", idx.sum())    
    idx = ~output_data["pep_match"] & output_data["proteome_match"]
    print("Proteome matches w/o GT peptide match:", idx.sum())
    print("\n", "=" * 100, "\n")


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


# Clean tmp folders
shutil.rmtree(search_tmp_dir)
