"""Shared data loading for ranking visualizations.

Reads all peptide_precision_plot_data.csv files, extracts AUC per
algorithm (best version) per dataset, and computes ranks.
"""

import csv
import glob
import os

import numpy as np
import pandas as pd


def load_auc_data(results_dir="results"):
    """Load AUC values: returns DataFrame with columns [dataset, algorithm, auc]."""
    rows = []
    csv_files = sorted(glob.glob(os.path.join(results_dir, "*/peptide_precision_plot_data.csv")))

    for path in csv_files:
        dataset = os.path.basename(os.path.dirname(path))
        with open(path) as f:
            for row in csv.DictReader(f):
                rows.append({
                    "dataset": dataset,
                    "algorithm": row["algorithm"],
                    "version": row["version"],
                    "auc": float(row["auc"]),
                })

    df = pd.DataFrame(rows)
    # Keep best version per algorithm per dataset
    df = df.loc[df.groupby(["dataset", "algorithm"])["auc"].idxmax()]
    return df[["dataset", "algorithm", "auc"]].reset_index(drop=True)


def compute_ranks(df):
    """Add a 'rank' column (1 = best) per dataset. Returns the augmented DataFrame."""
    df = df.copy()
    df["rank"] = df.groupby("dataset")["auc"].rank(ascending=False, method="min").astype(int)
    return df


def algo_order_by_median_rank(df):
    """Return list of algorithms sorted by median rank (best first)."""
    median_ranks = df.groupby("algorithm")["rank"].median().sort_values()
    return median_ranks.index.tolist()
