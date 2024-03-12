import argparse
import os
import depthcharge
import pandas as pd
from pyteomics import mgf

from ground_truth_mapper import convert_GT_to_output_format
from metrics import aa_match_metrics, aa_match_batch


parser = argparse.ArgumentParser()
parser.add_argument("output_path", help="")
parser.add_argument("data_dir", help="")  # change to ground truth paths?
args = parser.parse_args()


input_paths = [
    os.path.join(args.data_dir, x) for x in os.listdir(args.data_dir) if ".mgf" in x
]
input_paths = sorted(input_paths)

# convert predictions to the common (GT) format
sequences_true = {}
n_seq = 0
for mgf_path in input_paths:
    spectra = mgf.read(mgf_path)
    sequences_true.update(
        {n_seq + i: spec["params"]["seq"] for i, spec in enumerate(spectra)}
    )
    n_seq = len(sequences_true)
sequences_true = pd.Series(sequences_true)
sequences_true = sequences_true.apply(convert_GT_to_output_format)

psms = pd.read_csv(args.output_path)
psms["sequence_pred"] = psms["sequence"]
psms["sequence"] = sequences_true.values
psms = psms[["sequence", "sequence_pred", "search_engine_score[1]"]]

# calculate metrics
aa_matches_batch, n_aa1, n_aa2 = aa_match_batch(
    psms["sequence"],
    psms["sequence_pred"],
    depthcharge.masses.PeptideMass("massivekb").masses,
)

aa_precision, aa_recall, pep_precision = aa_match_metrics(
    aa_matches_batch, n_aa1, n_aa2
)

# log results
with open("./metrics.txt", "w") as f:
    f.write(f"AA precision: {round(aa_precision, 3)}\n")
    f.write(f"AA recall: {round(aa_recall, 3)}\n")
    f.write(f"Pep precision: {round(pep_precision, 3)}\n")
