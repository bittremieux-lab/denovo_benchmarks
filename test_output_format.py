import re
import pandas as pd
from pyteomics import proforma


def get_n_tokens(sequence: str) -> int:
    seq = proforma.parse(sequence)
    n_tokens = len(seq[0])
    if seq[1]["n_term"]:
        n_tokens += len(seq[1]["n_term"])
    if seq[1]["c_term"]:
        n_tokens += len(seq[1]["c_term"])
    return n_tokens

def validate_scan_index(scan_index: str) -> bool:
    SCAN_IDX_PATTERN = r"F\d+:\d+"
    return bool(re.fullmatch(SCAN_IDX_PATTERN, scan_index))

def validate_scan(scan: str) -> bool:
    SCAN_IDX_PATTERN = r"F\d+:\d+"
    return bool(re.fullmatch(SCAN_IDX_PATTERN, scan))

def validate_sequence(sequence: str) -> bool:
    try:
        seq = proforma.parse(sequence)
    except:
        return False
    return True

def validate_token_scores(scores: str, sequence: str) -> bool:
    n_tokens = get_n_tokens(sequence)
    return len(scores.split(",")) == n_tokens


output_data = pd.read_csv("test_outputs/test_output.csv")

for col in ["sequence", "score", "aa_scores"]:
    assert col in output_data, f"`{col}` must be presented."

assert (
    "scans" in output_data or "scan_indices" in output_data
), """
    At least one of {`scans`, `scan_indices`} columns 
    with the reference scan information must be presented.
"""

if "scans" in output_data:
    assert (
        output_data["scans"].apply(validate_scan).all()
    ), "`scans` do not have expected format."

if "scan_indices" in output_data:
    assert (
        output_data["scan_indices"].apply(validate_scan_index).all()
    ), "`scan_indices` do not have expected format."

assert (
    output_data["sequence"].apply(validate_sequence).all()
), """
    Predicted sequences are not in the output format.
"""

assert output_data.apply(
    lambda row: validate_token_scores(row["aa_scores"], row["sequence"]),
    axis=1,
).all(), """
    Number of per-token scores (','-separated scores in `aa_scores`) 
    must match number of individual tokens in a predicted sequence.
"""
