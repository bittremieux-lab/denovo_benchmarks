"""Script to convert predicted labels from the original algorithm 
output format to the common data format."""

import argparse
import re
from pyteomics.mztab import MzTab


def format_sequence_and_scores(sequence: str, aa_scores: str) -> str:
    """
    Convert peptide sequence to the common output data format.
    If it changes the number of tokens in a sequence, 
    adjust per-token algorithm scores accordingly.

    Parameters
    ----------
    sequence : str
        Peptide sequence in the original algorithm output format.
    aa_scores: str
        Algorithm confidence scores for each token in sequence.
        Stored as a string of float scores separated by ','.

    Returns
    -------
    transformed_sequence : str
        Peptide sequence in the common output data format.
    """

    REPLACEMENTS = [
        ("C+57.021", "C")  # C is written without Carbamidomethyl modification
    ]

    PTM_PATTERN = r"^([0-9.+-]+)([A-Z])"  # TODO: move inside format_sequence_and_scores function?

    def transform_match_ptm(match: re.Match) -> str:
        """
        Transform representation of amino acids substring matching
        the PTM pattern.

        Parameters
        ----------
        match : re.Match
            Substring matching the PTM pattern.

        Returns
        -------
        transformed_match : str
            Transformed PTM pattern representation.
        """
        ptm, aa = match.group(1), match.group(2)
        return aa + ptm

    def parse_scores(scores: str) -> list[float]:
        """
        Convert per-token scores from a string of float scores 
        separated by ',' to a list of float numbers.
        """
        scores = scores.split(",")
        scores = list(map(float, scores))
        return scores

    def format_scores(scores: list[float]) -> str:
        """
        Write a list of float per-token scores
        into a string of float scores separated by ','.
        """
        return ",".join(map(str, scores))
    

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    # transformation of PTM notation
    # move N-terminus modifications BEYOND 1st AA (if any)
    # TODO: check & replacement can be not optimal!
    if re.search(PTM_PATTERN, sequence):
        sequence = re.sub(PTM_PATTERN, transform_match_ptm, sequence)
        # the N-terminus modification and the next AA will be considered as a single AA+PTM token,
        # so their scores should also be aggregated
        aa_scores = parse_scores(aa_scores)
        aa_scores[1] = (aa_scores[0] + aa_scores[1]) / 2
        aa_scores = aa_scores[1:]
        aa_scores = format_scores(aa_scores)

    return sequence, aa_scores


def format_scan_index(scan_index: str) -> str:
    """TODO."""

    FILE_IDX_PATTERN = "\[(\d+)\]"

    def transform_match_file_idx(match: re.Match) -> str:
        """TODO."""

        file_idx = int(match.group(0)[1:-1])
        return f"F{file_idx - 1}"

    scan_index = re.sub("[a-z=_]", "", scan_index)
    scan_index = re.sub(FILE_IDX_PATTERN, transform_match_file_idx, scan_index)
    return scan_index


parser = argparse.ArgumentParser()
parser.add_argument(
    "output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# read predictions from output file
output_data = MzTab(args.output_path)

output_data = output_data.spectrum_match_table
output_data = output_data.rename(
    {
        "search_engine_score[1]": "score",
        "spectra_ref": "scan_indices",
        "opt_ms_run[1]_aa_scores": "aa_scores",
    },
    axis=1,
)

output_data[["sequence", "aa_scores"]] = output_data.apply(
    lambda row: format_sequence_and_scores(row["sequence"], row["aa_scores"]),
    axis=1,
    result_type="expand",
)

output_data["scan_indices"] = output_data["scan_indices"].apply(
    format_scan_index
)

# save processed predictions to the same file
output_data.to_csv("outputs.csv", index=False)
