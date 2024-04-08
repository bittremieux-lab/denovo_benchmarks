# python 2.7
"""Script to convert predicted labels from the original algorithm 
output format to the common data format."""

import argparse
import pandas as pd

REPLACEMENTS = [
    ("Cmod", "C")  # C is written without Carbamidomethyl modification
]


def convert_to_output_format(sequence):
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

    # remove separator
    sequence = sequence.replace(",", "")

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    return sequence


parser = argparse.ArgumentParser()
parser.add_argument(
    "output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# read predictions from output file
output_data = pd.read_csv(args.output_path, sep="\t")

output_data = output_data.rename(
    {
        "output_seq": "sequence",
        "output_score": "score",
        "scan": "scans",
    },
    axis=1,
)
output_data = output_data[output_data["sequence"].notnull()]
output_data["sequence"] = output_data["sequence"].apply(
    convert_to_output_format
)

# save processed predictions to the same file
output_data.to_csv("outputs.csv", index=False)
