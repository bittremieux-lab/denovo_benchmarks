# python 2.7

import argparse
import pandas as pd

REPLACEMENTS = [("Cmod", "C")]  # C is written without Carbamidomethyl modification


def convert_to_output_format(sequence):
    # remove separator
    sequence = sequence.replace(",", "")

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    return sequence


parser = argparse.ArgumentParser()
parser.add_argument("output_path", help="")
args = parser.parse_args()

# read predictions from output file
output_data = pd.read_csv(args.output_path, sep="\t")

output_data = output_data.rename({"output_seq": "sequence", "scan": "scans"}, axis=1)
output_data = output_data[output_data["sequence"].notnull()]
output_data["sequence"] = output_data["sequence"].apply(convert_to_output_format)

# save processed predictions to the same file # ? only use ["sequence", "scan"] ?
output_data.to_csv("outputs.csv", index=False)
