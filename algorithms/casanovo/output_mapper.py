import argparse
import re
from pyteomics.mztab import MzTab


REPLACEMENTS = [("C+57.021", "C")]  # C is written without Carbamidomethyl modification

PTM_PATTERN = r"^([0-9.+-]+)([A-Z])"


def transform_match(match: str) -> str:
    ptm, aa = match.group(1), match.group(2)
    return aa + ptm


def convert_to_output_format(sequence):

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    # transformation of PTM notation
    # move N-terminus modifications BEYOND 1st AA
    sequence = re.sub(PTM_PATTERN, transform_match, sequence)

    return sequence


parser = argparse.ArgumentParser()
parser.add_argument("output_path", help="")
args = parser.parse_args()

# read predictions from output file
output_data = MzTab(args.output_path)

output_data = output_data.spectrum_match_table
output_data["sequence"] = output_data["sequence"].apply(convert_to_output_format)

# save processed predictions to the same file
output_data.to_csv(args.output_path.replace(".mztab", ".csv"), index=False)
