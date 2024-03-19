# python 2.7

import argparse
import os
import re
from pyteomics import mgf


REPLACEMENTS = [("C", "C(+57.02)")]  # C always has fixed Carbamidomethyl modification

PTM_PATTERN = r"([A-Z])\[([0-9.-]+)\]"


def transform_match(match):
    aa, ptm = match.group(1), match.group(2)
    if not ptm.startswith("-"):
        ptm = "+" + ptm
    return "{}({})".format(aa, ptm)


def convert_sequence_to_input_format(sequence):
    # remove cleavage sites
    if re.match(r"[A-Z-].*.[A-Z-]", sequence) is not None:  # check not mandatory
        sequence = sequence[2:-2]

    # transformation of PTM notation
    # AA[ptm_mass] -> AA(+ptm_mass)
    sequence = re.sub(PTM_PATTERN, transform_match, sequence)

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    return sequence


def convert_to_input_format(spectrum, file_i):  # TODO naming
    spectrum["params"]["seq"] = convert_sequence_to_input_format(
        spectrum["params"]["seq"]
    )

    # add file id to scan data
    spectrum["params"]["scans"] = "F{}:{}".format(file_i, spectrum["params"]["scans"])
    return spectrum


parser = argparse.ArgumentParser()
parser.add_argument(
    "input_dir", help=""
)  # change to input_path={input_path or input_dir}?
parser.add_argument("output_path", help="")
args = parser.parse_args()

input_paths = [
    os.path.join(args.input_dir, x) for x in os.listdir(args.input_dir) if ".mgf" in x
]
input_paths = sorted(input_paths)

mapped_spectra = []
for file_i, input_path in enumerate(input_paths):
    spectra = mgf.read(input_path)
    mapped_spectra += [
        convert_to_input_format(spectra[i], file_i) for i in range(len(spectra))
    ]

mgf.write(
    mapped_spectra,
    args.output_path,
    key_order=["title", "pepmass", "charge", "scans", "rtinseconds"],
    file_mode="w",
)
print("{} spectra written to {}.".format(len(mapped_spectra), args.output_path))
