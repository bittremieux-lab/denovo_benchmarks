# python 2.7
"""Script to convert input .mgf files from the original format 
to the format expected by the algorithm."""

import argparse
import os
import re
from pyteomics import mgf


REPLACEMENTS = [
    ("C", "C(+57.02)")
]  # C always has fixed Carbamidomethyl modification

PTM_PATTERN = r"([A-Z])\[([0-9.-]+)\]"

from collections import namedtuple

PTM = namedtuple("DeepNovoPTM", ["amino_acid", "ptm_mass", "representation"])

REPLACEMENTS = [
    ("C", "C(+57.02)")
]  # C always has fixed Carbamidomethyl modification

PTM_PATTERN = r"([A-Z])\[([0-9.-]+)\]"
SUPPORTED_PTMS = [
    PTM("M", 15.99, "M(+15.99)"),
    PTM("N", 0.98, "N(+.98)"),
    PTM("Q", 0.98, "Q(+.98)"),
]
PTM_MASS_TOL = 0.01


def equal_with_tolerance(a, b, tolerance=1e-9):
    """
    Check if two float numbers are equal within a specified tolerance.
    """
    return abs(a - b) <= tolerance


def transform_match(match):
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
    aa, ptm = match.group(1), match.group(2)

    # transform PTMs supported by DeepNovo to the model's expected representation
    for supported_ptm in SUPPORTED_PTMS:
        if aa == supported_ptm.amino_acid and equal_with_tolerance(
            float(ptm), supported_ptm.ptm_mass, PTM_MASS_TOL
        ):
            return supported_ptm.representation

    # transform other PTMs
    if not ptm.startswith("-"):
        ptm = "+" + ptm
    return "{}({})".format(aa, ptm)


def convert_sequence_to_input_format(sequence):
    """
    Convert peptide sequence to the algorithm input format.

    Parameters
    ----------
    sequence : str
        Peptide sequence in the original format.

    Returns
    -------
    transformed_sequence : str
        Peptide sequence in the algorithm input format.
    """

    # remove cleavage sites
    if (
        re.match(r"[A-Z-].*.[A-Z-]", sequence) is not None
    ):  # check is not mandatory
        sequence = sequence[2:-2]

    # transformation of PTM notation
    # AA[ptm_mass] -> AA(+ptm_mass)
    sequence = re.sub(PTM_PATTERN, transform_match, sequence)

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    return sequence


def convert_to_input_format(spectrum, file_i=0):
    """
    Convert the spectrum (annotation sequence and params) to the
    input format expected by the algorithm.

    Parameters
    ----------
    spectrum : dict
        Peptide sequence in the original format.
    file_i: int
        Number of .mgf file being processed. Used to ensure a unique
        scan_id for each spectrum.

    Returns
    -------
    transformed_spectrum : dict
        Peptide sequence in the algorithm input format.
    """
    spectrum["params"]["seq"] = convert_sequence_to_input_format(
        spectrum["params"]["seq"]
    )

    # add file id to scan data
    spectrum["params"]["scans"] = "F{}:{}".format(
        file_i, spectrum["params"]["scans"]
    )
    return spectrum


parser = argparse.ArgumentParser()
parser.add_argument(
    "input_dir",
    help="The directory containing the input .mgf files. All files will be used.",
)
parser.add_argument(
    "output_path",
    help="The path to write prepared input data in the format expected by the algorithm.",
)
args = parser.parse_args()

input_paths = [
    os.path.join(args.input_dir, x)
    for x in os.listdir(args.input_dir)
    if ".mgf" in x
]
input_paths = sorted(input_paths)

mapped_spectra = []
for file_i, input_path in enumerate(input_paths):
    spectra = mgf.read(input_path)
    mapped_spectra += [
        convert_to_input_format(spectra[i], file_i)
        for i in range(len(spectra))
    ]

mgf.write(
    mapped_spectra,
    args.output_path,
    key_order=["title", "pepmass", "charge", "scans", "rtinseconds"],
    file_mode="w",
)
print(
    "{} spectra written to {}.".format(len(mapped_spectra), args.output_path)
)
