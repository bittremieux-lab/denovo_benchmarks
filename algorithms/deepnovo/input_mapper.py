# python 2.7
"""Script to convert input .mgf files from the original format 
to the format expected by the algorithm."""

import argparse
import os
import re
from pyteomics import mgf
from pyteomics.mass.unimod import Unimod
from tqdm import tqdm
from collections import namedtuple

PTM = namedtuple("DeepNovoPTM", ["amino_acid", "ptm_unimod_id", "representation"])

REPLACEMENTS = []

PTM_PATTERN = r"([A-Z])\[UNIMOD:([0-9]+)\]"
SUPPORTED_PTMS = [
    PTM("C", 4, "C(+57.02)"),
    PTM("M", 35, "M(+15.99)"),
    PTM("N", 7, "N(+.98)"),
    PTM("Q", 7, "Q(+.98)"),
]

UNIMOD_DB = Unimod()


def transform_match_ptm(match):
    """
    Transform representation of amino acids substring matching
    the PTM pattern.
    Expects PTMs in ProForma notation, e.g. 'M[UNIMOD:35]'.

    Parameters
    ----------
    match : re.Match
        Substring matching the PTM pattern.

    Returns
    -------
    transformed_match : str
        Transformed PTM pattern representation.
    """
    aa, ptm_id = match.group(1), int(match.group(2))
    
    # transform PTMs supported by DeepNovo to the model's expected representation
    for supported_ptm in SUPPORTED_PTMS:
        if aa == supported_ptm.amino_acid and ptm_id == supported_ptm.ptm_unimod_id:
            return supported_ptm.representation
    
    # transform other PTMs
    ptm = str(UNIMOD_DB.get(ptm_id).monoisotopic_mass)
    if not ptm.startswith("-"):
        ptm = "+" + ptm
    return "{}({})".format(aa, ptm)


def format_sequence(sequence):
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
        re.match(r"[A-Z-_].*.[A-Z-_]", sequence) is not None
    ):  # check is not mandatory
        sequence = sequence[2:-2]

    # transformation of PTM notation
    # AA[ptm_mass] -> AA(+ptm_mass)
    sequence = re.sub(PTM_PATTERN, transform_match_ptm, sequence)

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    return sequence


def format_input(spectrum, file_i=0):
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
    spectrum["params"]["seq"] = format_sequence(
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
# parser.add_argument(
#     "input_path",
#     help="The path to the input .mgf file.",
# )
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
        format_input(spectra[i], file_i)
        for i in tqdm(range(len(spectra)))
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
