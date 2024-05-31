"""Methods to convert ground truth labels to the common data format."""

import re
from pyteomics.mass.unimod import Unimod


REPLACEMENTS = [
    ("C[UNIMOD:4]", "C"), # C is written without Carbamidomethyl modification
]
PTM_PATTERN = r"([A-Z])\[UNIMOD:([0-9]+)\]"

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

    ptm = str(UNIMOD_DB.get(ptm_id).monoisotopic_mass)
    if not ptm.startswith("-"):
        ptm = "+" + ptm
    return f"{aa}[{ptm}]"

def format_sequence(sequence: str) -> str:
    """
    Convert peptide sequence to the common output data format.

    Parameters
    ----------
    sequence : str
        Peptide sequence in the original ground truth format.

    Returns
    -------
    transformed_sequence : str
        Peptide sequence in the common output data format.
    """

    # direct (token-to-token) replacements (if any)
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    # transformation of PTM notation
    # AA[ptm_mass] -> AA+ptm_mass
    sequence = re.sub(PTM_PATTERN, transform_match_ptm, sequence)

    return sequence
