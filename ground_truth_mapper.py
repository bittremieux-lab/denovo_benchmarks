import re

REPLACEMENTS = [
    # ("C", "C+57.021") # Consider C having a fixed Carbamidomethyl modification
]

PTM_PATTERN = r"([A-Z])\[([0-9.-]+)\]"


def transform_match(match: str) -> str:
    aa, ptm = match.group(1), match.group(2)
    if not ptm.startswith("-"):
        ptm = "+" + ptm
    return aa + ptm


def convert_GT_to_output_format(sequence):
    # remove cleavage sites
    if re.match(r"[A-Z].*.[A-Z]", sequence) is not None:  # check not mandatory
        sequence = sequence[2:-2]

    # direct (token-to-token) replacements (if any)
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    # transformation of PTM notation
    # AA[ptm_mass] -> AA+ptm_mass
    sequence = re.sub(PTM_PATTERN, transform_match, sequence)

    return sequence
