"""Methods to convert ground truth labels to the common data format."""

import re


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

    REPLACEMENTS = []

    PTM_PATTERN = r"([A-Z])\[([0-9.-]+)\]"

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
        # TODO: transformation with adding "+" is no longer required
        # since ptm masses are now compared as close enough FLOAT NUMBERS, 
        # not exactly matching STRINGS 
        aa, ptm = match.group(1), match.group(2)
        if not ptm.startswith("-"):
            ptm = "+" + ptm
        return aa + ptm
    

    # remove cleavage sites
    if re.match(r"[A-Z].*.[A-Z]", sequence) is not None:
        sequence = sequence[2:-2]

    # direct (token-to-token) replacements (if any)
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    # transformation of PTM notation
    # AA[ptm_mass] -> AA+ptm_mass
    sequence = re.sub(PTM_PATTERN, transform_match_ptm, sequence)

    return sequence
