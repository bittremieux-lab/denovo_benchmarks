"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""

import argparse
import os
import re
import numpy as np
from collections import OrderedDict
from tqdm import tqdm
from pyteomics import mgf
from pyteomics.mztab import MzTab
from base import OutputMapperBase


class OutputMapper(OutputMapperBase):
    REPLACEMENTS = [
        ("C+57.021", "C[UNIMOD:4]"),
        # Amino acid modifications.
        ("M+15.995", "M[UNIMOD:35]"),    # Met oxidation
        ("N+0.984", "N[UNIMOD:7]"),     # Asn deamidation
        ("Q+0.984", "Q[UNIMOD:7]"),     # Gln deamidation
        # N-terminal modifications.
        ("+42.011", "[UNIMOD:1]"),      # Acetylation
        ("+43.006", "[UNIMOD:5]"),      # Carbamylation
        ("-17.027", "[UNIMOD:385]"),     # NH3 loss
        # ("+43.006-17.027", "[UNIMOD:5][UNIMOD:385]"),   # Carbamylation and NH3 loss
    ]
    PEP_SPLIT_PATTERN = r"(?<=.)(?=[A-Z])"
    MOD_PATTERN = r"\[UNIMOD:[0-9]+\]"
    N_TERM_MOD_PATTERN = r"^((\[UNIMOD:[0-9]+\])+)" # find N-term modifications

    def __init__(self, metadata: OrderedDict) -> None:
        """TODO."""
        self.filename_mapping = {}
        for k, v in metadata.items():
            if "ms_run" in k:
                ms_run_key = k.split("-")[0]
                filename = v.split("/")[-1].split(".")[0]
                self.filename_mapping[ms_run_key] = filename
        return
    
    def _transform_match_n_term_mod(self, match: re.Match) -> str:
        """
        Transform representation of peptide substring matching
        the N-term modification pattern.
        `[n_mod]PEP` -> `[n_mod]-PEP`
        
        Parameters
        ----------
        match : re.Match
            Substring matching the N-term modification pattern.

        Returns
        -------
        transformed_match : str
            Transformed N-term modification pattern representation.
        """
        ptm = match.group(1)
        return f"{ptm}-"

    def _spectrum_id_to_filename_idx(self, spectrum_id: str) -> str:
        filename, spectrum_idx = spectrum_id.split(":")
        
        filename = self.filename_mapping[filename]
        spectrum_idx = spectrum_idx.split("=")[-1]      
        return filename, spectrum_idx
    
    def _parse_scores(self, scores: str) -> list[float]:
        """
        Convert per-token scores from a string of float scores 
        separated by ',' to a list of float numbers.
        """
        scores = scores.split(",")
        scores = list(map(float, scores))
        return scores

    def format_spectrum_id(self, spectrum_id: str) -> str:
        """
        TODO.
        Transform spectrum id generated by the algorithm to the common format
        `ms_run[i]:index=j` -> `filename:j`.
        """
        filename, spectrum_idx = self._spectrum_id_to_filename_idx(spectrum_id)
        spectrum_id = filename + ":" + spectrum_idx
        return spectrum_id

    def format_sequence_and_scores(self, sequence, aa_scores):
        """
        Convert peptide sequence to the common output data format
        (ProForma with modifications represented with 
        Unimod accession codes, e.g. M[UNIMOD:35])
        and modify per-token scores if needed.

        This method is only needed if per-token scores have to be modified 
        to correspond the transformed sequence in ProForma format.
        Otherwise use `format_sequence` method instead.

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original algorithm output format.
        aa_scores: str
            String of per-token scores for each token in the sequence.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the common output data format.
        transformed_aa_scores: str
            String of per-token scores corresponding to each token
            in the transformed sequence.
        """   
        # direct (token-to-token) replacements
        for repl_args in self.REPLACEMENTS:
            sequence = sequence.replace(*repl_args)

        # format sequence and scores for n-term modifications
        if re.search(self.N_TERM_MOD_PATTERN, sequence):
            # transform n-term modification notation
            # represent in ProForma delta mass notation [+n_term_mod]-PEP
            sequence = re.sub(self.N_TERM_MOD_PATTERN, self._transform_match_n_term_mod, sequence)
            
            # count non-terminal tokens
            # assume modification is always attached a token or is terminal;
            # then any non-terminal modification does not increase number of tokens
            n_non_term_tokens = len(re.split(
                self.PEP_SPLIT_PATTERN,
                re.sub(self.MOD_PATTERN, "", sequence).strip("-")
            ))
            # number of scores should match number of tokens (+1 for n-term mods token, if any).
            # merge together all scores for n-term modification tokens
            aa_scores = self._parse_scores(aa_scores)
            
            n_term_scores, non_term_scores = aa_scores[:-n_non_term_tokens], aa_scores[-n_non_term_tokens:]
            term_score = np.mean(n_term_scores)

            aa_scores = [term_score] + non_term_scores
            aa_scores = self._format_scores(aa_scores)
        
        return sequence, aa_scores


parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_path", required=True, help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# Read predictions from output file
output_data = MzTab(args.output_path)
metadata = output_data.metadata
output_data = output_data.spectrum_match_table

# Rename columns to the expected column names if needed
output_data = output_data.rename(
    {
        "search_engine_score[1]": "score",
        "spectra_ref": "spectrum_id",
        "opt_ms_run[1]_aa_scores": "aa_scores",
    },
    axis=1,
)

# Transform data to the common output format
output_mapper = OutputMapper(metadata=metadata)
output_data = output_mapper.format_output(output_data)

# Save processed predictions to outputs.csv
# (the expected name for the algorithm output file)
output_data.to_csv("outputs.csv", index=False)
