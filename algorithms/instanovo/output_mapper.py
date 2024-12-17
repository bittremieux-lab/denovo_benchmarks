"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""
import argparse
import pandas as pd
import numpy as np
from base import OutputMapperBase

import re

class OutputMapper(OutputMapperBase):
    REPLACEMENTS = {
        "C(+57.02)": "C[UNIMOD:4]",  # C Carboxyamidomethylation
        "N(+.98)":   "N[UNIMOD:7]",  # N Deamidation
        "Q(+.98)":   "Q[UNIMOD:7]",  # Q Deamidation
        "S(+79.97)": "S[UNIMOD:21]", # S Phosphorylation
        "T(+79.97)": "T[UNIMOD:21]", # T Phosphorylation
        "Y(+79.97)": "Y[UNIMOD:21]", # Y Phosphorylation
        "M(+15.99)": "M[UNIMOD:35]", # M Oxidation
         "(+42.01)": "[UNIMOD:1]",   # Acetylation
        "(+43.01)":  "[UNIMOD:5]",   # Carbamylation
        "(-17.03)":  "[UNIMOD:385]", # NH3 loss


    }

    def __init__(self):
        self.pattern = re.compile("|".join(map(re.escape, self.REPLACEMENTS.keys())))

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
        sequence = self.format_sequence(sequence)
        aa_scores = self.format_scores(aa_scores)
        return sequence, aa_scores
    
    def format_scores(self, scores):
        """
        Write a list of float per-token scores
        into a string of float scores separated by ','.
        """
        return ','.join(scores.strip('[]').split(', '))

    def format_sequence(self, sequence):
        """
        Convert peptide sequence to the common output data format 
        (ProForma with modifications represented with 
        Unimod accession codes, e.g. M[UNIMOD:35]).

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original algorithm output format.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the common output data format.  
        """
        sequence = str(sequence).replace(",", "")
        return self.pattern.sub(lambda m: self.REPLACEMENTS[m.group(0)], sequence)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_path", required=True, help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# Read predictions from output file
output_data = pd.read_csv(args.output_path)

# convert log probabilities to confidence score
output_data['score'] = output_data['log_probs'].apply(np.exp)

# Rename columns to the expected column names
output_data = output_data.rename(columns={'preds': 'sequence', 
                        'token_log_probs': 'aa_scores'})

# Select only the necessary columns
output_data = output_data[['sequence', 'score', 'aa_scores', 'spectrum_id']]


# Transform data to the common output format
output_mapper = OutputMapper()
output_data = output_mapper.format_output(output_data)

# Save processed predictions to outputs.csv
# (the expected name for the algorithm output file)
output_data.to_csv("outputs.csv", index=False)