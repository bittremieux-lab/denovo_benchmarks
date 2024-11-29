"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""

import argparse
import re
import pandas as pd
from base import OutputMapperBase


class OutputMapper(OutputMapperBase):
    REPLACEMENTS = [
        ("C", "C[UNIMOD:4]"), # C has Carbamidomethylation by default
    ]

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

        # direct (token-to-token) replacements
        for repl_args in self.REPLACEMENTS:
            sequence = sequence.replace(*repl_args)
        
        return sequence

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

        aa_scores = ",".join(aa_scores.strip("[]").split(" "))
        return sequence, aa_scores

    def format_spectrum_id(self, spectrum_id):
        """
        Represent spectrum spectrum id as {filename}:{index} string, 
        where
        - `filename` - name of the .mgf file in a dataset 
            (without .mgf extension)
        - `index` - index (0-based) of each spectrum in an .mgf file.
        """

        # Keep only the filename (without preceding path)
        spectrum_id = spectrum_id.split("/")[-1]

        # Remove .mgf extension
        spectrum_id = spectrum_id.replace(".mgf", "")
        return spectrum_id
    

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# Read predictions from output file
output_data = pd.read_csv(args.output_path, sep=",")

# Transform data to the common output format
# Modify OutputMapper to customize arguments and transformation.
output_mapper = OutputMapper()
output_data = output_mapper.format_output(output_data)

# Save processed predictions to outputs.csv
# (the expected name for the algorithm output file)
output_data.to_csv("outputs.csv", index=False)
