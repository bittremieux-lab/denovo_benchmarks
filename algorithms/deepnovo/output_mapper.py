# python 2.7
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
        ("Cmod", "C"),  # C is written without Carbamidomethyl modification
        ("Mmod", "M[+15.995]"),  # M Oxidation in ProForma delta mass notation
        ("Nmod", "N[+0.984]"),  # N Deamidation in ProForma delta mass notation
        ("Qmod", "Q[+0.984]"),  # Q Deamidation in ProForma delta mass notation
    ]

    def format_sequence(self, sequence):
        """
        Convert peptide sequence to the common output data format 
        (ProForma with modifications represented 
        in the delta mass notation).

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original algorithm output format.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the common output data format.  
        """

        # remove separator
        sequence = sequence.replace(",", "")

        # direct (token-to-token) replacements
        for repl_args in self.REPLACEMENTS:
            sequence = sequence.replace(*repl_args)

        return sequence
    

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_path", required=True, help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# Read predictions from output file
output_data = pd.read_csv(args.output_path, sep="\t")
output_data = output_data[(output_data != output_data.columns).all(axis=1)]

# Rename columns to the expected column names if needed
output_data = output_data.rename(
    {
        "output_seq": "sequence",
        "output_score": "score",
        "scan": "spectrum_id",
    },
    axis=1,
)
output_data = output_data[output_data["sequence"].notnull()]

# Transform data to the common output format
output_mapper = OutputMapper()
output_data = output_mapper.format_output(output_data)

# Save processed predictions to outputs.csv
# (the expected name for the algorithm output file)
output_data.to_csv("outputs.csv", index=False)
