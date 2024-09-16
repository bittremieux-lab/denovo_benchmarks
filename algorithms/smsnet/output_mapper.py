"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""

import argparse
import re
import pandas as pd
import numpy as np
from base import OutputMapperBase

class OutputMapper(OutputMapperBase):
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
        return sequence.replace('m', 'M[UNIMOD:35]').replace('s', 'S[UNIMOD:21]').replace('t', 'T[UNIMOD:21]').replace('y', 'Y[UNIMOD:21]').replace('C', 'C[UNIMOD:4]')

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# Read predictions from output file
output_data = pd.read_csv(args.output_path, sep="\t")

# Add new columns
output_data['score'] = [np.mean([float(x) for x in aa_scores.split(';')]) for aa_scores in output_data['Scores'].values]
output_data['aa_scores'] = [x.replace(';', ',') for x in output_data['Scores'].values]
output_data['spectrum_id'] = [output_data['MS File'].loc[x] + ':' + str(output_data['ScanNum'].loc[x]) for x in output_data.index]

# Rename columns to the expected column names if needed
output_data = output_data.rename(
    {
        "RawPrediction": "sequence",
    },
    axis=1,
)

# Transform data to the common output format
# Modify OutputMapper to customize arguments and transformation.
output_mapper = OutputMapper()
output_data = output_mapper.format_output(output_data)

# Save processed predictions to outputs.csv
# (the expected name for the algorithm output file)
output_data.to_csv("outputs.csv", index=False)
