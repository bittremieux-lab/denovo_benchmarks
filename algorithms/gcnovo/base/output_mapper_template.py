"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""

import argparse
import re
import pandas as pd
from base import OutputMapperBase


class OutputMapper(OutputMapperBase):
    pass
    # Redefine base class methods 
    # or implement new methods if needed.
    

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# Read predictions from output file
output_data = pd.read_csv(args.output_path, sep="\t")

# Rename columns to the expected column names if needed
output_data = output_data.rename(
    {
        # "output_sequence": "sequence",
        # "output_score": "score",
        # "output_spectrum_id": "spectrum_id",
        # "output_aa_scores": "aa_scores",
        # ...
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
