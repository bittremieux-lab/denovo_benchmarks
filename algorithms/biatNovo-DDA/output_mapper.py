"""
Script to convert predictions from the algorithm output format
to the common output format.
"""

import argparse
import pandas as pd
from base import OutputMapperBase


class OutputMapper(OutputMapperBase):
    REPLACEMENTS = [
        ("C(Carbamidomethylation)", "C[UNIMOD:4]"),
        ("M(Oxidation)", "M[UNIMOD:35]"),
        ("N(Deamidation)", "N[UNIMOD:7]"),
        ("Q(Deamidation)", "Q[UNIMOD:7]"),
    ]

    def format_sequence(self, sequence):
        sequence = str(sequence).replace(",", "")

        for repl_args in self.REPLACEMENTS:
            sequence = sequence.replace(*repl_args)

        return sequence


# setup the argument parser
parser = argparse.ArgumentParser()
parser.add_argument("--output_path", help="The path to the algorithm predictions file.")
args = parser.parse_args()

# read the algorithm predictions
output_data = pd.read_csv(args.output_path, sep="\t")
# print(output_data.dtypes)
# convert the algorithm predictions to the common output format
output_mapper = OutputMapper()
output_data = output_data[(output_data != output_data.columns).all(axis=1)]
output_data = output_data.rename(
    {
        "predicted_sequence": "sequence",
        "predicted_score": "score",
        "feature_id": "spectrum_id",
        "predicted_position_score": "aa_scores",
    },
    axis=1,
)
output_data = output_mapper.format_output(output_data)

output_data.to_csv("outputs.csv", index=False)
