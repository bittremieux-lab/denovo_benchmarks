"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""

import argparse
import re
import pandas as pd
from base import OutputMapperBase


class OutputMapper(OutputMapperBase):
    # Split before each AA and at "-" to separate terminal modifications
    # (Only works if PTM notation does not contain capital letters,
    # e.g. works with A[+1.234], but not with A[UNIMOD:X])
    PEP_SPLIT_PATTERN = r"(?<=[^-])(?=[A-Z])|-"

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

        # Fix for the case when multiple N-term modifications 
        # are predicted in a row. Removes redundant dashes. 
        sequence = sequence.replace("]-[", "][")

        aa_scores = aa_scores.split("|")

        # FIXME: Temporary fix for len(aa_scores) > len(sequence)
        # (probable reason: sequence tokenization pattern does not work
        # correctly with ProFroma Unimod notation for PTMs - [UNIMOD:X];
        # number of aa_scores is too large for sequences with PTMs). 
        # This aa_scores transformation is based on denovo/model/predict_step 
        # implementation logic, but can be incorrect or not 100% precise.
        n_tokens = len(
            re.split(self.PEP_SPLIT_PATTERN, sequence.replace("UNIMOD:", ""))
        )
        if len(aa_scores) != n_tokens:
            print(
                f"For sequence {sequence} expected {n_tokens} per-token scores, but got {len(aa_scores)} scores."
            )
            aa_scores = aa_scores[::-1][:n_tokens]
            aa_scores = aa_scores[::-1]

        aa_scores = ",".join(aa_scores)
        return sequence, aa_scores

    # TODO: add recalculation of peptide score (from corrected aa_scores) 
    

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
