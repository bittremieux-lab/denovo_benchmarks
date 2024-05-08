# python 2.7
"""Script to convert predicted labels from the original algorithm 
output format to the common data format."""

import argparse
import re
import pandas as pd
# from .base.output_mapper import OutputMapperBase


class OutputMapper:
    REPLACEMENTS = [
        ("Cmod", "C")  # C is written without Carbamidomethyl modification
    ]
    SPLIT_SEQ_PATTERN = r"(?<=.)(?=[A-Z])"
    
    def _format_scores(self, scores):
        """
        Write a list of float per-token scores
        into a string of float scores separated by ','.
        """
        return ",".join(map(str, scores))
    
    def format_scan(self, scan):
        return scan

    def format_sequence(self, sequence):
        """
        Convert peptide sequence to the common output data format.

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
    
    def simulate_token_scores(self, pep_score, sequence):
        """
        TODO.
        (if per-token scores are not provided by the model,
        define proxy per-token scores from the peptide score)
        """
        scores = [str(pep_score),] * len(re.split(self.SPLIT_SEQ_PATTERN, sequence))
        return self._format_scores(scores)
    
    def format_output(self, output_data):
        """TODO."""
        
        if "aa_scores" in output_data:
            output_data[["sequence", "aa_scores"]] = output_data.apply(
                lambda row: self.format_sequence_and_scores(row["sequence"], row["aa_scores"]),
                axis=1,
                result_type="expand",
            )
        
        else:
            output_data["sequence"] = output_data["sequence"].apply(
                self.format_sequence,
            )
            output_data["aa_scores"] = output_data.apply(
                lambda row: self.simulate_token_scores(row["score"], row["sequence"]), 
                axis=1,
            )
            
        if "scans" in output_data:
            output_data["scans"] = output_data["scans"].apply(self.format_scan)

        if "scan_indices" in output_data:
            output_data["scan_indices"] = output_data["scan_indices"].apply(
                self.format_scan_index
            )

        return output_data


parser = argparse.ArgumentParser()
parser.add_argument(
    "output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# read predictions from output file
output_data = pd.read_csv(args.output_path, sep="\t")

output_data = output_data.rename(
    {
        "output_seq": "sequence",
        "output_score": "score",
        "scan": "scans",
    },
    axis=1,
)
output_data = output_data[output_data["sequence"].notnull()]

output_mapper = OutputMapper()
output_data = output_mapper.format_output(output_data)

# save processed predictions to the same file
output_data.to_csv("outputs.csv", index=False)
