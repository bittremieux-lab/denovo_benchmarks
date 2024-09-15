# python 2.7
"""
Script to convert predictions from the algorithm output format 
to the common output format.
"""

import argparse
import os
import pandas as pd

class OutputMapper:
    REPLACEMENTS = [
        ("C(Cam)", "C[UNIMOD:4]"),
        ("K(Carbamyl)", "K[UNIMOD:5]"),
        ("M(O)", "M[UNIMOD:35]"),
        ("S(Phospho)", "S[UNIMOD:21]"),
        ("T(Phospho)", "T[UNIMOD:21]"),
        ("Y(Phospho)", "Y[UNIMOD:21]"),
        ("(N-term|Acetyl)", "[UNIMOD:1]-"),
        ("(N-term|Carbamyl)", "[UNIMOD:5]-")
    ]

    def __init__(self, output_dir: str) -> None:
        output_frame = None
        for fname in os.listdir(output_dir):
            if fname.endswith(".novorpt.csv"):
                mgf_file = fname.replace(".novorpt.csv", "")
                output_path = os.path.join(output_dir, fname)
                output_data = self.read_csv(output_path)
                output_data = output_data.rename(columns=lambda x: x.strip())
                # Rename columns to the expected column names if needed
                output_data = output_data.rename(
                    {
                        "peptide": "sequence",
                        "score": "score",
                        "scanNum": "spectrum_id",
                        "aaScore": "aa_scores"
                    },
                    axis=1,
                )
                output_data = output_data[["spectrum_id", "sequence", "score", "aa_scores"]]
                output_data["spectrum_id"] = output_data.spectrum_id.apply(lambda x: f"{mgf_file}:{x}")
                output_data["sequence"] = output_data.sequence.apply(lambda x: self.format_sequence(x))
                output_data["aa_scores"] = output_data.aa_scores.apply(lambda x: x.replace("-", ","))
                output_data = output_data.apply(lambda e: e.map(lambda x: x.strip() if isinstance(x, str) else x))
                if output_frame is None:
                    output_frame = output_data
                else:
                    output_frame = pd.concat([output_frame, output_data])

        output_frame.to_csv("outputs.csv", index=False)
        return

    def read_csv(self, csv_path):
        def get_last_comment_line(lines):
            # Find all lines that start with '#'
            comment_lines = [i for i, line in enumerate(lines) if line.startswith('#')]
            if comment_lines:
                return comment_lines[-1]  # Return the index of the last comment line
            return -1  # Return -1 if no comment lines are found

        # Read the CSV file
        with open(csv_path, 'r') as file:
            lines = file.readlines()

            # Find the last line index that starts with '#'
        last_comment_index = get_last_comment_line(lines)

        # Extract data starting from the line after the last comment line
        if last_comment_index != -1:
            relevant_lines = lines[last_comment_index:]  # Lines after last #
        else:
            relevant_lines = []

            # Join the relevant lines into a single string
        filtered_csv = ''.join(relevant_lines)

        # Read the data into a DataFrame if there are any relevant lines
        if filtered_csv.strip():
            from io import StringIO
            df = pd.read_csv(StringIO(filtered_csv))
        else:
            df = pd.DataFrame()  # Create an empty DataFrame if no relevant lines
        return df

    def format_sequence(self, sequence):
        for repl_args in self.REPLACEMENTS:
            sequence = sequence.replace(*repl_args)

        return sequence


parser = argparse.ArgumentParser()
parser.add_argument(
    "--output_dir", required=True, help="The path to the algorithm predictions file."
)
args = parser.parse_args()

output_mapper = OutputMapper(args.output_dir)

