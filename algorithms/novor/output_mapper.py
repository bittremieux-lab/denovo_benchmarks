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
        ("(N-term|Carbamyl)", "[UNIMOD:5]-"),
        ("(N-term|Acetyl)", "[UNIMOD:1]-"),
        ("K(Acetyl)", "K[UNIMOD:1]"),
        ("M(O)", "M[UNIMOD:35]"),
        ("Q(Pyro-Glu)", "Q[UNIMOD:28]"),
        ("E(Pyro-Glu)", "E[UNIMOD:27]"),
        ("K(Frm)", "K[UNIMOD:122]"),
        ("K(Frm)", "K[UNIMOD:122]"),
        ("S(Frm)", "S[UNIMOD:122]"),
        ("(N-term|Frm)", "[UNIMOD:122]-"),
        ("D(Methyl)", "D[UNIMOD:34]"),
        ("E(Methyl)", "E[UNIMOD:34]"),
        ("S(Phospho)", "S[UNIMOD:21]"),
        ("T(Phospho)", "T[UNIMOD:21]"),
        ("Y(Phospho)", "Y[UNIMOD:21]"),
        ("K(TMT6)", "K[UNIMOD:737]"),
        ("K(TMT10plex)", "K[UNIMOD:737]"),
        ("(N-term|TMT6)", "[UNIMOD:737]-"),
        ("(N-term|TMT10plex)", "[UNIMOD:737]-"),
        ("K(Lys4)", "K[UNIMOD:481]"),
        ("R(Arg6)", "R[UNIMOD:188]"),
        ("K(Lys8)", "K[UNIMOD:259]"),
        ("R(Arg10)", "R[UNIMOD:267]")
    ]

    def __init__(self, output_dir: str) -> None:
        output_frame = None
        for fname in os.listdir(output_dir):
            if fname.endswith(".novorai.csv"):
                mgf_file = fname.replace(".novorai.csv", "")
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
                output_data = output_data.apply(lambda e: e.map(lambda x: x.strip() if isinstance(x, str) else x))
                output_data['spectrum_id'] = output_data['spectrum_id'].astype(str)
                for index, row in output_data.iterrows():
                    spectrum_id = "%s:%s" % (mgf_file, row["spectrum_id"])

                    aa_scores = row["aa_scores"].strip().split("-")
                    if row['sequence'].startswith("(N-term"):
                        aa_scores.insert(0, aa_scores[0])
                    aa_scores = ",".join(aa_scores)

                    sequence = self.format_sequence(row["sequence"])

                    output_data.at[index, "spectrum_id"] = spectrum_id
                    output_data.at[index, "sequence"] = sequence
                    output_data.at[index, "aa_scores"] = aa_scores

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

