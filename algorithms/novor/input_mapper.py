"""
Script to convert input .mgf files from the common input format
to the algorithm expected format.
"""

import argparse


class InputMapper:
    def format_input(self, input_file, output_file):
        with open(input_file, 'r') as input_fd:
            with open(output_file, 'w') as output_fd:
                scan_idx = 0
                for line in input_fd:
                    output_fd.write(line)
                    if line.startswith("TITLE="):
                        output_fd.write(f"SCANS={scan_idx}\n")
                        scan_idx += 1
        return output_file


parser = argparse.ArgumentParser()
parser.add_argument(
    "--input_path",
    help="The path to the input .mgf file.",
)
parser.add_argument(
    "--output_path",
    help="The path to write prepared input data in the format expected by the algorithm.",
)
args = parser.parse_args()

# Transform data to the algorithm input format.
# Modify InputMapper to customize arguments and transformation.
input_mapper = InputMapper()
spectra = input_mapper.format_input(args.input_path, args.output_path)
