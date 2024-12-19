"""
Script to convert input .mgf files from the common input format 
to the algorithm expected format.
"""

import argparse
import os
import re
from pyteomics import mgf
from pyteomics.mass.unimod import Unimod
from tqdm import tqdm
from collections import namedtuple
from base import InputMapperBase

class InputMapper(InputMapperBase):
    def __init__(self,):
        pass

    def format_input(self, spectrum):
        """
        Convert the spectrum (annotation sequence and params) to the
        input format expected by the algorithm.

        Parameters
        ----------
        spectrum : dict
            Peptide sequence in the original format.

        Returns
        -------
        transformed_spectrum : dict
            Peptide sequence in the algorithm input format.
        """
        # No format changes
        
        return spectrum

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

spectra = mgf.read(args.input_path)
mapped_spectra = [
    input_mapper.format_input(spectra[i])
    for i in tqdm(range(len(spectra)))
]

# Save spectra in the algorithm input format.
# Modify the .mgf key order if needed.
mgf.write(
    mapped_spectra,
    args.output_path,
    key_order=["title", "rtinseconds", "pepmass", "charge"],
    file_mode="w",
)
print(
    "{} spectra written to {}.".format(len(mapped_spectra), args.output_path)
)
