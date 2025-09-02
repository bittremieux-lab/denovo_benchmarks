"""
Script to convert input .mgf files from the common input format 
to the algorithm expected format.
"""

import argparse
import numpy as np
import os
import yaml
from pyteomics import mgf
from tqdm import tqdm
from base import InputMapperBase


class InputMapper(InputMapperBase):
    
    def format_input(self, spectrum, valid_charge=None):
        """
        Convert the spectrum (annotation sequence and params) to the
        input format expected by the algorithm.

        For pi-HelixNovo, check if the spectrum charge
        falls into the valid range. If not, clip it to the valid range
        (to avoid runtime errors) and empty the spectrum.

        Parameters
        ----------
        spectrum : dict
            Peptide sequence in the original format.
        valid_charge: list of int, optional
            List of valid charge states. If None, 
            any charge state is allowed.
            Must contain all the allowed charge states
            (not just boundaries).

        Returns
        -------
        transformed_spectrum : dict
            Peptide sequence in the algorithm input format.
        """
        
        if valid_charge is not None:
            charge = int(spectrum["params"]["charge"][0])
            if charge not in valid_charge:
                print(
                    f"Spectrum {spectrum['params']['title']} has invalid charge {charge}. "
                    f"Clipping to the valid range [{valid_charge.min()}, {valid_charge.max()}] "
                    "and emptying the spectrum."
                )
                charge = int(np.clip(charge, valid_charge.min(), valid_charge.max()))
                spectrum["params"]["charge"] = charge

                spectrum["m/z array"] = np.empty(0)
                spectrum["intensity array"] = np.empty(0)
                spectrum["charge array"] = np.ma.empty(0, dtype=np.int64, fill_value=0)

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
parser.add_argument(
    "--config_path",
    help="The path to pi-HelixNovo config.yaml file with charge information.",
)
args = parser.parse_args()

# Transform data to the algorithm input format.
# Modify InputMapper to customize arguments and transformation.
input_mapper = InputMapper()

with open(args.config_path) as f:
    config = yaml.safe_load(f)
    max_charge = config["max_charge"]
    valid_charge = np.arange(1, config["max_charge"] + 1)

spectra = mgf.read(args.input_path)
mapped_spectra = [
    input_mapper.format_input(spectra[i], valid_charge)
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
