"""TODO"""

import os
import re
import numpy as np
import pandas as pd
from pyteomics import mgf
from tqdm import tqdm

from . import utils
from .metrics import spectral_angle
from token_masses import AA_MASSES
 

# import spectrum_utils.spectrum as sus
from koinapy import Koina
from sklearn.linear_model import LinearRegression


N_CALIBRATION_PSMS = 2000
WINDOW_SIZE = 2000
FRAGMENT_MASS_TOL = 0.02 # Da

# Methods to predict spectrum intensity and RT with Koina models

# - deeplc_hela_hf
# Modifications for this model include carbamidomethyl, oxidation of methionine, and N-terminal acetylatio (in train?). 
# Supported are all modification from UNIMOD. 
# Cysteine residues are assumed to be carbamidomethylated (C == C[UNIMOD:4]). 
# TODO: so how do we predict for immunopeptides then?

# - ms2pip_HCD2021
# Valid sequences lengths up to 30 AA
# There are no limitations to valid Precursor Charges.
# The model was trained on peptides with oxidation of methionine and fixed or variable carbamidomethylation of cysteine. 
# Supported are all modification from UNIMOD but modifications are only used to shift fragment mz they do not affect intensity.

model_rt = Koina("Deeplc_hela_hf", "koina.wilhelmlab.org:443")
supported_mods_rt = ["[UNIMOD:4]", "[UNIMOD:35]", "[UNIMOD:1]"] # deeplc, not used currently

models_I = {
    # TMT data
    "TMT": Koina("ms2pip_CID_TMT", "koina.wilhelmlab.org:443"),
    # TimsTOF data
    "TimsTOF": Koina("ms2pip_timsTOF2024", "koina.wilhelmlab.org:443"),
    # TOF data
    "TOF": Koina("ms2pip_TTOF5600", "koina.wilhelmlab.org:443"),
    # others (Orbitrap HCD data)
    "default": Koina("ms2pip_HCD2021", "koina.wilhelmlab.org:443"),
}
supported_mods_I = None # ms2pip 

def map_mods_delta_mass_to_unimod(sequence):
    """
    Map modifications in the sequence from delta mass format to UNIMOD format.
    Only modifications from SUPPORTED_MODS are mapped to their UNIMOD accession.
    Other modifications are replaced with generic [UNIMOD:] tag to be later filtered out.
    """
    SUPPORTED_MODS = [
        ("[+57.022]", "[UNIMOD:4]"),
        ("[+15.995]", "[UNIMOD:35]"),
        ("[+42.011]", "[UNIMOD:1]"),
    ]

    for replacement in SUPPORTED_MODS:
        sequence = sequence.replace(*replacement)
        
    sequence = re.sub("\[[^A-z]+\]", "[UNIMOD:]", sequence)
    return sequence

def check_supported_by_model(peptide, charge, supported_mods=None, min_seq_len=6):
    if not isinstance(peptide, str):
        return False

    if not (1 <= charge):
        return False

    # If supported_mods is given, check that peptide contains only these mods
    if supported_mods is not None:
        # Replace mods - only to test for presence of unsupported mods
        mod_replacements = {mod: mod.strip("[]").split(":")[-1] for mod in supported_mods}
        for mod, replacement in mod_replacements.items():
            peptide = peptide.replace(mod, replacement)
        
        if "UNIMOD:" in peptide:
            return False
    
    seq_len = len(peptide.replace("[^A-Z]", ""))
    if not (min_seq_len <= seq_len <= 30):
        return False
    
    return True

def predict_intensities(
    data: pd.DataFrame, 
    # TODO: should be DatasetTag, but we have a relative import problem!
    dataset_tags: tuple[str] = (), 
) -> tuple[pd.DataFrame, pd.DataFrame]:
    inputs = data[["sequence", "charge"]]
    inputs.columns = ["peptide_sequences", "precursor_charges"]

    # get model according to dataset_tags
    if "tmt" in dataset_tags:
        model_I = models_I["TMT"]
    elif "timstof" in dataset_tags:
        model_I = models_I["TimsTOF"]
    elif "sciex" in dataset_tags or "agilent" in dataset_tags:
        model_I = models_I["TOF"]
    else:
        model_I = models_I["default"]

    predictions = model_I.predict({col: inputs[col].values for col in inputs})
    predictions_I = pd.DataFrame(predictions["intensities"], index=inputs.index)
    predictions_mz = pd.DataFrame(predictions["mz"], index=inputs.index)
    return predictions_mz, predictions_I

def predict_RT(data: pd.DataFrame) -> np.array:
    inputs = data[["sequence"]]
    inputs.columns = ["peptide_sequences"]
    predictions_rt = model_rt.predict({col: inputs[col].values for col in inputs})
    return predictions_rt["irt"][:, 0]

def get_calibration_model(calib_psms):
    """
    Train RT calibration model to map true_RT to iRT units.
    """
    rt_calib_reg = LinearRegression()
    calib_rt_true = calib_psms["true_RT"].values
    calib_rt_pred = calib_psms["pred_RT"].values
    rt_calib_reg.fit(calib_rt_true[:, None], calib_rt_pred[:, None])
    print("DEBUG RT calibration:")
    print("Train:", rt_calib_reg.score(calib_rt_true[:, None], calib_rt_pred[:, None]))
    return rt_calib_reg

def calculate_SA(spec_idxs, predictions_mz, predictions_I, mgf_path):
    """
    Load mgf file, iterate through experimental spectra, calculate SA.

    Parameters
    ----------
    spec_idxs : pd.Series or dict
        spec_idx: index - psm idx in dataframe, value - 0-based spectrum idx in mgf file
    predictions_mz : pd.DataFrame
        DataFrame with predicted m/z arrays, indexed by psm idx
    predictions_I : pd.DataFrame
        DataFrame with predicted intensity arrays, indexed by psm idx
    mgf_path : str
        Path to the mgf file
    """
    with mgf.IndexedMGF(mgf_path) as spectra:
        SA_metrics = []
        for (psm_idx, spec_idx) in tqdm(spec_idxs.items(), total=len(spec_idxs)):
            spectrum = spectra[spec_idx]
            mz_true, I_true = utils.get_mz_I(spectrum)
            prec_mass = utils.get_precursor_mass(spectrum)
            charge = utils.get_charge(spectrum)

            mz_true, I_true = utils.remove_precursor_peak(
                np.array(mz_true), np.array(I_true), prec_mass, charge, FRAGMENT_MASS_TOL, isotope=3,
            )
            mz_pred, I_pred = predictions_mz.loc[psm_idx].values, predictions_I.loc[psm_idx].values

            SA = spectral_angle(mz_pred, I_pred, mz_true, I_true, FRAGMENT_MASS_TOL)
            SA_metrics.append(SA)
    return np.array(SA_metrics)
