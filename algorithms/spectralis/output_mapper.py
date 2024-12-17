import argparse
from pyteomics.mztab import MzTab
import urllib.parse
from base import OutputMapperBase
import pandas as pd
from glob import glob
import re

parser = argparse.ArgumentParser()
parser.add_argument(
    "--spectralis_output_dir", required=True, help="The path to the algorithm predictions file."
)
parser.add_argument(
    "--mgf_in_dir", required=True, help="The dir with the input dataset .mgf files.",
)
parser.add_argument(
    "--casanovo_output", required=True, help="The path to casanovo's predictions file."
)
args = parser.parse_args()

class OutputMapper(OutputMapperBase):
    REPLACEMENTS = [
        ("C+57.021", "C[UNIMOD:4]"),
        # Amino acid modifications.
        ("M+15.995", "M[UNIMOD:35]"),    # Met oxidation
        ("N+0.984", "N[UNIMOD:7]"),     # Asn deamidation
        ("Q+0.984", "Q[UNIMOD:7]"),     # Gln deamidation
        # N-terminal modifications.
        ("+42.011", "[UNIMOD:1]"),      # Acetylation
        ("+43.006", "[UNIMOD:5]"),      # Carbamylation
        ("-17.027", "[UNIMOD:385]"),     # NH3 loss
        # "+43.006-17.027": 25.980265      # Carbamylation and NH3 loss
    ]
    N_TERM_MOD_PATTERN = r"^(\[UNIMOD:[0-9]+\])" # find N-term modifications

    def __init__(self, spectralis_output_dir: str, input_dir: str, casanovo_output: str):
        output_MzTab = MzTab(casanovo_output)
        self.output_data = output_MzTab.spectrum_match_table[["spectra_ref", "sequence", "search_engine_score[1]"]]
        self.output_data['file_ref'] = self.output_data['spectra_ref'].apply(lambda x: x.split(":")[0])
        self.output_data['spectrum_idx'] = self.output_data['spectra_ref'].apply(lambda x: x.split(":index=")[1])

        # file://xyz%20abc.mgf => xyz abc.mgf
        file_parse = lambda path: urllib.parse.unquote(urllib.parse.urlparse(path).path)

        # Create a map from the spectra_ref to the file name
        meta_keys = list(output_MzTab.metadata.keys())
        self.id2file = {key.split("-")[0]: file_parse(output_MzTab.metadata[key]) for key in meta_keys if "ms_run[" in key}

        self.output_data['file_name'] = self.output_data['file_ref'].map(self.id2file).apply(lambda x: (x.split("/")[-1])[:-4])
        self.output_data['spectrum_id'] = self.output_data['file_name'].astype(str) + ":" + self.output_data['spectrum_idx'].astype(str)

        file_names = [(x.split("/")[-1])[:-4] for x in glob(spectralis_output_dir + "/*.csv")]
        dfs = [pd.read_csv(f, header=0, index_col=None) for f in glob(spectralis_output_dir + "/*.csv")]
        for file_name, df in zip(file_names, dfs):
            df['file_name'] = file_name

        self.spectralis_data = pd.concat(dfs, ignore_index=True)
        self.spectralis_data['spectrum_id'] = self.spectralis_data['file_name'].astype(str) + ":" + self.spectralis_data['scans'].astype(str)

        self.output_data = self.output_data.merge(self.spectralis_data, on='spectrum_id', how='left')

        self.output_data.rename(columns={'Spectralis_score': 'score'}, inplace=True)
        self.output_data = self.output_data[['sequence', 'score', 'spectrum_id']]

    def _transform_match_n_term_mod(self, match: re.Match) -> str:
        """
        Transform representation of peptide substring matching
        the N-term modification pattern.
        `[n_mod]PEP` -> `[n_mod]-PEP`
        
        Parameters
        ----------
        match : re.Match
            Substring matching the N-term modification pattern.

        Returns
        -------
        transformed_match : str
            Transformed N-term modification pattern representation.
        """
        ptm = match.group(1)
        return f"{ptm}-"

    def format_sequence(self, sequence: str) -> str:
        """
        Convert peptide sequence to the common output data format 
        (ProForma with modifications represented 
        in the delta mass notation).

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original algorithm output format.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the common output data format.   
        """

        # direct (token-to-token) replacements
        for repl_args in self.REPLACEMENTS:
            sequence = sequence.replace(*repl_args)

        # transform n-term modification notation
        # represent in ProForma delta mass notation [+n_term_mod]-PEP
        if re.search(self.N_TERM_MOD_PATTERN, sequence):
            sequence = re.sub(self.N_TERM_MOD_PATTERN, self._transform_match_n_term_mod, sequence)

        return sequence



output_mapper = OutputMapper(args.spectralis_output_dir, args.mgf_in_dir, args.casanovo_output)
output_mapper.format_output(output_mapper.output_data)
output_mapper.output_data.to_csv("outputs.csv", index=False)
