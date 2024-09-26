import argparse
from pyteomics.mztab import MzTab
import pandas as pd
from glob import glob
import urllib.parse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--mztab_path", required=True, help="The path to the predictions mztab file made by casanovo."
)
parser.add_argument(
    "--mgf_in_dir", required=True, help="The dir with the input dataset .mgf files.",
)
parser.add_argument(
    "--mgf_out_dir", required=True, help="The dir with the output dataset .mgf files(i.e. in mgfs + initial sequence).",
)
args = parser.parse_args()

class IntermediateMapper():
    def __init__(self, output_path: str, input_dir: str, output_dir: str):
        # Read predictions from output file
        output_MzTab = MzTab(output_path)
        self.output_data = output_MzTab.spectrum_match_table[["spectra_ref", "sequence", "search_engine_score[1]"]]
        self.output_data['file_ref'] = self.output_data['spectra_ref'].apply(lambda x: x.split(":")[0])
        self.output_data['spectrum_idx'] = self.output_data['spectra_ref'].apply(lambda x: x.split(":index=")[1])
        self.input_dir = input_dir
        self.output_dir = output_dir

        # file://xyz%20abc.mgf => xyz abc.mgf
        file_parse = lambda path: urllib.parse.unquote(urllib.parse.urlparse(path).path)

        # Create a map from the spectra_ref to the file name
        meta_keys = list(output_MzTab.metadata.keys())
        self.id2file = {key.split("-")[0]: file_parse(output_MzTab.metadata[key]) for key in meta_keys if "ms_run[" in key}

        self.output_data['file_name'] = self.output_data['file_ref'].map(self.id2file).apply(lambda x: x.split("/")[-1])

        self.file_names = self.output_data['file_name'].unique()

        self.seq_mapping = self.output_data.groupby('file_name').apply(lambda x: x.set_index('spectrum_idx')['sequence'].to_dict()).to_dict()
        self.score_mapping = self.output_data.groupby('file_name').apply(lambda x: x.set_index('spectrum_idx')['search_engine_score[1]'].to_dict()).to_dict()

    def write_initial_seq_to_mfg(self):
        for input_mgf_path in glob(self.input_dir + "/*.mgf"):
            file_name = input_mgf_path.split("/")[-1]
            output_mgf_path = self.output_dir + "/" + file_name
            current_mapping = self.seq_mapping[file_name]
            current_score_mapping = self.score_mapping[file_name]
            
            with open(input_mgf_path, 'r') as infile, open(output_mgf_path, 'w+') as outfile:
                current_spectrum = []
                spectrum_idex = 0
                
                for line in infile:
                    line = line.strip()
                    
                    if line == "BEGIN IONS":
                        current_spectrum = [line]
                    elif line == "END IONS":
                        if str(spectrum_idex) in current_mapping:
                            sequence = current_mapping[str(spectrum_idex)]
                            current_spectrum.insert(1, f"SEQ={sequence}")
                            current_spectrum.insert(2, f"SCORE={current_score_mapping[str(spectrum_idex)]}")
                        else:
                            current_spectrum.insert(1, f"SEQ=MISSING")
                            current_spectrum.insert(2, f"SCORE={0}")
                            with open("missing_sequences.csv", "a+") as missing_file:
                                missing_file.write(f"{file_name},{spectrum_idex}\n")
                        
                        outfile.write("\n".join(current_spectrum) + "\n\n")

                        current_spectrum = []
                        spectrum_idex += 1
                    else:
                        current_spectrum.append(line)

intermediate_mapper = IntermediateMapper(args.mztab_path, args.mgf_in_dir, args.mgf_out_dir)
intermediate_mapper.write_initial_seq_to_mfg()

# python3 intermediate_mapper.py --mztab_path outputs.mztab --mgf_in_dir ../../sample_data/SMALL/mgf --mgf_out_dir ./tmp