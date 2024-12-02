from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
use_lstm = True
import argparse
import csv
import re
from dataclasses import dataclass
# ==============================================================================
# FLAGS (options) for this app
# ==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--beam_size", type=int, default="5")
parser.add_argument("--knapsack_file", type=str, default="./models/knapsack.npy")
parser.add_argument("--forward_model", type=str, default="./models/forward_deepnovo.pth")
parser.add_argument("--backward_model", type=str, default="./models/backward_deepnovo.pth")
parser.add_argument("--init_model", type=str, default="./models/init_net.pth")
parser.add_argument("--spectrum", type=str, default="spectrum.mgf")
args = parser.parse_args()

denovo_mgf_file = args.spectrum

denovo_output_feature_file = 'features.csv'
denovo_spectrum_fw = open("spectrum.mgf", 'w')

@dataclass
class Feature:
    spec_id: str
    mz: str
    z: str
    rt_mean: str
    scan: str

    def to_list(self):
        return [self.spec_id, self.mz, self.z, self.rt_mean, self.scan, "0.0:1.0", "1.0"]


def transfer_mgf(old_mgf_file_name, output_feature_file_name, spectrum_fw=denovo_spectrum_fw):
    with open(old_mgf_file_name, 'r') as fr:
        with open(output_feature_file_name, 'w') as fw:
            writer = csv.writer(fw, delimiter=',')
            header = ["spec_group_id","m/z","z","rt_mean","scans","profile","feature area"]
            writer.writerow(header)
            flag = False
            for line in fr:
                if "BEGIN ION" in line:
                    flag = True
                    spectrum_fw.write(line)
                elif not flag:
                    spectrum_fw.write(line)
                elif line.startswith("TITLE="):
                    spectrum_fw.write(line)
                elif line.startswith("PEPMASS="):
                    mz = re.split("=|\r|\n", line)[1]
                    spectrum_fw.write(line)
                elif line.startswith("CHARGE="):
                    z = re.split("=|\r|\n|\+", line)[1]
                    spectrum_fw.write("CHARGE=" + z + '\n')
                elif line.startswith("SCANS="):
                    scan = re.split("=|\r|\n", line)[1]
                    spectrum_fw.write(line)
                elif line.startswith("RTINSECONDS="):
                    rt_mean = re.split("=|\r|\n", line)[1]
                    spectrum_fw.write(line)
                elif line.startswith("SEQ="):
                    seq = re.split("=|\r|\n", line)[1]
                elif line.startswith("END IONS"):
                    feature = Feature(spec_id=scan, mz=mz, z=z, rt_mean=rt_mean, scan=scan)
                    writer.writerow(feature.to_list())
                    flag = False
                    del scan
                    del mz
                    del z
                    del rt_mean
                    spectrum_fw.write(line)
                else:
                    spectrum_fw.write(line)

transfer_mgf(denovo_mgf_file, denovo_output_feature_file, denovo_spectrum_fw)

denovo_spectrum_fw.close()
# ==============================================================================
# GLOBAL VARIABLES for VOCABULARY
# ==============================================================================
_PAD = "_PAD"
_GO = "_GO"
_EOS = "_EOS"
_START_VOCAB = [_PAD, _GO, _EOS]
PAD_ID = 0
GO_ID = 1
EOS_ID = 2
assert PAD_ID == 0
vocab_reverse = ['A',
                 'R',
                 'N',
                 'N[UNIMOD:7]',
                 'D',
                 #~ 'C',
                 'C[UNIMOD:4]',
                 'E',
                 'Q',
                 'Q[UNIMOD:7]',
                 'G',
                 'H',
                 'I',
                 'L',
                 'K',
                 'M',
                 'M[UNIMOD:35]',
                 'F',
                 'P',
                 'S',
                 'T',
                 'W',
                 'Y',
                 'V',
                ]

vocab_reverse = _START_VOCAB + vocab_reverse
vocab = dict([(x, y) for (y, x) in enumerate(vocab_reverse)])
vocab_size = len(vocab_reverse)
# ==============================================================================
# GLOBAL VARIABLES for THEORETICAL MASS
# ==============================================================================
mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949
mass_AA = {'_PAD': 0.0,
           '_GO': mass_N_terminus-mass_H,
           '_EOS': mass_C_terminus+mass_H,
           'A': 71.03711, # 0
           'R': 156.10111, # 1
           'N': 114.04293, # 2
           'N[UNIMOD:7]': 115.02695,
           'D': 115.02694, # 3
           #~ 'C(Carbamidomethylation)': 103.00919, # 4
           'C[UNIMOD:4]': 160.03065, # C(+57.02)
           #~ 'C(Carbamidomethylation)': 161.01919, # C(+58.01) # orbi
           'E': 129.04259, # 5
           'Q': 128.05858, # 6
           'Q[UNIMOD:7]': 129.0426,
           'G': 57.02146, # 7
           'H': 137.05891, # 8
           'I': 113.08406, # 9
           'L': 113.08406, # 10
           'K': 128.09496, # 11
           'M': 131.04049, # 12
           'M[UNIMOD:35]': 147.0354,
           'F': 147.06841, # 13
           'P': 97.05276, # 14
           'S': 87.03203, # 15
           'T': 101.04768, # 16
           'W': 186.07931, # 17
           'Y': 163.06333, # 18
           'V': 99.06841, # 19
          }

mass_ID = [mass_AA[vocab_reverse[x]] for x in range(vocab_size)]
mass_ID_np = np.array(mass_ID, dtype=np.float32)
mass_AA_min = mass_AA["G"] # 57.02146

# ==============================================================================
# GLOBAL VARIABLES for PRECISION, RESOLUTION, temp-Limits of MASS & LEN
# ==============================================================================
WINDOW_SIZE = 10 # 10 bins
MZ_MAX = 3000.0
MAX_NUM_PEAK = 500
KNAPSACK_AA_RESOLUTION = 10000 # 0.0001 Da
mass_AA_min_round = int(round(mass_AA_min * KNAPSACK_AA_RESOLUTION)) # 57.02146
KNAPSACK_MASS_PRECISION_TOLERANCE = 100 # 0.01 Da
num_position = 0
PRECURSOR_MASS_PRECISION_TOLERANCE = 0.01
AA_MATCH_PRECISION = 0.1
MAX_LEN = 50

# ==============================================================================
# HYPER-PARAMETERS of the NEURAL NETWORKS
# ==============================================================================
num_ion = 58
weight_decay = 0.0  # no weight decay lead to better result.
embedding_size = 512
num_lstm_layers = 1
num_units = 64
lstm_hidden_units = 512
dropout_rate = 0.25
batch_size = 1
num_workers = 6
steps_per_validation = 300
max_gradient_norm = 5.0


# ==============================================================================
# DATASETS
# ==============================================================================
knapsack_file = args.knapsack_file
topk_output = 1
denovo_input_spectrum_file = "spectrum.mgf"
denovo_input_feature_file = "features.csv"
denovo_output_file = "output.tsv" # # fixed for denovo_benchmarks
predicted_format = "deepnovo"
target_file = denovo_input_feature_file
predicted_file = denovo_output_file

# feature file column format
col_feature_id = "spec_group_id"
col_precursor_mz = "m/z"
col_precursor_charge = "z"
col_rt_mean = "rt_mean"
col_raw_sequence = "seq"
col_scan_list = "scans"
col_feature_area = "feature area"

# predicted file column format
pcol_feature_id = 0
pcol_feature_area = 1
pcol_sequence = 2
pcol_score = 3
pcol_position_score = 4
pcol_precursor_mz = 5
pcol_precursor_charge = 6
pcol_protein_id = 7
pcol_scan_list_middle = 8
pcol_scan_list_original = 9
pcol_score_max = 10


distance_scale_factor = 100.
sinusoid_base = 30000.
spectrum_reso = 10
n_position = int(MZ_MAX) * spectrum_reso

