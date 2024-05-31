import json
import os
import ppx 
import subprocess
import pandas as pd
from tqdm import tqdm
from pyteomics import fasta, mgf
from oktoberfest.runner import run_job

VSC_DATA = "/data/antwerpen/209/vsc20960/"
VSC_SCRATCH = "/scratch/antwerpen/209/vsc20960/"
ROOT = os.path.join(VSC_SCRATCH, "benchmarking")

# Dataset id in MS repository 
DSET_ID = "MSV000086665" # "PXD004280"
# Number of raw files to download
N_FILES = 3
# Path to ThermoRawFileParser apptainer container
RAW_FILE_PARSER_PATH = os.path.join(VSC_SCRATCH, "benchmarking", "thermorawfileparser_latest.sif")
# Path to target database (fasta) for DB search  
DB_FILE = "uniprotkb_proteome_UP000005640_2024_05_16.fasta"
# Path to MSFragger executable file (jar)
MSFRAGGER_PATH = os.path.join(VSC_DATA, "easybuild", "build", "MSFragger-4.0", "MSFragger-4.0.jar")
# DB search params file (in SEARCH_PARAMS_DIR)
SEARCH_PARAMS_FILE = "human_closed_fragger.params"
# Max q-value threshold for rescored DB search results
Q_VAL_THRESHOLD = 0.01
# Spectrum params order for saving labeled mgf files
MGF_KEY_ORDER = ["title", "pepmass", "rtinseconds", "charge", "scans", "seq"]

PROTEOMES_DIR = os.path.join(ROOT, "proteomes")
SEARCH_PARAMS_DIR = os.path.join(ROOT, "search_params")
RESCORE_PARAMS_DIR = os.path.join(ROOT, "rescore_params")

RAW_DATA_DIR = os.path.join(ROOT, "raw")
MZML_DATA_DIR = os.path.join(ROOT, "mzml")
RESCORED_DATA_DIR = os.path.join(ROOT, "rescored")
MGF_DATA_DIR = os.path.join(ROOT, "mgf")

DATASET_STORAGE_DIR = os.path.join(VSC_DATA, "benchmarking", "datasets")

for data_dir in [
    RAW_DATA_DIR, 
    MZML_DATA_DIR, 
    RESCORED_DATA_DIR, 
    MGF_DATA_DIR, 
    RESCORE_PARAMS_DIR, 
    DATASET_STORAGE_DIR
]:
    os.makedirs(data_dir, exist_ok=True)


def get_psm_scan_id(psm_id):
    return psm_id.split("-")[1]


# Download raw files from repository
dset_dir = os.path.join(RAW_DATA_DIR, DSET_ID)
proj = ppx.find_project(
    DSET_ID, 
    local=dset_dir,
)
fnames = [fname for fname in proj.remote_files() if fname.lower().endswith(".raw")][:N_FILES]
proj.download(fnames)

raw_file_pathes = [file_path for file_path in proj.local_files() if file_path.suffix.lower() == ".raw"]
print(raw_file_pathes)

# Convert raw files to mzml (ThermoRawFileParser)
mzml_files_dir = os.path.join(MZML_DATA_DIR, DSET_ID)
os.makedirs(mzml_files_dir, exist_ok=True)

for file_path in tqdm(raw_file_pathes):
    cmd = [
        "apptainer",
        "exec",
        "--cleanenv",
        RAW_FILE_PARSER_PATH,
        "ThermoRawFileParser.sh",
        "-i",
        str(file_path),
        "-o",
        mzml_files_dir,
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)
print(os.listdir(mzml_files_dir), "\n")

# Convert raw files to mgf (ThermoRawFileParser)
mgf_files_dir = os.path.join(MGF_DATA_DIR, DSET_ID)
os.makedirs(mgf_files_dir, exist_ok=True)

for file_path in tqdm(raw_file_pathes):
    cmd = [
        "apptainer",
        "exec",
        "--cleanenv",
        RAW_FILE_PARSER_PATH,
        "ThermoRawFileParser.sh",
        "-i",
        str(file_path),
        "-o",
        mgf_files_dir,
        "-f",
        "0"
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)
print(os.listdir(mgf_files_dir), "\n")

# Generate decoys for DB search
db_path = os.path.join(PROTEOMES_DIR, DB_FILE)

name, ext = DB_FILE.split(".")
name = name + "_w_decoys"
db_w_decoys_file = ".".join([name, ext])

db_w_decoys_path = os.path.join(PROTEOMES_DIR, db_w_decoys_file)
fasta.write_decoy_db(
    db_path, 
    db_w_decoys_path,
    mode='reverse',
    prefix='rev_', # prefix Percolator expects to see
)

# Run DB search on mzml files (MSFragger)
search_params_file = os.path.join(SEARCH_PARAMS_DIR, SEARCH_PARAMS_FILE)
mzml_files = os.path.join(mzml_files_dir, "*.mzML")

cmd = [
    "java",
    "-Xmx64G",
    "-jar",
    MSFRAGGER_PATH,
    search_params_file,
    mzml_files,
]
subprocess.run(" ".join(cmd), shell=True, check=True)

# Rescore DB search results (Percolator with Prosit features)
rescored_files_dir = os.path.join(RESCORED_DATA_DIR, DSET_ID)
os.makedirs(rescored_files_dir, exist_ok=True)

spectra = mzml_files_dir  # the location of the mzML file containing the measured spectra
spectra_type = "mzml"  # the format the spectra are provided in ("mzml", "RAW", "d")
search_results = mzml_files_dir  # the location of the search engine output
search_results_type = "MSFragger" # the name of the search engine that produced the search results
intensity_model = "Prosit_2020_intensity_CID"  # the model used for fragment intensity prediction
retention_time_model = "Prosit_2019_irt"  # the model used for retention time prediction
prediction_server = "koina.wilhelmlab.org:443" # the Koina server that provides access to the specified models
output_directory = rescored_files_dir # the output folder for everything Oktoberfest produces during rescoring

# TODO: redefine other rescoring params, depending on data/search type?
rescoring_config = {
    "type": "Rescoring",
    "inputs":{
        "search_results": search_results,
        "search_results_type": search_results_type,
        "spectra": spectra,
        "spectra_type": spectra_type
    },
    "output": output_directory,
    "models": {
        "intensity": intensity_model,
        "irt": retention_time_model
    },
    "prediction_server": prediction_server,
    "ssl": True,
    "numThreads": 1,
    "fdr_estimation_method": "percolator",
    "massTolerance": 20,
    "unitMassTolerance": "ppm"
}
rescoring_config_path = os.path.join(RESCORE_PARAMS_DIR, f"{DSET_ID}_rescoring_config.json")
with open(rescoring_config_path, 'w') as fp:
    json.dump(rescoring_config, fp)

run_job(rescoring_config_path)

# Load DB search + rescoring results
results_path = os.path.join(rescored_files_dir, "results", "percolator", "rescore.percolator.psms.txt")
results_df = pd.read_csv(results_path, sep="\t")
results_df = results_df[results_df["q-value"] < Q_VAL_THRESHOLD][["PSMId", "filename", "proteinIds", "q-value"]]
results_df["scan_id"] = results_df["PSMId"].apply(get_psm_scan_id)

# Create dir for labeled mgf files
labeled_mgf_files_dir = os.path.join(DATASET_STORAGE_DIR, DSET_ID)
os.makedirs(labeled_mgf_files_dir, exist_ok=True)

# Annotated spectra with found PSMS, save as a dataset in DATASET_STORAGE_DIR
for mgf_file in os.listdir(mgf_files_dir):
    fname = mgf_file.split(".")[0]
    print(fname)
    
    file_labels_df = results_df[results_df["filename"] == fname]
    file_labels_df = file_labels_df.sort_values("scan_id", key=lambda x: x.apply(int))
    file_scan_ids = file_labels_df.scan_id.values.tolist()

    assert len(file_scan_ids) == len(set(file_scan_ids)), "Contains non-unique scan_ids."
    print("Number of PSMs:", len(file_scan_ids))

    file_labels_df = file_labels_df.set_index("scan_id")
    
    # Load original spectra (.mgf)
    unlabeled_mgf_path = os.path.join(mgf_files_dir, f"{fname}.mgf")
    spectra = mgf.IndexedMGF(unlabeled_mgf_path)
    print("Number of unlabeled spectra:", len(spectra))

    # Annotate spectra if possible, keep annotated only
    labeled_spectra = []
    for spectrum in tqdm(spectra):
        scan_id = spectrum["params"]["scans"]
        if scan_id in file_scan_ids:
            spectrum["params"]["seq"] = file_labels_df.loc[scan_id, "proteinIds"]
            labeled_spectra.append(spectrum)
    print("Number of labeled spectra:", len(labeled_spectra))
    
    # Write annotated spectra (to DATASET_STORAGE_DIR)
    labeled_mgf_path = os.path.join(labeled_mgf_files_dir, f"{fname}.mgf")
    mgf.write(
        labeled_spectra,
        labeled_mgf_path,
        key_order=MGF_KEY_ORDER,
    )
    print(f"{len(labeled_spectra)} spectra written to {labeled_mgf_path}.")
