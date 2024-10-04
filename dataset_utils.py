import json
import os
import ppx 
import re
import shutil
import subprocess
import pandas as pd
import alphatims.bruker
from tqdm import tqdm
from pyteomics import fasta, mgf

VSC_DATA =  os.environ['VSC_DATA']
VSC_SCRATCH = os.environ['VSC_DATA'] 
ROOT = os.environ['ROOT'] 
VSC_FRAGPIPE = os.environ['VSC_FRAGPIPE']

# Path to ThermoRawFileParser apptainer container
RAW_FILE_PARSER_PATH = os.path.join(VSC_SCRATCH, "benchmarking", "thermorawfileparser_latest.sif")
# Path to msconvert apptainer container
MSCONVERT_PATH = os.path.join(VSC_SCRATCH, "benchmarking", "pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif")
# Path to MSFragger DB split script (for running DB search with large DB)
SPLIT_SCRIPT_PATH = os.path.join(VSC_FRAGPIPE, "FragPipe/21.1-Java-11/tools/msfragger_pep_split.py")
# Path to MSFragger executable file (jar)
MSFRAGGER_PATH = os.path.join(VSC_DATA, "easybuild", "build", "MSFragger-4.0", "MSFragger-4.0.jar")
# Path to MSBooster executable file (jar)
MSBOOSTER_PATH = os.path.join(VSC_FRAGPIPE, "MSBooster", "1.2.31-Java-11", "MSBooster-1.2.31.jar")
DIANN_PATH = os.path.join(VSC_FRAGPIPE, "FragPipe/21.1-Java-11/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8")
KOINA_URL = "https://koina.wilhelmlab.org:443/v2/models/"
MSBOOSTER_BASE_PARAMS = os.path.join(VSC_SCRATCH, "benchmarking", "params_rescore", "msbooster_base.params")
# Spectrum params order for saving labeled mgf files
MGF_KEY_ORDER = ["title", "pepmass", "rtinseconds", "charge", "scans", "seq"]
# Path to the file with datasets properties (tags)
DATASET_TAGS_PATH = os.path.join(ROOT, "denovo_benchmarks", "dataset_tags.tsv")

PROTEOMES_DIR = os.path.join(ROOT, "proteomes")

RAW_DATA_DIR = os.path.join(ROOT, "raw")
MZML_DATA_DIR = os.path.join(ROOT, "mzml")
RESCORED_DATA_DIR = os.path.join(ROOT, "rescored")
MGF_DATA_DIR = os.path.join(ROOT, "mgf")

DATASET_STORAGE_DIR = os.path.join(VSC_DATA, "benchmarking", "datasets")

# Spectra smoothing for .d to .mgf conversion with alphatims
CENTROID_WINDOW = 5
MAX_SPECTRA_PER_FILE = 20000


def get_files_list(download_config):
    """
    Select files for dataset based on 
    selection rules defined in the download_config.
    
    TODO.
    """
    def check_file(file_path):
        for keyword in download_config.keywords:
            if keyword not in file_path:
                return False
        if not file_path.lower().endswith(ext):
            return False
        return True

    dset_id = download_config.dset_id
    dset_dir = os.path.join(RAW_DATA_DIR, dset_id)
    ext = download_config.ext

    if "PXD" in dset_id or "MSV" in dset_id:
        proj = ppx.find_project(
            dset_id, 
            local=dset_dir,
        )

        files_list = [
            file_path
            for file_path 
            in proj.remote_files() 
            if check_file(file_path)
        ]
    else:
        # if load via wget, file_path is just fname.ext
        files_list = [
            os.path.basename(file_link) 
            for file_link 
            in download_config.links
            # if check_file(file_link)
        ]

    files_list = files_list[:download_config.n_files]
    files_list = {
        os.path.basename(file_path)[:-len(ext)]: file_path
        for file_path
        in files_list
    }
    return files_list


def download_files(download_config, files_list):#, unpack_dir=None, raw_ext=None):
    # TODO: now RAW_DATA_DIR can also contain mzml and mgf -- if specified? 
    # (just all the source data for a DSET_ID always is stored at RAW_DATA_DIR?)
    dset_id = download_config.dset_id
    dset_dir = os.path.join(RAW_DATA_DIR, dset_id)
    print(f"Loading dataset {dset_id} to the folder {dset_dir}")

    if "PXD" in dset_id or "MSV" in dset_id:
        proj = ppx.find_project(
            dset_id, 
            local=dset_dir,
        )
        print("Local files:", proj.local_files())
    
        # select files to download
        # TODO: change: not skipping existing files so far = download all from fnames
        fnames = list(files_list.values())
        if download_config.ext == ".wiff":
            fnames += [fname + ".scan" for fname in fnames]
        proj.download(fnames)
        print("Loaded files:\n", proj.local_files()) # TODO: remove

    else:
        # if load via wget, need to take download_links from downloag_config for each file
        file_links = {
            os.path.basename(file_link)[:-len(download_config.ext)]: file_link
            for file_link 
            in download_config.links
        }
        print("Local files:", os.listdir(dset_dir))
        for fname, file_path in files_list.items():
            if not os.path.exists(os.path.join(dset_dir, file_path)):
                cmd = [
                    "wget",
                    "-P",
                    dset_dir,
                    file_links[fname]
                ]
                subprocess.run(" ".join(cmd), shell=True, check=True)
        print("Loaded files:", os.listdir(dset_dir))


def convert_raw(dset_id, files_list, target_dir, target_ext=".mzml"):
    os.makedirs(target_dir, exist_ok=True)    

    dset_dir = os.path.join(RAW_DATA_DIR, dset_id)
    raw_file_pathes = {
        fname: os.path.join(dset_dir, file_path)
        for fname, file_path
        in files_list.items()
    }
    print("Files:\n", list(raw_file_pathes.values()))
    print(f"Converting to {target_ext}. Storing to {target_dir}")

    ext_flag = "--mgf" if target_ext == ".mgf" else "--mzML"
    filter_ms_level = 2 if target_ext == ".mgf" else 1
    for fname, file_path in tqdm(raw_file_pathes.items()):
        # TODO в идеале тоже надо пропускать файлы, которые уже существуют
        out_fname = fname + target_ext
        cmd = [
            "apptainer",
            "exec",
            "--cleanenv",
            MSCONVERT_PATH,
            "wine msconvert",
            ext_flag,
            "-z",
            "-o",
            target_dir,
            "--outfile",
            out_fname,
            "--filter",
            '"peakPicking vendor"',
            "--filter",
            '"precursorRefine"',
            # "--filter",
            # '"chargeStatePredictor singleChargeFractionTIC=0.9 maxKnownCharge=4"',
            "--filter",
            f'"msLevel {filter_ms_level}-"',
            "--filter",
            '"titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>"'
            # '"titleMaker <RunId>.<Index>.<ChargeState>.<MsLevel>"'
        ]
        cmd += [file_path]
        subprocess.run(" ".join(cmd), shell=True, check=True)
    print(os.listdir(target_dir))
    
    
# should we run it only if there is no prepared decoys file? 
def generate_decoys_fasta(dset_name, db_file):
    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    db_path = os.path.join(PROTEOMES_DIR, db_file)

    cmd = [
        "cd",
        mzml_files_dir,
        "&&",
        "philosopher workspace --init",
        "&&",
        "philosopher database",
        "--custom",
        db_path,
        "--contam"
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)

    db_w_decoys_path = [fname for fname in os.listdir(mzml_files_dir) if fname.endswith(".fas")][0]
    db_w_decoys_path = os.path.join(mzml_files_dir, db_w_decoys_path)
    return db_w_decoys_path


def run_database_search(dset_name, db_w_decoys_path, db_search_config):
    search_ext = db_search_config.ext
    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    mzml_files = [f for f in os.listdir(mzml_files_dir) if os.path.splitext(f)[1].lower() == search_ext]
    mzml_files = [os.path.join(mzml_files_dir, f) for f in mzml_files]
    
    options = [
        "--database_name",
        db_w_decoys_path,
        "--decoy_prefix",
        "rev_",
        "--output_format",
        "pepxml_pin", # .pin outputs for MSBooster
    ]
    if db_search_config.ext == ".d":
        options += ["--write_uncalibrated_mgf", "1"]
    
    # Parse additional search params from config if provided
    for arg in [*db_search_config.search_params.items()]:
        options += list(map(str, arg))
    
    cmd = [
        "java",
        "-Xmx160G",
        "-jar",
        MSFRAGGER_PATH,
        *options,
        *mzml_files,
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)
    print("DB search results (.pepXML, .pin):\n", os.listdir(mzml_files_dir))

    if search_ext == ".d":
        # Use uncalibrared mzml files to get "proxy" mzml files
        for file_path in mzml_files:
            fname = os.path.splitext(file_path)[0]
            src_fname = fname + "_uncalibrated.mzML"
            dst_fname = fname + ".mzML"
            shutil.move(src_fname, dst_fname) # TODO: replace copy to rename (move)
        print("Created mzML files:\n", os.listdir(mzml_files_dir))

        # Use uncalibrared mgf files to get unlabeled mgf files from .d
        mgf_files_dir = os.path.join(MGF_DATA_DIR, dset_name)
        for file_path in mzml_files:
            fname = os.path.splitext(file_path)[0]
            src_fname = fname + "_uncalibrated.mgf"
            dst_fname = os.path.join(mgf_files_dir, os.path.basename(fname + ".mgf"))
            shutil.move(src_fname, dst_fname) # TODO: replace copy to move
        print("Created unlabeled mgf files\n", os.listdir(mgf_files_dir))


def run_database_search_split(dset_name, db_w_decoys_path, db_search_config):
    # Search for mzml files in the directory
    search_ext = db_search_config.ext
    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    mzml_files = [f for f in os.listdir(mzml_files_dir) if os.path.splitext(f)[1].lower() == search_ext]
    mzml_files = [os.path.join(mzml_files_dir, f) for f in mzml_files]

    # TODO: better add to closed.params
    # db_w_decoys_path -- path to DB file with decoys and contams

    n_db_splits = db_search_config.n_db_splits
    params_file = os.path.join(mzml_files_dir, "closed.params")
    # Iterate over each mzML file and run the msfragger command
    for mzml_file in mzml_files:
        print(f"PROCESSING {mzml_file}")
        cmd = [
            "python",
            SPLIT_SCRIPT_PATH,
            f"{n_db_splits}",
            '"java -jar -Dfile.encoding=UTF-8 -Xmx160G"',
            MSFRAGGER_PATH,
            params_file,
            mzml_file,
        ]

        # Run the command
        subprocess.run(" ".join(cmd), shell=True, check=True)
        print(f"PROCESSED {mzml_file}\n\n")

    # Print results
    print("DB search results (.pepXML, .pin):\n", os.listdir(mzml_files_dir))


def get_psm_rescoring_features(dset_name, rescoring_config):
    """Create PSMs rescoring features with MSBooster."""
    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    rescored_files_dir = os.path.join(RESCORED_DATA_DIR, dset_name)
 
    # select all the .mzml files from mzml_files_dir (# only select with fnames in files_list?)
    mzml_files = [f for f in os.listdir(mzml_files_dir) if os.path.splitext(f)[1].lower() == ".mzml"]
    mzml_files = [os.path.join(mzml_files_dir, f) for f in mzml_files]
    print(".mzML files available for rescoring:\n", mzml_files)

    # select .pin files with fnames in files_list
    # TODO: check if there are no problems with existing _rescore.pin files
    pin_files = [f for f in os.listdir(mzml_files_dir) if os.path.splitext(f)[1].lower() == ".pin"]
    pin_files = [os.path.join(mzml_files_dir, f) for f in pin_files]
    print(".pin files available for rescoring:\n", pin_files)

    file_prefix = "rescore" # TODO: move outside?

    options = [
        "--DiaNN",
        DIANN_PATH,
        "--KoinaURL",
        KOINA_URL,
        "--editedPin",
        file_prefix,
        "--paramsList",
        MSBOOSTER_BASE_PARAMS,
        "--mzmlDirectory",
        *mzml_files,
        "--pinPepXMLDirectory",
        *pin_files,
        "--outputDirectory",
        rescored_files_dir,
    ]
    # TODO: Parse additional params from config if provided
    for arg in [*rescoring_config.feat_pred_params.items()]:
        options += list(map(str, arg))
    
    cmd = [
        "java",
        "-Xmx64G",
        "-jar",
        MSBOOSTER_PATH,
        *options,
    ]
    print("MSBOOSTER DEBUG:\n") # TODO for debug only, remove
    print(" ".join(cmd)) # TODO
    subprocess.run(" ".join(cmd), shell=True, check=True)
    print("Created PSMs features (_rescore.pin):\n", os.listdir(mzml_files_dir))


def run_psm_rescoring(dset_name, rescoring_config, files_list):
    """Run Percolator for PSMs rescoring (using MSBooster features)."""
    # TODO: move outside (to constants?)
    num_threads = 3
    test_fdr = 0.01
    train_fdr = 0.01
    file_prefix = "rescore"

    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    rescored_files_dir = os.path.join(RESCORED_DATA_DIR, dset_name)

    # Merge together PSMs features for all the _rescore.pin files in files_list
    print(files_list, "\n")

    # Get number of columns
    fname = list(files_list.keys())[0]
    file_path = os.path.join(mzml_files_dir, f"{fname}_rescore.pin")
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    # Split the first line into column names
    column_names = first_line.split("\t")
    n_cols = len(column_names)

    dfs = [
        pd.read_csv(
            os.path.join(mzml_files_dir, f"{fname}_{file_prefix}.pin"),
            usecols=list(range(n_cols)),
            sep="\t"
        ) for fname in files_list
    ]
    print(len(dfs), "\n")
    df = pd.concat(dfs, axis=0).reset_index(drop=True)
    print(df.shape)
    # Save merged PSMs features df to be used by Percolator
    df.to_csv(
        os.path.join(rescored_files_dir, f"{file_prefix}.pin"), 
        sep="\t", 
        index=False
    )

    input_file = os.path.join(rescored_files_dir, f"{file_prefix}.pin")
    weights_file = os.path.join(rescored_files_dir, f"{file_prefix}.percolator.weights.csv")
    target_psms = os.path.join(rescored_files_dir, f"{file_prefix}.percolator.psms.txt")
    decoy_psms = os.path.join(rescored_files_dir, f"{file_prefix}.percolator.decoy.psms.txt")
    target_peptides = os.path.join(rescored_files_dir, f"{file_prefix}.percolator.peptides.txt")
    decoy_peptides = os.path.join(rescored_files_dir, f"{file_prefix}.percolator.decoy.peptides.txt")
    # log_file = os.path.join(rescored_files_dir, f"{file_prefix}.log")

    cmd = f"percolator --weights {weights_file} \
            --num-threads {num_threads} \
            --subset-max-train 500000 \
            --post-processing-tdc \
            --testFDR {test_fdr} \
            --trainFDR {train_fdr} \
            --results-psms {target_psms} \
            --decoy-results-psms {decoy_psms} \
            --results-peptides {target_peptides} \
            --decoy-results-peptides {decoy_peptides} \
            {input_file}"
    subprocess.run(cmd, shell=True, check=True)
    print(
        "PSMs rescoring results (.percolator.psms.txt):\n", 
        os.listdir(rescored_files_dir)
    )


def get_filename(psm_id: str):
    """Assumes that there are no `.` in the file name."""
    return psm_id.split(".")[0]


def format_peptide_notation(sequence: str):
    """TODO: PTMs may need conversion to ProForma notation."""
    # remove cleavage sites
    if (
        re.match(r"[A-Z-_].*.[A-Z-_]", sequence) is not None
    ):  # check is not mandatory
        sequence = sequence[2:-2]
    return sequence


def collect_dataset_tags(config):
    if os.path.exists(DATASET_TAGS_PATH):
        tags_df = pd.read_csv(DATASET_TAGS_PATH, sep="\t").to_dict('records')
    else:
        tags_df = []
        
    # add reference proteome information (for each dataset)
    dset_tags = {"proteome": config.db_search.database_path}
    # add dataset property tags
    dset_tags.update({tag.name: 1 for tag in config.tags})
    print(f"Dataset {config.name} tags:\n", dset_tags)
    
    tags_df.append({"dataset": config.name, **dset_tags})
    tags_df = pd.DataFrame(tags_df).fillna(0)
    tags_df = tags_df.drop_duplicates(subset="dataset", keep="last")
    tags_df.to_csv(DATASET_TAGS_PATH, index=False, sep="\t")
    print(f"Written to {DATASET_TAGS_PATH}")




# # Set maximum number of spectra per MGF file
# MAX_SPECTRA_PER_FILE = 20000

# # Store the 0-based indices of spectra with respect to new split files
# spectra_idxs_0 = []

# for fname in tqdm(files_list):
#     input_path = os.path.join(mgf_files_dir, fname + ".mgf")
#     spectra = mgf.read(input_path)
    
#     # Filter spectra based on "charge"
#     spectra_filtered = []
#     for spectrum in tqdm(spectra):
#         if "charge" in spectrum["params"]:
#             spectra_filtered.append(spectrum)
    
#     # Split filtered spectra into multiple files if needed
#     num_filtered_spectra = len(spectra_filtered)
#     num_splits = (num_filtered_spectra // MAX_SPECTRA_PER_FILE) + (1 if num_filtered_spectra % MAX_SPECTRA_PER_FILE > 0 else 0)
    
#     split_idx_0 = []
    
#     for k in range(num_splits):
#         # Define the range of spectra to include in the current split file
#         start_idx = k * MAX_SPECTRA_PER_FILE
#         end_idx = min((k + 1) * MAX_SPECTRA_PER_FILE, num_filtered_spectra)
#         spectra_chunk = spectra_filtered[start_idx:end_idx]
        
#         # Create the new MGF file with split spectra
#         output_filename = f"{fname}_{k}.mgf"
#         output_path = os.path.join(mgf_files_dir, output_filename)
        
#         mgf.write(spectra_chunk, output_path)
#         print(f"Written {len(spectra_chunk)} spectra to {output_path}.")
        
#         # Store the 0-based index and spectrum title for this chunk
#         idxs_0 = {idx: spectrum["params"]["title"] for idx, spectrum in enumerate(spectra_chunk)}
#         idxs_0 = pd.Series(idxs_0).reset_index()
#         idxs_0.columns = ["idx_0", "title"]
#         idxs_0["filename"] = output_filename
        
#         # Append indices for this split
#         split_idx_0.append(idxs_0)
    
#     # Concatenate indices from all splits for this file
#     if split_idx_0:
#         spectra_idxs_0.append(pd.concat(split_idx_0, axis=0).reset_index(drop=True))
    
#     # Remove the original MGF file after splitting
#     os.remove(input_path)
#     print(f"Original file {input_path} removed.")
    
# # Combine indices from all files
# spectra_idxs_0 = pd.concat(spectra_idxs_0, axis=0).reset_index(drop=True)

# # Update the results dataframe with the new idx_0 and filename fields
# results_df = pd.merge(results_df, spectra_idxs_0, on="title")
# results_df["spectrum_id"] = results_df["filename"] + ":" + results_df["idx_0"].astype(str)