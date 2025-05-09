import argparse
import os
import re
import shutil
from tqdm import tqdm
from pyteomics import mzml
from dataset_utils import *
from dataset_config import get_config, DatasetTag


# Here I describe the whole logic/structure of the database searching workflow
# that this script should be able to perform.

# For each dataset, for each file, we want to eventually get:
# - PSM labels from database search with MSFragger
# (workflow: run MSFragger, create features with MSBooster, 
# optionally with Prosit to predict RT & Intensity features, run Percolator rescoring)
# - PSM labels from database search with MSGF+
# (workflow: run MSGF+ (separately for targets and decoys), 
# combine results into featrures with msgf2pin, run Percolator rescoring)
# - PSM labels from database search with Comet
# (workflow: run Comet, collect result files (.pep.xml, .txt) to a separate folder, rename to .txt to .pin,
# run Percolator rescoring) 
# - Full spectra files (both labeled and unlabeled) in mgf format,
# saved in mgf_files_dir in chunks of MAX_SPECTRA_PER_FILE spectra.

# We would like to avoid any unnecessary re-running of the database search or rescoring steps,
# and especially we don't want to re-write the mgf files if they already exist.
# But note that original file has name fname.raw/fname.mzml/etc, 
# and the corresponding mgf files with have names fname_0.mgf, fname_1.mgf, etc. for each chunk.

# Previously, the script was only running MSFragger DB search & rescoring, 
# so it was always performing the same pipeline:
# 0. Load dataset config file, get list of files to process
# 1. Download raw files from repository if needed
# 2. Unpack if unpacking needed (raw files come as .zip) and not unpacked yet
# 3. Convert raw files to mzml (msconvert) if mzml files not exist yet in mzml_files_dir
# 4. Run DB search with MSFragger if MSFragger results (.pepXML, .pin) not exist yet in mzml_files_dir
# 5. If want to use deep learning-based features, run MSBooster to get _rescore.pin
# 6. Run Percolator rescoring on the created features: 
# mzml_file_dir/*_rescore.pin if use deep learning-based features, else mzml_file_dir/*.pin
# 7. Collect PSM labels from the rescoring results, 
# Filter them based on q-value,
# Give each PSM a spectrum_id based on:
# filename: name of the chunk(!) mgf file containing this spectrum
# idx_0: 0-based index of this spectrum in the chunk mgf file
# 8. Save chunk mgf files to the mgf_files_dir: 
# if mgf files not exist yet,get them from the raw files, 
# then iterate over mgf files and make sure that only spectra containing charge are kept,
# split them into chunks of MAX_SPECTRA_PER_FILE spectra,
# and save them to the mgf_files_dir with names fname_0.mgf, fname_1.mgf, etc.
# In the original pipeline, spectrum_ids for labels are also collected during this process. 
# However, this needs to be adjusted, since we don't want to re-save mgf files after each DB search.

# Now, I need to adjust and restructure this script workflow to be able to handle the following scenarios:
# - 1. MSFragger DB search + rescoring: 
# get files list, make sure mzml files exist (if not, download them and convert from raw), run MSFragger,
# run MSBooster to get _rescore.pin (consider it mandatory, unless corresponding lines are not commented out), 
# run Percolator rescoring on the created features, collect Percolator results and save them as labels.csv to labels_path.
# If chunk mgf files already exist for all the files in the files_list, they should be used to get spectrum_ids.
# - 2. MSGF+ DB search + rescoring:
# get files list, make sure mzml files exist (if not, download them and convert from raw), run MSGF+,
# run msgf2pin to get .pin files, run Percolator rescoring on the created features, 
# collect Percolator results and save them as labels.csv to labels_path.
# If chunk mgf files already exist for all the files in the files_list, they should be used to get spectrum_ids.
# - 3. Comet DB search + rescoring:
# get files list, make sure mzml files exist (if not, download them and convert from raw), run Comet,
# collect result files (.pep.xml, .txt) to a separate folder, rename to .txt to .pin,
# run Percolator rescoring on the created features,
# collect Percolator results and save them as labels.csv to labels_path.
# - 4. Save mgf files to the mgf_files_dir:
# If chunk mgf files already exist for all the files in the files_list, nothing to do.
# If chunk mgf files don't exist yet, get them from full mgf files. Make sure that only spectra containing charge are kept.
# If mgf files not exist yet, get them from the raw files, and then split into chunks of MAX_SPECTRA_PER_FILE spectra.
# Original full mgf files should be removed after splitting.
# This mgf creation pipeline should be run if chunk mgf files were not found during any of the database searching workflows. 

# Important note:
# Looks like MSBooster is only needed to get deep learning-based features for rescoring:
# Based on this part:
# # [! uncomment to use deep learning-based features in rescoring]
#     file_path = os.path.join(mzml_files_dir, f"{fname}_{file_prefix}.pin")
#     # file_path = os.path.join(mzml_files_dir, f"{fname}.pin")
# So essentialy what we do for MSFragger is:
# - get DB search results as .pin
# - (!) if we want deep learning-based features, run MSBooster to get also _rescore.pin
# - run Percolator rescoring:
    # - (!) if we want deep learning-based features, use _rescore.pin features
    # - else, just use original .pin features from MSFragger!
# Then for all the other tools, the equivalent workflow will be:
# (taking into account that we DON'T need deep learning-based features for them!)
# - run DB search with MSGF+ or Comet
# - collect their outputs in separate dirs mzml_files_dir/{tool_name}_features in .pin format
# - run Percolator rescoring on the collected features - .pin files from mzml_files_dir/{tool_name}_features
# - Percolator outputs should be saved in separate dirs rescored_files_dir/{tool_name}_rescored
# And only for MSFragger we will keep the option to use deep learning-based features
# by running MSBooster on the original .pin files (optionally)

Q_VAL_THRESHOLD = 0.01
Q_VAL_THRESHOLD_SYNTHETIC = 0.001
MAX_SPECTRA_PER_FILE = 20000

def check_mzml_files_exist(files_list, mzml_files_dir):
    # Check if mzML files exist
    mzml_files = [fname.lower() + ".mzml" for fname in files_list]
    return set(mzml_files).issubset(set(fname.lower() for fname in os.listdir(mzml_files_dir)))

def check_raw_files_exist(files_list, raw_files_dir, ext=".raw"):
    # Check if raw files exist
    raw_files = [fname.lower() + ext for fname in files_list]
    return set(raw_files).issubset(set(fname.lower() for fname in os.listdir(raw_files_dir)))

def check_unpacked_files_exist(files_list, unpack_dir, base_file_ext):
    # Check if unpacked files exist
    unpacked_files = [fname.lower() + base_file_ext for fname in files_list]
    return set(unpacked_files).issubset(set(fname.lower() for fname in os.listdir(unpack_dir)))

def prepare_mzml_files(dset_id, files_list, raw_files_dir, mzml_files_dir, download_config):
    """Ensure mzML files exist by downloading and converting raw files if needed."""
    
    if not check_mzml_files_exist(files_list, mzml_files_dir):
        # Download raw files if needed
        if not check_raw_files_exist(files_list, raw_files_dir, download_config.ext):
            download_files(download_config, files_list)

        # Unpack raw files if needed
        if download_config.ext.endswith(".zip"):
            base_file_ext = download_config.ext[:-len(".zip")]
            # If files can be directly processed by MSFragger (.d, .mzml), unpac to mzml_files_dir
            unpack_dir = mzml_files_dir if base_file_ext == ".d" else raw_files_dir
            if not check_unpacked_files_exist(files_list, unpack_dir, base_file_ext):
                for fname in files_list.values():
                    file_path = os.path.join(raw_files_dir, fname)
                    shutil.unpack_archive(filename=file_path, extract_dir=unpack_dir)

        # Convert raw files to mzML
        if download_config.ext in [".raw", ".wiff"]:
            convert_raw(dset_id, files_list, mzml_files_dir, target_ext=".mzml")
    
    mzml_files = [fname for fname in os.listdir(mzml_files_dir) if fname.lower().endswith(".mzml")]
    print("mzML files are ready in", mzml_files_dir, ":\n", mzml_files)

def check_db_search_results_exist(search_tool, files_list, mzml_files_dir):
    """Check if database search results already exist."""
    if search_tool == "msfragger":
        expected_files = [os.path.join(mzml_files_dir, f"{fname}.pin") for fname in files_list]
    elif search_tool == "msgf":
        expected_files = [os.path.join(mzml_files_dir, "msgf_features", f"{fname}.pin") for fname in files_list]
    elif search_tool == "comet":
        expected_files = [os.path.join(mzml_files_dir, "comet_features", f"{fname}.pin") for fname in files_list]
    else:
        raise ValueError(f"Unsupported search tool: {search_tool}")
    
    return all(os.path.exists(file) for file in expected_files)

def run_db_search(search_tool, dset_name, db_search_config, files_list, mzml_files_dir):
    """Run the appropriate DB search tool if results do not already exist."""
    if check_db_search_results_exist(search_tool, files_list, mzml_files_dir):
        print(f"Database search results already exist for {search_tool}. Skipping search.\n")
        return

    if search_tool == "msfragger":
        run_msfragger_search(dset_name, db_search_config)
    elif search_tool == "msgf":
        run_msgf_search(dset_name, db_search_config)
    elif search_tool == "comet":
        run_comet_search(dset_name, db_search_config)
    else:
        raise ValueError(f"Unsupported search tool: {search_tool}")

def map_psm_id_index_to_scan_id(results_df):
    
    def get_filename(psm_id):
        return psm_id.split("_SII_")[0]
    
    def map_mzml_index_to_scan_id(psm_id):
        filename = psm_id.split("_SII_")[0]

        psm_id = psm_id.split("_")
        mzml_index = int(psm_id[-5]) - 1
        scan_id = mzml_index_to_scan_id[filename][mzml_index]

        psm_id[-5] = str(scan_id)
        psm_id = "_".join(psm_id)
        return psm_id
    
    results_df["filename"] = results_df["PSMId"].apply(get_filename)
    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    mzml_files = results_df.filename.unique().tolist()
    mzml_index_to_scan_id = {}
    for fname in mzml_files:
        file_path = os.path.join(mzml_files_dir, fname + ".mzML")
        mzml_spectra = mzml.MzML(file_path)
        print(fname, len(mzml_spectra))

        index_to_scan_id_mapping = {}
        for spectrum in tqdm(mzml_spectra):
            idx = spectrum["index"]
            scan_id = int(spectrum["id"].split("=")[-1]) # ! Assume id="scanId=X" -- id doesn't contain any other information
            index_to_scan_id_mapping[idx] = scan_id
        mzml_index_to_scan_id[fname] = index_to_scan_id_mapping
        
    results_df["PSMId"] = results_df["PSMId"].apply(map_mzml_index_to_scan_id)
    results_df = results_df.drop("filename", axis=1)
    return results_df

def check_chunked_mgf_files_exist(files_list, mgf_files_dir):
    """Check if chunked MGF files exist."""
    chunked_mgf_files = [fname.lower() + "_0.mgf" for fname in files_list]
    return set(chunked_mgf_files).issubset(set(fname.lower() for fname in os.listdir(mgf_files_dir)))

def check_full_mgf_files_exist(files_list, mgf_files_dir):
    """Check if full MGF files exist."""
    full_mgf_files = [fname.lower() + ".mgf" for fname in files_list]
    return set(full_mgf_files).issubset(set(fname.lower() for fname in os.listdir(mgf_files_dir)))

def get_mgf_files_spectra_idxs(files_list, mgf_files_dir, raw_files_dir, dset_id, download_config):
    """Ensure chunked MGF files exist and extract spectra indices."""
    if not check_chunked_mgf_files_exist(files_list, mgf_files_dir):
        # Convert raw files to MGF if needed
        if not check_full_mgf_files_exist(files_list, mgf_files_dir):
            # Download raw files if needed
            if not check_raw_files_exist(files_list, raw_files_dir, download_config.ext):
                download_files(download_config, files_list)
            convert_raw(dset_id, files_list, mgf_files_dir, target_ext=".mgf")

        # Filter and chunk MGF files
        spectra_idxs_0 = []
        for fname in tqdm(files_list):
            input_path = os.path.join(mgf_files_dir, fname + ".mgf")
            spectra = mgf.read(input_path)
            spectra_filtered = [spectrum for spectrum in spectra if "charge" in spectrum["params"]]
            n_splits = len(spectra_filtered) // MAX_SPECTRA_PER_FILE + int(len(spectra_filtered) % MAX_SPECTRA_PER_FILE > 0)
            for k in range(n_splits):
                chunk = spectra_filtered[k * MAX_SPECTRA_PER_FILE:(k + 1) * MAX_SPECTRA_PER_FILE]
                output_fname = f"{fname}_{k}"
                output_path = os.path.join(mgf_files_dir, output_fname + ".mgf")
                if not os.path.exists(output_path):
                    mgf.write(chunk, output_path)
                idxs_0 = {idx: spectrum["params"]["title"] for idx, spectrum in enumerate(chunk)}
                idxs_0 = pd.DataFrame.from_dict(idxs_0, orient="index", columns=["title"]).reset_index()
                idxs_0 = idxs_0.rename(columns={"index": "idx_0"})
                idxs_0["filename"] = output_fname
                spectra_idxs_0.append(idxs_0)
        return pd.concat(spectra_idxs_0, axis=0).reset_index(drop=True)

    # If chunked MGF files already exist, just read them & extract spectra idxs
    spectra_idxs_0 = []
    for fname in tqdm(files_list):
        chunk_files = [f for f in os.listdir(mgf_files_dir) if re.fullmatch(fr"{fname}_[\d]+.mgf", f)]
        # chunk_files = [f for f in os.listdir(mgf_files_dir) if f.startswith(fname + "_") and f.endswith(".mgf")]
        for chunk_file in chunk_files:
            chunk_path = os.path.join(mgf_files_dir, chunk_file)
            spectra = mgf.read(chunk_path)
            idxs_0 = {idx: spectrum["params"]["title"] for idx, spectrum in enumerate(spectra)}
            idxs_0 = pd.DataFrame.from_dict(idxs_0, orient="index", columns=["title"]).reset_index()
            idxs_0 = idxs_0.rename(columns={"index": "idx_0"})
            idxs_0["filename"] = chunk_file[:-len(".mgf")]
            spectra_idxs_0.append(idxs_0)
    return pd.concat(spectra_idxs_0, axis=0).reset_index(drop=True)

def run_rescoring(search_tool, dset_name, rescoring_config, files_list, rescored_files_dir, rescore_file_prefix="rescore"):
    """Run Percolator rescoring for the specified search tool."""
    # TODO: check_rescoring_results_exist
    if f"{rescore_file_prefix}.percolator.psms.txt" in os.listdir(rescored_files_dir):
        print(f"Rescoring results already exist for {search_tool}. Skipping rescoring.\n")
        return

    tool_rescoring_features_dirs = {
        "msfragger": os.path.join(MZML_DATA_DIR, dset_name),
        "msgf": os.path.join(MZML_DATA_DIR, dset_name, "msgf_features"),
        "comet": os.path.join(MZML_DATA_DIR, dset_name, "comet_features"),
    }
    features_dir = tool_rescoring_features_dirs[search_tool]
    if not os.path.exists(features_dir):
        raise FileNotFoundError(f"Rescoring features directory not found for {search_tool}: {features_dir}")

    # Create a single merged rescoring features file here (from files in features_dir)
    file_paths = [os.path.join(features_dir, f"{fname}.pin") for fname in files_list]
    # [! uncomment to use deep learning-based features in rescoring]
    if search_tool == "msfragger":
        get_psm_rescoring_features(dset_name, rescoring_config)
        file_paths = [os.path.join(features_dir, f"{fname}_rescore.pin") for fname in files_list]
    
    skiprows = [1] if search_tool == "msgf" else None
    droprows = [i - 1 for i in skiprows] if skiprows else []
    dfs = []
    for file_path in file_paths:
        print(file_path)
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
        # Split the first line into column names
        column_names = first_line.split("\t")
        try:
            df = pd.read_csv(file_path, sep="\t", usecols=column_names)
        except:
            df = pd.read_csv(file_path, sep="\t", usecols=column_names, skipfooter=1)
        df = df.drop(droprows, axis=0).reset_index(drop=True)
    #     df = pd.read_csv(file_path, sep="\t", usecols=column_names, skiprows=skiprows) 
        dfs.append(df)
        print()
    print("Found rescoring features files:", len(dfs))

    df = pd.concat(dfs, axis=0).reset_index(drop=True)
    df = df.fillna(0)
    print("Merged rescoring features dataframe:", df.shape)
    # Save merged PSMs features df to be used by Percolator
    df.to_csv(
        os.path.join(rescored_files_dir, f"{rescore_file_prefix}.pin"), 
        sep="\t", 
        index=False
    )
    # Run rescoring
    run_psm_rescoring(dset_name, rescoring_config, rescored_files_dir, rescore_file_prefix)

def get_spectrum_title(psm_id, search_tool="msfragger"):
    if search_tool == "comet":
        title = psm_id.split("/")[-1]
        title = title.split("_")
        filename, scan_id, charge, rank = title[:-3], title[-3], title[-2], title[-1]
        filename = "_".join(filename)
        title = ".".join([filename, scan_id, scan_id, charge])
        return title
    
    elif search_tool == "msgf":
        title = psm_id.split("_")
        # PSM "index": SII_13752_1 = title[-6:-3], 
        filename, scan_id, charge, rank = title[:-6], title[-5], title[-2], title[-1]
        filename = "_".join(filename)
        title = ".".join([filename, scan_id, scan_id, charge])
        return title
    
    # else (msfragger, msgf):
    title = "_".join(psm_id.split("_")[:-1])
    return title

def collect_labels(
    search_tool, rescored_files_dir, rescore_file_prefix, spectra_idxs_0, labels_path, 
    q_val_threshold=Q_VAL_THRESHOLD,
    ):
    """Collect PSM labels from Percolator results and save them."""
    results_path = os.path.join(rescored_files_dir, f"{rescore_file_prefix}.percolator.psms.txt")
    # results_df = pd.read_csv(results_path, sep="\t")
    results_df = pd.read_csv(
        results_path, sep="\t", 
        usecols=['PSMId', 'score', 'q-value', 'posterior_error_prob', 'peptide', 'proteinIds']
    )
    results_df = results_df[results_df["q-value"] < q_val_threshold][["PSMId", "peptide", "q-value"]]

    if search_tool == "msgf" and DatasetTag.agilent in config.tags:
        results_df = map_psm_id_index_to_scan_id(results_df)
    results_df["title"] = results_df["PSMId"].apply(lambda x: get_spectrum_title(x, search_tool=search_tool))
    results_df = pd.merge(results_df, spectra_idxs_0, on="title")

    results_df["spectrum_id"] = results_df["filename"] + ":" + results_df["idx_0"].astype(str)
    results_df["peptide"] = results_df["peptide"].apply(format_peptide_notation)
    sequences_true = results_df[["peptide", "spectrum_id"]].rename(columns={"peptide": "seq"})
    sequences_true.to_csv(labels_path, index=False)

# Main script logic
if __name__ == "__main__":
    # Paths parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_path", type=str, help="path to the dataset config file")
    parser.add_argument("--search_tool", type=str, default="msfragger", help="database search tool to use (e.g., msfragger, msgf, comet)")
    args = parser.parse_args()

    # Config setup
    config = get_config(args.config_path)
    dset_id = config.download.dset_id
    dset_name = config.name
    search_tool = args.search_tool
    rescore_file_prefix = f"{search_tool}_rescore" if search_tool != "msfragger" else "rescore"

    # Main output dirs setup
    for data_dir in [RAW_DATA_DIR, MZML_DATA_DIR, RESCORED_DATA_DIR, DATASET_STORAGE_DIR]:
        os.makedirs(data_dir, exist_ok=True)

    # Create dirs for intermediate & final dataset files
    raw_files_dir = os.path.join(RAW_DATA_DIR, dset_id)
    mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name)
    rescored_files_dir = os.path.join(RESCORED_DATA_DIR, dset_name, rescore_file_prefix)
    mgf_files_dir = os.path.join(DATASET_STORAGE_DIR, dset_name, "mgf")
    for data_dir in [raw_files_dir, mzml_files_dir, rescored_files_dir, mgf_files_dir]:
        os.makedirs(data_dir, exist_ok=True)

    print("Creating dataset at:", mgf_files_dir)
    files_list = get_files_list(dset_name, config.download)
    print("Processing files:\n", files_list)#list(files_list.keys()))

    # Prepare mzML files
    # TODO: make optional: we don't need mzml files if DB search outputs for search_tool already exist
    # if config.db_search.ext == ".mzml":
    prepare_mzml_files(dset_id, files_list, raw_files_dir, mzml_files_dir, config.download)

    # Run DB search
    run_db_search(search_tool, dset_name, config.db_search, files_list, mzml_files_dir)

    # Run rescoring
    run_rescoring(search_tool, dset_name, config.rescoring, files_list, rescored_files_dir, rescore_file_prefix)

    # Prepare MGF files and extract spectra indices
    spectra_idxs_0 = get_mgf_files_spectra_idxs(files_list, mgf_files_dir, raw_files_dir, dset_id, config.download)

    # Collect labels
    q_val_threshold = Q_VAL_THRESHOLD_SYNTHETIC if (DatasetTag.synthetic in config.tags) else Q_VAL_THRESHOLD
    labels_fname = "labels.csv" if search_tool == "msfragger" else f"{search_tool}_labels.csv"
    labels_path = os.path.join(DATASET_STORAGE_DIR, dset_name, labels_fname)
    collect_labels(search_tool, rescored_files_dir, rescore_file_prefix, spectra_idxs_0, labels_path, q_val_threshold)

    # Add dataset tags
    collect_dataset_tags(config)
