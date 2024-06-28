import argparse
import os
from oktoberfest.runner import run_job
from dataset_utils import *
from dataset_config import get_config


# Paths parsing
parser = argparse.ArgumentParser()
parser.add_argument("--config_path", type=str, help="path to the dataset config file")
args = parser.parse_args()

# Config setup
config = get_config(args.config_path)
dset_id = config.download.dset_id
dset_name = config.name

# Output dirs setup
for data_dir in [
    RAW_DATA_DIR, 
    MZML_DATA_DIR, 
    RESCORED_DATA_DIR, 
    MGF_DATA_DIR, 
    RESCORE_PARAMS_DIR, 
    DATASET_STORAGE_DIR
]:
    os.makedirs(data_dir, exist_ok=True)

# Create dirs for intermediate dataset files
raw_files_dir = os.path.join(RAW_DATA_DIR, dset_id)
mzml_files_dir = os.path.join(MZML_DATA_DIR, dset_name) # TODO: renamve to smth search_files_dir? 
rescored_files_dir = os.path.join(RESCORED_DATA_DIR, dset_name)
mgf_files_dir = os.path.join(MGF_DATA_DIR, dset_name)
labeled_mgf_files_dir = os.path.join(DATASET_STORAGE_DIR, dset_name)

for data_dir in [
    raw_files_dir, mzml_files_dir, rescored_files_dir, mgf_files_dir, labeled_mgf_files_dir
]:
    os.makedirs(data_dir, exist_ok=True)

print("Creating dataset at:", labeled_mgf_files_dir)
labeled_fnames = [os.path.splitext(f)[0] for f in os.listdir(labeled_mgf_files_dir)]
print("Existing labeled files:")
print(labeled_fnames)

# Get list of files to process 
files_list = get_files_list(config.download)
print("Processing files:\n", list(files_list.keys()))

# If already have pepxml, no need to run DB search (TODO: add flag for FORCED RE-RUN)
if not set(fname.lower() + ".pin" for fname in files_list).issubset(
        set(fname.lower() for fname in os.listdir(mzml_files_dir))
    ):
    # Prepare mzml (or other DB search input format) files if they don't exist
    if not set(fname.lower() + config.db_search.ext for fname in files_list).issubset(
            set(fname.lower() for fname in os.listdir(mzml_files_dir))
        ):
        # Download raw files from repository if needed      
        if not set(fname.lower() + config.download.ext for fname in files_list).issubset(
                set(fname.lower() for fname in os.listdir(raw_files_dir))
            ):
            download_files(config.download, files_list)

        # распаковать ЕСЛИ они еще не распакованы
        if config.download.ext.endswith(".zip"):
            base_file_ext = config.download.ext[:-len(".zip")]
            # If files can be directly processed by MSFragger (.d, .mzml), unpack to mzml_files_dir
            if base_file_ext == ".d":
                unpack_dir = mzml_files_dir
            else:
                unpack_dir = raw_files_dir

            if not set(fname.lower() + base_file_ext for fname in files_list).issubset(
                set(fname.lower() for fname in os.listdir(unpack_dir))
            ):
                print(f"Files will be unpacked to {unpack_dir}")
                fnames = list(files_list.values())
                for fname in fnames:
                    file_path = os.path.join(raw_files_dir, fname)
                    shutil.unpack_archive(filename=file_path, extract_dir=unpack_dir)

        if config.download.ext == ".raw":
            # Convert raw files to mzml (ThermoRawFileParser) # TODO: check for files from files_list only?
            convert_raw(dset_id, files_list, mzml_files_dir, target_ext=config.db_search.ext) # also always mzml??

        elif config.download.ext == ".wiff":
            # Convert raw files to mzml (msconvert) # TODO: check for files from files_list only?
            convert_wiff(dset_id, files_list, mzml_files_dir, target_ext=config.db_search.ext) # always mzml?
        
        else:
            print(f"Unknown file extension {config.download.ext}")
    
    # TODO: add contaminants (before or after decoys?)
    # Generate decoys for DB search # TODO: add decoys generation only if doesn't exist
    db_w_decoys_path = generate_decoys_fasta(db_file=config.db_search.database_path)

    # Run DB search (MSFragger)
    run_database_search(dset_name, db_w_decoys_path, config.db_search)

# Rescore DB search results (MSBooster + Percolator)
# Create features with MSBooster
# If already have _rescore.pin, no need to run feature prediction (TODO: add flag for FORCED RE-RUN)
file_prefix = "rescore"
if not set(fname.lower() + f"_{file_prefix}.pin" for fname in files_list).issubset(
        set(fname.lower() for fname in os.listdir(mzml_files_dir))
    ):
    # TODO: ideally should also check whether mzml files exist and add them otherwise
    get_psm_rescoring_features(dset_name, config.rescoring)

# Run Percolator on created features
# TODO: any checks whether PSMs rescoring results already available?
if not f"{file_prefix}.percolator.psms.txt" in os.listdir(rescored_files_dir):
    run_psm_rescoring(dset_name, config.rescoring, files_list)

# Prepare unlabeled mgf files if they don't exist
if not set(fname.lower() + ".mgf" for fname in files_list).issubset(
        set(fname.lower() for fname in os.listdir(mgf_files_dir))
    ):
    # Download raw files from repository if needed  
    if not set(fname.lower() + config.download.ext for fname in files_list).issubset(
            set(fname.lower() for fname in os.listdir(raw_files_dir))
        ):
        download_files(config.download, files_list)

    # распаковать ЕСЛИ они еще не распакованы
    if config.download.ext.endswith(".zip"):
        base_file_ext = config.download.ext[:-len(".zip")]
        # If files can be directly processed by MSFragger (.d, .mzml), unpack to mzml_files_dir
        if base_file_ext == ".d":
            unpack_dir = mzml_files_dir
        else:
            unpack_dir = raw_files_dir

        if not set(fname.lower() + base_file_ext for fname in files_list).issubset(
            set(fname.lower() for fname in os.listdir(unpack_dir))
        ):
            print(f"Files will be unpacked to {unpack_dir}")
            fnames = list(files_list.values())
            for fname in fnames:
                file_path = os.path.join(raw_files_dir, fname)
                shutil.unpack_archive(filename=file_path, extract_dir=unpack_dir)

    if config.download.ext == ".raw":
        # Convert raw files to mgf (ThermoRawFileParser)
        convert_raw(dset_id, files_list, mgf_files_dir, target_ext=".mgf")

    elif config.download.ext == ".wiff":
        # Convert wiff files to mgf (msconvert)
        convert_wiff(dset_id, files_list, mgf_files_dir, target_ext=".mgf")

    elif config.download.ext == ".d.zip":
        # After DB search files should be already created 
        # from uncalibrated.mgf produced by MSFagger
        print("No mgf files found for .d data")
    
    else:
        print(f"Unknown file extension {config.download.ext}")

# Load DB search + rescoring results
create_labeled_mgf(dset_name, labeled_mgf_files_dir, config.rescoring.q_val_threshold)
