#!/bin/bash
algorithm_name="$1"
test_spectra_dir="datasets/sample_data"
overlay_size=512
test_output_dir="./test_outputs"

# Create the output directory if it doesn't exist
mkdir -p "${test_output_dir}"

# List input files
echo "Processing test dataset"
ls "$test_spectra_dir"/*.mgf

# Check if the output file does not exist
echo "Processing algorithm: $algorithm_name"

# [Optional] Create a test container (from container.def)
# apptainer build --fakeroot \
#     "algorithms/${algorithm_name}/test_container.sif" \
#     "algorithms/${algorithm_name}/container.def"

# Create writable overlay for the container
apptainer overlay create --fakeroot --size $overlay_size \
    --sparse "algorithms/${algorithm_name}/test_overlay.img"

# Calculate predicitons
echo "RUN ALGORITHM"
apptainer exec --fakeroot --nv \
    --overlay "algorithms/${algorithm_name}/test_overlay.img" \
    -B "${test_spectra_dir}":/algo/data \
    "algorithms/${algorithm_name}/container.sif" \
    bash -c "cd /algo && ./make_predictions.sh data"

# Collect predictions in output_dir
echo "EXPORT PREDICTIONS"
apptainer exec --fakeroot \
    --overlay "algorithms/${algorithm_name}/test_overlay.img" \
    -B "${test_output_dir}":/algo/outputs \
    "algorithms/${algorithm_name}/container.sif" \
    bash -c "cp /algo/outputs.csv /algo/outputs/test_output.csv"

# Check predictions output format
# Use evaluation.sif container for running output format tests
echo "VALIDATE PREDICTIONS OUTPUT FORMAT"
apptainer exec --fakeroot "evaluation.sif" \
    bash -c "python test_output_format.py"

echo "OUTPUT FORMAT VALIDATED."

# Remove test container image and overlay
# TODO: make a flag to not remove container if needed
# rm -rf "algorithms/${algorithm_name}/test_container.sif"
rm -rf "algorithms/${algorithm_name}/test_overlay.img"
rm -rf "${test_output_dir}"
