#!/bin/bash
spectra_dir="$1"
output_dir="./outputs"
overlay_size=1024

# Clean output dir 
rm -rf "$output_dir"
# Create the output directory if it doesn't exist
mkdir "$output_dir"

# List input files
echo "Spectra data dir: $spectra_dir"
ls "$spectra_dir"/*.mgf

# Loop through each algorithm in the algorithms directory
for algorithm_dir in algorithms/*; do
    if [ -d "$algorithm_dir" ]; then
        algorithm_name=$(basename "$algorithm_dir")
        
        echo "$algorithm_dir"
        echo "Running $algorithm_name"
        
        # Remove an existing container overlay, if any
        rm -rf "algorithms/${algorithm_name}/overlay.img"
        # Create writable overlay for the container
        apptainer overlay create --fakeroot --size $overlay_size --sparse "algorithms/${algorithm_name}/overlay.img"

        # Calculate predicitons
        echo "RUN ALGORITHM"
        apptainer exec --fakeroot --nv \
            --overlay "algorithms/${algorithm_name}/overlay.img" \
            -B "${spectra_dir}":/algo/data \
            "algorithms/${algorithm_name}/container.sif" \
            bash -c "cd /algo && ls && ./make_predictions.sh data" # TODO: remove ls (for debug only)
        
        # Collect predictions in output_dir
        echo "EXPORT PREDICTIONS"
        apptainer exec --fakeroot \
            --overlay "algorithms/${algorithm_name}/overlay.img" \
            -B "${output_dir}"/:/algo/outputs \
            "algorithms/${algorithm_name}/container.sif" \
            bash -c "cp /algo/outputs.csv /algo/outputs/${algorithm_name}_output.csv"

    fi
done

# Setup cluster environment for metric calculation # TODO: remove / move outside
# module restore pytorch-modules 
module load calcua/all
module load Python/3.10.4-GCCcore-11.3.0
module load SciPy-bundle/2022.05-foss-2022a
module load PyYAML/6.0-GCCcore-11.3.0
module load PyTorch/1.13.1-foss-2022a-CUDA-11.7.0
# TODO: resolve the issue that full module list is required to run on cluster
source $VSC_DATA/venvs/bm_test/bin/activate

# Evaluate predictions
echo "EVALUATE PREDICTIONS"
python evaluate.py "$output_dir/" "$spectra_dir"

cat metrics.csv # 1) remove 2) replace / add building plots

# Clean environment # TODO: also preferably outside
deactivate 
module purge
