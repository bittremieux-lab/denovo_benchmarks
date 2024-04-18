#!/bin/bash
spectra_dir="$1"
output_dir="./outputs"
overlay_size=1024

# TODO: make optional, only if recalculation is needed (flag)
# # Clean output dir 
# rm -rf "$output_dir"
# # Create the output directory if it doesn't exist
# mkdir "$output_dir"

# List input files
echo "Spectra data dir: $spectra_dir"
ls "$spectra_dir"/*.mgf

# Loop through each algorithm in the algorithms directory
for algorithm_dir in algorithms/*; do
    if [ -d "$algorithm_dir" ]; then
        algorithm_name=$(basename "$algorithm_dir")
        output_file="$output_dir/${algorithm_name}_output.csv"
        echo "$output_file"
        
        # Check if the output file does not exist
        if [ ! -e "$output_file" ]; then
            echo "Processing algorithm: $algorithm_name"

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

        else
            echo "Skipping algorithm: $algorithm_name. Output file already exists."
        fi
    fi
done

# Evaluate predictions
echo "EVALUATE PREDICTIONS"
apptainer exec --fakeroot "evaluation.sif" \
    bash -c "python evaluate.py ${output_dir}/ ${spectra_dir}"
