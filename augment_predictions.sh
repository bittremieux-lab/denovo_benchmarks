#!/bin/bash
dset_dir="$1"
algorithm="$2"
spectra_dir="$dset_dir/mgf"
output_root_dir="./outputs"

dset_name=$(basename "$dset_dir")

# Echo message based on whether an algorithm is provided
if [ -z "$algorithm" ]; then
    echo "Augment predictions for all algorithms on dataset $dset_name."
else
    echo "Augment predictions for $algorithm on dataset $dset_name."
fi

# List input files
echo "Processing dataset: $dset_name ($dset_dir)"
ls "$spectra_dir"/*.mgf

# Loop through each algorithm in the algorithms directory
for algorithm_dir in algorithms/*; do

    if [ -d "$algorithm_dir" ] && [ $(basename "$algorithm_dir") != "base" ]; then
        algorithm_name=$(basename "$algorithm_dir")

        # If an algorithm is specified, only continue if algorithm_name matches
        if [ -z "$algorithm" ] || [ "$algorithm_name" == "$algorithm" ]; then

            # Extract latest container version
            algorithm_version=$(grep -m 1 "container_version:" "algorithms/${algorithm_name}/versions.log" | awk -F'"' '{print $2}')
            
            # Validate version extraction worked
            if [ -z "$algorithm_version" ]; then
                echo "Warning: Could not extract container_version from algorithms/${algorithm_name}/versions.log. Skipping." >&2
                continue
            fi

            echo "Processing algorithm: $algorithm_name (version $algorithm_version)"
            
            # Build output directory path: outputs/{algorithm}/{version}/{dataset}/
            output_dir="$output_root_dir/$algorithm_name/$algorithm_version/$dset_name"
            output_file="$output_dir/output.csv"
            echo "Output file: $output_file"

            # Check if output file exists
            if [ ! -e "$output_file" ]; then
                echo "Warning: Output file not found for $algorithm_name. Skipping augmentation."
                continue
            fi

            # Augment algorithm predictions with RT and SA (if not already present)
            echo "AUGMENT PREDICTIONS"
            apptainer exec --fakeroot --env-file .env "evaluation.sif" \
                bash -c "python -m evaluation.augment_predictions --output_dir ${output_dir} --data_dir ${dset_dir}"
                
            # bash -c "python -m evaluation.augment_predictions --output_dir ${output_dir} --data_dir ${dset_dir} --force"

        fi

    fi
done
