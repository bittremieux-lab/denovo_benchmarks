#!/bin/bash
dset_dir="$1"
algorithm="$2"
spectra_dir="$dset_dir/mgf"
output_root_dir="./outputs"

dset_name=$(basename "$dset_dir")
output_dir="$output_root_dir/$dset_name"

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

            output_file="$output_dir/${algorithm_name}_output.csv"
            echo "Output file: $output_file"

            # Augment algorithm predictions with RT and SA (if not already present)
            echo "AUGMENT PREDICTIONS"
            apptainer exec --fakeroot --env-file .env "evaluation.sif" \
                bash -c "python -m evaluation.augment_predictions --output_dir ${output_dir} --data_dir ${dset_dir} --algo_name ${algorithm_name} --force"

        fi

    fi
done
