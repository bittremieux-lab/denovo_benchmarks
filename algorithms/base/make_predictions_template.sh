#!/bin/bash

# List and sort files in the dataset
read -a sorted_files -d "\n" < <(ls -1 "$@"/*.mgf | sort)

# Iterate through sorted files in the dataset
for file_idx in "${!sorted_files[@]}"; do

    input_file=${sorted_files[file_idx]}
    echo "Processing file $file_idx: $input_file"

    # Convert input data to model format
    python input_mapper.py \
        --input_path "$input_file" \
        --file_i "$file_idx" \
        --output_path ./input_data.mgf

    # Run de novo algorithm on the input data
    python ...
done

# Convert predictions to the general output format
python output_mapper.py --output_path=...
