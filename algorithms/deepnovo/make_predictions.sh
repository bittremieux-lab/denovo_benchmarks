#!/bin/bash

# Replace deepnovo_config.py with custom config file
cp config.py DeepNovo/deepnovo_config.py
rm -f outputs.tab

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
# Parse tags and set individual environment variables for each of them
# (variable names are identical to tag names
#  -- check DatasetTag values in config.py)
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Iterate through files in the dataset
for input_file in "$@"/*.mgf; do

    echo "Processing file: $input_file"

    # Convert input data to model format
    python input_mapper.py \
        --input_path "$input_file" \
        --output_path ./input_data.mgf

    # Run de novo algorithm on the input data
    python DeepNovo/deepnovo_main.py \
        --train_dir . --decode \
        --beam_search --beam_size 5

    # Collect predictions
    cat decode_output.tab >> outputs.tab
done

# Convert predictions to the general output format
python output_mapper.py --output_path=outputs.tab
