#!/bin/bash

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
# Parse tags and set individual environment variables for each of them
# (variable names are identical to tag names
#  -- check DatasetTag values in dataset_config.py)
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Output file header
echo -e "sequence\tscore\taa_scores\tspectrum_id" > outputs.tsv

# Iterate through files in the dataset
for input_file in "$@"/*.mgf; do

    echo "Processing file: $input_file"

    # Convert input data to model format
    python input_mapper.py \
        --input_path "$input_file" \
        --output_path ./input_data.mgf

    cd DePS4DenovoBenchmarks
    python main.py --spectrum ../input_data.mgf

    # Collect predictions (from algorithm output file output.tsv)
    tail -n+2 output.tsv >> /algo/outputs.tsv
    cd /algo
done

# Convert predictions to the general output format
python ./output_mapper.py --output_path=outputs.tsv
