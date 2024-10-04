#!/bin/bash

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Output file head
echo "sequence,score,aa_scores,spectrum_id" > outputs.csv

# Iterate through files in the dataset
for input_file in "$@"/*.mgf; do

    echo "Processing file: $input_file"

    python denovo.py  --input "$input_file" --output decode.rst

    # Collect predictions
    cat decode.rst >> outputs.csv
done