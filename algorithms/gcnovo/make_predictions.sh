#!/bin/bash

# Get dataset property tags
#DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
## Parse tags and set individual environment variables for each of them
## (variable names are identical to tag names
##  -- check DatasetTag values in dataset_config.py)
#while IFS='=' read -r key value; do
#    export "$key"="$value"
#done <<< "$DSET_TAGS"

# Iterate through files in the dataset
mkdir -p /algo/input_data
mkdir -p /algo/outputs

# Output file header
echo -e "spectrum_id,feature_area,sequence,score,predicted_position_score,precursor_mz,precursor_charge,protein_access_id,scan_list_middle,scan_list_original,predicted_score_max,aa_scores" > /algo/outputs.csv

for input_file in "$@"/*.mgf; do
    echo "Processing file: $input_file"

    # cp $input_file /algo/input_data/input_data.mgf
    # Convert input data to model format
    python input_mapper.py \
        --input_path "$input_file" \
        --output_path /algo/input_data/input_data.mgf

    python gcnovo_main.py \
        --denovo_input_spectrum_file="/algo/input_data/input_data.mgf" \
        --denovo_input_feature_file="/algo/input_data/input_data.mgf.csv" \
        --denovo_output_file="/algo/outputs/output"

    # Collect predictions
    tail -n+2 "/algo/outputs/output[benchmark].csv" >> /algo/outputs.csv
done
