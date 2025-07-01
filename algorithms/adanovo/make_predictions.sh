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
echo -e "sequence,PSM_ID,accession,unique,database,database_version,search_engine,score,modifications,retention_time,charge,exp_mass_to_charge,calc_mass_to_charge,spectrum_id,pre,post,start,end,aa_scores" > all_outputs.csv

# Iterate through files in the dataset
for input_file in "$@"/*.mgf; do

    echo "Processing file: $input_file"

    # Convert input data to model format
    python ./adanovo_v1/adanovo.py --mode=denovo --model=./ckpt_1.ckpt --peak_path="$input_file" --config=./adanovo_v1/config.yaml --output=./demo.mztab

    # Convert predictions to the general output format
    python ./output_mapper.py --output_path=./demo.mztab --input_path="$input_file"

    # Collect predictions
    tail -n+2 outputs.csv >> /algo/all_outputs.csv
done

# Rename final output
mv /algo/all_outputs.csv /algo/outputs.csv
