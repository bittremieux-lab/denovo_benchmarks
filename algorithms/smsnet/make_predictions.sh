#!/bin/bash

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
# Parse tags and set individual environment variables for each of them
# (variable names are identical to tag names
#  -- check DatasetTag values in dataset_config.py)
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Create output folder
mkdir /algo/data
mkdir /algo/data/smsnet_temp
mkdir /algo/data/smsnet_temp/mgf
mkdir /algo/data/smsnet_temp/mgf_output

rm '/algo/data/smsnet_temp/mgf/*'
rm '/algo/data/smsnet_temp/mgf_output/*'

# Iterate through files in the dataset
for input_file in "$@"/*.mgf; do

    echo "Processing file: $input_file"

    FNAME="${input_file##*/}"
    FPREFIX="${FNAME%.*}"

    # Convert input data to model format
    python input_mapper.py \
        --input_path "$input_file" \
        --output_path /algo/data/smsnet_temp/mgf/$FNAME

    # [Optionally] use tag variables to specify de novo algorithm
    # for the particular dataset properties
    if  [[ -v phosphorylation && $phosphorylation -eq 1 ]]; then
        echo "Using phospho model (phosphotylation of STY, oxidation of M)."
        cd /algo/phospho/
        python run.py --model_dir ./model --inference_input_file /algo/data/smsnet_temp/mgf/$FNAME --rescore
        MODE='p-mod'
    # use standard mode
    else
        echo "Using standard model (oxidation of M only)."
        cd /algo/standard/
        python run.py --model_dir ./model --inference_input_file /algo/data/smsnet_temp/mgf/$FNAME --rescore
        MODE='m-mod'
    fi

done

# Go back to writable location
cd /algo/data

# Compile predictions 
python /algo/create_denovo_report.py /algo/data/smsnet_temp/mgf_output /algo/data/smsnet_temp/mgf $MODE

# Convert predictions to the general output format
python /algo/output_mapper.py --output_path=/algo/data/smsnet_temp/mgf_$MODE"_fdr5.tsv"

# Store predictions in /algo (predictions are expected to be written to /algo/outputs.csv)
cp /algo/data/outputs.csv /algo/outputs.csv
