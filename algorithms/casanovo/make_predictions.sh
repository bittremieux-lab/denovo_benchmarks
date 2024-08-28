#!/bin/bash

# Get dataset property tags
DSET_TAGS=$(python base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export $key=$value
done <<< "$DSET_TAGS"


# Use tag variables to specify de novo algorithm
# for the particular dataset properties
if  [[ -v nontryptic && $nontryptic -eq 1 ]]; then
    # Run de novo algorithm on the input data
    echo "Using non-tryptic model."
    casanovo sequence -c config.yml -o outputs.mztab "$@"/*.mgf --model ./casanovo_nontryptic.ckpt
else
    # Run de novo algorithm on the input data
    echo "Using general model."
    casanovo sequence -c config.yml -o outputs.mztab "$@"/*.mgf
fi

# Convert predictions to the general output format
python ./output_mapper.py --output_path=outputs.mztab --input_dir="$@"
