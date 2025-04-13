#!/bin/bash

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"


# Run de novo algorithm on the input data
echo "Using general model."
casanovo sequence -c config.yml -o outputs.mztab "$@"/*.mgf

# Convert predictions to the general output format
python ./output_mapper.py --output_path=outputs.mztab
