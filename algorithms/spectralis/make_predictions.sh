#!/bin/bash

# Make pyenvs available
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# Run Casanovo to obtain initial starting sequences
pyenv activate casanovo
casanovo_output="outputs.mztab"

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Use tag variables to specify de novo algorithm
# for the particular dataset properties
if  [[ -v nontryptic && $nontryptic -eq 1 ]]; then
    # Run de novo algorithm on the input data
    echo "Using non-tryptic model."
    casanovo sequence -c casanovo_config.yml -o $casanovo_output "$@"/*.mgf --model ./casanovo_nontryptic.ckpt
else
    # Run de novo algorithm on the input data
    echo "Using general model."
    casanovo sequence -c casanovo_config.yml -o $casanovo_output "$@"/*.mgf
fi

# Change to spectralis pyenv
pyenv deactivate
pyenv activate spectralis

# Write the initial starting sequences to the input MGFs
spectralis_data_dir="./seq_data"
mkdir -p $spectralis_data_dir
python /algo/intermediate_mapper.py --mztab_path $casanovo_output --mgf_in_dir "$@" --mgf_out_dir $spectralis_data_dir