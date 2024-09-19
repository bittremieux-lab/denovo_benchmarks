#!/bin/bash
. /home/novor/.bashrc

# Start novor server
/algo/start-novor-server.sh

# Get dataset property tags
DSET_TAGS=$(python3 /algo/base/dataset_tags_parser.py --dataset "$@")
# Parse tags and set individual environment variables for each of them
# (variable names are identical to tag names
#  -- check DatasetTag values in dataset_config.py)
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Iterate through files in the dataset
for input_file in "$@"/*.mgf; do
    echo "Processing file: $input_file"
    input_basename=$(basename "$input_file")
    output_novor=${input_basename/.mgf/.denovo.csv}
    output_novorpt=${input_basename/.mgf/.novorpt.csv}
    
    # Convert input data to model format
    python3 input_mapper.py \
        --input_path "$input_file" \
        --output_path "$input_basename"

    # for the particular dataset properties
    if  [[ -v phosphorylation && $phosphorylation -eq 1 ]]; then
        echo "Using phosphorylation tag"
        pho="-phospho"
    fi

    if  [[ -v nontryptic && $nontryptic -eq 1 ]]; then
        echo "Using non-tryptic model."
        param_file="param-nonspecific${pho}.txt"
    else
        echo "Using general model."
        param_file="param-trypsin${pho}.txt"
    fi

    echo "Using parameter file: $param_file"
    #run novor
    /home/novor/bin/novor -f --topnout 12 -p "$param_file" -o "$output_novor" "$input_basename"
    #run novorpt
    /home/novor/bin/novorpt -f -d "$output_novor"
done

# Stop novor server
/algo/stop-novor-server.sh

# Convert predictions to the general output format
python3 output_mapper.py --output_dir="."
