#!/bin/bash
#. /home/novor/.bashrc

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
    output_novor=${input_basename/.mgf/.novorai.csv}

    # Convert input data to model format
    python3 input_mapper.py \
        --input_path "$input_file" \
        --output_path "$input_basename"

    # for the particular dataset properties
    PTM_LIST=()
    if  [[ -v tmt && $tmt -eq 1 ]]; then
    	PTM_LIST+=("TMT6 (K)")
    	PTM_LIST+=("TMT6 (N-term)")
    fi

    if  [[ -v silac && $silac -eq 1 ]]; then
        PTM_LIST+=("Silac-Lys4")
        PTM_LIST+=("Silac-Arg6")
        PTM_LIST+=("Silac-Lys8")
        PTM_LIST+=("Silac-Arg10")
    fi

    if  [[ -v phosphorylation && $phosphorylation -eq 1 ]]; then
    	PTM_LIST+=("Phospho (STY)")
    fi

    if  [[ -v oxidation && $oxidation -eq 1 ]]; then
    	PTM_LIST+=("Oxidation (M)")
    fi

    if  [[ -v formylation && $formylation -eq 1 ]]; then
    	PTM_LIST+=("Formyl (N-term)")
    	PTM_LIST+=("Formyl (KST)")
    fi

    if  [[ -v acetylation && $acetylation -eq 1 ]]; then
    	PTM_LIST+=("Acetyl (K)")
    	PTM_LIST+=("Acetyl (N-term)")
    fi

    if  [[ -v methylation && $methylation -eq 1 ]]; then
    	PTM_LIST+=("Methyl (DE)")
    fi

    if  [[ -v carbamidomethylation && $carbamidomethylation -eq 1 ]]; then
    	PTM_LIST+=("Carbamidomethyl (C)")
    fi

    if [ ${#PTM_LIST[@]} -eq 0 ]; then
    	PTM_LIST+=("Carbamidomethyl (C)")
    	PTM_LIST+=("Oxidation (M)")
    fi

    PTM_STRING=$(IFS=', '; echo "${PTM_LIST[*]}")
    PTM_STRING="variableModifications = "$PTM_STRING

    param_file="param-input.txt"
    cp param.txt $param_file
    echo $PTM_STRING >> $param_file

    echo "Using parameter file: $param_file"
    echo "---"
    cat $param_file
    echo "---"

    NOVOR_BIN=/home/novor/novorai/novoraidenovo
    if echo "$input_file" | grep -q -e "mAb" -e "herceptin"; then
       NOVOR_BIN=/home/novor/novorai/novoraidenovo-ab
    fi
    if  [[ -v timstof && $timstof -eq 1 ]]; then
    	NOVOR_BIN=/home/novor/novorai/novoraidenovo-timstof
    fi

    #run novor
    $NOVOR_BIN -p "$param_file" -o "$output_novor" -i "$input_basename"
done

# Convert predictions to the general output format
python3 output_mapper.py --output_dir="."
