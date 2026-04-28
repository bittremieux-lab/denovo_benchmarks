#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Get dataset property tags
DSET_TAGS=$(python3 /algo/base/dataset_tags_parser.py --dataset "$@")
# Parse tags and set individual environment variables for each of them
# (variable names are identical to tag names
#  -- check DatasetTag values in dataset_config.py)
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

add_ptms_for_tag() {
    local tag="$1"
    shift
    if [[ -v $tag && ${!tag} -eq 1 ]]; then
        PTM_LIST+=("$@")
    fi
}

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

    add_ptms_for_tag tmt \
        "TMT6 (K)" \
        "TMT6 (N-term)"

    add_ptms_for_tag silac \
        "Silac-Lys4" \
        "Silac-Arg6" \
        "Silac-Lys8" \
        "Silac-Arg10"

    add_ptms_for_tag phosphorylation \
        "Phospho (STY)"

    add_ptms_for_tag oxidation \
        "Oxidation (M)"

    add_ptms_for_tag formylation \
        "Formyl (N-term)" \
        "Formyl (KST)"

    add_ptms_for_tag acetylation \
        "Acetyl (K)" \
        "Acetyl (N-term)"

    add_ptms_for_tag methylation \
        "Methyl (DE)"

    add_ptms_for_tag carbamidomethylation \
        "Carbamidomethyl (C)"

    # 21 PTMs
    add_ptms_for_tag dimethylation \
        "Dimethyl (K)" \
        "Dimethyl (R)"

    add_ptms_for_tag trimethylation \
        "Trimethyl (K)"

    add_ptms_for_tag crotonylation \
        "Crotonyl (K)"

    add_ptms_for_tag hydroxyisobutyryl \
        "Hydroxyisobutyryl (K)"

    add_ptms_for_tag biotinylation \
        "Biotin (K)"

    add_ptms_for_tag propionyl \
        "Propionyl (K)"

    add_ptms_for_tag ubiquitination \
        "GlyGly (K)"

    add_ptms_for_tag butyryl \
        "Butyryl (K)"

    add_ptms_for_tag succinyl \
        "Succinyl (K)"

    add_ptms_for_tag malonylation \
        "Malonyl (K)"

    add_ptms_for_tag glutarylation \
        "Glutaryl (K)"

    add_ptms_for_tag citrullination \
        "Citrullination (R)"

    add_ptms_for_tag nitro \
        "Nitro (Y)"

    add_ptms_for_tag acetaldehyde \
        "Acetaldehyde_26 (HK)" \
        "Acetaldehyde_28 (HK)"

    add_ptms_for_tag iron_3h_substitution \
        "Iron3 (DE)"

    # Set default PTMs if empty
    if [ ${#PTM_LIST[@]} -eq 0 ]; then
    	PTM_LIST+=("Carbamidomethyl (C)")
    	PTM_LIST+=("Oxidation (M)")
    fi

    PTM_STRING=$(IFS=', '; echo "${PTM_LIST[*]}")
    PTM_STRING="variableModifications = "$PTM_STRING

    param_file="param-input.txt"
    cp param.txt $param_file
    echo $PTM_STRING >> $param_file

    ptm_file="ptms.txt"
    echo "Using ptm file: $ptm_file"
    echo "---"
    cat $ptm_file
    echo "---"

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
    $NOVOR_BIN -m "$ptm_file" -p "$param_file" -o "$output_novor" -i "$input_basename"
done

# Convert predictions to the general output format
python3 output_mapper.py --output_dir="."
