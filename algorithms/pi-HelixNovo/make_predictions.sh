#!/bin/bash

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Choose the device  
if command -v nvidia-smi &> /dev/null
then
    nvidia-smi > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        device=0 # GPU:0
    else
        device=-1 # CPU
    fi
else
    device=-1 # CPU
fi

# Use tag variables to specify de novo algorithm
# for the particular dataset properties
conda activate main_env
cd /algo
echo "Processing mgf files:"
python pi-HelixNovo/main.py --mode=denovo --config=pi-HelixNovo/config.yaml --gpu=$device --output=outputs.csv --peak_path="$@"/*.mgf --model=pi-HelixNovo/MSV000081142-epoch-5-step-800000.ckpt
