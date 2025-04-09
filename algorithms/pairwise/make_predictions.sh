#!/bin/bash
conda init
source /opt/conda/etc/profile.d/conda.sh
conda activate pairwise_env

if command -v nvidia-smi &> /dev/null && nvidia-smi > /dev/null 2>&1; then
    device=0
    accelerator=cuda
else
    device=-1
    accelerator=cpu
fi
echo "Using device $device with accelerator $accelerator."

cd /algo/pairwise

python src/data/create_lance.py --input /algo/"$@"/mgf --output /algo/data.lance

#TODO: add checkpoint
python src/main.py --config=configs/master_bm.yaml --num_workers=0 --accelerator=$accelerator --downstream_root_dir=/algo/data.lance #--downstream_weights=checkpoint.ckpt 

# Placeholder for output mapper:
# python output_mapper.py --input predictions_table.mzTab --output final_outputs.csv
