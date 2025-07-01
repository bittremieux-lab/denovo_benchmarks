#!/bin/bash
conda init
source /opt/conda/etc/profile.d/conda.sh
conda activate pairwise_env

if command -v nvidia-smi &> /dev/null && nvidia-smi > /dev/null 2>&1; then
    device=0
    accelerator=gpu
else
    device=-1
    accelerator=cpu
fi
echo "Using device $device with accelerator $accelerator."

cd /algo/pairwise

python src/data/create_lance.py --input /algo/"$@"/ --output /algo/data.lance

#TODO: add checkpoint
python src/main.py --config=configs/master_bm.yaml --accelerator=$accelerator --downstream_root_dir=/algo/data.lance --downstream_weights=/algo/pairwise_mskb.ckpt 

cd /algo

# Placeholder for output mapper:
python output_mapper.py --input_dir /algo/"$@"/ --output_path pairwise/outs/logs/log/predictions_table.mzTab
