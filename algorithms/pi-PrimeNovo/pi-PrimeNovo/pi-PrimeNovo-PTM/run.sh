#!/bin/bash

#SBATCH --job-name=xiang
#SBATCH --partition=beam
#SBATCH --quotatype=reserved
#SBATCH --gres=gpu:8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8




PTM_phos=../NON790/epoch=134-step=54400.ckpt
data_PTM=~/test_LUAD_NF

srun python -m PrimeNovo.PrimeNovo --mode=eval --model=$PTM_phos --peak_path=$data_PTM