# PTM Model

This repository contains the code and weights for a PrimeNovo model finetuned on "Phosphorylation (+79.97)" data. The model was adapted from the original PrimeNovo (trained on MassiveKB) by finetuning it on the 2020-Cell-LUAD dataset. The training data was enriched with phosphorylation data, so the model is designed for detecting phosphorylation PTMs.

## Environment Setup

Set up the environment exactly as before to ensure compatibility.

## Model Overview

The finetuned model has been modified to predict one additional token compared to the original PrimeNovo. For convenience, we added an extra PTM token `"B"` in the `config.yaml` located in this directory. When the model outputs a `"B"`, it represents "Phosphorylation (+79.97)".

You can also train or finetune your own PrimeNovo-PTM model by reassigning the token `"B"` to a different PTM of your choice. **Important:** If you change the token `"B"` to another value, ensure that your training data reflects this change so that the target PTM matches the modified token during training.

## Modifying the PTM Token

If you decide to change the PTM from "Phosphorylation (+79.97)" to a different modification, update the PTM molecular weight in the following files:

1. in `mass_con.py`: change variable `mass_b =79.9663` to whatever new PTM mole weight it has.

2. in `model.py`: change variable `mass_b =79.9663` to whatever new PTM mole weight it has.

3. in `transformers.py`: change variable `mass_b =79.9663` to whatever new PTM mole weight it has.

4. in `config.yaml`: change "B": 79.9663 to whatever new PTM mole weight it has.

And that's it, you can then finetune on your own dataset with additional `1` PTM, just replace that PTM with letter "B" anywhere it appears in your training/testing data.

## Run Instructions

**Note!!!!!!!!!!!!!!!!!!:** All the following steps should be performed under the main directory: `pi-PrimeNovo-PTM`. Do **not** use `cd ..` !!!!!!!!!!!!!!!!!!!

### Step 1: Download Required Files

To evaluate the mgf containing Phosphated PTM:

1. **Model Checkpoint**: [PTM_Phosphorylation.ckpt](https://drive.google.com/file/d/1YcF9VNE1gFF8T0EfwcFb7v1tiKw25-ai/view?usp=share_link)


**Note:** If you are using a remote server, you can use the `gdown` package to easily download the content from Google Drive to your server disk.

### Step 2: Choose the Mode

The `--mode` argument can be set to either:

- `eval`: Use this mode when evaluating data with a labeled dataset.
- `denovo`: Use this mode for de novo analysis on unlabeled data.

**Important**: Select `eval` only if your data is labeled.

### Step 3: Run the Commands

Execute the following command in the terminal:

```bash
python -m PrimeNovo.PrimeNovo --mode=eval --peak_path=./PTM_test.mgf --model=./PTM_Phosphorylation.ckpt
```

This automatically uses all GPUs available in the current machine.

### Step 4: analyze the output

We include a sample running output ```./output.txt```. The performance for evaluation will be reported at the end of the output file.

If you are using ```denovo``` mode, you will get a ```denovo.tsv``` file under the current directory. The file has the following structure:

| label | prediction | charge | score |
| --- | --- | --- | --- |
| Title in MGF document | Sequence in [ProForma](https://doi.org/10.1021/acs.jproteome.1c00771) notation| Charge, as a number | Confidence score as number in range 0 and 1 using scientific notation |

The example below contains two peptides predicted based on some given spectrum:

```tsv
label	prediction	charge	score
MS_19321_2024_02_DDA	ATTALP	2	0.99
MS_19326_2024_02_DDA	TAM[+15.995]TR	2	0.87
```

NOTE: "B" amino acid will stands for the additional PTM added to this code, in this case it is Phosphorylation.

## Citation

Bibtex：

```bibtex
@article{zhang2025pi,
  title={$\pi$-PrimeNovo: an accurate and efficient non-autoregressive deep learning model for de novo peptide sequencing},
  author={Zhang, Xiang and Ling, Tianze and Jin, Zhi and Xu, Sheng and Gao, Zhiqiang and Sun, Boyan and Qiu, Zijie and Wei, Jiaqi and Dong, Nanqing and Wang, Guangshuai and others},
  journal={Nature Communications},
  volume={16},
  number={1},
  pages={267},
  year={2025},
  publisher={Nature Publishing Group UK London}
}
```

Standard Citation:

```
Zhang, X., Ling, T., Jin, Z. et al. π-PrimeNovo: an accurate and efficient non-autoregressive deep learning model for de novo peptide sequencing. Nat Commun 16, 267 (2025). https://doi.org/10.1038/s41467-024-55021-3
```
