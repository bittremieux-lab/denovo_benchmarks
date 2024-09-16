# SMSNet

## Sequence-Mask-Search framework and SMSNet model for de novo peptide sequencing.

The pre-trained model can be downloaded from [FigShare](https://figshare.com/articles/SMSNet_s_trained_models/8259122).

SMSNet's predicted amino acid sequences for a public HLA peptidome dataset (MassIVE accession MSV000080527) and phosphoproteome dataset (PRIDE accession PXD009227) can be found on [FigShare](https://figshare.com/articles/SMSNet_s_predictions_for_HLA_peptidome_and_human_phosphoproteome_datasets/8259134).

The preprint can be found on [bioRxiv](http://biorxiv.org/cgi/content/short/667527v1).

## Dependencies
This project uses Python 3.5.2, with the following lib dependencies:
* [Tensorflow 1.4](https://www.tensorflow.org/) (compatible up to 1.11)
* [Keras 2.2.4](https://keras.io/)

A list of all python packages can be found in ```requirement.txt```


## Instructions
### Decode
```
python run.py --model_dir <model_directory> --inference_input_file <input_file_directory/input_file.mgf> --rescore
```
```<model_directory>``` is the directory of the model (can be downloaded from the link above).

```<input_file_directory/input_file.mgf>``` is the path to the input file.

Using --rescore flag will generate another probability file with suffix “_rescore” in the same directory. The output will be in “<input_file_directory>_output/”.

Other options can be found in "run.py".

Model parameters (including possible amino acids) can be found in "nmt/input_config.py".

### Note
* In order to generate the report file, the ```TITLE``` lines in .mgf file must end with "scan=<number>".
* To switch between m-mod and p-mod model, the following changes are needed (default: p-mod):
    * ```tgt_vocab_size``` (24 for m_mod / 27 for p-mod) and ```tgt_vocab_file``` in ```run.py``` line 61-62.
    * Comment/uncomment the possible vocab in ```inverse_vocab``` in ```nmt/input_config.py``` accordingly ('s', 't', 'y' at line 65, 67, 71).
    * Select the corresponding ```AA_weight_table``` in function ```create_aa_tables()``` in ```nmt/input_config.py``` (by comment/uncomment line 169-174 or 176-180).


### Outputs
* For each input file, three output files will be generated in the output directory: ```<input_file>```, ```<input_file>_prob``` and, ```<input_file>_rescore```. They are the output sequences, probabilities of each amino acid, and probabilities of each amino acid after rescoring, respectively.
* The report summarizing the outputs in .tsv format will be in the same parent directory as the input.

## Example
Decoding “test_decode/test_file.mgf” with a model in "model/m_mod/".
```

python run.py --model_dir model/m_mod/ --inference_input_file test_decode/test_file.mgf --rescore
```
The report will be in test_decode_output.

The model will also produce three files:
```
test_decode_output/test_file (sequence)
test_decode_output/test_file_prob (score)
test_decode_output/test_file_rescore (score after rescoring with post-processing model)
```

## Database search
For database searching, change the file name in ```utils_masking/SMSNet_final_database_search.py```, then run ```python utils_masking/SMSNet_final_database_search.py```
    
    
## Acknowledgement
This code is based on TensorFlow Neural Machine Translation [GNMT](https://github.com/tensorflow/nmt). 
