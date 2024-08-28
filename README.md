TODO: update for Apptainer (Singularity) containers.

# Benchmarking de novo peptide sequencing algorithms

## Adding a new algorithm

Make a pull request to add your algorithm to the benchmarking system.

Add your algorithm in the `denovo_benchmarks/algorithms/algorithm_name` folder by providing  
`container.def`, `make_predictions.sh`, `input_mapper.py`, `output_mapper.py` files.  
Detailed files descriptions are given below.  

Templates for each file implementation can be found in the 
`algorithms/base/` [folder](https://github.com/PominovaMS/denovo_benchmarks/tree/main/algorithms/base).  
It also includes the `InputMapperBase` and `OutputMapperBase` base classes for implementing input and output mappers.  
For examples, you can check 
[Casanovo](https://github.com/PominovaMS/denovo_benchmarks/tree/main/algorithms/casanovo) 
and [DeepNovo](https://github.com/PominovaMS/denovo_benchmarks/tree/main/algorithms/deepnovo) implementations. 

- **`container.def`** — definition file of the [Apptainer](https://apptainer.org/docs/user/main/definition_files.html) 
container image that creates environment and installs dependencies required for running the algorithm.
    
- **`make_predictions.sh`** — bash script to run the *de novo* algorithm on the input dataset 
(folder with MS spectra in .mgf files) and generate an output file with per-spectrum peptide predictions.  
    **Input**: path to a dataset folder containing .mgf files with spectra data  
    **Output**: output file (in a common output format) containing predictions for all spectra in the dataset

- **`input_mapper.py`** — python script to convert input data 
from its original representation (**input format**) to the format expected by the algorithm.

    **Input format**
    - Input: a dataset folder with separate .mgf files containing MS spectra with annotations.
    - Keys order for a spectrum in .mgf file:  
    `[TITLE, PEPMASS, RTINSECONDS, CHARGE, SCANS, SEQ]`
    - Annotation (peptide sequence) format:
        - `SEQ=PEPTIDE` string
        - 20 amino acid tokens:  
        `G, A, S, P, V, T, C, L, I, N, D, Q, K, E, M, H, F, R, Y, W`
        - PTMs represented in the [ProForma](https://github.com/HUPO-PSI/ProForma/tree/master) notation:  
        `C[UNIMOD:4]`, `M[UNIMOD:35]` , etc.
    - If the algorithm uses the ground truth sequence from the input file (e.g. for evaluation), 
    tokens in the annotation sequence may need conversion to the expected format  
    (e.g. `'M[UNIMOD:35]' → 'Mmod'` or `'M[+15.999]' → 'Mmod'` )  
` `  
- **`output_mapper.py`** — python script to convert the algorithm output to the common **output format**.

    **Output format**
    - .csv file (with `sep=","`)
    - must contain columns:
        - `"sequence"` — predicted peptide sequence, written in the predefined **output sequence format**
        - `"score"` — *de novo* algorithm "confidence" score for a predicted sequence
        - `"aa_scores"` — per-amino acid scores, if available. If not available, the whole peptide `score` will be used as a score for each amino acid.
        - `"spectrum_id"` — information to match each prediction with its ground truth sequence.  
            `F{filename}:{scan_id}` string, where  
            `filename` — name of the .mgf file in a dataset,  
            `scan_id` — scan id of each spectrum (SCANS= in an .mgf file)  
        ` `  
    
    - **Output sequence format**
        - 20 amino acid tokens:  
        `G, A, S, P, V, T, C, L, I, N, D, Q, K, E, M, H, F, R, Y, W`
        - `C` is written **without** Carbamidomethyl modification 
        (even if it was considered as modified on the prediction stage).
        - Amino acids with post-translational modifications (PTMs) are written in 
        **[ProForma](https://github.com/HUPO-PSI/ProForma/tree/master) format** **Delta mass notation**:  
        `M[+15.999]`, etc.
        - N-terminus and C-terminus modifications, if supported by the algorithm, are also written in **ProForma Delta mass notation**:  
        `[+n_mod]-PEPTIDE-[+c_mod]`
