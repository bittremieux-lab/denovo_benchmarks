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

    To configure the model for specific data properties (e.g. non-tryptic data, data from a particular instrument, etc.), please use **dataset tags**. 
    Current set of tags can be found in the `DatasetTag` in [dataset_config.py](https://github.com/PominovaMS/denovo_benchmarks/blob/main/dataset_config.py) and includes `nontryptic`, `timstof`, `waters`, `sciex`.
    Example usage can be found in `algorithms/base/make_predictions_template.sh`.

- **`input_mapper.py`** — python script to convert input data 
from its original representation (**input format**) to the format expected by the algorithm.

    **Input format**
    - Input: a dataset folder with separate .mgf files containing MS spectra.
    - Keys order for a spectrum in .mgf file:  
    `[TITLE, RTINSECONDS, PEPMASS, CHARGE]`


- **`output_mapper.py`** — python script to convert the algorithm output to the common **output format**.

    **Output format**
    - .csv file (with `sep=","`)
    - must contain columns:
        - `"sequence"` — predicted peptide sequence, written in the predefined **output sequence format**
        - `"score"` — *de novo* algorithm "confidence" score for a predicted sequence
        - `"aa_scores"` — per-amino acid scores, if available. If not available, the whole peptide `score` will be used as a score for each amino acid.
        - `"spectrum_id"` — information to match each prediction with its ground truth sequence.  
            `{filename}:{index}` string, where  
            `filename` — name of the .mgf file in a dataset,  
            `index` —  index (0-based) of each spectrum in an .mgf file.
        
    
    - **Output sequence format**
        - 20 amino acid tokens:  
        `G, A, S, P, V, T, C, L, I, N, D, Q, K, E, M, H, F, R, Y, W`
        - Amino acids with post-translational modifications (PTMs) are written in 
        **[ProForma](https://github.com/HUPO-PSI/ProForma/tree/master) format** with **Unimod accession codes** for PTMs:  
        `C[UNIMOD:4]` for Cysteine Carbamidomethylation, `M[UNIMOD:35]` for Methionine Oxidation, etc.
        - N-terminus and C-terminus modifications, if supported by the algorithm, are also written in **ProForma notation** with **Unimod accession codes**:  
        `[UNIMOD:xx]-PEPTIDE-[UNIMOD:yy]`


## Running the benchmark

To run the benchmark locally:

1. **Clone the repository**:
    ```bash
    git clone https://github.com/PominovaMS/denovo_benchmarks.git
    cd denovo_benchmarks
    ```

2. **Build containers for algorithms and evaluation**:
    To build all apptainer images, make sure you have [apptainer installed](https://apptainer.org/docs/user/main/quick_start.html#installation). Then run:

    ```bash
    chmod +x build_apptainer_images.sh
    ./build_apptainer_images.sh
    ```

    This will build the apptainer images for all algorithms and the evaluation apptainer image.

    If an apptainer image already exists, the script will ask if you want to rebuild it.

    ```bash
    A .sif image for casanovo already exists. Force rebuild? (y/N) 
    ```

    If a container is missing, that algorithm will be skipped during benchmarking. We don't share or store containers publicly yet due to ongoing development and their large size.

3. **Run benchmark on a dataset**:
    ```bash
    ./run.sh /path/to/dataset/dir
    ```
    Example:
    ```bash
    ./run.sh sample_data/PXD004424
    ```


## Input data structure

The benchmark expects input data to follow a specific folder structure. 

- Each dataset is stored in a separate folder with unique name.
- **Spectra** are stored as `.mgf` files inside the `mgf/` subfolder.
- **Ground truth labels** (PSMs found via database search) are contained in `labels.csv` file within each dataset folder.

Below is an example layout for our evaluation datasets stored on the HPC:

```
datasets/
    PXD004424/
        labels.csv
        mgf/
            151009_exo3_1.mgf
            151009_exo3_2.mgf
            151009_exo3_3.mgf
            ...
    PXD004947/
        labels.csv
        mgf/...
    PXD004948/
        labels.csv
        mgf/...
    PXD004325/
        labels.csv
        mgf/...
    ...
```

Note that algorithm containers only get as input the `/mgf` subfolder with spectra files and **do not** have access to the `labels.csv` file. 
Only the evaluation container accesses the `labels.csv` file to evaluate algorithm predictions.


## Running Streamlit dashboard locally:
To view the Streamlit dashboard for the benchmark locally, run:
```bash
# If Streamlit is not installed
pip install streamlit

streamlit run dashboard.py
```

The dashboard reads the benchmark results stored in the `results/` folder.