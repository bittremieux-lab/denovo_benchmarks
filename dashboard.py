import os
import streamlit as st
import streamlit.components.v1 as components
from datasets_info import DATASETS


RESULTS_DIR = "results"
PLOT_HEIGHT = 600
PLOT_WIDTH = int(PLOT_HEIGHT * 1.2)

st.set_page_config(layout="wide")

tab1, tab2 = st.tabs(["Main", "Adding an algorithm"])

with tab1:
    st.title("*De novo* Benchmarks")
    # st.divider()

    st.header("Info")
    st.markdown(
        """
        *De novo* peptide sequencing identifies peptides from mass spectrometry data
        without a reference proteome database. 
        It is essential for discovering novel peptides and in applications 
        where database construction is complicated (immunopeptidomics, metabolomics).
        However, evaluation and comparison of multiple existing methods is challenging
        due to the lack of standardized metrics and universal test datasets.
       
        This project aims to provide a unified framework for the comprehensive benchmarking 
        of *de novo* peptide sequencing algorithms.
        It performs evaluation on both publicly available and internal datasets, 
        encompassing different species, variations in cleavage enzymes (tryptic and non-tryptic data), 
        and various post-translational modifications (PTMs), 
        as well as specific proteomics subfields such as immunopeptidomics. 
        Each algorithm is run in an isolated environment using Apptainer containers 
        to ensure reproducibility and consistency.

        Current benchmarking results are represented below. 
        """
    )
    # st.divider()

    st.header("Benchmarking results")
    st.markdown(
        """
        Results are split by the evaluation dataset. 

        Ground truth labels were obtained by running database search 
        with [MSFragger](https://msfragger.nesvilab.org/) and subsequence rescoring 
        with [Percolator](http://percolator.ms/) using features from Prosit-predicted spectra.  
        Unless otherwise specified, PSMs with estimated FDR (Percolator q-value) **<1%** 
        were selected as ground truth peptides.

        The plots demonstrate performance of *de novo* peptide sequencing algorithms according to **2** metrics:
        - Peptide prediction precision and coverage
        - Amino acids prediction precision and coverage

        Click on algorithms on the plots to select a subset of algorithms to display. 
        """
    )

    datasets = os.listdir(RESULTS_DIR)
    for dataset_name in datasets:
        st.subheader(dataset_name)
        if dataset_name in DATASETS:
            dataset_descr = DATASETS[dataset_name]
            st.text(dataset_descr)

        col1, col2 = st.columns(2, gap="medium")

        with col1:
            path_to_html = os.path.join(
                RESULTS_DIR, dataset_name, "peptide_precision_coverage.html"
            )
            with open(path_to_html, "r") as f:
                html_data = f.read()
            components.html(
                html_data, width=PLOT_WIDTH, height=PLOT_HEIGHT, scrolling=False
            )

        with col2:
            path_to_html = os.path.join(
                RESULTS_DIR, dataset_name, "AA_precision_coverage.html"
            )
            with open(path_to_html, "r") as f:
                html_data = f.read()
            components.html(
                html_data, width=PLOT_WIDTH, height=PLOT_HEIGHT, scrolling=False
            )


with tab2:
    st.header("Adding a new algorithm")
    # st.divider()

    st.markdown(
        """
        Make a pull request on [Github](https://github.com/PominovaMS/denovo_benchmarks) 
        to add your algorithm to the benchmarking system.
        """
    )

    main, sidebar = st.columns([2, 1])

    with main:
        st.markdown(
            """
            Add your algorithm in the `denovo_benchmarks/algorithms/algorithm_name` folder by providing  
            `container.def`, `make_predictions.sh`, `input_mapper.py`, `output_mapper.py` files.  
            Detailed files descriptions are given below.  
            
            Templates for each file implementation can be found in the 
            `algorithms/base/` [folder](https://github.com/PominovaMS/denovo_benchmarks/tree/main/algorithms/base).  
            It also includes the `InputMapperBase` and `OutputMapperBase` base classes for implementing input and output mappers.  
            For examples, you can check 
            [Casanovo](https://github.com/PominovaMS/denovo_benchmarks/tree/main/algorithms/casanovo) 
            and [DeepNovo](https://github.com/PominovaMS/denovo_benchmarks/tree/main/algorithms/deepnovo) implementations. 
            """
        )

        st.subheader("Files description")
        st.markdown(
            """
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
                    - `"scans"` — information to match each prediction with its ground truth sequence.  
                        `F{file_i}:{scan_id}` string, where  
                        `file_i` — number of the .mgf file in a dataset (lexicographically sorted),  
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
            """
        )

        # expander = st.expander("**`input_mapper.py`** — python script to convert input data from its original representation (**input format**) to the format expected by the algorithm.")
        # expander.markdown(
        #     """
        #     **Input format**
        #     - Input: a dataset folder with separate .mgf files containing MS spectra with annotations.
        #     - Keys order for a spectrum in .mgf file:  
        #     `[TITLE, PEPMASS, RTINSECONDS, CHARGE, SCANS, SEQ]`
        #     - Annotation (peptide sequence) format:
        #         - `SEQ=PEPTIDE` string
        #         - 20 amino acid tokens:  
        #         `G, A, S, P, V, T, C, L, I, N, D, Q, K, E, M, H, F, R, Y, W`
        #         - PTMs represented in the [ProForma](https://github.com/HUPO-PSI/ProForma/tree/master) notation:  
        #         `C[UNIMOD:4]`, `M[UNIMOD:35]` , etc.
        #     - If the algorithm uses the ground truth sequence from the input file (e.g. for evaluation), 
        #     tokens in the annotation sequence may need conversion to the expected format  
        #     (e.g. `'M[UNIMOD:35]' → 'Mmod'` or `'M[+15.999]' → 'Mmod'` )  
        #     """
        # )

        # expander = st.expander("**`output_mapper.py`** — python script to convert the algorithm output to the common **output format**.")
        # expander.markdown(
        #     """
        #     **Output format**
        #     - .csv file (with `sep=","`)
        #     - must contain columns:
        #         - `"sequence"` — predicted peptide sequence, written in the predefined **output sequence format**
        #         - `"score"` — *de novo* algorithm "confidence" score for a predicted sequence
        #         - `"aa_scores"` — per-amino acid scores, if available. If not available, the whole peptide `score` will be used as a score for each amino acid.
        #         - `"scans"` — information to match each prediction with its ground truth sequence. One of two columns must be provided.
        #             - `"scans"` — `F{file_i}:{scan_id}` string, where  
        #             `file_i` — number of the .mgf file in a dataset (lexicographically sorted),  
        #             `scan_id` — scan id of each spectrum (SCANS= in an .mgf file)  
        #         ` `  
            
        #     - **Output sequence format**
        #         - 20 amino acid tokens:  
        #         `G, A, S, P, V, T, C, L, I, N, D, Q, K, E, M, H, F, R, Y, W`
        #         - `C` is written **without** Carbamidomethyl modification 
        #         (even if it was considered as modified on the prediction stage).
        #         - Amino acids with post-translational modifications (PTMs) are written in 
        #         **[ProForma](https://github.com/HUPO-PSI/ProForma/tree/master) format** **Delta mass notation**:  
        #         `M[+15.999]`, etc.
        #         - N-terminus and C-terminus modifications, if supported by the algorithm, are also written in **ProForma Delta mass notation**:  
        #         `[+n_mod]-PEPTIDE-[+c_mod]`
        #     """
        # )