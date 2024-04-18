import os
import streamlit as st

RESULTS_DIR = "./results"

st.set_page_config(layout="wide")

st.title("*De novo* Benchmarks")
st.divider()

st.header("Info")
st.divider()

st.header("Benchmarking results")

dataset_name = "9 species dataset: Apis-mellifera"
st.subheader(dataset_name)

col1, col2 = st.columns(2, gap="medium")

with col1:
    # st.write("Peptide precision & coverage")
    st.image(os.path.join(RESULTS_DIR, "Peptide_precision_coverage.png"))

with col2:
    # st.write("Amino acid precision & coverage")
    st.image(os.path.join(RESULTS_DIR, "AA_precision_coverage.png"))

