import os
import streamlit as st
import streamlit.components.v1 as components

RESULTS_DIR = "./results"
PLOT_HEIGHT = 600
PLOT_WIDTH = int(600 * 1.2)

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
    path_to_html = os.path.join(RESULTS_DIR, "Peptide_precision_coverage.html")
    with open(path_to_html, 'r') as f: 
        html_data = f.read()
    components.html(html_data, width=PLOT_WIDTH, height=PLOT_HEIGHT, scrolling=False)

with col2:
    path_to_html = os.path.join(RESULTS_DIR, "AA_precision_coverage.html")
    with open(path_to_html, 'r') as f: 
        html_data = f.read()
    components.html(html_data, width=PLOT_WIDTH, height=PLOT_HEIGHT, scrolling=False)
