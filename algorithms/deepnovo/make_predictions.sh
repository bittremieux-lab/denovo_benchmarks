#!/bin/bash
# Build DeepNovo from source
cd DeepNovo/HighResolution/ && python setup.py build_ext --inplace && cd ../..
# Replace data_utils.py with config file
mv config.py DeepNovo/HighResolution/data_utils.py

# Convert input data to model format
python input_mapper.py "$@" ./input_data.mgf
# Run de novo algorithm on the input data
python DeepNovo/HighResolution/main.py --train_dir . --decode --beam_search
# Convert predictions to the general output format
python output_mapper.py decode_output.tab 
