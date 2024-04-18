#!/bin/bash
# Replace data_utils.py with config file
cp config.py DeepNovo/deepnovo_config.py

# Convert input data to model format
python input_mapper.py "$@" ./input_data.mgf
# Run de novo algorithm on the input data
python DeepNovo/deepnovo_main.py --train_dir . --decode --beam_search --beam_size 5
# Convert predictions to the general output format
python output_mapper.py decode_output.tab 
