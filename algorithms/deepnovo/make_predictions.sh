#!/bin/bash
# Replace deepnovo_config.py with config file
cp config.py DeepNovo/deepnovo_config.py

# Convert input data to model format
python input_mapper.py "$@" ./input_data.mgf
# Run de novo algorithm on the input data
python DeepNovo/deepnovo_main.py --train_dir . --decode --beam_search --beam_size 5
# Convert predictions to the general output format
python output_mapper.py decode_output.tab 

# # Iterate through each file in the input directory
# for input_file in "$@"/*; do
#     if [ -f "$input_file" ]; then
#         echo "Processing file: $input_file"
#         # Convert input data to model format
#         python input_mapper.py "$input_file" ./input_data.mgf
#         # Run de novo algorithm on the input data
#         python DeepNovo/deepnovo_main.py --train_dir . --decode --beam_search --beam_size 5
#         # Collect predictions
#         cat decode_output.tab >> outputs.tab
#     fi
# done

# # Convert predictions to the general output format
# python output_mapper.py outputs.tab 