#!/bin/bash
# Run de novo algorithm on the input data
casanovo sequence -c config.yml -o outputs.mztab "$@"/*.mgf
# Convert predictions to the general output format
python ./output_mapper.py --output_path=outputs.mztab
