#!/bin/bash
echo "Your container args are: $@"
# Run de novo algorithm on the input data
casanovo sequence -c config.yml -o outputs.mztab "$@"/*.mgf
# Convert predictions to the general output format
python3 ./output_mapper.py outputs.mztab
