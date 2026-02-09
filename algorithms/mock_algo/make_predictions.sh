#!/bin/bash
# Mock prediction script that mimics the container behavior
dataset_name="$1"

echo "Mock algorithm running for dataset: $dataset_name"

# Create outputs.csv in the mock_algo directory
# (Real container will write in internal /algo directory)
echo "spectrum_id,sequence,score,aa_scores" > outputs.csv
echo 'test_1:0,PEPTIDE,0.95,"0.95,0.95,0.95,0.95,0.95,0.95,0.95"' >> outputs.csv
echo 'test_1:1,PEPTIDE,0.87,"0.87,0.87,0.87,0.87,0.87,0.87,0.87"' >> outputs.csv
echo 'test_1:2,PTEPIDE,0.87,"0.87,0.87,0.87,0.87,0.87,0.87,0.87"' >> outputs.csv

realpath outputs.csv

exit 0