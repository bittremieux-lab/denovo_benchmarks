#!/bin/bash
spectra_dir="$1"
output_dir="./outputs"

# List input files
echo "Spectra data dir: $spectra_dir"
ls "$spectra_dir"/*.mgf

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each algorithm in the algorithms directory
for algorithm_dir in algorithms/*; do
    if [ -d "$algorithm_dir" ]; then
        algorithm_name=$(basename "$algorithm_dir")
        
        echo "$algorithm_dir"
        echo "Running $algorithm_name"
        
        # Build the Docker image for the current algorithm
        echo "BUILD ALGORITHM DOCKER"
        docker buildx build --platform linux/amd64 -f "$algorithm_dir/Dockerfile" -t "${algorithm_name}-docker" .
        
        # Run the algorithm Docker container
        echo "RUN ALGORITHM DOCKER"
        docker run --name "${algorithm_name}-container" "${algorithm_name}-docker" "$spectra_dir"
        
        # Export predictions from the Docker container
        echo "EXPORT PREDICTIONS"
        docker cp "${algorithm_name}-container:/app/outputs.csv" "$output_dir/${algorithm_name}_outputs.csv"

    fi
done

# Evaluate predictions
echo "EVALUATE PREDICTIONS"
python3 evaluate.py "$output_dir/" "$spectra_dir"

docker container prune -f
