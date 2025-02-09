#!/bin/bash
dset_dir="$1"
algorithm="$2"
spectra_dir="$dset_dir/mgf"
output_root_dir="./outputs"
time_log_root_dir="./times"
overlay_size=1024

# TODO maybe now we need separate output dir for each dataset
dset_name=$(basename "$dset_dir")
output_dir="$output_root_dir/$dset_name"
time_log_dir="$time_log_root_dir/$dset_name"

# Echo message based on whether an algorithm is provided
if [ -z "$algorithm" ]; then
    echo "Running benchmark with all algorithms on dataset $dset_name."
else
    echo "Running benchmark with $algorithm on dataset $dset_name."
fi

recalculate=false

while getopts ":r" opt; do
  case $opt in
    r) recalculate=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
echo "Recalculate all algorithm outputs: $recalculate"

if "$recalculate"; then
    # Clean output dir 
    rm -rf "$output_dir"
    rm -rf "$time_log_dir"
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
mkdir -p "$time_log_dir"

# List input files
echo "Processing dataset: $dset_name ($dset_dir)"
ls "$spectra_dir"/*.mgf

# Loop through each algorithm in the algorithms directory
for algorithm_dir in algorithms/*; do

    if [ -d "$algorithm_dir" ] && [ $(basename "$algorithm_dir") != "base" ]; then
        algorithm_name=$(basename "$algorithm_dir")

        # If an algorithm is specified, only continue if algorithm_name matches
        if [ -z "$algorithm" ] || [ "$algorithm_name" == "$algorithm" ]; then

            time_log_file="$time_log_dir/${algorithm_name}_time.log"
            output_file="$output_dir/${algorithm_name}_output.csv"
            echo "Output file: $output_file"
            
            # Check if the output file does not exist
            if [ ! -e "$output_file" ]; then
                echo "Processing algorithm: $algorithm_name"

                # Remove an existing container overlay, if any
                rm -rf "algorithms/${algorithm_name}/overlay.img"
                # Create writable overlay for the container
                apptainer overlay create --fakeroot --size $overlay_size --sparse "algorithms/${algorithm_name}/overlay.img"

                # Calculate predictions
                echo "RUN ALGORITHM $algorithm_name"
                { time ( apptainer exec --fakeroot --nv \
                    --overlay "algorithms/${algorithm_name}/overlay.img" \
                    -B "${spectra_dir}":"/algo/${dset_name}" \
                    --env-file .env \
                    "algorithms/${algorithm_name}/container.sif" \
                    bash -c "cd /algo && ./make_predictions.sh ${dset_name}" 2>&1 ); } 2> "$time_log_file"
                
                # Collect predictions in output_dir
                echo "EXPORT PREDICTIONS"
                apptainer exec --fakeroot \
                    --overlay "algorithms/${algorithm_name}/overlay.img" \
                    -B "${output_dir}":/algo/outputs \
                    --env-file .env \
                    "algorithms/${algorithm_name}/container.sif" \
                    bash -c "cp /algo/outputs.csv /algo/outputs/${algorithm_name}_output.csv"

            else
                echo "Skipping algorithm: $algorithm_name. Output file already exists."
            fi

        fi

    fi
done

# Evaluate predictions
# TODO: add results_dir explicit definition
echo "EVALUATE PREDICTIONS"
apptainer exec --fakeroot --env-file .env "evaluation.sif" \
    bash -c "python evaluate.py ${output_dir}/ ${dset_dir}"
# TODO change to dset_dir/labels? 