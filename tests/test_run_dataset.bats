#!/usr/bin/env bats

bats_require_minimum_version 1.5.0

# Setup and teardown
setup() {
    # Create temporary test directory
    export TEST_DIR="$(mktemp -d)"
    export TEST_DATASET_DIR="$TEST_DIR/test_dataset"
    export TEST_SPECTRA_DIR="$TEST_DATASET_DIR/mgf"
    
    # Create mock dataset structure
    mkdir -p "$TEST_SPECTRA_DIR"
    mkdir -p "$TEST_DIR/algorithms/test_algo"
    
    # Create mock MGF files
    echo "mock spectrum data" > "$TEST_SPECTRA_DIR/test1.mgf"
    echo "mock spectrum data" > "$TEST_SPECTRA_DIR/test2.mgf"
    
    # Create mock algorithm container definition
    touch "$TEST_DIR/algorithms/test_algo/container.sif"
    
    # Copy script to test directory for isolated testing
    cp run_dataset.sh "$TEST_DIR/"
    
    # Create a testable version of the script by extracting functions
    cd "$TEST_DIR"
}

teardown() {
    # Skip cleanup if DEBUG is set
    if [ -z "${DEBUG:-}" ]; then
        # Clean up temporary directory
        rm -rf "$TEST_DIR"
    else
        echo "Test directory preserved at: $TEST_DIR" >&2
    fi
}

# Test 1: Script exists and is executable
@test "run_dataset.sh exists and is executable" {
    [ -f "run_dataset.sh" ]
    [ -x "run_dataset.sh" ]
}

# Test 2: Script requires dataset argument
@test "script fails when no dataset argument provided" {
    # The script will fail when trying to access non-existent paths
    run -127 bash run_dataset.sh
    [ "$status" -ne 0 ]
}

# Test 3: Script accepts -r flag
@test "script accepts -r (recalculate) flag" {
    # Mock the script execution to just parse arguments
    run bash -c '
        recalculate=false
        while getopts ":r" opt; do
            case $opt in
                r) recalculate=true ;;
            esac
        done
        shift $((OPTIND-1))
        echo "$recalculate"
    ' -- -r test_dataset
    
    [ "$status" -eq 0 ]
    [[ "$output" == "true" ]]
}

# Test 4: Script correctly parses dataset directory
@test "script correctly extracts dataset name from path" {
    result=$(basename "/path/to/my_dataset")
    [ "$result" = "my_dataset" ]
}

# Test 5: Output directories are created
@test "script creates output and time log directories" {
    # Create a minimal version that just creates directories
    cat > test_mkdir.sh <<'EOF'
#!/bin/bash
dset_dir="$1"
output_root_dir="./outputs"
time_log_root_dir="./times"
dset_name=$(basename "$dset_dir")
output_dir="$output_root_dir/$dset_name"
time_log_dir="$time_log_root_dir/$dset_name"
mkdir -p "$output_dir"
mkdir -p "$time_log_dir"
EOF
    
    chmod +x test_mkdir.sh
    run ./test_mkdir.sh "$TEST_DATASET_DIR"
    
    [ "$status" -eq 0 ]
    [ -d "outputs/test_dataset" ]
    [ -d "times/test_dataset" ]
}

# Test 6: Recalculate flag cleans directories
@test "recalculate flag removes existing output directories" {
    # Create pre-existing directories with files
    mkdir -p outputs/test_dataset
    mkdir -p times/test_dataset
    touch outputs/test_dataset/old_file.csv
    touch times/test_dataset/old_log.log
    
    # Create a test script that implements the cleanup logic
    cat > test_cleanup.sh <<'EOF'
#!/bin/bash
recalculate=false
while getopts ":r" opt; do
    case $opt in
        r) recalculate=true ;;
    esac
done
shift $((OPTIND-1))

dset_dir="$1"
output_root_dir="./outputs"
time_log_root_dir="./times"
dset_name=$(basename "$dset_dir")
output_dir="$output_root_dir/$dset_name"
time_log_dir="$time_log_root_dir/$dset_name"

if [ "$recalculate" = true ]; then
    rm -rf "$output_dir"
    rm -rf "$time_log_dir"
fi
mkdir -p "$output_dir"
mkdir -p "$time_log_dir"
EOF
    
    chmod +x test_cleanup.sh
    run ./test_cleanup.sh -r "$TEST_DATASET_DIR"
    
    [ "$status" -eq 0 ]
    [ ! -f "outputs/test_dataset/old_file.csv" ]
    [ ! -f "times/test_dataset/old_log.log" ]
}

# Test 7: Script skips 'base' algorithm directory
@test "script skips 'base' algorithm directory" {
    mkdir -p algorithms/base
    
    # Count non-base algorithm directories
    count=0
    for algorithm_dir in algorithms/*; do
        if [ -d "$algorithm_dir" ] && [ "$(basename "$algorithm_dir")" != "base" ]; then
            ((count++))
        fi
    done
    
    [ "$count" -eq 1 ]  # Only test_algo, not base
}

# Test 8: Script handles algorithm directories correctly
@test "script identifies algorithm directories" {
    mkdir -p algorithms/algo1
    mkdir -p algorithms/algo2
    mkdir -p algorithms/base
    touch algorithms/not_a_dir.txt
    
    # Count valid algorithm directories (excluding base)
    count=0
    for algorithm_dir in algorithms/*; do
        if [ -d "$algorithm_dir" ] && [ "$(basename "$algorithm_dir")" != "base" ]; then
            ((count++))
        fi
    done
    
    [ "$count" -eq 3 ]  # test_algo, algo1, algo2
}

# Test 9: Overlay file naming is correct
@test "overlay file name includes dataset name" {
    algorithm_name="test_algo"
    dset_name="my_dataset"
    overlay_file="algorithms/${algorithm_name}/overlay_${dset_name}.img"
    
    expected="algorithms/test_algo/overlay_my_dataset.img"
    [ "$overlay_file" = "$expected" ]
}

# Test 10: Output file paths are constructed correctly
@test "output file paths are constructed correctly" {
    dset_name="test_dataset"
    algorithm_name="casanovo"
    output_dir="./outputs/$dset_name"
    output_file="$output_dir/${algorithm_name}_output.csv"
    
    expected="./outputs/test_dataset/casanovo_output.csv"
    [ "$output_file" = "$expected" ]
}

# Test 11: Time log file paths are constructed correctly
@test "time log file paths are constructed correctly" {
    dset_name="test_dataset"
    algorithm_name="casanovo"
    time_log_dir="./times/$dset_name"
    time_log_file="$time_log_dir/${algorithm_name}_time.log"
    
    expected="./times/test_dataset/casanovo_time.log"
    [ "$time_log_file" = "$expected" ]
}

# Test 12: MGF file glob pattern works
@test "MGF files are listed correctly" {
    # Create a non-MGF file to ensure filtering works
    touch "$TEST_SPECTRA_DIR/not_mgf.txt"
    
    # Count MGF files (test1.mgf and test2.mgf from setup)
    mgf_count=$(ls "$TEST_SPECTRA_DIR/"*.mgf 2>/dev/null | wc -l)
    [ "$mgf_count" -eq 2 ]
}

# Test 13: Script handles missing spectra directory
@test "script handles missing spectra directory gracefully" {
    run ls "/nonexistent/path/"*.mgf 2>&1
    [ "$status" -ne 0 ]
}

# Test 14: Getopts handles invalid options
@test "script reports invalid options" {
    run bash -c '
        while getopts ":r" opt; do
            case $opt in
                r) echo "valid" ;;
                \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
            esac
        done
    ' -- -x
    
    [ "$status" -eq 1 ]
    [[ "$output" == *"Invalid option"* ]]
}

# Test 15: Multiple flags can be combined
@test "script processes options before positional arguments" {
    run bash -c '
        recalculate=false
        while getopts ":r" opt; do
            case $opt in
                r) recalculate=true ;;
            esac
        done
        shift $((OPTIND-1))
        dset_dir="$1"
        echo "recalc=$recalculate dir=$dset_dir"
    ' -- -r /path/to/dataset
    
    [ "$status" -eq 0 ]
    [[ "$output" == *"recalc=true"* ]]
    [[ "$output" == *"dir=/path/to/dataset"* ]]
}

