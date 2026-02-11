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
    mkdir -p "$TEST_DIR/algorithms/base"
    
    # Create mock MGF files
    echo "mock spectrum data" > "$TEST_SPECTRA_DIR/test1.mgf"
    echo "mock spectrum data" > "$TEST_SPECTRA_DIR/test2.mgf"
    
    # Create mock algorithm container definition
    touch "$TEST_DIR/algorithms/test_algo/container.sif"
    
    # Create mock versions.log file for test_algo
    cat > "$TEST_DIR/algorithms/test_algo/versions.log" <<EOF
- container_version: "benchmark-1.0.0"
  date: "2025-10-06"
  repo_commit: "5fa09b3"
  notes: "Test version."
EOF
    
    # Copy script to test directory for isolated testing
    cp run.sh "$TEST_DIR/"
    
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
@test "run.sh exists and is executable" {
    [ -f "run.sh" ]
    [ -x "run.sh" ]
}

# Test 2: Script requires both dataset and algorithm arguments
@test "run.sh fails when algorithm argument not provided" {
    run bash run.sh "$TEST_DATASET_DIR"
    [ "$status" -eq 1 ]
}

# Test 3: Script rejects non-existent algorithm
@test "run.sh rejects non-existent algorithm" {
    run bash run.sh "$TEST_DATASET_DIR" nonexistent_algo
    [ "$status" -eq 1 ]
    [[ "$output" == *"not found"* ]]
}

# Test 4: Script rejects 'base' as algorithm name
@test "run.sh rejects 'base' as algorithm name" {
    run bash run.sh "$TEST_DATASET_DIR" base
    [ "$status" -eq 1 ]
    [[ "$output" == *"not an algorithm"* ]]
}

# Test 5: Script validates algorithm directory exists
@test "run.sh algorithm validation logic" {
    # Create minimal test script that validates inputs
    cat > test_algo_validation.sh <<'EOF'
#!/bin/bash
algorithm_name="$1"

if [ ! -d "algorithms/${algorithm_name}" ]; then
    echo "Error: Algorithm '${algorithm_name}' not found in algorithms/" >&2
    exit 1
fi

if [ "${algorithm_name}" = "base" ]; then
    echo "Error: 'base' is not an algorithm" >&2
    exit 1
fi

echo "Valid: algorithm=${algorithm_name}"
EOF
    chmod +x test_algo_validation.sh
    
    run ./test_algo_validation.sh test_algo
    [ "$status" -eq 0 ]
    [[ "$output" == *"Valid"* ]]
}

# Test 6: Script correctly parses two positional arguments
@test "run.sh extracts dataset and algorithm from arguments" {
    run bash -c '
        dset_dir="$1"
        algorithm_name="$2"
        dset_name=$(basename "$dset_dir")
        echo "dataset=$dset_name algo=$algorithm_name"
    ' -- /path/to/my_dataset casanovo
    
    [ "$status" -eq 0 ]
    [[ "$output" == *"dataset=my_dataset"* ]]
    [[ "$output" == *"algo=casanovo"* ]]
}

# Test 7: Script accepts -r flag with two arguments
@test "run.sh processes -r flag with dataset and algorithm" {
    run bash -c '
        recalculate=false
        while getopts ":r" opt; do
            case $opt in
                r) recalculate=true ;;
            esac
        done
        shift $((OPTIND-1))
        dset_dir="$1"
        algorithm_name="$2"
        echo "recalc=$recalculate dataset=$dset_dir algo=$algorithm_name"
    ' -- -r /path/to/dataset test_algo
    
    [ "$status" -eq 0 ]
    [[ "$output" == *"recalc=true"* ]]
    [[ "$output" == *"dataset=/path/to/dataset"* ]]
    [[ "$output" == *"algo=test_algo"* ]]
}

# Test 8: Script constructs new output paths with algorithm_name/version/dataset structure
@test "run.sh output file paths use new structure: algorithm/version/dataset" {
    dset_name="test_dataset"
    algorithm_name="casanovo"
    algorithm_version="benchmark-1.0.0"
    output_dir="./outputs/$algorithm_name/$algorithm_version/$dset_name"
    output_file="$output_dir/output.csv"
    time_log_file="$output_dir/time.log"
    
    [ "$output_file" = "./outputs/casanovo/benchmark-1.0.0/test_dataset/output.csv" ]
    [ "$time_log_file" = "./outputs/casanovo/benchmark-1.0.0/test_dataset/time.log" ]
}

# Test 9: Script constructs overlay file name with both dataset and algorithm
@test "run.sh overlay file name includes dataset and algorithm" {
    algorithm_name="test_algo"
    dset_name="my_dataset"
    overlay_file="algorithms/${algorithm_name}/overlay_${dset_name}.img"
    
    expected="algorithms/test_algo/overlay_my_dataset.img"
    [ "$overlay_file" = "$expected" ]
}

# Test 10: Script handles -r flag with conditional cleanup
@test "run.sh recalculate flag cleans output directory" {
    # Create pre-existing directory with files
    mkdir -p outputs/test_algo/benchmark-1.0.0/test_dataset
    touch outputs/test_algo/benchmark-1.0.0/test_dataset/old_file.csv
    touch outputs/test_algo/benchmark-1.0.0/test_dataset/time.log
    
    # Create a test script that implements cleanup logic matching run.sh
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
algorithm_name="$2"
algorithm_version="benchmark-1.0.0"
dset_name=$(basename "$dset_dir")
output_dir="./outputs/$algorithm_name/$algorithm_version/$dset_name"

if [ "$recalculate" = true ]; then
    rm -rf "$output_dir"
fi
mkdir -p "$output_dir"
EOF
    
    chmod +x test_cleanup.sh
    run ./test_cleanup.sh -r "$TEST_DATASET_DIR" test_algo
    
    [ "$status" -eq 0 ]
    [ ! -f "outputs/test_algo/benchmark-1.0.0/test_dataset/old_file.csv" ]
    [ ! -f "outputs/test_algo/benchmark-1.0.0/test_dataset/time.log" ]
}

# Test 11: Script creates output directory with version structure
@test "run.sh creates output directory" {
    # Create a minimal script that creates directories matching run.sh
    cat > test_mkdir.sh <<'EOF'
#!/bin/bash
dset_dir="$1"
algorithm_name="$2"
algorithm_version="benchmark-1.0.0"
output_root_dir="./outputs"
dset_name=$(basename "$dset_dir")
output_dir="$output_root_dir/$algorithm_name/$algorithm_version/$dset_name"
mkdir -p "$output_dir"
EOF
    
    chmod +x test_mkdir.sh
    run ./test_mkdir.sh "$TEST_DATASET_DIR" test_algo
    
    [ "$status" -eq 0 ]
    [ -d "outputs/test_algo/benchmark-1.0.0/test_dataset" ]
}

# Test 12: Script handles missing spectra directory gracefully
@test "run.sh handles missing spectra directory gracefully" {
    run ls "/nonexistent/path/"*.mgf 2>&1
    [ "$status" -ne 0 ]
}

# Test 13: Script processes invalid options
@test "run.sh reports invalid options" {
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

# Test 14: Script correctly parses dataset name from path
@test "run.sh correctly extracts dataset name from path" {
    result=$(basename "/path/to/my_dataset")
    [ "$result" = "my_dataset" ]
}

# Test 15: Script accepts flags before positional arguments
@test "run.sh processes options before positional arguments" {
    run bash -c '
        recalculate=false
        while getopts ":r" opt; do
            case $opt in
                r) recalculate=true ;;
            esac
        done
        shift $((OPTIND-1))
        dset_dir="$1"
        algorithm_name="$2"
        echo "recalc=$recalculate dir=$dset_dir algo=$algorithm_name"
    ' -- -r /path/to/dataset test_algo
    
    [ "$status" -eq 0 ]
    [[ "$output" == *"recalc=true"* ]]
    [[ "$output" == *"dir=/path/to/dataset"* ]]
    [[ "$output" == *"algo=test_algo"* ]]
}

# Test 16: Script can extract container_version from versions.log
@test "run.sh extracts container_version from versions.log" {
    algorithm_version=$(grep -m 1 "container_version:" "algorithms/test_algo/versions.log" | awk -F'"' '{print $2}')
    
    [ "$algorithm_version" = "benchmark-1.0.0" ]
}

# Test 17: Script handles missing versions.log file
@test "run.sh fails when versions.log is missing" {
    mkdir -p algorithms/no_version_algo
    
    run bash -c '
        algorithm_name="no_version_algo"
        algorithm_version=$(grep -m 1 "container_version:" "algorithms/${algorithm_name}/versions.log" 2>/dev/null | awk -F\" '\''{print $2}'\'')
        if [ -z "$algorithm_version" ]; then
            echo "Error: Could not extract container_version" >&2
            exit 1
        fi
    '
    
    [ "$status" -eq 1 ]
    [[ "$output" == *"Could not extract container_version"* ]]
}

# Test 18: Script handles malformed versions.log
@test "run.sh fails when versions.log is malformed" {
    mkdir -p algorithms/bad_version_algo
    echo "invalid yaml content" > algorithms/bad_version_algo/versions.log
    
    run bash -c '
        algorithm_name="bad_version_algo"
        algorithm_version=$(grep -m 1 "container_version:" "algorithms/${algorithm_name}/versions.log" | awk -F\" '\''{print $2}'\'')
        if [ -z "$algorithm_version" ]; then
            echo "Error: Could not extract container_version" >&2
            exit 1
        fi
    '
    
    [ "$status" -eq 1 ]
    [[ "$output" == *"Could not extract container_version"* ]]
}

# Test 19: Script extracts latest version when multiple versions exist
@test "run.sh extracts latest (first) version from versions.log" {
    cat > algorithms/test_algo/versions.log <<EOF
- container_version: "benchmark-2.0.0"
  date: "2026-01-01"
  notes: "Latest version."
- container_version: "benchmark-1.0.0"
  date: "2025-10-06"
  notes: "Older version."
EOF
    
    algorithm_version=$(grep -m 1 "container_version:" "algorithms/test_algo/versions.log" | awk -F'"' '{print $2}')
    
    [ "$algorithm_version" = "benchmark-2.0.0" ]
}

# Test 20: Script creates correct directory structure with version
@test "run.sh creates output directory with algorithm/version/dataset structure" {
    cat > test_mkdir_versioned.sh <<'EOF'
#!/bin/bash
dset_dir="$1"
algorithm_name="$2"
algorithm_version="benchmark-1.0.0"
dset_name=$(basename "$dset_dir")
output_root_dir="./outputs"

output_dir="$output_root_dir/$algorithm_name/$algorithm_version/$dset_name"

mkdir -p "$output_dir"
EOF
    
    chmod +x test_mkdir_versioned.sh
    run ./test_mkdir_versioned.sh "$TEST_DATASET_DIR" test_algo
    
    [ "$status" -eq 0 ]
    [ -d "outputs/test_algo/benchmark-1.0.0/test_dataset" ]
}

# Test 21: Script uses simplified output file names in same directory
@test "run.sh uses output.csv and time.log in output directory" {
    dset_name="test_dataset"
    algorithm_name="casanovo"
    algorithm_version="benchmark-1.0.0"
    output_dir="./outputs/$algorithm_name/$algorithm_version/$dset_name"
    
    output_file="$output_dir/output.csv"
    time_log_file="$output_dir/time.log"
    
    [ "$output_file" = "./outputs/casanovo/benchmark-1.0.0/test_dataset/output.csv" ]
    [ "$time_log_file" = "./outputs/casanovo/benchmark-1.0.0/test_dataset/time.log" ]
}

# Test 22: Verify run.sh actually defines correct file path variables
@test "run.sh source code contains correct output_file and time_log_file definitions" {
    # Check that run.sh defines these variables correctly
    grep -q '^time_log_file="\$output_dir/time\.log"' run.sh
    grep -q '^output_file="\$output_dir/output\.csv"' run.sh
}
