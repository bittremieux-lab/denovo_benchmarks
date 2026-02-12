#!/usr/bin/env bats

# Integration tests for run_dataset.sh
# These tests verify the script works end-to-end with mocked external dependencies

load 'test_helper/bats-support/load'
load 'test_helper/bats-assert/load'

setup() {
    export TEST_DIR="$(mktemp -d)"
    export ORIG_DIR="$PWD"
    
    # Create mock project structure
    mkdir -p "$TEST_DIR/sample_data/test_dataset/mgf"
    mkdir -p "$TEST_DIR/algorithms/mock_algo"
    mkdir -p "$TEST_DIR/algorithms/base"
    
    # Create mock MGF files
    cat > "$TEST_DIR/sample_data/test_dataset/mgf/spectrum1.mgf" <<EOF
BEGIN IONS
TITLE=Spectrum 1
PEPMASS=500.25
CHARGE=2+
100.0 1000.0
200.0 2000.0
END IONS
EOF
    
    # Create versions.log for mock algorithm
    cat > "$TEST_DIR/algorithms/mock_algo/versions.log" <<EOF
- container_version: "mock-1.0.0"
  notes: "Mock algorithm for testing"
EOF
    
    # Create mock make_predictions.sh inside algorithm
    cat > "$TEST_DIR/algorithms/mock_algo/make_predictions.sh" <<'SCRIPT'
#!/bin/bash
mgf_dir="$1"
output_csv="$2"
time_log="$3"

echo "scan,sequence,score" > "$output_csv"
echo "1,PEPTIDE,0.95" >> "$output_csv"
echo "real 0.100s" > "$time_log"
SCRIPT
    chmod +x "$TEST_DIR/algorithms/mock_algo/make_predictions.sh"
    
    # Create mock container (dummy file)
    touch "$TEST_DIR/algorithms/mock_algo/container.sif"
    
    # Create mock .env file
    touch "$TEST_DIR/.env"
    
    # Copy scripts to test environment
    cp "$ORIG_DIR/run_dataset.sh" "$TEST_DIR/"
    cp -r "$ORIG_DIR/mock_env" "$TEST_DIR/"
    
    cd "$TEST_DIR"
    
    # Load mock environment
    source mock_env/setup_mock.sh
}

teardown() {
    cd "$ORIG_DIR"
    
    # Restore original environment if in test directory
    if [ -f "$TEST_DIR/mock_env/teardown_mock.sh" ]; then
        source "$TEST_DIR/mock_env/teardown_mock.sh"
    fi
    
    rm -rf "$TEST_DIR"
}

@test "integration: script runs successfully with valid dataset" {
    skip "Requires full apptainer environment"
    
    run bash run_dataset.sh sample_data/test_dataset
    
    [ "$status" -eq 0 ]
    [ -d "outputs/test_dataset" ]
    [ -d "times/test_dataset" ]
}

@test "integration: script creates expected output files" {
    skip "Requires full apptainer environment"
    
    bash run_dataset.sh sample_data/test_dataset
    
    [ -f "outputs/test_dataset/mock_algo_output.csv" ]
    [ -f "times/test_dataset/mock_algo_time.log" ]
}

@test "integration: script skips existing outputs without -r flag" {
    skip "Requires full apptainer environment"
    
    # First run
    bash run_dataset.sh sample_data/test_dataset
    
    # Create marker file
    echo "marker" > outputs/test_dataset/mock_algo_output.csv
    
    # Second run without -r
    bash run_dataset.sh sample_data/test_dataset
    
    # Marker file should still exist
    grep -q "marker" outputs/test_dataset/mock_algo_output.csv
}

@test "integration: script recalculates with -r flag" {
    skip "Requires full apptainer environment"
    
    # First run
    bash run_dataset.sh sample_data/test_dataset
    echo "marker" > outputs/test_dataset/mock_algo_output.csv
    
    # Second run with -r
    bash run_dataset.sh -r sample_data/test_dataset
    
    # Marker should be gone
    ! grep -q "marker" outputs/test_dataset/mock_algo_output.csv || [ ! -f outputs/test_dataset/mock_algo_output.csv ]
}
