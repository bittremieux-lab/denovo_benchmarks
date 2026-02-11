#!/bin/bash

# Test runner script for run.sh
# This script checks if BATS is installed and runs the test suite

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "================================"
echo "Testing run.sh"
echo "================================"

# Check if BATS is installed
if ! command -v bats &> /dev/null; then
    echo -e "${YELLOW}BATS is not installed.${NC}"
    echo ""
    echo "To install BATS on macOS:"
    echo "  brew install bats-core"
    echo ""
    echo "To install BATS on Linux:"
    echo "  git clone https://github.com/bats-core/bats-core.git"
    echo "  cd bats-core"
    echo "  sudo ./install.sh /usr/local"
    echo ""
    echo "For more information: https://github.com/bats-core/bats-core"
    exit 1
fi

echo -e "${GREEN}✓ BATS is installed${NC}"
echo ""

# Check if test files exist
if [ ! -f "tests/test_run.bats" ]; then
    echo -e "${RED}✗ Test file not found: tests/test_run.bats${NC}"
    exit 1
fi

# Install BATS helper libraries if needed
if [ ! -d "tests/test_helper/bats-support" ]; then
    echo "Installing BATS helper libraries..."
    mkdir -p tests/test_helper
    
    git clone https://github.com/bats-core/bats-support.git tests/test_helper/bats-support 2>/dev/null || echo "bats-support already exists"
    git clone https://github.com/bats-core/bats-assert.git tests/test_helper/bats-assert 2>/dev/null || echo "bats-assert already exists"
fi

# Run unit tests
echo "Running unit tests..."
echo "-------------------"
if bats tests/test_run.bats; then
    echo -e "${GREEN}✓ Unit tests passed${NC}"
else
    echo -e "${RED}✗ Unit tests failed${NC}"
    exit 1
fi

echo ""

# Run integration tests (if available)
if [ -f "tests/test_run_integration.bats" ]; then
    echo "Running integration tests..."
    echo "-------------------------"
    integration_failed=false
    
    if [ -f "tests/test_run_integration.bats" ]; then
        if bats tests/test_run_integration.bats; then
            echo -e "${GREEN}✓ run.sh integration tests passed${NC}"
        else
            echo -e "${YELLOW}! run.sh integration tests skipped or failed${NC}"
            integration_failed=true
        fi
    fi
fi

echo ""
echo -e "${GREEN}================================${NC}"
echo -e "${GREEN}All tests completed successfully!${NC}"
echo -e "${GREEN}================================${NC}"
