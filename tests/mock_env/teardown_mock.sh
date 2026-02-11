#!/bin/bash
# Restore original environment after mock setup
# Undoes changes made by setup_mock.sh
# Usage: source tests/mock_env/teardown_mock.sh

if [ -z "$BASH_SOURCE" ]; then
    echo "Error: This script must be sourced, not executed" >&2
    return 1
fi

if [ -z "$ORIGINAL_PATH" ]; then
    echo "Error: ORIGINAL_PATH not set. Did you source tests/mock_env/setup_mock.sh?" >&2
    return 1
fi

# Restore original PATH
export PATH="$ORIGINAL_PATH"

# Clean up mock environment variables
unset MOCK_ENV_DIR
unset ORIGINAL_PATH

echo "✓ Original environment restored"