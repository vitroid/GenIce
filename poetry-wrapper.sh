#!/bin/bash
# Poetry wrapper script that sets POETRY_PYTHON and creates python -> python3.13 symlink in PATH
# Poetry sometimes looks for 'python' command directly, so we need to make it available

# Determine the correct Python path
PYTHON_PATH=""
PYTHON_BIN_DIR=""

# Try Homebrew Python 3.13 first (if available and in range >=3.11,<3.14)
if [ -f "/opt/homebrew/bin/python3.13" ]; then
    PYTHON_PATH=/opt/homebrew/bin/python3.13
    PYTHON_BIN_DIR=/opt/homebrew/bin
# Try Homebrew Python 3.11 (original .venv Python)
elif [ -f "/opt/homebrew/Cellar/python@3.11/3.11.4_1/Frameworks/Python.framework/Versions/3.11/bin/python3.11" ]; then
    PYTHON_PATH=/opt/homebrew/Cellar/python@3.11/3.11.4_1/Frameworks/Python.framework/Versions/3.11/bin/python3.11
    PYTHON_BIN_DIR=/opt/homebrew/Cellar/python@3.11/3.11.4_1/Frameworks/Python.framework/Versions/3.11/bin
# Try .venv/bin/python if it exists
elif [ -f ".venv/bin/python" ]; then
    PYTHON_PATH="$(pwd)/.venv/bin/python"
    PYTHON_BIN_DIR="$(pwd)/.venv/bin"
elif [ -f ".venv/bin/python3" ]; then
    PYTHON_PATH="$(pwd)/.venv/bin/python3"
    PYTHON_BIN_DIR="$(pwd)/.venv/bin"
# Fallback: try to find python3.13 or python3.11 in PATH
elif command -v python3.13 >/dev/null 2>&1; then
    PYTHON_PATH=$(command -v python3.13)
    PYTHON_BIN_DIR=$(dirname "$PYTHON_PATH")
elif command -v python3.11 >/dev/null 2>&1; then
    PYTHON_PATH=$(command -v python3.11)
    PYTHON_BIN_DIR=$(dirname "$PYTHON_PATH")
else
    echo "Warning: Could not find Python 3.11 or 3.13. Using system python3 (may be too old)." >&2
    PYTHON_PATH=/usr/bin/python3
    PYTHON_BIN_DIR=/usr/bin
fi

# Set POETRY_PYTHON environment variable (Poetry does reference this)
export POETRY_PYTHON="$PYTHON_PATH"

# Create a temporary directory with a 'python' symlink that points to the correct Python
# This is needed because Poetry sometimes looks for 'python' command directly
TEMP_BIN_DIR=$(mktemp -d)
trap "rm -rf $TEMP_BIN_DIR" EXIT

# Create python -> python3.13 (or appropriate version) symlink
PYTHON_NAME=$(basename "$PYTHON_PATH")
ln -s "$PYTHON_PATH" "$TEMP_BIN_DIR/python"
ln -s "$PYTHON_PATH" "$TEMP_BIN_DIR/$PYTHON_NAME"

# Add temp directory to PATH (at the beginning so it takes precedence)
export PATH="$TEMP_BIN_DIR:$PYTHON_BIN_DIR:$PATH"

# Also try to set the env explicitly if poetry env use works
# But first check if we're in a poetry project
if [ -f "pyproject.toml" ]; then
    # Try to set the env (this may fail if Poetry has permission issues, but we try anyway)
    poetry env use "$PYTHON_PATH" 2>/dev/null || true
fi

# Execute poetry with all arguments
exec poetry "$@"

