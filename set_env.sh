#!/bin/bash
# Setup environment for Moments2Amplitudes

# Determine repository root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MOMS2AMPS="${SCRIPT_DIR}"

# Source .env if present
if [ -f "${MOMS2AMPS}/.env" ]; then
    set -a
    source "${MOMS2AMPS}/.env"
    set +a
fi

# Set default BRUFIT if not already set in environment or .env
export BRUFIT="${BRUFIT:-/d/grid17/septian/brufit_dev/brufit}"

# Configure ROOT and Python include paths
export ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${BRUFIT}/include:${ROOT_INCLUDE_PATH}"
export PYTHONPATH="${MOMS2AMPS}/python:${PYTHONPATH}"
export MPLCONFIGDIR="/tmp"
