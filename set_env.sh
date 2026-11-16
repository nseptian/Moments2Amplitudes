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

# Set default BRUFIT and ROOTSYS if not already set
export BRUFIT="${BRUFIT:-/d/grid17/septian/brufit_dev/brufit}"
export ROOTSYS="${ROOTSYS:-/d/grid17/septian/root_v6-36-08_cuda_install}"

if [ -f "${ROOTSYS}/bin/thisroot.sh" ]; then
    source "${ROOTSYS}/bin/thisroot.sh"
    export LD_LIBRARY_PATH="${BRUFIT}/lib:${LD_LIBRARY_PATH}"
    export ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${MOMS2AMPS}/macros:${BRUFIT}/include:${ROOTSYS}/include:${ROOT_INCLUDE_PATH}"
elif [ -d "${ROOTSYS}" ]; then
    export PATH="${ROOTSYS}/bin:${PATH}"
    export LD_LIBRARY_PATH="${ROOTSYS}/lib:${BRUFIT}/lib:${LD_LIBRARY_PATH}"
    export ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${MOMS2AMPS}/macros:${BRUFIT}/include:${ROOTSYS}/include:${ROOT_INCLUDE_PATH}"
fi

# Pixi conda environment for Python runtime and packages (numpy, scipy, matplotlib, pyyaml, etc.)
if [ -d "${MOMS2AMPS}/.pixi/envs/default" ]; then
    export CONDA_PREFIX="${MOMS2AMPS}/.pixi/envs/default"
    export PATH="${PATH}:${CONDA_PREFIX}/bin"
    export PYTHONPATH="${MOMS2AMPS}/python:${CONDA_PREFIX}/lib/python3.9/site-packages:${PYTHONPATH}"
else
    export PYTHONPATH="${MOMS2AMPS}/python:${PYTHONPATH}"
fi

# Ensure ROOTSYS/bin takes precedence over any conda/pixi binaries for root and rootcling
if [ -d "${ROOTSYS}/bin" ]; then
    export PATH="${ROOTSYS}/bin:${PATH}"
fi

export MPLCONFIGDIR="/tmp"
