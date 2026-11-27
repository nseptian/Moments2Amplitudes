#!/bin/bash
set -e

# Setup ROOT 6.36 environment
if [ -f "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh" ]; then
    source /d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh
elif [ -f "/d/grid17/septian/root_v6-36-08_cuda_build/bin/thisroot.sh" ]; then
    source /d/grid17/septian/root_v6-36-08_cuda_build/bin/thisroot.sh
fi

if [ -z "$BRUFIT" ]; then
    if [ -d "/d/grid17/septian/brufit_dev/brufit_negative" ]; then
        export BRUFIT="/d/grid17/septian/brufit_dev/brufit_negative"
    elif [ -d "/d/grid17/septian/brufit_dev/brufit_cov" ]; then
        export BRUFIT="/d/grid17/septian/brufit_dev/brufit_cov"
    fi
fi

dir="$1"
savePath="$2"
outName="$3"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
macroPath="${SCRIPT_DIR}/ExportCovarianceC.C"
if [ ! -f "$macroPath" ]; then
    if [ -f "${SCRIPT_DIR}/legacy/ExportCovarianceC.C" ]; then
        macroPath="${SCRIPT_DIR}/legacy/ExportCovarianceC.C"
    elif [ -f "${SCRIPT_DIR}/../macros/ExportCovarianceC.C" ]; then
        macroPath="${SCRIPT_DIR}/../macros/ExportCovarianceC.C"
    fi
fi

if [ -n "$outName" ]; then
    root -l -b -q "${macroPath}(\"${dir}\",\"${savePath}\",\"${outName}\")"
else
    root -l -b -q "${macroPath}(\"${dir}\")"
fi
