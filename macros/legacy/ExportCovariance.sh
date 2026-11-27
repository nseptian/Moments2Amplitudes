#!/bin/bash
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

dir=$1
savePath=$2
outName=$3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
macroPath="${SCRIPT_DIR}/"
if [ ! -f "${macroPath}ExportCovarianceC.C" ]; then
    if [ -f "${SCRIPT_DIR}/legacy/ExportCovarianceC.C" ]; then
        macroPath="${SCRIPT_DIR}/legacy/"
    elif [ -f "${SCRIPT_DIR}/../legacy/ExportCovarianceC.C" ]; then
        macroPath="${SCRIPT_DIR}/../legacy/"
    fi
fi

root -l -b -q "${macroPath}ExportCovarianceC.C(\"${dir}\",\"${savePath}\",\"${outName}\")"
