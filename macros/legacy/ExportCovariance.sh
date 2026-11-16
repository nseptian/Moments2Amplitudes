source /d/grid17/septian/root_v6-36-08_cuda_build/bin/thisroot.sh

export BRUFIT="/d/grid17/septian/brufit_dev/brufit_cov/"

dir=$1
savePath=$2
outName=$3

macroPath=/d/grid17/septian/Moments2Amplitudes/macros/

root -l -b -q $BRUFIT/macros/LoadBru.C ${macroPath}'ExportCovarianceC.C("'${dir}'","'${savePath}'","'${outName}'")'