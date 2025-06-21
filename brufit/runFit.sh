#!/bin/sh

# This script is used to run the BruFit analysis on a ROOT file.
rootFile='TestMassDependentMoments.C'
loadFile='Load.C'

source /d/grid17/septian/brufit_dev/brufit/set_env.sh

for i in {0..9}
do
    echo "Running BruFit for file $i"
    brufit -l -b -q Load.C "$rootFile($i)" > "output_$i.txt"
    if [ $? -ne 0 ]; then
        echo "Error running BruFit for file $i"
        exit 1
    fi
done