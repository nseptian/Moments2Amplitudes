#!/bin/bash

# This script is used to run the BruFit analysis on a ROOT file with parallel processing.
rootFile='TestMassDependentMoments.C'
loadFile='Load.C'

source /d/grid17/septian/brufit_dev/brufit/set_env.sh

# Function to run a single fit
run_fit() {
    local i=$1
    echo "Starting BruFit for file $i"
    root -l -b -q $BRUFIT/macros/LoadBru.C Load.C "$rootFile($i)" > "output_$i.txt" 2>&1
    if [ $? -eq 0 ]; then
        echo "Completed BruFit for file $i"
    else
        echo "Error running BruFit for file $i" >&2
    fi
}

# Maximum number of simultaneous jobs
MAX_JOBS=100

# Run fits in parallel
for i in {1..10000}
do
    # Wait if we've reached the maximum number of background jobs
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 0.1
    done
    
    # Start the fit in background (call function directly, not exported)
    run_fit $i &
    
    echo "Started job $i ($(jobs -r | wc -l) jobs running)"
done

# Wait for all background jobs to complete
wait

echo "All fits completed!"