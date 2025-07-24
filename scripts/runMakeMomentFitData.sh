#!/usr/bin/bash

# This script runs the makeMomentFitData.C in root

# Usage: runMakeMomentFitData.sh <inputFile> <outputFile> <treeName>
# Example: runMakeMomentFitData.sh input.root output.root treeName

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <inputFile> <outputFile>"
    exit 1
fi

makeMomentFitDataCPath='/d/grid17/septian/Moments2Amplitudes/brufit/MakeMomentFitData.C'

filePath='/d/grid17/septian/Moments2Amplitudes/samples/'
inputFile=${filePath}$1
outputFile=${filePath}$2

# Check if the input file exists
if [ ! -f "$inputFile" ]; then
    echo "Input file $inputFile does not exist."
    exit 1
fi

# Load ROOT environment
source /d/grid17/septian/brufit_dev/brufit/set_env.sh

# Run the ROOT command
root -l -b -q "$makeMomentFitDataCPath(\"$inputFile\", \"$outputFile\")"
if [ $? -ne 0 ]; then
    echo "Error running ROOT command."
    exit 1
fi

# Check if the output file was created
if [ ! -f "$outputFile" ]; then
    echo "Output file $outputFile was not created."
    exit 1
fi
echo "Output file $outputFile created successfully."

echo "Adding ID branch to $outputFile".
root -l -b -q '$BRUFIT/macros/AddIDBranch.C("kin","'$outputFile'")'