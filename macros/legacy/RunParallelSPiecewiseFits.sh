#!/bin/bash
# ============================================================================
# RunParallelSPiecewiseFits.sh
# Runs randomized S-piecewise fits in parallel using available CPU cores
# ============================================================================

set -e

N_SEEDS=${1:-100}
MAX_JOBS=${2:-40}
D_WAVE_FILE=${3:-"result_tree_DBreitWigner_best.root"}
OUT_DIR=${4:-"logs_spiecewise"}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

mkdir -p "$OUT_DIR"

echo "========================================================================"
echo " Launching $N_SEEDS S-Piecewise Randomized Fits in Parallel"
echo " Max Parallel Workers: $MAX_JOBS"
echo " D-Wave Initial File:  $D_WAVE_FILE"
echo " Output Directory:     $OUT_DIR"
echo " Project Directory:    $PROJECT_DIR"
echo "========================================================================"

ROOT_CMD="root"
if command -v brufit >/dev/null 2>&1; then
    ROOT_CMD="brufit"
fi

# Pre-compile ACLiC libraries once upfront to avoid compilation race conditions
echo "Ensuring ACLiC libraries are loaded using $ROOT_CMD..."
$ROOT_CMD -l -b -q "${PROJECT_DIR}/macros/Load.C" -e 'return;' >/dev/null 2>&1 || true

run_single_seed() {
    local seed=$1
    local dwave_file=$2
    local out_dir=$3
    local proj_dir=$4
    local cmd_bin=$5
    local log_file="$out_dir/log_seed_${seed}.txt"

    $cmd_bin -l -b -q \
        "${proj_dir}/macros/Load.C" \
        "${proj_dir}/macros/RunSingleSeedSPiecewise.C(${seed}, \"${dwave_file}\", \"${out_dir}\")" \
        > "$log_file" 2>&1
    local status=$?
    if [[ $status -eq 0 ]]; then
        echo -n "."
    else
        echo -n "x"
    fi
}

export -f run_single_seed

echo "Running fits:"
for seed in $(seq 1 $N_SEEDS); do
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 0.05
    done
    run_single_seed "$seed" "$D_WAVE_FILE" "$OUT_DIR" "$PROJECT_DIR" "$ROOT_CMD" &
done

wait
echo ""
echo "========================================================================"
echo " All $N_SEEDS parallel S-piecewise fits completed."
echo "========================================================================"

