#!/bin/bash
# ============================================================================
# RunParallelMCMCFits.sh
# Runs MCMC uncertainty propagation fits in parallel across CPU cores
# ============================================================================

set -e

N_WORKERS=${1:-40}
START_ENTRY=${2:-0}
N_ENTRIES=${3:-200}
STRIDE=${4:-1}
CENTRAL_FILE=${5:-"logs/central_best_params.root"}
OUT_DIR=${6:-"logs/mcmc_fits"}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

mkdir -p "$OUT_DIR"

echo "========================================================================"
echo " Launching Parallel MCMC Uncertainty Propagation"
echo " Parallel Workers:    $N_WORKERS"
echo " Total MCMC Entries:  $N_ENTRIES (start: $START_ENTRY, stride: $STRIDE)"
echo " Central Param File:  $CENTRAL_FILE"
echo " Output Directory:    $OUT_DIR"
echo " Project Directory:   $PROJECT_DIR"
echo "========================================================================"

if [ ! -f "$CENTRAL_FILE" ]; then
    echo "Error: Central best-fit parameter file not found: $CENTRAL_FILE"
    exit 1
fi

ROOT_CMD="root"
if command -v brufit >/dev/null 2>&1; then
    ROOT_CMD="brufit"
fi

# Pre-compile ACLiC libraries once upfront to avoid compilation race conditions
echo "Ensuring ACLiC libraries are loaded using $ROOT_CMD..."
$ROOT_CMD -l -b -q "${PROJECT_DIR}/macros/Load.C" -e 'return;' >/dev/null 2>&1 || true

# Calculate entries per worker chunk
ENTRIES_PER_WORKER=$(( (N_ENTRIES + N_WORKERS - 1) / N_WORKERS ))

run_single_worker() {
    local worker_id=$1
    local start=$2
    local count=$3
    local stride=$4
    local central_file=$5
    local out_dir=$6
    local proj_dir=$7
    local cmd_bin=$8
    local log_file="$out_dir/log_worker_${worker_id}.txt"

    $cmd_bin -l -b -q \
        "${proj_dir}/macros/Load.C" \
        "${proj_dir}/macros/RunSingleMCMCWorker.C(${worker_id}, ${start}, ${count}, ${stride}, \"${central_file}\", \"${out_dir}\")" \
        > "$log_file" 2>&1
    local status=$?
    if [[ $status -eq 0 ]]; then
        echo -n "."
    else
        echo -n "x"
    fi
}

export -f run_single_worker

echo "Running MCMC worker jobs:"
ACTIVE_WORKERS=0
for w in $(seq 0 $((N_WORKERS - 1))); do
    OFFSET=$(( w * ENTRIES_PER_WORKER ))
    if [ $OFFSET -ge $N_ENTRIES ]; then
        break
    fi
    COUNT=$(( N_ENTRIES - OFFSET ))
    if [ $COUNT -gt $ENTRIES_PER_WORKER ]; then
        COUNT=$ENTRIES_PER_WORKER
    fi
    CURR_START=$(( START_ENTRY + OFFSET * STRIDE ))

    while [ $(jobs -r | wc -l) -ge $N_WORKERS ]; do
        sleep 0.05
    done

    run_single_worker "$w" "$CURR_START" "$COUNT" "$STRIDE" "$CENTRAL_FILE" "$OUT_DIR" "$PROJECT_DIR" "$ROOT_CMD" &
    ACTIVE_WORKERS=$(( ACTIVE_WORKERS + 1 ))
done

wait
echo ""
echo "========================================================================"
echo " All $ACTIVE_WORKERS parallel MCMC workers completed."
echo " Worker logs saved in: $OUT_DIR/log_worker_*.txt"
echo " Entry logs saved in:  $OUT_DIR/log_entry_*.txt"
echo "========================================================================"
