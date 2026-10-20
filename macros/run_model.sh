#!/bin/bash

# Organized Model Run Script
# This script runs a modeling analysis and organizes all outputs in a timestamped directory

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================

# Default values
MACRO=""
SEED_START=1
SEED_END=100
MAX_JOBS=80
RUN_NAME=""
OUTPUT_BASE="runs"

# ============================================================================
# Functions
# ============================================================================

usage() {
    cat << EOF
Usage: $(basename $0) [OPTIONS] MACRO

Organized model run script - saves macro, logs, and results in one directory

Arguments:
  MACRO                 ROOT macro to run (e.g., H_odd_model.C)

Options:
  -n, --name NAME       Run name (default: MACRO_TIMESTAMP)
  -s, --start N         Starting seed number (default: 1)
  -e, --end N           Ending seed number (default: 100)
  -j, --jobs N          Max parallel jobs (default: 80)
  -o, --output DIR      Output base directory (default: runs)
  -h, --help            Show this help

Example:
  $(basename $0) -n test_run_001 -s 1 -e 500 -j 40 H_odd_model.C

Output structure:
  runs/
    test_run_001_20251119_143022/
      ├── macro/          # Copy of macro and dependencies
      ├── logs/           # Execution logs
      ├── results/        # Output ROOT files
      └── metadata.txt    # Run information

EOF
    exit 0
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n|--name)
                RUN_NAME="$2"
                shift 2
                ;;
            -s|--start)
                SEED_START="$2"
                shift 2
                ;;
            -e|--end)
                SEED_END="$2"
                shift 2
                ;;
            -j|--jobs)
                MAX_JOBS="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_BASE="$2"
                shift 2
                ;;
            -h|--help)
                usage
                ;;
            -*)
                echo "Error: Unknown option: $1" >&2
                usage
                ;;
            *)
                if [[ -z "$MACRO" ]]; then
                    MACRO="$1"
                else
                    echo "Error: Multiple macros specified" >&2
                    usage
                fi
                shift
                ;;
        esac
    done

    # Validate required arguments
    if [[ -z "$MACRO" ]]; then
        echo "Error: MACRO is required" >&2
        usage
    fi

}

# Setup run directory
setup_run_dir() {
    local macro_base=$(basename "$MACRO" .C)
    local timestamp=$(date +%Y%m%d_%H%M%S)
    
    # Set run name if not provided
    if [[ -z "$RUN_NAME" ]]; then
        RUN_NAME="${macro_base}_${timestamp}"
    else
        RUN_NAME="${RUN_NAME}_${timestamp}"
    fi
    
    # If OUTPUT_BASE is relative, make it relative to current directory
    # If we're in runs/, create directory here; otherwise in runs/ under project
    if [[ "$OUTPUT_BASE" = "runs" ]]; then
        if [[ "$(basename $(pwd))" = "runs" ]]; then
            # Already in runs directory, create here
            RUN_DIR="$(pwd)/${RUN_NAME}"
        else
            # Not in runs, create runs directory
            RUN_DIR="$(pwd)/${OUTPUT_BASE}/${RUN_NAME}"
        fi
    else
        RUN_DIR="$(cd "$(dirname "$OUTPUT_BASE")" 2>/dev/null && pwd)/$(basename "$OUTPUT_BASE")/${RUN_NAME}"
    fi
    
    # Create directory structure
    mkdir -p "$RUN_DIR"/{macro,logs,results}
    
    echo "Run directory created: $RUN_DIR"
}

# Copy macro and dependencies
copy_macros() {
    echo "Copying macros..."
    
    # Copy main macro
    cp "$PROJECT_DIR/macros/$MACRO" "$RUN_DIR/macro/"
    
    # Copy Load.C
    cp "$PROJECT_DIR/macros/Load.C" "$RUN_DIR/macro/"
    
    # Copy ConfigureAmpsNoValues.C if it exists and is used
    if [[ -f "$PROJECT_DIR/macros/ConfigureAmpsNoValues.C" ]] && grep -q "ConfigureAmpsNoValues" "$PROJECT_DIR/macros/$MACRO"; then
        cp "$PROJECT_DIR/macros/ConfigureAmpsNoValues.C" "$RUN_DIR/macro/"
    fi
    
    # Copy other included macros if needed
    for included_macro in $(grep -oP '#include\s+"\K[^"]+\.C' "$PROJECT_DIR/macros/$MACRO" 2>/dev/null || true); do
        if [[ -f "$PROJECT_DIR/macros/$included_macro" ]]; then
            cp "$PROJECT_DIR/macros/$included_macro" "$RUN_DIR/macro/"
        fi
    done
    
    echo "Macros copied to $RUN_DIR/macro/"
}

# Create metadata file
create_metadata() {
    cat > "$RUN_DIR/metadata.txt" << EOF
================================================================================
Model Run Metadata
================================================================================
Run Name:           $RUN_NAME
Date:               $(date '+%Y-%m-%d %H:%M:%S')
Host:               $(hostname)
User:               $(whoami)
Working Directory:  $(pwd)

Model Configuration:
  Macro:            $MACRO
  Seed Range:       $SEED_START to $SEED_END
  Total Jobs:       $((SEED_END - SEED_START + 1))
  Max Parallel:     $MAX_JOBS

Environment:
  BRUFIT:           ${BRUFIT:-Not set}
  ROOT Version:     $(root-config --version 2>/dev/null || echo "Unknown")
  ROOT Path:        $(which root 2>/dev/null || echo "Unknown")

Git Information:
  Branch:           $(git branch --show-current 2>/dev/null || echo "Not a git repo")
  Commit:           $(git rev-parse --short HEAD 2>/dev/null || echo "N/A")
  Status:           $(git status --short 2>/dev/null | wc -l) modified files

Output Structure:
  Run Directory:    $RUN_DIR
  Logs:             $RUN_DIR/logs/
  Results:          $RUN_DIR/results/
  Macros:           $RUN_DIR/macro/

================================================================================
EOF
    
    echo "Metadata saved to $RUN_DIR/metadata.txt"
}

# Run a single fit
run_fit() {
    local seed=$1
    local log_file="$RUN_DIR/logs/log_seed_${seed}.txt"
    
    # Run the fit with ROOT_INCLUDE_PATH set for ACLiC compilation
    # Output files will go to current directory (results/)
    # ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${ROOT_INCLUDE_PATH}" \
    root -l -b -q \
        "${BRUFIT}/macros/LoadBru.C" \
        "${PROJECT_DIR}/macros/Load.C" \
        "${PROJECT_DIR}/macros/${MACRO}(${seed})" \
        > "$log_file" 2>&1
    
    local status=$?
    
    if [[ $status -eq 0 ]]; then
        echo "[$(date '+%H:%M:%S')] ✓ Seed $seed completed"
    else
        echo "[$(date '+%H:%M:%S')] ✗ Seed $seed failed (exit code: $status)" >&2
    fi
}

# Run all fits
run_all_fits() {
    echo ""
    echo "========================================================================"
    echo "Starting Model Run: $RUN_NAME"
    echo "========================================================================"
    echo "Macro:        $MACRO"
    echo "Seeds:        $SEED_START to $SEED_END"
    echo "Total Jobs:   $((SEED_END - SEED_START + 1))"
    echo "Max Parallel: $MAX_JOBS"
    echo "Output:       $RUN_DIR"
    echo "========================================================================"
    echo ""
    
    local start_time=$(date +%s)
    local completed=0
    local failed=0
    
    # Change to results directory so all output files go there
    cd "$RUN_DIR/results"
    
    # Run fits in parallel
    for seed in $(seq $SEED_START $SEED_END); do
        # Wait if we've reached the maximum number of background jobs
        while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
            sleep 0.1
        done
        
        # Start the fit in background
        run_fit $seed &
        
        echo "[$(date '+%H:%M:%S')] Started seed $seed ($(jobs -r | wc -l) jobs running)"
    done
    
    # Wait for all background jobs to complete
    echo ""
    echo "All jobs launched. Waiting for completion..."
    wait
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # Count successes and failures
    local total_jobs=$((SEED_END - SEED_START + 1))
    failed=$(grep -l "Error\|error\|ERROR\|failed" "$RUN_DIR/logs"/*.txt 2>/dev/null | wc -l)
    completed=$((total_jobs - failed))
    
    echo ""
    echo "========================================================================"
    echo "Run Completed: $RUN_NAME"
    echo "========================================================================"
    echo "Duration:     $((duration / 60))m $((duration % 60))s"
    echo "Completed:    $completed/$total_jobs jobs"
    echo "Failed:       $failed/$total_jobs jobs"
    echo "Output:       $RUN_DIR"
    echo "========================================================================"
    echo ""
    
    # Save completion info
    cat >> "$RUN_DIR/metadata.txt" << EOF

Run Completion:
  Completion Time:  $(date '+%Y-%m-%d %H:%M:%S')
  Duration:         $((duration / 60))m $((duration % 60))s
  Successful:       $completed
  Failed:           $failed
  Total:            $total_jobs

================================================================================
EOF
}

# ============================================================================
# Main
# ============================================================================

main() {
    # Determine project directory (script might be run from runs/ or project root)
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
    export PROJECT_DIR
    
    # Parse arguments
    parse_args "$@"
    
    # Check if macro exists
    if [[ ! -f "$PROJECT_DIR/macros/$MACRO" ]]; then
        echo "Error: Macro not found: $PROJECT_DIR/macros/$MACRO" >&2
        exit 1
    fi
    
    # Source BruFit environment
    if [[ -f "/d/grid17/septian/brufit_dev/brufit/set_env.sh" ]]; then
        source /d/grid17/septian/brufit_dev/brufit/set_env.sh
    else
        echo "Warning: BruFit environment not found" >&2
    fi
    
    # Check BRUFIT is set
    if [[ -z "$BRUFIT" ]]; then
        echo "Error: BRUFIT environment variable not set" >&2
        exit 1
    fi
    
    # Setup run directory
    setup_run_dir
    
    # Copy macros
    copy_macros
    
    # Create metadata
    create_metadata
    
    # Export variables for subshells
    export BRUFIT RUN_DIR PROJECT_DIR MACRO
    
    # Run all fits
    run_all_fits
    
    echo "Results are in: $RUN_DIR"
}

# Run main function
main "$@"
