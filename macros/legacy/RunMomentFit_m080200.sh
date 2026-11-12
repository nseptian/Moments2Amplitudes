#!/bin/bash
set -e

# ==============================================================================
# MomentFit Pipeline: Fit -> Monitor -> Calculate Acceptance -> Plot Moments
# ==============================================================================

# Configuration parameters with sensible defaults
MAX_L=${1:-4}
INIT_SCALE=${2:-0.2}
STEP=${3:-0.001}
TBIN_ARG=${4:-"all"}
MODE=${5:-"all"}   # "all" (default: fit + monitor) or "monitor_only" / "no_fit"

# 1. Ensure ROOT is available
if ! command -v root &>/dev/null; then
    if [ -f "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh" ]; then
        source "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh"
    fi
fi

# 2. Source BruFit environment
if [ -f "/d/grid17/septian/brufit_dev/brufit_cov/set_env.sh" ]; then
    source "/d/grid17/septian/brufit_dev/brufit_cov/set_env.sh"
elif [ -z "$BRUFIT" ]; then
    if [ -d "/d/home/septian/brufit" ]; then
        export BRUFIT="/d/home/septian/brufit"
    elif [ -d "/d/grid17/septian/brufit" ]; then
        export BRUFIT="/d/grid17/septian/brufit"
    fi
fi

# 3. Select Python environment with ROOT bindings
PYTHON_BIN="/d/grid17/septian/EtaPi0Analysis/.pixi/envs/default/bin/python"
if [ ! -f "$PYTHON_BIN" ]; then
    PYTHON_BIN="/d/home/septian/.pixi/envs/default/bin/python"
fi
if [ ! -f "$PYTHON_BIN" ]; then
    PYTHON_BIN="python3"
fi

# 4. Resolve Macro & Script paths
MACRO_PATH="/d/grid17/septian/EtaPi0Analysis/macros/MomentFit_m080200.C"
if [ ! -f "$MACRO_PATH" ]; then
    MACRO_PATH="/d/home/septian/Moments2Amplitudes/macros/MomentFit_m080200.C"
fi

MONITOR_SCRIPT="/d/grid17/septian/EtaPi0Analysis/macros/MonitorMomentFit.py"
if [ ! -f "$MONITOR_SCRIPT" ]; then
    MONITOR_SCRIPT="/d/home/septian/Moments2Amplitudes/macros/MonitorMomentFit.py"
fi

RUN_DIR="/d/grid17/septian/EtaPi0Analysis/run/MomentFit_m080200_FixedCov_Lmax${MAX_L}_InitValueScale${INIT_SCALE}_step${STEP}/"
LOGS_DIR="${RUN_DIR}logs/"

mkdir -p "$LOGS_DIR"

echo "=========================================================="
echo " Starting MomentFit_m080200 Pipeline (All-in-One)"
echo " Max L:            $MAX_L"
echo " Init Scale:       $INIT_SCALE"
echo " Step:             $STEP"
echo " t-bin selection:  $TBIN_ARG"
echo " Mode:             $MODE"
echo " BruFit Path:      $BRUFIT"
echo " Python Binary:    $PYTHON_BIN"
echo " Output Directory: $RUN_DIR"
echo " Logs Directory:   $LOGS_DIR"
echo "=========================================================="

if [ "$TBIN_ARG" == "all" ]; then
    TBIN_LIST=(0 1 2 3 4)
else
    TBIN_LIST=($TBIN_ARG)
fi

# Step 1: Run Fit Jobs (unless in monitor_only mode)
if [ "$MODE" != "monitor_only" ] && [ "$MODE" != "no_fit" ]; then
    PIDS=()
    for i in "${TBIN_LIST[@]}"; do
        LOG_FILE="${LOGS_DIR}MomentFit_m080200_tbin_${i}.log"
        echo "Launching MomentFit_m080200 for tbin_index = $i (MaxL = $MAX_L, Scale = $INIT_SCALE, Step = $STEP)"
        echo "  -> Log: $LOG_FILE"
        root -l -b -q "${BRUFIT}/macros/LoadBru.C" "${MACRO_PATH}($i, $MAX_L, $INIT_SCALE, $STEP)" > "$LOG_FILE" 2>&1 &
        PIDS+=($!)
    done

    echo ""
    echo "Fit job(s) running in background (PIDs: ${PIDS[*]})."
    echo "To monitor real-time fit logs in another terminal:"
    echo "  tail -f ${LOGS_DIR}MomentFit_m080200_tbin_*.log"
    echo ""
    echo "Waiting for all fit jobs to finish..."

    FAILED_JOBS=0
    for pid in "${PIDS[@]}"; do
        if ! wait "$pid"; then
            ((FAILED_JOBS++))
        fi
    done

    echo ""
    if [ "$FAILED_JOBS" -gt 0 ]; then
        echo "⚠️  $FAILED_JOBS fit job(s) encountered non-zero return code. Inspect logs in: $LOGS_DIR"
    else
        echo "✅ All MomentFit jobs completed!"
    fi
else
    echo "Skipping fit execution (MODE=$MODE). Running monitoring and downstream steps directly."
fi

# Step 2: Run Monitor, Acceptance Calculation & Plotting in One Go
echo ""
echo "=========================================================="
echo " Running Post-Fit Monitoring & Acceptance Pipeline"
echo "=========================================================="

"$PYTHON_BIN" "$MONITOR_SCRIPT" \
    --run-dir "$RUN_DIR" \
    --tbin "$TBIN_ARG" \
    --max-l "$MAX_L"

echo ""
echo "=========================================================="
echo " Pipeline Complete for: $RUN_DIR"
echo "=========================================================="
