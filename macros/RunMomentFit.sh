#!/bin/bash
set -e

# ==============================================================================
# RunMomentFit.sh
#
# Generalized MomentFit Pipeline (arbitrary mass ranges via --config <yaml> or CLI)
# powered by MomentFit.C and MomentFitCollector.py:
#   1. Runs BruFit MomentFit (MomentFit.C) across t-bins in parallel.
#   2. Process control: cancel running jobs, detach/background, live watch.
#   3. Runs post-fit monitoring and diagnostics (MonitorMomentFit.py).
#   4. Collects unsuccessful bins (fit failure, stuck chains, or ESS < threshold).
#   5. Diagnoses failure mechanisms (zero acceptance Ph1, stuck Ph2, or low ESS).
#   6. Generates ready-to-run targeted single-bin retry script (RunMomentFit_retry.sh).
#   7. Optional --auto-retry loop across configurable retry rounds.
# ==============================================================================

# Default Configuration Parameters (overridden by --config <yaml> or CLI flags)
MAX_L=8
INIT_SCALE=1.0
STEP=1.0
TBIN_ARG="all"
MODE="all"          # "all" (fit + wait + monitor + collect), "monitor_only", "collect_only", "cancel", "fit_only"
M_MIN=0.8
M_MAX=2.0
N_BINS=30
MCMC_ITERS=200000
MCMC_BURNIN=2000
ESS_THRESHOLD=20.0
ESS_METRIC="mean"
MIN_ESS_THRESHOLD="10.0"
AUTO_RETRY=false
MAX_RETRIES=2
NUM_CPUS=$(nproc 2>/dev/null || grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 2)
DEFAULT_MAX_PARALLEL=$(( NUM_CPUS / 2 ))
[ "$DEFAULT_MAX_PARALLEL" -lt 1 ] && DEFAULT_MAX_PARALLEL=1
MAX_PARALLEL="${MAX_PARALLEL:-$DEFAULT_MAX_PARALLEL}"
# New Likelihood & PDF options
APPLY_MC_VAR=false
ALLOW_NEG_PDF=true
USE_TUNING=true
GIBBS_BLOCK_SIZE=0
COLD_WARMUP_ACC=0.15
MIN_ACCEPTANCE=0.05
MAX_ACCEPTANCE=0.15
TARGET_ACCEPTANCE=""
TUNE_MODE="MappedRhat"
RHAT_WINDOW=1000
RHAT_TARGET=1.05
ESS_FRACTION=0.02
RHAT_MAX_RETRIES=5
BAILOUT_THRESHOLD=10000
BAILOUT_MIN_ACC=0.0001
MC_INT_BINS=50
FIXED_COV=true
USE_COV=true
SET_INIT_VALUE=true
USE_GJ_FRAME=true
NUM_CPU=1
MCMC_SEED=111
MAX_RETRY_ITERS=200000
BURNIN_FRACTION=0.10
LOW_ACC_THRESH=0.03
GIBBS_LOW_ACC_THRESH=0.12
HIGH_ACC_THRESH=0.35
GIBBS_ACT_ACC=0.05
RECONDITION_COV=false
RECONDITION_SHRINKAGE=0.2
RECONDITION_MIN_EIGEN=1.0e-4

NO_WAIT=false
WATCH=false
FORCE=false
OVERWRITE=false
DRY_RUN=false
INTERVAL=10.0
CUSTOM_RUN_DIR=""
CUSTOM_COV_DIR=""
CUSTOM_DATA_DIR=""
CUSTOM_DATA_FILE=""
CUSTOM_MCRECON_FILE=""
CUSTOM_MCGEN_FILE=""
CUSTOM_WEIGHT_FILE=""
CONFIG_FILE=""
ACCEPTANCE_MACRO=""
SKIP_ACCEPTANCE=false
FORCE_ACCEPTANCE=false
EXTRA_MONITOR_ARGS=()

EXPLICIT_STEP=false
EXPLICIT_SCALE=false
EXPLICIT_MAX_L=false
EXPLICIT_M_MIN=false
EXPLICIT_M_MAX=false
EXPLICIT_N_BINS=false
EXPLICIT_ESS_THRESHOLD=false
EXPLICIT_ESS_METRIC=false
EXPLICIT_MIN_ESS_THRESHOLD=false
EXPLICIT_MAX_RETRIES=false
EXPLICIT_MAX_PARALLEL=false

declare -A TBIN_INDEX_MAP=(
    [0]="t010020"
    [1]="t020032"
    [2]="t032050"
    [3]="t050075"
    [4]="t075100"
)

print_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] [MAX_L] [INIT_SCALE] [STEP] [TBIN] [MODE]

Options:
  --monitor-only, --monitor, -m   Monitor progress only without launching fits or waiting.
  --collect-only, collect         Run diagnostics collection & retry generation only.
  --cancel, --stop                Cancel running fit jobs (for specified t-bin or all).
  --force, -f                     Force terminate immediately using SIGKILL (with --cancel).
  --dry-run                       Preview which processes would be canceled without killing.
  --no-wait, --detach, -b         Launch fit jobs in background and return immediately.
  --watch, -w                     Continuously refresh progress table (Ctrl+C to exit).
  --interval <sec>                Refresh interval for --watch (default: 10s).
  --run-dir, -d <path>            Specify run directory directly.
  --config, -c <yaml>             Route execution through a MomentFit YAML configuration file.
  --cov-dir <path>                Specify linear algebra covariance directory.
  --data-dir <path>               Specify input data directory.
  --data-pattern <pattern>        Specify data file path or pattern (e.g. {data_dir}/data_merged_{tbin}.root).
  --mcrecon-pattern <pattern>     Specify mcrecon file path or pattern (e.g. {data_dir}/reconmc_merged_{tbin}.root).
  --mcgen-pattern <pattern>       Specify mcgen file path or pattern (e.g. {data_dir}/genmc_merged_{tbin}.root).
  --weight-pattern <pattern>      Specify weight file path or pattern (e.g. {data_dir}/weights_data_merged_{tbin}.root).
  --tbin, -t <bin>                t-bin selection (0..4, or 'all', default: all).
  --max-l, -l <L>                 Maximum L multipole (default: 4).
  --m-min <float>                 Minimum mass range in GeV (default: 0.8).
  --m-max <float>                 Maximum mass range in GeV (default: 2.0).
  --n-bins <int>                  Number of mass bins (default: 30).
  --init-scale, -s <scale>        Initial value scale (default: 1.0).
  --step <step>                   MCMC step size (default: 1.0).
  --iters <int>                   MCMC iterations (default: 10000).
  --burn-in <int>                 MCMC burn-in (default: 2000).
  --ess-threshold <float>         Effective Sample Size threshold (default: 20.0).
  --ess-metric <mean|min>         ESS metric to evaluate against threshold (default: mean).
  --auto-retry                    Automatically execute generated retries until ESS >= threshold.
  --max-retries <int>             Maximum automatic retry rounds (default: 2).
  --max-parallel, -j <int>        Concurrency limit for single-bin refits (default: half of CPUs, currently $DEFAULT_MAX_PARALLEL).
  --plot, --force-plot, --plot-only Draw trace & posterior plots without running fits.
  --downsample-trace [N]          Downsample MCMC trace plots to N points (default: disabled; N=300 if flag given without value).
  --replot-from-root, --from-root Re-render all PDF plots directly from saved ROOT files (fast cosmetics mode).
  --overwrite, -y                 Overwrite existing fit logs and results without prompting.
  --mc-variance                   Enable MC variance scaling (default: true).
  --no-mc-variance                Disable MC variance scaling.
  --negative-pdf                  Allow negative PDF / negative intensities (default: true).
  --no-negative-pdf               Disable negative PDF.
  --tuning                        Activate MCMC covariance tuning engine (default: true).
  --no-tuning                     Disable MCMC covariance tuning engine.
  --gibbs-block-size <int>        Gibbs block size (default: 0).
  --cold-warmup-acc <float>       Cold warm-up acceptance target (default: 0.15).
  --min-acceptance <float>        Minimum acceptance bound for tuning (default: 0.05).
  --max-acceptance <float>        Maximum acceptance bound for tuning (default: 0.15).
  --target-acceptance <float>     Target acceptance rate (default: midpoint of [min, max]).
  --tune-mode <mode>              Tuning mode: MappedRhat or Acceptance (default: MappedRhat).
  --rhat-window <int>             R-hat tuning window size (default: 1000).
  --rhat-target <float>           R-hat tuning target (default: 1.05).
  --ess-fraction <float>          ESS fraction for MappedRhat tuning (default: 0.02).
  --rhat-max-retries <int>        Max adaptive R-hat retries (default: 5).
  --bailout-threshold <int>       Bailout check proposal count (default: 10000).
  --bailout-min-acceptance <flt>  Bailout min acceptance fraction (default: 0.0001).
  --mc-integration-bins <int>     Angular MC integration bins (default: 50).
  --fixed-covariance / --no-fixed-covariance  Use linear-algebra fixed covariance (default: true).
  --use-covariance / --no-use-covariance      Use covariance proposal matrix (default: true).
  --set-init-value / --no-set-init-value      Initialize from Stage 2 moments (default: true).
  --gj-frame / --hel-frame        Decay frame (default: --gj-frame).
  --num-cpu <int>                 RooFit NumCPU per job (default: 1).
  --seed <int>                    gRandom RNG seed (default: 111).
  --max-retry-iters <int>         Retry MCMC iteration cap (default: 200000).
  --burnin-fraction <float>       Retry burn-in fraction of iterations (default: 0.10).
  --recondition-covariance / --no-recondition-covariance  Opt-in empirical covariance reconditioning (default: false).
  --calculate-acceptance [path]   Run acceptance calculation directly (with optional macro path).
  --skip-acceptance               Skip running acceptance calculation in MonitorMomentFit.py.
  --force-acceptance              Force running acceptance calculation even if not all bins converged.
  --progress-only, --status-only  Print status table without full plot generation.
  --help, -h                      Show this help message.

Examples:
  $(basename "$0")                               # Full fit, monitor, and diagnostic collection
  $(basename "$0") --force-plot --run-dir <dir>  # Draw diagnostic plots only (no fits)
  $(basename "$0") --monitor-only                # Inspect current progress anytime
  $(basename "$0") --watch                       # Continuously watch progress in real-time
  $(basename "$0") --collect-only                # Check diagnostics & generate retries for bad bins
  $(basename "$0") --auto-retry                  # Run fit and auto-retry until all bins meet ESS >= 100
  $(basename "$0") --cancel                      # Cancel all running fit jobs
  $(basename "$0") --cancel --tbin 0             # Cancel running fit for t-bin 0 only
  $(basename "$0") 4 0.2 0.1 all --no-wait       # Launch fits in background without waiting
EOF
}

# ------------------------------------------------------------------------------
# Flexible Argument Parsing (supports flags, positional args, and shortcuts)
# ------------------------------------------------------------------------------
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cancel|--stop|cancel|stop|--kill|kill)
            MODE="cancel"
            shift
            ;;
        --force|-f)
            FORCE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --monitor-only|--monitor|-m|monitor_only|monitor|no_fit|--no-fit|status|progress)
            MODE="monitor_only"
            shift
            ;;
        --collect-only|collect_only|collect)
            MODE="collect_only"
            shift
            ;;
        --no-wait|--detach|-b|no_wait|detach|fit_only|--fit-only)
            NO_WAIT=true
            shift
            ;;
        --watch|-w)
            WATCH=true
            MODE="monitor_only"
            shift
            ;;
        --interval)
            INTERVAL="$2"
            shift 2
            ;;
        --run-dir|-d)
            CUSTOM_RUN_DIR="$2"
            shift 2
            ;;
        --cov-dir)
            CUSTOM_COV_DIR="$2"
            shift 2
            ;;
        --data-dir)
            CUSTOM_DATA_DIR="$2"
            shift 2
            ;;
        --data-pattern|--data-file)
            CUSTOM_DATA_FILE="$2"
            shift 2
            ;;
        --mcrecon-pattern|--mcrecon-file|--reconmc-pattern|--reconmc-file)
            CUSTOM_MCRECON_FILE="$2"
            shift 2
            ;;
        --mcgen-pattern|--mcgen-file|--genmc-pattern|--genmc-file)
            CUSTOM_MCGEN_FILE="$2"
            shift 2
            ;;
        --weight-pattern|--weight-file|--weights-pattern|--weights-file)
            CUSTOM_WEIGHT_FILE="$2"
            shift 2
            ;;
        --tbin|-t)
            TBIN_ARG="$2"
            shift 2
            ;;
        --max-l|-l)
            MAX_L="$2"
            EXPLICIT_MAX_L=true
            shift 2
            ;;
        --m-min|--mmin)
            M_MIN="$2"
            EXPLICIT_M_MIN=true
            shift 2
            ;;
        --m-max|--mmax)
            M_MAX="$2"
            EXPLICIT_M_MAX=true
            shift 2
            ;;
        --n-bins|--nbins)
            N_BINS="$2"
            EXPLICIT_N_BINS=true
            shift 2
            ;;
        --init-scale|-s)
            INIT_SCALE="$2"
            EXPLICIT_SCALE=true
            shift 2
            ;;
        --step)
            STEP="$2"
            EXPLICIT_STEP=true
            shift 2
            ;;
        --iters)
            MCMC_ITERS="$2"
            shift 2
            ;;
        --burn-in)
            MCMC_BURNIN="$2"
            shift 2
            ;;
        --ess-threshold)
            ESS_THRESHOLD="$2"
            EXPLICIT_ESS_THRESHOLD=true
            shift 2
            ;;
        --ess-metric)
            ESS_METRIC="$2"
            EXPLICIT_ESS_METRIC=true
            shift 2
            ;;
        --min-ess-threshold)
            MIN_ESS_THRESHOLD="$2"
            EXPLICIT_MIN_ESS_THRESHOLD=true
            shift 2
            ;;
        --mc-variance)
            APPLY_MC_VAR=true
            shift
            ;;
        --no-mc-variance)
            APPLY_MC_VAR=false
            shift
            ;;
        --negative-pdf)
            ALLOW_NEG_PDF=true
            shift
            ;;
        --no-negative-pdf)
            ALLOW_NEG_PDF=false
            shift
            ;;
        --tuning)
            USE_TUNING=true
            shift
            ;;
        --no-tuning)
            USE_TUNING=false
            shift
            ;;
        --gibbs-block-size)
            GIBBS_BLOCK_SIZE="$2"
            shift 2
            ;;
        --cold-warmup-acc)
            COLD_WARMUP_ACC="$2"
            shift 2
            ;;
        --min-acceptance)
            MIN_ACCEPTANCE="$2"
            shift 2
            ;;
        --max-acceptance)
            MAX_ACCEPTANCE="$2"
            shift 2
            ;;
        --target-acceptance)
            TARGET_ACCEPTANCE="$2"
            shift 2
            ;;
        --tune-mode)
            TUNE_MODE="$2"
            shift 2
            ;;
        --rhat-window)
            RHAT_WINDOW="$2"
            shift 2
            ;;
        --rhat-target)
            RHAT_TARGET="$2"
            shift 2
            ;;
        --ess-fraction)
            ESS_FRACTION="$2"
            shift 2
            ;;
        --rhat-max-retries)
            RHAT_MAX_RETRIES="$2"
            shift 2
            ;;
        --bailout-threshold)
            BAILOUT_THRESHOLD="$2"
            shift 2
            ;;
        --bailout-min-acceptance)
            BAILOUT_MIN_ACC="$2"
            shift 2
            ;;
        --mc-integration-bins)
            MC_INT_BINS="$2"
            shift 2
            ;;
        --fixed-covariance)
            FIXED_COV=true
            shift
            ;;
        --no-fixed-covariance)
            FIXED_COV=false
            shift
            ;;
        --use-covariance)
            USE_COV=true
            shift
            ;;
        --no-use-covariance)
            USE_COV=false
            shift
            ;;
        --set-init-value)
            SET_INIT_VALUE=true
            shift
            ;;
        --no-set-init-value)
            SET_INIT_VALUE=false
            shift
            ;;
        --gj-frame)
            USE_GJ_FRAME=true
            shift
            ;;
        --hel-frame)
            USE_GJ_FRAME=false
            shift
            ;;
        --num-cpu)
            NUM_CPU="$2"
            shift 2
            ;;
        --seed)
            MCMC_SEED="$2"
            shift 2
            ;;
        --max-retry-iters)
            MAX_RETRY_ITERS="$2"
            shift 2
            ;;
        --burnin-fraction)
            BURNIN_FRACTION="$2"
            shift 2
            ;;
        --low-acceptance-threshold)
            LOW_ACC_THRESH="$2"
            shift 2
            ;;
        --gibbs-low-acceptance-threshold)
            GIBBS_LOW_ACC_THRESH="$2"
            shift 2
            ;;
        --high-acceptance-threshold)
            HIGH_ACC_THRESH="$2"
            shift 2
            ;;
        --gibbs-activation-acceptance)
            GIBBS_ACT_ACC="$2"
            shift 2
            ;;
        --recondition-covariance)
            RECONDITION_COV=true
            shift
            ;;
        --no-recondition-covariance)
            RECONDITION_COV=false
            shift
            ;;
        --recondition-shrinkage)
            RECONDITION_SHRINKAGE="$2"
            shift 2
            ;;
        --recondition-min-eigen-floor)
            RECONDITION_MIN_EIGEN="$2"
            shift 2
            ;;
        --auto-retry)
            AUTO_RETRY=true
            shift
            ;;
        --max-retries)
            MAX_RETRIES="$2"
            EXPLICIT_MAX_RETRIES=true
            shift 2
            ;;
        --max-parallel|-j)
            MAX_PARALLEL="$2"
            EXPLICIT_MAX_PARALLEL=true
            shift 2
            ;;
        --overwrite|--force-overwrite|-y|--yes)
            OVERWRITE=true
            shift
            ;;
        --plot|--force-plot|--force-plots|--plot-only)
            EXTRA_MONITOR_ARGS+=("--plot")
            if [ "$MODE" = "all" ]; then
                MODE="plot_only"
            fi
            shift
            ;;
        --skip-trace|--skip-posterior|--skip-moments)
            EXTRA_MONITOR_ARGS+=("$1")
            shift
            ;;
        --downsample-trace|--max-trace-points)
            if [ -n "$2" ] && [[ "$2" =~ ^[0-9]+$ ]]; then
                EXTRA_MONITOR_ARGS+=("--downsample-trace" "$2")
                shift 2
            else
                EXTRA_MONITOR_ARGS+=("--downsample-trace")
                shift 1
            fi
            ;;
        --replot-from-root|--from-root)
            EXTRA_MONITOR_ARGS+=("--replot-from-root")
            if [ "$MODE" = "all" ]; then
                MODE="plot_only"
            fi
            shift
            ;;
        --progress-only|--status-only)
            EXTRA_MONITOR_ARGS+=("--progress-only")
            shift
            ;;
        --config|-c)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --calculate-acceptance)
            MODE="calculate_acceptance"
            if [ -n "${2:-}" ] && [[ "${2:-}" != -* ]] && [[ "${2:-}" != *.yaml ]] && [[ "${2:-}" != *.yml ]]; then
                if [[ "${2:-}" == *.C ]] || [ -f "${2:-}" ]; then
                    ACCEPTANCE_MACRO="$2"
                    shift 2
                    continue
                fi
            fi
            shift
            ;;
        --skip-acceptance)
            SKIP_ACCEPTANCE=true
            EXTRA_MONITOR_ARGS+=("--skip-acceptance")
            shift
            ;;
        --force-acceptance)
            FORCE_ACCEPTANCE=true
            EXTRA_MONITOR_ARGS+=("--force-acceptance")
            shift
            ;;
        --retry|retry)
            MODE="retry_only"
            shift
            ;;
        --help|-h)
            print_usage
            exit 0
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Parse any remaining positional arguments
NUM_POS=${#POSITIONAL_ARGS[@]}
if [ "$NUM_POS" -ge 1 ]; then
    case "${POSITIONAL_ARGS[0]}" in
        cancel|stop|kill)
            MODE="cancel"
            if [ "$NUM_POS" -ge 2 ]; then
                TBIN_ARG="${POSITIONAL_ARGS[1]}"
            fi
            ;;
        plot|plot_only|plots)
            MODE="plot_only"
            EXTRA_MONITOR_ARGS+=("--plot")
            if [ "$NUM_POS" -ge 2 ]; then
                TBIN_ARG="${POSITIONAL_ARGS[1]}"
            fi
            ;;
        *)
            if [[ "${POSITIONAL_ARGS[0]}" =~ ^([0-4]|all|t[0-9]{6})$ ]] && [ "$NUM_POS" -le 2 ]; then
                TBIN_ARG="${POSITIONAL_ARGS[0]}"
                if [ "$NUM_POS" -ge 2 ]; then
                    case "${POSITIONAL_ARGS[1]}" in
                        cancel|stop|kill) MODE="cancel" ;;
                        monitor_only|monitor|no_fit) MODE="monitor_only" ;;
                        collect_only|collect) MODE="collect_only" ;;
                        plot_only|plot|plots) MODE="plot_only"; EXTRA_MONITOR_ARGS+=("--plot") ;;
                        *) MODE="${POSITIONAL_ARGS[1]}" ;;
                    esac
                fi
            elif [ "$NUM_POS" -ge 3 ] && [[ "${POSITIONAL_ARGS[2]}" =~ ^([0-4]|all|t[0-9]{6})$ ]]; then
                MAX_L="${POSITIONAL_ARGS[0]}"; EXPLICIT_MAX_L=true
                INIT_SCALE="${POSITIONAL_ARGS[1]}"; EXPLICIT_SCALE=true
                STEP="1.0"
                TBIN_ARG="${POSITIONAL_ARGS[2]}"
                if [ "$NUM_POS" -ge 4 ]; then
                    case "${POSITIONAL_ARGS[3]}" in
                        cancel|stop|kill) MODE="cancel" ;;
                        monitor_only|monitor|no_fit) MODE="monitor_only" ;;
                        collect_only|collect) MODE="collect_only" ;;
                        plot_only|plot|plots) MODE="plot_only"; EXTRA_MONITOR_ARGS+=("--plot") ;;
                        *) MODE="${POSITIONAL_ARGS[3]}" ;;
                    esac
                fi
            else
                [ "$NUM_POS" -ge 1 ] && MAX_L="${POSITIONAL_ARGS[0]}" && EXPLICIT_MAX_L=true
                [ "$NUM_POS" -ge 2 ] && INIT_SCALE="${POSITIONAL_ARGS[1]}" && EXPLICIT_SCALE=true
                [ "$NUM_POS" -ge 3 ] && STEP="${POSITIONAL_ARGS[2]}" && EXPLICIT_STEP=true
                [ "$NUM_POS" -ge 4 ] && TBIN_ARG="${POSITIONAL_ARGS[3]}"
                if [ "$NUM_POS" -ge 5 ]; then
                    case "${POSITIONAL_ARGS[4]}" in
                        cancel|stop|kill) MODE="cancel" ;;
                        monitor_only|monitor|no_fit) MODE="monitor_only" ;;
                        collect_only|collect) MODE="collect_only" ;;
                        plot_only|plot|plots) MODE="plot_only"; EXTRA_MONITOR_ARGS+=("--plot") ;;
                        *) MODE="${POSITIONAL_ARGS[4]}" ;;
                    esac
                fi
            fi
            ;;
    esac
fi

# 1. Setup BruFit & Analysis Environment (skip if already set up)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# If --config <yaml> was supplied, route execution cleanly through the YAML-driven Python runner
if [ -n "$CONFIG_FILE" ]; then
    PYTHON_BIN="${SCRIPT_DIR}/../.pixi/envs/default/bin/python"
    [ ! -x "$PYTHON_BIN" ] && PYTHON_BIN="$(command -v python3 || command -v python)"
    RUNNER_ARGS=(--config "$CONFIG_FILE")
    case "$MODE" in
        monitor_only) RUNNER_ARGS+=(--monitor) ;;
        collect_only) RUNNER_ARGS+=(--collect) ;;
        retry_only) RUNNER_ARGS+=(--retry) ;;
        plot_only) RUNNER_ARGS+=(--plot) ;;
        calculate_acceptance)
            RUNNER_ARGS+=(--calculate-acceptance)
            [ -n "$ACCEPTANCE_MACRO" ] && RUNNER_ARGS+=("$ACCEPTANCE_MACRO")
            ;;
        cancel) RUNNER_ARGS+=(--cancel) ;;
        *) RUNNER_ARGS+=(--stage moment_fit) ;;
    esac
    [ "$EXPLICIT_MAX_L" = true ] && RUNNER_ARGS+=(--max-l "$MAX_L")
    [ "$TBIN_ARG" != "all" ] && RUNNER_ARGS+=(--tbin "$TBIN_ARG")
    [ -n "$CUSTOM_RUN_DIR" ] && RUNNER_ARGS+=(--run-dir "$CUSTOM_RUN_DIR")
    [ "$DRY_RUN" = true ] && RUNNER_ARGS+=(--dry-run)
    [ "$OVERWRITE" = true ] && RUNNER_ARGS+=(--overwrite)
    [ "$AUTO_RETRY" = true ] && RUNNER_ARGS+=(--auto-retry)
    [ "$EXPLICIT_ESS_THRESHOLD" = true ] && RUNNER_ARGS+=(--ess-threshold "$ESS_THRESHOLD")
    [ "$EXPLICIT_ESS_METRIC" = true ] && RUNNER_ARGS+=(--ess-metric "$ESS_METRIC")
    [ "$EXPLICIT_MIN_ESS_THRESHOLD" = true ] && RUNNER_ARGS+=(--min-ess-threshold "$MIN_ESS_THRESHOLD")
    [ "$EXPLICIT_MAX_RETRIES" = true ] && RUNNER_ARGS+=(--max-retries "$MAX_RETRIES")
    [ "$EXPLICIT_MAX_PARALLEL" = true ] && RUNNER_ARGS+=(--max-parallel "$MAX_PARALLEL")
    [ "$SKIP_ACCEPTANCE" = true ] && RUNNER_ARGS+=(--skip-acceptance)
    [ "$FORCE_ACCEPTANCE" = true ] && RUNNER_ARGS+=(--force-acceptance)
    [ "$WATCH" = true ] && RUNNER_ARGS+=(--watch --interval "$INTERVAL")
    [ "$FORCE" = true ] && RUNNER_ARGS+=(--force)
    exec env PYTHONPATH="${SCRIPT_DIR}/../python:${PYTHONPATH}" "$PYTHON_BIN" -m moments2amplitudes.cli momentfit "${RUNNER_ARGS[@]}"
fi

if [ -z "$BRUFIT" ]; then
    CURRENT_DIR="$(pwd)"
    if [ -f "$SCRIPT_DIR/../set_env.sh" ]; then
        source "$SCRIPT_DIR/../set_env.sh"
    elif [ -f "$SCRIPT_DIR/../set_ana_env.sh" ]; then
        source "$SCRIPT_DIR/../set_ana_env.sh"
    elif [ -f "/d/home/septian/EtaPi0Analysis/set_ana_env.sh" ]; then
        source "/d/home/septian/EtaPi0Analysis/set_ana_env.sh"
    elif [ -f "/d/grid17/septian/EtaPi0Analysis/set_ana_env.sh" ]; then
        source "/d/grid17/septian/EtaPi0Analysis/set_ana_env.sh"
    fi
    cd "$CURRENT_DIR" 2>/dev/null || true
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 2. Ensure ROOT is available
if ! command -v root &>/dev/null; then
    if [ -d "/d/grid17/septian/root_v6-36-08_cuda_install/bin" ]; then
        export PATH="/d/grid17/septian/root_v6-36-08_cuda_install/bin:$PATH"
    fi
    if [ -f "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh" ]; then
        source "/d/grid17/septian/root_v6-36-08_cuda_install/bin/thisroot.sh"
    fi
fi

# 3. Select Python 3.12 environment matching ROOT 6.36 / BruFit
if [ -f "$SCRIPT_DIR/../.pixi/envs/default/bin/python" ] && [ "$("$SCRIPT_DIR/../.pixi/envs/default/bin/python" -c "import sys; print(sys.version_info[0:2] >= (3, 12))" 2>/dev/null)" = "True" ]; then
    PYTHON_BIN="$SCRIPT_DIR/../.pixi/envs/default/bin/python"
elif [ -f "$SCRIPT_DIR/../.pixi/envs/py312/bin/python" ]; then
    PYTHON_BIN="$SCRIPT_DIR/../.pixi/envs/py312/bin/python"
elif [ -f "/d/home/septian/EtaPi0Analysis/.pixi/envs/default/bin/python" ]; then
    PYTHON_BIN="/d/home/septian/EtaPi0Analysis/.pixi/envs/default/bin/python"
elif [ -f "/d/grid17/septian/EtaPi0Analysis/.pixi/envs/default/bin/python" ]; then
    PYTHON_BIN="/d/grid17/septian/EtaPi0Analysis/.pixi/envs/default/bin/python"
elif command -v python3.12 &>/dev/null; then
    PYTHON_BIN="$(command -v python3.12)"
elif [ -f "$SCRIPT_DIR/../.pixi/envs/default/bin/python" ]; then
    PYTHON_BIN="$SCRIPT_DIR/../.pixi/envs/default/bin/python"
else
    PYTHON_BIN="python3"
fi

# Ensure ROOT 6.36 PyROOT bindings and local python libraries are in PYTHONPATH
CLEAN_PYTHONPATH=$(echo "$PYTHONPATH" | tr ':' '\n' | grep -v 'lib/python3\.' | tr '\n' ':' | sed 's/:$//')
export PYTHONPATH="${ROOTSYS}/lib:${SCRIPT_DIR}/../python:${CLEAN_PYTHONPATH}"

# 4. Resolve Macro & Script paths (using local macros in Moments2Amplitudes)
MACRO_PATH="$SCRIPT_DIR/MomentFit.C"
MONITOR_SCRIPT="$SCRIPT_DIR/MonitorMomentFit.py"
COLLECTOR_SCRIPT="$SCRIPT_DIR/MomentFitCollector.py"
if [ -z "$ACCEPTANCE_MACRO" ]; then
    ACCEPTANCE_MACRO="${MOM2AMP_CALCULATE_ACCEPTANCE:-$SCRIPT_DIR/CalculateAcceptance.C}"
fi

# 5. Resolve Output Directory & Mass Binning Auto-Detection
DIR_FOR_DETECTION="${CUSTOM_RUN_DIR:-}"
if [ -n "$DIR_FOR_DETECTION" ]; then
    if [ "$EXPLICIT_M_MIN" = false ] || [ "$EXPLICIT_M_MAX" = false ]; then
        if [[ "$DIR_FOR_DETECTION" =~ m([0-9]{3})([0-9]{3}) ]]; then
            TAG_MIN="${BASH_REMATCH[1]}"
            TAG_MAX="${BASH_REMATCH[2]}"
            [ "$EXPLICIT_M_MIN" = false ] && M_MIN=$(awk "BEGIN {printf \"%.4f\", $TAG_MIN / 100.0}")
            [ "$EXPLICIT_M_MAX" = false ] && M_MAX=$(awk "BEGIN {printf \"%.4f\", $TAG_MAX / 100.0}")
            if [ "$EXPLICIT_N_BINS" = false ]; then
                N_BINS=$(awk "BEGIN {printf \"%d\", int(($TAG_MAX - $TAG_MIN) / 4.0 + 0.5)}")
            fi
        elif [[ "$DIR_FOR_DETECTION" =~ m([0-9]{2})([0-9]{2}) ]]; then
            TAG_MIN="${BASH_REMATCH[1]}"
            TAG_MAX="${BASH_REMATCH[2]}"
            [ "$EXPLICIT_M_MIN" = false ] && M_MIN=$(awk "BEGIN {printf \"%.4f\", $TAG_MIN / 10.0}")
            [ "$EXPLICIT_M_MAX" = false ] && M_MAX=$(awk "BEGIN {printf \"%.4f\", $TAG_MAX / 10.0}")
            if [ "$EXPLICIT_N_BINS" = false ]; then
                N_BINS=$(awk "BEGIN {printf \"%d\", int(($M_MAX - $M_MIN) / 0.04 + 0.5)}")
            fi
        fi
    fi
fi

TAG_MIN_INT=$(awk "BEGIN {print int($M_MIN * 100 + 0.5)}")
TAG_MAX_INT=$(awk "BEGIN {print int($M_MAX * 100 + 0.5)}")
MASS_TAG=$(printf "m%03d%03d" "$TAG_MIN_INT" "$TAG_MAX_INT")

if [ -n "$CUSTOM_RUN_DIR" ]; then
    RUN_DIR="$CUSTOM_RUN_DIR"
else
    REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
    DEFAULT_RUN_DIR="${REPO_DIR}/run/MomentFit_${MASS_TAG}_FixedCov_Lmax${MAX_L}_InitValueScale${INIT_SCALE}_step${STEP}/"
    if [ "$MODE" != "all" ] && [ ! -d "$DEFAULT_RUN_DIR" ]; then
        LATEST_CANDIDATE=$(ls -td "${REPO_DIR}/run/MomentFit_${MASS_TAG}_FixedCov_*" "${REPO_DIR}/run/MomentFit_"*"_FixedCov_*" 2>/dev/null | grep -v "_stepall" | head -n 1 || true)
        if [ -n "$LATEST_CANDIDATE" ] && [ -d "$LATEST_CANDIDATE" ]; then
            RUN_DIR="$LATEST_CANDIDATE"
            echo "Auto-detected newest run directory: $RUN_DIR"
        else
            RUN_DIR="$DEFAULT_RUN_DIR"
        fi
    else
        RUN_DIR="$DEFAULT_RUN_DIR"
    fi
fi
[ "${RUN_DIR: -1}" != "/" ] && RUN_DIR="${RUN_DIR}/"
LOGS_DIR="${RUN_DIR}logs/"

if [ "$TBIN_ARG" == "all" ]; then
    TBIN_LIST=(0 1 2 3 4)
else
    TBIN_LIST=($TBIN_ARG)
fi

# ==============================================================================
# Cancellation Helper Function
# ==============================================================================
kill_tbin_fit() {
    local tbin_idx="$1"
    local pid_file="${LOGS_DIR}tbin_${tbin_idx}.pid"
    local log_file="${LOGS_DIR}MomentFit_${MASS_TAG}_tbin_${tbin_idx}.log"
    if [ ! -f "$log_file" ]; then
        local alt_log
        alt_log=$(ls -1 "${LOGS_DIR}"MomentFit*_tbin_"${tbin_idx}".log 2>/dev/null | head -n 1 || true)
        [ -n "$alt_log" ] && log_file="$alt_log"
    fi
    local tbin_name="${TBIN_INDEX_MAP[$tbin_idx]:-tbin_${tbin_idx}}"
    local target_pids=()

    # 1. Read PID file
    if [ -f "$pid_file" ]; then
        local spid
        spid=$(cat "$pid_file" 2>/dev/null || true)
        if [ -n "$spid" ] && kill -0 "$spid" 2>/dev/null; then
            target_pids+=("$spid")
            local children
            children=$(pgrep -P "$spid" 2>/dev/null || true)
            for c in $children; do
                target_pids+=("$c")
                local grandchildren
                grandchildren=$(pgrep -P "$c" 2>/dev/null || true)
                for gc in $grandchildren; do
                    target_pids+=("$gc")
                done
            done
        fi
    fi

    # 2. Search processes matching MomentFit macro invocation
    local pattern="MomentFit\.C\(${tbin_idx},"
    local matched
    matched=$(pgrep -u "$USER" -f "$pattern" 2>/dev/null || true)
    for p in $matched; do
        target_pids+=("$p")
        local children
        children=$(pgrep -P "$p" 2>/dev/null || true)
        for c in $children; do
            target_pids+=("$c")
        done
    done

    # Exclude our own PID
    local unique_pids=()
    for p in $(printf "%s\n" "${target_pids[@]}" 2>/dev/null | sort -u); do
        if [ "$p" != "$$" ]; then
            unique_pids+=("$p")
        fi
    done

    if [ ${#unique_pids[@]} -eq 0 ]; then
        echo "  ℹ️  No active processes found for tbin ${tbin_idx} (${tbin_name})."
        rm -f "$pid_file"
        return 0
    fi

    echo "  Found ${#unique_pids[@]} active process(es) for tbin ${tbin_idx} (${tbin_name}): ${unique_pids[*]}"

    if [ "$DRY_RUN" = true ]; then
        echo "  [DRY-RUN] Would terminate: ${unique_pids[*]}"
        return 0
    fi

    if [ "$FORCE" = true ]; then
        echo "  Force terminating with SIGKILL..."
        kill -9 "${unique_pids[@]}" 2>/dev/null || true
    else
        echo "  Sending SIGTERM for graceful shutdown..."
        kill -TERM "${unique_pids[@]}" 2>/dev/null || true
        sleep 1.5

        local alive=()
        for p in "${unique_pids[@]}"; do
            if kill -0 "$p" 2>/dev/null; then
                alive+=("$p")
            fi
        done
        if [ ${#alive[@]} -gt 0 ]; then
            echo "  Processes still alive after SIGTERM (${alive[*]}), sending SIGKILL..."
            kill -9 "${alive[@]}" 2>/dev/null || true
        fi
    fi

    local now_str
    now_str=$(date '+%Y-%m-%d %H:%M:%S')
    if [ -f "$log_file" ]; then
        echo "=== Job canceled by user at ${now_str} ===" >> "$log_file"
    fi

    local tbin_dir="${RUN_DIR}${tbin_name}/"
    if [ -d "$tbin_dir" ]; then
        for rlog in "${tbin_dir}"logRooFit*.txt; do
            if [ -f "$rlog" ] && [[ "$rlog" != *".bak"* ]]; then
                # Check if this bin already completed and has ResultsBruMcmcFixedCovariance.root
                already_done=false
                for mdir in "${tbin_dir}"Mpi0eta*; do
                    if [ -d "$mdir" ] && [[ "$mdir" != *".bak"* ]]; then
                        cand_root="${mdir}/ResultsBruMcmcFixedCovariance.root"
                        if [ -f "$cand_root" ] && [ "$(stat -c%s "$cand_root" 2>/dev/null || echo 0)" -gt 10240 ]; then
                            if grep -q "--> Target Bin.*$(basename "$mdir")" "$rlog" 2>/dev/null && grep -q "PosteriorCovariance" "$rlog" 2>/dev/null; then
                                already_done=true
                                break
                            fi
                        fi
                    fi
                done
                if [ "$already_done" = true ]; then
                    continue
                fi

                if ! grep -q "PosteriorCovariance" "$rlog" 2>/dev/null && ! grep -q "canceled by user" "$rlog" 2>/dev/null; then
                    echo "=== MCMC fit canceled by user at ${now_str} ===" >> "$rlog"
                fi
            fi
        done
    fi

    rm -f "$pid_file"
    echo "  ✅ Successfully canceled tbin ${tbin_idx} (${tbin_name})."
}

# ==============================================================================
# MODE: CANCEL
# ==============================================================================
if [ "$MODE" == "cancel" ]; then
    echo "=========================================================="
    echo " Cancelling MomentFit (${MASS_TAG}) Jobs"
    echo " Target t-bin(s):  $TBIN_ARG"
    echo " Run Directory:    $RUN_DIR"
    echo " Logs Directory:   $LOGS_DIR"
    echo " Mode:             $([ "$DRY_RUN" = true ] && echo "DRY-RUN" || ([ "$FORCE" = true ] && echo "FORCE (SIGKILL)" || echo "GRACEFUL (SIGTERM + SIGKILL)"))"
    echo "=========================================================="

    for i in "${TBIN_LIST[@]}"; do
        kill_tbin_fit "$i"
    done

    echo ""
    echo "Updated Pipeline Status:"
    "$PYTHON_BIN" "$MONITOR_SCRIPT" --run-dir "$RUN_DIR" --progress-only
    exit 0
fi

echo "=========================================================="
if [ "$MODE" == "monitor_only" ]; then
    echo " MomentFit (${MASS_TAG}) Monitoring Pipeline"
elif [ "$MODE" == "collect_only" ]; then
    echo " MomentFit (${MASS_TAG}) Diagnostic Collection & Retry Generation"
elif [ "$MODE" == "plot_only" ]; then
    echo " MomentFit (${MASS_TAG}) Diagnostic Plot Generation (No Fits)"
else
    echo " Starting MomentFit (${MASS_TAG}) Pipeline"
fi
echo " Max L:            $MAX_L"
echo " Init Scale:       $INIT_SCALE"
echo " Step:             $STEP"
echo " Mass Range:       [$M_MIN, $M_MAX] GeV ($N_BINS bins)"
echo " MCMC Samples:     $MCMC_ITERS (burn-in: $MCMC_BURNIN)"
echo " ESS Threshold:    $ESS_THRESHOLD"
echo " t-bin selection:  $TBIN_ARG"
echo " Mode:             $MODE"
echo " Auto-Retry:       $AUTO_RETRY (max retries: $MAX_RETRIES)"
echo " Max Parallel:     $MAX_PARALLEL jobs"
echo " BruFit Path:      ${BRUFIT:-Not set}"
echo " Python Binary:    $PYTHON_BIN"
echo " Macro Path:       $MACRO_PATH"
if [ -n "$CUSTOM_COV_DIR" ]; then
    echo " Cov Directory:    $CUSTOM_COV_DIR"
fi
if [ -n "$CUSTOM_DATA_DIR" ]; then
    echo " Data Directory:   $CUSTOM_DATA_DIR"
fi
echo " Output Directory: $RUN_DIR"
echo " Logs Directory:   $LOGS_DIR"
echo "=========================================================="

# ==============================================================================
# Step 1: Run Fit Jobs (ONLY if not in monitor_only, collect_only, or plot_only mode)
# ==============================================================================
if [ "$MODE" != "monitor_only" ] && [ "$MODE" != "no_fit" ] && [ "$MODE" != "collect_only" ] && [ "$MODE" != "plot_only" ] && [ "$MODE" != "calculate_acceptance" ]; then
    mkdir -p "$LOGS_DIR"

    # Safety check: Detect if fit logs or results already exist for target t-bins
    EXISTING_LOGS=()
    for i in "${TBIN_LIST[@]}"; do
        LOG_FILE="${LOGS_DIR}MomentFit_${MASS_TAG}_tbin_${i}.log"
        if [ ! -f "$LOG_FILE" ]; then
            ALT_LOG=$(ls -1 "${LOGS_DIR}"MomentFit*_tbin_"${i}".log 2>/dev/null | head -n 1 || true)
            [ -n "$ALT_LOG" ] && LOG_FILE="$ALT_LOG"
        fi
        tbin_name="${TBIN_INDEX_MAP[$i]:-tbin_${i}}"
        if [ -s "$LOG_FILE" ]; then
            EXISTING_LOGS+=("$LOG_FILE")
        elif [ -d "${RUN_DIR}${tbin_name}" ] && [ "$(ls -A "${RUN_DIR}${tbin_name}" 2>/dev/null)" ]; then
            EXISTING_LOGS+=("${RUN_DIR}${tbin_name}")
        fi
    done

    if [ "$OVERWRITE" != true ] && [ ${#EXISTING_LOGS[@]} -gt 0 ]; then
        echo ""
        echo "⚠️  Existing fit log(s) and results were detected in:"
        echo "    $RUN_DIR"
        for exf in "${EXISTING_LOGS[@]}"; do
            if [ -f "$exf" ]; then
                echo "    - $(basename "$exf") ($(du -h "$exf" 2>/dev/null | cut -f1))"
            else
                echo "    - $(basename "$exf")/ (result directory exists)"
            fi
        done
        echo ""

        CONFIRM_OVERWRITE="n"
        if [ -t 0 ]; then
            read -r -p "Do you want to overwrite the fit results? [y/N]: " CONFIRM_OVERWRITE
        else
            echo "Non-interactive session detected and --overwrite not specified."
            echo "Aborting fit launch to protect existing results."
            exit 1
        fi

        case "$CONFIRM_OVERWRITE" in
            [yY]|[yY][eE][sS])
                echo "Backing up existing log(s) and proceeding with overwrite..."
                TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
                for exf in "${EXISTING_LOGS[@]}"; do
                    if [ -f "$exf" ]; then
                        cp "$exf" "${exf}.bak_${TIMESTAMP}" 2>/dev/null || true
                    fi
                done
                ;;
            *)
                echo ""
                echo "Aborting fit execution to protect existing logs and results."
                echo "💡 Helpful commands:"
                echo "  • To draw plots without re-running fits: $0 --run-dir \"$RUN_DIR\" --force-plot"
                echo "  • To inspect progress table:             $0 --run-dir \"$RUN_DIR\" --monitor-only"
                echo "  • To force overwrite without prompt:     $0 --run-dir \"$RUN_DIR\" --overwrite"
                exit 0
                ;;
        esac
    fi

    PIDS=()

    cleanup_on_interrupt() {
        echo ""
        echo "⚠️  Received interrupt signal (Ctrl+C). Terminating running fit jobs..."
        for i in "${TBIN_LIST[@]}"; do
            kill_tbin_fit "$i"
        done
        exit 130
    }
    trap cleanup_on_interrupt SIGINT SIGTERM

    export_mom2amp_env() {
        export MOM2AMP_MIN_ACCEPTANCE="$MIN_ACCEPTANCE"
        export MOM2AMP_MAX_ACCEPTANCE="$MAX_ACCEPTANCE"
        export MOM2AMP_TARGET_ACCEPTANCE="$TARGET_ACCEPTANCE"
        export MOM2AMP_TUNE_MODE="$TUNE_MODE"
        export MOM2AMP_RHAT_WINDOW="$RHAT_WINDOW"
        export MOM2AMP_RHAT_TARGET="$RHAT_TARGET"
        export MOM2AMP_ESS_FRACTION="$ESS_FRACTION"
        export MOM2AMP_RHAT_MAX_RETRIES="$RHAT_MAX_RETRIES"
        export MOM2AMP_BAILOUT_THRESHOLD="$BAILOUT_THRESHOLD"
        export MOM2AMP_BAILOUT_MIN_ACC="$BAILOUT_MIN_ACC"
        export MOM2AMP_MC_INT_BINS="$MC_INT_BINS"
        export MOM2AMP_FIXED_COV="$FIXED_COV"
        export MOM2AMP_USE_COV="$USE_COV"
        export MOM2AMP_SET_INIT_VALUE="$SET_INIT_VALUE"
        export MOM2AMP_USE_GJ_FRAME="$USE_GJ_FRAME"
        export MOM2AMP_NUM_CPU="$NUM_CPU"
        export MOM2AMP_SEED="$MCMC_SEED"
    }

    run_tbin_worker() {
        local i="$1"
        local tbin_name="${TBIN_INDEX_MAP[$i]:-tbin_${i}}"
        local log_file="${LOGS_DIR}MomentFit_${MASS_TAG}_tbin_${i}.log"
        local pid_file="${LOGS_DIR}tbin_${i}.pid"

        local data_f="$CUSTOM_DATA_FILE"
        if [ -n "$data_f" ]; then
            data_f="${data_f//\{tbin\}/$tbin_name}"
            data_f="${data_f//\{data_dir\}/$CUSTOM_DATA_DIR}"
            data_f="${data_f//\{base_dir\}/$CUSTOM_DATA_DIR}"
        fi

        local mcrecon_f="$CUSTOM_MCRECON_FILE"
        if [ -n "$mcrecon_f" ]; then
            mcrecon_f="${mcrecon_f//\{tbin\}/$tbin_name}"
            mcrecon_f="${mcrecon_f//\{data_dir\}/$CUSTOM_DATA_DIR}"
            mcrecon_f="${mcrecon_f//\{base_dir\}/$CUSTOM_DATA_DIR}"
        fi

        local mcgen_f="$CUSTOM_MCGEN_FILE"
        if [ -n "$mcgen_f" ]; then
            mcgen_f="${mcgen_f//\{tbin\}/$tbin_name}"
            mcgen_f="${mcgen_f//\{data_dir\}/$CUSTOM_DATA_DIR}"
            mcgen_f="${mcgen_f//\{base_dir\}/$CUSTOM_DATA_DIR}"
        fi

        local weight_f="$CUSTOM_WEIGHT_FILE"
        if [ -n "$weight_f" ]; then
            weight_f="${weight_f//\{tbin\}/$tbin_name}"
            weight_f="${weight_f//\{data_dir\}/$CUSTOM_DATA_DIR}"
            weight_f="${weight_f//\{base_dir\}/$CUSTOM_DATA_DIR}"
        fi

        export_mom2amp_env
        root -l -b -q "${BRUFIT}/macros/LoadBru.C" \
            "${MACRO_PATH}($i, $MAX_L, $INIT_SCALE, $STEP, -1, $MCMC_ITERS, $MCMC_BURNIN, $M_MIN, $M_MAX, $N_BINS, \"$RUN_DIR\", \"$CUSTOM_COV_DIR\", \"$CUSTOM_DATA_DIR\", $APPLY_MC_VAR, $ALLOW_NEG_PDF, $USE_TUNING, $GIBBS_BLOCK_SIZE, $COLD_WARMUP_ACC, \"$data_f\", \"$mcrecon_f\", \"$mcgen_f\", \"$weight_f\")" > "$log_file" 2>&1
        local fit_status=$?
        rm -f "$pid_file"

        if [ $fit_status -eq 0 ]; then
            echo "  [$tbin_name] ✅ Fit completed! Running post-fit monitoring and plot generation in parallel..."
        else
            echo "  [$tbin_name] ⚠️  Fit exited with code $fit_status. Running monitoring..."
        fi

        # Immediately monitor and generate diagnostic plots for this tbin as soon as its fit finishes
        local WORKER_MON_ARGS=(
            --run-dir "$RUN_DIR"
            --tbin "$i"
            --max-l "$MAX_L"
            --m-min "$M_MIN"
            --m-max "$M_MAX"
            --n-bins "$N_BINS"
            --ess-threshold "$ESS_THRESHOLD"
            --ess-metric "$ESS_METRIC"
            --burn-in "$MCMC_BURNIN"
        )
        [ -n "$CUSTOM_DATA_DIR" ] && WORKER_MON_ARGS+=(--data-dir "$CUSTOM_DATA_DIR")
        MOM2AMP_CALCULATE_ACCEPTANCE="$ACCEPTANCE_MACRO" "$PYTHON_BIN" "$MONITOR_SCRIPT" "${WORKER_MON_ARGS[@]}" "${EXTRA_MONITOR_ARGS[@]}" > /dev/null 2>&1 || true

        echo "  [$tbin_name] 🏁 Finished pipeline for $tbin_name."
        return $fit_status
    }

    for i in "${TBIN_LIST[@]}"; do
        tbin_name="${TBIN_INDEX_MAP[$i]:-tbin_${i}}"
        LOG_FILE="${LOGS_DIR}MomentFit_${MASS_TAG}_tbin_${i}.log"
        PID_FILE="${LOGS_DIR}tbin_${i}.pid"
        echo "Launching parallel MomentFit pipeline for tbin_index = $i ($tbin_name, MaxL = $MAX_L, Scale = $INIT_SCALE, Step = $STEP)"
        echo "  -> Log: $LOG_FILE"
        run_tbin_worker "$i" &
        PID=$!
        PIDS+=($PID)
        echo "$PID" > "$PID_FILE"
    done

    echo ""
    echo "Parallel t-bin fit pipelines running in background (PIDs: ${PIDS[*]})."
    echo "To monitor real-time fit logs in another terminal:"
    echo "  tail -f ${LOGS_DIR}MomentFit_${MASS_TAG}_tbin_*.log"
    echo "To inspect structured progress table anytime:"
    echo "  $0 --run-dir \"$RUN_DIR\" --monitor-only"
    echo "To cancel running job(s) anytime:"
    echo "  $0 --run-dir \"$RUN_DIR\" --cancel"

    if [ "$NO_WAIT" = true ]; then
        echo ""
        echo "🚀 Fits launched in background (--no-wait enabled). Exiting immediately."
        exit 0
    fi

    echo ""
    echo "Waiting for all parallel t-bin pipelines to finish..."

    FAILED_JOBS=0
    for pid in "${PIDS[@]}"; do
        if ! wait "$pid"; then
            FAILED_JOBS=$((FAILED_JOBS + 1))
        fi
    done

    for i in "${TBIN_LIST[@]}"; do
        rm -f "${LOGS_DIR}tbin_${i}.pid"
    done

    echo ""
    if [ "$FAILED_JOBS" -gt 0 ]; then
        echo "⚠️  $FAILED_JOBS fit job(s) encountered non-zero return code. Inspect logs in: $LOGS_DIR"
    else
        echo "✅ All parallel t-bin MomentFit pipelines completed!"
    fi
else
    if [ "$MODE" == "plot_only" ]; then
        echo "Skipping fit execution (MODE=plot_only). Generating diagnostic plots directly."
    else
        echo "Skipping fit execution (MODE=$MODE). Running monitoring and downstream steps directly."
    fi
fi

# ==============================================================================
# Step 2: Run Monitor & Progress Pipeline
# ==============================================================================
run_monitoring() {
    echo ""
    echo "=========================================================="
    echo " Running Post-Fit Monitoring & Progress Pipeline"
    echo "=========================================================="

    local MON_ARGS=(
        --run-dir "$RUN_DIR"
        --tbin "$TBIN_ARG"
        --max-l "$MAX_L"
        --m-min "$M_MIN"
        --m-max "$M_MAX"
        --n-bins "$N_BINS"
        --ess-threshold "$ESS_THRESHOLD"
        --ess-metric "$ESS_METRIC"
        --burn-in "$MCMC_BURNIN"
    )
    [ -n "$CUSTOM_DATA_DIR" ] && MON_ARGS+=(--data-dir "$CUSTOM_DATA_DIR")

    if [ "$MODE" = "calculate_acceptance" ]; then
        if [ -n "$ACCEPTANCE_MACRO" ]; then
            MON_ARGS+=(--calculate-acceptance "$ACCEPTANCE_MACRO")
        else
            MON_ARGS+=(--calculate-acceptance)
        fi
        "$PYTHON_BIN" "$MONITOR_SCRIPT" "${MON_ARGS[@]}" "${EXTRA_MONITOR_ARGS[@]}"
        return $?
    fi

    if [ "$WATCH" = true ]; then
        MON_ARGS+=(--watch --interval "$INTERVAL")
    fi

    MOM2AMP_CALCULATE_ACCEPTANCE="$ACCEPTANCE_MACRO" "$PYTHON_BIN" "$MONITOR_SCRIPT" "${MON_ARGS[@]}" "${EXTRA_MONITOR_ARGS[@]}" || true
}

# ==============================================================================
# Step 3: Diagnostic Collector & Next-Fit Generation
# ==============================================================================
run_collection_and_retry() {
    echo ""
    echo "=========================================================="
    echo " Running Diagnostic Collection & Next-Fit Generation"
    echo "=========================================================="

    local COLLECTOR_CMD=(
        "$PYTHON_BIN" "$COLLECTOR_SCRIPT"
        --run-dir "$RUN_DIR"
        --tbin "$TBIN_ARG"
        --m-min "$M_MIN"
        --m-max "$M_MAX"
        --n-bins "$N_BINS"
        --ess-threshold "$ESS_THRESHOLD"
        --ess-metric "$ESS_METRIC"
        --min-ess-threshold "$MIN_ESS_THRESHOLD"
        --burn-in "$MCMC_BURNIN"
        --max-parallel "$MAX_PARALLEL"
        --macro-path "$MACRO_PATH"
        --gibbs-block-size "$GIBBS_BLOCK_SIZE"
        --cold-warmup-acc "$COLD_WARMUP_ACC"
        --min-acceptance "$MIN_ACCEPTANCE"
        --max-acceptance "$MAX_ACCEPTANCE"
        --tune-mode "$TUNE_MODE"
        --rhat-window "$RHAT_WINDOW"
        --rhat-target "$RHAT_TARGET"
        --ess-fraction "$ESS_FRACTION"
        --rhat-max-retries "$RHAT_MAX_RETRIES"
        --bailout-threshold "$BAILOUT_THRESHOLD"
        --bailout-min-acceptance "$BAILOUT_MIN_ACC"
        --mc-integration-bins "$MC_INT_BINS"
        --num-cpu "$NUM_CPU"
        --seed "$MCMC_SEED"
        --max-retry-iters "$MAX_RETRY_ITERS"
        --burnin-fraction "$BURNIN_FRACTION"
        --low-acceptance-threshold "$LOW_ACC_THRESH"
        --gibbs-low-acceptance-threshold "$GIBBS_LOW_ACC_THRESH"
        --high-acceptance-threshold "$HIGH_ACC_THRESH"
        --gibbs-activation-acceptance "$GIBBS_ACT_ACC"
        --recondition-shrinkage "$RECONDITION_SHRINKAGE"
        --recondition-min-eigen-floor "$RECONDITION_MIN_EIGEN"
    )

    if [ -n "$TARGET_ACCEPTANCE" ]; then COLLECTOR_CMD+=(--target-acceptance "$TARGET_ACCEPTANCE"); fi
    if [ -n "$CUSTOM_COV_DIR" ]; then COLLECTOR_CMD+=(--cov-dir "$CUSTOM_COV_DIR"); fi
    if [ -n "$CUSTOM_DATA_DIR" ]; then COLLECTOR_CMD+=(--data-dir "$CUSTOM_DATA_DIR"); fi
    if [ -n "$CUSTOM_DATA_FILE" ]; then COLLECTOR_CMD+=(--data-pattern "$CUSTOM_DATA_FILE"); fi
    if [ -n "$CUSTOM_MCRECON_FILE" ]; then COLLECTOR_CMD+=(--mcrecon-pattern "$CUSTOM_MCRECON_FILE"); fi
    if [ -n "$CUSTOM_MCGEN_FILE" ]; then COLLECTOR_CMD+=(--mcgen-pattern "$CUSTOM_MCGEN_FILE"); fi
    if [ -n "$CUSTOM_WEIGHT_FILE" ]; then COLLECTOR_CMD+=(--weight-pattern "$CUSTOM_WEIGHT_FILE"); fi
    if [ "$APPLY_MC_VAR" = true ]; then COLLECTOR_CMD+=(--mc-variance); else COLLECTOR_CMD+=(--no-mc-variance); fi
    if [ "$ALLOW_NEG_PDF" = true ]; then COLLECTOR_CMD+=(--negative-pdf); else COLLECTOR_CMD+=(--no-negative-pdf); fi
    if [ "$USE_TUNING" = true ]; then COLLECTOR_CMD+=(--tuning); else COLLECTOR_CMD+=(--no-tuning); fi
    if [ "$FIXED_COV" = true ]; then COLLECTOR_CMD+=(--fixed-covariance); else COLLECTOR_CMD+=(--no-fixed-covariance); fi
    if [ "$USE_COV" = true ]; then COLLECTOR_CMD+=(--use-covariance); else COLLECTOR_CMD+=(--no-use-covariance); fi
    if [ "$SET_INIT_VALUE" = true ]; then COLLECTOR_CMD+=(--set-init-value); else COLLECTOR_CMD+=(--no-set-init-value); fi
    if [ "$USE_GJ_FRAME" = true ]; then COLLECTOR_CMD+=(--gj-frame); else COLLECTOR_CMD+=(--hel-frame); fi
    if [ "$RECONDITION_COV" = true ]; then COLLECTOR_CMD+=(--recondition-covariance); else COLLECTOR_CMD+=(--no-recondition-covariance); fi

    if [ "$EXPLICIT_MAX_L" = true ]; then COLLECTOR_CMD+=(--max-l "$MAX_L"); fi
    if [ "$EXPLICIT_SCALE" = true ]; then COLLECTOR_CMD+=(--init-scale "$INIT_SCALE"); fi
    if [ "$EXPLICIT_STEP" = true ]; then COLLECTOR_CMD+=(--step "$STEP"); fi
    COLLECTOR_CMD+=(--iters "$MCMC_ITERS")

    # Pass --no-exec so RunMomentFit.sh orchestrates the retry execution loop
    COLLECTOR_CMD+=(--no-exec)

    if [ "$AUTO_RETRY" = true ]; then
        COLLECTOR_CMD+=(--non-interactive)
    fi

    "${COLLECTOR_CMD[@]}"
}

# Run monitoring step unless in collect_only mode
if [ "$MODE" != "collect_only" ] && [ "$MODE" != "fit_only" ]; then
    run_monitoring
fi

# Run collector and optional auto-retry loop
if [ "$MODE" != "fit_only" ] && [ "$MODE" != "monitor_only" ] && [ "$MODE" != "plot_only" ] && [ "$MODE" != "calculate_acceptance" ] && [ "$WATCH" != true ]; then
    CURRENT_ROUND=1
    while true; do
        run_collection_and_retry

        RETRY_SCRIPT="${RUN_DIR}RunMomentFit_retry.sh"
        SUMMARY_JSON="${RUN_DIR}unsuccessful_bins_summary.json"

        if [ ! -f "$SUMMARY_JSON" ]; then
            break
        fi

        UNSUCC_COUNT=$(python3 -c "import json; d=json.load(open('$SUMMARY_JSON')); print(d.get('unsuccessful_count', 0))" 2>/dev/null || echo 0)

        if [ "$UNSUCC_COUNT" -eq 0 ]; then
            echo ""
            echo "🎉 ALL mass bins across all t-bins have achieved $ESS_METRIC ESS >= $ESS_THRESHOLD!"
            echo "   Convergence criterion fully satisfied."
            break
        fi

        if [ "$AUTO_RETRY" = true ]; then
            if [ "$CURRENT_ROUND" -le "$MAX_RETRIES" ] && [ -f "$RETRY_SCRIPT" ]; then
                echo ""
                echo "=========================================================="
                echo " [Auto-Retry Round $CURRENT_ROUND / $MAX_RETRIES]"
                echo " Automatically executing targeted refits for $UNSUCC_COUNT bin(s)..."
                echo "=========================================================="
                bash "$RETRY_SCRIPT" --non-interactive --no-exec

                ((CURRENT_ROUND++))
            else
                echo ""
                echo "⚠️  Maximum auto-retry limit ($MAX_RETRIES) reached or retry script not found."
                echo "    Remaining unsuccessful bins: $UNSUCC_COUNT. Check diagnostics in: $SUMMARY_JSON"
                break
            fi
        else
            # Check if running interactively on a terminal
            if [ -t 0 ] && [ -f "$RETRY_SCRIPT" ]; then
                echo ""
                echo "----------------------------------------------------------"
                read -r -p "Would you like to execute the retry script now for these $UNSUCC_COUNT unsuccessful bin(s)? [Y/n]: " USER_PROMPT_RETRY
                case "$USER_PROMPT_RETRY" in
                    [nN][oO]|[nN])
                        echo "Skipping retry execution. You can manually launch it later:"
                        echo "  bash $RETRY_SCRIPT"
                        break
                        ;;
                    *)
                        echo "=========================================================="
                        echo " [Interactive Retry Round $CURRENT_ROUND]"
                        echo " Executing targeted refits for $UNSUCC_COUNT bin(s)..."
                        echo "=========================================================="
                        bash "$RETRY_SCRIPT" --non-interactive --no-exec
                        ((CURRENT_ROUND++))
                        ;;
                esac
            else
                break
            fi
        fi
    done
fi

echo ""
echo "=========================================================="
echo " Pipeline Complete for: $RUN_DIR"
echo "=========================================================="
