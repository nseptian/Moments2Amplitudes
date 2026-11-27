# Moments2Amplitudes

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/nseptian/Moments2Amplitudes)

Tools for **mass-dependent modeling of angular moments of two-pseudoscalar mesons produced in polarized photon-proton collisions** (such as $\eta\pi^0$). 

The main goal is to extract partial wave amplitudes by fitting mass-dependent, physically-motivated parameterizations directly to the estimated physical moments and covariance matrices produced by the [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) MCMC sampler — going from **Moments to Amplitudes**.

---

## Physics Overview

The angular distribution of two pseudoscalar mesons produced in photoproduction can be expanded in terms of angular moments $H(L, M)$, which are bilinear sums of the underlying partial wave amplitudes $T_{J, m}^\epsilon$ (spin $J$, projection $m$, reflectivity $\epsilon = \pm 1$). 

Rather than estimating amplitudes independently bin-by-bin in invariant mass, this framework imposes mass-dependent parameterizations across the entire mass spectrum simultaneously. It reads the moment values and covariance matrices from brufit MCMC fit outputs (`ResultsBruMcmcCovariance.root`) and performs a global $\chi^2$ minimization.

### Amplitude Parameterizations

Supported parameterizations in `MassDependentFunction`:

| Model | Description |
|---|---|
| Breit-Wigner | Resonances such as $a_2(1320)$ ($D$-wave) |
| Flatté | Coupled-channel threshold effects, e.g., $a_0(980) \to \eta\pi, K\bar{K}$ ($S$-wave) |
| Breit-Wigner + Polynomial | Resonance on top of a smooth background |
| Flatté + Polynomial | Flatté amplitude with polynomial background term |
| Two Breit-Wigners + Polynomial | Coherent sum of two Breit-Wigner states plus background |
| Flatté + Breit-Wigner + Polynomial | Combined $S$-wave model |

---

## Configuration-Driven Framework & Pixi Integration

Moments2Amplitudes provides a human-readable, declarative YAML configuration system and a command-line interface (`mom2amp`) integrated with **Pixi**:

```bash
# 1. Install dependencies with Pixi
pixi install

# 2. Configure BruFit and ROOT (copy .env.example -> .env or set $BRUFIT)
cp .env.example .env

# 3. Validate a YAML fit configuration
pixi run validate configs/d_breit_wigner.yaml

# 4. Run multi-start randomized fit in parallel (e.g., 20 starts on 8 cores)
pixi run fit configs/d_breit_wigner.yaml

# 5. Display formatted physics parameter summary table and inspect results
pixi run summary configs/d_breit_wigner.yaml

# 6. Generate diagnostic distribution, amplitude plots, and moment comparisons with MCMC
pixi run plot configs/d_breit_wigner.yaml
# Optional: overlay N lowest chi2 fit curves and filter specific moments
pixi run plot configs/s_piecewise_d_bw.yaml --n-lowest-chi2 10 --moments L=2

# 7. Calculate reflectivity-separated resonance cross sections
pixi run cs configs/s_piecewise_d_bw.yaml
# Multi-t bin evaluation:
# pixi run cs configs/s_piecewise_d_bw.yaml --runs logs_t010020 logs_t020032 logs_t032050 logs_t050075

# 8. Propagate MCMC experimental uncertainty across posterior samples
pixi run propagate configs/s_piecewise_d_bw.yaml
```

See **[`docs/CONFIG_FRAMEWORK.md`](docs/CONFIG_FRAMEWORK.md)** for complete configuration options, two-stage fit workflows, and schema specifications, and **[`docs/MOMENTFIT.md`](docs/MOMENTFIT.md)** for the end-to-end linear algebra -> SPD covariance export -> BruFit MCMC (`pixi run momentfit`) workflow.

---

- [ROOT](https://root.cern/) (v6.28+, with Minuit2 and ACLiC enabled)
- [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) (PWA-moment relations and MCMC fitting)
- [Moments](https://github.com/bgrube/Moments) (linear algebra method for calculating angular moments)

---

## Repository Structure

```
Moments2Amplitudes/
├── configs/                          # Declarative YAML fit configuration files
│   ├── d_breit_wigner.yaml           # Stage 1: D-wave Breit-Wigner fit to H4 moments
│   ├── s_piecewise_d_bw.yaml         # Stage 2: S-wave piecewise + D-wave Breit-Wigner to H0, H2, H4
│   └── ...
│
├── include/                          # C++ headers for mass-dependent fitting
│   ├── MassDependentFitter.h         # Main coordinator and Minuit2 fitting interface
│   ├── Chi2Engine.h                  # High-performance chi2 computation and covariance inversion
│   ├── ParameterManager.h            # Parameter mapping, limits, and wave initialization rules
│   ├── MCMCFitter.h                  # Multi-threaded MCMC posterior uncertainty propagation
│   ├── MassDependentFunction.h       # Parameterized line shapes (Breit-Wigner, Flatté, polynomials)
│   ├── MassDependentMoments.h        # Loads moment results and covariances across mass bins
│   ├── MomentHelper.h                # Per-bin moment data and covariance matrix handler
│   ├── ParameterHelper.h             # Parameter naming and mapping utilities
│   ├── Equation.h / EquationSolver.h # Bilinear equations relating moments to amplitudes
│
├── src/                              # C++ implementation files (compiled via ROOT ACLiC)
│   ├── MassDependentFitter.cpp
│   ├── Chi2Engine.cpp
│   ├── ParameterManager.cpp
│   ├── MCMCFitter.cpp
│   ├── ParameterHelper.cpp
│   ├── Equation.cpp / EquationSolver.cpp
│   └── MomentHelper.cpp
│
├── macros/                           # Fit execution and helper macros
│   ├── Load.C                        # Loads required libraries, headers, and ACLiC modules
│   ├── ConfigureAmpsNoValues.C       # Sets up amplitude wave definitions
│   └── ...
│
├── python/                           # Python analysis tools, driver, CLI, and post-processing
│   ├── moments2amplitudes/           # Modular package
│   │   ├── cli.py                    # Unified CLI (`mom2amp validate/fit/summary/plot/propagate/momentfit`)
│   │   ├── driver.py                 # Macro generation, worker dispatching, and run orchestration
│   │   ├── momentfit/                # End-to-end MomentFit workflow (linear algebra, cov export, BruFit MCMC, diagnostics)
│   │   ├── summary.py                # Formatted ASCII and Markdown physics summary tables
│   │   ├── plotter.py                # Chi2 distribution and diagnostic plotting
│   │   ├── config.py                 # Configuration loader and parser
│   │   ├── schema.py                 # Pydantic schema validation
│   │   └── rundir.py                 # Standardized run directory hierarchy manager
│   ├── BruFitResults.py              # BruFit fit-status verification and bounded-memory MCMC chain diagnostics
│   ├── moments_to_amplitudes.py      # Core Moments2Amplitudes analysis class
│   └── MomentsToAmplitudes.py        # Module alias / backwards-compatibility wrapper
│
└── set_env.sh                        # Environment setup script
```

---

## Setup

### 1. Build brufit

Follow the [brufit build instructions](https://github.com/dglazier/brufit/tree/R6.34Test). Make sure the BruFit environment script is sourced before running any macros.

### 2. Clone and set up Moments

```bash
git clone https://github.com/bgrube/Moments.git
cd Moments
pip install -r requirements.txt
```

> If running on a cluster other than ifarm, edit `OpenMpUtilities.py` line 47 and replace `"ifarm"` with your hostname.

### 3. Set environment variables

Edit [`set_env.sh`](set_env.sh) to point to your local installations:

```bash
export Moments=/path/to/Moments/
export MOMS2AMPS=/path/to/Moments2Amplitudes/
export ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${ROOT_INCLUDE_PATH}"
```

Then source it:

```bash
source set_env.sh
```

---

## Running Example Fits

The scripts in `macros/` provide test and example setups for mass-dependent fits using moments obtained from brufit MCMC sampling.

### Single Fit Execution

To run a single fit macro with ROOT:

```bash
root -l -b -q ${BRUFIT}/macros/LoadBru.C macros/Load.C 'macros/Mom2Amp_SFlatteDBreitWigner.C(0)'
```

Here `0` is the random seed passed to initialize free parameters.

### Multi-Seed Parallel Execution

To explore parameter space and check fit stability across different initial seeds, use `run_model.sh`:

```bash
./macros/run_model.sh -n fit_run_01 -s 1 -e 100 -j 40 Mom2Amp_SFlatteDBreitWigner.C
```

Options:
- `-n, --name`: Name of the run (creates a timestamped directory under `run/`).
- `-s, --start`: Starting seed number (default: `1`).
- `-e, --end`: Ending seed number (default: `100`).
- `-j, --jobs`: Maximum concurrent background jobs (default: `80`).

Each run creates an organized output directory containing:
```
run/fit_run_01_YYYYMMDD_HHMMSS/
├── macro/       # Snapshot of macro files used for the run
├── logs/        # Individual seed logs
├── results/     # Output ROOT trees and fit parameters
└── metadata.txt # Run environment and summary
```

### Multi-Start Global Minimum Search & MCMC Uncertainty Propagation

To avoid local minima traps and propagate experimental uncertainties cleanly, use the randomized multi-start and parallel MCMC propagation workflows:

1. **Pure D-wave fit to $H_4$ ($a_2(1320)$)**:
   ```bash
   brufit -l -b -q macros/Load.C 'macros/Mom2Amp_DBreitWigner_MCMC.C(100, 0, 200, 1, "mcmc_propagation_results.root", "logs", 40)'
   ```

2. **S-wave Piecewise + D-wave Breit-Wigner to $H_0, H_2, H_4$**:
   ```bash
   brufit -l -b -q macros/Load.C 'macros/Mom2Amp_SPiecewiseDBreitWigner.C(100, 0, 200, 1, "logs/central_best_params.root", "mcmc_propagation_spiecewise_results.root", "logs_spiecewise", 40)'
   ```

See [`docs/MCMC_PROPAGATION.md`](docs/MCMC_PROPAGATION.md) for full documentation of options, initialization rules, and output tree schemas.

### Configuring Amplitude Models

In the macros, model parameterizations are assigned per partial wave using `MassDependenceConfig`:

```cpp
using MDModel = m2pw::MassDependentFitter::MassDependenceConfig::ModelType;
using WaveModelConfig = m2pw::MassDependentFitter::MassDependenceConfig::WaveModelConfig;

// S-wave (L=0) parameterized with Flatté + polynomial; D-wave (L=2) with Breit-Wigner; P-wave (L=1) mass-independent
auto config = m2pw::MassDependentFitter::CreateCustomConfig(
    {
        {0, {WaveModelConfig{"a0_980", MDModel::FlattePlusPolynomial, 2}}},
        {2, {WaveModelConfig{"a2_1320", MDModel::BreitWigner}}}
    },
    {1} // Mass-independent waves
);
```

---

## Python Analysis Tools (`Moments2Amplitudes`)

The `python/` directory provides Python utilities for post-processing and visualizing fit results from `MassDependentFitter::MakeResultTree`:

```python
from moments_to_amplitudes import Moments2Amplitudes

# Initialize analyzer with fit results directory
analyzer = Moments2Amplitudes(
    fit_result_dir="path/to/fit/results",
    output_dir="./plots",
    filename_pattern="result_tree_[]_nominal.root",
    n_mass_bins=30,
    first_mass_bin_center=0.82,
    mass_bin_width=0.04
)

# 1. Load multi-seed fit files and filter valid minima (isValid == 1)
analyzer.load_fit_files(max_files=100)

# 2. Inspect chi2 distribution and rank best fits
analyzer.print_chi2_values()
analyzer.plot_chi2_distribution(output_filename="chi2_distribution.pdf")

# 3. Load the best N fits
analyzer.load_trees_from_best_fits(n_best=10)

# 4. Extract mass-dependent resonance parameters (e.g. mass, width)
params = analyzer.get_mass_dependent_parameters(
    parameters=["posRefl_a2_1320_mass", "posRefl_a2_1320_width"],
    lowest_chi2_tree_index=0
)

# 5. Plot fitted amplitudes vs invariant mass overlaid with H0(0,0)
analyzer.plot_amplitudes_vs_mass(
    unnormalized_moments_file="path/to/UnnormalizedMoments_vs_Mass.root",
    output_filename="Mom2Amps_a2_vs_mass.pdf",
    wave_type="a_2"
)
```

---

## Key Classes

| Class | Header | Description |
|---|---|---|
| `m2pw::MassDependentFitter` | [`include/MassDependentFitter.h`](include/MassDependentFitter.h) | High-level coordinator, wave configuration, and Minuit2 interface |
| `m2pw::Chi2Engine` | [`include/Chi2Engine.h`](include/Chi2Engine.h) | High-performance $\chi^2$ calculation, pre-compiled routing, covariance matrix inversion, and singular fallback |
| `m2pw::ParameterManager` | [`include/ParameterManager.h`](include/ParameterManager.h) | Parameter lifecycle, naming conventions, bounds, and wave initialization rules |
| `m2pw::MCMCFitter` | [`include/MCMCFitter.h`](include/MCMCFitter.h) | Multi-threaded MCMC posterior uncertainty propagation worker management |
| `m2pw::MassDependentFunction` | [`include/MassDependentFunction.h`](include/MassDependentFunction.h) | Computes complex partial wave amplitudes for Breit-Wigner, Flatté, and polynomial models |
| `MassDependentMoments` | [`include/MassDependentMoments.h`](include/MassDependentMoments.h) | Loads brufit MCMC moment values and covariance matrices across mass bins |
| `MomentHelper` | [`include/MomentHelper.h`](include/MomentHelper.h) | Interface for extracting moments and covariance matrices per mass bin from MCMC trees |
| `Equation` / `EquationSolver` | [`include/Equation.h`](include/Equation.h) | Encodes symbolic/numerical relationships between angular moments and partial wave bilinear products |

---

## References

- [brufit](https://github.com/dglazier/brufit) — ROOT-based framework for partial wave analysis and moment estimation
- B. Grube, [Moments framework](https://github.com/bgrube/Moments) — Linear algebra method for calculating angular moments from data
