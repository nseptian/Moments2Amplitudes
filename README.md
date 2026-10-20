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

## Requirements

- [ROOT](https://root.cern/) (v6.28+, with Minuit2 and ACLiC enabled)
- [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) (PWA-moment relations and MCMC fitting)
- [Moments](https://github.com/bgrube/Moments) (linear algebra method for calculating angular moments)

---

## Repository Structure

```
Moments2Amplitudes/
├── include/                          # C++ headers for mass-dependent fitting
│   ├── MassDependentFitter.h         # Core fitter: chi2 calculation, Minuit2 interface, parameter management
│   ├── MassDependentFunction.h       # Parameterized line shapes (Breit-Wigner, Flatté, polynomials)
│   ├── MassDependentMoments.h        # Loads moment results and covariances across mass bins
│   ├── MomentHelper.h                # Per-bin moment data and covariance matrix handler
│   ├── ParameterHelper.h             # Parameter naming and mapping utilities
│   ├── Equation.h / EquationSolver.h # Bilinear equations relating moments to amplitudes
│
├── src/                              # C++ implementation files (compiled via ROOT ACLiC)
│   ├── MassDependentFitter.cpp
│   ├── ParameterHelper.cpp
│   ├── Equation.cpp / EquationSolver.cpp
│   └── MomentHelper.cpp
│
├── macros/                           # Example fit macros and execution scripts
│   ├── Load.C                        # Loads required libraries and includes
│   ├── ConfigureAmpsNoValues.C       # Sets up amplitude wave definitions
│   ├── run_model.sh                  # Shell script for parallel multi-seed fits with organized output
│   ├── H_two_model.C                 # Example: S + D waves fitting H_0 and H_2 moments
│   ├── H_four_model.C                # Example: waves up to L=4 (H_0, H_2, H_4)
│   ├── H_odd_model.C                 # Example: odd moments (P-wave contribution)
│   ├── Mom2Amp_DBreitWigner.C        # Example: D-wave Breit-Wigner fit
│   ├── Mom2Amp_SFlatteDBreitWigner.C # Example: S-wave Flatté + D-wave Breit-Wigner
│   ├── Mom2Amp_SFlattepolyDBreitWigner.C  # Example: S-wave (Flatté + poly) + D-wave (BW)
│   └── RunGivenBruMoments.C         # Example: Fit using pre-existing moments
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

## Key Classes

| Class | Header | Description |
|---|---|---|
| `m2pw::MassDependentFitter` | [`include/MassDependentFitter.h`](include/MassDependentFitter.h) | Configures waves, sets up the $\chi^2$ objective function, and executes Minuit2 minimization |
| `m2pw::MassDependentFunction` | [`include/MassDependentFunction.h`](include/MassDependentFunction.h) | Computes complex partial wave amplitudes for Breit-Wigner, Flatté, and polynomial models |
| `MassDependentMoments` | [`include/MassDependentMoments.h`](include/MassDependentMoments.h) | Loads brufit MCMC moment values and full covariance matrices across all mass bins |
| `MomentHelper` | [`include/MomentHelper.h`](include/MomentHelper.h) | Interface for extracting moments and covariance matrices per mass bin |
| `Equation` / `EquationSolver` | [`include/Equation.h`](include/Equation.h) | Encodes symbolic/numerical relationships between angular moments and partial wave bilinear products |

---

## References

- [brufit](https://github.com/dglazier/brufit) — ROOT-based framework for partial wave analysis and moment estimation
- B. Grube, [Moments framework](https://github.com/bgrube/Moments) — Linear algebra method for calculating angular moments from data
