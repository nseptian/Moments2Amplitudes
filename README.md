# Moments2Amplitudes

[![DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/nseptian/Moments2Amplitudes)

Tools for **mass-dependent partial wave analysis (PWA) of two-pseudoscalar meson photoproduction** — specifically $\eta\pi^0$ produced in polarized photon-proton collisions. The name reflects the full pipeline: starting from measured angular **Moments** (obtained via brufit MCMC), then fitting mass-dependent **Amplitudes** to the output.

Inspired by and compatible with [brufit](https://github.com/dglazier/brufit/tree/R6.34Test/tutorials/PhotoAmps/TwoSpin0AmpsFromMoments) and the [Moments](https://github.com/bgrube/Moments) framework.

---

## Physics Overview

The angular distribution of $\gamma p \to \eta\pi^0 p$ events can be decomposed into a set of angular moments $H(L, M)$ — bilinear combinations of the underlying partial wave amplitudes $T_{J,m}^{\epsilon}$ (where $J$ is spin, $m$ is helicity, $\epsilon$ is reflectivity). This project takes the posterior moment distributions from brufit's MCMC sampler and fits physically-motivated, mass-continuous amplitude models across all mass bins simultaneously to extract resonance parameters.

### Supported Amplitude Parameterizations

| Model | Use case |
|---|---|
| Breit-Wigner | Narrow resonances (e.g. $a_2(1320)$, D-wave) |
| Flatté | Coupled-channel threshold effects (e.g. $a_0(980) \to \eta\pi, K\bar{K}$) |
| Breit-Wigner + polynomial | Resonance on a smooth background |
| Flatté + polynomial | Coupled-channel + smooth background |
| Flatté + Breit-Wigner + polynomial | Full S-wave model |
| Conformal polynomial | Model-independent, analytic S-wave expansion |

The conformal polynomial maps the physical mass region $s \in [s_\text{min}, \infty)$ to the unit disk $|z| < 1$ and expands the amplitude as $T(m) = \sum_{i=0}^{N}(a_i + ib_i)\,z^i$, respecting analyticity without imposing a specific resonance hypothesis.

---

## Requirements

- [ROOT](https://root.cern/) 6.28+ (with ACLiC, RDataFrame, Minuit2)
- [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) — PWA/moment relations and MCMC fitting
- Python 3.9 with the packages listed in [`pixi.toml`](pixi.toml) (managed via [pixi](https://pixi.sh))

---

## Repository Structure

```
Moments2Amplitudes/
├── include/                          # C++ headers
│   ├── MassDependentFitter.h         # Main fitter: chi2, minimizer, parameter management
│   ├── MassDependentFunction.h       # Line-shape library (BW, Flatté, conformal poly, ...)
│   ├── MassDependentMoments.h        # Loads brufit MCMC results across mass bins
│   ├── MomentHelper.h                # Per-bin moment I/O and covariance access
│   ├── ParameterHelper.h             # Parameter naming conventions and helpers
│   ├── Equation.h / EquationSolver.h # System of equations relating moments to amplitudes
│
├── src/                              # C++ implementation files (compiled by ROOT ACLiC)
│   ├── MassDependentFitter.cpp
│   ├── ParameterHelper.cpp
│   ├── Equation.cpp / EquationSolver.cpp
│   └── MomentHelper.cpp
│
├── macros/                           # Example ROOT macros for mass-dependent fits
│   ├── Load.C                        # ROOT library loader for brufit + this project
│   ├── ConfigureAmpsNoValues.C       # Set up amplitude basis (Lmax, Mmax, Nref)
│   ├── run_model.sh                  # Parallel fit runner (multi-seed, timestamped output)
│   ├── H_two_model.C                 # Example: S+D waves (Flatté a0 + BW a2), L=0,2
│   ├── H_four_model.C                # Example: L=0,2,4 waves
│   ├── H_odd_model.C                 # Example: odd-L (P-wave) model
│   ├── Mom2Amp_DBreitWigner.C        # Example: D-wave only BW
│   ├── Mom2Amp_SFlatteDBreitWigner.C # Example: S(Flatté) + D(BW)
│   ├── Mom2Amp_SFlattepolyDBreitWigner.C   # Example: S(Flatté+poly) + D(BW)
│   ├── Mom2Amp_SFlattepolyPDBreitWigner.C  # Example: S(Flatté+poly) + P + D(BW)
│   ├── Mom2Amp_SPiecewiseDBreitWigner.C    # Example: S(piecewise) + D(BW)
│   └── RunGivenBruMoments.C         # Run fit from previously stored brufit moments
│
├── set_env.sh                        # Set MOMS2AMPS, Moments, ROOT_INCLUDE_PATH
└── pixi.toml                         # Python environment specification
```

---

## Setup

### 1. Build brufit

Follow the [brufit build instructions](https://github.com/dglazier/brufit/tree/R6.34Test). Make sure the BruFit environment script is sourced before running any macros.

### 2. Set environment variables

Edit [`set_env.sh`](set_env.sh) to point to your local installations:

```bash
export MOMS2AMPS=/path/to/Moments2Amplitudes/
export ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${ROOT_INCLUDE_PATH}"
```

Then source it:

```bash
source set_env.sh
```

---

## Mass-Dependent Amplitude Fit

After obtaining `ResultsBruMcmcCovariance.root` per mass bin from brufit's MCMC sampler, run a mass-dependent model macro. For example:

```bash
# Single seed
root -l -b -q ${BRUFIT}/macros/LoadBru.C macros/Load.C 'macros/Mom2Amp_SFlatteDBreitWigner.C(0)'

# Run many seeds in parallel with organized output
./macros/run_model.sh -n my_fit -s 1 -e 100 -j 40 Mom2Amp_SFlatteDBreitWigner.C
```

The `run_model.sh` script creates a timestamped output directory:
```
run/my_fit_YYYYMMDD_HHMMSS/
  ├── macro/      # snapshot of macros used
  ├── logs/       # per-seed logs
  ├── results/    # result ROOT trees
  └── metadata.txt
```

### Configuring the amplitude model

```cpp
using MDModel = m2pw::MassDependentFitter::MassDependenceConfig::ModelType;
using WaveModelConfig = m2pw::MassDependentFitter::MassDependenceConfig::WaveModelConfig;

auto config = m2pw::MassDependentFitter::CreateCustomConfig(
    {
        {0, {WaveModelConfig{"a0_980", MDModel::FlattePlusPolynomial, 2}}},  // S-wave
        {2, {WaveModelConfig{"a2_1320", MDModel::BreitWigner}}}              // D-wave
    },
    {1}  // P-wave mass-independent
);
```

---

## Key Classes

| Class | File | Description |
|---|---|---|
| `m2pw::MassDependentFitter` | [`include/MassDependentFitter.h`](include/MassDependentFitter.h) | Top-level fitter; sets up chi2, runs Minuit2, manages parameters |
| `m2pw::MassDependentFunction` | [`include/MassDependentFunction.h`](include/MassDependentFunction.h) | Evaluates complex amplitude for a given line-shape model |
| `MassDependentMoments` | [`include/MassDependentMoments.h`](include/MassDependentMoments.h) | Loads brufit MCMC results and covariances across mass bins |
| `MomentHelper` | [`include/MomentHelper.h`](include/MomentHelper.h) | Per-bin moment access: values, errors, covariance matrix |
| `Equation` / `EquationSolver` | [`include/Equation.h`](include/Equation.h) | Encodes the system $H_{L,M} = f(\{T_{J,m}^\epsilon\})$ |

---

## References

- B. Grube, [Moments framework](https://github.com/bgrube/Moments)
- D. Glazier, [brufit](https://github.com/dglazier/brufit)
