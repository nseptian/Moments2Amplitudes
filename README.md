# Moments2Amplitudes

[![DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/nseptian/Moments2Amplitudes)

Tools for **mass-dependent partial wave analysis (PWA) of two-pseudoscalar meson photoproduction** — specifically $\eta\pi^0$ and $\pi^+\pi^-$ produced in polarized photon-proton collisions (GlueX experiment). The framework extracts partial wave amplitudes from angular moments ($H_{L,M}$) using physically-motivated, mass-dependent parameterizations, going beyond mass-independent bin-by-bin fitting to impose analyticity and continuity across the full mass spectrum.

Inspired by and compatible with [brufit](https://github.com/dglazier/brufit/tree/R6.34Test/tutorials/PhotoAmps/TwoSpin0AmpsFromMoments) and the [Moments](https://github.com/bgrube/Moments) framework.

---

## Physics Overview

The angular distribution of $\gamma p \to \eta\pi^0 p$ events can be decomposed into a set of angular moments $H(L, M)$ — bilinear combinations of the underlying partial wave amplitudes $T_{J,m}^{\epsilon}$ (where $J$ is spin, $m$ is helicity, $\epsilon$ is reflectivity). This project addresses two complementary problems:

1. **Linear Algebra method** — direct extraction of mass-independent amplitudes from moments in each mass bin via matrix inversion, with full covariance propagation.
2. **Mass-dependent fit** — parameterize amplitudes with physical line-shapes (Breit-Wigner, Flatté, conformal polynomials) and fit them simultaneously across all mass bins by minimizing a global $\chi^2$ against the measured moments.

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
- [Moments](https://github.com/bgrube/Moments) — moment calculation and linear algebra pipeline
- Python 3.9 with the packages listed in [`pixi.toml`](pixi.toml) (managed via [pixi](https://pixi.sh))

---

## Repository Structure

```
Moments2Amplitudes/
├── include/                          # C++ headers
│   ├── MassDependentFitter.h         # Main fitter: chi2, minimizer, parameter management
│   ├── MassDependentFunction.h       # Line-shape library (BW, Flatté, conformal poly, ...)
│   ├── MassDependentMoments.h        # Loads brufit moment results across mass bins
│   ├── MomentHelper.h                # Moment I/O and covariance access
│   ├── ParameterHelper.h             # Parameter naming conventions and helpers
│   ├── Equation.h / EquationSolver.h # System of equations relating moments to amplitudes
│
├── src/                              # C++ implementation files (compiled by ROOT ACLiC)
│   ├── MassDependentFitter.cpp
│   ├── ParameterHelper.cpp
│   ├── Equation.cpp / EquationSolver.cpp
│   └── MomentHelper.cpp
│
├── macros/                           # Scripts and ROOT macros
│   ├── WeightPhSp.py                 # Generate weighted phase-space MC from amplitude model
│   ├── CalcMoments_IO_test.py        # Calculate angular moments (linear algebra method)
│   ├── PlotMoments_IO_test.py        # Plot moment results and comparisons
│   ├── ExportCovariance_IO_test.py   # Export moment covariance matrices for brufit
│   ├── ExportCovarianceC.C           # ROOT macro: write covariance to brufit format
│   ├── ExportCovariance.sh           # Shell wrapper for ExportCovarianceC
│   ├── photoProdCalcMoments_IO_test.py   # Photo-production moment calculation variant
│   ├── photoProdPlotMoments_IO_test.py   # Photo-production plot variant
│   ├── run_model.sh                  # Parallel fit runner (multi-seed, timestamped output)
│   ├── Load.C                        # ROOT library loader for brufit + this project
│   ├── ConfigureAmpsNoValues.C       # Set up amplitude basis (Lmax, Mmax, Nref)
│   ├── H_two_model.C                 # Mass-dep fit: S+D waves (Flatté a0 + BW a2)
│   ├── H_four_model.C                # Mass-dep fit: L=0,2,4 waves
│   ├── H_odd_model.C                 # Mass-dep fit: odd-L (P-wave) model
│   ├── Mom2Amp_DBreitWigner.C        # D-wave only BW model
│   ├── Mom2Amp_SFlatteDBreitWigner.C # S(Flatté) + D(BW) model
│   ├── Mom2Amp_SFlattepolyDBreitWigner.C   # S(Flatté+poly) + D(BW) model
│   ├── Mom2Amp_SFlattepolyPDBreitWigner.C  # S(Flatté+poly) + P + D(BW) model
│   ├── Mom2Amp_SPiecewiseDBreitWigner.C    # S(piecewise) + D(BW) model
│   └── RunGivenBruMoments.C         # Run fit from previously stored brufit moments
│
├── Moments.config/                   # Analysis configuration
│   └── Config_IO_test.py             # AnalysisConfig objects (mass binning, t-bins, data paths)
│
├── notebook/                         # Jupyter notebooks for analysis and visualization
│   ├── Moments_IO_test.ipynb         # Inspect and compare moment analysis results
│   ├── H_two_model.ipynb             # Visualize H2 mass-dependent fit results
│   ├── H_four_model.ipynb            # Visualize H4 mass-dependent fit results
│   └── ConformalPolynomial_study.ipynb  # Explore conformal polynomial parameterization
│
├── set_env.sh                        # Set MOMS2AMPS, Moments, ROOT_INCLUDE_PATH
└── pixi.toml                         # Python environment specification
```

---

## Setup

### 1. Clone and set up the Moments package

```bash
git clone https://github.com/bgrube/Moments.git
cd Moments
pip install -r requirements.txt
```

> If running on a cluster other than ifarm, edit `OpenMpUtilities.py` line 47 and replace `"ifarm"` with your hostname.

### 2. Build brufit

Follow the [brufit build instructions](https://github.com/dglazier/brufit/tree/R6.34Test). Make sure the BruFit environment script (`set_env.sh`) is sourced before running any macros.

### 3. Set environment variables

Edit [`set_env.sh`](set_env.sh) to point to your local installations:

```bash
export Moments=/path/to/Moments/
export MOMS2AMPS=/path/to/Moments2Amplitudes/
export ROOT_INCLUDE_PATH="${MOMS2AMPS}/include:${ROOT_INCLUDE_PATH}"
```

Then:

```bash
source set_env.sh
```

### 4. Set up the Python environment (pixi)

```bash
# Install pixi if needed: https://pixi.sh
pixi install
pixi shell
```

---

## Workflow

### Path A — Linear Algebra Method (mass-independent, per-bin)

This path extracts amplitudes bin-by-bin via matrix inversion of the moment system.

**Step 1 — Generate weighted phase-space MC** (for acceptance correction):

```bash
python macros/WeightPhSp.py --seed 12345 --tag test_run
```

**Step 2 — Calculate moments** from data:

```bash
python macros/CalcMoments_IO_test.py
```

Output: `Moments.out/<config>/` with moment pickle files and integral matrices.

**Step 3 — Plot moments**:

```bash
python macros/PlotMoments_IO_test.py
```

**Step 4 — Export covariance matrices** for brufit:

```bash
python macros/ExportCovariance_IO_test.py
# or directly:
./macros/ExportCovariance.sh <dir> <matrixFile> <outputName>
```

---

### Path B — Mass-Dependent Fit (global, continuous in mass)

This path fits a physically-motivated model to all mass bins simultaneously.

**Step 1 — Obtain moment fit results from brufit.**

Run brufit's MCMC moment fitter to get `ResultsBruMcmcCovariance.root` per mass bin.

**Step 2 — Run a mass-dependent model macro.**

Choose a model macro from `macros/`:

```bash
# Example: S-wave Flatté + D-wave Breit-Wigner, single seed
root -l -b -q ${BRUFIT}/macros/LoadBru.C macros/Load.C 'macros/Mom2Amp_SFlatteDBreitWigner.C(0)'

# Example: run 100 seeds in parallel with organized output
./macros/run_model.sh -n my_fit -s 1 -e 100 -j 40 Mom2Amp_SFlatteDBreitWigner.C
```

The `run_model.sh` script creates a timestamped output directory under `run/` with:
```
run/my_fit_YYYYMMDD_HHMMSS/
  ├── macro/      # snapshot of macros used
  ├── logs/       # per-seed logs
  ├── results/    # result ROOT trees
  └── metadata.txt
```

**Step 3 — Visualize results** in Jupyter:

```bash
pixi run jupyter notebook notebook/H_two_model.ipynb
```

---

## Key Classes

| Class | File | Description |
|---|---|---|
| `m2pw::MassDependentFitter` | [`include/MassDependentFitter.h`](include/MassDependentFitter.h) | Top-level fitter; sets up chi2, runs Minuit2, manages parameters |
| `m2pw::MassDependentFunction` | [`include/MassDependentFunction.h`](include/MassDependentFunction.h) | Evaluates complex amplitude for a given line-shape model |
| `MassDependentMoments` | [`include/MassDependentMoments.h`](include/MassDependentMoments.h) | Loads brufit moment results and covariances across mass bins |
| `MomentHelper` | [`include/MomentHelper.h`](include/MomentHelper.h) | Per-bin moment access: values, errors, covariance matrix |
| `Equation` / `EquationSolver` | [`include/Equation.h`](include/Equation.h) | Encodes the system $H_{L,M} = f(\{T_{J,m}^\epsilon\})$ |

### `MassDependenceConfig` model types

```cpp
// Configure which partial waves get which model
m2pw::MassDependentFitter::MassDependenceConfig config;
config.massDependentWaves[0] = {{"a0_980", MDModel::FlattePlusPolynomial, 2}};  // S-wave
config.massDependentWaves[2] = {{"a2_1320", MDModel::BreitWigner}};             // D-wave
config.massIndependentL = {1};  // P-wave treated mass-independently
```

---

## Analysis Configurations

Pre-defined analysis configurations live in [`Moments.config/Config_IO_test.py`](Moments.config/Config_IO_test.py):

| Config | Description |
|---|---|
| `CFG_IO_TEST_LMAX4_ZEROPWAVES` | $L_\text{max}=4$, no P-waves, single bin $m \in [1.28, 1.32]$ GeV |
| `CFG_IO_TEST_LMAX4_WITHPWAVES` | $L_\text{max}=4$, with P-waves, single bin $m \in [1.28, 1.32]$ GeV |
| `CFG_IO_TEST_LMAX6` | $L_\text{max}=6$, single bin $m \in [1.28, 1.32]$ GeV |

---

## References

- B. Grube, [Moments framework](https://github.com/bgrube/Moments)
- D. Glazier, [brufit](https://github.com/dglazier/brufit)
- GlueX Collaboration, [Physics book](https://arxiv.org/abs/1112.1734)
- Caprini et al., *Nucl. Phys. B* **530** (1998) 153 — conformal mapping for form factors
- COMPASS, [arXiv:1509.00992](https://arxiv.org/abs/1509.00992) — conformal polynomials in PWA
