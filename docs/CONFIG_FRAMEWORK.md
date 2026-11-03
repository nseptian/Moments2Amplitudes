# Moments2Amplitudes Configuration Framework & User Guide

The **Moments2Amplitudes** framework provides a declarative, configuration-driven system for mass-dependent partial-wave amplitude analysis of angular moments, complete with parallel randomized multi-start fitting and MCMC posterior uncertainty propagation.

---

## 1. Quick Start with Pixi

### Installation
Clone the repository and enter the directory:
```bash
git clone https://github.com/nseptian/Moments2Amplitudes.git
cd Moments2Amplitudes
pixi install
```

### Configuring BruFit and ROOT
The framework requires [BruFit](https://github.com/dglazier/brufit) (for spherical harmonic models and moment extraction) and CERN ROOT.

You can configure the BruFit and ROOT locations through several mechanisms:

1. **Via `.env` file (Recommended for local developer setups)**:
   Copy `.env.example` to `.env` and specify your paths:
   ```bash
   cp .env.example .env
   ```
   Edit `.env`:
   ```bash
   BRUFIT=/path/to/your/brufit
   # ROOTSYS=/path/to/custom/root  # Optional: leave unset to use pixi conda ROOT
   ```
   Pixi automatically loads `.env` upon environment activation.

2. **Via `pixi.toml`**:
   In `pixi.toml` under `[activation.env]`:
   ```toml
   [activation.env]
   BRUFIT = "${BRUFIT:-/d/grid17/septian/brufit_dev/brufit}"
   ```

3. **Via Shell Environment**:
   ```bash
   export BRUFIT=/path/to/your/brufit
   ```

4. **Via the YAML Configuration File**:
   Under the `environment` block in your fit config:
   ```yaml
   environment:
     brufit_dir: "/path/to/your/brufit"
     root_dir: ""  # Optional: uses conda ROOT if empty
   ```

---

## 2. Running Fits with Pixi Tasks

Pre-configured tasks are available in `pixi.toml`:

### 1. Validate a Configuration
Check syntax, file paths, wave definitions, and mass binning:
```bash
pixi run validate configs/s_piecewise_d_bw.yaml
```

### 2. Multi-Start Global Minimum Fit
Run 100 randomized starts in parallel to find the central global minimum and plot the \(\chi^2\) distribution:
```bash
pixi run fit configs/s_piecewise_d_bw.yaml
```
Optional CLI overrides:
```bash
pixi run fit configs/s_piecewise_d_bw.yaml --n-starts 50 --n-threads 12 --output-dir my_run
```

### 3. MCMC Posterior Uncertainty Propagation
Propagate experimental uncertainty across MCMC posterior samples using the global minimum best parameters:
```bash
pixi run propagate configs/s_piecewise_d_bw.yaml
```
Optional CLI overrides:
```bash
pixi run propagate configs/s_piecewise_d_bw.yaml --n-workers 8 --max-entries 200
```

### 4. End-to-End Pipeline
Run multi-start fitting followed directly by MCMC propagation:
```bash
pixi run mom2amp run configs/s_piecewise_d_bw.yaml
```

---

## 3. Configuration File Reference (`config.yaml`)

### Example: Piecewise S-waves + D-wave Breit-Wigner
```yaml
project:
  name: "EtaPi0_SPiecewise_DBreitWigner"
  description: "Stage 2: Piecewise S-wave + D-wave Breit-Wigner fit"

# Environment configuration (BruFit and ROOT)
environment:
  brufit_dir: "${BRUFIT:-/d/grid17/septian/brufit_dev/brufit}"
  root_dir: "${ROOTSYS:-}"

# Dataset moments and covariance settings
dataset:
  moments_dir: "/path/to/Moments.out/nominal"
  results_filename: "ResultsBruMcmcFixedCovariance.root"
  bin_prefix: "Mpi0eta"
  var_bin_prefix: "t010030"
  mass_bins:
    min: 1.04
    max: 1.72
    step: 0.04
  use_covariance: true

# Moments included in chi2 calculation
chi2:
  included_l: [0, 2, 4]  # H0, H2, H4 moments
  excluded_l: []
  max_l: 4

# Wavesets and parameter initialization
waves:
  # L=0: Piecewise (mass-independent per mass bin)
  - l: 0
    name: "S-wave"
    type: "piecewise"
    init:
      strategy: "random"
      magnitude:
        strategy: "random"
        min: 0.0
        max: 200.0
      phase:
        strategy: "fixed"
        value: 0.0

  # L=1: P-wave fixed to 0
  - l: 1
    name: "P-wave"
    type: "piecewise"
    init:
      fixed_zero: true

  # L=2: D-wave Breit-Wigner
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "breit_wigner"
    formula: "k / (M^2 - s - i * M * width)"  # Explicit formula configured directly in YAML
    parameters: ["k", "M", "width"]
    poly_order: -1  # Disabled polynomial background
    resonances:
      - name: "a2_1320"
        mass: 1.318
        width: 0.105
        mass_fixed: false
        width_fixed: false
        mass_limits: [1.20, 1.45]
        width_limits: [0.05, 0.25]
    init:
      strategy: "from_file"
      file_path: "logs_dbw/central_best_params.root"
      option: 0

# Numerical Minimizer Settings (Minuit2)
minimizer:
  algorithm: "Minuit2"
  strategy: 1
  tolerance: 1.0e-3
  max_calls: 1000000
  print_level: 0

# Multi-Start Randomized Central Minimum Search
multi_start:
  enabled: true
  n_starts: 100
  n_threads: 8
  output_dir: "logs_spiecewise"
  plateau_delta_chi2: 0.05
  seed_offset: 1

# MCMC Uncertainty Propagation
mcmc_propagation:
  enabled: true
  mcmc_tree_path: "/path/to/Moments.out/nominal/Mpi0eta1.320000_t010030_/ResultsBruMcmcFixedCovariance.root"
  burn_in: 0
  max_entries: 200
  stride: 1
  n_workers: 8
  output_file: "mcmc_propagation_spiecewise_results.root"
  log_dir: "logs_spiecewise/mcmc_fits"
```

---

## 4. Explicit & Custom Mass-Dependent Model Formulas

The framework allows physicists to define the **explicit mathematical formula** directly in the YAML configuration. The framework compiles these expressions directly into high-speed C++ lambdas evaluated during Minuit2 minimization.

### Supported Built-in Formulas
* **Breit-Wigner**:
  ```yaml
  model: "breit_wigner"
  formula: "k / (M^2 - s - i * M * width)"
  parameters: ["k", "M", "width"]
  ```
* **Coherent Sum of Two Breit-Wigners**:
  Available in `configs/d_two_breit_wigner.yaml`:
  ```yaml
  model: "two_breit_wigner"
  formula: "k1 / (Mass1^2 - s - i * Mass1 * Width1) + (k2 / (Mass2^2 - s - i * Mass2 * Width2)) * exp(i * phase)"
  parameters: ["k1", "Mass1", "Width1", "k2", "Mass2", "Width2", "phase"]
  ```
* **Flatté (Coupled-Channel)**:
  Available in `configs/s_flatte_d_bw.yaml`:
  ```yaml
  model: "flatte"
  formula: "k / (Mass^2 - s - i * (g_etapi^2 * rho_etapi + g_KK^2 * rho_KK + g_etaprimepi^2 * rho_etaprimepi))"
  parameters: ["k", "Mass", "g_etapi", "g_KK", "g_etaprimepi"]
  ```
  *(where `rho_etapi`, `rho_KK`, `rho_etaprimepi` automatically map to the analytic phase-space factors and handle negative values as CBAR fixed ratios).*

### Custom User Formulas
You can also write custom mathematical amplitude expressions with arbitrary parameters:
```yaml
waves:
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "custom"
    formula: "k * pow(s - 0.7, 0.5) / (M^2 - s - i * M * width)"
    parameters: ["k", "M", "width"]
    resonances:
      - name: "a2_custom"
        mass: 1.318
        width: 0.105
        mass_limits: [1.20, 1.45]
        width_limits: [0.05, 0.25]
```

---

## 5. C++ API Modernization & Parameter Conventions

The underlying C++ library provides clean, standardized parameter management:

- **Piecewise parameters**: `MI_<mass_bin>_<wave>` (e.g. `MI_1.320000_a_0_0` for magnitude, `MI_1.320000_aphi_0_0` for phase).
- **Shared resonance parameters**: `MD_shared_<resonance>_<param>` (e.g. `MD_shared_a2_1320_Mass`, `MD_shared_a2_1320_Width`).
- **Wave couplings**: `MD_<wave>_<resonance>_<param>` (e.g. `MD_a_2_1_a2_1320_k`).
- **Configurable Wave Initialization**: `WaveInitRule` structs allow setting custom initialization strategies (`RandomPositiveMagnitudeZeroPhase`, `FixedZero`, `FixedValue`, `FromFile`) per angular momentum \(L\).

---

## 6. High-Performance Execution Engine & Optimization

The fitting core has been optimized for high-throughput multi-start fitting and massive MCMC posterior evaluation:

1. **Pre-Compiled Static Parameter Routing**:
   - In standard Minuit2 fits, `Chi2Function::operator()` is called tens of thousands of times per minimization.
   - The engine pre-compiles all parameter destinations (mass-independent direct routing, mass-dependent magnitude/phase mappings, and Gaussian constraint penalties) once during initialization into flat integer arrays.
   - All inner-loop string operations (`Tokenize`, `BeginsWith`, `Index`, `TString::Format`) and red-black tree searches (`std::map<TString, ...>`) are completely eliminated from the minimization loop.

2. **\(O(1)\) Contiguous Covariance Matrix Cache**:
   - Covariance matrices and residuals for each mass bin are cached in contiguous `std::vector<BinCovarianceCache>` and `std::vector<equation_t*>`, replacing tree lookups by mass bin `double`.
   - Residual vectors \(\delta_k\) and inverse covariance similarity operations run with direct array indexing and zero heap allocations.

3. **Pre-Assembled Wave Parameter Buffers**:
   - Resonance parameters (e.g. coupling \(k\), mass \(M\), width \(\Gamma\), polynomial coefficients, and global phases) are assembled once per Minuit step into contiguous memory buffers and evaluated across mass bins.

4. **Multi-Threaded Parallel Execution**:
   - Evaluation across mass bins is thread-safe and can be dynamically parallelized using `ROOT::TThreadExecutor` for large datasets.
   - Multi-start fits and MCMC chains run concurrently on multi-core workstations with near-linear scaling.


