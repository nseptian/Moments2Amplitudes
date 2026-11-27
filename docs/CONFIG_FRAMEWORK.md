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

> [!NOTE]
> The `environment:` block in YAML config files is completely optional and no longer required. BruFit and ROOT paths are resolved automatically from your active Pixi environment, `.env`, or `$BRUFIT`/`$ROOTSYS`.

---

## 2. CLI Commands & Pixi Tasks

Pre-configured tasks are available in `pixi.toml` (and directly via `mom2amp`):

### 1. Validate a Configuration
Check YAML syntax, file paths, wave definitions, and mass binning without running:
```bash
pixi run validate configs/d_breit_wigner.yaml
```

### 2. Multi-Start Global Minimum Fit
Run randomized starts in parallel to find the central global minimum:
```bash
pixi run fit configs/d_breit_wigner.yaml
```
Optional CLI overrides:
```bash
pixi run fit configs/d_breit_wigner.yaml --n-starts 50 --n-threads 12 --output-dir my_run

# Parallel execution across all or specific t-bins:
pixi run fit configs/a2CS_D_BreitWigner.yaml --tbins all
pixi run fit configs/a2CS_D_BreitWigner.yaml --tbins t010020 t020032 t032050 --tbin-workers 3
```

### 3. Display Physics Summary Report
Format and print a physics summary table (separated by resonance parameters, couplings, piecewise waves, and fit statistics) and write reports to `reports/fit_summary.txt` and `reports/fit_summary.md`:
```bash
pixi run summary configs/d_breit_wigner.yaml
```
Optional output specification:
```bash
pixi run summary configs/d_breit_wigner.yaml --fit-file logs_dbw/fit/nominal_best.root
```

### 4. Generate Diagnostic Plots & Moment Comparison
Plot the \(\chi^2\) distribution landscape, fitted amplitude magnitudes, and overlay mass-dependent fit moments against experimental unnormalized MCMC moments:
```bash
pixi run plot configs/d_breit_wigner.yaml
```
Optional CLI overrides:
```bash
pixi run plot configs/s_piecewise_d_bw.yaml --n-lowest-chi2 10 --moments L=2
```
Saves publication-ready PNG and PDF figures in `<output_dir>/plots/` and individual moments in `<output_dir>/plots/moments/`.

### 5. Calculate Reflectivity-Separated Cross-Section
Calculate acceptance-corrected yields and differential cross section \(\frac{d\sigma}{dt}\) separated by reflectivity for a fitted resonance (e.g. \(a_2(1320)\)):
```bash
pixi run cs configs/s_piecewise_d_bw.yaml
```
Optional CLI overrides:
```bash
pixi run cs configs/s_piecewise_d_bw.yaml --resonance a2_1320 --t-min 0.1 --t-max 0.2 --luminosity 257.4
```
For multi-\(t\)-bin analysis, pass multiple run directories:
```bash
pixi run cs configs/s_piecewise_d_bw.yaml --runs logs_t010020 logs_t020032 logs_t032050 logs_t050075
```
Saves reports to `<output_dir>/reports/cross_section_{resonance}.txt`, `.md`, `.json`, and publication plot `<output_dir>/plots/{resonance}_CrossSection_vs_t.png` and `.pdf`.

### 6. MCMC Posterior Uncertainty Propagation
Propagate experimental uncertainty across MCMC posterior samples using the global minimum best parameters:
```bash
pixi run propagate configs/s_piecewise_d_bw.yaml
```
Optional CLI overrides:
```bash
pixi run propagate configs/s_piecewise_d_bw.yaml --n-workers 8 --max-entries 200
```

### 7. End-to-End Pipeline
Run multi-start fitting followed directly by MCMC propagation:
```bash
pixi run mom2amp run configs/s_piecewise_d_bw.yaml
```

---

## 3. Standard Run Directory Layout

All outputs from fits, reports, plots, and MCMC propagation are organized into a clean, predictable hierarchy managed by `rundir.py`:

```
<output_dir>/                      # Configured via YAML `multi_start.output_dir` or `--output-dir`
├── fit/                           # Fit outputs and parameters
│   ├── starts/                    # Individual randomized start ROOT files and logs
│   │   ├── start_000001.root
│   │   ├── start_000001.log
│   │   ├── start_000002.root
│   │   └── ...
│   ├── nominal_best.root          # Lowest-chi2 best fit file (used for summaries and downstream stages)
│   └── nominal_best_params.root   # Extracted best-fit parameters (used for MCMC propagation)
│
├── reports/                       # Human-readable and Markdown reports
│   ├── fit_summary.txt            # Plain ASCII summary table
│   ├── fit_summary.md             # GitHub-flavored Markdown summary report
│   └── fit.log                    # Combined driver execution log
│
├── plots/                         # Generated figures
│   ├── chi2_distribution.png     # High-resolution PNG of chi2 distribution across starts
│   └── chi2_distribution.pdf     # Vector PDF for publication/presentation
│
├── mcmc_fits/                     # Individual worker fit outputs during MCMC propagation
│   └── worker_*.root
│
└── config.yaml                    # Frozen copy of configuration used for the run
```

---

## 4. Two-Stage Fitting Workflow

In photoproduction of two pseudoscalar mesons (e.g. \(\gamma p \to \eta\pi^0 p\)), higher moments like \(H_4\) are purely dominated by spin \(J \ge 2\) amplitudes (\(D\)-wave), with zero contribution from \(S\)-wave (\(J=0\)) or \(P\)-wave (\(J=1\)). This allows a robust, decoupled **Two-Stage Analysis**:

```
 ┌────────────────────────────────────────────────────────┐
 │ Stage 1: Fit Pure D-wave Breit-Wigner to H4            │
 │ • Config: configs/d_breit_wigner.yaml                 │
 │ • Extract clean resonance parameters: M(a2), Gamma(a2) │
 │ • Output: logs_dbw/fit/nominal_best.root              │
 └──────────────────────────┬─────────────────────────────┘
                            │
                            │ Initialize D-wave (init.strategy: "from_file")
                            ▼
 ┌────────────────────────────────────────────────────────┐
 │ Stage 2: Fit S-wave Piecewise + D-wave to H0, H2, H4   │
 │ • Config: configs/s_piecewise_d_bw.yaml               │
 │ • S-wave: Mass-independent magnitude/phase per bin     │
 │ • D-wave: Initialized or fixed from Stage 1 best fit   │
 │ • Output: logs_spiecewise/fit/nominal_best.root        │
 └────────────────────────────────────────────────────────┘
```

### Stage 1: D-Wave Breit-Wigner Fit (`configs/d_breit_wigner.yaml`)
Fits only \(H_4\) moments (`included_L: [4]`) with \(L=2\) Breit-Wigner and `poly_order: -1` (disabling polynomial background):
```bash
pixi run fit configs/d_breit_wigner.yaml
```

### Stage 2: S-Piecewise + D-Breit-Wigner Fit (`configs/s_piecewise_d_bw.yaml`)
Fits \(H_0, H_2, H_4\) moments (`included_L: [0, 2, 4]`). The \(D\)-wave resonance parameters are initialized directly from Stage 1:
```yaml
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "breit_wigner"
    formula: "k / (M^2 - s - i * M * width)"
    parameters: ["k", "M", "width"]
    poly_order: -1
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
      file_path: "logs_dbw/fit/nominal_best.root"
      option: 0
```
Run Stage 2 with:
```bash
pixi run fit configs/s_piecewise_d_bw.yaml
```

---

## 5. Configuration File Reference (`config.yaml`)

### Example: Piecewise S-waves + D-wave Breit-Wigner
```yaml
project:
  name: "EtaPi0_SPiecewise_DBreitWigner"
  description: "Stage 2: Piecewise S-wave + D-wave Breit-Wigner fit"

# Dataset moments and covariance settings
dataset:
  moments_dir: "/path/to/MomentFit"
  t_bin: "t010020"       # Configurable t-bin directory (e.g. t010020, t020032, or index 0..4)
  results_filename: "ResultsBruMcmcFixedCovariance.root"
  bin_prefix: "Mpi0eta"
  var_bin_prefix: ""
  mass_bins:
    min: 1.04
    max: 1.72
    step: 0.04
  use_covariance: true

# Moments included in chi2 calculation
chi2:
  included_L: [0, 2, 4]  # H0, H2, H4 moments
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
        max: 350.0           # Calibrated above S-wave a0(980) peak (~314)
      phase:
        strategy: "random"   # Initial draw in [-pi, pi], kept unbounded in Minuit (no hard wall)
        min: -3.14159265
        max: 3.14159265
        fixed: false

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
      file_path: "logs_dbw/fit/nominal_best.root"
      spread: 0.25           # Perturb loaded couplings/phases across multi-start seeds (optional)
      coupling:
        strategy: "random"   # Defaults for random/missing couplings: [-80, 80] (or poly_scale/poly_growth)
        min: -80.0
        max: 80.0
      option: 0

# Numerical Minimizer Settings (Minuit2)
minimizer:
  algorithm: "Minuit2"
  strategy: 1
  tolerance: 1.0e-3
  max_calls: 1000000
  print_level: 0
  reject_at_limit: true      # Mark fit as not converged if any free parameter stops at its limit
  limit_tolerance: 1.0e-5    # Relative distance to [lower, upper] considered "at limit"

# Multi-Start Randomized Central Minimum Search
multi_start:
  enabled: true
  n_starts: 100              # Total trials (used if seed_range not set)
  seed_range: [1, 100]       # Optional explicit [start, end] seed range (e.g. [5001, 10000] or "5001:10000")
  seed_start: 1              # Or specify starting seed (1-based, inclusive)
  seed_end: 100              # Or specify ending seed (inclusive)
  skip_existing: false       # Skip seeds whose fit output already exists in fit/starts/
  n_threads: 8
  output_dir: "logs_spiecewise"
  plateau_delta_chi2: 0.05
  seed_offset: 1
  n_lowest_chi2: 10

# Moment & Amplitude Plotting Settings
plotting:
  enabled: true
  n_lowest_chi2: 10          # Number of lowest chi2 fit curves to draw
  draw_best_fits: true      # Overlay best fit curves on data moments
  unnormalized_moments_file: "/path/to/UnnormalizedMoments_vs_Mass.root" # (optional, auto-detected)
  moments: ["all"]          # "all", "L=2", "L=4", or list: ["H_0_0_0", "H_0_2_0"]
  formats: ["png", "pdf"]
  moments_vs_t: true        # Side-by-side 1 x N_tbins moment comparison plots
  moments_vs_t_sharey: false
  amplitudes_vs_t: true     # Side-by-side 1 x N_tbins partial wave intensity & phase comparison plots
  amplitudes_vs_t_sharey: false
  distributions: true       # 1D posterior distributions of MCMC parameters and yields
  distribution_bins: 50     # Histogram bins for 1D MCMC posterior distributions
  chi2_landscape_bins: 60   # Histogram bins for multi-start Chi2 landscape plots

# Reflectivity-Separated Resonance Cross-Section Settings
cross_section:
  enabled: true
  resonance: "a2_1320"      # Target resonance name from waves[].resonances
  branching_ratio: 0.0565   # Total BR: BR(res -> final state) * BR(final state decay)
  luminosity: 257.4         # Integrated luminosity in pb^-1
  luminosity_unit: "pb^-1"  # "pb^-1" (factor 1e6 to mub^-1), "mub^-1", or "nb^-1"
  t_min: 0.1                # Minimum -t in GeV^2 (auto-detected from run directory if omitted)
  t_max: 0.2                # Maximum -t in GeV^2
  mass_limits: [0.8, 2.0]   # Integration limits [m_min, m_max] in GeV
  phase1_comparison: true   # Overlay published Phase-1 PWA values

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

## 6. Covariance Matrix Extraction, Inversion, & Robust Singular Fallback

### Moment Covariance Calculation from MCMC Trees
Covariance matrices are computed directly from the BruFit `MCMCTree` saved in `ResultsBruMcmcCovariance.root` using `m2pw::MomentHelper::ComputeMCMCCovariance`:
1. Physical unnormalized moment samples are calculated per MCMC step $j$:
   $$H_{i, j} = \frac{Y_j \cdot h_{i, j}}{\eta}, \quad H_{0, j} = \frac{2 Y_j}{\eta}$$
   where $Y_j$ is the acceptance-corrected event yield, $h_{i, j}$ are normalized spherical harmonic moments, and $\eta$ is the phase-space acceptance integral.
2. The empirical sample covariance matrix across $N$ posterior samples is:
   $$V_{ik} = \frac{1}{N - 1} \sum_{j=1}^N (H_{i, j} - \bar{H}_i)(H_{k, j} - \bar{H}_k)$$

### Full Covariance vs. Diagonal Weighting
In the YAML configuration under `dataset`:
- `use_covariance: true` (default): Uses the full inverse covariance matrix $W = V^{-1}$ per mass bin:
  $$\chi^2 = \sum_{\text{bins}} (\vec{H}^{\text{meas}} - \vec{H}^{\text{theo}})^T V^{-1} (\vec{H}^{\text{meas}} - \vec{H}^{\text{theo}})$$
- `use_covariance: false`: Uses only the diagonal variance terms ($W_{ii} = 1/\sigma_i^2$), ignoring correlations.

### Handling Singular Matrices & Automatic Fallback
When fitting a subset of angular moments (such as fitting only $H_4$ in Stage 1 without $H_0, H_2$), the corresponding covariance submatrix may be rank-deficient or ill-conditioned ($det(V) \approx 0$).
`m2pw::Chi2Engine` handles this seamlessly:
- LU decomposition inversion (`TDecompLU`) errors are intercepted with ROOT's error suppression (`gErrorIgnoreLevel = kBreak`) to prevent stderr clutter.
- If matrix inversion fails or $det(V) \le 0$, `Chi2Engine` prints an informative log notice:
  ```
  [FIT][CHI2] mass_bin=1.06 covariance_matrix_singular (det=0) fallback_to_diagonal
  ```
- It automatically falls back to diagonal weighting $W_{ii} = 1/\sigma_i^2$, allowing Minuit2 minimization to proceed cleanly to convergence without crashing or stalling.

---

## 7. Polynomial Background Deactivation (`poly_order: -1`)

For clean resonance parameter extraction (e.g., measuring pure $a_2(1320)$ mass and width in Stage 1), an empirical polynomial background is unnecessary and can introduce artificial parameter correlations.

Setting `poly_order: -1` in the YAML configuration cleanly deactivates the polynomial contribution across all layers:
1. **Model Line Shape**: In `MassDependentFunction`, `hasPoly = (_PolyOrder >= 0)` evaluates to `false`, fixing the polynomial amplitude contribution strictly to $A_{\text{poly}} \equiv 0$.
2. **Minuit Parameter Registry**: `ParameterManager::BuildParamTypesForWave` skips allocating polynomial parameters (`m_threshold`, `m_expansion`, `re_*`, `im_*`). Minuit2 minimizes only the physical resonance parameters ($k, M, \Gamma$).
3. **Summary Filtering**: `mom2amp summary` automatically identifies `poly_order: -1` and suppresses non-existent polynomial parameters from the printed ASCII tables and Markdown reports.

---

## 8. Explicit & Custom Mass-Dependent Model Formulas

The framework allows physicists to define the **explicit mathematical formula** directly in the YAML configuration. The framework compiles these expressions directly into high-speed C++ lambdas evaluated during Minuit2 minimization.

### Supported Built-in Formulas
* **Breit-Wigner**:
  ```yaml
  model: "breit_wigner"
  formula: "k / (M^2 - s - i * M * width)"
  parameters: ["k", "M", "width"]
  ```
* **Realistic Breit-Wigner with Kinematic Phase-Space Factor** (JPAC PRD 100, 054017):
  Multiplies the Breit-Wigner lineshape by $\sqrt{\rho_{\eta\pi}(m)}$ so the physical amplitude vanishes smoothly at the $\eta\pi$ reaction threshold:
  ```yaml
  model: "breit_wigner"
  formula: "sqrt(rho_etapi) * k / (M^2 - s - i * M * width)"
  parameters: ["k", "M", "width"]
  ```
  *(Available in `configs/a2CS_D_BreitWigner_Kinematic.yaml` and `configs/a2CS_SPiecewise_DBreitWigner_Kinematic.yaml`)*
* **Relativistic Blatt-Weisskopf Breit-Wigner with Centrifugal Barrier Factors**:
  Incorporates the Blatt-Weisskopf centrifugal barrier penetration factor ratio $B_L(q, R)/B_L(q_0, R)$ in the numerator and the energy-dependent dynamic width $\Gamma_L(m) = \Gamma_0 (M_0/m)(q/q_0)[B_L(q, R)/B_L(q_0, R)]^2$ in the denominator:
  ```yaml
  model: "barrier_breit_wigner"
  barrier: true
  barrier_l: 2
  barrier_radius: 5.0677
  formula: "sqrt(rho_etapi) * barrier_bw_etapi(k, M, width, 2)"
  parameters: ["k", "M", "width"]
  ```
  *(Available in `configs/a2CS_D_BreitWigner_Barrier.yaml`)*
  Additional barrier BW helper variants supported inside `formula`:
  - `peak_barrier_bw_etapi(k, M, width, L)`: multiplies numerator by $M_0 \Gamma_0$ (peak-normalized, dynamic $\eta\pi$ width)
  - `const_barrier_bw_etapi(k, M, width, L)`: uses constant total width $\Gamma_0$ in the denominator while keeping the $\eta\pi$ barrier ratio $B_L(q)/B_L(q_0)$ in the numerator
  - `peak_const_barrier_bw_etapi(k, M, width, L)`: combines peak normalization ($M_0 \Gamma_0$ in numerator) with constant total width $\Gamma_0$ in the denominator and $\eta\pi$ barrier ratio $B_L(q)/B_L(q_0)$ in the numerator (used for broad multi-channel states like $a_2(1700)$ in `configs/SPiecewise_DTwoBreitWigner_FitA21700.yaml`)
* **Coherent Sum of Two Barrier Breit-Wigners** ($a_2(1320) + a_2(1700)$):
  ```yaml
  model: "barrier_two_breit_wigner"
  barrier: true
  barrier_l: 2
  barrier_radius: 5.0677
  formula: "sqrt(rho_etapi) * (barrier_bw_etapi(k1, Mass1, Width1, 2) * exp(i * phase1) + barrier_bw_etapi(k2, Mass2, Width2, 2) * exp(i * phase2))"
  parameters: ["k1", "Mass1", "Width1", "phase1", "k2", "Mass2", "Width2", "phase2"]
  ```
  *(Available in `configs/a2CS_published_model_Barrier.yaml`)*
* **Coherent Sum of Two Breit-Wigners**:
  Available in `configs/d_two_breit_wigner.yaml` (and with kinematic factor in `configs/a2CS_D_TwoBreitWigner_Kinematic.yaml`):
  ```yaml
  model: "two_breit_wigner"
  formula: "sqrt(rho_etapi) * (k1 / (Mass1^2 - s - i * Mass1 * Width1) + (k2 / (Mass2^2 - s - i * Mass2 * Width2)) * exp(i * phase))"
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

### Gaussian Parameter Constraints (`gaussian_constraints`)
You can apply Gaussian priors on any resonance mass, width, coupling, or global phase directly under a `resonances:` entry in YAML. When constrained with `(mean, sigma)`, the minimizer adds $(\theta - \mu)^2 / \sigma^2$ to total $\chi^2$.

Make sure the parameter is floated (`*_fixed: false`) so Minuit can vary it under the constraint:
```yaml
waves:
  - l: 2
    name: "D-wave"
    type: "mass_dependent"
    model: "barrier_two_breit_wigner"
    resonances:
      - name: "a2_1320"
        mass: 1.3182
        mass_fixed: false
        mass_limits: [1.20, 1.45]
        width: 0.107
        width_fixed: false
        width_limits: [0.05, 0.20]
        gaussian_constraints:
          mass: [1.3182, 0.0006]       # [mean, sigma] in GeV
          width: {mean: 0.107, sigma: 0.005}
      - name: "a2_1700"
        mass: 1.706
        width: 0.38
        gaussian_constraints:
          mass: [1.706, 0.020]
          width: [0.380, 0.050]
```
Supported keys in `gaussian_constraints`:
- Single Breit-Wigner / Flatté: `mass` (`m`), `width` (`w`), `mass_flatte`, `mass_bw`, `width_bw`, `g_s` (`g_etapi`), `g_kk`, `g_etaprimepi`
- Two-Breit-Wigner: `mass`/`width` automatically map to `Mass1`/`Width1` (resonance index 1) or `Mass2`/`Width2` (resonance index 2), or specify `mass1`/`mass2`/`width1`/`width2` explicitly
- Any exact internal Minuit parameter name (e.g. `MD_shared_a2_1320_Mass`, `MD_shared_a2_1320_1700_Width2`, `MD_shareda_a2_1320_1700_globalphase2`, `MD_sharedb_a2_1320_1700_globalphase2`)

---

## 9. Physics Summary Reports & Diagnostic Plotting

### Formatted Physics Summary Table (`mom2amp summary`)
The `summary` command inspects the nominal best fit ROOT file (`<output_dir>/fit/nominal_best.root`), extracts Minuit parameters, calculates goodness-of-fit statistics, and generates both plain text and Markdown reports:

```bash
pixi run summary configs/d_breit_wigner.yaml
```

The report provides:
1. **Fit Quality Overview**:
   - Total \(\chi^2\), Degrees of Freedom (NDF), and reduced \(\chi^2 / \text{NDF}\).
   - Minuit status, EDM (estimated distance to minimum), and convergence state.
2. **Resonance Parameters**:
   - Mass \(M\) and width \(\Gamma\) (in GeV) with Minuit symmetric error, initial value, bounds, and status (Free/Fixed).
3. **Couplings & Piecewise Amplitudes**:
   - Coupling constants \(k\) and per-mass-bin piecewise magnitudes and phases.
4. **Exported Reports**:
   - Plain ASCII: `<output_dir>/reports/fit_summary.txt`
   - GitHub-flavored Markdown: `<output_dir>/reports/fit_summary.md`

### Diagnostic Plots & Moment Comparison (`mom2amp plot`)
The `plot` command evaluates the nominal best fit, all randomized starts in `<output_dir>/fit/starts/`, and compares reconstructed moments with experimental MCMC moments:

```bash
pixi run plot configs/s_piecewise_d_bw.yaml
```

CLI options:
- `--n-lowest-chi2 <N>`: Draw the top \(N\) lowest \(\chi^2\) fit curves (Rank 1 solid red, ranks 2..N semi-transparent red).
- `--moments <filter>`: Filter moments to plot (`"all"`, `"L=2"`, `"L=4"`, or comma-separated list `"H_0_0_0,H_0_2_0"`).
- `--unnormalized-moments-file <path>`: Path to `UnnormalizedMoments_vs_Mass.root` (auto-detected if omitted).
- `--formats <png,pdf>`: Output graphic formats (default: `png,pdf`).

Generated Figures:
1. **Multi-Start \(\chi^2\) Landscape**:
   - Histogram of minimized \(\chi^2\) values across starts.
   - Highlights the global minimum plateau (all starts within \(\Delta\chi^2 \le \text{plateau\_delta\_chi2}\) of the best fit).
   - Saved to `<output_dir>/plots/multistart_chi2_landscape.png` and `.pdf`.
2. **Individual Moment Comparison Figures**:
   - Experimental MCMC moments plotted as black markers with error bars.
   - Top \(N\) lowest \(\chi^2\) mass-dependent fit curves overlaid: Rank 1 (solid red line) and Ranks 2..N (semi-transparent red lines).
   - Saved to `<output_dir>/plots/moments/{moment_name}.png` and `.pdf` (e.g. `H_0_0_0.png`, `H_0_2_0.png`).
3. **Multi-Panel Moments Overview Grid**:
   - Comprehensive multi-panel overview grid of all requested angular moments.
   - Saved to `<output_dir>/plots/moments_overview.png` and `.pdf`.
4. **Fitted Amplitude Magnitudes & Phases**:
   - Saved to `<output_dir>/plots/amplitudes_magnitude_vs_mass.png` and `.pdf` as a two-panel publication figure sharing the mass axis:
     - **Top panel**: Amplitude intensities $|V_\alpha|^2$ and total fit intensity $\sum |V|^2$ overlaid on acceptance-corrected yield data $H_0(0,0)$ with error bars.
     - **Bottom panel**: Wave phases $\phi_\alpha$ as a function of mass with MCMC 68% credible intervals.
   - Also saved to `amplitudes_intensity_vs_mass.png` and standalone `amplitudes_phase_vs_mass.png`, along with individual wave two-panel plots in `<output_dir>/plots/amplitudes/individual/{wave}_intensity_and_phase.png`.
5. **Grouped Wave Overview Panels (Shared Y-Axis & Aligned Phases)**:
   - Saved to `<output_dir>/plots/amplitudes/overview_{D_waves_positive,D_waves_negative,S_waves,all_waves}.png` and `<output_dir>/plots/amplitudes_overview_vs_mass.png`.
   - All wave intensity panels share the exact same y-axis with the background $H_0(0,0)$ data, removing redundant y-ticks and labels across columns to use space efficiently.
   - Corresponding wave phases are aligned below sharing the phase y-axis and the mass x-axis.
6. **Multi-$t$ Resonance Parameters & $S,P,D$ Wave Intensities (Total & Reflectivity-Separated)**:
   - Saved under `<base_dir>/plots/amplitudes_vs_t/` (when multiple `-t` bins are present):
     - `spd_waves_refl_separated_vs_mass_all_tbins.{png,pdf}`: $2 \times N_t$ reflectivity-separated mass spectra showing natural-parity ($\epsilon = +1$: $S^+, P^+, D^+$) on top and unnatural-parity ($\epsilon = -1$: $S^-, P^-, D^-$) on bottom, overlaid with $H_0(0,0)$ data and Total $\sum |V|^2$ on a shared y-axis.
     - `spd_waves_refl_separated_vs_t.{png,pdf}`: Integrated reflectivity-separated wave yields ($S^+, P^+, D^+$ vs $S^-, P^-, D^-$) and shares of total intensity vs $-t$.
     - `spd_waves_total_reflectivity_vs_t.{png,pdf}`: Total $|V|^2$ across both reflectivities ($\epsilon = \pm 1$) and $M$ projections for $S$, $P$, $D$ waves plus Total $\sum |V|^2$ and $H_0(0,0)$ data vs $-t$, with fractional intensity breakdown below.
     - `spd_waves_vs_mass_all_tbins.{png,pdf}`: Mass spectrum comparison of total-reflectivity $S$, $P$, $D$ waves across all $-t$ bins.
     - `resonance_parameters_vs_t.{png,pdf}`, `resonance_<name>_params_vs_t.{png,pdf}`, and `parameters_vs_t/<param>_vs_t.{png,pdf}`: Fitted resonance masses, widths, couplings, and global phases as a function of $-t$ with Minuit / MCMC error bars and configured bounds/initial lines.

---

## 10. Reflectivity-Separated Resonance Cross-Section Calculation

The framework extracts the physical acceptance-corrected yields and differential cross section \(\frac{d\sigma}{dt}\) separated by natural parity exchange (positive reflectivity \(\epsilon = +1\)) and unnatural parity exchange (negative reflectivity \(\epsilon = -1\)) from the fitted mass-dependent resonance parameters.

### 1. Mathematical Formalism
For each partial wave \(w\) carrying the resonance (e.g. \(D\)-waves \(a_2^m\) for \(\epsilon=+1\) and \(b_2^m\) for \(\epsilon=-1\)), the resonance coupling \(k_w\) is extracted from `mass_dependent_params` (or `mcmc_params`).

1. **Acceptance-Corrected Resonance Yield per Wave**:
   The acceptance-corrected yield attributed to the resonance in wave \(w\) is obtained by integrating the continuous Breit-Wigner intensity over the mass spectrum and normalizing by the invariant mass bin width \(\Delta m\):
   $$Y_w = \frac{1}{\Delta m} \int_{m_{\text{min}}}^{m_{\text{max}}} \frac{k_w^2}{(m^2 - M^2)^2 + (M \Gamma)^2} \, dm$$
   Because \(Y_w \propto k_w^2\), linear error propagation from the Minuit symmetric parameter error \(\sigma(k_w)\) yields:
   $$\sigma(Y_w) = 2 \frac{Y_w}{|k_w|} \sigma(k_w)$$

2. **Reflectivity Sums**:
   Yields are summed independently over positive (\(a_L^m\)) and negative (\(b_L^m\)) reflectivity waves:
   $$Y^{(+)} = \sum_{w \in \text{PosRefl}} Y_w, \quad \sigma(Y^{(+)}) = \sqrt{\sum_{w \in \text{PosRefl}} \sigma^2(Y_w)}$$
   $$Y^{(-)} = \sum_{w \in \text{NegRefl}} Y_w, \quad \sigma(Y^{(-)}) = \sqrt{\sum_{w \in \text{NegRefl}} \sigma^2(Y_w)}$$
   $$Y^{(\text{tot})} = Y^{(+)} + Y^{(-)}, \quad \sigma(Y^{(\text{tot})}) = \sqrt{\sigma^2(Y^{(+)}) + \sigma^2(Y^{(-)})}$$

3. **Differential Cross Section \(\frac{d\sigma}{dt}\)**:
   The differential cross section in \(\mu\text{b} / \text{GeV}^2\) is given by:
   $$\frac{d\sigma^{(\pm, \text{tot})}}{dt} = \frac{Y^{(\pm, \text{tot})}}{\text{den\_factor} \cdot \mathcal{L} \cdot \Delta t \cdot \mathcal{B}}$$
   where:
   - \(\mathcal{L}\) is the integrated luminosity (in \(\text{pb}^{-1}\), with \(\text{den\_factor} = 10^6\ \text{pb}/\mu\text{b}\)).
   - \(\Delta t = t_{\text{max}} - t_{\text{min}}\) is the \(-t\) bin width (in \(\text{GeV}^2\)).
   - \(\mathcal{B}\) is the total product branching fraction (e.g. \(\mathcal{B}(a_2 \to \eta\pi^0) \times \mathcal{B}(\eta \to 2\gamma) \times \mathcal{B}(\pi^0 \to 2\gamma) \approx 0.0565\)).

4. **Production Reflectivity Fraction**:
   $$F^{(+)} = \frac{Y^{(+)}}{Y^{(\text{tot})}} = \frac{d\sigma^{(+)}/dt}{d\sigma^{(\text{tot})}/dt}$$

5. **MCMC Posterior Uncertainty Propagation**:
   While Minuit error propagation assumes linearity and uncorrelated waves around the minimum, MCMC uncertainty propagation evaluates the full non-linear integral for each posterior sample \(j \in [1, N]\) from the `mcmc_params` tree:
   $$Y_{w, j} = \frac{1}{\Delta m} \int_{m_{\text{min}}}^{m_{\text{max}}} \frac{k_{w, j}^2}{(m^2 - M_j^2)^2 + (M_j \Gamma_j)^2} \, dm$$
   $$Y^{(+)}_{j} = \sum_{w \in \text{PosRefl}} Y_{w, j}, \quad Y^{(-)}_{j} = \sum_{w \in \text{NegRefl}} Y_{w, j}, \quad \left(\frac{d\sigma}{dt}\right)_j = \frac{Y_j}{\text{den\_factor} \cdot \mathcal{L} \cdot \Delta t \cdot \mathcal{B}}$$
   From the resulting \(N\) posterior samples, the empirical standard deviation and asymmetric 68% credible intervals (16th and 84th percentiles) are calculated:
   $$\sigma_{\text{MCMC}}\left(\frac{d\sigma}{dt}\right) = \text{std}\left(\left\{\left(\frac{d\sigma}{dt}\right)_j\right\}\right), \quad \Delta_{-\text{mcmc}} = \text{median} - p_{16}, \quad \Delta_{+\text{mcmc}} = p_{84} - \text{median}$$
   This captures the exact correlations and covariances between partial waves and resonance mass/width without Gaussian approximations.

### 2. CLI Execution & Multi-\(t\)-bin Evaluation
Run on a single fit result (automatically auto-detecting MCMC propagated files):
```bash
pixi run cs configs/s_piecewise_d_bw.yaml
```
Optionally provide an explicit MCMC file or disable MCMC propagation:
```bash
pixi run cs configs/s_piecewise_d_bw.yaml --mcmc-file mcmc_propagation_spiecewise_results.root
pixi run cs configs/s_piecewise_d_bw.yaml --no-mcmc
```
Run across multiple \(-t\) bins simultaneously to generate the cross section vs \(-t\) distribution:
```bash
pixi run cs configs/s_piecewise_d_bw.yaml --runs logs_t010020 logs_t020032 logs_t032050 logs_t050075
```
Automatically overlays published Phase-1 PWA benchmarks (stat and syst error boxes) and theoretical model bands (`Phys. Rev. D 102, 014003`), plotting MCMC propagated statistical uncertainties on current points and saving figures to `<output_dir>/plots/{resonance}_CrossSection_vs_t.png` and `.pdf`.

---

## 11. Modular C++ Architecture & Optimization

The core C++ engine is organized into focused, modular components:

| Component | Files | Role |
|---|---|---|
| **Chi2Engine** | [`include/Chi2Engine.h`](include/Chi2Engine.h), [`src/Chi2Engine.cpp`](src/Chi2Engine.cpp) | Compiles parameter routing arrays, manages the contiguous \(O(1)\) covariance matrix cache, performs LU matrix inversion with singular fallback, and evaluates \(\chi^2\). |
| **ParameterManager** | [`include/ParameterManager.h`](include/ParameterManager.h), [`src/ParameterManager.cpp`](src/ParameterManager.cpp) | Registers Minuit parameters, applies wave initialization strategies (`WaveInitRule`), handles parameter bounds/steps, and serializes results to ROOT trees. |
| **MCMCFitter** | [`include/MCMCFitter.h`](include/MCMCFitter.h), [`src/MCMCFitter.cpp`](src/MCMCFitter.cpp) | Multi-threaded worker pool orchestrator for evaluating MCMC posterior chains and propagating experimental uncertainties. |
| **MassDependentFitter** | [`include/MassDependentFitter.h`](include/MassDependentFitter.h), [`src/MassDependentFitter.cpp`](src/MassDependentFitter.cpp) | High-level coordinator interfacing with `ROOT::Math::Minimizer` (Minuit2), managing datasets, and controlling fit execution. |

### Parameter Conventions
- **Piecewise parameters**: `MI_<mass_bin>_<wave>` (e.g. `MI_1.320000_a_0_0` for magnitude, `MI_1.320000_aphi_0_0` for phase).
- **Shared resonance parameters**: `MD_shared_<resonance>_<param>` (e.g. `MD_shared_a2_1320_Mass`, `MD_shared_a2_1320_Width`).
- **Wave couplings**: `MD_<wave>_<resonance>_<param>` (e.g. `MD_a_2_1_a2_1320_k`).
- **Configurable Wave Initialization**: `WaveInitRule` structs allow setting custom initialization strategies (`RandomPositiveMagnitudeZeroPhase`, `FixedZero`, `FixedValue`, `FromFile`) per angular momentum \(L\).

### Minimizer Performance Features
1. **Pre-Compiled Static Parameter Routing**:
   - `Chi2Engine` pre-compiles all parameter destinations (mass-independent direct routing, mass-dependent magnitude/phase mappings, and Gaussian constraint penalties) once during initialization into flat integer arrays.
   - All inner-loop string operations (`Tokenize`, `BeginsWith`, `Index`, `TString::Format`) and red-black tree searches (`std::map<TString, ...>`) are completely eliminated from the minimization loop.

2. **\(O(1)\) Contiguous Covariance Matrix Cache**:
   - Covariance matrices and residuals for each mass bin are cached in contiguous memory (`std::vector<BinCovarianceCache>` and `std::vector<equation_t*>`), replacing tree lookups by mass bin `double`.
   - Residual vectors \(\delta_k\) and inverse covariance similarity operations run with direct array indexing and zero heap allocations.

3. **Pre-Assembled Wave Parameter Buffers**:
   - Resonance parameters (e.g. coupling \(k\), mass \(M\), width \(\Gamma\), polynomial coefficients, and global phases) are assembled once per Minuit step into contiguous memory buffers and evaluated across mass bins.

---

## 12. Multi-Start Parallel Execution

Non-linear partial-wave amplitude fits to angular moments possess multiple local minima and saddle points. A multi-start randomized strategy is critical to ensure discovery of the true global minimum:

- The Python driver (`mom2amp fit`) automatically parallelizes the randomized starts across available CPU cores specified by `--n-threads` or `multi_start.n_threads`.
- Each worker executes independently with an isolated random seed, saving individual fit outputs and logs to `<output_dir>/fit/starts/MassDependentResults_<seed>.root` and `<output_dir>/fit/starts/fit_seed_<seed>.log`.
- The driver tracks real-time progress, aggregates results, detects the global minimum plateau (all starts within \(\Delta\chi^2 \le \text{plateau\_delta\_chi2}\)), and links/copies the best start to `<output_dir>/fit/nominal_best.root`.
- Downstream tools (`mom2amp summary`, `mom2amp plot`, and Stage 2 fits) read `<output_dir>/fit/nominal_best.root` directly.

### Continuation Fits & Seed Ranges

If an initial multi-start fit (e.g. 5,000 trials) yields an inadequate plateau or poor clustering, you can continue sampling additional seeds without re-running earlier fits:

1. **Via CLI flags**:
   ```bash
   # Continue fit from seed 5001 to 10000
   mom2amp fit config.yaml --seed-range 5001:10000

   # Or using separate start and end flags
   mom2amp fit config.yaml --seed-start 5001 --seed-end 10000

   # Skip seeds that have already completed in the output directory
   mom2amp fit config.yaml --seed-range 1:10000 --skip-existing
   ```

2. **Via YAML configuration**:
   ```yaml
   multi_start:
     seed_range: [5001, 10000]   # Or "5001:10000"
     # seed_start: 5001
     # seed_end: 10000
     skip_existing: true
   ```

3. **Seamless Aggregation**:
   When running a continuation range like `5001..10000`, the driver executes only the requested seeds, but post-fit analysis scans `<output_dir>/fit/starts/` and combines existing completed fits (e.g. seeds `1..5000`) with newly completed fits (`5001..10000`). The global minimum, plateau counts, summary reports, and \(\chi^2\) landscape plots automatically reflect the entire accumulated sample of 10,000 trials.




