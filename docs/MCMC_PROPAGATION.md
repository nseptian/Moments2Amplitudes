# MCMC Uncertainty Propagation & Multi-Start Global Minimum Search

This document explains the multi-start global minimum search and full MCMC uncertainty propagation workflow implemented in **Moments2Amplitudes**.

---

## 1. Motivation

Fitting mass-dependent partial-wave amplitudes directly to measured angular moments involves non-linear bilinear relations:

$$H(L, M) = \sum_{\ell, m, \ell', m'} C_{\ell, m, \ell', m'}^{L, M} \, T_{\ell, m} T_{\ell', m'}^*$$

These relations define a non-convex parameter space with multiple local minima. To ensure reliable physical results:
1. **Multi-Start Search on Central Moments**: Run multiple fits with randomized initial seeds to avoid local minima traps. A flat plateau where many independent starts converge to identical $\chi^2$ ($\Delta \chi^2 \le 0.05$) verifies that the fit reached the genuine global minimum.
2. **MCMC Propagation from Global Minimum**: Once the true central global minimum is found, use it as the initial parameter point for fitting each posterior sample in the brufit `MCMCTree`. This propagates parameter correlations and experimental uncertainties cleanly into physical amplitude intensities and moments.

---

## 2. Fit Workflows

### Workflow A: Pure D-Wave Breit-Wigner (`Mom2Amp_DBreitWigner_MCMC.C`)

Isolates the $D$-wave ($a_2(1320)$) using $L=4$ moments ($H_4$), which only receive contributions from $D \times D$.

```bash
brufit -l -b -q macros/Load.C 'macros/Mom2Amp_DBreitWigner_MCMC.C(100, 0, 200, 1, "mcmc_propagation_results.root", "logs", 40)'
```

* **Step 1**: Runs 100 randomized starts in parallel across 40 CPU cores. Generates `logs/central_fits/log_seed_<seed>.txt`, `logs/central_fits_summary.txt`, and `logs/central_chi2_distribution.png` (plateau + histogram). Saves the best central parameters to `logs/central_best_params.root`.
* **Step 2**: Evaluates 200 MCMC entries in parallel across 40 cores (`RunParallelMCMCFits.sh`). Each entry starts from the best central parameters and saves its log to `logs/mcmc_fits/log_entry_<entry>.txt`. Outputs `mcmc_propagation_results.root`.

---

### Workflow B: S-Wave Piecewise + D-Wave Breit-Wigner (`Mom2Amp_SPiecewiseDBreitWigner.C`)

Fits piecewise mass-independent $S$-waves across 20 mass bins ($0.82\text{ to }1.58\text{ GeV}$) along with mass-dependent $D$-wave ($a_2(1320)$) to moments $H_0, H_2, H_4$.

```bash
brufit -l -b -q macros/Load.C 'macros/Mom2Amp_SPiecewiseDBreitWigner.C(100, 0, 200, 1, "logs/central_best_params.root", "mcmc_propagation_spiecewise_results.root", "logs_spiecewise", 40)'
```

#### Initialization Rules
* **$D$-waves ($\ell=2$)**: Initialized **only** from the best D-wave fit file (`logs/central_best_params.root`). Mass ($1.318\text{ GeV}$) and width ($0.107\text{ GeV}$) are fixed to the best-fit values.
* **Piecewise $S$-waves ($\ell=0$)**:
  * **Magnitudes** ($a_{0,0}$, $b_{0,0}$): Randomized uniformly from positive values $\in [0.0, 200.0]$ for each seed.
  * **Phases** ($\phi_{a_{0,0}}$, $\phi_{b_{0,0}}$): Initialized strictly to **$0.0$**.
* **$P$-waves ($\ell=1$)**: Fixed strictly to **$0.0$**.

---

## 3. Output ROOT Trees

Output files (`mcmc_propagation_results.root` and `mcmc_propagation_spiecewise_results.root`) contain:

| Tree Name | Rows | Contents |
| :--- | :---: | :--- |
| `central_params` | 1 | Best-fit global parameters ($M, \Gamma$, couplings, polynomial coefficients) and errors |
| `central_chi2` | 1 | Best-fit total $\chi^2$, Minuit status code, validity flag |
| `central_result` | $N_{\text{bins}}$ | Nominal amplitudes ($a_{\ell,m}, b_{\ell,m}$) and moments $H(L,M)$ per mass bin |
| `mcmc_params` | $N_{\text{MCMC}}$ | Fitted global parameters for each MCMC posterior sample |
| `mcmc_chi2` | $N_{\text{MCMC}}$ | Total $\chi^2$, status, and validity flag for each MCMC sample |
| `mcmc_result` | $N_{\text{MCMC}} \times N_{\text{bins}}$ | Fitted amplitudes and moments per mass bin for each MCMC sample |
| `c_chi2_dist` | 1 (Canvas) | Diagnostic sorted $\chi^2$ plateau graph and distribution histogram |

---

## 4. Parallelization Architecture

```
Mom2Amp_SPiecewiseDBreitWigner.C
 ├── Step 1: RunParallelSPiecewiseFits.sh (Central Multi-Start)
 │    ├── Worker 1  --> RunSingleSeedSPiecewise.C(seed=1)  --> logs/log_seed_1.txt
 │    ├── ...
 │    └── Worker 40 --> RunSingleSeedSPiecewise.C(seed=100)--> logs/log_seed_100.txt
 │
 └── Step 2: RunParallelMCMCFitsSPiecewise.sh (MCMC Uncertainty Propagation)
      ├── Worker 0  --> RunSingleMCMCWorkerSPiecewise (entries 0-4)   --> chunk_0.root
      ├── ...
      └── Worker 39 --> RunSingleMCMCWorkerSPiecewise (entries 195-199)--> chunk_39.root
           └── Merged into final ROOT file (mcmc_propagation_spiecewise_results.root)
```

Because separate ROOT child processes run concurrently, each process has full thread isolation without lock contention on `TFormula` or `RooWorkspace`, and stdout/stderr is cleanly redirected to dedicated per-fit log files.
