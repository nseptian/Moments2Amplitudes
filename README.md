# Moments2Amplitudes

This repository provides tools for **mass-dependent modeling of angular moments of two-pseudoscalar mesons produced in polarized photon-proton collisions** by imposing mass-dependent models on partial wave amplitudes. The main goal is to fit experimental or simulated moment data using physically-motivated, mass-dependent parameterizations of the underlying amplitudes.

The code and workflow are inspired by and compatible with [brufit](https://github.com/dglazier/brufit/tree/R6.34Test/tutorials/PhotoAmps/TwoSpin0AmpsFromMoments).

## Requirements

- [ROOT](https://root.cern/) (with modules required by brufit installed)
- [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) (PWA-Moment relations)
- [Moments](https://github.com/bgrube/Moments) (for weighting phase space, install python requirements from this repository)

## Repository Structure

```
.
├── WeightPhSp.py                # Main script for phase-space weighting
├── src/                         # C++ source files for mass-dependent fitting
├── include/                     # C++ headers for mass-dependent fitting
├── macros/                      # ROOT macros for running fits and analyses
├── jupyter_notebook/            # Jupyter notebooks for analysis and plotting
├── brufit/                      # brufit scripts
├── samples/                     # Sample data files
└── set_env.sh                   # Environment setup script
```

## How to Use

1. Install and build [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) as described in its documentation.

2. **Clone the [Moments](https://github.com/bgrube/Moments) repository and install its requirements:**
    ```bash
    git clone https://github.com/bgrube/Moments.git
    cd Moments
    pip install -r requirements.txt
    ```
    Change "ifarm" in OpenMpUtilities.py line 47 to your computer hostname

3. Edit set_env.sh and set the path of your Moments and BruFit paths.
4. Source set_env.sh

---

## Step-by-Step Tutorial (Mass-independent Partial Wave Amplitudes estimation from Moments)

Replace file names and parameters as needed in WeightPhSp.py

### 1. Generate Weighted Phase-Space MC Data for BruFit

```bash
python WeightPhSp.py --seed 12345 --tag test_run
```

- This will generate phase-space MC samples and weight them according to a model amplitude set.
- Setup and use brufit/runMakeMomentFitData.sh to produce BruFit formatted data.

### 2. Fit Moments

- Setup and use brufit/MomentFit.C for moments fitting.

### 3. Estimate Partial Wave Amplitudes from Moments

- Set partial wave amplitudes parametrization and path to the moments fit results in brufit/RunGivenBruMoments.C
- In the brufit directory, run:
  
```bash
brufit Load.C RunGivenMoments.C
```

### 4. Analyze and Visualize Results

Use `jupyter_notebook/Moments_IO_test.ipynb` to analyze the results.

---

## Step-by-Step Tutorial of Mass-dependent Modelling of the Moments
Work in progress