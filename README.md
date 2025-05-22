# Moments2Amplitudes

This repository contains code and tools to develop and test algorithms for extracting Partial Wave Amplitudes from Moments of Angular Distributions in two-pseudoscalar mesons produced in polarized photon-proton collisions.

## Requirements

- [ROOT](https://root.cern/) (with modules required by brufit installed)
- [brufit](https://github.com/dglazier/brufit/tree/R6.34Test) (main inversion algorithm)
- [Moments](https://github.com/bgrube/Moments) (for weighting phase space, install python requirements from this repository)

## Repository Structure

```
.
├── WeightPhSp.py                # Main script for phase-space weighting and amplitude extraction
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

## Step-by-Step Tutorial

Below is a template for running a typical analysis workflow. Replace file names and parameters as needed in WeightPhSp.py

### 1. Generate Weighted Phase-Space MC Data for BruFit

```bash
python WeightPhSp.py --seed 12345 --tag test_run
```

- This will generate phase-space MC samples from  and weight them according to a model amplitude set.
- Setup and use brufit/runMakeMomentFitData.sh to produce BruFit formatted data.

### 2. Fit Moments

- Setup and use brufit/MomentFit.C for moments fitting.

### 3. Estimate Partial Wave Amplitudes from Moments
- Set partial wave amplitudes parametrization and path to the moments fit results in brufit/RunGivenBruMoments.C
- In brufit directory, run:
  
```bash
brufit Load.C RunGivenMoments.C
```

### 3. Analyze and Visualize Results

Use jupyter_notebook/Moments_IO_test.ipynb to analyze the results.
