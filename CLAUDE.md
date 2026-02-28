# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MATLAB-based 2-photon calcium imaging data analysis pipeline for neuroscience research. Processes multi-Z volumetric imaging data through dark correction, motion correction, AI cell segmentation (Cellpose), and fluorescence (dF/F) analysis. Integrates physiological signals (breathing, sensory stimuli) for event-triggered analysis.

## Running Scripts

There is no build system. Scripts are run directly in MATLAB:

```matlab
run run_stackZ_dFF.m                       % Single file: full pipeline
run Batch_BG_MC_SAM_260224.m              % Batch: dark→MC→segmentation
run Batch_breathing_triggered_avg_260226.m % Batch: breathing-triggered averaging
run Batch_sensory_STIM_260226.m            % Batch: stimulus-triggered averaging
run Breathing_triggered_avg_260226.m       % Single: breathing analysis
run select_ROI_dFF_plot_260117_v2.m        % Interactive ROI viewer
```

Parameters are hardcoded as user-editable variables at the top of each script — there are no config files.

## Requirements

- MATLAB R2021a+ with Image Processing, Statistics, and Parallel Computing toolboxes
- Python conda environment for Cellpose and OASIS:
  ```bash
  conda create -n cellpose-gpu -c conda-forge cellpose cuda-version=11.8
  conda activate cellpose-gpu
  conda install -c conda-forge oasis-deconv
  ```
- Python is invoked from MATLAB via `pyenv` (configured at top of pipeline scripts)

## Architecture

### Pipeline Stages (run_stackZ_dFF.m is the canonical reference)

```
ScanImage TIFF → Dark Correction → Motion Correction → Cellpose SAM → dF/F → Analysis
```

1. **Dark correction** (`minusDark/`): Estimates per-pixel noise floor, subtracts dark current, optionally masks dead pixels. Output: `*_minusDark.tif`.
2. **Motion correction** (`NoRMCorre-master/`): Piecewise rigid registration via FFT cross-correlation. Output: `*_minusDark_MC.tif` + `*_MCinfo.mat`.
3. **Cell segmentation** (`CellPose/`): Calls Python Cellpose CLI with SAM backend (GPU-accelerated). Extracts fluorescence traces F [T×N]. Output: `*_cpSAM_output.mat`.
4. **dF/F computation** (`2p_breathing_coherence/+helper/dFF_RZ.m`): Sliding median baseline normalization. Optional OASIS AR(1) spike deconvolution. Output: `*_dFF.mat` + `*_stackDFF.png`.

### Key Conventions

- All fluorescence data is stored as `[T × N]` matrices (time × ROIs/neurons).
- Output files follow the pattern: `[stem]_minusDark.tif`, `[stem]_minusDark_MC.tif`, `[stem]_cpSAM_output.mat`, `[stem]_dFF.mat`.
- Z-stack outputs go into `[filename]_ZSTACK_PIPELINE/z0/`, `z10/`, etc.
- Batch scripts skip pipeline stages if output files already exist (resumable processing).
- Functions return info structs (e.g., `darkInfo`, `MCinfo`) containing parameters and diagnostics.

### Helper Package: `2p_breathing_coherence/+helper/`

Called as `helper.functionName()`. Key functions:
- `dFF_RZ.m` — core dF/F computation with optional OASIS deconvolution
- `stackDFF.m` — stacked trace visualization (publication-quality plots)
- `oasis_deconv_and_dff_AR1.m` — OASIS spike inference via Python
- `phase_jackknife.m` — jackknife resampling for coherence statistics

### Third-Party Libraries

| Directory | Purpose |
|-----------|---------|
| `NoRMCorre-master/` | Non-rigid motion correction (piecewise rigid, FFT-based) |
| `chronux_2_12/` | Spectral analysis: multitaper FFT, coherence, spectrograms |
| `minusDark/` | Dark frame estimation and per-pixel masking |
| `DLC_tracking/` | DeepLabCut pose/breathing marker tracking |
| `Bassler/` | Basler secondary camera controller |

### Key Parameters (defaults used across scripts)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FPS` | 30 | Imaging frame rate (Hz) |
| `TossFrames` | 30 | Drop first N frames (stabilization) |
| `BaselineWinSec` | 20 | Sliding median baseline window (s) |
| `DarkSigma` | 3 | Sigma threshold for dark pixel masking |
| `Diameter` | 30 | Expected cell diameter (pixels) for Cellpose |
| `TauDecay` | 0.45 | GCaMP8s decay constant for OASIS (s) |
| `UseGPU` | true | CUDA acceleration for Cellpose |
