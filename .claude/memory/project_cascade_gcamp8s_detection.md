---
name: project_cascade_gcamp8s_detection
description: "CASCADE deep-net spike inference for GCaMP8s somata — the label-free detector for Ventral_surface_summary; setup, env routing, and the pitfalls that produced wrong numbers"
metadata: 
  node_type: memory
  type: project
  originSessionId: a1713ea0-2ee8-4943-8d90-ba971a8be969
  modified: 2026-07-23T06:06:49.336Z
---

**CASCADE is the automatic GCaMP8s event detector** (Rupprecht et al., Nat Neurosci 2021). A deep net trained on ground truth from SIMULTANEOUS ephys + 2P calcium in other animals, so it needs **no labels from this dataset** and no threshold tuning — it measures each trace's noise and routes it to the ensemble member trained at that level. Set up 2026-07-22 for [[project_ca_spike_detector_260722]].

## Why the Torch fork, and which env
Published CASCADE pins `python<3.8` + `tensorflow==2.3` — unusable here. Use **CascadeTorch** (`github.com/PTRRupprecht/CascadeTorch`, torch port released 2026-02-11, same converted weights, identical output). Cloned to `%TEMP%\CascadeTorch`.

Env routing, chosen to install NOTHING and touch no working env:
- **dlc310** runs CascadeTorch: torch 2.9 + CUDA, and its LAPACK **works** (unlike the `oasis`/`cellpose-gpu` envs — see [[reference_oasis_on_this_machine]]). Models are native `.pth`, so no h5py needed.
- **cellpose-gpu** does the HDF5 reading (only env with h5py). `*_dFF.mat` are v7.3, which `scipy.io` cannot read. Do NOT add h5py to dlc310 — a second HDF5 binary risks DeepLabCut.

## Scripts (repo root)
- `cascade_export_dff.py` (cellpose-gpu) — v7.3 → `.npz`. Note h5py returns MATLAB arrays transposed, so dFF arrives as [N x T] = exactly the [neurons x time] CASCADE wants.
- `cascade_run_gc8s.py` (dlc310) — inference, model chosen per cell type.
- `cascade_fix_discrete.py` (dlc310) — recompute discrete spikes from saved rates.
- `cascade_plot_gallery.py` (dlc310) — 100-ROI gallery, 5 pages, PNG + PDF.
- `cascade_vs_labels.py` (cellpose-gpu) — external validation vs the hand clicks.

Outputs under `Ventral_surface_summary\_cascade_260722\{traces,predictions,gallery}`.

## Model: use `GC8s_EXC_30Hz_smoothing25ms_high_noise` for ALL groups here
Including Sst and Vgat, despite those being GABAergic.

**`Interneurons_GC8+_30Hz` was tried on Sst/Vgat and REJECTED.** It looked like a fix (rate jumped from 3-10 to ~246 spikes/min) but the output is a flat prior, not a measurement: all six sessions landed within 9% of each other (243-265 spikes/min), CV(rate) 0.20-0.99, correlation with the actual dF/F ~0 (-0.05 to 0.16), and it fires for **56-92% of the recording**. The EXC model on the same traces gives sparse, onset-locked bursts confined to <1% of frames.

**Diagnose this failure by CV + duty cycle, NOT by corr(rate, dF/F)** — a sparse spike train cannot correlate instantaneously with a slow decaying integral, so low correlation is normal even when the detection is good.

## THE noise-level caveat that bounds all of this
CASCADE's nu = median|diff(dF/F)| / sqrt(fps) * 100 (it scales to percent internally, so **ratio input is correct**). On this data nu is 0.35-2.35, **median 1.24**, and **90.3% of ROIs fall BELOW the model's lowest trained level (nu = 2)**. So nearly every trace routes to the same nu=2 network and the per-trace noise matching does not actually happen. Extrapolation is in the benign direction (data cleaner than training) but the calibration is not guaranteed. Sst/Vgat are the cleanest of all (nu 0.35-0.69), furthest below the floor, which is likely why the interneuron model collapsed to its prior there. Cascade issue #61 covers requesting models at other noise levels.

## Pitfalls that each produced silently wrong numbers
- `infer_discrete_spikes` returns **(approximations, spikes)** — the SECOND is spike times. Taking the first gives an array the length of the trace, which looks like one spike per frame (324k "spikes" over 100 ROIs).
- Discrete spike positions come back as **floats**; round before using as indices.
- **MATLAB v7.3 stores an empty array as a dataset holding its dimensions**, tagged `MATLAB_empty=1`. Reading naively yields 2 fake labels on every uncurated ROI — inflated 5337 real labels to 6096 and destroyed the precision estimate. Always check the attribute.
- CASCADE expects dF/F as a **ratio, not percent**; percent input silently routes to the wrong noise model. `cascade_run_gc8s.py` guards on p99 > 20.
- **CASCADE infers ACTION POTENTIALS; the hand labels mark CALCIUM EVENTS.** One transient holds several spikes, so a 1:1 comparison charges a false positive per extra spike. Burst-merging at 0.35 s lifts precision 0.46→0.54 but costs recall 0.57→0.34 (it merges genuinely separate events too) — neither view is definitive.

## Agreement with the hand labels (EXC model everywhere, i.e. before the interneuron fix)
Pooled P 0.46 / R 0.57 / F1 0.51; per-ROI median F1 0.55. Far below the 0.96 the hand-tuned OASIS detector reached — **but that detector was tuned to these labels and CASCADE has never seen them.**

**The disagreement is diagnostic, not just error.** CASCADE finds 4 events where 219 were labelled in `ChAT/0523/cell1`, and 9 vs 89 in `ChAT/0522/IO` — the SAME two sessions that show sign-flipped FDR ~1.0 and contributed ZERO usable events to the supervised kernel fit. Three independent lines now say those labels are mostly noise.

## Output shape
Primary output is a **spike rate in spikes/s**, calibrated in absolute units — for rhythmic breathing-locked ROIs that is more useful than discrete onsets, since it gives a continuous firing-rate estimate to regress on breath phase. Discrete times are secondary. In the gallery the rate tracks the **rising phase** of each transient (spikes drive the rise; the decay is indicator kinetics), and silent ROIs get ~0 with no tuning.
