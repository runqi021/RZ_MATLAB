---
name: project_ca_spike_detector_260722
description: "Calibrated spike detector for Ventral_surface_summary — progress tracker: kernel fit, OASIS benchmark, QC gallery, open next steps"
metadata: 
  node_type: memory
  type: project
  originSessionId: a1713ea0-2ee8-4943-8d90-ba971a8be969
  modified: 2026-07-23T03:49:46.077Z
---

**Goal (started 2026-07-22):** reliable dF/F spike detection that finds real Ca events and not motion/noise, calibrated against the user's supervised labels in `C:\Users\Admin\Desktop\Ventral_surface_summary`. Two-stage design: **temporal model (deconvolution) = detector**, **spatial test = artifact veto**. They are orthogonal; deconvolution alone cannot reject artifacts because it only sees one trace.

**Ground truth:** 25 sessions with `ca_spike_data.mat` (from [[project_dffQC_pipeline]]'s sibling `calcium_spike_gui.m`), 868 ROIs, 5337 hand-clicked events, 4208 at >=4 sigma. Only **362 ROIs carry >=1 label**. Groups: ChAT / Sst / Vgat / Vglut2 (`Vglut2_test` folded into Vglut2).

## Scripts written (all NEW standalone, nothing existing modified — [[feedback_no_inplace_edits]])
- `ca_spike_kernel_fit_260722.m` — step 1: measures the calcium kernel from isolated >=4σ labels. Out: `Ventral_surface_summary\_spike_kernel_260722\` (`spike_kernel_fit.mat` has `K` per group and `E` per event with `.session`/`.tau`).
- `ca_spike_deconv_bench_260722.m` — step 2: benchmarks raw dF/F vs matched filter (4τ and 1.5τ) vs OASIS. Out: `_spike_deconv_bench_260722\`.
- `ca_spike_detect_gallery_260722.m` — QC gallery, 100 random curated ROIs, 5 pages, PNG + vector PDF. Out: `_spike_gallery_260722\`.
- `ca_oasis_run.py` — OASIS as a subprocess (see [[reference_oasis_on_this_machine]]).

## Measured kernel (step 1, 1804 isolated events)
| group | n | rise 10-90% | tau | g@30fps |
|---|---|---|---|---|
| ChAT | 270 | 67 ms | 0.49 s | 0.934 |
| Sst | 29 | 133 ms | 0.64 s | 0.949 |
| Vgat | 25 | 100 ms | 0.48 s | 0.933 |
| Vglut2 | 1480 | 133 ms | 0.91 s | 0.964 |

`calcium_spike_gui.m:761` hardcodes `g = 0.93` (tau 0.46 s): fine for ChAT/Vgat, **wrong for Vglut2 (~2x)**. **Per-session tau varies as much as between groups** (Vglut2 0.36 -> 1.20; ChAT 0.47 -> 0.80) and per-session g is what decides whether OASIS works on a session — group-level g is what made OASIS fail on Vglut2/0224 (true tau 0.36 s). Rise is 2-4 frames for every group but ChAT, so AR(1)'s instantaneous-rise assumption is violated. Sst/Vgat tau rest on only 25-29 events — provisional.

## Hard-won gotchas (each cost a failed run)
- **Restrict scoring to ROIs with >=1 label.** `ifSpike=false` cannot distinguish "inspected and silent" from "never inspected". Including uncurated ROIs let raw dF/F fire at 117x the label rate while appearing clean.
- **Alignment**: OASIS marks the transient ONSET, labels sit on the PEAK. Measured shift = **+4 frames**. The biexponential kernel argmax predicts 5-9 frames and OVER-corrects — use the measured value, do not add to the model value.
- **Match tolerance ±4 frames is too strict** (catches 0.67 of labels); use **±6** (0.78). Rise is 2-4 frames and the GUI snaps clicks to a local max within ±3 frames.
- Compare methods at matched false-alarm rate or by F1, **never at equal sigma threshold** — each statistic has its own noise scale.
- `conv(...,'same')` for a matched filter recentres by floor(K/2) = 30-55 frames here. Use explicit 'full' indexing.
- Sparse spike trains make MAD degenerate (median 0) — scale OASIS by `sn`, not MAD.
- The `-0.2/+2.5 s` post-label exclusion mask eats nearly all of a dense session's trace, so `det/min` and any per-minute rate is nonsense there. F1/precision/recall are unaffected.

## Results so far
- Step 2 benchmark (curated ROIs, group g, **pre-alignment-fix so all numbers are depressed**): best F1 raw 0.45, mf 0.53, mfs 0.52, **OASIS 0.66**. Per group OASIS: ChAT 0.95, Vgat 0.72, Vglut2 0.54, Sst 0.47.
- QC gallery with per-session g + measured lag + ±6 tol. **Threshold sweep on the 100-ROI sample settles the operating point:**

| th (*sn) | precision | recall | F1 |
|---|---|---|---|
| **1.50** | **0.87** | **0.76** | **0.81** |
| 2.00 | 0.96 | 0.69 | 0.80 |
| 3.50 | 0.98 | 0.44 | 0.61 |
| 4.00 | 0.98 | 0.39 | 0.56 |

  **Use 1.5*sn.** Below 1.5 precision falls faster than recall rises. The old 2.0 came from the stale broken-alignment benchmark and was far too conservative (P 0.96 with 11 false positives in 100 ROIs = obviously over-tight).
- **Failure mode is specific**: dense, low-amplitude ROIs (e.g. `Vglut2/pFN_roi4` ROI7, ~2 Hz oscillation at 0.2-0.3 dF/F, 19 of 51 found). Sparse clean ROIs are near-perfect.

## Motion (the original question)
`sx`/`sy` from `*_MCinfo.mat` do **NOT** predict labelled spikes: motion-speed AUC median **0.500** (range 0.41-0.65) across 25 sessions. So an xy-shift-based motion veto buys nothing. Residual artifacts must be z-motion / vessel / flow, which need a **spatial** test from the MC TIFFs (core vs ring, and projection on the image gradient, since a translation changes F by grad(F).delta). Only `Vglut2/1124/cell1` (0.65) and `Vglut2/0224/cell1` (0.62) show real motion-locking. NOTE: breathing drives both motion and the activity of interest, so motion must be a covariate, never a hard veto.

## Next
1. Re-sweep threshold on the gallery sample (running) -> real operating point; expect well below 2.0.
2. Re-run step 2 benchmark with per-session g + lag +4 + tol 6 for honest numbers.
3. Then step 3: spatial veto (needs one pass over `*_preproc_MC*.tif`).

## Caveats that bound everything
- "Precision" = agreement with the user's clicks. A real unlabelled event scores as a false positive.
- 12.8% of labels are <3 sigma, 5.2% <2 sigma — the labels themselves have a soft low-confidence tail, which caps achievable F1.
- `ChAT/0521/cell1` and `ChAT/0523/cell1` contributed ZERO usable events to the kernel fit and have the worst negative-going rates (1.80, 1.12) — their labels look dominated by low-confidence detections.
