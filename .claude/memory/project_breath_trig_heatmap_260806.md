---
name: project-breath-trig-heatmap-260806
description: "analysis_260806 breath-triggered dF/F heatmap scripts (single FOV, pooled session, cross-genotype) and their design rules"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2c1b7eed-e822-44e4-86fb-1002c555d522
  modified: 2026-08-08T07:04:50.948Z
---

`analysis_260806\` — deliberately simple breath-triggered dF/F figures, built 2026-08-06 for
`C:\260806_sst-soma-g8s\phys` (SST soma-G8s). Three standalone scripts, same panel design:

- `breath_trig_heatmap_260806.m` — ONE recording folder (edit `folderPath` at top)
- `breath_trig_heatmap_pooled_260806.m` — all recordings in one session folder
- `breath_trig_heatmap_genotypes_260806.m` — the 260806 session + `D:\Ventral_surface_summary`
  Vglut2/Vgat/Sst/ChAT/Sert (63 recordings, 1913 ROIs, 156871 single cycles)

Layout: 2 rows, row 1 = 1/4 the height of row 2. Row 1 = breath trace | breath power.
Row 2 = dF/F traces (or single-cycle heatmap) | mean dF/F heatmap | z heatmap | z power.
Power panels always in the RIGHT column.

Settled choices (asked for explicitly, don't re-litigate):
- **NO sorting** of ROIs, ever — native cpSAM order, every ROI included, no significance test
- z = per-ROI over the WHOLE trace, computed on dF/F (not raw F), then trigger-averaged
- clim = `prctile(...,[0.5 99.5])`; dF/F panel = `gray(256)`, z panel = `parula(256)`
- **windows in breath CYCLES, not seconds** — `ibi = median(diff(trig))/fps` per recording;
  top trace = 10 cycles, each panel below = 10/3 cycles, so the trace spans the 3 panels
- multitaper knob is `W_hz` (half-bandwidth, default 0.25) -> `TW = W_hz*T_sec`, `K = 2TW-1`,
  `err = 0`. NOT a bare TW: this keeps smoothing constant in Hz across recording lengths
- spectra `fpass = [0.1 14]`, log-log, NOTHING normalised (PC1 and dF/F are both z-scored,
  so both are already z^2/Hz — normalising throws that away)
- dF/F trace stack is on ONE common scale with a 5% dF/F bar; per-trace auto-scaling was
  explicitly rejected
- `clear; clc;` but NOT `close all` — figures accumulate so FOVs can be compared

Cross-session pooling: breath rate (IBI 0.5-1.5 s) and frame rate (30/42/47 fps) both vary, so
each PETH is built on its own axis then resampled onto a common axis in CYCLES, labelled in
seconds at the pooled median IBI. Pooling in raw seconds would smear the +/-1 cycle peaks.

Data gotcha on 260806: `breath_pc1.mat` came from the fixed-metric route (`src_metric='fb'`),
which truncated every run to `T=3000` (the shortest run), so the breath trace covers only the
first 100 s of the 200 s recordings — `fb_full` in `phys\breath_fixedmetric.mat` holds all 6000
frames if the rest is ever needed. See [[project_breath_fixedmetric]] and [[project_breath_svd_pc1]].

Observation, NOT yet tested: a one-frame dF/F DROP at exactly t=0, simultaneous in all ROIs,
present in raw F, appears in most recordings. Sign/timing identical across FOVs and frame rates
is the signature of breathing-locked motion, not calcium — needs the MC-shift or core/ring test
before anyone calls it physiology. See [[project_ca_spike_detector_260722]].
