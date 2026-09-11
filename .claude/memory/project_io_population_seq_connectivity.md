---
name: project_io_population_seq_connectivity
description: IO population ensemble-event viewer (IO_population_seq_connectivity_260607.m) — spike-prob histogram + per-ensemble latency map on avg proj
metadata: 
  node_type: memory
  type: project
  originSessionId: 3c68632d-e187-42f7-b413-da2dc5adcb4b
---

`IO_population_activity/IO_population_seq_connectivity_260607.m` (2026-06-07) — standalone
inferior-olive population tool. NEW file (NOT a rework of `breathing_population_manifold.m`, which was
left untouched per [[feedback_no_inplace_edits]]; that script is pure dF/F PCA, no spikes).

## Signal = GLOBAL dF/F + calcium threshold, NO spike sorting (2026-06-07, latest)
User realized the `ca_spike_data` "spikes" are just raw-dF/F threshold peaks (spike_params method
`'Raw dFF'`, oasis_run=0) → spike sorting adds nothing here. Switched the whole analysis to dF/F:
`Z=zscore(dFF,0,1)`, `seqSig=smooth(Z, seq_smooth_sec=0.1s)`. (a) Detection: `P_pop=mean(seqSig,2)`,
peaks > mean+2SD. (b) Per-cell activation timing = first **rising-edge crossing** of z-dF/F through
`caThreshSD=2` SD in the window (uniform threshold since z-scored; captures onset, and uniform GCaMP lag
preserves order — user "we can assume they all have the same lag"). (c) Raster = z-dF/F threshold
crossings. `ca_spike_data` no longer consumed (loader still reads spikeProb but it's unused). Single
knob `caThreshSD`. Verdict from comparison: spike-based vs dF/F-threshold gave ~identical results here
(active sets identical, 99% of crossings same frame) because the spikes were just smoothed-dF/F peaks.
The PRIOR spike-based design (first detected spike from `roi_spikes.spike_train`) is superseded.

## CURRENT SCOPE (drastically reduced 2026-06-07 — user "take out all other plots")
The script was first built big (PCA, connectivity, modules, sequence/latency stats, split-half,
lead-lag, wave fit, time tuning, shuffle) then **stripped to 3 things only**. Filename still says
`seq_connectivity` but content is now an ensemble-event viewer (rename not requested). outTag =
`IO_ensemble_events_260607`. Produces **2 figures + results.mat**:
1. **Exact population activity** (Fig `01_population_activity_raster`, 4-row subplot): TOP = population
   spike probability `P_spk(t)=mean(seqSig,2)` filled area + threshold (mean+2SD or pct) + ensemble-colored
   event triangles; BOTTOM (rows 2-4) = **full spike raster** of all ROIs (exact spike times from binary
   `DATA.spikeProb` as dots; falls back to z-dFF imagesc if no spikes) with faint ensemble-colored event
   xlines, x-linked to the top. User wanted to "see what is really going on", not just the mean.
2. **Ensemble events** — peaks in `P_spk` (MinPeakDist=1s) within a ±`winSec` window (default ±0.3 s),
   drop events with <`minActiveNeurons` active cells, then **k-means** on per-event across-ROI activity
   vectors `EvVec[nEv×N]` (z-scored per event). **k is data-driven**: `nEnsembles=[]` → silhouette
   auto-select via `evalclusters` over `2:maxEns` (maxEns=8); set a number to force. (Was hard-coded 4.)
3. **Per-ensemble activation SEQUENCE on avg projection** — Fig `02_ensemble_activation_sequence`, one
   panel per ensemble: active cells ordered by mean **first-activation time** and connected by a line
   through their centroids over the avg proj; segments+nodes colored by activation time (parula,
   clim=[-winSec winSec]), open circle = first-activated cell. (2026-06-07 user: replaced the latency
   HEATMAP with connecting LINES, and switched COM→**first-spike**.) Helper `draw_sequence_on_fov`; old
   `draw_latency_on_fov`/`roi_boundaries` removed. Tangled lines here = near-synchronous, not a clean wave.
   - **Activation time = FIRST DETECTED SPIKE** (2026-06-07 user "are we not using detected spikes?"):
     read straight from the binary `DATA.spikeProb` (roi_spikes.spike_train) — `find(Wspk(:,i)>0,1)` in
     the window — NOT a crossing on the smoothed trace. Active = ≥1 detected spike in window (no SD
     threshold). Smoothed `seqSig` is still used for P_spk + EvVec clustering only. Fallback (no spike
     data): threshold crossing on seqSig at `roiThr=mean+actThreshSD·SD`. Feeds `minActiveNeurons`.
   - **eventAmp** = `mean(EvVec,2)` saved separately (clustering on `zscore(EvVec,0,2)` is amplitude-blind
     by design); `results.ens_summary` = per-ensemble nEvents/meanAmp/meanActive; titles show amp=.

## 2P epoch stacks of ensemble activity — 16-bit TIFF (added 2026-06-07)
`saveEpochTiffs=true` writes **one multipage 16-bit TIFF per ensemble** into
`<outDir>/ensemble_2p_epochs/` (`ensemble#_2p_epochs.tif`). Source = `*_MC_MC.tif` (uint16, already
trimmed to T=2970, SAME as F/dFF → **offset=0**, event frames map directly; code still computes
`offset=max(0,nMovie-T)` defensively). Each stack = that ensemble's event peaks cut to
`[-vidPreSec..+vidPostSec]` (default ±0.5 s), concatenated with `vidGapFrames` black separators.
**Frames are RAW uint16, ORIGINAL dynamic range — NO scaling, NO clip, NO gamma** (LZW lossless via
`imwrite`). Open in ImageJ/Fiji; playback speed is a viewer setting (no baked-in fps).
`tiffreadVolume` loads whole movie (512×256×2970 ~778MB; fine).
**Iteration history (all 2026-06-07, same session):** started MPEG-4 mp4 2×-slow with colored border →
user "terrible quality" → uncompressed Grayscale AVI, border removed → user "u gamma the video?" →
gamma removed (linear) → user "save original dynamic range, why 256?" (AVI is 8-bit) → user "i mean tif"
→ final = 16-bit LZW TIFF, untouched. See [[feedback_never_gamma_video]] (NEVER gamma video frames).
User-tuned detection params: evMinDist=0.5s, winSec=0.25, minActiveNeurons=10 → fewer/cleaner events.

## Reused infra (loader + display helpers kept)
- `seqSig` = deconvolved spike prob smoothed `seq_smooth_sec=0.1 s` (else positive z-dFF fallback).
- Loader supports session FOLDER (`*cpSAM_output.mat`→F→`helper.dFF_RZ`; `*_AVG_for_CP.tif` avg proj;
  `detect_session_fps`) or direct `.mat` (`dFF`+optional `spikeProb`/`maskL`).
- **`ca_spike_data.mat` here = per-ROI struct array `roi_spikes(1×N).spike_train`** (NOT a T×N matrix);
  loader does `[SP.roi_spikes.spike_train]`. Sparse BINARY.
- Avg-proj overlay trick: show grayscale as TRUECOLOR `repmat(img,[1 1 3])` so axes colormap stays free
  for the latency colorbar; ROI fills via `val2rgb` patches. clip [0.5 99.9] + gamma 0.6.
- pixelSize_override param (e.g. 0.4444551) else metadata/zoom.

## Verified
Runs on `Ventral_surface_summary/Vglut2/1124/IO/roi3_-1000-300-60...` (T=2970, N=127, 30 Hz): 47 events,
k=4 ensembles (counts 13/17/5/12), both figures render correctly.
