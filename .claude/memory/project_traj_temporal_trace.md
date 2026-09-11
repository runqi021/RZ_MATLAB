---
name: project_traj_temporal_trace
description: Temporal-trace trajectory (delay-embed PCA + Isomap) for calcium vs breathing motion — two svd_breath_motion scripts
metadata: 
  node_type: memory
  type: project
  originSessionId: 2a7a33fd-578b-4aff-90c1-792d0024ee79
---

Two new standalone scripts in `svd_breath_motion/` (2026-06-06). **Current design (after iteration):**
NOT delay-embedding/SSA anymore, NOT dF/F anymore. Each sig ROI's trajectory = PCA / Isomap over its
**50 best-breath-correlated cycles** of the **spike probability** (Gaussian-smoothed deconvolved spike
train), embedding the 24 phase bins → a clean per-ROI phase LOOP.

- `Ventral_surface_traj_pca_svd_260606.m` — PCA; per-ROI 24-pt loop, top-3 ≈ 89–100%.
- `Ventral_surface_traj_isomap_svd_260606.m` — Isomap; per-ROI res3 ≈ 0.00–0.04.

## Trace = spike probability (not dF/F)
`useSpikes=true` (toggle): trace = `smooth_gauss(spike_train, round(spk_smooth_sec*fps))`,
`spk_smooth_sec=0.1`. The stored `spike_train` is sparse BINARY (`ca_spike_data.mat`, ~1.6% nnz; also
has `spike_amp`). Reason: dF/F carries a breath-coupled MOTION artifact (the FOV moves with breathing)
— likely the dF/F PC2 — so deconvolved spikes are the motion-free signal. Still shifted earlier by the
7-frame lead (`ca_lag_sec=0.2333`).

## Per-ROI = `roi_cycle_matrix` (shared helper)
Per ROI, pick the `popCycles` breath cycles whose lead-shifted spike-prob is best `corr` with breath
PC1, phase-resample each to `nPhaseBins=24` via `interp1` on cumulative `F.phi` → `cycMat
=[popCycles×24]`. Per-ROI fig = `pca_traj`/`isomap_traj`(cycMat) (cycles = ensemble dims, embed the 24
phase obs). Sparse ROIs: no-spike cycles have NaN corr → ranked LAST (relaxed guard: need only
`nCyc>=popCycles`, not 50 valid-corr), so all 10 ROIs render.

## Per-significant-ROI design
- Loops over **each significant ROI** (`PP.r >= confC`) from `coherence_polar_data.mat` — **10 sig
  ROIs across 7 FOVs** (confC=0.8269, 10/297; the `roi5_1400-1230-0` FOV has 4, others 1 each). One
  figure per ROI (`..._roi##`), folder load cached per FOV.
- **Only the CALCIUM trace is embedded** (single-panel figure). Breathing motion was REMOVED from
  the dimension reduction — its cycle is already defined by the 0/π breath phase, so embedding it is
  circular (user call). The breath phase is still used, only to COLOR the calcium loop.
- Calcium trace = that ROI's OWN dF/F column (was per-FOV mean over sig ROIs earlier; switched to
  per-ROI), shifted EARLIER by the data-driven lead `ca_lag_sec=0.2333` (7 fr @30 Hz,
  [[project_ventral_temporal_phase]]).
- Each trace z-scored; window `L = round(win_cycles · medianBreathPeriod)` (one cycle, from
  `median(diff(footIdx))`). Trajectory colored by `piecewise_phase_local` breath phase (onset=0 red,
  peak=π). Vglut2/1124 `+1` event shift applied in the loader.
- One figure per FOV (1×2 tiles) + `.mat`, into `traj_pca_svd_260606/` and `traj_isomap_svd_260606/`.
- Ran on **7 FOVs** with sig ROIs (mostly `pFN` fast-breath L≈14; two slow-breath `roi5` L=40/85).
  Breath loops are clean phase-colored rings; calcium manifolds are noisier but phase-organized.

## Population trajectory (pooled across the 10 sig ROIs)
- Added a `traj_{pca,isomap}_POPULATION` figure: pool all 10 sig ROIs into one trajectory. Since they
  are from 7 NON-simultaneous recordings, individual breaths can't truly align across FOVs.
- **Multi-cycle (current, `popCycles=50`)**: per ROI, select the `popCycles` breath cycles whose
  lead-shifted calcium is **best correlated with breath PC1** (`corr(calU,brU)` per cycle, take top
  N), resample each to `nPhaseBins=24` phase points (via `interp1` on cumulative `F.phi`), concatenate
  → `M=[10 ROIs × popCycles*24]`; embed the (cycle,phase) observations, color by within-cycle phase.
  All sig FOVs have ≥70 cycles (min `roi5_1400`=70; pFN 135–408) so popCycles up to ~70 fits without
  repeats. Within a FOV cycles are real; ACROSS FOVs the cycle pairing is **arbitrary** (accepted).
  Result: dense phase-organized cloud (50 cyc: 1200 obs, PCA top3≈58%, Isomap res3≈0.31).
  `popCycles` is a one-line knob at the top of each script (tried 10/20/50).
- Helpers `pca_traj` / `isomap_traj` (embed observations, z-score per ROI) added, copied from
  `phase_pca_svd`.

Reused `detect_session_fps`, and the `piecewise_phase_local` / `isomap_traj`/`knn_adj` / `draw_loop`
patterns copied from `Ventral_surface_phase_pca_svd_260606.m` (which stays untouched).
