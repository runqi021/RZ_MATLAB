---
name: project_breathing_is_motion
description: "breathing_is_motion module — shift-based breathing analysis (GUI + analysis + phase), data alignment and Chronux gotchas"
metadata: 
  node_type: memory
  type: project
  originSessionId: f0ab7c50-89c3-4140-932f-26e3de194f8b
---

`breathing_is_motion/` (under RZ_MATLAB) = breathing analysis WITHOUT DLC, using the
2-photon motion-correction `|shift|` as a breathing surrogate (brain moves with each breath).
Built 2026-05-22 for sessions where running the full DLC pipeline isn't worth it.

Three scripts:
- `breathing_is_motion_gui.m` — clone of [[breathing_peak_gui]] (`breathing_peak_gui.m`) but
  fed from `*_preproc_MC_output.mat` instead of DLC CSV. Browses a master folder, builds the
  `|shift|` (µm) trace, interactive peak detect/edit, saves `*_breath_peak_data.mat` with the
  SAME fields as the DLC GUI (`breath`, `t_breath`, `insp_onset_idx`, `insp_onsets_train`,
  `findpeak_params`). `findpeak_params` adds `shift_source` ('MC1'|'MC1+MC2'), `axis`,
  `pixelSize_um`, `mc_file`. Shift Source panel: MC1 (default) or MC1+MC2 sum; axis magnitude
  (default) / signed x / y.
- `breathing_is_motion_analysis.m` — single folder (`folderPath` at top). Loads breath peaks +
  dFF + curated calcium spikes. Produces: |shift| µm trace + stacked dFF (20% dFF scale bar,
  Ca peaks overlaid red), inspiration-triggered dFF overlay, breath-train↔Ca-train coherence,
  breath spectrum. Saves figures + `shift_breath_master.mat`.
- `breathing_is_motion_phase.m` — phase modulation for one ROI from `shift_breath_master.mat`.
  Phase interpolated linearly 0→2π between consecutive inspiration onsets. 3 panels: dFF by
  breath cycle sorted by cycle duration dt (in time, in phase) + cosine fit `a+b*cos(phi-mu)`.
  μ = preferred/modulated phase. Spiking ROI → P(Ca spike|phase); non-spiking → mean dFF/phase.

Key non-obvious facts / gotchas:
- **Frame alignment**: use the FIRST MC pass `*_preproc_MC_output.mat` (NOT `*_MC_MC_output`).
  Its `mcOut.shifts` is a struct array, `shifts(k).shifts` = rigid `[y x]` pixels. nFrames ==
  cpSAM `F` rows == dFF rows (all post-toss, e.g. 2970 for a 3000-frame 30 fps run). No offset.
- `shift_um = sqrt(sum((yx*pixelSize_um).^2,2))`. pixelSize_um from `_meta.mat` (=1.7778/zoom).
- **Calcium train = curated spikes** from `calcium_spike_gui.m` → `ca_spike_data.mat`
  (`roi_spikes(i).spike_idx`, top-level `ifSpike` logical per ROI). islocalmax is FALLBACK only.
- dFF loaded from `*_dFF.mat` (the trace spikes were detected on), not recomputed, so indices align.
- **Chronux `coherencypt`**: pass point times as structs (`d.times = ...`) or it length-checks
  them as continuous data ("inconsistent dimensions"). Must request ≥9 outputs
  `[C,phi,~,~,~,f,~,confC,~]` or confC isn't computed. Pass full `t_img` grid as 5th arg.
- **Phase modulation null**: circularly shifting the *tuning curve* can't change its cosine
  amplitude (rotation preserves 1st Fourier magnitude) → trivial p=1. Correct null circularly
  shifts the dFF *time series* against phase (or random spike phases for the spike case).
- `|shift|` magnitude rectifies breathing motion (frequency-doubles); signed x/y axis is cleaner
  for phase but user chose magnitude.

Related: [[reference_cosine_phase_tuning]] (the planned master_cosine_fit_phase.m was never built).
