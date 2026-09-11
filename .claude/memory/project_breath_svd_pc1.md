---
name: project_breath_svd_pc1
description: "breath_svd_pc1.m — DLC-free breathing trace from 2P-triggered Basler video via SVD PC1, fps auto-read from calcium metadata"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2d9a96e6-14a8-4f27-8d60-c0c91c42e3bd
  modified: 2026-07-29T19:10:43.329Z
---

**MOVED 2026-07-29: no longer at the repo root.** The whole breath front end now lives in `analysis_260727\breath_svd\` (with `legacy_dlc\` for the DLC-dot GUIs and `diagnostics\` for the June method-development scripts). It has its own README. `run_analysis_260727.m` and `analysis_preflight_260727.m` both addpath it. Moved files needed `repoRoot = fileparts(fileparts(here))` — two levels, not one — or `detect_session_fps`/Chronux/`+helper` are not found; `orofacial_crop_extract.py` deliberately stayed at the repo root because the orofacial analysis shares it. Verified after the move: 16 files parse, repoRoot resolves to the repo root, all shared deps reachable.

**Which script is the current per-video breath trace: `breath_svd_pc1.m`.** Every analysed dataset was built with it and `analysis_preflight_260727.m` names it. `breath_fixedmetric_*` (see [[project_breath_fixedmetric]]) is a SECOND front end writing the SAME `breath_pc1.mat`, for when breath *amplitude* must be comparable across runs — one ROI per session, frozen basis, `disp` in pixels. Its `pv` metric IS per-video PC1, kept only as the benchmark. Per-video / different views / timing-only → `breath_svd_pc1`. Many runs of the same view + amplitude → fixedmetric.

`breath_svd_pc1.m` — extract a breathing trace from the orofacial Basler video **without DLC**, using SVD of facial motion (Facemap-style). Created 2026-06-05 for the Ventral_surface_summary dataset.

**Idea:** the breath cam is 2P-frame-triggered, so its true fps = calcium imaging fps (the AVI header is bogus). SVD of mean-subtracted intensity → **PC1 = dominant motion = breathing waveform** (signed; replaces a DLC nostril marker). On the ChAT 0522 test, breathing held ~66% of variance → PC1 is clean.

**AVI naming (updated 2026-07-22):** discovery now uses a pattern list `AVI_PATTERNS = {'Basler_*.avi','cam*.avi'}` (local helper `find_avis`), not a single `Basler_*.avi` glob. Legacy Pylon-Viewer name = `Basler_*.avi`; new dual-cam GUI (basler_dual_acq / [[project_basler_acquisition_gui]]) writes `cam1_YYYYMMDD_HHMMSS_runNNN.avi`. In `260721_Sert_soma_G8s\phys` only the first folder is legacy Basler, the other 26 are cam1 — the old single-glob would have silently processed 1/27. One video per recording folder; first sorted match used. **No `timestamps.csv` in these folders** — fine, `orofacial_crop_extract.py` prioritizes the `--fps` override (imaging fps) over timestamps, so ts absence is harmless.

**Per folder under ROOT_DIR (auto-finds AVIs via `AVI_PATTERNS`, or set FOLDERS subset):**
1. `fps = detect_session_fps(folder)` — reuses the repo fps detector (reads ScanImage `scanFrameRate` from `_meta.mat`/TIFF/`_dFF.mat`; returns 30 for these, scanFrameRate_raw≈29.99, numSlices=1). See [[feedback_detect_session_fps]].
2. interactive crop via `orofacial_crop_extract.py --fps <fps>` (no `--fullframe`); cube → `breath_crop.mat`.
3. SVD (mean-sub intensity, auto spatial-bin to <1.5GB) → top 10 PCs.
4. PC1 figure (eigenimage + trace + PSD) + saves `breath_pc1.mat` (`breathTrace`=PC1, t, fps, eigImg, sv, varExp) + `breath_pc1.png`. Summary table at end.

**Diff map in peak GUI (2026-06-06):** under the trace, `breathing_peak_gui_pc1.m` shows the PC1 eigenimage AND a `top100 - bottom100` diff map (inspiration - baseline; frames ranked by oriented breath signal) side by side, BOTH parula. The diff map is the full-res `diffImg` cached in `breath_pc1.mat` (`breath_svd_pc1.m` now computes diffImg+camImg from the cube and saves them; `breath_cache_diffmap.m` backfills old files; `breath_diff_image.m` is the standalone viewer). GUI `compute_diff_map` prefers cached `diffImg`, else reconstructs a coarser one from U*sv*(meanV_top-meanV_bot). eigImg is binned (SVD res), diffImg is full crop res — same FOV. The SVD eigenimage was always CACHED (not real-time).

**Peak/trough GUIs on PC1 (added 2026-06-05):** new standalone copies of `breathing_peak_gui.m` / `breathing_trough_gui.m` that read PC1 instead of DLC (originals untouched, per [[feedback_no_inplace_edits]]):
- `breathing_peak_gui_pc1.m` — sessions = folders with `breath_pc1.mat`; signal = `breathTrace` detrended+z-scored (Invert toggle since SVD sign arbitrary); shows PC1 eigenimage panel to confirm snout; same `find_peaks_local` + click-edit; saves `breath_peak_pc1.mat` with the SAME field names as the DLC peak output (insp_onset_idx, insp_onsets_t/amp/train, breath, t_breath, findpeak_params{camera_fps,inverted,source='svd_pc1'}).
- `breathing_trough_gui_pc1.m` — near-identical to the DLC trough GUI (loads peak waveform verbatim, `find_feet_lineFit`); only the session scan changed → finds `breath_peak_pc1.mat`; saves `breath_insp_start_pc1.mat`.
- Pipeline: breath_svd_pc1.m → breath_pc1.mat → peak_gui_pc1 → breath_peak_pc1.mat → trough_gui_pc1 → breath_insp_start_pc1.mat. Same-field outputs so existing downstream breathing analyses ingest them unchanged.

**DLC-vs-PC1 validation + heartbeat (2026-06-05, `breath_dlc_vs_pc1.m`, 17 Vglut2+ChAT traces):**
- PC1 vs DLC: **median |r|=0.97** (0.81–0.99), zero lag. PC1 is a validated DLC replacement and is as clean or CLEANER (SVD pushes heartbeat into higher PCs — PC2 in fast breathers).
- **Heartbeat ~4.1 Hz (3.9–4.6), steady across all sessions** (independent of breath rate → it's the heart, not a harmonic). Mouse HR ~250 bpm under anesthesia.
- Contamination is breath-rate-dependent: **slow breathers (≤1.2 Hz) = 1–10% heart (trivially clean); fast breathers (~2.1 Hz) = 18–28% heart in BOTH DLC and PC1** because heartbeat (4.2) ≈ 2× breath rate (overlaps the 2nd harmonic, SVD can't fully separate).
- Outputs in `Ventral_surface_summary\_dlc_vs_pc1\` (per-folder PNG + SUMMARY.png + summary.csv).
- **Recommended peak/trough detection:** PC1 → adaptive breathing **low-pass at ~1.8×breath-rate** (kills 4 Hz heartbeat; THE key fix for noisy traces) → events via Hilbert phase (peak=crest, insp-onset=rising zero-cross; one/cycle) OR existing prominence-peak + line-fit-foot run on the low-passed signal. See [[project_phase_time_analysis]].

**Finalized onset (inspiration-start) detector (2026-06-05) — `find_feet_deriv` in `breathing_trough_gui_pc1.m`:** anchored to the user's already-detected peaks (search ONLY each inter-peak gap). Per peak: take the steepest rising sample, walk back to the baseline foot, CONTINUING through any local min still above `baseline+20%*amp` (so doublet notches don't trap it), bounded by the previous peak. RAW signal, NO amplitude filtering (filtering rounds the sharp foot = the feature). Signal must be inspiration-UP (auto-invert by skewness). This is the derivative sign-flip / local-min idea the user already uses; the baseline-walkback is the only addition (for notches/doublets).
- **Doublets** (2 stacked inspiratory pulses, rare — only 3 across 17 traces) resolve by PEAK count: user marks 2 peaks → 2 feet (2nd=notch=2nd inspiration); marks 1 peak → 1 foot at true baseline (walkback skips notch). No special doublet logic — it falls out of the inter-peak-gap bounding.
- **Truncated-first-inspiration guard:** if a recording starts MID-inspiration the true onset is off-screen — don't fabricate it. Test: first peak gets an onset only if `min(signal[1:firstPeak]) <= baseline+20%amp` (trace actually returned to baseline before it). Else NaN/drop. Don't toss ALL firsts — only 2/17 sessions (ChAT roi3_R, roi5_7x) are truncated; the other 15 start at rest and keep their real first onset. In `find_feet_deriv`.

**Duration guard [0.05, 1.0] s**: upper 1.0 verified safe (0% exceed); lower must be 0.05 NOT 0.1 (13.6% of inspirations <0.1s, all fast Vglut2 pFN 2Hz breathers where inspiration is genuinely 1-3 frames at 30fps). Fact-checked on 2379 DLC inspirations.
- DLC peak files often have bogus `camera_fps=50.7` → compute durations from frame-idx ÷ true fps (`detect_session_fps`=30), NOT stored `_t`.
- `breath_detect_all.m` = run detector on all traces → per-session full-trace PNGs in `Ventral_surface_summary\_detect_derivative\`. `breath_dlc_vs_pc1.m` / `breath_onset_methods.m` / `breath_doublet_diag.m` = analysis/diagnostics.
- DLC onsets+peaks are the **working ground truth** (user used them for all prior analysis). Heartbeat ~4.1 Hz contaminates fast breathers (see [[project_orofacial_crop_dff_analysis]] sibling note).

**`svd_breath_motion/` folder (2026-06-06):** home for downstream analyses driven by the SVD breath signal. First: `Ventral_surface_coherence_polar_svd_260606.m` = SVD-breath version of `Ventral_surface_coherence_polar_260528.m` (breath-phase × Ca-spike coherence polar). Only change: reads `breath_peak_pc1.mat` (insp_onset_idx=peaks) + `breath_insp_start_pc1.mat` (insp_start_idx=feet) + BP.breath(=PC1 waveform) instead of DLC `*breath_peak_data.mat`/`*breath_insp_start_data.mat`; adds repoRoot to path (subfolder, [[feedback_path_setup]]); out→`coherence_polar_svd_260606/`. Same field names so it's a drop-in swap. Verified: 23 recordings, 297 ROIs, 10 sig. Per-group N matches DLC version (ROI inclusion gated by Ca spikes not breath) → directly comparable to DLC `coherence_polar_260528`. Whole Ventral_surface phase/coherence SERIES ported (10 scripts, all `*_svd_260606.m` in `svd_breath_motion/`): coherence_polar, temporal_phase (single/batch/population/population_groups), phase_pca (+combined), phase_polar_allROI, spike_phase_polar, coherence_polar_phase. **Port recipe (sed):** breath globs `'*DLC*breath_peak_data.mat'`/`'*breath_insp_start_data.mat'` (and bare `'breath_*_data.mat'` for find_data_file) → `'breath_peak_pc1.mat'`/`'breath_insp_start_pc1.mat'`; `addpath(scriptDir)` → `repoRoot=fileparts(scriptDir); addpath(repoRoot)` + explicit 2p_breathing_coherence + genpath(chronux_2_12) (5 of them only had `addpath(scriptDir)`, assuming chronux pre-loaded — must add it); `coherence_polar_260528`→`coherence_polar_svd_260606` (sig-list source); outdir/output stamps → `_svd_260606`. Outputs go to `*_svd_*` dirs so DLC results aren't overwritten → side-by-side comparable. Pre-existing cosmetic warnings (true()-shape, unreachable-stmt) are in the originals, not from the swap.

**Key facts / caveats:**
- Signal = mean-sub intensity, NOT dF/F (dF/F divides by ~0 baselines on dark pixels → 1e15 blowup; see [[project_orofacial_crop_dff_analysis]]).
- Can't just add PC1+..+PC5: sign is arbitrary + modes orthogonal (different motions). For a waveform use ONE signed PC; for "total movement" use RMS energy (but that's rectified → 2× freq).
- PC1 isn't ALWAYS breathing — script computes which PC has most breath-band power (`pick_breath_pc`) and prints a NOTE if it's not PC1. `SHOW_PC` overrides.
- `RE_CROP=true` re-crops every run (24 videos = crop each); `false` reuses cubes, redoes SVD only. MUST run in a real MATLAB (interactive cv2 crop window), not headless/MCP.
- 24 breath videos across ChAT/Sst/Vgat/Vglut2 under Ventral_surface_summary, each beside its calcium data. Cameras: acA2040-90umNIR (ChAT) / acA1300-200um (others); rawvideo/pal8.
