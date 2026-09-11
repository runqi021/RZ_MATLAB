---
name: project_bilateral_whisk_breath_pipeline
description: "Bilateral whisking synchronization + breathing-coordination analysis pipeline (260615 thermal+Basler); scripts, event/peak conventions, camera-sync facts"
metadata: 
  node_type: memory
  type: project
  originSessionId: df3ac351-d090-46b7-9550-12bbcefb3f43
---

Event-timing analysis of **left/right whisker synchrony** and **whisking↔breathing coordination** for the 260615 thermal+Basler dataset. Whisk = Basler camera DLC; breath = thermal nose camera `_breath.mat`. See [[project_thermal_basler_sync_260613]], [[project_260615_whisk_breath_coordination]], [[project_thermal_ats_breathing]].

## Camera sync (CRITICAL — settled with user 2026-06-20)
- Basler (whisk) + thermal (nose) are driven by the **same 400 Hz WFG trigger → synchronized FRAME-TO-FRAME**. Frame i ↔ frame i, always. One camera may **overrun** (thermal is usually shorter); just truncate to common length. **No lag/circshift correction belongs in analysis** — `(0:N-1)/fps` with both starting at t=0 is correct.
- `whisk_breath_sync_report.m` cross-correlates the two breathing proxies to a `LAG` and saves `sync_result.lag_frames` — this is a **QC confirmation of frame sync, NOT a correction to apply downstream**.
- Sub-frame exposure centroid is **negligible**: real values in `thermal_breathing/thermal_basler_sync_video.py` (`FLIR_EXP=0.97648e-3`): **Basler 0.25 ms (center +0.125), FLIR/thermal 0.97648 ms (center +0.488)** → center diff ≈ 0.36 ms. (The `INTEG_MS=2.3754` in `whisk_breath_sync_report.m` is OUTDATED — fix if touched.) Do NOT add any ½-integration shift. **Why:** spent a long detour debating exposure-centroid + LAG shifts; user clarified frame-to-frame trigger makes all of it moot.

## Event convention — STANDARD (settled 2026-06-20, final)
- **Breath-triggered / coordination scripts use ONSETS, and positive/right = AFTER the trigger.**
- **Protraction onset** = local MINIMUM (trough) of BP whisk angle (`findpeaks(-x)`) = start of forward sweep. Markers `v`.
- **Inspiration onset** = local MINIMUM (trough) of inhale-up breath (`findpeaks(-brf)`) = start of inhalation. (NOT the peak/peak-inhalation.)
- Sign: x = `event − trigger`, so `x>0` = event after trigger; L-R lag = `tL − tR`, `+` = L after R.
- History: I briefly flipped everything to PEAKS mid-session, then user reverted to ONSETS for breath-triggered scripts (raster, coord_6panel). `whisk_breath_trig_avg/amp_freq/freq_hist/detect_diag/breath_freq_whisk_psd` were always onset-based. **EXCEPTION — still PEAKS by explicit user request:** `whisk_LR_sync_RZ` (L-R protraction-peak synchrony), `whisk_LR_5s_random_RZ`, `whisk_LR_breath_random_RZ` (L-R only, no breath trigger).
- Whisk angle base→tip (deg), **LEFT x mirrored so protraction is +up on both sides**: L `rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),-(M(:,5)-M(:,2)))))`, R `rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),M(:,11)-M(:,8))))`. DLC csv via `pick_csv` (best-N snapshot), `readmatrix(...,'NumHeaderLines',3)`.

## DLC likelihood gating (standard, 2026-06-20)
- **Tracking points with likelihood < 0.6 → NaN → linear interpolation over frames** (endpoints held). Shared MATLAB helper `dlc_gate_interp.m` (repo root): drop-in for `readmatrix(csv,'NumHeaderLines',3)`, cleans every bodypart's (x,y), leaves likelihood cols.
- **Nasal already done in Python**: `thermal_breathing/thermal_nostril_breath_single.py` `gate_interp()` with `LIK_THRESH=0.6` (and `_batch`) — nostril centers gated+interp at crop-build time. MATLAB helper mirrors it.
- **Wired into ALL whisker-csv readers (2026-06-20):** whisk_breath_raster, whisk_breath_coord_6panel, whisk_LR_sync, whisk_LR_phase, whisk_LR_5s_random, whisk_LR_breath_random, whisk_protraction_onset_methods, breath_whisk_overlay, whisk_detect_diag, onset_hilbert_diag, trig_detect_overlay, whisker_DLC_hilbert, whisk_epoch_analyze/_auto/_quiet, whisk_epoch_gui, whisk_epoch_single_overlay, whisk_ITI_dist/_epoch, breath_freq_whisk_psd, raster_L_trig_R. (whisk_breath_trig_avg reads breath only → N/A.) `dlc_clean_whisker.m` is a SEPARATE kinematics-based cleaner (NOT likelihood) writing _clean.mat that nobody loads — left alone.
- Whisker csv layout: col1=frame, then (x,y,lik) triplets; bodyparts vL-base/vL-tip/vR-base/vR-tip → lik cols 4,7,10,13.

## Detection params (shared across scripts; keep consistent)
- `fpsW=400`; whisk `BP=[6 30]` (6-panel still `[5 30]`); breath `BR_BP=[2 15]`, `BR_PROM=0.5`.
- `MINDIST` = min peak spacing. **0.02 s was too loose (allows 50 Hz → too many peaks per protraction)**; raise toward 0.03–0.05. Raster currently 0.03.
- `PROT_PROM` ≈ 0.8–1.0 ×std. Whisk-epoch gating: envelope `(|hilbert(xL)|+|hilbert(xR)|)/2`, `detect()` with `THR_FRAC*prctile(env,95)`, `MIN_DUR`, `MERGE_GAP`, `ENV_SMOOTH`; peaks restricted to epochs (`iL=iL(m(iL))`). `THR_FRAC=0` = pass-through (no gating).

## Active working folder: `breathing_whisking_coord/` (2026-06-22)
Current scripts under active development. Data import chain shared by all of them:

**Session discovery:** numeric-named subdirs in `dataRoot`; cross-referenced with `cam1_*` run folders (sorted chronologically) for session index `kk`; whisk CSV exists in `whiskDir` → session included.  
**Whisk import:** `pick_csv(whiskDir, '{animal}_whisk_n{kk}')` → `dlc_gate_interp(csv, 0.6)` → `atan2` angle (L: col 2-6, R: col 8-12) → `filtfilt(butter(3,[5 30]/(400/2),'bandpass'))`. Hilbert envelope from (|hilbert(xL)|+|hilbert(xR)|)/2. `fpsW=400`.  
**Breath import:** `thermal_resolve_paths(pick_csv(noseDir, '{animal}_nose_n{kk}'))` → `Pn.breath` → `load()` → `Bs.breath` + `Bs.fps`.  
**Epoch gating:** `detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP)` with defaults THR_FRAC=0.2, MIN_DUR=0.5, MERGE_GAP=0.2.  
**EXCLUDE = "5840027"** across all scripts.  
**Key threshold params:** BASAL_HZ=5 (some scripts use 6), SNIFF_HZ=8, SLOW_HZ=10, FAST_HZ=20, SYNC_THRESH_MS=7.7, MINDIST=0.025 s, PROT_PROM=1×std, BR_PROM=0.5×std.

| Script | Purpose |
|--------|---------|
| `whisk_LR_sync_combined_10panels_RZ.m` | 10-panel L/R protraction-peak synchrony (\|lag\| histogram); rows: ITI/all / slow+fast / basal+sniff / B×slow+B×fast / S×slow+S×fast |
| `whisk_epoch_analyze_auto_RZ.m` | Epoch-gated scatter plots (xL vs xR, breath vs xL/R) split by basal/sniff; + 3-panel triggered averages (R-trig-L, insp-trig-L/R); both angle and velocity |
| `breath_freq_whisk_psd_RZ.m` | 4-panel: breath freq histogram / whisking PSD (basal vs sniff) / whisk freq from sub-frame onset timing / amplitude vs freq |
| `whisk_LR_lag_vs_breath_whisk_freq_RZ.m` | Continuous \|lag\| vs breath freq (scatter+bin), vs whisk freq, + 2D heatmap (breath-rate × whisk-freq) |
| `whisk_meta_shuffle_controls_RZ.m` | Circular-shift shuffle controls for all 3 metric types (L/R lag / breath phase / peri-insp PSTH); 500-rep null; figures 1-3 |
| `drawROI_N_lpsub_gui.m` | Interactive nostril ROI GUI (thermal crops → `_nostrilROI.mat` + `_breath.mat`) |
| `thermal_nostril_breath_single.py` | Python (flir env) single-session thermal `.ats` → nostril breathing signal `.mat` |

**Analysis method change incoming** — user stated the analysis approach will be replaced; data import structure above stays stable.

## Scripts (older, outside breathing_whisking_coord/)
- `whisk_LR_sync_RZ.m` — L-R protraction-peak synchrony, pooled all sessions. Pair each L peak to nearest R peak within `MAXLAG=0.07`; lag = L−R (ms); sync window = `SYNC_FRAC`×mean whisk cycle. **Lag histogram bin CENTERED on 0** (edges at ±bwid/2, ±3bwid/2…; `bwid=2000*MAXLAG/NBINS`).
- `whisk_breath_raster_RZ.m` — breath-cycle-triggered raster of protraction peaks vs inspiration peak. **Two side-by-side panels (LEFT | RIGHT whisker, not pooled)**, same breath order (sorted by duration). Window `tInsp(i)-0.2 → tInsp(i+2)`; `xlim([-0.2 0.4])`; red curve = breath duration.
- `whisk_breath_coord_6panel_RZ.m` — single-session 2×3 figure: (1) L/R whisk+peaks, (2) L-R peak-lag hist, (3) L-R lag over time, (4) insp-triggered protraction raster, (5) first-protraction latency hist, (6) whisks/breath hist. Circular PLV stats (respiratory phase + L-R Hilbert) printed to console only.
- `whisk_LR_5s_random_RZ.m` — random `WIN`-s window (user set 3 s) from a whisking epoch (length-weighted pick); L/R BP traces + protraction-peak `^` markers. `RNG_SEED` to reproduce.
- `whisk_LR_breath_random_RZ.m` — random window L/R whisk (left y-axis) + breathing (right y-axis); `FIX_T0` pins window start. **Peak markers removed — traces only** (user request).

## Nostril breath extraction GUI (drawROI_N_lpsub_gui.m, thermal_breathing/)
- Interactive replacement for drawROI_N_lpsub_batch: two 24×24 nostril crops (WIN_NATIVE=24, ±12 px around DLC nostril point) with draggable ellipse/circle ROIs, L/R/avg detrended traces, ellipse⇄circle toggle, numeric Size + −/+ buttons (circle handle is unusable on a 24-px image), auto-load saved ROI, Prev/Next, Save/Save&Next (same `_nostrilROI.mat` + `_breath.mat` format).
- **Trace updates ONLY on ROI release (ROIMoved), NOT live during drag** — live MovingROI updates were too laggy (filtfilt+render of a long trace per mouse-move). Draw/move ROI → release → plots. Display trace is decimated to ~3000 pts; Save uses full-res. Auto ROI / load / size-change also trigger a replot.
- **Breathing band is ~2–10 Hz.** "Auto ROI" button + "SNR view": per-pixel map = **breath SNR = power(2–10 Hz) / power(10–40 Hz noise band)**, absolute-signal floor (pixels <10% of max in-band power zeroed). **CRITICAL: SNR peaks at the crop EDGES — high-contrast nose-boundary motion makes huge 2-10 Hz swings (artifact), NOT breathing.** So Auto ROI must **search ONLY within MAX_OFF=4 px of the crop center (the tracked dot)** then take the brightest SNR there (3×3 refine); searching the whole crop grabs an edge and the leash just clamps it to the ring → "absolutely failing". Verified on real data: global argmax→(16,20) edge; within-4px→(13,16) real breathing pixel (agrees with in-band-power peak). Size fixed = user Size, ≤3 px; no covariance/rotation (unstable on few-px blob).
- **CONTRA/MIDPOINT ANCHORING TESTED → IT HURTS, DO NOT USE for extracting a nostril's own breath.** Geometry: the own dot best tracks its own nostril; anchoring that crop to the contralateral dot (or midpoint) injects the inter-dot RELATIVE jitter (~1.4 px) into the nostril → smears the breathing pixel. Same-pixel test on 5916300 R: SNR 6.44 (self) → 1.71 (contra); near-center self L 8.77/R 6.44 vs contra ~0. The earlier "contra is better" impression + the 15 s anchor_compare (5916296 L 2.82→3.19) were small/noisy/confounded; R was always best with self. **Keep ANCHOR='own' / GUI anchor='self'.**
- GUI still has an `anchor` dropdown self/contra/midpoint (per-frame `stabil` re-shift from stored L_center/R_center; circshift with NaN-nulled wrapped border) and the extractor has matching ANCHOR — both DEFAULT to self/own and are kept only for A/B judging. `thermal_anchor_compare.py` does the offline 3-way comparison. `breathPowerMap()` actually returns the SNR map. The breathing TRACE itself is still the LP-subtracted ROI mean (SNR only drives pixel selection). FINE_BP [1 20]→[2 10]. Bigger FOV needs larger WIN_NATIVE re-extract (Python). uifigure gotcha: set per-axes colormap with a MATRIX (`hot(256)`) AFTER CLim + `drawnow`, else the swap silently no-ops and power renders through bluewhitered. See [[project_thermal_ats_breathing]] [[project_breath_svd_pc1]].

## Full file inventory (2026-06-20 sweep)
- **Core analysis:** whisk_LR_sync_RZ, whisk_breath_raster_RZ, whisk_breath_coord_6panel_RZ, whisk_LR_5s_random_RZ, whisk_LR_breath_random_RZ, breath_whisk_overlay_RZ, whisk_LR_phase_RZ, raster_L_trig_R_RZ.
- **Whisk detect/epoch/diag:** whisk_detect_diag_RZ, whisk_epoch_gui_RZ, whisk_epoch_analyze_RZ / _auto / _quiet / _single_overlay, whisk_ITI_dist_RZ, whisk_ITI_epoch_RZ, whisker_DLC_hilbert_RZ, onset_hilbert_diag_RZ, trig_detect_overlay_RZ.
- **Whisk+breath combined:** whisk_breath_trig_avg_RZ, whisk_breath_detect_diag_RZ, whisk_breath_freq_hist_RZ, whisk_breath_amp_freq_RZ, breath_freq_whisk_psd_RZ.
- **Camera sync:** whisk_breath_sync_extract.py + whisk_breath_sync_report.m (QC); thermal_breathing/thermal_basler_sync_video.py (CORRECT exposures) + _sync_batch.py + _epoch_video.py + check_sync.py; sync_motion_extract.py / sync_motion_report.m / sync_combined_video(.py/_batch.py).
- **Breath signal gen (thermal nose→_breath.mat):** thermal_breathing/drawROI_N_lpsub.m(+_batch); thermal_nostril_breath_{single,batch}.py/_lpsub.m/_analyze.m/_view.m; thermal_nostril_roi_draw.m/_linescan.m; thermal_breath_svd.m; thermal_ats_to_mat.py; thermal_resolve_paths.m (path helper used everywhere); thermal_make_video*/thermal_dlc_video.py/thermal_inspect_videos.py.
- **Whisker DLC tracking + crop (→whisk csv):** DLC_tracking/whisker_{LR_batch,L_3w_batch,phase_single,pc1_velocity,dlc_mtspec(_test)}.m; whisker_{label_video_batch,label_video_test,select_crop}.py; whisker_crop_{extract,batch}.py / whisker_crop_track.m / whisker_roi_check.py.
- **Raw acquisition:** Bassler/basler_gui.m, basler_dual_acq.py, run_basler_n_runs_n_cam.m.
- **NOT part of this pipeline** (separate calcium-breath work): Breath_summary_timeNphase_260331/*, svd_breath_motion/*, 2p_breathing_coherence/*, breathing_is_motion/*, chat/vglut_breath_coherence*, IO_population_*.
- Exact exposures fixed in whisk_breath_sync_report.m: INTEG_MS 2.3754→0.97648 ms (the sub-frame correction is cosmetic; changes NO results — see camera-sync section).

## Working rules honored
New standalone `*_RZ.m` files for new analyses; only edited existing scripts in place when user explicitly asked (raster, sync). See [[feedback_no_inplace_edits]], [[feedback_detect_session_fps]]. Plot recommendation came from user: event-timing (raster/lag-hist/latency/whisks-per-breath), NOT continuous correlation/linear fit.
