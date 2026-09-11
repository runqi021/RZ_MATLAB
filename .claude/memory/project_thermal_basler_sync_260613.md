---
name: project-thermal-basler-sync-260613
description: "260613 thermal+basler whisking/breathing dataset - sync assessment, scripts, frame-timing facts"
metadata: 
  node_type: memory
  type: project
  originSessionId: 18f3e613-c0b7-4c78-8b44-95de195de636
---

**Dataset:** `C:\260613_breathing_thermalNbasler\<scope>\cam1_<ts>_run001\` — each run folder has a Basler `.avi` (FFV1, 992x504 BGR), a FLIR thermal `.ats` (56x176), and `timestamps.csv`. Goal: bilateral whisking synchrony + breathing-whisking coordination. First pass = run `cam1_20260613_181906_run001`.

**Acquisition:** Basler + thermal both hardware-triggered by ONE 400 Hz WFG. Basler = `basler_dual_acq.py` hardware FrameStart on Line3 falling edge (exposure STARTS at sync-in; `camera_timestamp_ns` latched there; exposure was 2200 us this session). Thermal integration 2.3754 ms. Per-run we save ONLY `timestamps.csv` (frame_idx, camera_timestamp_ns, wall_time_s) — NO exposure/trigger sidecar (consider adding `run_meta.json`).

**Frame timing:** both cameras START exposure on the same pulse; effective frame CENTER = sync-in + exposure/2 (Basler +1.10 ms, thermal +1.19 ms) → differ ~0.09 ms, cancels. Absolute clocks NOT comparable (Basler PC-epoch/counter vs FLIR internal IRIG ~3 h off) → frame alignment must come from a shared physical signal, NOT timestamps.

**Sync findings (181906):** both clean 400 Hz, ZERO internal drops; Basler 53894 fr, thermal 52002 fr → thermal stopped 1892 fr (4.73 s) earlier. Best-guess offset LAG=0 (both caught first pulse, thermal ended early). NOT yet hard-confirmed: breathing waveforms (thermal airflow-temp vs Basler marker-motion) don't cross-correlate (corr ~0); thermal frame-diff motion is sensor-noise dominated. → verify offset by EYE from the combined video (shared movement simultaneous in both panels) and set `--lag N` if needed.

**Channels:** thermal = breathing only (nostril cool-spot bottom-center; whiskers NOT thermally resolved). Basler = bilateral whisker fans (each cheek a comb of bright marker streaks) + nostril dots + both eyes. Whisker tracking plan (MATLAB, no DLC needed to start): per-side ROI → threshold → `regionprops` WeightedCentroid/Orientation (or `vision.PointTracker` KLT / `normxcorr2`).

**New scripts (this session, all parameterized, save outputs into the run folder):**
- `thermal_breathing/thermal_inspect_videos.py` — .ats → `_track.mp4` (gray) + `_heat.mp4` (inferno); per-frame timing/drop report.
- `basler_preview.py` — FFV1 .avi → decimated `_preview.mp4` + `_proj.png` (mean|max; dots brightest in max).
- `sync_combined_video.py` — synchronized `_combined_lag{N}.mp4` (Basler top / thermal bottom, paired basler=thermal+lag); the deliverable + visual sync check.
- `whisk_breath_sync_extract.py` + `whisk_breath_sync_report.m` — breathing-proxy cross-corr (gave weak corr; superseded by motion + eyeball).
- `sync_motion_extract.py` + `sync_motion_report.m` — motion-energy cross-corr.
- `sync_combined_video_batch.py` — drives sync_combined_video.py over ALL run folders (recursive `**/*.ats` catches archived/arhived; atomic `.combining.lock` lets 2 instances run concurrently; `--filter hived` = archived only). 260613 root has 31 recordings (17 active + 14 archived).
- `whisker_crop_extract.py` — INTERACTIVE: draw 2 whisker ROIs on a max/blend projection → crops both across all frames (auto size-budget bin) → `<run>_whiskcrop.mat` (mov_<name> uint8 [H,W,T], roi/bin, t_s, fps) + per-ROI `_crop.mp4` + `_whisk_rois.png`. Run via the flir env python; prints "cropped N frames" BEFORE writing the big mat (don't load mid-write).
- `whisker_crop_track.m` — whisking trace per side from the crop. The bright marker is a PAINTED WHISKER STREAK rotating about its follicle (NOT a compact dot). Robust readout = project the MOVING-marker centroid `max(frame-mean,0)` weighted onto the sweep PC axis (bounded, never flips); band-pass [4 25] Hz, L-R xcorr (tight ±60 ms lag) + mscohere. Rejected: absolute-brightness orientation/angle (locks on static follicle, ±180 flips, unwrap runaway). 190730 pilot: whisking is bursty; bilateral near-synchronous (corr 0.50 @ -38 ms, coherence 0.60 @ 16.5 Hz). Outputs `<run>_whiskcrop_track.{png,mat}`.

- `whisker_crop_batch.py` — batch the 2-ROI whisker crop over ALL runs (active+archived). ROIs are PER-SCOPE (camera frame differs per scope: 5840027 976x406, 5916296 960x460, 5916297 880x349, 5916300 992x504). Phase 1 `--draw` (interactive, caches `<scope>/whisker_rois.json`; a scope already cropped is seeded from its `_whiskcrop.mat` — 5840027 auto-seeds from 190730). Phase 2 (no flag) crops every run unattended (~8-10 min each, skip done). Run in flir env.
- `DLC_tracking/dlc_setup_5916297.py` + `dlc_setup_per_scope.py` — create ONE DLC project per animal/scope (`<scope>_experimental-RZ-<date>`), full raw Basler videos incl. archived, bodyparts vL0/vL1/vR0/vR1, net resnet_50, numframes2pick 20, FRAMES NOT EXTRACTED yet. Done for all 4 scopes (5840027/5916296/5916297/5916300). Run in dlc310 env (`C:\Users\Admin\.conda\envs\dlc310\python.exe`). DLC can't symlink (no admin/exFAT) so it copies+deletes a seed video (smallest picked to minimize). Next: extract_frames + label_frames per project.

**Next:** breath-whisk coordination (thermal/Basler nostril breathing vs the whisk traces in `_whiskcrop_track.mat`); label+train the DLC projects.

**Gotcha:** FLIR fnv `im.final` is a reused buffer view — must `.copy()` or stored frames alias → static video / zero motion. See [[reference-fnv-buffer-copy]]. Related: [[project_thermal_ats_breathing]], [[project_basler_acquisition_gui]], [[feedback_ffv1_video_codec]].
