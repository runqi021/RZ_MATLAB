---
name: pair-behavior-to-phys-m-script-reference
description: "Script that pairs behavior AVIs + DLC CSVs with calcium FOV folders by timestamp, copies, and trims"
metadata: 
  node_type: memory
  type: reference
  originSessionId: d0d6f257-9f5d-4554-aa44-3425c396829e
  modified: 2026-07-23T01:20:31.375Z
---

`C:\Users\Admin\Desktop\RZ_MATLAB\pair_behavior_to_phys.m`

- Discovers FOV folders recursively under `physDir` (handles intermediate category subfolders)
- Matches by nearest `|TIF_epoch - AVI_filename_timestamp|` within 120 sec tolerance
- Shows pairing table for user confirmation before copying (`dryRun = true` = safe preview)
- Uses `_meta.mat` `framesPerSlice * numSlices` for raw TIF frame count (not cpSAM which has frame-dropping)

**Dual-format update (2026-07-22):** now handles BOTH Basler naming/formats automatically (verified on 260721_Sert_soma_G8s: 26 new + 1 legacy, 0 unparsed).
- **AVI discovery is now recursive** (`behavDir/**/*.avi`), because the new dual-cam GUI nests videos as `<camName>/<camName>_<YYYYMMDD_HHMMSS>_run<NNN>/*.avi`. Legacy AVIs sit flat.
- **Two filename→timestamp regexes:** new `_(YYYY)(MM)(DD)_(HH)(MM)(SS)_run\d+\.avi$` (no ms, run suffix) checked first; legacy `_(YYYY)(MM)(DD)_(HH)(MM)(SS)(mmm)\.avi$` (3-digit ms). Unambiguous (run-suffix vs 9-digit-time). Cam start time ≈ TIF ScanImage epoch within ~2s (confirmed: cam1_..._175021 vs epoch 17:50:22.8).
- **CODEC: new videos are FFV1 (lossless), NOT MJPG** — the MATLAB GUI `run_basler_n_runs_n_cam.m` is just a wrapper that calls the SAME `basler_dual_acq.py` FFV1 encoder (see [[project_basler_acquisition_gui]], [[feedback_ffv1_video_codec]]). MATLAB `VideoReader` CANNOT decode FFV1. So new-format trim goes through **ffmpeg `-frames:v N -c copy`** (lossless, frame-exact since FFV1 is all-intra); `ffmpegExe` auto-detected (PATH / `dlc310` env / `cellpose-gpu` env; found at `...\.conda\envs\dlc310\Library\bin\ffmpeg.exe`). Legacy trim still uses VideoReader→VideoWriter.
- **New-format frame count** comes from the sibling `timestamps.csv` row count (`count_ts_frames`, no video decode), NOT VideoReader. If timestamps.csv absent → NaN → copy untrimmed + warn. `timestamps.csv` is copied into the FOV folder as `<aviBase>_timestamps.csv` (trimmed to header+nKeep rows if the video was trimmed).
- Still greedy 1:1 AVI→FOV (single-camera assumption). If a session records cam1 AND cam2, only one would pair per FOV — revisit if dual-cam is used.
- DLC CSV matching (`behavDir/**/*DLC*.csv`) retained but harmless when absent (new sessions have no DLC).
