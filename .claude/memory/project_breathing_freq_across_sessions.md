---
name: project_breathing_freq_across_sessions
description: breathing_freq_across_sessions_260629.m — per-session breath rate (1/ITI) from DLC across a long multi-run recording; fps=50.68 not 90
metadata: 
  node_type: memory
  type: project
  originSessionId: 5352112a-72f5-4c22-93bb-da0378863e7d
---

`breathing_freq_across_sessions_260629.m` (repo root): per-session breathing frequency across a long multi-run behavior recording, straight from DLC csvs already on disk (no video read, no SVD). Loops run folders under `folderPath`, takes the breathing dot (default `dot4`, **inverted `-x`** so findpeaks lands on breath events), high-pass (movmedian) + z-score, findpeaks → ITI → freq=1./ITI, per-run mean±std. Plots mean±std + median vs wall-clock recording time, plus #breaths/run QC bar. Outputs `breathing_freq_across_sessions.png/.pdf/.csv/.mat` into folderPath.

First run on `Z:\Runqi\behavior\260505_breathing_wt\cam1` (112 runs, 5-min interval, ~60s each): aroused ~5 Hz start → stable ~3 Hz plateau ~6.5 h → drop to ~1.5-2 Hz at ~7 h (animal quieting). Last 2 runs unreliable (few breaths, huge error bars).

**fps gotcha:** the old `mtspec_dot4_x_per_run.mat` in that folder stored `fps=90`, which is WRONG (nominal header). Real rate from `timestamps.csv` `camera_timestamp_ns` median-diff is **50.68 Hz** — its PSD freqs were inflated ~1.78×. Script derives fps per-run from timestamps and reads wall-clock recording time from `wall_time_s` (posixtime). See [[feedback_detect_session_fps]].

Refractory = `MINDIST_SEC` (0.10 s default = 10 Hz ceiling). Z: network reads are flaky — script has a 4× readmatrix retry. Related: [[project_basler_acquisition_gui]] (timestamps.csv format), [[project_chat_analysis_breath_alignment]].
