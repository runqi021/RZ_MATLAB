---
name: project_fix_avi_timing
description: fix_avi_timing.py — losslessly re-time Basler FFV1 AVIs to the true fps from timestamps.csv (fixes VLC playing too fast / wrong duration)
metadata: 
  node_type: memory
  type: project
  originSessionId: 2d9a96e6-14a8-4f27-8d60-c0c91c42e3bd
---

`fix_avi_timing.py` (repo root) — fixes the WRONG playback timing of Basler recordings.

**Root cause:** `basler_dual_acq.py` (see [[project_basler_acquisition_gui]]) stamps the FFV1 AVI with a NOMINAL fps (e.g. 200), but the camera runs slower (~111 fps). So VLC/players play ~1.8× too fast and show the wrong duration. The real per-frame times are only in the sibling `timestamps.csv`. This affects EVERY Basler run.

**Fix:** rewrite the video PTS at the true constant fps = `(N-1)/(t_last-t_first)` from timestamps.csv, via ffmpeg `setts` bitstream filter — **lossless STREAM COPY** (FFV1 frames untouched, ~15-20× realtime, bit-faithful). Output `<stem>_fixedfps.mkv` next to source; original `.avi` never touched.

```
ffmpeg -y -i in.avi -map 0 -c copy -bsf:v "setts=pts=N/(FPS*TB):dts=N/(FPS*TB)" out.mkv
```

**Gotchas learned:**
- `-r F -i in.avi -c copy` does NOT work — only relabels the header; actual PTS still encode 200 fps so players ignore it.
- `setts=pts=N` (no TB) also fails — N is in the original 1/200 timebase. Must scale by TB: `pts=N/(FPS*TB)`.
- `mkvmerge` (gold-standard `--timestamps`) is NOT installed on this machine; ffmpeg `setts` is the fallback and works.
- Verified on cam1_20260602_180811_run001: 33396 frames → 111.313 fps, output duration 300.0 s (was 166.98 s).

**Usage:** `python fix_avi_timing.py <run_folder | file.avi | parent_dir>` (batch finds every `.avi` with a sibling `timestamps.csv`; `--fps` override; skips existing unless `--overwrite`). This is a VIEWING fix only — the dF/F analysis already uses timestamps.csv (see [[project_orofacial_crop_dff_analysis]], [[feedback_detect_session_fps]]).
