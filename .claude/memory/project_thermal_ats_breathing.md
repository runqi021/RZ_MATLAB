---
name: project_thermal_ats_breathing
description: "FLIR .ats thermal breathing pipeline — fnv reader (conda env flir), thermal_ats_to_mat.py + thermal_breath_svd.m, nostril power-map QC"
metadata: 
  node_type: memory
  type: project
  originSessionId: 9757af76-e3c6-4f7f-93db-480b30c15ecc
---

Thermal breathing readout from FLIR A6753sc `.ats` files (data: `D:\260611_thermal_breathing\Active\<session>\Rec-*.ats`, 59 files, ~890 MB each).

**Reading `.ats` (the hard part):** `.ats` = proprietary FLIR radiometric sequence; nothing in the stack reads it natively (not VideoReader, not cv2). Solved via the **FLIR Science File SDK** (`FileSDK 2025.10.0`, installed to `C:\Users\Admin\AppData\Local\Programs\FLIR Systems\sdks\file\`). Its Python `fnv` wheel was pip-installed into a dedicated conda env **`flir`** (`C:\Users\Admin\.conda\envs\flir\python.exe`, Python 3.13 + numpy + scipy). Import: `import fnv, fnv.file`. Key API: `im=fnv.file.ImagerFile(path); im.unit=fnv.Unit.TEMPERATURE_FACTORY; im.temp_type=fnv.TempType.CELSIUS; im.get_frame(i); arr=im.final` (reshape to `im.height x im.width`); `im.frame_info.time` is a datetime per frame.

**Why a separate Python env, not MATLAB pyenv:** MATLAB here is R2021b → only drives Python ≤3.9 via pyenv, so it CANNOT load the 3.13 fnv env. Pipeline is therefore split: Python reads → compact `.mat` → MATLAB does the SVD.

**Two scripts** (in `thermal_breathing/`):
- `thermal_ats_to_mat.py` — runs in `flir` env. Reads `.ats` as °C, derives TRUE fps from timestamp span (native ≈396 Hz, header's 400 is nominal), anti-alias decimates to ~50 Hz, saves `*_thermal.mat` (`stack [T×H×W] single °C`, `fps`, `fps_orig`, `t`, `H`,`W`,`T`). ~142 MB/file.
- `thermal_breath_svd.m` — pure MATLAB. Loads `*_thermal.mat`, SVD of per-pixel mean-subtracted temperature, sign-stabilized breath PC, + per-pixel **breath-band power map** (FFT band-integral) = the QC image. Mirrors `breath_svd_pc1.m` conventions (mt_psd, pick_breath_pc); `breath_svd_pc1.m` NOT touched. Writes `*_thermal_svd.mat` + `.png`.

**Reader VERIFIED correct** against Research Studio: `fnv` (with `unit=TEMPERATURE_FACTORY, temp_type=CELSIUS`) vs the user's CSV export of the same frame (matched by timestamp to 0.00 ms) = mean abs diff **0.014 °C**, max 0.13, corr 0.9998. GOTCHA: Research Studio CSV export UPSAMPLES 3× per axis (192×80 native → 576×240 csv); collapse 3× before comparing. SECOND GOTCHA: if you forget to set `im.unit`, fnv returns raw COUNTS (~3000), not °C.

**CRITICAL — dropped frames / timestamps:** camera set to 400 Hz (24000 frames/60 s) but DROPS frames (200–2000 per file; -0015 lost 235 = 1.0%, scattered single-frame drops). Frames present are therefore NOT uniformly spaced. The reader MUST resample onto a uniform grid from the real per-frame `frame_info.time` (linear interp fills gaps) BEFORE decimating — naive integer decimation assuming uniform spacing corrupts the time axis. `thermal_ats_to_mat.py` now does this (saves `fps_nominal`, `fps_mean`, `n_dropped`, `n_expected`). frame_info.time is a placeholder-year (1976) IRIG datetime; fnv index 0-based = camera FrameNumber−1.

**Viewing the movie (no FLIR GUI):** `thermal_make_video.m` (MATLAB, RELIABLE via MCP) reads a `*_thermal.mat` → `*_view.mp4`, FIXED linear temperature window (no gamma/AGC, per [[feedback_never_gamma_video]]), hot LUT, 6× nearest zoom. `thermal_make_video_native.py` (cv2, full 400 Hz no-decimation) is correct code but the long invocation hit a harness empty-output exit-127 glitch (cv2/fnv/imports all work in isolation) — prefer the MATLAB path. KEY framing the user insisted on: the timestamp-resampled video IS the correct one; a "raw frames played uniformly" video would be WRONG (dropped frames → compressed gaps). 50 Hz decimation loses only >25 Hz, so breathing is fully intact.

**Open analysis issues (NOT reading):** (1) decimation leaves big edge-transients at first/last sample → trim ~1 s each end before SVD/PSD. (2) slow whole-face thermal drift + sniff-bout envelope (<1 Hz) carry ~4× the breath-band power, so a naive SVD PC1 / peak-pick reports a bogus slow rate — need per-pixel high-pass (~1 Hz) before SVD. (3) signal is real bursty awake sniffing (active bouts ~3–5 Hz ripples on ~1 Hz envelope, separated by quiet/flat periods); verify flat periods aren't head MOTION moving the nostril off fixed pixels (may need motion-tracked nostril ROI). Power map (nostril hotspot) is robust regardless.

Related: [[project_breath_svd_pc1]] (Basler-video analog), [[feedback_detect_session_fps]], [[feedback_never_gamma_video]] (radiometric = no AGC/gamma, concern doesn't arise).
