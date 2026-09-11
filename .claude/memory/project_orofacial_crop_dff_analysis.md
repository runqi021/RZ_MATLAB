---
name: project_orofacial_crop_dff_analysis
description: Cropped-ROI spatial/temporal/SVD analysis of Basler orofacial video via per-pixel dF/F (two-script Python→MATLAB)
metadata: 
  node_type: memory
  type: project
  originSessionId: 2d9a96e6-14a8-4f27-8d60-c0c91c42e3bd
---

Two-script tool to explore a cropped region of a Basler orofacial behavior video in space + time.
Created 2026-06-05 for `260602_orofacial_ChAT_ChroME`.

**Single MATLAB entry point** `orofacial_crop_analyze.m` (repo root): calls the Python extractor via `system()` if `crop_for_analysis.mat` is missing (delete to re-extract — simple delete+redo, [[feedback_rerun_logic]]), then does all analysis. `PYTHON_EXE = C:\Program Files\Python314\python.exe` (has cv2 4.13 + scipy; decodes FFV1 fine). Extraction params (time clip `EXTRACT_T_START/END`, `TEMPORAL_STRIDE`, `SPATIAL_BIN`, `ROI_XYWH=[]`→interactive) all set at top of the `.m`.
- `orofacial_crop_extract.py` — argparse-driven (MATLAB passes `--avi/--ts/--out/--t_start/--t_end/--temporal_stride/--spatial_bin/--roi/--fps/--fullframe`); OpenCV reads FFV1 (MATLAB VideoReader CAN'T), interactive `cv2.selectROI` on mean-of-60 preview when `--roi` omitted, auto spatial-bin to fit MAT5 2GB, saves uint8 cube `[H,W,T]` + `t_s` + `fps_eff`. Also runs standalone.

**Generalized 2026-06-05** to handle non-pipeline videos: time-axis source priority = `--fps` override → `timestamps.csv` → AVI-header fps (`build_timeaxis`). `--fullframe` skips crop/ROI entirely. Handles `rawvideo`/`pal8` Pylon-Viewer AVIs too (not just FFV1). Verified on a ChAT breath cam (`Basler_acA2040...avi`, 512x512, 3000 frames): no timestamps, header 98.62 Hz is BOGUS (2P-frame-triggered breath cam → true rate = imaging 30 Hz, see [[project_chat_analysis_breath_alignment]]), so ran `--fullframe --fps 30`. In the `.m`: `FULL_FRAME=true`, `FPS_OVERRIDE=30`.

**SIGNAL param (added 2026-06-05): `'raw'` (default) | `'dF'` | `'dff'`.** dF/F is WRONG for raw behavior video — it divides by baseline F0, and dark/background pixels (F0~0) explode to ~1e15, swamping every map/PSD/SVD (caught on the ChAT breath cam: std(dFF) max 8e15). Use `'raw'` = intensity − temporal mean (no division; robust), or `'dF'` = intensity − sliding-median baseline (linear, removes drift). `'dff'` = `(I-F0)/F0` is ONLY for calcium ROIs with real baseline F. Result on breath cam with 'raw': std map median 1.4 / max 83 (sane), breathing peak ~0.83 Hz, SVD mode1 = 66% var. Downstream code uses generic `X` matrix + `sigLabel`.

**Original framing — per-pixel dF/F, NOT motion energy (for CALCIUM data):** User correction: motion energy `|I(t)-I(t-1)|` is rectified → nonlinear → frequency-doubles/harmonics ("double counts"; ok for whiskers only). dF/F `(I-F0)/F0` is linear → preserves true breathing/whisking frequencies + normalizes baseline brightness. Reshape cube → `F [T×Npix]` → `helper.dFF_RZ(F,'FPS',fps_eff,'BaselineWinSec',20,'UseOASIS',false)` (see [[feedback_helper_package_prefix]]).

**Video facts:** FFV1, 944×920, 33396 frames, **300.0 s**, AVI header says **200 fps but that's WRONG** — real ~**111.3 fps** from `timestamps.csv` (camera_timestamp_ns col). Always use timestamps for time axis (cf [[feedback_detect_session_fps]], [[feedback_tiff_metadata]]).

**Figures:** (1) mean frame + std(dFF) + peak|dFF| heatmaps; (2) 2D spatial FFT log-power + radial spatial PSD; (3) global trace + temporal PSD (Chronux `mtspectrumc` TW=3, fpass→min(25,fps/2), breathing/whisking bands shaded) + kymographs; (4) `svds` SVD (auto spatial-bin to <2GB double) eigenimages + temporal modes + var-explained. Reuses SVD idiom from `SVD/Run_ROI_QC.m`, PSD idiom from `breathing_is_motion_analysis.m`. See [[project_chronux_multitaper_conventions]].
