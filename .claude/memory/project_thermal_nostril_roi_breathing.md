---
name: project_thermal_nostril_roi_breathing
description: Planned task — DLC L/R nostril tracking → circular ROI on thermal video → breathing signal in actual °C
metadata: 
  node_type: memory
  type: project
  originSessionId: 5b72e4e4-eb69-4350-aa72-7dd8848e4124
---

Upcoming task (raised 2026-06-17, DLC still training): use DLC nostril tracking (LEFT + RIGHT nostril positions) to place a circular ROI around each nostril and extract the **breathing signal in actual temperature (°C)** from the FLIR thermal camera video. Breathing modulates nostril temp (inhale cools, exhale warms), so the per-frame ROI temperature oscillates at the breathing rate.

**Good tracker (2026-06-17):** project `D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17`, shuffle2, snapshot_best-100. Bodyparts **L1 = LEFT nostril, L2 = RIGHT nostril** (named L1/L2 by mistake but that's the mapping). QC on 5840027_nose_n1: median ~0.92, 0% <0.6, longest low-conf run = 1 frame → interpolation effectively a no-op. NOTE: the earlier `260615_nose-RZ-2026-06-15` project regressed at iter-1/2 (corrupted refine labels: same video 0%→99% <0.6); abandoned, this redone Jun17 project replaces it.

**Build on existing thermal pipeline:** [[project_thermal_ats_breathing]] (.ats → fnv → thermal_ats_to_mat.py → MATLAB) and [[project_thermal_basler_sync_260613]]. Validate against SVD method [[project_breath_svd_pc1]] / thermal_breath_svd.m.

**ROI-definition method (user-specified 2026-06-17):** per animal, per nostril dot (L/R): for every frame crop a 30x30 px window CENTERED ON THAT FRAME'S DLC DOT (tracking-aligned), average the aligned crops over time → motion-stabilized average-projection image of the nostril; user draws an ellipse/circle on it. Store ROI RELATIVE TO DOT CENTER, then apply each frame re-centered on the moving dot to pull per-frame °C → breathing trace. For ROI-placement guidance use a BREATH-BAND POWER MAP (per-pixel FFT/PSD power at the breathing frequency, or correlation with a reference breath trace) computed on the ALIGNED stack — NOT raw temporal std. Reason (user pushed back 2026-06-17): std on RAW fixed-grid frames = a MOTION map (nose-edge sweep), useless. Alignment (dot-centered crop) removes gross translation first; even then plain std still picks up residual sub-pixel/edge jitter at the nostril rim. Breath-band power isolates pixels modulated at the breathing rate from motion/noise. This is the existing "nostril power-map" idea in thermal_breath_svd.m. Order: (1) DLC align, (2) breath-band power map, (3) draw/auto-place ellipse, (4) ROI-mean degC. MATLAB drawellipse (rotation; circle = special case); capture center-offset + semi-axes + angle.

**Before coding, confirm with user:**
- DLC ran on the THERMAL frames directly (so coords map 1:1)? If DLC tracked the Basler video instead, need thermal↔Basler registration (different camera/FOV) — much harder.
- Thermal .mat must hold **radiometric Celsius**, NOT display-normalized 8-bit frames (the whole point is °C). Per [[project_thermal_ats_breathing]] the Python→.mat split exists because R2021b can't pyenv py3.13.
- fps from timestamps, never the .ats header (it's wrong) — [[feedback_detect_session_fps]].
- ROI radius (px) — fixed, or auto from inter-nostril distance? Per-frame stat: mean vs max temp in ROI (max may track the warm exhale plume better).
- Output: per-nostril (L, R) °C trace over time, follow the moving nostril each frame.

**Scripts built (2026-06-17, single-video v1) — SEPARATED concerns (user: "draw tool is just for drawing"):**
- `thermal_nostril_breath_single.py` (flir env): .ats °C + DLC csv (numpy) → `_nostrilC.mat` (tracking-aligned crops, avg proj, raw disk trace, 100Hz stack).
- `thermal_nostril_roi_draw.m`: DRAWING ONLY — draw ellipse/circle on avg proj (bluewhitered), extract RAW °C ROI trace, save `_nostrilROI.mat`. No analysis here.
- `thermal_nostril_breath_analyze.m`: ALL analysis — bandpass+invert (inhale up), raw/bp figs, Hilbert instantaneous phase L vs R overlay + phase-diff + PLV, L+R averaged raw/bp. Consumes `_nostrilROI.mat`.
- `thermal_nostril_breath_view.m`: default-disk QC (pre-draw "is it good"), raw+bp.
Defaults: BP=[2 12] Hz, INVERT=true (inhale cools nostril → negate so inhale=rise). bandpass needs fps/2>BP(2) → stack at 100Hz.

**CANONICAL breathing writer — updated to `drawROI_N_lpsub_gui.m` (2026-06-20):** The old two-step (`thermal_nostril_roi_draw.m` → `thermal_nostril_breath_lpsub.m`) was replaced by a combined interactive GUI. The **actual extraction workflow** is: user runs `drawROI_N_lpsub_gui.m`, draws L/R nostril ROIs on the avg projection (live-updating traces), then Save → writes `<stem>_nostrilROI.mat` + canonical `<stem>_breath.mat`. Same LP-subtraction logic (LP_CUT=1 Hz, INVERT=true, FINE_BP=[2 10]). The `_breath.mat` struct is identical (B.breath = L+R avg, LP-subtracted, inhale-up °C; B.breath_bp; fps; t; method='lpsub'). `thermal_nostril_breath_lpsub.m` is now superseded. Other variants: `drawROI_N_lpsub.m` (single-run non-GUI), `drawROI_N_lpsub_batch.m` (batch dialog), `drawROI_N_lpsub_powermap.m` (diagnostic). Resolver `thermal_resolve_paths.m` returns P.breath. flir env has fnv+scipy+numpy but NO matplotlib/pandas → must read CSV with numpy, do all viz in MATLAB.

**Alignment chain (verified from _dlc.json next to each .ats):** thermal_dlc_video.py made the DLC video = 1 frame per .ats frame, decim=1, **upscale=6**. So `ats_frame = video_frame*decim` and `ats_px = dlc_px / upscale`. .ats native res 96x60, 32002 frames. .ats files at `D:\260615_thermalNbasler\<animal>\cam1_*\Rec-*.ats`; each has sibling `Rec-*_dlc.json` + `Rec-*_dlc.mp4` (the mp4 was copied+renamed into the nose project as `<animal>_nose_n#.mp4`). v1 uses a default disk ROI (r=2.5 native px) for an immediate °C trace; interactive ellipse drawing on the avg/power map is the NEXT step. Build on [[project_thermal_ats_breathing]] reader conventions (CELSIUS, TEMPERATURE_FACTORY, np.array copy of im.final).

Per [[feedback_no_inplace_edits]] create a NEW standalone script.
