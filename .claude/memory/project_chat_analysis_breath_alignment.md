---
name: project_chat_analysis_breath_alignment
description: ChAT_analysis dataset — breath/calcium are frame-triggered 1:1; stored camera_fps is a bogus placeholder; toss-30-then-truncate alignment
metadata: 
  node_type: memory
  type: project
  originSessionId: 8681f366-ed2b-4752-b32e-0201447c7c20
---

`C:\Users\Admin\Desktop\ChAT_analysis\` holds ChAT-soma and IO neuron recordings organized as `<date>/<cell|IO>/<recording>/`. Each recording folder has calcium pipeline outputs (`*_ch1_dFF.mat`, `*cpSAM_output.mat`, `*_ch1_meta.mat`, `ca_spike_data.mat`) and DLC breath outputs (`*DLC_*_breath_peak_data.mat` with fields `breath`, `t_breath`, `insp_onset_idx`, `findpeak_params`).

**Breath/calcium alignment convention (confirmed by user 2026-05-26):** the Basler breath camera is hardware frame-triggered 1:1 with the 2P imaging, so effective breath fps = imaging fps (~30 Hz). The `camera_fps = 50.70` stored in EVERY `breath_peak_data.mat` is a hardcoded placeholder, NOT measured — ignore it. The AVI container `FrameRate` (~98.6 Hz) is also bogus. `t_breath` is built from the bad 50.7 and is therefore on a wrong axis — do not use it for timing.

Alignment recipe: toss first 30 breath frames (same as calcium toss), align frame-by-frame, and if breath is longer (camera left running), truncate the tail. `T = min(len(breath)-30, len(calcium))`; clip both; shift `insp_onset_idx` by -30. No resampling. Evidence: 4/5 recordings have breath frames exactly equal to calcium frames (e.g. 3000=3000, 6000=6000); the one mismatch (0522/cell2/roi12: breath 8011 vs calcium 5970) is a forgot-to-stop tail.

**Why:** getting the time base right sets the entire coherence frequency axis. **How to apply:** when computing breath↔calcium coherence here, use [[project_breathing_is_motion]]/[[breath_combined_summary]] math (Chronux coherencyc, TW=4) but set fs_breath = detect_session_fps imaging rate, never the stored camera_fps. See [[feedback_detect_session_fps]].

**Verified ChAT neurons (user-designated, 2026-05-27):** only THREE ROIs are true ChAT; everything else is IO/non-ChAT. They are recording-folder-substring + ROI index: `roi5_7x_x-1200y200z-30_3000f_23lp` ROI1, `roi3_8x_x-1070y730z0_3000f_15lp` ROI1, `roi1_4x_x-900y700z-15_6000f_13lp` ROI1. Coherence inclusion = ROIs with >1 detected spike (minSpikes=2), not just the ifSpike flag.

**Analysis scripts (repo root, RZ_MATLAB):** `chat_breath_coherence_polar_260526.m` (breath-peak × Ca-spike coherence polar; ChAT colored w/ mag+phase CI, IO black circles no CI; auto band from breath PSD peak ± fwhm_factor×FWHM), `chat_dFF_vs_breath_260526.m` (per-recording dFF stack + breath overlay, meta-driven 3-line title), `chat_single_overlay_quick_260526.m` (single-folder dFF+breath and dFF+rawF overlays).
