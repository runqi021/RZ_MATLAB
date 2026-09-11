---
name: project_representative_traces_260527
description: "representative_traces_avgproj_260527.m — hand-picked SST/Vgat representative ROIs, breath+dFF overlay + gamma avg proj"
metadata: 
  node_type: memory
  type: project
  originSessionId: 959165bd-a94a-4536-be5e-f8543c8230bf
---

`representative_traces_avgproj_260527.m` (created 2026-05-27) — figure script for hand-picked representative ROIs across many dffQC FOVs. Per FOV: PLAIN avg projection (done like `avg_proj_scalebar_260527.m`: crop crop_um, clip [0.5 99.9], gamma 0.6, burned 50µm scalebar — NO ROI outline, NO text label) + that ROI's FULL-recording dF/F trace. Combined summary = one row per FOV, PLUS N standalone per-input PNGs (`repr_##_<stem>_ROI##.png`) collected in the out folder. Outputs to `D:\batch_dffQC_test_260325\representative_260527\`. NOTE: breath overlay, ROI outline, and ROI text label were ALL explicitly REMOVED on user request (none of these cells are breath-tuned) — keep avg proj plain + traces dFF-only.

**Curated selections (user's "representative" cells, 2026-05-27):**
- SST (260322_sst_soma_g8s): `maybe breathing\6x_-820-1070-z20_00001` ROI4; `maybe breathing\6x_-743-840-z30_00001` ROI2; `random\7x_-850-1032-z30_3000f_00001` ROI2 (NO breath data — dFF-only).
- Vgat (260323_vgat_g8s): `730-930-z25_9x_3000f_00001` ROI3; `760-750-z30_6x_3000f_00001` ROI3; `-500-990-z-30_8x_3000f_00001` ROI1.

**Gotcha baked in:** dF/F in dffQC came from the SINGLE-MC cpSAM (verified F_roi==cpSAM F, maxdiff 0), so maskL ← `*_preproc_MC_cpSAM_output.mat` and avg ← `*_preproc_MC_AVG_for_CP.tif` (NOT the `_MC_MC` double-MC versions — those would misalign outlines).
