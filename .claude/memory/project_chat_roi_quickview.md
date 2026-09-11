---
name: project_chat_roi_quickview
description: chat_roi_quickview_260527.m — single-ROI breath/calcium quick view; heatmap sort modes
metadata: 
  node_type: memory
  type: project
  originSessionId: 959165bd-a94a-4536-be5e-f8543c8230bf
---

`chat_roi_quickview_260527.m` (created 2026-05-27) — single-FOV, single-ROI quick view distilled from [[project_chat_analysis_breath_alignment]]'s `chat_breath_coherence_polar_260526.m`. Set `folderPath` + `roi`, run. One figure (tiledlayout 3x2): dF/F+breath overlay, PSD (breath + dF/F'), breath×dF/F coherence, breath-triggered dF/F heatmap, triggered average ±SD. Saves `quickview_ROI##.png/.pdf` into the folder. Reuses Chronux + detect_session_fps; same toss-nDrop-then-truncate breath alignment.

**Heatmap sort (the parent script's Fig 4 sorts the SAME way — NOT dt-nearest):** default `sortMode='postmean'` = mean dF/F in post-inspiration window (tau>=0), descending. Added options: `'dt'` (interval to NEXT breath onset, ascending = "dt nearest"), `'none'` (chronological). The parent `chat_breath_coherence_polar_260526.m` only ever did postmean (`post=mean(E(:,tau>=0),2)`).
