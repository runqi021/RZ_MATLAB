---
name: TIFF metadata only in raw files
description: Only the original raw ScanImage TIFF has metadata (zoom, fps, etc). Processed TIFFs (MC, preproc) do NOT. Always read metadata from _meta.mat or the raw TIFF, never from processed TIFFs.
type: feedback
---

Only the most raw/original ScanImage TIFF has embedded metadata (Software/ImageDescription fields with zoom factor, frame rate, etc). All processed TIFFs (_preproc.tif, _MC.tif, _MC_MC.tif) do NOT carry this metadata — they are written by saveastiff/run_rigid_mc without ScanImage headers.

**Why:** The QC pipeline (Batch_dffQC_260325.m) saves a `_meta.mat` file per session specifically because processed TIFFs lose metadata. This was a recurring source of bugs — code trying to parse zoom/fps from MC TIFFs gets nothing and falls back to wrong defaults (e.g., 1x zoom instead of 10x, making ROI crops way too small).

**How to apply:** When you need zoom factor, fps, or any ScanImage metadata:
1. First try `*_meta.mat` in the session folder (has `pixelSize_um`, `zoomFactor`, `fps`, etc.)
2. If no _meta.mat, look for the raw source TIFF (filename matches parent folder, no _preproc/_MC suffix)
3. NEVER try to parse metadata from processed TIFFs — it won't be there
