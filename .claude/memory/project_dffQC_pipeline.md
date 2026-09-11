---
name: dffqc-pipeline-naming-and-architecture
description: "The single-Z pipeline script is called \"dffQC\" — runs both rigid and non-rigid MC, then cpSAM + dFF for both, with QC plots saved."
metadata: 
  node_type: memory
  type: project
  originSessionId: d0d6f257-9f5d-4554-aa44-3425c396829e
  modified: 2026-07-23T05:27:29.433Z
---

"dffQC" is the shorthand name for the single-Z quality-control pipeline script.

**Architecture:**
- `Batch_dffQC_260325.m` (repo root) is the current batch driver. It runs 2-pass **rigid** MC (`run_rigid_mc`) -> cpSAM -> dFF per file, resumable (skips steps whose outputs exist). Input `folderPath`/`masterFolder`; outputs land beside each TIFF.
- **Active MC folder = `NoRMCorre-master/`** — `Batch_dffQC_260325.m` line ~32 does `addpath(NoRMCorre-master)`, so `NoRMCorre-master\run_rigid_mc.m` is the one that runs. (The old `MC_spnr_RZ_v1/` note is stale; today only two `run_rigid_mc.m` exist: `NoRMCorre-master\` [active] and `NoRMCorre_modified_RZ_v1\` [off-path spare].)

**MC big-file / OOM fixes (2026-07-22, applied to BOTH run_rigid_mc copies):** long recordings (12000-18000 frames) failed two ways — see [[feedback_tiff_metadata]] sibling context.
- **Classic-TIFF 4 GB cap:** the MC output was written with a frame-by-frame `imwrite(...,'Compression','none')` loop = classic TIFF only -> `TIFFAppendToStrip: Maximum TIFF file size exceeded` past 4 GB. Fixed: write via `saveastiff` (which auto-BigTIFFs >4GB; `saveastiff.m` line ~224 `s.bytes>2^32-1`). Preproc writes were already safe (Batch uses saveastiff). Downstream `imread`/cpSAM read BigTIFF fine.
- **Out-of-memory during rigid MC — the REAL culprit was `motion_metrics`, not the save.** Traceback: `run_rigid_mc` -> `motion_metrics` line 79 `cY = corr(Yr,mYr)` where `Yr=reshape(Y,[],T)` is `[262144 x T]`; `corr` builds several whole-movie copies (crop + z-score internals) -> ~19GB x several -> OOM/hang on long movies (12000/18000 fr). The QC metrics are diagnostic-only. Fixes in `run_rigid_mc` metrics block: **subsample frames when T large** (`mstride=ceil(T/3000)`, feed <=3000 frames into corr, QC plot x-axis = `tcorr=1:mstride:T`), compute display range from a cheap linear-index subsample (not `quantile(Y(:))` which full-copies), and **`clear Y` BEFORE corr()**. Saved MC movie unaffected (full-res, all frames). Also (pre-existing this session): removed the wasted auto-started 12-worker parpool (rigid uses `use_parallel=false` so it never used it), and write output in 1000-frame chunks via `saveastiff(append,big)` + local `clip_cast_block` (no whole-movie cast duplicate). Peak after all fixes ~28GB@18000fr. Hard floor is still Y+M1 (~38GB@18000fr) coexisting at normcorre return -> that survives via pagefile swap (18000fr registration ~3800s = slow but completes). If #10 is too slow, next lever = stream input from the preproc TIFF via NoRMCorre `output_type='tiff'` (Y+M1 never fully load) — not yet done (breaks in-memory toss/baseline). Machine: 64GB RAM, 20 cores; free RAM was ~30-35GB during batch (not fully restarted).

**How to apply:** When user says "dffQC", they mean this pipeline. MC functions come from `NoRMCorre-master/`.
