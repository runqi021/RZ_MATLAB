---
name: project_roi_crossfov_matching
description: Cross-FOV ROI cell-registration pipeline — roi_pair_morph_match_260724.m (match+group) + roi_match_review_260724.m (manual prune)
metadata: 
  node_type: memory
  type: project
  originSessionId: 1e48fc28-81f5-4afc-8886-ddee3a1a2d43
  modified: 2026-07-31T05:30:48.818Z
---

Cross-FOV ROI MATCHING = cell registration: the same patch is imaged REPEATEDLY as several FOVs (different zoom/O2-vs-roomAir/run), and repeats of one physical cell should collapse to "one cell seen N times". Built 2026-07-24, evolved from the older `roi_pair_morph_quant.m` (which only did an NCC-vs-IoU scatter; the matched-ROI gallery lived in its `.asv` autosave).

**Scripts (in subfolder `match_roi_morph/`; repoRoot = fileparts(fileparts(mfilename)) so detect_session_fps in repo root still resolves):**
- `roi_pair_morph_match_260724.m` — match + group + per-cell figures + save.
- MANUAL CURATION = 4 fast click-only round GUIs (the monolithic `roi_match_review_260724.m` was too slow/unusable and was DELETED). Display patch size = compare size (`viewHalf_um = patchHalf_um`) so figures/GUIs show exactly what NCC scored. All share ONE assignment vector `grpOf` per observation in `roi_match_curated.mat`: **>0 = cell id, 0 = ungrouped, -1 = tossed** (non-destructive: NA excludes from downstream, data kept). Each round is a FUNCTION file, loads the patch cache for instant startup:
  - `roi_review1_toss_260724.m` — per-cell montage, click a patch to toss/un-toss.
  - `roi_review2_assign_260724.m` — pool = TOSSED. candidate on top + its NEAREST existing cells (coord-shortlisted) shown by MEDOID below; click a cell to add, g=new group, x=toss, k=skip.
  - `roi_review3_merge_260724.m` — NCC-SUGGESTED merge: 2nd round of NCC BETWEEN cells (medoid vs medoid, only DISJOINT-FOV pairs within maxMergeDist_um=100), plots cell-cell NCC-vs-IoU scatter + NCC histogram, then steps through suggestions with NCC>mergeSuggestThr(0.5) high→low; y=merge(2nd→1st) n=skip. Re-run to catch transitive merges. (Replaced the old click-2-tiles grid.)
  - `roi_review4_ungrouped_260724.m` — never-grouped singletons, TWO phases: phase 1 = triage grid, click junk to toss ('d'=done); phase 2 = AUTO-SUGGEST (like round 3): each ungrouped ROI's best-matching existing cell (medoid NCC, disjoint-FOV, <=100um) proposed if NCC>assignThr(0.5); scatter+hist then step high→low, y=assign n=skip(stays 1x) x=toss. ~300 clicks → a few confirmations.
  - COMPLETENESS: save_state (all rounds) writes every NON-tossed ROI — grpOf>0 as their cells, and each remaining ungrouped grpOf==0 ROI as its OWN 1x cell (so skipped/unmatched ROIs are KEPT, per user's keep-all rule). grpOf itself keeps 0='pending' in the .mat so rounds re-run. Tossed grpOf<0 excluded.
  - assign GUIs show Kshow=20 nearest cells (was 12).
  - DISPLAY look: PER-ROI-PATCH normalization — each display patch stretched to its OWN dispClip=[1 99.9] pct then dispGamma=0.6 (via disp_norm_patch on the raw patch), so every cell pops regardless of FOV brightness. Bicubic resample (imresize default, antialias on downscale). Matching stays LINEAR (imgClip [5 99.5] per-FOV, NO gamma — would break NCC invariance); Umat from fov_img, display pC from fov_raw+disp_norm_patch. Baked into cache; get_patches fallback replicates. (Rejected per-FOV and global-across-FOV display norm.)
  - Cell "face" = MEDOID (most-central real member, highest mean NCC to rest), NOT a pixel-average (averaging blurs across plane drift = bias). Each writes roi_match_curated.{mat,csv}.
  - SPEED: matcher writes `roi_match_patches.mat` (per-ROI display patch+mask+Umat, single); all 3 GUIs load it via get_patches() for instant startup (fallback = reload FOVs). Per-cell PNGs are OFF by default (`savePerCellPng=false`) because exportgraphics-per-cell was the hang (hundreds × ~0.5s, NOT infinite); the `.mat`/CSV save BEFORE that loop.

**How matching works (settled):**
- Each ROI patch is cut CENTERED on its own Cellpose centroid, resampled to a common um grid, so a true match is aligned at ZERO shift. Score = `sum(unit_i .* unit_j)` == corr2 == NCC peak. Fully VECTORIZED in chunks (no normxcorr2, no parfor) — minutes→seconds. See [[feedback_vectorize_over_loops]].
- Coordinate (motorPosition) is ONLY a prefilter (`maxPairDist_um`,`maxPairDz_um`) to cut ~N^2 pairs; morphology (NCC) is THE metric, never coordinate. Within-FOV pairs excluded (`fov(ia)~=fov(ib)`); one-ROI-per-FOV enforced in grouping.
- Grouping = greedy SINGLE-LINK agglomeration, best NCC first, with one-per-FOV must-not-link. Kept single-link on purpose: repeats sit at slightly different PLANES so a cell's appearance is a TRAJECTORY (adjacent-plane repeats match strong, distant weak) — average/complete-link would wrongly fragment true drifting members. Group-average NCC is a REVIEW FLAG, not an auto-reject.
- Confidence per cell = all-pairs NCC inside the finished group (`Umat(mem,:)*Umat(mem,:)'`): cellMinNCC (weakest internal pair = honest score), cellMedNCC, memberNCC (each obs' mean NCC to rest). Weak members (r<nccMergeThr) flagged red.

**Standard settings:** maxPairDist_um=100, maxPairDz_um=40, patchHalf_um~18 (context size = the knob that fixes dense-cluster confusion, e.g. roi6/7), targetUmPerPx=1.5 (optical floor; DOWN-samples all FOVs → denoise + comparable; 0.5 upsamples→noise), nccMergeThr=0.60 (set off the NCC-vs-IoU scatter), cellSortBy='confidence' (shakiest first).

## Threshold retune 2026-07-30 (Vglut2/1124/pFN) — and why NCC alone cannot gate this

Params in `roi_pair_morph_match_260727.m` were loosened after a visibly-correct pair was
being rejected: `maxPairDist_um` 300→**400**, `maxPairDz_um` 50→**100**, `nccMergeThr`
0.60→**0.40** (`patchHalf_um`=30, `targetUmPerPx`=0.5). Result on that dataset: 48 ROIs →
35 cells, 11 → **13** cross-recording merges. **One of the two new merges is right and one is
wrong**, and NCC cannot tell them apart:

| new merge | NCC | position residual | verdict |
|---|---|---|---|
| recA#5 ↔ recB#12 | 0.573 | **3.0 µm** | correct — the pair that motivated the change |
| recA#23 ↔ recB#3 | 0.449 | **236 µm** | **FALSE** — two similar silhouettes, different neighbourhoods |

**The position residual is a clean discriminator where NCC is not.** The two recordings sit
at a fixed frame offset (measured from the confirmed pairs: dx +3.6, dy +175.9 µm); subtract
it and every one of the 12 true pairs lands at **≤ 3.4 µm** while the false merge is at 236 µm
— a 70× gap with nothing in between. NCC meanwhile puts a true pair (0.573) *below* a false
one is not quite right, but they overlap badly: 0.449 false vs 0.573 true, against accepted
pairs spanning 0.449–0.892. **Recommended: gate on offset-corrected residual (<~10 µm) IN
ADDITION to NCC**, rather than lowering nccMergeThr further. Coordinates are currently
"prefilter only" — this is the one place they carry more information than morphology.

**Why a true pair scores low: zero-shift NCC + centroid drag.** Patches are cut centred on the
Cellpose centroid and scored at zero shift (no translation search). When Cellpose segments the
same neuron with a different extent — recA#5 257 px vs recB#12 324 px, the recB mask swallowing
more of the process — the centroid moves, the two patches end up misregistered by several px,
and the correlation collapses. So *segmentation disagreement, not appearance change*, is what
sinks these scores. The principled fix is a ±3–5 px translation search in the NCC; lowering the
threshold treats the symptom and admits false merges (as it just did).

**Diagnostic recipe** (used above, worth repeating): estimate the inter-recording frame offset
from the *confirmed* pairs only (median dx/dy), then report every merge's residual after removing
it. True pairs cluster at the ~2 µm noise floor; anything ≫ that is a false merge regardless of
its NCC. `match.pairs` (i, j, ncc, iou, dist) holds every scored pair, so a rejected pair can
always be re-examined without re-running the matcher.

**Beware the truncated `fov_name`.** `match.fov_name` for a folder containing a dot comes back
cut at it — `roi5_1400-1230-0_x4.4_15lp_6000f_00001` → `roi5_1400-1230-0_x4` (`fileparts` on a
folder path; same bug class as [[project_260728_vglut2_vessel_depth]]). `match.fov_folder` keeps
the full path and is correct. Matching is unaffected, but the `(rec_name, roi_index)` key that
curation portability depends on is NOT unique-safe under truncation — see
[[feedback_roi_curation_port]].

**Intensity:** per-FOV percentile stretch [5 99.5] for display only; NCC's own mean-sub+normalize makes matching invariant to LINEAR brightness/contrast (handles O2/roomAir/power). Do NOT add gamma — nonlinear, breaks that invariance + amplifies the noise the 1.5um grid suppressed. More context = bigger patchHalf, not gamma.

**Downstream: cell_activity_summary_260726.m** (in match_roi_morph/) — pools activity per matched CELL into one figure: LEFT = medoid avg-proj patch (per-patch gamma + outline), RIGHT = one stacked row per recording = that recording's dF/F (black, column roi.roi(r) of *_ch1_dFF.mat) + breathing overlaid (blue, from breath_peak_pc1.mat `breath`, length-aligned to dFF by dropping front toss frames = numel(breath)-numel(dFF)). fps via detect_session_fps. Cells rebuilt from the per-ROI group vector (curated grpOf: >0=cell, 0=ungrouped-untossed→own 1x cell, <0=tossed excluded; else results match.grp) — NOT from curatedCellMembers, so ALL non-tossed ROIs included regardless of how curated was saved. minCellSize=1 (include 1x singletons = full inventory; 2 = pooled multi only). Saves to <datasetPath>\cell_id_test\ (wiped each run). Params: maxCells (test cap 30; Inf=all), minCellSize, leftHalf_um. dFF drawn raw (noisy); breath scaled into each row's dFF range. Verified: 455 total cells (166 were size>=2-ish; ~289 singletons now included).

**Rerun-safety:** matcher clean-slate deletes only ITS OWN products (roi_match_results.mat, cells.csv, patches.mat, cell_*.png, pairs_top_*.png, scatter) — it explicitly KEEPS roi_match_curated.{mat,csv} (human curation). GUIs resume from curated only if numel(grpOf)==nRoiTot (else warn + start fresh). To restart curation, delete roi_match_curated.mat by hand.

**Outputs (in datasetPath/roi_match_out_260724/, own products regenerated each run):** roi_match_results.mat (`match` struct: params, roi table, pairs, grp, cellMembers, memberNCC, cellMin/MedNCC, fov_folder), roi_match_cells.csv, cell_###_NNx_minMM.png (one figure per cell), pairs_top_*.png, ncc_vs_iou_scatter.png. Reviewer writes roi_match_curated.{mat,csv} (non-destructive: keeps drops + original).

**Data note (260721_Sert_soma_G8s):** 27 cpSAM FOVs, 1 each; roi7 imaged ~8x (1.6/1.7/1.8x × O2/roomAir) so 12x cells are REAL not dup. cpSAM_output.mat has F/maskL/avgPath/t/opt — NO cellpose flow/gradient field (would need re-running cellpose; and flows are derived from mask ≈ redundant with IoU anyway).

Possible next: z-consistency check (weak link OK if large plane gap, suspicious if same-z); constellation/neighbour-arrangement term for dense fields. Related: [[feedback_no_inplace_edits]].
