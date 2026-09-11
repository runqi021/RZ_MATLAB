---
name: feedback_vectorize_over_loops
description: Performance lesson — vectorize batched math instead of per-item loops / parfor; question whether an expensive op is even doing work
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1e48fc28-81f5-4afc-8886-ddee3a1a2d43
  modified: 2026-07-25T01:58:07.332Z
---

The user values code that runs "faster and smarter and better" — treat obvious slowness as a bug to fix, not a given.

Concrete win (roi_pair_morph_match_260724.m, cross-FOV ROI matching): pairing 83k candidate pairs went from **minutes → seconds** by replacing ~83k FFT-based `normxcorr2` calls (wrapped in `parfor`) with a **vectorized zero-shift correlation** — precompute one unit-norm, mean-subtracted patch vector per ROI (`Umat`), then score all pairs in memory-capped chunks as `sum(Umat(A,:).*Umat(B,:),2)`. Identity used: `sum(unit_i .* unit_j) == corr2(patch_i, patch_j) == normxcorr2 peak at zero shift`.

**Why:** the heavy tool was doing almost no work — patches are extracted CENTERED on each ROI's own centroid, so a true match is already aligned in both patches; `searchHalf`≈`patchHalf` meant normxcorr2 searched only ±1 µm. A giant FFT machine for a ±1 µm shift is pure waste. Also `parfor` cost a 30–60 s `parpool` startup + broadcasting a ~50 MB cache to 12 workers — more overhead than the light per-iteration compute saved. Vectorization beat parallelism outright.

**How to apply (general playbook):**
1. `tic/toc` around each phase FIRST to find the real bottleneck — don't optimize by guessing (here the user confirmed "it's the pairing part", not saving).
2. Ask "is this expensive op actually doing work?" If a search window / tolerance is tiny, the general tool (normxcorr2, imregister, optimization) may collapse to a one-liner (corr2, a dot product).
3. Precompute per-item ONCE, never per-pair; then batch the pairwise step as chunked matrix ops (guard memory with a chunk size like `CH=10000`).
4. Reach for `parfor` only when per-iteration work is heavy; for light work, vectorize instead — pool startup + broadcast can dominate.
5. Reusing one hidden figure (`clf`) beats creating N figures when exporting many PNGs.

Related: [[feedback_no_inplace_edits]] (new analyses get a NEW .m), [[feedback_detect_session_fps]].
