---
name: feedback-roi-review4-save-bugs
description: roi_review4_ungrouped round 4 silently dropped ungrouped ROIs on save and never remembered skips — fixed only in the analysis_260727 copy
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-07-28T01:13:25.998Z
---

`roi_review4_ungrouped_260724.m` (the ORIGINAL, in `match_roi_morph\`) has two
work-losing bugs. They are fixed only in the `analysis_260727\cell_pair_morph\`
copy — the original still has them.

1. Its `save_state` built `curatedCellMembers` from `grpOf>0` only, so every
   ungrouped ROI vanished from the cell list and `roi_match_curated.csv` the moment
   round 4 saved. Rounds 2 and 3 already appended them as 1x cells; round 4 — the
   round that is *about* ungrouped ROIs — deleted them.
2. A phase-2 "skip" left `grpOf==0`, indistinguishable from never-looked-at, so every
   reviewed-and-kept ROI reappeared in the triage grid and suggestion queue next run.
   Fixed via a separate `roi_match_r4_reviewed.mat` (separate because rounds 1-3
   rewrite the curated struct wholesale and would drop an extra field).

**Why:** the user asked "is this not reloading my saved results?" — the reload was
fine; the save was destroying the work, so it *looked* like a load bug.

**How to apply:** when a resumable review GUI seems not to resume, check what its
save writes before checking what its load reads. And diagnose which script last
wrote a .mat by comparing its FIELD SET against each script's save code — that is
how it emerged that no round-4 save had ever landed.

See [[project_roi_crossfov_matching]], [[project_analysis_260727_cell_pooled_coherence]].
