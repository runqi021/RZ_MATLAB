---
name: project-cell-link-unapplied-curation
description: Two sessions had curated cross-FOV matches that cell_link never joined; merges now live in cell_merge_overrides_260814.m; matcher fov_name is fileparts-truncated
metadata: 
  node_type: memory
  type: project
  originSessionId: 7abf9dd7-af19-40fb-a912-10c9692c4000
  modified: 2026-08-15T05:42:38.657Z
---

Three cell pairs spotted by eye in the per-cell summaries (2026-08-14) turned out to have **two different causes**, and only one was a matcher failure.

**Cause 1 — curation done, never joined.** `Vglut2/0224/pFN` and `Vglut2/1124/pFN` both have `roi_match_curated.mat` but **no `cell_link.mat`**. The pairing existed all along (the Vglut2/0224 pair was curated group 6) and was simply never joined into the registry. Their matcher output sits one level deeper, under a `pFN/` subfolder, which is the likely reason `cell_link` did not find them. Joining those two curations recovers **5 real cell-cell merges** (0224 groups 6 and 14; 1124 groups 2, 3, 4).

**Cause 2 — genuine matcher miss.** `Vglut2/0810` (cells 48+52) and `Sst/0807` (218+220) DO have `cell_link.mat`. Both ROIs of each pair were present as observations and were left ungrouped. Position rules them out: the 0807 pair sits **161 µm apart**, inside the range where a loose gate previously made a 236 µm FALSE merge. Those merges rest on the visual match, not on position.

**Where the corrections live:** `analysis_260806\cell_merge_overrides_260814.m` (7 groups) + `apply_cell_merges_260814.m` (union-find) + `pooled_obs_260814.m`. Keyed on **(recName, roi)** parsed from `OBS.label`, never on a cell index — cell numbering is regenerated on every rebuild. Transcribed rather than fixed upstream because re-running `cell_link` renumbers every cell and would invalidate all figure filenames. A merged cell takes the **lowest** of its constituent numbers, so only the merged ones change name. Re-running `cell_link` later makes these no-ops ("already merged"), not double counts. Result: 276 → 269 cells, 326 observations conserved.

**TRAP — `fov_name` in `roi_match_results.mat` is `fileparts`-truncated.** Everything after the first dot is eaten as a file extension:

```
roi1_z5_2.5x_3000f_lp15_00001  ->  roi1_z5_2
fov1_1.7x_y1650x1320_z30_...   ->  fov1_1
```

This **collides**: in `Vglut2/0728` three distinct z-planes (z30/z35/z40) all become `fov1_1`; in `Vgat/0730` a `cell1` and a `quiet` recording both become `roi1_1`. The matcher itself is fine (it indexes by FOV number and keeps `fov_folder`), but **never join on `fov_name`** — use `match.fov_folder`, or `OBS.label` from the registry, which keeps full names.

**Still open:** four curated groups pair an active cell with an ROI that failed the active gate in its other recording (0224 grps 9/16/17 → cells 26, 28, 30; 1124 grp 9 → cell 177). Pooling them would not change the cell count but WOULD shift those cells' events, logZ and permutation p. Left out pending a decision.

Related: [[project-roi-crossfov-matching]], [[feedback-roi-curation-port]], [[project-analysis-260727-cell-pooled-coherence]]
