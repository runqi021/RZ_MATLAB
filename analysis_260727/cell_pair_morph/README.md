# cell_pair_morph — cross-FOV ROI matching (cell registration)

Copy of `..\..\match_roi_morph\`, retargeted to `cell_cfg_260727.m` and with two
work-losing bugs fixed in round 4. The originals are untouched.

The problem this solves: the same patch of tissue is often imaged repeatedly as
several FOVs, so one neuron shows up as N unrelated ROIs. This pipeline decides which
ROIs across recordings are the same cell.

## The only file you edit: `cell_cfg_260727.m`

```matlab
cfg.datasetPath = 'D:\260721_Sert_soma_G8s\phys';
cfg.outDirName  = 'roi_match_out_260724';   % where the curation already lives
```

`outDirName` is deliberately still `_260724`: that folder holds the curation you have
already done. Renaming it would orphan that work. Change it only for a new experiment.
It must match `cfg.cellOutDirName` in `..\coh_ca_breath\coh_cfg_260727.m`.

## Run order

| step | script | what you do |
|---|---|---|
| 1 | `roi_pair_morph_match_260727.m` | automatic. Zero-shift NCC between centroid-centred patches on a common µm grid, single-link grouping. Writes `roi_match_results.mat` + the patch cache. |
| 2 | `roi_review1_toss_260727.m` | click junk observations to toss |
| 3 | `roi_review2_assign_260727.m` | confirm or deny each proposed cell member |
| 4 | `roi_review3_merge_260727.m` | merge cells that are the same cell |
| 5 | `roi_review4_ungrouped_260727.m` | place the never-grouped leftovers |
| — | `cell_activity_summary_260727.m` | one figure per cell: medoid patch + stacked dF/F with breathing |

Rounds 2–5 all resume from `roi_match_curated.mat`; you can quit and come back.

## The shared state

Everything is one per-ROI vector, `grpOf`:

```
grpOf(r) > 0   member of that cell
grpOf(r) == 0  ungrouped, not tossed -> saved as its own 1x cell
grpOf(r) < 0   tossed
```

`homeOf` remembers which cell a tossed ROI came from, so un-tossing restores it.
Saving writes `roi_match_curated.mat` (the state) and `roi_match_curated.csv` (one row
per observation). Re-running the matcher (step 1) deliberately does **not** delete the
curation.

## Two bugs fixed in round 4

**1. Round 4 deleted the ungrouped ROIs it exists to handle.** Its `save_state` built
`curatedCellMembers` from `grpOf>0` only, so every ungrouped observation disappeared
from the cell list and the CSV as soon as it saved — 274 of 1290 on this dataset, i.e.
the round that is *about* ungrouped ROIs was the one that dropped them. Rounds 2 and 3
already appended them as 1× cells; round 4 now does the same. (The header comment had
always claimed this behaviour; only the code disagreed.)

**2. A skip in phase 2 was not remembered.** Skipping left `grpOf` at 0, which is
indistinguishable from never having been looked at, so every ROI you reviewed and chose
to keep reappeared in the phase-1 triage grid and the phase-2 queue on the next run.
Only tosses and assignments survived. Skips are now recorded in
`roi_match_r4_reviewed.mat`, and phases 1 and 2 both skip already-reviewed ROIs.

That file is kept **separate** from `roi_match_curated.mat` on purpose: rounds 1–3
rewrite the curated struct wholesale and would silently drop an extra field. Delete
`roi_match_r4_reviewed.mat` to re-review everything from scratch.

Round 4 also now prints its state on entry and a summary on exit, so it is obvious
whether it resumed or started fresh.

## Note on `repoRoot`

These copies sit two levels below the repo root instead of one, so every
`fileparts(fileparts(mfilename('fullpath')))` became a triple `fileparts`. Without that
they resolve `repoRoot` to `analysis_260727` and fail to find `detect_session_fps`.
