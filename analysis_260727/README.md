# analysis_260727 — breathing modulation of calcium events

**Read `PROGRESS_REPORT.md` first** — it has the current results, every bug found,
and what is still open.

Three pipelines, all self-contained copies; the originals in `..\match_roi_morph\`
and `..\coh_ca_breath_260722\` are untouched and still work.

```
analysis_260727\
    run_analysis_260727.m   driver: runs the whole chain on any dataset
    breath_svd\             UPSTREAM: .avi -> breathing trace + peak/onset triggers
    cell_pair_morph\        which ROIs in different recordings are the SAME neuron
    coh_ca_breath\          coherence + the cell-identity join
    breath_time\            PRIMARY: inspiration-triggered PETH in absolute time
    phase_rayleigh\         RETIRED — kept and runnable, not part of the chain
```

`breath_svd\` sits before everything else: it produces the `breath_pc1.mat`,
`breath_peak_pc1.mat` and `breath_insp_start_pc1.mat` that
`analysis_preflight_260727.m` looks for. It is a manual step (you draw an ROI and
curate events), so it is not in the driver — but nothing downstream runs without
it. **`breath_svd_pc1.m` is the current per-video trace producer**; see that
folder's README for when to use `breath_fixedmetric_*` instead.

The point of joining cell identity to the physiology: a neuron imaged in four
recordings used to appear as four unrelated ROIs, each with its own weak estimate.
Now it is one cell with one estimate built from all four — and the four are kept
alongside it, so pooling can always be checked against its parts.

---

## Just run it

```matlab
run_analysis_260727()                                        % default dataset
run_analysis_260727('C:\...\Ventral_surface_summary','Ventral')
```

Everything lands in `<dataset>\analysis_260727\`; nothing outside it is touched.
Cell matching is **optional** — with no `roi_match_results.mat` the chain falls back
to identity grouping (every ROI its own cell) and still completes, with pooling off.

The driver runs, in order:

```
1  coherence_polar_general_260727   per-ROI coherence (secondary measure)
2  cell_link_260727                 join cell identity to recordings + AUDIT
3  cell_pool_260727                 lossless per-cell container
4  breath_time_overlay_260727       breath diagnostics: window + trigger QC
5  breath_time_peth_260727          PRIMARY: per-cell PETH + shuffle test
```

Step 2 computes nothing and only checks the join — read its output before trusting
anything after it.

## Manual steps (curation, not in the driver)

```
breath_svd\                                       (required; see its README)
  breath_svd_pc1.m                    avi -> breath_pc1.mat        (ROI per video)
  breathing_peak_gui_pc1.m            -> breath_peak_pc1.mat       TRIGGER
  breathing_trough_gui_pc1.m          -> breath_insp_start_pc1.mat TRIGGER

cell_pair_morph\                                  (skip if already curated)
  roi_pair_morph_match_260727.m       auto-match ROIs across FOVs
  roi_review1_toss_260727.m           drop junk observations
  roi_review2_assign_260727.m         confirm/deny proposed members
  roi_review3_merge_260727.m          merge cells that are the same cell
  roi_review4_ungrouped_260727.m      place the never-grouped leftovers
     -> roi_match_curated.mat
```

## Why the primary analysis is time, not phase or coherence

- **Coherence** asks whether respiration and calcium share power at a stable
  *frequency*. It weakens whenever breathing rate drifts, even for a cell that fires
  at exactly the same point of every breath.
- **Cycle-interpolated phase** fixes that but is not uniformly occupied in time —
  only 15% of each cycle falls in `[0, π)`, a 12× occupancy imbalance, so an
  unmodulated cell looks modulated until an ECDF correction is applied.
- **Absolute time from inspiration onset** has neither problem: every trigger
  contributes to every time bin, so the null is flat with no correction at all.

All three were run and they agree (see `PROGRESS_REPORT.md` §2), which is why the
simplest one is the one reported.

## The two config files

Paths live in exactly two places, and nowhere else:

| file | sets |
|---|---|
| `cell_pair_morph\cell_cfg_260727.m` | `datasetPath`, the matcher's output folder |
| `coh_ca_breath\coh_cfg_260727.m` | `rootPath`, `genotype`, and where the curation is |

`cfg.cellOutDirName` in the second must equal `cfg.outDirName` in the first. Both
currently point at `roi_match_out_260724`, which is where your existing curation
(1290 ROIs → 165 multi-recording cells + 274 singletons, 16 tossed) already lives.

## Renamed on purpose

`coh_cfg` → `coh_cfg_260727`, and every copied script calls the new name. Two files
called `coh_cfg.m` on the MATLAB path shadow each other depending on which folder was
added last, which would silently point this pipeline at the 260722 config. The output
folder is likewise `coherence_polar_sert_260727`, so nothing here can overwrite the
260722 results.

`repoRoot` was also corrected in every copy: these files sit two levels below the repo
root, not one, so the original `fileparts(fileparts(...))` would have resolved to
`analysis_260727` and failed to find `detect_session_fps`, Chronux and `+helper`.

## What gets written where

Into `<rootPath>\coherence_polar_sert_260727\`:

| file | from | holds |
|---|---|---|
| `coherence_polar_data.mat` | step 6 | per-ROI coherence, now with an explicit `(recName, roiIdx)` join key |
| `cell_link.mat` / `.csv` | step 7 | every ROI observation ↔ cell id ↔ recording, plus the audit |
| `cell_link_recordings.csv` | step 7 | per-recording ROI-count cross-checks |
| `cell_pool.mat` | step 8 | the lossless container: breath per recording, spikes per observation |
| `cell_coherence_pooled.mat` / `.csv` / `.png` / `.pdf` | step 9 | one coherence per cell + its parts |

## Two bugs fixed in the copies

Both were losing work, both are in `roi_review4_ungrouped_260727.m`:

1. **Round 4 deleted the ungrouped ROIs it exists to handle.** Its `save_state` looped
   over `grpOf>0` only, so every ungrouped observation vanished from
   `curatedCellMembers` and from `roi_match_curated.csv` the moment it saved — 274 of
   1290 observations, on this dataset. Rounds 2 and 3 already emitted them as 1× cells;
   round 4 now matches.
2. **A skip in phase 2 was not remembered.** It left `grpOf` at 0, so every ROI you
   looked at and chose to keep came back in the triage grid and in the suggestion queue
   on the next run. Skips are now recorded in `roi_match_r4_reviewed.mat` — a separate
   file, so re-running rounds 1–3 cannot clobber it. Delete that file to re-review
   everything.

Also fixed, in `coherence_polar_general_260727.m`: `sig_rois.csv` split its label on
`#`, a character the label never contains, so the `fov` and `roi` columns both came out
as the entire label string.
