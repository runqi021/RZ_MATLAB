---
name: project_roi_stitch_coordinate_first_matching
description: "roi_stitch_260801/ — standalone stitch-first ROI grouping; position is the metric, appearance only reports"
metadata: 
  node_type: memory
  type: project
  originSessionId: b2f03943-7192-4cdf-a370-82cf33218d5f
  modified: 2026-08-02T04:15:37.162Z
---

`C:\Users\Admin\Desktop\RZ_MATLAB\roi_stitch_260801\` — a second, **independent** way to
decide which ROIs across recordings are the same neuron. Built 2026-08-01 at the user's
request: *"i will work on the fov and try to stitch them first, then we have the new
coordinate, then match roi based on that. make sure this is standing alone so no cross
contamination."*

**Inverts the existing matcher.** `roi_pair_morph_match_260724/260727` uses coordinates as a
prefilter and decides by patch NCC. This one stitches first, then treats "same place" as
"same neuron". Appearance is computed and written to the CSV but **never gates** — which is
what made it a usable independent check (see [[reference_stage_axis_convention]]).

**Standalone by construction:** own config (`stitch_cfg_260801`), own output dir
(`roi_stitch_out_260801/`), own filenames. The scan explicitly skips `roi_match_out*`
folders. Neither pipeline can damage the other's curation.

## Files
- `stitch_cfg_260801.m` — config; `('set', path)` to retarget
- `roi_stitch_scan_260801.m` — dataset → `fov_table.mat`
- `roi_stitch_map_260801.m` — **the aligner to use**: ONE XY map, all tiles, depth ignored;
  click a tile, arrows move 1 µm, `g` auto-registers it by phase correlation, `shift+g` all
- `roi_stitch_align_260801.m` — pairwise fallback for one ambiguous pair
- `roi_stitch_match_260801.m` — group by position
- `roi_stitch_overview_260801.m` — overview figures
- `roi_stitch_check_axes_260801.m` — **run this first on any new rig/session**

## The three gates
- lateral ≤12 µm — set from the *registration residual* (~2 µm), NOT soma size
- |dz| ≤25 µm — Z spans 330 µm over 21 planes on Vgat/0730; without it, cells stacked in
  depth coincide laterally and merge
- unambiguous — no match if the runner-up is nearly as close

Plus a structural rule: **a group may not contain two ROIs from the same recording.** True by
construction, and it is what blocks the single-link chaining that produced Vgat cell 2 in the
appearance matcher.

## Two coordinate transforms, applied in order (`roi_stitch_scan_260801`)
1. **rig**: `cfg.stageYSign` — applies to the ROI OFFSET ONLY, never to the FOV centre.
   Flipping both is a global reflection: relative geometry unchanged, so it silently does
   nothing. See [[reference_stage_axis_convention]].
2. **prep**: `cfg.AXMAP_ROT90 = [0 1; -1 0]` — lateral = +stage y, rostral = −stage x, for
   sessions from 260721 on. Rigid, so matching is invariant; half-extents SWAP (a 360×512 px
   FOV reports 457×650 µm, long axis rostral).

`roi_stitch_orient_260801` carries image PIXELS through the same transform
(`rot90(A.',2)` rotated / `flipud(A)` not). Every consumer calls it, so map, overview and
auto-register cannot drift apart.

**The invariance check is the test that catches all of this**: a pure rotation must leave
matching *bit-identical*. 691 cells / 143 seen >1× / max 6× before and after. When a
"rotation" changes the numbers, the transform is wrong.

## Traps hit while building it
- Applying the y-sign to centre AND offset → reverted matching to the wrong-axis numbers
  (787/78/3×). Caught only by re-running the matcher, not by reading the code.
- Rotating coordinates but not the images → tiles drawn in the wrong orientation while ROI
  markers sat correctly; auto-registration fell from 15/23 tiles to 7/23.
- Auto-register with no depth gate composited 15 FOVs spanning 330 µm into one reference;
  phase correlation returned confident wrong peaks of 40–82 µm. `AUTOREG_DZ = 15` →
  median 13.8 → 2.7 µm.
- Overlap saturates at 1.00 for many FOVs, so ordering by overlap alone falls through to FOV
  index — a dz=19 µm pair got presented before a dz=2 µm one. Score by `ov × (1 − dz/max)`.
- Single-linkage on z chains a finely sampled depth series into one slab (merged 80 µm / 12
  FOVs). Use lateral site × fixed-width depth band instead.
- `shortname({x}){1}` is a syntax error in R2021b — cannot index a call result.
- Rebuilding the axes per keypress destroys zoom/pan; build once, then only update `CData`
  and `XData`/`YData`.
