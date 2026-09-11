---
name: project-structural-stack-260910
description: "structural_stack_260910/ two-stage pipeline: raw SI z-stack -> group-averaged 1 um isotropic volume -> falloff-style projection + true 3-D figure, one per channel"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0c9f4688-9d18-4f93-9896-d7dc4ce25892
  modified: 2026-09-11T00:14:51.003Z
---

**`structural_stack_260910\`** — generic structural 3-D stack processing + visualisation,
built 2026-09-10 on `C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\map_FN_fov1_00001`.
Group-agnostic and channel-count-agnostic: **everything geometric is read from the SI
header**, so a 1-channel and a 2-channel stack take the same path with no edit.

| file | role |
|---|---|
| `struct_stack_prep_260910.m` | read once, group-average frames per z, resample to 1 um isotropic, write one TIFF per channel + `_stackmeta.mat` |
| `struct_stack_render_260910.m` | one figure per channel: 100 um depth-bin MIPs, coronal, sagittal, **true 3-D oblique ray-cast** |
| `run_map_FN_fov1_260910.m` | driver; edit `TIF` to point at another stack |

Stage 1 is the only expensive step and it CACHES, so replotting is free. Delete
`*_avg1um.tif` + `*_stackmeta.mat` to force a re-read.

## 260909 numbers (map_FN_fov1_00001)
512x512 zoom 2 -> **0.8889 um/px**; 531 slices, **1 um step**, 30 frames/slice;
`channelSave [1 3]` -> **ch1 = ChAT-G8m somata, ch3 = vessel dye**; motor z **+30 -> -500**;
`pzAdjust` 15->68 %, Lz -350.655, zPowerReference +30. 31 860 pages = 15.63 GiB.
Read **288 s**; output **[455 455 531]** per channel, 0.23 GB each. Ray-cast 36-89 s.
FN motoneurons resolve to ~350 um; ch3 penetrating vessels to ~400.

## Decisions that are easy to get wrong on a re-run
- **Pages interleave channels PER FRAME**: `page = ((slice-1)*fpsl + frame-1)*nCh + chIdx`.
  A per-slice stride still returns correctly-sized images and never errors — it silently
  averages the OTHER label. `chIdx` is the position in `channelSave`; outputs are named by
  the SI channel NUMBER.
- **Pixel size from `imagingFovUm`, not `PIX_BASE/zoom`** — the FOV form already carries
  `scanAngleMultiplier` and rotation.
- **Walk pages with `nextDirectory`**, never `setDirectory(t,k)`: quadratic over 31 860 IFDs.
- **Do NOT sanity-check page count as `fileBytes/frameBytes`** — each page carries ~2.4 kB
  of IFD/tags (72 MB total here), so the naive division invents 144 phantom pages and cries
  truncation on a complete file. Compare against expected pixel bytes instead.
- **int16 output, not uint16.** ScanImage digitises signed and the dark level sits near 0;
  clamping negatives away floors the noise and biases every later background estimate.
- **No depth correction.** The stack is `pzAdjust`-ramped, so attenuation is already
  compensated during acquisition — `exp(z/l)` on top double-counts it. `'Norm','perplane'`
  exists if the deep half still reads dark; off by default.
- **`ZSurface` is an ASSUMPTION, not a measurement.** depth = ZSurface - motor z, default 0.
  Stage 1 *reports* the (top-bot) half-rise so it can be checked and shifts nothing:
  ch1 **+29 um** (meaningless — a sparse soma label has no surface edge), ch3 **-12 um**.
- Mount `'auto'` = `rot90` for 260909 (>= 2026-07-21), so coronal = max along cols and the
  ML axis is labelled `lateral -> medial`. That DIRECTION is inherited from 260728/260804 —
  see [[project_rotated_mount_sessions]], it is not to be carried across sessions unasked.
- Layout: **3 tile rows**, side views span the lower two. Equal rows leave a broad empty
  band because the x-y panels are square and a side view is as deep as the stack.

## Running it
`matlab -batch` in a DETACHED background process, not the MATLAB MCP:
`run_matlab_file` is killed at the 120 s mark when the call is backgrounded, which
terminated the first read at slice 250 ("Operation terminated by user").

All three FOVs (`map_FN_fov1/2/3_00001`) are filed in their own subfolders and fully
prepped + rendered as of 2026-09-10 (reads 288 / 157 / 114 s). The session also holds a
separate `map\` tile folder, which is a different job.

See [[project_falloff_optical_penetration]], [[project_falloff_render_and_box_260901]],
[[reference_matlab_opengl_export_stubs]], [[feedback_detect_session_fps]].

## Manual stitcher for the 260909 strip (added 2026-09-10)
`stack_stitch_gui_260910.m` + `stack_stitch_apply_260910.m`, same folder. Built because
`260818-map-registration\stitch_manual_gui_fast_260818.m` assumes a col/row LATTICE of
avg-projection images, and this is 3 VOLUMES with no lattice where z matters.

**The strip:** fov1/fov2/fov3 sit at stage x **-1727.12 / -1330.11 / -933.11**, identical
stage y (1164.43) and z (30). Step **397.0 um**, FOV 455.13 -> **58.1 um overlap (12.8%)**.
At 1 um isotropic that is 397 px apart, so 1 px = 1 um everywhere in the GUI.

- **Seam window** (`v`, or the button; opens on start): for the selected tile, per
  overlapping neighbour, the overlap region from each tile in grey plus a RED/GREEN overlay
  in X-Y and again in X-Z. The overlay is normalised PER IMAGE while the grey panels keep
  the shared window — on a shared window the brighter tile paints the whole overlap its own
  colour and says nothing about alignment. Yellow = registered.
- Typed **nudge step box, floored at 1 um** (the voxel is a micron; a smaller step is
  discarded by round() and the tile appears stuck). Focus is parked on a hidden control
  after typing so the arrow keys go back to the tiles, not the caret.
- GUI shows an **X-Y mosaic AND an X-Z mosaic**; the X-Z one is the point — two tiles can
  look registered in X-Y and still be a slab apart in depth. Arrows nudge x/y, `[`/`]`
  nudge z, click or 1..9 selects, r / R resets, scroll+right-drag zoom/pan the panel under
  the cursor, s saves.
- **No flip control, on purpose.** Stage mode owns the sign
  ([[reference_stage_axis_convention]]); a flip here mirrors the whole mosaic.
- **Z nudges are WHOLE PLANES** — a plane is a micron, and interpolating to chase a
  sub-micron offset blurs every tile to fix a shift under the axial PSF.
- Save -> `<root>\matlab_stitch\stack_stitch_coords.mat`; apply renders any channel
  headless with a linear feather blend. **Weights are 2-D (x,y only)**, so the denominator
  is one plane and the numerator single [455 x 1249 x 531] ~ 1.2 GB fits in RAM.

## LATERAL SCALE MISMATCH, stage vs ScanImage FOV -- **+6.4%** (2026-09-10)
Two independent hand-curated stitches on the same rig/session, both saying the header
UNDERSTATES the pixel size. **Trust the map number: 144 pairs beats 2.**

**`map\` 9x9 mosaic, 144 adjacent pairs (`map_curated_step_260910.m`) -- the good estimate:**
image step **418.75 +/- 3.33 px in x**, **426.88 +/- 4.31 px in y**, for 400 um of stage
travel -> **0.9550 / 0.9374 um/px**, pooled **0.9462 +/- 0.0126** against the header's
0.8889 = **+6.4%**. Implied `objectiveResolution` **41.4**, not 38.9.
Curated overlap: **93 +/- 3 px in x, 85 +/- 4 px in y** (88 / 80 um at the measured scale).
**x and y differ by 1.9%** (0.955 vs 0.937) -- plausibly a real resonant-vs-galvo axis
difference (`forceSquarePixels` fixes the pixel COUNT, not the optical scale), but curation
error is correlated across pairs so do not quote the SEM as if the pairs were independent.

**The fov1/2/3 strip gave ~+10%, from 2 pairs** (steps 358, 360 resampled px for 397 um,
i.e. ~0.985 um/px in original pixels). It is 4% off the map value and is the weaker
measurement -- 2 pairs, judged on a nominal-58-px overlap that is really ~97 px.
**An earlier note here claimed "~11%, ACCUMULATES" as if settled; it was not.**

**Consequences if the map number holds:** `_avg1um.tif` volumes are ~6% coarser laterally
than 1 um/px, every 100 um lateral scale bar is ~6% short, and z (from the motor) is
unaffected -- so the "1 um isotropic" volumes are really ~0.95 x 0.95 x 1 um. Still
UNRESOLVED as a rig calibration; a stage-step test would settle it independently.

**Two traps that cost real numbers here:**
- **`stitch_manual_gui_fast_260818` stores `x = base + manX`** -- the nudge is ALREADY in
  `x`, and `manX` is only a record. Adding them double-counts and inflated the step
  432.5 vs 418.75 px, i.e. the answer itself. `stack_stitch_gui_260910` uses the OPPOSITE
  convention (base and nudge separate, must be added). Guard with the stored `outW`/`outH`:
  the placement must span exactly that canvas.
- **That dataset has `flipY = true`**, so the signed row step is NEGATIVE and
  `tile - step` produced an overlap of 960 px on a 512 px tile. Use magnitudes.
- **Nudges are quantised to the GUI's move step (5 px here)**, so the per-pair spread is
  mostly quantisation, not measurement scatter.

**A phase-correlation estimate of +19 px was WRONG and the failure is instructive:** it
searched inside the NOMINAL 58 px overlap while the true overlap is 98 px, so the real peak
was partly outside the window (peaks were 0.06 / 0.03). Do not size a correlation window
from the placement you are trying to verify.

**The brain surface is tilted across the strip:** ch3 (top-bot) half-rise runs
**-12 / -19 / -27 um** for fov1/2/3, ~15 um of drop over 794 um. That is what the per-tile
z nudge is for.

**Two bugs worth remembering:** `S.tileH(i)=` builds a ROW while offsets are COLUMNS, and
`ox + ww` on mixed orientations broadcasts to n-by-n instead of erroring, so the canvas size
comes out a vector; and `axis image` on the 531-deep X-Z panel makes it height-limited, so
it renders as a thumbnail unless the axes gets nearly as much canvas as the X-Y view.

## Shift-null validation of the manual stitch (`stitch_shift_null_260910.m`, 2026-09-10)
Asks whether RZ's hand alignment beats a random one. Per seam: displace tile B by a random
(dx, dy), score, repeat 500x; read the curated placement against that null. 3 tiles = 2 seams.

**Metric, designed here: YF = sum(min(a,b)) / sum(max(a,b))** — the seam view's "yellow",
summed: a soft IoU of the two overlap images. **Background subtraction is what makes it
work**: flat background sits in min() and max() equally, so on raw data YF climbs toward 1
for ANY shift. Each patch gets its own 25th pct removed, /(99th-25th), clipped to [0,1] —
per patch, so a brightness difference between tiles cannot masquerade as structure.
Pearson NCC is reported alongside purely to catch YF misbehaving.

**Two design points that matter more than the metric:**
- **The comparison window must be FIXED.** Recomputing the overlap per shift changes its
  size AND content, so the score partly measures how many pixels were compared. Erode the
  curated overlap by ShiftMax, hold tile A fixed, move only B's sampling. Erosion by exactly
  ShiftMax guarantees B always has real data.
- **Shifts are drawn from an ANNULUS, |d| in [3,30] um per axis with random sign** (RZ's
  call). A 1 um draw is not a misalignment, it is the curated placement plus noise, and such
  draws pile up in the null's upper tail and make a good alignment look unremarkable.
  Integer microns: the voxel is 1 um, and sub-pixel would need interpolation that blurs the
  patch and confounds the score with a smoothing. Note both axes are forced past ShiftMin,
  so a PURE x or y misalignment is never sampled — set ShiftMin 0 if that is the question.

**Result: the manual stitch is unambiguous — 0/500 draws beat it on either seam, either
channel, either metric** (p < 1/501 = 0.002 each).

| ch | seam | YF curated | YF null | z | NCC curated |
|---|---|---|---|---|---|
| 1 | 1-2 | 0.749 | 0.532 +/- 0.029 | +7.6 | +0.912 |
| 1 | 2-3 | 0.587 | 0.334 +/- 0.026 | +9.8 | +0.747 |
| 3 | 1-2 | 0.811 | 0.551 +/- 0.096 | +2.7 | +0.950 |
| 3 | 2-3 | 0.474 | 0.233 +/- 0.071 | +3.4 | +0.799 |

**ch3's z is far lower than ch1's despite a BIGGER absolute gap** — its null SD is 3x wider
(0.096 vs 0.029) because the vessel channel is sparse and high-contrast, so a random shift
sometimes lands vessel-on-vessel by luck. z is the wrong summary for a sparse label; the
0/500 separation is the honest one.
**Do not read absolute YF across seams as alignment quality** — seam 2-3 scores lower than
1-2 on both channels, but its z is HIGHER on ch1. Absolute YF tracks how much shared
structure the overlap contains, not how well it is registered.

Output: `<coordDir>\shift_null_ch<N>\shift_null_ch<N>.png/.pdf/.csv/.mat`.
