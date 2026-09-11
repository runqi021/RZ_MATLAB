---
name: project-falloff-render-and-box-260901
description: "WT-vs-shiverer vessel figures (cube view, oblique ray-cast render, boxed falloff) and the traps found building them - above-pia extension, saturation in MIPs, per-plane normalisation, and the E-drive session the scatter was silently missing"
metadata: 
  node_type: memory
  type: project
  originSessionId: 41ccab4b-4965-4c62-809e-6991e34f1bb8
  modified: 2026-09-01T07:43:09.852Z
---

Four scripts in `falloff-analysis-260805\`, all comparing **Vglut2 roi1 (WT)** against
**shiverer_dbh roi2**, all on the SAME 200x200 um sub-volume so the figures describe one
object:

| script | output in `C:\fall-off\` |
|---|---|
| `falloff_cubeview_260831.m` | `falloff_cube_<tag>` — single block, two MIP faces |
| `falloff_cube_pair_260831.m` | `falloff_cube_pair_260831` — WT + shiverer side by side |
| `falloff_render3d_260831.m` | `falloff_render3d_260831` — oblique ray-cast, real parallax |
| `falloff_metrics_box_260901.m` | `falloff_metrics_box_260901` (+ `.mat`) — boxed top5-bot5 |

Boxes, read off the axes of the existing `falloff_proj_<tag>.png`:
WT `ML 100-300, AP 50-250`; shiverer `ML 0-200, AP 100-300`.

## The above-pia extension — do NOT re-read 13 GB
`falloff_projections_260805.m:364` builds its grid as `(0 : dzNat : max)`, so the cached
volume STARTS AT THE PIA and every above-surface plane is discarded.
`falloff_extend_above_260831.m` prepends them by reading only the SHALLOWEST stack's first
slices — **300 pages in 2.4 s (WT), 1500 in 3.6 s (shiverer)** — and writes a separate
`falloff_proj_<tag>_vol_from<D>.mat` so the original cache is untouched.
**How far up is set by acquisition:** WT `roi1_z-10-100` reaches only **-10 um**, shiverer
`roi2_2x_z-50-150` reaches **-50**. +20 was never acquired for the WT at all.
`pick_vol` ERRORS if the extended cache is missing rather than falling back — a figure
captioned "from -10" that silently started at 0 would be wrong with nothing to show it.

## Saturation is real, and a MIP amplifies it
`satFrac` reaches **7.5%** (WT `lp23`) and **9.6%** (WT `lp21`) of pixels pinned at the
32767 digitizer ceiling — a clot or dye pooling. A MAX projection lets those voxels win
EVERY ray they lie on, so a few bad voxels smear into bright bands. Fix: keep the K
brightest samples per ray and average (`TOP_K = 4`).
**Specify K as a COUNT, not a percentile:** there are only ~175-210 ray steps, so 99.5 pct
rounds to K=1 and silently gives back the plain MAX the accumulator exists to avoid.

## Brightness normalisation — `NORM_MODE` in the render
The depth correction makes the two incomparable: at 200 um it multiplies the WT by **446x**
(l=24.6) and the shiverer by **26x** (l=61.4), so the WT's deep half becomes bright haze.
- `'depthcorr'` the pipeline default; per-dataset display window, NOT comparable
- `'perplane'` **RZ's choice** — no depth correction, each z plane stretched to its own
  percentiles, ONE window shared across datasets. Removes the brightness ramp but NOT the
  SNR collapse (a noise plane normalises to grain, a vessel plane to vessels), so
  penetration still reads, uncoupled from exposure.
- `'snr'` per plane `(I-median)/(1.4826*MAD)`, units of noise SD — the quantitative version
Caveat on `'perplane'`: deep WT grain is stretched to full contrast and can read as
structure. `'snr'` with a fixed window keeps noise near 0 sigma.

**Frangi/vesselness was NOT used** and should not be: already rejected as "looks fake"
([[project_subvolume_render]]), and it is contrast-normalised by construction, so it would
enhance the WT's deep grain into tubes and make the genotypes CONVERGE.

## Boxed falloff — the metrics figure
Metrics are read off the ASSEMBLED volume, so each depth appears exactly once from the
round that owns it (build_volume already split overlaps at their midpoint) — that is what
removes the overlapping dots the whole-field plot has. `pctPlane` carries that round's
power. Definitions unchanged from `falloff_260804.m` (`top`/`bot` = mean of brightest/
dimmest 5%, `diff = top-bot`, `diff_n = diff/mW^2`).
Two fits: GLOBAL on `diff_n` (red, log panel, drawn only over the fit range), and PER POWER
ROUND with **l HELD FIXED, amplitude free** (orange, raw linear panel, joined into a
sawtooth by vertical risers at the round boundaries). The per-round R^2 is a TEST of the
global l, not a fit. Legend is in **mW, not laser %** — mW ~ pct^1.9, so 11%->48% is
14.1->238.5 mW, a 17x change from a 4.4x setpoint change.

**Boxed l is shorter than whole-field** because the box excludes surface vessels that keep
`top` high and flatten the early decay: WT **23.5** vs 24.6 (fit 10-200), shiverer **60.1**
vs 61.4 (fit 10-500). Sensitivity: shiverer 58.2 fitting to 450.
**WT fit range matters.** Its stored `FIT_DEPTH` is [10 150] and that looks deliberate —
extending to 200 pulls in the 238.5 mW round, whose per-round R^2 is **-0.34** (worse than
a flat line) because it has reached the noise floor. That is what moves l 21.5 -> 23.5.

## Real z0 is NOT auto-detected anywhere
`DEPTH_MODE = 'filename'`; depth comes from the filename z token. Measured pia (peak of the
7 um smoothed top-bot): **WT +3 um, shiverer -2 to -3 um** — a few um off in OPPOSITE
directions, so the two blocks are not referenced to the same landmark and "-10 um above the
pia" is really -13 for the WT and -7 for the shiverer. `PIA_MODE` exists in
`falloff_metrics_box_260901.m` but is **'off'**; nothing is shifted. A constant shift moves
the fitted AMPLITUDE, not the slope, so l is unaffected either way.
The shiverer's top-bot PLATEAUS over ~15 um near the surface (bright vessel plexus pins the
top 5%), so its peak locates the pia only to about +/-5 um; the WT's rise is sharp.

## The scatter was silently missing a session
`falloff_summary_scatter_260805.m` listed `C:\260810_vglut2-soma-g8s\...` but the data is on
**`E:\260810_vglut2-soma-g8s\roi1_stack_z-5-70_30f_00001`**. It warned and dropped the row,
giving **n=4** WT. Fixed to E:. Its real fit was already on disk: **l = 20.9, R^2 0.997,
fit 10-70** — the shortest WT FOV, at least partly because a 70 um stack only supports a
10-65 um fit, a third of the range behind every other WT FOV.
The scatter now READS `falloff_metrics_box_260901.mat` and overrides the whole-field l for
those two FOVs (`BOX_MAT`), drawing them with a black marker edge because **the pool is then
mixed** — 2 boxed FOVs among whole-field ones. Set `BOX_MAT = ''` to revert.

**Current pooled result:** wildtype **n=5, l = 23.8 +/- 1.8 um (EAL 47.7)**;
shiverer **n=6, l = 74.5 +/- 10.5 um (EAL 148.9)**. No overlap (WT max 25.5, shi min 60.1).
Supersedes the n=4 numbers in [[project_falloff_optical_penetration]].

## MATLAB geometry traps hit while building these
- **A figure canvas is capped at the screen.** Asking for 1700x950 on a 1536x960 screen
  returns 1537.8x893, so a design applied unscaled makes "squares" 271x282 and puts anything
  positioned in raw pixels ~90 px out. Read the size BACK and scale by one factor.
- **`colorbar` shrinks its peer axes** even with `PositionConstraint='innerposition'` —
  re-apply the axes Position AFTER placing the colourbar.
- **`axis image` resizes the AXES** to the image aspect; pad the image to square instead if
  the box must stay square.
- A **Windows path in a TeX-interpreted title fails to render at all** (backslashes are
  control sequences, `%` starts a comment) — use `Interpreter','none'`.

See [[project_falloff_optical_penetration]], [[reference_laser_power_lut]],
[[feedback_falloff_analysis_style]], [[reference_matlab_opengl_export_stubs]].
