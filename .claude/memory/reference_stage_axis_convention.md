---
name: reference_stage_axis_convention
description: "Image ROW maps to NEGATIVE stage y when placing FOVs by motorPosition — measured, not assumed; wrong sign silently costs half your cross-FOV cell matches"
metadata: 
  node_type: memory
  type: reference
  originSessionId: b2f03943-7192-4cdf-a370-82cf33218d5f
  modified: 2026-08-02T04:03:28.258Z
---

When placing FOVs into a common frame from `motorPosition`, the correct mapping is:

```
X = stage x        + (col - W/2) * px_um
Y = -1 * stage y   + (row - H/2) * px_um     <-- note the sign
```

**Measured, not assumed.** `roi_stitch_260801/roi_stitch_check_axes_260801.m` phase-correlates
overlapping FOV pairs and scores **all eight** symmetries of a square — the four sign
combinations of (col→x, row→y) *and* those four transposed (col→y, row→x). Testing only
the flips cannot detect a rotated frame; a 4-way version of this test gave a confident
wrong answer for exactly that reason. On `Vgat/0730` (15 pairs, overlap ≥0.5, |dz| ≤15 µm):

| rank | convention | median | pairs <10 µm |
|---|---|---|---|
| 1 | **col=+x row=−y (flipud)** | **9.1 µm** | 53% |
| 2 | col=+y row=−x (rot90cw) | 44.0 µm | 20% |
| 3 | col=+x row=+y (identity) | 44.7 µm | 13% |
| … | all four transposed | 44–76 µm | ≤20% |

Clean separation (4.8× to the runner-up), best pairs at 1–5 µm. **All transposed mappings
lose badly — there is no 90° rotation at the rig level.**

## Two distinct layers — do not conflate

| layer | transform | affects matching? |
|---|---|---|
| **rig** (scan mirrors → stage) | col→+x, **row→−y** | **YES** — sets where ROIs sit inside their own FOV |
| **prep** (stage → anatomy) | `AXMAP_ROT90 = [0 1; −1 0]`, lateral=+stage y, rostral=−stage x, sessions from 260721 on | **NO** — global rotation, distances invariant |

The second is [[project_rotated_mount_sessions]] and is handled in
`Ventral_surface_ROI_on_cartoon_260729.m`. Stitching and matching should be done in STAGE
coordinates and anatomy applied only for display.

**Why it matters — the failure is silent.** With the wrong sign, FOVs are placed tens
to hundreds of µm apart and nothing errors. On `Vgat/0730`, fixing it alone (no manual
stitching at all) moved cross-FOV matching from 787 cells / 78 seen >1× / max 3× to
691 cells / **143 seen >1× / max 6×**. The 6× matches the number of FOVs overlapping
in the z≈0 slab, i.e. the physically predicted ceiling was finally reached. Median
patch NCC also rose 0.50 → 0.64 *even though appearance never gates the matching* —
an independent metric confirming the fix.

The tell that something is wrong: you find yourself hand-nudging tiles by 50+ µm to
make a map close. That is compensating for an axis error, not stage error — real
stage error on this rig is ~2 µm (tie-point registration).

**Distinct from the rotated mount.** [[project_rotated_mount_sessions]] (everything from
260721 onward: −stage x = rostral, stage y = lateral) is a *global* transform for
reading anatomy. It cannot change which ROIs match, because distances are invariant
under it. This y-sign issue is separate and DOES change matching.

## Open: does the cartoon map need the same fix?

`Ventral_surface_ROI_on_cartoon_260729.m:303` places per-ROI positions with
`my_r = motor(2) + (c(2) - Hf/2) * px_um_fov` — row → **+**stage y, the convention this
measurement rules out. If it applies there too, every ROI is mirrored about its FOV's
horizontal centre line: ±150 µm for a 512 px FOV at 0.59 µm/px, ±379 µm at 1.48 µm/px.
FOV *centres* use `motor` directly and are unaffected, so site-level layout is fine — it
is the fine structure within each patch that would be mirrored. NOT YET VERIFIED on the
older sessions; the y-sign is presumably a rig property and so constant, but that is an
assumption until measured. Worth running the 8-way check on a non-rotated session
(e.g. Vglut2/0224 or 1124) before changing anything.

See [[project_roi_stitch_coordinate_first_matching]].

## CONFIRMED AGAIN on 260829, and the matcher was fixed (2026-08-30)

Re-measured independently on `260829_Sert-soma-g8s_vagotomized` (24 FOVs, 515 ROIs,
26 overlapping same-site FOV pairs with |dz| <= 15 um). Metric: fraction of ROIs whose
nearest cross-FOV neighbour is within 10 um.

| ySign | ROIs with a partner <10 um | median NN |
|---|---|---|
| +1 | 15.0 % | 27.1 um |
| **-1** | **43.7 %** | **11.7 um** |

2.92x, the same size of effect as the Vgat/0730 finding. Individual pairs are
unambiguous: `roi14_3x_z0` x `roi7_3x_z0` matched **96.6 %** of ROIs at a 6.3 um median.

**The matcher had it wrong and has been corrected.**
`roi_pair_morph_match_260727.m:181` places ROIs with
`y = motor(2) + ySign*(row - H/2)*px_um`, and `ySign` was **+1**. Correct value is
**-1** (the memory formula `Y = -stage_y + (row-H/2)*px` rearranged; a global Y flip
leaves distances invariant). Changed to -1 on 2026-08-30 in
`roi_pair_morph_match_260727.m` and `roi_pair_morph_quant.m`. The
`match_roi_morph\roi_pair_morph_*_260724.m` pair still carries +1.

**Why it silently destroys matches:** at 3x zoom `pixelSize_um = 0.5926`, so a 512-row
FOV has a 152 um half-height and the wrong sign displaces an ROI by up to **304 um** --
past `maxPairDist_um = 200`. Real repeats are dropped at the PREFILTER, before any
appearance scoring, so **no amount of curation can recover them**: curation can only
reject pairs the matcher proposed.

**Archive-wide implication (OPEN).** Every session matched with these scripts ran at
+1, including Sert/0828, Vglut2/0824 and every `roi_match_out_260724` output. Vgat/0730
is clean only because it went through the `roi_stitch_260801` path, which measures the
axes. So cross-FOV under-matching probably caps the whole archive. Sert/0829's curation
was kept as-is (RZ decision, 2026-08-30) rather than re-run -- it still found 108 cells
spanning >1 recording, but that is a floor, not the true number.

The diagnostic that produced the table is a scratchpad one-off (`ysign_check_260830.m`):
cluster FOVs into sites by `motorPosition`, then for each sign score cross-FOV
nearest-neighbour distances. It only discriminates where FOVs differ in stage y or in
zoom -- an identical mirror applied to two identical FOVs cancels out. Easy to rebuild.
