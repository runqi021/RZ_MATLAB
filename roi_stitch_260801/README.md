# roi_stitch_260801 — stitch-first ROI grouping

A second, **independent** way to decide which ROIs across recordings are the same
neuron. Stitch the FOVs first, then group ROIs that land in the same place.

This is the inverse of the existing matcher (`roi_pair_morph_match_260724/260727`),
which uses coordinates only as a prefilter and decides identity by patch NCC.

## Standalone — no cross-contamination

| | appearance matcher | this pipeline |
|---|---|---|
| config | `cell_cfg_260727` | `stitch_cfg_260801` |
| output dir | `roi_match_out_260727/` | `roi_stitch_out_260801/` |
| products | `roi_match_results.mat`, `roi_match_curated.mat`, `cell_link.mat` | `fov_table.mat`, `fov_offsets.mat`, `roi_stitch_results.mat` |

Separate config functions, so retargeting one cannot move the other. This pipeline
never reads `roi_match_*` (the scan explicitly excludes those folders) and never
writes outside its own output dir. Both can run on the same dataset and be
compared; neither can damage the other's curation.

## Run order

```matlab
addpath('C:\Users\Admin\Desktop\RZ_MATLAB\roi_stitch_260801');
stitch_cfg_260801('set', 'D:\Ventral_surface_summary\Vgat\0730');

roi_stitch_check_axes_260801;  % 0. FIRST on any new rig/session -- see below
roi_stitch_scan_260801();      % 1. read dataset -> fov_table.mat
roi_stitch_map_260801();       % 2. stitch: ONE XY map, click a tile, arrow it
roi_stitch_match_260801();     % 3. group ROIs by position -> results + CSV
roi_stitch_overview_260801();  % 4. overview figures (optional, any time)
```

## Step 0 is not optional: the axis convention

Placing a FOV needs to know where an image column and row land in STAGE space.
There are eight possibilities (the symmetries of a square) and the wrong one fails
**silently** — FOVs land tens to hundreds of µm apart and nothing errors. The tell
is finding yourself hand-nudging tiles by 50+ µm to make a map close.

`roi_stitch_check_axes_260801` phase-correlates overlapping pairs and scores all
eight. Measured on `Vgat/0730` (15 pairs):

| rank | convention | median | pairs <10 µm |
|---|---|---|---|
| 1 | **col=+x row=−y** | **9.1 µm** | 53% |
| 2 | col=+y row=−x | 44.0 µm | 20% |
| 3 | col=+x row=+y | 44.7 µm | 13% |

Set the winner via `cfg.stageYSign`. Testing only the four *flips* cannot detect a
rotated frame — an earlier 4-way version of this test gave a confident wrong answer.

### Two coordinate layers — do not conflate them

| layer | transform | affects matching? |
|---|---|---|
| **rig** (scan mirrors → stage) | col→+x, row→−y | **YES** — sets where ROIs sit inside their own FOV |
| **prep** (stage → anatomy) | `[0 1; −1 0]`, lateral=+stage y, rostral=−stage x, sessions from 260721 on | **NO** — global rotation, distances invariant |

The rotated-mount convention is the *second* one and belongs to
`Ventral_surface_ROI_on_cartoon_260729.m`. Stitch and match in **stage** coordinates;
apply anatomy only for display.

Step 3 runs without step 2, on raw stage coordinates. That is a useful baseline
but not the intended result — see below.

### Two aligners

| | `roi_stitch_map_260801` **(default)** | `roi_stitch_align_260801` |
|---|---|---|
| view | all FOVs in one XY projection | one pair at a time |
| depth | not used at all | dz gates + cell/vessel regimes |
| step | 1 µm (shift ×10) | 0.5 µm (shift ×10) |
| use it when | normal case — you want to see the whole map | a pair is ambiguous and you want it isolated with a phase-correlation seed and a live match count |

Both write the same `fov_offsets.mat`, so they are interchangeable and can be
mixed. The map tool is the one to reach for.

## Results so far on `Vgat/0730` (23 FOVs, 871 ROIs)

| stage | cells | seen >1× | largest group | median min-NCC |
|---|---|---|---|---|
| wrong axis convention, no stitch | 787 | 78 | 3× | 0.50 |
| **axis fixed**, no stitch | 691 | 143 | **6×** | 0.64 |
| axis fixed + auto-registration | 688 | **145** | 6× | 0.64 |

Six FOVs overlap at z ≈ −12…+15 µm, so a neuron there should be recoverable from
up to 6 recordings. Under the wrong axis convention the ceiling was 3×; with it
fixed, 6× is reached — the physically predicted maximum.

Note that **median patch NCC rose 0.50 → 0.64 even though appearance never gates
the matching**. An independent metric improving is the strongest confirmation
available that the fix is real and not just a change in bookkeeping.

Auto-registration adds little here (2 cells) because the corrections it finds are
already only ~2.7 µm median — the frame was the whole problem, not the stage.

## The three gates

| gate | default | why |
|---|---|---|
| lateral distance | 12 µm | the identity test. Set from the **registration residual** (~2 µm achieved by tie-point registration on this rig), not from soma size. 30+ µm starts absorbing genuine neighbours. |
| \|dz\| | 25 µm | Z spans 330 µm over 21 distinct planes here. Without this, cells stacked in depth coincide laterally and merge. |
| unambiguous | ratio 0.7 | if the runner-up candidate is nearly as close as the winner, no match. In a dense field the nearest neighbour is not by itself evidence of identity. |

Plus a hard structural constraint in the merge step: **a group may not contain two
ROIs from the same recording.** One neuron is seen at most once per FOV, so this
is true by construction — and it is what prevents the single-link chaining that
walked the appearance matcher from small round blobs to large elongated somata
(`Vgat` cell 2).

## Appearance is reported, never used to gate

`cell_minNCC` is computed for every accepted group and written to the CSV, on
zoom-normalised 41×41 patches. It is a **check on the result**, not the metric.
Low values are flagged for review, not dropped — NCC has a known failure mode
here (segmentation disagreement moves the centroid, patches misregister, a true
pair scored 0.57 against a false pair at 0.45), which is the reason this pipeline
exists.

## Coordinate frame

Raw `motorPosition`, no anatomical axis map and no vessel calibration. Within one
animal the frame only has to be **internally consistent**, and every extra
transform is another thing that can be wrong. Whatever mapping the cartoon applies
is a separate, later, global step.

## Outputs

- `fov_table.mat` — per-FOV motor position, pixel size, ROI centroids in µm
- `fov_offsets.mat` — `fovTable(fov, tx_um, ty_um, placed, z_um)`; resumable
- `roi_stitch_results.mat` — per-ROI `grp`, coordinates; per-cell size / minNCC / spread
- `roi_stitch_cells.csv` — flat table, sorted by group size
