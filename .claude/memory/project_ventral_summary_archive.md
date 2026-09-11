---
name: project-ventral-summary-archive
description: "Ventral_surface_summary layout, its 'active cell' criterion, and the cartoon map's hard dependency on a hand-clicked vessel calibration"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-08-09T06:35:08.298Z
---

**MOVED 2026-08: the archive is now `D:\Ventral_surface_summary`** — the `C:\Users\Admin\Desktop` copy no longer exists (checked 2026-08-17). All `_260816` scripts point at D:. It is the curated cross-animal archive. See [[project_ventral_selected_cells_260816]].

**THE ARCHIVE IS SELF-CONTAINED as of 2026-07-31 — keep it that way.** Acquisition
sessions moved D: → E:, and E: is an external drive that gets unplugged mid-session. Each
session's grouping was therefore copied INTO the archive beside its recordings
(`<Genotype>\<MMDD>\roi_match_out_*\` and `<Genotype>\<MMDD>\cell_pooled\`), and
`cell_link_sources` in `Ventral_surface_polar_coh_vs_rayleigh_260729.m` now resolves with
`fullfile(rootPath, ...)` instead of absolute drive paths. **If a summary script ever needs
D: or E: again that is a regression** — a summary that silently loses cell identities when a
drive is absent is worse than one that refuses to run. Work log: `analysis_260727\PROGRESS_260731.md`.

**Contents 2026-07-31: 13 genotype/date blocks, 52 recordings.** ChAT 0124/0521/0522/0523/0524,
Sert 0721 (6 sites/10), Sst 0322/0330, Vgat 0323 (3/3) and **Vgat 0730 (3 sites/9, added
2026-07-31)**, Vglut2 0224 (8) / 0728 (3) / 1124 (8).

**`archive_into_summary_260729('template')` guesses site from the FILENAME PREFIX and that is
routinely wrong** — on Vgat/0730 it put all ten spike-bearing `roi1_*` into one `cell1`.
Site must come from `motorPosition` clustered at `mergeSiteUm = 100`; that gave 3 sites
~230 µm apart. Note `cell1` there was formed by single-link chaining and its extremes are
145 µm apart (defensible: those five span zoom 1.4×–3.2× and a 1.4× FOV contains all of them).
Vgat/0730 also has 7 `roi2_*` recordings with **no breath data at all** — `cpSAM` done, but
`breath_svd_pc1` + both GUIs never run, so they cannot enter until that is done.

**Vgat/0730 breathes at ~0.4 Hz — 5x slower than every other session in the archive.**
Median cycle **2.567 s** vs ~0.45 s for a 2.2 Hz Vglut2. Verified independently of the
pipeline from the hand-curated peaks (fps / median inter-peak interval = 0.37/0.38/0.38/0.68 Hz),
so it is the animal, not a bug. Phase pooling across genotypes is still valid (phase is
normalised per cycle) but any ABSOLUTE-TIME cross-genotype comparison is now comparing 2.6 s
cycles against 0.45 s ones — and check the TW=4 coherence bandwidth is adequate at 0.35 Hz
before quoting Vgat coherence beside the others. `roi1_3x_x1350y850_z265` runs at 0.68 Hz,
double the rest, and supplies 3 of the session's 4 significant coherence ROIs.
**Also real, do not "fix":** two Vgat/0730 recordings ran at `scanFrameRate_raw` 42.084 and
47.147 Hz against 30 for the rest; TIF page counts match the filenames.

Chain run 2026-07-31: 871 obs → 709 cells (62 multi-recording, 93 tossed), JOIN IS CLEAN,
0 excitation-significant of 12 testable cells. `cell_link_sources` in the polar script now
lists Sert/0721 + Vglut2/0728 + Vgat/0730; **0224 and 1124 have curated matchers but no
`cell_link`** (their chains were never run), so pooling is off for them.

**Layout**: `<Genotype>\<MMDD>\<site>\<recording>\`. Site is `cell1..cellN` (one
imaging site each) or an anatomical label (`IO`, `7N`). pFN recordings are filed as
`cellN`, NOT a `pFN` folder — see `Vglut2\0224\cell1..cell5`. MAO is excluded (part of
IO, and its cells are not active). Genotypes: ChAT, Sst, Vgat, Vglut2, and Sert added
2026-07-29 (`Sert\0721`, baseline condition only).

**ACTIVE CELL = `nnz(spike_train>0) > 5`** i.e. >=6 events, raw count. Defined in
`spike_trigger_dFF.m:34` and `temporal_phase_perROI.m:56`. The coherence foundation's
`minSpikes = 2` is a separate, weaker gate. Match this for any POOLED/population
analysis; the per-cell PETH test needs >=20 events instead, because a 6-event cell has
no power in a shuffle test and only inflates the BH denominator.

**The cartoon map assumes ONE stage zero for every animal.**
`Ventral_surface_ROI_on_cartoon_*.m` applies a single `stage_to_canvas` affine across
all genotypes, so every session must have zeroed its stage on the same anatomical
landmark. That is why `ro1_2x_x1018y1130_18lp_3000f_00001` (260721_Sert, first FOV, a
different zero) is in `coh_cfg_260727.cfg.excludeRecordings`. It does NOT assume one
stage *orientation* any more — see [[project_rotated_mount_sessions]].

**Cartoon map calibration — RECOVERED, no longer a blocker.** It needs
`ROI_on_vessel_data.mat` (`ax`,`ay` stage→canvas affine; `bx`,`by` hand-clicked basilar
branch = anatomical origin; `theta_rad` vessel-axis rotation; `um_per_px_stitch`).
That file was missing, and its producer `Ventral_surface_ROI_on_vessel_260528.m` is
INTERACTIVE. But a prior successful run had saved all six variables into its own output
`Desktop\ROI_on_cartoon_data.mat`, so `recover_vessel_calib_260729.m` reads them back and
writes `ChAT_MAP\matlab_stitch\ROI_on_vessel_data.mat` — recovery, not re-derivation, so
the map stays in the same anatomical frame as before. Recovered values: `ax=[1.03797
1629.81]`, `ay=[-1.06443 3581.34]`, `bx,by=1634.53,3633.19`, `theta_rad=0.083686`
(4.79°), `0.8889 µm/px`. It refuses to overwrite an existing file; re-click only if the
stitch changes.

**One patch per SITE, not per recording.** Several sites were re-imaged at a second Z, and
`ro7`/`roi7` were recorded twice under near-identical names. Each stamped its own
0.85-alpha patch at the same place, which read as one dark smear. The map now merges FOV
centres within `mergeSiteUm = 100` (same session + group), draws the recording with the
most ROIs, and keeps every ROI in `R` — two Z planes are two sets of cells, only the
rendering was duplicated. Merged 5: Sert cell1 (42 µm), Sert cell7 (10 µm),
Vglut2/0728 cell2 (28 µm), Vglut2/1124 IO (0 µm), and Sert `ro3_4x_x700y-1233` (73 µm).
That last one is a **mislabelled archive site**: by stage position it is a 4× zoom-in on
the cell1 site, while its supposed sibling `ro3_2.5x` is 312 µm away. `cell7`'s pair being
10 µm apart confirms filing both as one site was right.

**A MASK IS NOT A CELL — never report summed Cellpose masks as a cell count.** The same
neuron re-imaged at a second Z or zoom yields one mask per recording. Resolve identity
through the matcher's `cell_link.mat`, keyed on the portable `(rec_name, maskL_label)` pair
(`obsT`), and honour its `status`: **`tossed` → `cell_id` NaN, rejected in curation, not a
cell, drop it**; `ungrouped` → a real cell seen once, counts as one. Beware
`numel(unique(cell_id))` — MATLAB `unique` returns every NaN separately, so tossed rows
each inflate the count by one. Sert on the 10 archive recordings:
**402 masks → 13 tossed → 389 → 296 cells → 48 active → 34 testable → 18 modulated**
(the 34 independently matches what the PETH reported).

**CORRECTED 2026-07-30: only Sert/0721 has a matcher run.** An earlier note here claimed
Vglut2/0728 did too — it does not. There is no `roi_match*` anywhere under
`E:\260728_vglut2_soma-g8s`; its `analysis_260727` chain ran 2026-07-29 on the identity
fallback (one cell per ROI, pooling off), which the driver does silently by design.
Matching would only affect Vglut2 `cell2`, whose two recordings (`roi2_4x` z−15 and z−20)
are 30 µm apart laterally and 5.6 µm in Z, so they almost certainly share neurons — 6 ROIs
each. Everything else in the archive is one recording per site, where matching is a no-op.

**Map significance source**: `coherence_polar_svd_260729/` (not `_260606`, which predated
Sert and Vglut2/0728). Produced by `Ventral_surface_coherence_polar_svd_260729.m` — same
analysis and parameters, `scan_dirs`/`groups` extended. 368 ROIs / 14 significant
(9 Vglut2, 4 Sert, 1 ChAT) vs 297 / 10 before; `confC = 0.827` unchanged. The 4 Sert cells
sit at phase +1.68…+2.50 rad, i.e. pre-peak inspiration, same range as the Vglut2 pFN
cells. Note the newer PETH result cannot be substituted here: the map's right hemisphere is
coloured by coherence *phase*, and the PETH gives *latency*.

**Tools (2026-07-29)**: `archive_into_summary_260729.m` (template → edit csv → copy;
copy not move, skips existing), `population_hist_260729.m` (fixed 3 s window, per
genotype, onset+peak), `percell_8panel_260729.m` (8 panels/cell),
`recover_vessel_calib_260729.m`, `Ventral_surface_coherence_polar_svd_260729.m`,
`Ventral_surface_ROI_on_cartoon_260729.m`.

**`Ventral_surface_breath_time_summary_260731.m` → `breath_time_summary\` (2026-07-31)**
— population + per-cell histograms in time AND linearized phase, one file set per group,
IO standalone. Standalone: scans the archive itself, does NOT use `ensure_pool_260727`.
Run: **1506 obs → 1384 cells → 235 ACTIVE → 41 sig** (IO 140/4, Sert 48/17, Vglut2 24/13,
Vgat 16/5, ChAT 3/2, Sst 4/0). Sert's 48 active independently reproduces the count above.
Layout: row 1 = pooled `ALL | SIG | NON-SIG`, then cells 3/row left→right top→bottom,
sig block first sorted by logZ; every slot is a `[time | phase]` pair.
Decisions that are easy to get wrong on a re-run:
- **time is PEAK-centred** (t=0 = peak, sky blue); red marks the onset one median
  onset-to-peak interval earlier, and the thin red curve is the onset-triggered
  histogram placed on that same time base. Window ±2.0 IBI **of that genotype**.
- **SIG = Rayleigh `p_analytic < 0.001`, not the shuffle** — a 500-shuffle p floors at
  1/501 = 0.002 and *cannot express* α=0.001. Shuffle p is still written to the CSV.
- **units are spk/cyc %**, and for phase the order is load-bearing: occupancy-normalise
  to events/s FIRST, then rescale to %. A raw % per phase bin makes every cell look
  expiratory, because inspiration is only ~15-21% of the cycle in time.
- red = onset and sky = peak throughout, so a significant cell is marked by
  **bold title + heavy frame, NOT a red frame** (red is already taken).
- `nShuffle = 0` → no grey band (bars only); set it >0 to restore the band, the
  machinery and the correct %-scaling are still in place.
- `coh_cfg_260727.m` lives in `analysis_260727\coh_ca_breath\`, NOT `analysis_260727\` —
  the polar script's preamble omits that addpath and only works if it is already on the
  saved MATLAB path.

**SOURCES MOVED D: → E: on 2026-07-30.** `E:\260721_Sert_soma_G8s\phys\baseline` and
`E:\260728_vglut2_soma-g8s\phys`. Both `SRC` in `archive_into_summary_260729.m` and all 31
`source_path` rows of `analysis_260727\site_mapping_260729.csv` were repointed; the archive
copies themselves needed nothing, because portability rides on the recording folder NAME,
not the absolute path (see [[feedback_roi_curation_port]]).

**Vglut2/0728 breath curation is now complete (19/19 recordings) and it changes nothing in
the archive.** The binding constraint is `ca_spike_data` (6 of 19) and then ACTIVE cells: of
those 6, `MAO_×3` is excluded by policy, and `fov1_1.7x_z35` (50 ROIs) and `roi2_2.5x_z−25`
(30 ROIs) both have **0 active ROIs / 0 events** and were already excluded in the mapping csv
on 2026-07-29. Only the 3 already-archived recordings qualify. Do not re-derive this — the csv
records each exclusion and its reason. Verified by stage position too: those two sit 278 µm and
977 µm from the nearest archived site, so they were never additions to `cell1`/`cell2` anyway.

**120 µm DEPTH BAR on Vgat/0730, added 2026-08-08.** This archive is a VENTRAL SURFACE
analysis, and 0730 (Vgat-G8m shiverer) is the only session imaged deliberately deep, so
every 0730 FOV deeper than 120 µm is now in `coh_cfg_260727.cfg.excludeRecordings` (7 entries
total). **Depth is `|motorPosition(3)|`, NOT the z in the filename** — they disagree once and
it matters: `roi1_1.5x_x1050y1000_z100_3000f_26.5lp` is named z100 but its motor reads
201.9 µm, so it is excluded. Excluded six: z250 (244.7), z265 (264.7), that z100 (201.9),
z320 (314.7), z120 (**121.6 — 1.6 µm over the bar**), z160 (161.6).
**Cost, stated up front: this removes Vgat's four strongest phase-locked cells** (logZ 3.79,
3.75, 3.58, 3.44, all from z265 and z250). Vgat's best becomes 2.43 and the top of the polar
plot becomes Vglut2 3.66. Deltas: polar 419→415 cells, Vgat 16→12; map 99→93 FOVs,
2888→2598 ROIs, Vgat 722→434 cells. Only 2 of the 6 had spike data, which is why the map
loses far more ROIs than the polar loses cells.
`Ventral_surface_ROI_on_cartoon_*.m` does NOT read `excludeRecordings` on its own — the
`_260808` copy applies it explicitly. Its skip message still prints "different stage zero"
for depth-barred recordings (hardcoded string).

**Sst 0806/0807 added to the archive 2026-08-08** as `Sst\0806` and `Sst\0807`, site =
`roiN` prefix, small files only (dFF, meta, cpSAM, AVG, breath, spikes; ~0.1 GB, no TIFFs)
via `analysis_260806\archive_sst_260806_260807.m`. Both are on the NEW mount axis. **No
matcher run, so one ROI = one cell — the Sst count is an overcount.** Sst goes 4 → 45 cells
in the polar and 4 → 40 active in the per-cell summary. Once archived, do NOT also list the
`C:\` session paths as separate groups anywhere — that double-counts every Sst ROI.

Re-runs live in **sibling** folders so the reference sets survive: `polar_coh_vs_rayleigh_260808\`,
`breath_time_summary_260808\`, `breath_trig_heatmap_260806\`. The map writes to the Desktop
under `ROI_on_cartoon_260808_depthbar` for the same reason (the first run overwrote the Aug-1
`ROI_on_cartoon_stitch260801.*`).

See [[project_analysis_260727_cell_pooled_coherence]], [[feedback_roi_curation_port]],
[[project_breath_trig_heatmap_260806]].
