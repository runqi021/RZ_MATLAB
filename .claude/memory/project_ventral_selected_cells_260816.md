---
name: project-ventral-selected-cells-260816
description: Per-cell summary figures -> popsel include/exclude GUI -> selected-cell polar + ventral map; the selection provenance chain in analysis_260806
metadata: 
  node_type: memory
  type: project
  originSessionId: f7925dbb-4fd1-49ae-800d-a7a640f2e0a7
  modified: 2026-08-17T23:37:24.295Z
---

**THE ARCHIVE IS NOW `D:\Ventral_surface_summary`** — the `C:\Users\Admin\Desktop`
copy is GONE (checked 2026-08-17). Every `_260816` script hardcodes `sumRoot = 'D:\...'`.

The chain that produces "the selected cells", in order:

1. **`analysis_260806\per_cell_summary_260812.m`** — one figure PER CELL (not per ROI),
   pooled over that cell's recordings. `selMode='active'` (pooled rate >= 2 events/min)
   or `'logz'`. Renders through `analysis_260727\coh_ca_breath\temporal_phase_cell_fig_260812.m`,
   the SAME function the GUI uses, so batch and GUI figures are identical. **The cut is
   applied to the logZ the figure itself computes**, never the archive CSV's `ray_logZ` —
   those are different quantities (archive pools with a 30 ms ca-lag; here there is NO
   ca-lag anywhere). Out: `per-cell-summary_active_260812\`, with side-bins
   `_below_active_gate_2evmin\` and `_superseded_by_merge_260814\`.
2. **`popsel_precompute_260816.m` → `popsel_gui_260816.m` → `popsel_launch_260816.m`** —
   include/exclude curation with a LIVE population readout (mean dF/F + event histogram +
   pooled Rayleigh over the included set). Keys i/e/u + arrows. One cache and one CSV per
   group in `popsel_260816\`. Two things make the click instant and exact:
   **everything is on a normalised tau/IBI axis (-1..+1 cycles, not seconds** — IBI spans
   0.43–2.83 s across the archive), and **the pooled statistic is a SUM** (`evSum(b)`,
   `occ(b)` cached per cell; each event keeps 1/occ of ITS OWN recording, so pooling =
   addition and Kish n_eff falls out of the same sums).
   State 2026-08-16: **IO 107 include / 20 exclude / 13 undecided; Sert 16 / 10 / 22.**
3. **`cell_classes_260816.m`** — hand class assignment (pre-I / I / post-I / null),
   checked against the registry. Sert is a CRITERION not a list, and the readings differ:
   `'included'` (default, popsel-included AND p<0.01) = 10; `'notexcluded'` = 17; `'all'` = 22.
   The gap is 7 significant Sert cells left UNDECIDED in the GUI (242 244 248 252 255 257 259).
   p = circular-shift permutation p on the PSTH, 1200 shuffles → floor 1/1201, so p<0.01 is
   real but p<0.001 is NOT expressible.
4. **`polar_selected_260816.m`** — dot on the weighted circular **MEDIAN** (the 260808
   archive figure moved the dot to the MEAN but kept a median-centred IQR arc — dot and arc
   described different centres; >30° disagreement for 61 of 445 cells). Recomputed from RAW
   event phases, not the 36-bin profile. Occupancy weighting mandatory (~7x uneven dwell).
   Self-checks against the archive CSV. Out: `polar_selected_260816\` (**147 cells**),
   plus `group_overlays_{onset,peak}`.
5. **`ventral_map_selected_260816.m`** — 260808 map's left/anatomy side untouched; right
   side = one filled dot per CELL (not per mask), coloured by `group_display_colors_260816`,
   over the same 147 cells. IO not drawn (107 of 147 in one clump); Sert significant only.
   Deep Vgat cells 277/278/279 DO get dots even though the 120 µm depth bar suppresses their
   anatomy patch.

**`trace_rawF_panel_260817.m`** (`analysis_260727\coh_ca_breath\`) — standalone single-panel
RAW F trace, wired into the GUI as a "Raw F trace" button; saves with `_rawF` in the stem so
it can never overwrite a summary render. Deliberately a separate file so no raw-F switch ever
lands inside `temporal_phase_cell_fig_260812`. **Raw F cannot be pooled** — counts depend on
zoom/laser/PMT gain — so ONE recording in, one panel out. `F_roi_raw`, `F_roi` and the cpSAM
`F` were verified identical (2026-08-17); alignment comes from the same file and the same
truncation as the dF/F trace, so the two panels cannot drift by a frame.

`ppt_260815_picked_cells.csv` (42 rows) = hand-picked cells for slides, with the exact
recording, ROI, trace window and trigger window each figure used.

**MANUAL F OVERRIDE, 2026-08-17 — cell 277 / ROI 32 only.** `manual_F_override_260817.m`
(analysis_260806) replaced the cpSAM fluorescence of ROI 32 in
`Vgat\0730\deep\roi1_2.4x_x1300y900_z250_...` with a hand-drawn ImageJ trace. That
recording's dFF.mat now carries a `manualF` provenance struct and the original is beside
it as `*_ch1_dFF_preManualF_260817_171221.mat`. Three things to know:
- **The csv is named `cell279_manual_F.csv` but cell 279 is a DIFFERENT recording**
  (z265/ROI29). The z250 recording's only registered cell is 277 = ROI 32. RZ confirmed
  277 is meant; the filename is a slip.
- Per-column patching is exact: dF/F is `(F - movmedian(F,win))/movmedian` with no
  cross-ROI term, so one column recomputes bit-identically. The script PROVES this by
  reproducing the stored dFF from the stored F before writing (matched to 0.0) and aborts
  if it cannot.
- **`ca_spike_data.mat` was NOT touched** — the 46 events, logZ 3.74, phase and class for
  cell 277 all still describe the CELLPOSE trace. The avg-projection panel also still
  outlines the cpSAM mask, not the hand-drawn ROI. Manual F changed the trace panels only.
Raw F mean fell 12170 -> 7047 counts, CV 0.178 -> 0.081, dF/F max 0.744 -> 0.301.

**Heatmap rows are CHRONOLOGICAL as of 2026-08-17** (`P.sortMode='none'`). They used to be
`'dt'` = sorted by cycle duration, while the y-axis said only "breath #" — which reads as
chronological and was wrong. Changed in all 7 callers of `temporal_phase_cell_fig_260812`
and the axis label now names the order. `sortMode` affects ONLY the heatmap; no statistic
depends on it. For a pooled multi-recording cell "chronological" means within each
recording, recordings stacked in registry order — there is no clock across acquisitions.

**Cell 192 removed from the Vgat `I` class** (2026-08-17, `cell_classes_260816.m`). It was
the group's only SURFACE cell (Vgat/0730 cell1, `roi1_2x_x1100y950_z13`, 81 events); the
remaining three (277/278/279) are all from the deep z250/z265 FOVs. Reason not recorded.

See [[project_ventral_summary_archive]], [[project_analysis_260727_cell_pooled_coherence]],
[[project_breath_trig_heatmap_260806]], [[reference_rayleigh_on_linear_phase]].
