---
name: project-stack-dff-logz
description: "stack_dff_logz_260812.m — diagonal-wave dF/F stack of Rayleigh logZ>=2 cells, coloured by genotype from the polar palette"
metadata: 
  node_type: memory
  type: project
  originSessionId: f145e055-fef7-4364-bc1d-f95f7ce506ec
  modified: 2026-08-12T08:46:17.382Z
---

`analysis_260806\stack_dff_logz_260812.m` — sheared "diagonal wave" stack of
dF/F traces for every phase-locked cell, N chunks per cell, coloured by genotype.

Descendant of `Breath_summary_timeNphase_260331\stackDFF_chunkALL_260413.m`
(same parallelogram geometry), but the cell list now comes from the CURRENT
Rayleigh analysis instead of the old coherence `fov_map.mat`.
See [[reference-stackdff-chunkall]].

## Selection
- `rayZ >= zMin` (zMin = 2) out of
  `breath_trig_heatmap_260806\breath_trig_heatmap_GENOTYPES_peak.mat`.
  `rayZ` there is the per-ROI Rayleigh logZ, already joined onto `roiLbl`.
- At zMin=2 that is **58 cells**: Vglut2 27, Sert 19, Sst 7, IO 2, ChAT 2, Vgat 1.
  (456 of 2761 ROIs have a finite logZ at all; 133 clear logZ>=1.)

## roiLbl has no site level — and that is where IO lives
`roiLbl` is `MMDD/FOV/roi`. The archive is `Genotype/MMDD/<site>/<FOV>`, so
resolve with `dir(fullfile(root, MMDD, '*', FOV, '*_cpSAM_output.mat'))` and
read the site back off `h(1).folder`. **IO is a SITE inside the genotype
folders** (`Vglut2/1124/IO`, `ChAT/0521/IO`), not a genotype — that recovered
site is the only way to label the IO cells. All 58 resolve this way.

## Colours
Read `groups` + `group_colors` from
`polar_coh_vs_rayleigh_260808\polar_coh_vs_rayleigh_data.mat` — never a local
literal — so this figure and the polar figures cannot drift. That file already
carries all six groups including IO:
IO black, ChAT #d91a1a, Vglut2 #1aa633, Vgat #1a4dd9, Sst #8c33bf, Sert #e6731a.
Rows are ordered in genotype blocks; with six colours at once, colour alone
cannot carry identity, so block labels + legend repeat it.

## Open question — Vgat coverage
Only **14 of 521** Vgat ROIs have a finite logZ, and 1 clears 2. Also 26 Vgat
cpSAM FOVs exist on disk but only 19 appear in the dataset. Vgat sites are
cell1/2/3 + quiet (0323, 0730) — all ventral surface, no IO. If "all Vgat"
is wanted, the upstream Rayleigh/heatmap build has to be re-run over the
missing FOVs; it is not a plotting-side fix.
