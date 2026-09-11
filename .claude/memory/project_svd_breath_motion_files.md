---
name: project_svd_breath_motion_files
description: Current inventory of the svd_breath_motion/ scripts (13) after the 2026-06-06 rename/de-dup
metadata: 
  node_type: memory
  type: project
  originSessionId: 2a7a33fd-578b-4aff-90c1-792d0024ee79
---

`svd_breath_motion/` was renamed/de-duplicated (2026-06-06) so every script's filename matches its
line-1 header and the confusing "phase/pca/traj/polar" names are self-explanatory. **13 scripts** (all
`Ventral_surface_*_svd_260606.m` except the STA one):

- **coherence_polar** — THE coherence analysis; writes `coherence_polar_data.mat` (sig=`PP.r>=confC`)
  that every other script reads. Source/foundation.
- **temporal_phase_singleROI** — 12-panel temporal-phase fig for ONE chosen ROI (was `*_SingleView`).
- **temporal_phase_activeROI** — batch 12-panel per ACTIVE ROI (`nnz>5`).
- **temporal_phase_sigROI** — batch 12-panel per SIGNIFICANT ROI.
- **temporal_phase_population** — population grid, all sig ROIs, [time|phase|polar] triples.
- **temporal_phase_population_groups** — same grid for IO / nonIO-nonsig sets.
- **phase_polar_allROI** — all-ROI spike-phase POLAR overlay (one line/ROI, by group).
- **pop_polar_overlay** — ONE polar axes overlaying the POPULATION-MEAN spike-phase histogram for 3
  sets: sig (n=10) / IO (n=272) / nonIO_nonsig (n=15); same method as population_groups; `normEach`
  toggle for raw spk/cyc % vs shape-only. (added 2026-06-06)
- **phase_popstate_traj** — population STATE-SPACE trajectory: `[24 bins × N neurons]`, embed bins via
  PCA+Isomap (2×3 fig). NEURON-space (was `phase_pca`). [[project_ventral_temporal_phase]]
- **phase_neuron_manifold** — neuron MANIFOLD PCA, each neuron a dot by pref phase (was `phase_pca_combined`).
- **cycletraj_pca** / **cycletraj_isomap** — TEMPORAL-TRACE per-cycle trajectory (per sig ROI's 50
  best-breath-corr cycles of spike-prob + pooled population loop). Were `traj_pca`/`traj_isomap`.
  [[project_traj_temporal_trace]]
- **dFF_heatmaps** — 6 population dF/F+spike figures (time/phase heatmaps, hists, single-cycle).
- **spike_triggered_dFF** — spike-triggered dF/F STA + per-ROI & individual-spike heatmaps.

DELETED: `temporal_phase_polarHistFit_allROI_svd_260606` (byte-identical to phase_polar_allROI except a
crossed header line). Earlier-deleted: `coherence_polar_phase_svd`, plain `temporal_phase_svd_260606`,
`spike_phase_polar_svd`. `ca_lag=0.1 s` across all (see [[project_ventral_temporal_phase]]).
