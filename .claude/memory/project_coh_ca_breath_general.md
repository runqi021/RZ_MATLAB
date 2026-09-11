---
name: project_coh_ca_breath_general
description: "coh_ca_breath_260722/ — general (flat single-genotype) version of the breath×Ca coherence pipeline, for any new experiment folder"
metadata: 
  node_type: memory
  type: project
  originSessionId: d0d6f257-9f5d-4554-aa44-3425c396829e
  modified: 2026-07-23T07:49:37.735Z
---

`coh_ca_breath_260722/` (repo root subfolder, created 2026-07-22) = the **general, reusable** breath×Ca-spike coherence pipeline for a SINGLE-genotype experiment whose recordings sit FLAT in one folder (e.g. `D:\260721_Sert_soma_G8s\phys\<recording>\`). Generalizes the `svd_breath_motion\Ventral_surface_*_svd_260606` set, which stays as-is for the multi-line `Ventral_surface_summary` meta-folder layout (`<root>\<Group>\<Date>\<cell|IO>\<rec>\`). See [[project_svd_breath_motion_files]].

**Why:** new experiments (Sert here; more transgenic lines coming) are one genotype per session folder, flat — NOT the nested cell-type layout the `_260606` foundation hardcodes (`scan_dirs={ChAT,Vglut2,Vgat,Sst}` → finds 0 recordings on a flat folder).

**Foundation script:** `coherence_polar_general_260722.m` (run FIRST — everything downstream reads its output).
- Config at top: `rootPath` (the flat experiment folder), `genotype` (label for all recordings, e.g. 'Sert'), `genotype_color`. Output → `<rootPath>\coherence_polar_<genotype>_260722\`.
- Discovery = flat recursive `dir(fullfile(rootPath,'**','ca_spike_data.mat'))`; every recording labeled `genotype`, except folders whose name contains 'IO' → IO group (hollow markers, group 1). `groups={'IO',genotype}`.
- Coherence math / output IDENTICAL to `_260606`: writes `coherence_polar_data.mat` (PP, labels, confC, groups, group_colors, params) + `sig_rois.csv` + `coherence_polar_all.png/.pdf`. Same struct → downstream read unchanged.
- Auto-alignment (same as original): `fps=detect_session_fps`, `nDrop=30` front-toss on breath, `T=min(numel(bw),nCa)` truncates breath to Ca length. Assumes breath cam is 2P-frame-triggered (1:1, same fps) — trims tail, does NOT resample. `ca_lag_sec=0.1` (spikes lead-shifted before coherence).
- **Verified 2026-07-22** on `D:\260721_Sert_soma_G8s\phys`: 9 recordings (those with ca_spike_data of 27), 43 ROIs, 1 sig (r=0.839). Partial because Ca spike detection only done on 9/27 + 4 dFF pending MC.

**Full folder contents (2026-07-22, renamed 2026-07-23):** `README.md` (run order) + `coh_cfg.m` (shared config — the ONLY per-experiment edit: rootPath/genotype/genotype_color; derives cohDir/cohData) + `coherence_polar_general_260722.m` (foundation) + **16 downstream scripts** copied from `svd_breath_motion\` and RENAMED to clean functional names (dropped `Ventral_surface_` prefix + `_svd_260606/07/09` suffix; line-1 header rewritten to match filename per convention): `temporal_phase_{singleROI,activeROI,sigROI,population,population_groups,window_gui}.m`, `phase_polar_allROI.m`, `phase_hist_sigROI_overlay.m`, `pop_polar_overlay.m`, `pop_group_means_spike_phase_polar.m`, `phase_popstate_traj.m`, `phase_neuron_manifold.m`, `cycletraj_pca.m`, `cycletraj_isomap.m`, `dFF_heatmaps.m`, `spike_triggered_dFF.m`. Every script's path-setup got `addpath(scriptDir)` (so `coh_cfg` resolves) and hardcoded `rootPath`/`cohData` replaced with `cfg=coh_cfg(); rootPath=cfg.rootPath; cohData=cfg.cohData;`. Analysis code otherwise byte-identical. `svd_breath_motion\` originals untouched (nested summary layout). Scripts don't call each other by name → renaming safe.
- Downstream were ALREADY flat-compatible (recursive `dir(rootPath,'**','ca_spike_data.mat')` + group labels from cohData); the copy just centralized config.
- **Two single-recording pickers** take one `folderPath` (not a scan): `temporal_phase_singleROI` → `folderPath=fullfile(cfg.rootPath, recName)` (edit recName+roi; default = sig recording ro3_4x roi 24); `temporal_phase_window_gui` → `folderPath_manual=''` + `extraROIs={}` cleared. NOTE: singleROI's header still lists DLC breath filenames (`*DLC*breath_peak_data.mat`) — verify it reads pc1 files before relying on it (pre-existing).
- **Workflow (no master script; each .m run individually):** edit `coh_cfg.m` → run `coherence_polar_general_260722.m` FIRST (auto-discovers all recordings, writes coherence_polar_data.mat) → run any downstream. New experiment = change only coh_cfg.m rootPath+genotype. See README.md in the folder.

**Consolidation (2026-07-23) — current 11-file set** (user pruned + merged + renamed; coh_cfg + coherence_polar_general_260722 = "fine, don't touch"). FINAL names (line-1 header == filename for all):
- `temporal_phase_perROI.m` — MERGED old activeROI+sigROI; top param `select='sig'|'active'|'all'` (sig=PP.r>=confC, active=nnz>minEvents, all=every ROI); each fig saved BOTH into the recording folder (`temporal_phase_svd_ROI##`) AND a per-select batch dir `temporal_phase_perROI_<select>\`.
- `temporal_phase_hist_population.m` — MERGED old population+population_groups; per-cell paginated, 4 groups; entity-1 of each group = that group's MEAN.
- `temporal_phase_hist_group_means.m` — group MEANS only (compact, one row/group); renamed from pop_group_means_spike_phase_polar.
- `temporal_phase_single.m` (static one ROI) + `temporal_phase_window_gui.m` (interactive) — per-ROI 12-panel inspectors (kept both).
- `phase_polar_allROI.m` (polar) + `phase_hist_sigROI_overlay.m` (linear) — all-ROI overlay twins (NOT merged, user declined).
- `dFF_heatmaps.m` + `spike_trigger_dFF.m` (note: "trigger" not "triggered" per user rename) — calcium.
- Naming NOT fully family-consistent (user OK with it): `hist` appears in both the population summaries and the linear overlay; overlay/calcium pairs named divergently. A prefix=family scheme was proposed and DECLINED for now.
- **4-GROUP scheme everywhere** (population + group_means): `sig / nonsig / IO_sig / IO_nonsig` = (IO? = colorIdx==ioGroup) x (sig? = r>=confC); empty groups auto-skipped (`ioGroup=-1` guard if no 'IO' in S.groups). On a flat single-genotype experiment only sig+nonsig appear.
- USER REMOVED (extra for now, originals safe in svd_breath_motion): `cycletraj_pca`, `cycletraj_isomap`, `phase_popstate_traj`, `phase_neuron_manifold`, `pop_polar_overlay`.
- KEPT distinct: `temporal_phase_singleROI` (static one ROI) + `temporal_phase_window_gui` (interactive) [user chose keep both]; `phase_polar_allROI` (polar per-ROI) + `phase_hist_sigROI_overlay` (linear per-ROI) [NOT merged — user declined]; `dFF_heatmaps`, `spike_triggered_dFF`. Names kept as-is (user declined consistency renames).

**GOTCHA — ChAT-specific group splits crash on single-genotype experiments (2026-07-23):** several downstream scripts hardcode the ChAT dataset's 3 groups `{ALL sig, IO, nonIO_nonsig}` (IO = `PP.colorIdx==ioGroup`). On a flat single-genotype experiment there are NO IO recordings → the IO group is empty → the script errors. First hit: `pop_group_means_spike_phase_polar.m` — `group_mean_dist` had `assert(nC>0,...)` at old line 132. FIXED in the general copy: empty group returns `nC=0` (no assert), the group loop `continue`s past empty groups, and the figure `tiledlayout` sizes to the number of NON-empty groups (`nRows=numel(G)`). On 260721 Sert (78 ROIs, confC=0.827, only 1 sig): plots 2 rows (ALL sig N=1, nonIO_nonsig N=77), IO skipped. **Same empty-IO-group pattern likely affects `temporal_phase_population_groups.m` and `pop_polar_overlay.m`** (both split IO/nonIO) — harden them the same way if used. Originals in `svd_breath_motion\` NOT touched (they run on the ChAT-populated summary set). Also: population-mean of N=1 is just that one cell (not meaningful) — low sig count is because Ca spike detection is only partly done.

**Inputs per recording** (produced upstream): `ca_spike_data.mat` (Ca spike detector — see [[project_ca_spike_detector_260722]]), `breath_peak_pc1.mat` + `breath_insp_start_pc1.mat` (breath_svd_pc1 → peak/trough GUIs, see [[project_breath_svd_pc1]]), `*_dFF.mat`, `*_meta.mat`. For 260721 all breath files done (27/27); ca_spike_data only 9/27.
