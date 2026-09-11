# coh_ca_breath_260722 — breath × Ca-spike coherence pipeline (general)

General, reusable version of the ventral-surface breath×calcium analysis for a
**single-genotype experiment** whose recordings sit **flat** in one folder, e.g.

```
D:\260721_Sert_soma_G8s\phys\
    ro1_2x_...\        ca_spike_data.mat, breath_peak_pc1.mat, breath_insp_start_pc1.mat, *_dFF.mat, *_meta.mat
    roi7_1.8x_...\     (same)
    ...
```

The original nested pipeline in `..\svd_breath_motion\` (layout
`<root>\<Group>\<Date>\<cell|IO>\<rec>\`) is kept as-is for the multi-line
`Ventral_surface_summary` dataset. Use THIS folder for each new experiment.

---

## The only file you edit: `coh_cfg.m`

All scripts read the experiment path from `coh_cfg.m`. Edit these three lines:

```matlab
cfg.rootPath       = 'D:\260721_Sert_soma_G8s\phys';   % the experiment's phys folder
cfg.genotype       = 'Sert';                           % the transgenic line
cfg.genotype_color = [0.90 0.45 0.10];                 % marker color
```

`cfg.cohDir` / `cfg.cohData` are derived automatically. Do NOT put paths in any
other script.

## Run order

1. Edit `coh_cfg.m`.
2. Run **`coherence_polar_general_260722.m`** (the foundation) — run this FIRST.
   It recursively scans `rootPath` for every recording that has
   `ca_spike_data.mat` + `breath_peak_pc1.mat` + `breath_insp_start_pc1.mat`,
   computes breath×spike coherence, and writes
   `…\phys\coherence_polar_<genotype>_260722\coherence_polar_data.mat`
   (+ `sig_rois.csv`, `coherence_polar_all.png/.pdf`).
3. Run **any downstream script** you want (below). Each reads `coh_cfg()` and the
   foundation's `coherence_polar_data.mat`.

There is no master script — each `.m` is run on its own (open it, press Run).
Every script `addpath`s itself + the repo root, so run order/current folder
doesn't matter as long as the foundation ran first.

## How recordings are chosen

Automatic. The foundation (and each downstream script) does
`dir(fullfile(rootPath,'**','ca_spike_data.mat'))` — it finds and processes
**all** recordings under `rootPath` that have the required inputs. You never list
files. A recording with a folder name containing `IO` is split into an `IO`
group (hollow markers); everything else is `genotype`.

## Downstream scripts

Per-ROI batch (one 12-panel figure per ROI; scans all recordings):
- `temporal_phase_perROI.m`         — set `select = 'sig' | 'active' | 'all'` at top

Population summaries (4 groups: sig / nonsig / IO_sig / IO_nonsig; empty groups skipped):
- `temporal_phase_hist_population.m`  — every cell, paginated, [time|phase|polar] per cell
- `temporal_phase_hist_group_means.m` — group MEANS only (one compact row per group)

Overlays (all ROIs on one axes, colored by group):
- `phase_polar_allROI.m`            — per-ROI spike-phase, POLAR
- `phase_hist_sigROI_overlay.m`     — per-ROI spike-phase, LINEAR, sig ROIs

Calcium:
- `dFF_heatmaps.m`                  — population dF/F + spike heatmaps (sig vs non-sig)
- `spike_trigger_dFF.m`            — spike-triggered dF/F STA (GCaMP rise time)

Single-recording inspectors (you also pick WHICH recording/ROI):
- `temporal_phase_single.m`         — set `recName` (a folder under rootPath) + `roi`
- `temporal_phase_window_gui.m`     — interactive; `selMode='manual'` needs a folder

## New experiment (e.g. a ChAT line)

Change only `coh_cfg.m`:
```matlab
cfg.rootPath = 'D:\2607xx_ChAT_soma_G8s\phys';
cfg.genotype = 'ChAT';
```
Re-run the foundation, then the downstream. Nothing else changes.

## Upstream inputs (per recording, produced before this pipeline)

- `*_dFF.mat`         — dffQC pipeline (MC → cpSAM → dFF)
- `breath_pc1.mat`, `breath_peak_pc1.mat`, `breath_insp_start_pc1.mat`
                       — `breath_svd_pc1.m` → peak/trough GUIs (DLC-free breathing)
- `ca_spike_data.mat` — Ca spike detector
- `*_meta.mat`        — for `detect_session_fps`

## Alignment (automatic)

- `fps = detect_session_fps(folder)` per recording.
- `nDrop = 30` breath frames tossed up front to match the calcium's tossed frames.
- `T = min(numel(breath), nCa)` — breath truncated to the Ca length.
- Assumes the breath camera is 2P-frame-triggered (1:1, same fps): it trims the
  tail, it does NOT resample. Valid on this rig.
