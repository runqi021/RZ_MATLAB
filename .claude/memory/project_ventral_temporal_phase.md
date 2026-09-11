---
name: project_ventral_temporal_phase
description: "Ventral surface temporal-phase analysis — single-file + batch over significant ROIs, significance criteria"
metadata: 
  node_type: memory
  type: project
  originSessionId: bf702b06-f693-419e-9d92-9f3239b43cc4
---

Temporal-domain companion to the Ventral_surface coherence-polar (phase) analysis. Uses the EXACT piecewise breath phase (onset=0, peak=π; `piecewise_phase_local`) — NOT a cosine fit, NOT Hilbert — to bin Ca-spikes.

## Scripts
- `Ventral_surface_temporal_phase_260601.m` — interactive single-file (set folderPath+roi). 7×4 tiled figure: dF/F+breath trace (rows1-2), phase sawtooth (row3), PSD/coherence/peak-trig avg+heatmap, peri-peak hist, avg-proj crop, event-phase polar, coherence polar, linear 0→4π, spike-triggered avg. Tile map block at top of FIGURE section; edit one tile number to move a panel. User edits this file directly — do NOT modify it for batch needs.
- `Ventral_surface_temporal_phase_batch_260602.m` — batch over significant ROIs; standalone PORT of the single-file figure as a local function (so the interactive file is untouched). One figure per ROI, saved in each ROI folder + collected in `Ventral_surface_summary/temporal_phase_batch_260602/`.

## Significance source
Significant ROI = `PP.r >= confC` from `Ventral_surface_summary/coherence_polar_260528/coherence_polar_data.mat` (confC=0.8269, TW=4, alpha=0.001). 11 significant ROIs (1 ChAT/0521, 1 ChAT/0522, 5 Vglut2/0224, 4 Vglut2/1124). Batch parses `labels{i}='group/date/recName/roi'` and resolves the folder by recursive search for recName (handles the `cell#` level).

## Coherence dot must replicate the polar analysis EXACTLY (else borderline ROIs disagree with the threshold)
The coherence-polar panel recomputes coherence; to match the stored PP.r it MUST use: TW=4 PSD for the FWHM band, band clamped to `f_breath_search=[0.2 4]`, the **Vglut2/1124 +1-frame trigger fix** (shift events +1, delay breath 1 frame), and `Tc=min(numel(bw),numel(spike_train))`. Using the quickview TW=6 band or skipping the 1124 fix shifts r at the 3rd decimal — enough to drop ChAT/0522 (r=0.828, confC=0.827) under the line. After the replica, all 11 recomputed r match PP.r to 3 decimals. See [[project_dffQC_pipeline]].

## `_svd_260606` set: onset-aligned histograms + Vglut2/1124 shift now hard-coded EVERYWHERE
The current iteration is the `svd_breath_motion/*_svd_260606.m` family (source = `coherence_polar_svd_260606/coherence_polar_data.mat`). Two cross-script conventions were unified (2026-06-06):
- **Histogram style**: time hist is ONSET-aligned (0 = insp onset, spikes spk/cycle-% on left axis, breath-PEAK distribution sky on right axis, single cycle tiled to 1.5 period via `time_from_prev`); linear phase hist is single-cycle `0..2π` (no `0..4π` duplicate). Applied to SingleView, batch, population, groups. (Plain `temporal_phase_svd_260606.m` + `spike_phase_polar_svd` still carry the OLD peak-aligned / `0..4π` panels.)
- **Vglut2/1124 rising-edge +1 shift**: previously only in the coherence source + batch's internal coherence recompute. Now **hard-coded in the breath-event load path of ALL 10 consumers** (SingleView, plain single, spike_phase_polar, batch main path, population, groups, polarHistFit_allROI, phase_polar_allROI, phase_pca, phase_pca_combined) so the HISTOGRAMS shift too, not just the coherence dot. Mechanism: delay breath by 1 frame — `ev=[0;ev(1:end-1)]; ef=[0;ef(1:end-1)]` (+ `bw=[bw(1);bw(1:end-1)]` where breath feeds coherence). Detection: single-folder/loader scripts use `contains(folderPath, fullfile('Vglut2','1124'))`; batch uses `strcmpi(group,'Vglut2') && strcmp(recDate,'1124')`.
- **GCaMP lead-comp: 0.015 → 0 → 0.2333 → now 0.1 s (user-set 2026-06-06)**: final value `ca_lag = 0.1 s`
  = **3 frames @30 Hz**, set in ALL svd_breath_motion scripts that define it (16 incl. coherence source,
  SingleView, plain single, batch→activeROI/sigROI, population, groups, allROI×2, phase_pca×2,
  spike_phase_polar, traj_pca/isomap, dFF_heatmaps; STA script has unused `ca_lag=0.10`). The 0.2333
  (7-frame, data-driven STA peak→onset) value below is the PRIOR setting — user reverted to 0.1. Source of 0.2333: the spike-triggered dF/F STA (`Ventral_surface_spike_triggered_dFF_260606.m`) showed the `ca_spike_data` events sit ON the dF/F PEAK (time-to-peak=0 for 100% of ROIs), so they lag the true spike by the rise (~233 ms); shifting spikes EARLIER by 7 frames re-anchors event→rise-onset. Applied in TWO ways: (a) histogram/population/traj spike shift `find(spk>0)-round(ca_lag·fps)` (=3 frames @0.1s); (b) **coherence now lead-SHIFTS the spike train before `coherencyc`** — `stL=[st(1+lag:end);zeros(lag,1)]` — NOT a post-hoc θ rotation anymore (changed 2026-06-06 per user "shift spike train then go into coherence"). Done in the source `coherence_polar_svd` AND the 3 redraws (`make_temporal_phase_fig` in activeROI+sigROI, and singleROI); the old `wrapToPi(th-2π·f_pk·ca_lag)` lines were deleted. Because a time-shift leaves coherence MAGNITUDE unchanged, **PP.r/confC/the 10 sig ROIs are identical** — only `PP.th` moves, and now for ALL 297 ROIs (the old rotation was sig-only). `coherence_polar_data.mat` was regenerated (confC=0.827, 10 sig). dF/F-TRIGGERED views (STA `spike_triggered_dFF`, the spike-triggered-avg + raster panels in single/active/sig) stay on RAW spikes by design. fps-safe via `round(0.1·fps)`; all data 30 Hz so =3.
- **`Ventral_surface_coherence_polar_phase_svd_260606.m` DELETED (2026-06-06)** — the sibling that excluded Vglut2/1124/IO; user judged it pointless. Nothing in svd_breath_motion loaded its output. (The older root `Ventral_surface_coherence_polar_phase_260529.m` + its consumer `cycle_dFF_insp_onset_260529.m` are separate and untouched.)
- The batch script was SPLIT + renamed (2026-06-06): `Ventral_surface_temporal_phase_batch_svd_260606.m`
  **DELETED**, replaced by two scripts that share the same `make_temporal_phase_fig` 12-panel figure:
  **`Ventral_surface_temporal_phase_activeROI_svd_260606.m`** (active = `nnz(spike_train>0) > minEvents`,
  minEvents=5; iterates all `S.labels`; out `temporal_phase_activeROI_svd_260606/active##_...`) and
  **`Ventral_surface_temporal_phase_sigROI_svd_260606.m`** (significant only, `find(PP.r>=confC)`,
  the 10 sig ROIs; out `temporal_phase_sigROI_svd_260606/sig##_...`; drops the `count_events` gate).
- Why the shift: breath cam is 2P-frame-triggered (1:1); 1124 used the rising edge so breath leads calcium by 1 frame. the `ca_lag_sec` GCaMP lead-comp is now 0 (see below). The standard `nDrop=30` breath trim exists to match dFF's already-removed `TossFrames=30`.

## Population heatmaps — SVD version
`svd_breath_motion/Ventral_surface_dFF_heatmaps_svd_260606.m` (2026-06-06) = current-method port of the
old root `Ventral_surface_dFF_heatmaps_260528.m`. **6 figures** (SIG | NON-SIG, out
`dFF_heatmaps_svd_260606/`): ROI-mean dF/F heatmap (time + phase), pooled spike-count hist (time),
spike-phase hist+cosine, duplicating PSTH (time), per-ROI nearest-onset spike-count heatmap. User
DROPPED the old spike-rate heatmaps + single-cycle dF/F stacks. Conventions: source
`coherence_polar_svd_260606` (SIG=`PP.r>=confC`, n=10/297, sorted by `PP.th`); breath PC1 files; phase
onset=0/peak=π via `piecewise_phase_local` (per-cycle resample on cumulative phi, NOT uniform);
time figs ONSET-triggered; spikes lead-shifted 7 fr (`ca_lag=0.2333`); Vglut2/1124 +1 fix; phase panels
single-cycle 0..2π (no 4π tile). folderMap by recName from `**/ca_spike_data.mat`.

## Phase histograms are now OCCUPANCY-COMPENSATED (2026-06-07)
All phase spike histograms plot `spk/cyc %` reweighted by phase dwell: `pctC = (100·N_b/nCyc)·(Ō/O_b)`
where `O_b = histcounts(mod(phi,2π), edgesP)` = frames the breath spends in phase bin b (WRAPPED phase —
all cycles), `Ō = mean(O_b over occupied bins)`. Keeps spk/cyc % UNITS (NOT a Hz rate — user explicitly
rejected the Hz version). Rationale: the piecewise onset=0/peak=π phase maps insp→[0,π] and exp→[π,2π]
regardless of duration, so equal-phase bins have unequal dwell; raw `N_b/nCyc` over-weights the
longer-dwell half (usually expiration). A uniform-firing neuron → flat after compensation. Applied to
ALL 10 phase scripts: pop_polar_overlay (2-panel raw vs compensated), population, population_groups,
phase_polar_allROI, temporal_phase_singleROI, activeROI+sigROI (make_temporal_phase_fig), neuron_manifold
+ popstate_traj (dim-reduction inputs), dFF_heatmaps Fig4 (pooled: sum per-ROI `O_b`, reweight, cosine
fit from reweighted dist). NOT applied to coherence (freq-domain) or cycletraj (continuous spike-prob).
TIME histograms NOT compensated (equal-duration bins already). Known leftover: `singleROI`'s `occ`/`rate`/
`muT` arrow still uses UNWRAPPED `phi_pw` (first-cycle only) — separate from the corrected histogram `occW`.

## Gotcha: Vglut2/0224 breathing files have NO `DLC` prefix
The 0224 Vglut2 recordings use DLC-free (MC-shift) breathing → file is `breath_peak_data.mat` (not `*DLC*...`). Globs must fall back from `*DLC*breath_peak_data.mat` to `*breath_peak_data.mat` (and same for insp_start). [[project_breathing_is_motion]]
