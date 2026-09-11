---
name: reference-rayleigh-on-linear-phase
description: "Rayleigh test on cycle-interpolated breath phase: occupancy is ~7x uneven so it MUST be corrected, and the weighted Rayleigh needs Kish n_eff not raw n"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-07-30T06:53:00.782Z
---

`analysis_260727\Ventral_surface_polar_coh_vs_rayleigh_260729.m` — two polar panels,
same ROIs, same phase, scored two ways. Panel A coherence (r = band-mean |C|, dashed
confC). Panel B Rayleigh (r = **log Z**, dashed log(−ln α)).

**The phase is piecewise-LINEAR between landmarks** — `piecewise_phase_local`: feet at
0/2π, peaks at π, linear ramps in TIME between consecutive landmarks, strict alternation
enforced, NaN outside the first/last landmark. Interpolating peak↔trough warps every
cycle onto its own [0, 2π] regardless of duration, which is what makes phase usable when
breathing is irregular. This is the phase-domain half; `breath_time_peth_260727.m` is the
absolute-time half.

**Occupancy is badly uneven and MUST be corrected.** Measured on Vglut2/1124: 27–47
frames per 10° bin across the fast inspiratory half vs **281–304 across the slow
expiratory half — ~7×**. A cell firing uniformly in time would look strongly
expiration-locked. Correction = weight each event by 1/(frames in its bin), i.e. take
the resultant of the events/frames rate profile. Keep the ANGLE in raw linear-phase
units so 0 = onset and π = peak for every recording and angles stay poolable; remapping
phase by its own ECDF also debiases but then the peak lands at a different angle in
every recording.

**THE TRAP: a weighted Rbar with raw n breaks the analytic Rayleigh Z.** `Z = n·Rbar²`
assumes n independent, EQUALLY weighted angles. With ~7× weights a couple of events in a
sparse inspiratory bin get huge leverage — this produced **p = 5e-79 where the
circular-shift shuffle said p = 0.05**. Fix is Kish's effective sample size:

    nEff = (Σw)² / Σw²        Z = nEff·Rbar²

which reduces to n under uniform occupancy. Here median `nEff/n = 0.66`. After the fix,
analytic vs shuffle log10 p agree to **−0.13 decades, Spearman ρ = 0.957**. Analytic
stays mildly liberal (83 vs 52 at α=0.05), so the dashed circle is a guide and
`ray_p_shuffle` in `polar_scores.csv` is the authority. Also: `spike_train` is binary
with single-frame events (run length always 1), so `nnz` is a true event count — event
inflation was NOT the cause.

**Panel B marks are boxplot-style, not the resultant vector**: dot = occupancy-WEIGHTED
circular median of the event phases, arc = weighted circular IQR — i.e. a direct summary
of the per-ROI phase histogram. Median via grid search (no closed form on the circle);
quartiles from `wprctile` of `angle(exp(1i*(a-med)))`. Both use the CONTINUOUS event
phase, not the bin centre. Arcs are drawn **only for shuffle-confirmed ROIs** — arcs for
all 368, or even the 78 above the liberal analytic circle, are mostly near-uniform cells
whose IQR approaches the full circle and they pile into rings of ink. Median IQR width:
Vglut2 61°, ChAT 50°, IO 68°, Sert 75°. The weighted vs unweighted median differs by only
~18° in practice, so the correction is not distorting the angle.

**Why log Z works as a radius**: p ≈ exp(−Z) so Z_crit = −ln α is n-INDEPENDENT — one
circle valid for every ROI, unlike Rbar whose critical value depends on n. α=0.05 →
logZ = 1.10; α=0.001 → 1.93. Z ≤ nEff, so an ROI needs ≥4 events to reach α=0.05 at
all — self-limiting, no extra event-count gate needed. Clamp r at 0: `polarplot` with a
negative radius silently reflects to the opposite angle, i.e. the WRONG phase.

**POOLED PER CELL, not per (recording, ROI).** A cell re-imaged at another Z appears once
per recording; scoring each separately double-counts it and makes significance
recording-dependent — Sert cell 57 gave logZ [0.98 2.24 2.83], p_shuffle [0.166 0.006
0.020]. Phase is the coordinate where pooling is legitimate (0 and π mean the same
landmark in every recording; absolute time does not), so event phases are concatenated
across recordings, each keeping **its own recording's occupancy weights**, and n_eff
accumulates. The null shifts each recording's train independently and re-pools. Identity
comes from the same `cell_link_sources` registry as the cartoon map, and `tossed` rows are
dropped. 367 observations → **358 cells**; only Sert actually merges (67 → 58, 7 cells
span >1 recording), since the other sessions have no matcher run. Pooling can *create*
significance from power: cell 159 was n.s. in all three recordings (p 0.08/0.22/0.07) but
p = 0.022 pooled. Two CSVs: `polar_scores_percell.csv` (what is plotted) and
`polar_scores_perobs.csv` (unpooled coherence, for reference).

**Panel layout**: the coherence panel is HIDDEN by default (`showCoherencePanel = false`)
— log Z on the linear phase is the measure. Coherence is still computed so its columns
survive in the CSV. Dashed circles at α = 0.05 (logZ 1.10) and α = 0.001 (logZ 1.93),
labelled on the figure with a white background box because they cross the data cloud.
Note the panel-A arc (jackknife phase CI, ~21° median) and the panel-B arc (IQR, ~61°) are
NOT comparable — a confidence interval versus a dispersion, differing ~3×.

**Results (archive, 368 ROIs, minSpikes=2).** Rayleigh is far more sensitive than
coherence — 65 Rayleigh-only vs **1** coherence-only, 13 both. They agree on phase:
circular-mean offset −0.33 rad (sd 0.56) over the 13. Shuffle-confirmed circular-mean
phase per group:

| group | n | phase | concentration |
|---|---|---|---|
| Vglut2 | 8 | **+1.66 rad (95°)** mid-inspiration | 0.96 |
| Sert | 26 | **−3.10 rad (178°)** at the inspiratory peak | 0.72 |
| ChAT | 2 | +3.10 rad (177°) | 0.74 |
| IO | 16 | −1.88 rad (−108°) expiratory | 0.56 |

So Sert lags Vglut2 pFN by ~85° of the breath cycle.

See [[project_analysis_260727_cell_pooled_coherence]], [[project_ventral_summary_archive]].
