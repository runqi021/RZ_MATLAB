# phase_rayleigh — breathing modulation by Rayleigh statistic, not coherence

Replaces coherence as the **primary** breathing-modulation score, after
Sirota / Karalis, adapted where their assumptions don't transfer to this preparation.

Coherence asks whether respiration and calcium share power at a stable *frequency*.
It weakens whenever breathing rate drifts — even for a cell that fires at exactly the
same point of every breath. This asks that second question directly.

Runs on `cell_pool.mat`, so phase, events, breath boundaries and cell identity are the
same objects the coherence pipeline used. Nothing is re-derived and the two cannot drift.

## Run order

```
1. phase_rayleigh_260727.m           statistics + nulls + occupancy QC figure
2. phase_rayleigh_polar_260727.m     the polar replot (radius = logZ)
3. phase_rayleigh_percell_260727.m   one QC panel per cell
```

Output: `<phys>\analysis_260727\phase_rayleigh\`. Config is `..\coh_ca_breath\coh_cfg_260727.m`
(`cfg.phaseDir`) — no separate config file.

## The occupancy problem, measured

Phase is cycle-interpolated, not Hilbert: foot = 0, peak = π, next foot = 2π. That is
robust to breathing-rate change, which is why it's used — but inspiration and expiration
have very unequal durations, so phase is **not uniformly occupied in time**.

On this dataset: only **15% of each cycle** falls in `[0, π)`, which is half the phase
axis. The 24-bin occupancy max/min ratio averages **12×**, reaching 22×.

An entirely unmodulated cell firing uniformly *in time* therefore produces a large and
highly significant *raw* Rayleigh vector pointing into expiration. Uncorrected phase
statistics on this data measure breathing asymmetry, not neural tuning. The ECDF
correction `ψ = 2π·F(φ)` is load-bearing, not a refinement — `phase_occupancy_qc.png`
is the evidence, and panel 1 of every per-cell figure shows it for that cell
individually.

`F` is built **per recording** (time in `[0, π)` ranges 0.10–0.24 across the 26
recordings) and applied before pooling across a cell's recordings.

## Two preferred phases, and why

Sirota reports preferred phase from **raw** phases, since the ECDF distorts the
biological axis. That holds at mild asymmetry; at 12× it does not — the raw resultant is
dragged toward expiration for *every* cell, modulated or not. Measured here: the two
differ by a **median of 40°**.

| column | meaning |
|---|---|
| `pref_phase_rate_deg` | **primary.** Circular mean of the occupancy-normalised rate profile (events per second per bin). A true rate maximum on an undistorted axis. |
| `pref_phase_raw_deg` | Sirota's convention, kept for comparison. |

## Significance

Sirota excluded point processes with < 200 events because Z is sample-size biased. Only
**15 cells** here reach 200, so that threshold cannot be copied and the asymptotic
Rayleigh p-value cannot be trusted. Significance comes **entirely** from a
cycle-preserving shuffle, valid at any n because the null is built from the cell's own
event train. Two nulls, because they fail in opposite directions:

- **`breath`** (primary) — circular shift by a whole number of complete breaths.
  Preserves autocorrelation *and* alignment to cycle boundaries. **Caveat:** with very
  regular breathing a whole-breath shift barely changes phase, making it conservative.
  `null_phase_shift_rad` reports how far it actually moved phases — it measured
  **1.31 rad** here (π/2 = 1.57 would be full decorrelation), so it is not degenerate on
  this data. If that number is ever small, trust the uniform null instead.
- **`uniform`** — uniform random circular shift. Fully decorrelates phase; cannot be
  degenerate.

1000 iterations each, Benjamini–Hochberg FDR at q < 0.05 across tested cells. The null
distributions themselves are stored, and plotted per cell, rather than summarised.

**Inclusion:** pooled event rate ≥ 1/min → 79 cells. Every cell with ≥1 event is still
written to the CSV with its counts and `tested = false`.

## Modulation depth — use `mod_depth_cos`, not `mod_depth`

`mod_depth = (max−min)/(max+min)` of the 24-bin rate profile **saturates at exactly 1
as soon as any bin is empty**, which at a median of 70 events over 24 bins is nearly
every cell. It was uninformative on the first run — every point piled on 1.00. It is
kept in the CSV only for continuity.

The usable effect size is `mod_depth_cos`: a cosine fit to the occupancy-normalised
rate profile, `rate(θ) ≈ a + b·cos(θ − μ)`, with **depth = b/a**. Robust to empty bins,
and the same `a + b·cos(phase − μ)` convention already used elsewhere in this project.
It spans 0.4–1.4 here with R² 0.2–0.7, and separates significant from non-significant
cells cleanly.

The fit is **weighted by occupancy**: for a Poisson rate `N/T`, variance scales as
`λ/T`, so bins holding less time are less precise. Inspiration bins hold ~1/7 of the
time and are weighted down accordingly. This is a precision weighting, not a preference
for expiration — the fitted `μ` remains free to point anywhere, and
`pref_phase_cos_deg` is reported so it can be checked against `pref_phase_rate_deg`.

## Result on 260721_Sert_soma_G8s

| method | significant |
|---|---|
| per-ROI coherence | 4 / 247 |
| pooled per-cell coherence | 14 / 99 |
| **Rayleigh, occupancy-corrected, shuffle-tested** | **40 / 79** |

Both nulls returned 40. Median logZ = 1.18, max 4.16. Top cells prefer ~170–250°, i.e.
just after the inspiratory peak — early expiration.

The gain over coherence is the expected consequence of asking a cycle-locked rather than
a frequency-locked question in a preparation where breathing rate drifts.

## Reading the polar

Angle = preferred phase (rate-normalised). Radius = logZ. Fill = significant after
shuffle + FDR. **logZ goes negative** when n·r̄² < 1, which is a genuine "no modulation"
result; a polar radius cannot be negative, so those are clamped to the origin and marked
with a grey ×, with true values in the CSV. The radial origin is logZ = 0, i.e. Z = 1,
the null expectation — so distance from the centre reads as evidence.

## Kept, not replaced

Coherence stays as a secondary measure — the two answer different questions:

| measure | captures |
|---|---|
| Rayleigh phase modulation | activity at a consistent point of each breath, even when rate varies |
| coherence | stable shared oscillatory power and phase at a frequency |
| modulation depth | how much phase-binned rate changes across the cycle |
| cycle reliability | fraction of breaths the cell is recruited near its preferred phase |

Cycle reliability is reported with **both denominators**: `recruitment` (÷ all breaths)
and `precision` (÷ event-bearing breaths). Low recruitment with high precision = fires
rarely, but always at the same point of the breath.
