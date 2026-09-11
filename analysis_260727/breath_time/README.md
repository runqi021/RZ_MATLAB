# breath_time — inspiration-triggered PETH in absolute time

**The primary breathing-modulation analysis.** Replaces both the coherence score and
the phase/Rayleigh score.

## Why time, not phase

Cycle-interpolated phase (foot = 0, peak = π) stretches every breath onto a common
0–2π axis. That is robust to breathing-rate drift — the reason it was tried — but it
is **not uniformly occupied in time**. Measured on 260721_Sert: only **15% of each
cycle** falls in `[0, π)`, which is half the phase axis, and the 24-bin occupancy
max/min ratio averages **12×**.

So an unmodulated cell firing uniformly *in time* produces a large, highly
significant *raw* Rayleigh vector pointing into expiration. It needs an ECDF
correction before it can be read at all.

In absolute time from inspiration onset there is no stretch and no such bias: every
trigger contributes to every time bin, so **the null is flat with no correction**.

The cost is the mirror image — a fixed window only stays aligned while breath
duration is reproducible — and `breath_time_overlay_260727.m` measures exactly how
long that is instead of assuming it.

## Run order

```
1. breath_time_overlay_260727.m   diagnostics: window choice + trigger QC
2. breath_time_peth_260727.m      per-cell PETH, scores, shuffle test
```

Both read `cell_pool.mat` only. Config: `..\coh_ca_breath\coh_cfg_260727.m`
(`cfg.outRoot`). Or just run `..\run_analysis_260727.m`.

## The score

| column | meaning |
|---|---|
| `peak_over_baseline` | peak rate ÷ the cell's own overall rate — the modulation index |
| `latency_ms` | when that peak occurs, relative to inspiration onset |
| `trough_over_baseline`, `trough_latency_ms` | the same for suppression, so cells that go *down* are not scored as unmodulated |

**Baseline = the cell's overall rate** (total events ÷ total exposure) across the
window. That is exactly the flat-null expectation, so the ratio is the classic
modulation index. A median baseline was tried first and is wrong here: at these
event counts most bins are empty, the median is 0, and the ratio explodes — it
reached 3.7e15 on one cell before this was caught.

## Three things that are easy to get wrong, and are handled

**1. Exposure must be counted in FRAMES, not triggers.** At 30 fps a frame arrives
every 33.3 ms, so a 50 ms bin holds one frame sometimes and two others. Dividing by
trigger count leaves a regular sawtooth in the PETH that is pure aliasing. Panel 7
of the overlay figure plots the exposure ripple that gets divided out. Recordings
here span 30–47 fps, so exposure is accumulated per recording at its own rate.

**2. The window must come from the data.** No fixed window serves two datasets: the
breath cycle is **1.833 s** in 260721_Sert but **0.567 s** in
Ventral_surface_summary. A −0.5…2.0 s window is 1.1 cycles for the first and **4.4
cycles** for the second, where the pre-onset "baseline" then sits inside the
*previous* inspiration. `AUTO_WINDOW` derives it as −0.25…+1.25 × the measured
median cycle.

**3. The shuffle count must scale with the number of tests.** The smallest p a
shuffle can return is `1/(nShuffle+1)`. Benjamini–Hochberg needs `p ≤ q/nTests` at
rank 1, so 1000 shuffles against 296 tested cells puts the floor (9.99e−4) *above*
the threshold (1.7e−4) and **no cell can be called at all** — which is exactly what
happened on the first Ventral run. `nShuffle = 40 × nTested`, bounded 1000–20000.

## Significance

Circular shift of the event train, which preserves event count and autocorrelation
exactly. Implemented as an exact circular cross-correlation — the full
trigger-by-event correlogram is computed once by FFT, after which every shuffle is
an index shift into it.

- **`uniform`** (primary) — shift by a uniform random number of frames. Fully
  decorrelates from the breath.
- **`breath`** — shift by a whole number of complete breaths. For a *phase*
  statistic this is the stricter null; for a PETH in absolute time it is close to
  degenerate, since shifting by whole breaths lands events at nearly the same lag
  relative to some other onset. Reported alongside, not relied on.

**When the two disagree, do not pick one — the statistic is not behaving.** On
Ventral they returned 0 vs 19 of 296, which is a power problem, not a result.

Bins of a PETH with overlapping windows are correlated (each event falls inside
several trigger windows), which is why significance comes from this shuffle and not
from a bin-wise Poisson test.

## Trigger QC

Every detected inspiration onset becomes a trigger, and a trigger that is not really
an inspiration dilutes every PETH. Cycles whose peak-minus-foot breath amplitude
falls below `ampFrac` × that **recording's own** median are flagged.

On Sert this catches **6 of 2690 (0.2%)** and moves the pooled peak/baseline from
2.28 to 2.27 — i.e. detection is clean. (Testing against a *global* median instead
of a per-recording one wrongly reports ~10%, because it conflates between-recording
brightness with within-recording pauses.)

## Result on 260721_Sert_soma_G8s

**41 of 89 cells significant** under both nulls. Peak/baseline median **3.15**
(1.64–8.00), latency median **401 ms** (IQR 401–470 ms) — early expiration, just
after the inspiratory peak at 267 ms.

The pooled PETH rises from 0.20 to 0.46 events/s and is back to baseline by ~1 s.

Independently reproduces the retired phase analysis (preferred phase 196° ≈ 410 ms)
and the pooled coherence, neither of which shares any machinery with this.

## Measured breath facts (Sert, 2690 cycles)

| | median | CV |
|---|---|---|
| inspiration (foot→peak) | 0.267 s | 0.21 |
| expiration (peak→foot) | 1.567 s | 0.38 |
| full cycle | 1.833 s | 0.33 |

**6 × inspiration = 1.60 s = 0.87 cycles** — less than one breath, not the 2–3 it
intuitively looks like; three cycles would need 21×. And the median breath waveform
is flat after ~0.5 s, so breath 2 is unresolvable at any window length: 5% of cycles
have restarted by 0.90 s, half by 1.83 s. Long pauses exist (11.71 s, 6.70 s,
5.17 s) but are 0.6% of cycles.
