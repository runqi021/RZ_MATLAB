---
name: Whisker signal choice history — what was tried, what was kept
description: Record of signal-type options considered for DLC whisker tracking (PC1, Y only, radial, speed, complex, sum-of-PSDs) and which the user explicitly accepted or rejected
type: project
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
When analyzing DLC-tracked whisker point motion, the user has strong opinions about what signal feeds into traces and spectra. Summary of what was tried in the 2026-04-12 session:

| Option | Status | Reason |
|---|---|---|
| **Y coordinate, demeaned** | ✅ baseline | Signed, no rectification, simple |
| **X coordinate, demeaned** | ✅ kept as complement to Y | Same reasoning |
| **PC1 of (x,y)** | ❌ explicitly rejected | User said "not pc1"; data-driven axis obscures anatomical interpretation |
| **z-score** | ❌ explicitly rejected | User wants raw pixels, not normalized |
| **Radial displacement `sqrt(x²+y²)`** | ❌ rejected | Rectified → frequency doubling ("that doesnt capture negative") |
| **Speed `|dx+dy|`** | ❌ same | Rectified, same issue |
| **Fixed-axis projection (`xcosθ+ysinθ`)** | considered, not used | Requires anatomical θ |
| **Complex signal `z=x+iy`, two-sided FFT** | considered, not used | Interpretation complex for user |
| **Hilbert transform / analytic signal** | discussed theoretically | Requires narrowband input (bandpass first) |
| **X and Y as separate signed traces** | ✅ | Cleanest bipolar, nothing lost |
| **PSD = PSD(X) + PSD(Y) summed** | ✅ kept | Combined axis info without rectification |

**Current default in `whisker_L_3w_batch.m`**: per point, compute demeaned X and demeaned Y separately. PSD computed on each axis (Chronux multitaper), then summed for plotting. Trajectory figure plots X vs Y directly (trace in 2D).

Bandpass is OPTIONAL and OFF by default (`bpBand = []`) — user preferred unfiltered signals to not mask any frequency content.

**Why:** User explicitly rejected any transformation that loses directional/signed info. Rectified signals (radial, speed, envelope, rectified traces) cause frequency doubling which ruins the spectrum. The "sum of per-axis PSDs" captures total motion power without rectification.

**How to apply:** When adding new whisker analyses, default to demeaned Y + demeaned X (no z-score, no PC1). For spectrum, compute per-axis and sum if a single curve is needed. Never rectify before FFT.
