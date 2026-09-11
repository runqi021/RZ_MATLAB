---
name: Cosine-fit phase tuning — script and formula
description: Where the per-ROI P(spike | breath phase) cosine fit lives and the exact formula used, so new phase-tuning analyses stay consistent
type: reference
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
**Script:** `Breath_summary_timeNphase_260331/master_cosine_fit_phase.m`

Per-ROI cosine fit of spike probability as a function of breathing phase:

```
y = a + b * cos(phase - mu)
```

Parameter definitions (from the script header):
- `mu` = circular mean of that ROI's spike phases
- `a`  = mean of the phase-binned P(spike) (the DC level)
- `b`  = 2 × mean( P(spike) · cos(ctrs − mu) ) (modulation amplitude, closed form)
- `ctrs` = bin centers of the phase histogram

Related conventions (match across the phase-tuning / breathing scripts):
- **Phase range**: `[0, 2π]` for display; `[0, 4π]` for two-cycle plots
- **Bins**: `nBinsHist = 24`
- **Dropped frames**: `nDrop = 30` (must match `breath_combined_summary_260329.m`)
- **Significance**: `all_is_sig_spk` from `coherence_data.mat` — matched by iteration order, used for optional filter (`sigOnly`)
- **Min spikes**: 5 per ROI minimum for a valid fit

Companion scripts in the same folder:
- `breath_combined_summary_260329.m` — per-session cosine fit + phase plots (source of the `mu`, binning conventions)
- `master_coherency_summary.m` — coherence-based significance gating
- `breath_polar_phase_diagnostic_260330.m` — polar visualization
- `breath_coherence_dFF_vs_spike_260410.m` — dFF vs spike coherence comparison

**How to apply:** When writing a new phase-tuning analysis (e.g., for whisking phase instead of breathing), reuse the `a`, `b`, `mu` closed-form formula and the 24-bin phase binning. Pair with coherence-based significance from the existing coherence pipeline rather than inventing a new one.
