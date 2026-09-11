---
name: project_whisk_breath_xcorr_epoch_merge
description: Whisk-breath pooled-xcorr masking gotcha — epoch MERGE_GAP>0 bridges quiet gaps and injects spurious long-lag peaks; use MERGE_GAP=0
metadata: 
  node_type: memory
  type: project
  originSessionId: 084a3a2a-e423-4593-9601-e540d173f351
---

In the whisk imbalance/global vs breath cross-correlation scripts
(`whisk_imbalance_vs_envelope_breath_xcorr_RZ.m` and its diag companions
`whisk_imbalance_breath_lag_outliers_RZ.m`,
`whisk_imbalance_breath_xcorr_diag_persession_RZ.m`), the xcorr is a **pooled
masked** estimate: one curve per session, summing over all whisking-masked
samples (`r = xcorr(xm,b)./sqrt(Sxx*xcorr(m,b2))`), with a **full-trace
circular-shift** null. Per-session curves are then averaged across sessions for
the summary; epochs are only a sample-selection mask, NOT an averaging unit.

**Gotcha:** `bool_to_epochs` merges epochs by extending the boundary across the
gap, so `MERGE_GAP > 0` makes the whisking mask span the quiet gap between
bouts. Those non-whisking samples are noise, and pooling them into the xcorr
produced spurious large peak lags (e.g. asym peak at -0.745 s, |r|~0.1, ns).
Setting **`MERGE_GAP = 0`** (mask follows threshold crossings exactly) removed
the weird long-lag peaks. `MIN_DUR = 1` drops tiny fragments. Pooling is robust
to bout fragmentation, so no merge is the safer default here.

**Result with clean masking:** global whisking |r| (~0.18-0.29) > asymmetry |r|
(~0.06-0.15) in every session; global lag tightly near 0; significant asymmetry
peaks all within +/-0.06 s; large off-zero asym lags were all non-significant
(flat-curve argmax artifacts). => Case 1: breath couples to total whisking
drive near-synchronously, only weakly/near-zero-lag to L-R asymmetry.

Signals: fast band `WHISK_BP=[5 50]`; global = 0.5*(env_L+env_R),
asymmetry = env(Lbp-Rbp); breath band-pass [1 20], already inhale-up
(INVERT_BREATH=false); thermal-vs-Basler exposure offset
`EXPOSURE_OFFSET_S=0.36324e-3` (0.48824-0.125 ms). See
[[project_bilateral_whisk_breath_pipeline]] and
[[project_whisk_breath_coordination_goal]].
