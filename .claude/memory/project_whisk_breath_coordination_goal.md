---
name: project-whisk-breath-coordination-goal
description: Overall aim of the 260615 thermal+Basler whisking/breathing analysis and why intervening-whisk coherence is infeasible
metadata: 
  node_type: memory
  type: project
  originSessionId: 46e13b9d-a1c8-421e-b6fb-05c853569c3b
---

Dataset 260615_thermalNbasler (mouse, animals 5916296/97/300): bilateral Basler whisker video (~400 fps, DLC L/R angle) + thermal nostril breathing (~400 fps). Goal = replicate the Moore 2013 Nature / Deschênes 2016 Neuron "breathing is the master clock for whisking" framework on our own data.

Three target analyses:
1. **Left/Right whisking synchrony** — L/R cross-correlation (canonical: [[reference-coherence-TW-design]] style + `whisk_LR_xcorr_all_vs_intervening_RZ.m`); shows L and R whisk in phase (peak at 0 ms).
2. **Phase resetting of whisking by breathing** — replicate Moore SI Fig S4/S5, see [[reference-phase-resetting-method]].
3. **Whisk↔breath coherence split by breathing state** (sniffing vs basal) — two figures; see [[reference-whisk-breath-coherence-params]]. We already have a population coherence+phase-vs-frequency script.

**Key constraint (confirmed against literature):** our recordings have NO long apneas/sighs. Intervening whisks (whisks >~100 ms from any inspiration during basal breathing) are SINGLE or few events, never an epoch. Deschênes 2016 got intervening-whisk L/R coherence ONLY by selecting whisks after a sigh when the expiratory pause lasted ≥500 ms — we don't have that. So **per-event coherence on intervening whisks is impossible** here (a 1 s window pollutes with neighbouring breath-locked whisks; <200 ms gives no spectral resolution). For intervening whisks use event-based / sample-mask measures (L/R xcorr over an intervening sample-mask) NOT spectral coherence. **Why:** breathing is ~6 Hz here so a whisk is at most ~180 ms from a breath. **How to apply:** reserve spectral coherence for sniffing-vs-basal epoch contrasts; characterize intervening whisks with the per-sample L/R xcorr already built.
