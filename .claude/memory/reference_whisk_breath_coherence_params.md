---
name: reference-whisk-breath-coherence-params
description: "Chronux coherence parameters from Moore 2013 / Deschênes 2016 for whisk-breath and L/R-whisker coherence, segmented by breathing state"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 46e13b9d-a1c8-421e-b6fb-05c853569c3b
---

Spectral coherence params for replicating Moore 2013 Nature / Deschênes 2016 Neuron (Chronux `coherencyc`/`coherencysegc`), phase 0 = protraction onset (Moore) or peak (Deschênes):

**Moore 2013 (awake behaving):** segment data by behavior BEFORE coherence, average over non-overlapping segments, **TW = 1**:
- basal respiration → **1 s** segments (peak ~2 Hz)
- sniffing & inspiratory whisks → **500 ms** segments (peak ~6 Hz)
- intervening whisks → **300 ms** segments (peak ~8 Hz)
- inspiratory whisk = onset within **100 ms** of an inspiration; intervening = not within 100 ms.
- Event cross-corr significance: KS test vs uniform; "modulation depth" = the KS statistic (0.08 basal, 0.26 sniff).

**Deschênes 2016:** Chronux coherence over **2 s** segments, **TW = 2** (K=3 tapers, W≈1 Hz). Bilateral L/R coherence: sniff-locked whisks high coherence + ~0 L/R phase; intervening whisks (only obtainable after a sigh, expiratory pause ≥500 ms) lower coherence + phase drift. 95% CI shaded bands. basal <4 Hz / sniff >4 Hz.

For OUR goal: make TWO coherence figures (sniffing epochs vs basal epochs) for L↔R whisker (and/or whisk↔breath). Segment each breathing-state epoch into fixed windows (1 s basal / 500 ms sniff per Moore, or uniform 2 s/TW=2 per Deschênes), stack as Chronux trials, report C(f) + phase(f). Our existing population coherence script does the trial-averaged `coherencyc` already — split its input by per-cycle breathing rate. NOTE intervening-whisk coherence not feasible in our data ([[project-whisk-breath-coordination-goal]]). Linked: [[reference-phase-resetting-method]], [[project-chronux-multitaper-conventions]].
