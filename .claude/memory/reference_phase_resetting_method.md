---
name: reference-phase-resetting-method
description: Moore 2013 Nature SI Fig S4/S5 phase-resetting-of-whisking-by-breathing definitions and how to build the PRC + period-shift plots
metadata: 
  node_type: memory
  type: reference
  originSessionId: 46e13b9d-a1c8-421e-b6fb-05c853569c3b
---

Phase resetting of whisking by breathing — Moore, Deschênes, Kleinfeld 2013 Nature 497:205, SI Fig S4/S5. Inspiration = external perturbation to the whisk oscillator. Per perturbed whisk:
- **t_reset** = latency from protraction ONSET to the inspiration inside that whisk cycle.
- **T0** = expected (unperturbed) whisk period = duration of the PRECEDING whisk (whisking freq stable within a bout).
- **T1** = duration of the perturbed whisk (the cycle containing the inspiration).
- **φ_reset = 2π·(t_reset / T0)** — breath phase within the whisk cycle.
- **Δφ_whisk = 2π·(T1 − T0) / T1** — resultant phase change (>0 next whisk delayed, <0 advanced).

Plots: **Fig S4b** = PDF of φ_reset during BASAL only (flat → inspiration hits all whisk phases). **Fig S4c** = PRC: Δφ_whisk vs φ_reset, basal only; result = inspirations late in the cycle (retraction phase) ADVANCE the next protraction. **Fig S5** = raw period shift (T1−T0, ms) vs t_reset (ms), basal+sniff pooled, done BOTH directions: (a) whisk-reset-by-breath → intervening whisks shortened by a breath; (b) breath-reset-by-whisk → no bias. Asymmetry ⇒ unidirectional breath→whisk coupling.

Apply a **~30 ms correction** (inspiratory drive reaches face ~30 ms before measured nasal inspiration) to inspiration onset times before S4c, per Fig 1d.

Phase/detection conventions (Moore): whisker BP **3–25 Hz** 3-pole Butterworth filtfilt → Hilbert; whisk cycle = Hilbert phase reset; accept if peak-to-peak >5° AND duration <250 ms; **protraction onset = crossing of 10% of that cycle's min-to-max amplitude**, phase 0 = protraction onset. Breath BP **1–15 Hz** same method, inspiration onset = 10% crossing. Rat thresholds basal <3 Hz / sniff >5 Hz; MOUSE sniff >6 Hz / basal <5 Hz. (Deschênes 2016 uses phase 0 = protraction PEAK and basal<4/sniff>4 — pick one convention and state it.) Linked: [[project-whisk-breath-coordination-goal]], [[reference-whisk-breath-coherence-params]].
