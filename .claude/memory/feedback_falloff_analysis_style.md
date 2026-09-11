---
name: feedback_falloff_analysis_style
description: "Runqi's corrections during the falloff work: no saturation exclusion, no stitching, no unrequested additions, verify metadata before iterating"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: f3c40add-59ae-4001-987d-927c5eb434ad
  modified: 2026-08-06T00:31:25.197Z
---

Corrections given repeatedly during the optical-penetration work ([[project_falloff_optical_penetration]]), 2026-07-29 → 2026-08-05. These generalise beyond that project.

**Do not exclude or fade data for "saturation".** I built saturation detection, faded those points, and dropped them from fits. Runqi: *"i dont care abut the saturation, what is your point?"* and *"stop fading and excluding what u defined as saturation"*. It also created invisible gaps between panels (some panels showed clipped slices, others silently dropped them). **Why:** it is his data and his call which points are valid; silently removing points is worse than showing them. **How to apply:** plot everything; if a subset is genuinely unusable, say so in text, don't delete it from the figure.

**Do not stitch or rescale stacks to make them agree.** I added per-round scale factors measured in the depth overlap so the rounds joined into one curve. Runqi: *"why would the stitch be real? what the fuck is this? take it out."* **Why:** the continuity is then guaranteed by construction, not evidence. **How to apply:** fit each stack on its own; show the disagreement.

**Do not add analyses that were not requested.** Over one task I added an offset correction, a mode diagnostic, a separate LUT-comparison script, local sliding-window fits, and commentary — none asked for. Runqi: *"i ask for a simple task ... you keep adding bs stuff"*. **How to apply:** deliver exactly the requested scope; offer extras in one sentence, don't build them.

**Verify metadata BEFORE iterating on a figure.** I re-ran the same 25 GB analysis four times chasing the depth axis (surface at shallowest slice → 0 → 117/134 → 100 → filename z) instead of reading the metadata properly once. Runqi: *"why are you runnning it some many times? what the fuck is your issue?"* **How to apply:** on a new dataset, dump the geometry fields and reconcile them against the filenames FIRST, in one cheap metadata-only pass, then run.

**Parameters tuned for one session are traps in the next.** `Z_TOP = 10` (written for 260728) silently deleted 115 of 151 slices on 260804. **How to apply:** a hard truncation constant must default to Inf/off and be justified per session.

**Say it plainly when a plot is a plain question.** "Is the attenuation just F(0)/F(z)?" wants yes/no first, then the caveat — not a lecture.
