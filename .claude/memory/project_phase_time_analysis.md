---
name: Phase and time analysis plan
description: Overall plan for breathing-calcium analysis at three levels (per-ROI, per-FOV, per-experiment) in both phase and time domains
type: project
---

User wants phase AND time analysis at three levels:
- Per ROI
- Per FOV (field of view / recording)
- Per experiment (cross-session summary)

**Why:** To characterize breathing-calcium coupling at multiple scales for publication figures.

**How to apply:** When building new analysis scripts, always consider which level (ROI/FOV/experiment) and which domain (phase [0,4π] or time [seconds around breath peak]).

Current scripts:
- `breath_phase_summary_260328.m` — experiment-level PHASE analysis (8-panel 2×4 grid)
- `breath_spike_raster_260328.m` — per-FOV PHASE analysis (4 separate figures)
- `Batch_breath_trig_combined_260328.m` — per-FOV TIME analysis (giant combined figure)
- Next: experiment-level TIME analysis script (breath_time_summary)
