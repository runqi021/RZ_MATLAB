---
name: File reorganization into subfolders (2026-04)
description: Breathing analysis scripts moved into Breath_summary_timeNphase_260331/; IO scripts into IO_population_activity/; requires path fixes
type: project
---

As of 2026-04-01, the user reorganized scripts from the repo root into subfolders:
- `Breath_summary_timeNphase_260331/` — breathing summary/diagnostic scripts (breath_combined_summary, breath_polar_phase_diagnostic, Batch_breath_trig_combined, breath_cycle_video)
- `IO_population_activity/` — population PSTH and IO analysis scripts

**Why:** Cleaner repo organization, grouping related analysis scripts together.

**How to apply:** Any new scripts placed in subfolders must include `addpath(fileparts(fileparts(mfilename('fullpath'))))` to add the repo root. Don't rely on `contains(path, ...)` guards for addpath — use unconditional addpath instead.
