---
name: feedback-cell-ids-not-stable
description: Ventral archive cell IDs are POSITIONAL and get renumbered on every registry rebuild; the "new sessions append at the end" behaviour is a dir() artifact, not a property
metadata:
  type: feedback
---

`Ventral_surface_event_latency_260811.m` numbers cells with
`unique({OBS.cellKey},'stable')` over `dir(gdir,'**','ca_spike_data.mat')`.
The ID is therefore a POSITION in scan order, not an identity.

**The append is luck.** Vgat/0730 -> 277-281, Vglut2/0824 -> 282-297,
ChAT/0826 -> 298-299 all landed at the end only because `dir('**')` on Windows
enumerates freshly-created folders last. Once those folders aged, the 2026-08-28
rebuild moved ChAT/0826 from 298 -> 23 and Vglut2/0824 from 282-297 -> 57-78.
Measured: of 299 old cells, **22 kept their ID, 256 moved, 21 vanished**, and
**125 of 288 rendered figures** then carried an ID resolving to a different cell
(`ChAT_0826_cell298.png` -> cell 298 is now `Sert/0721#c27`).

**What is safe:** the GUI (`temporal_phase_cell_gui_260812`) and the batch render
(`per_cell_summary_260812`) both read the live `event_latency_data.mat` and both
number through `apply_cell_merges_260814` keyed on `min(origCells)`, so they always
agree WITH EACH OTHER. `cell_merge_overrides_260814` matches on `(recName, roi)`,
so hand merges survive a renumber. Only figures/CSVs rendered against an OLDER
registry go stale.

**How to apply:**
- Never assume an ID from an existing figure or CSV survives a rebuild. Translate
  via `cell_key`, which IS stable, before trusting any saved selection
  (`cell_classes_260816.csv`, `ppt_260815_picked_cells.csv`, `popsel_260816\`).
- After any rebuild, diff old vs new `{CELL.key}` before rendering.
- Back up `event_latency_data.mat` first; rebuilds overwrite it in place.
- Only `polar_event_latency_260811.m` reads `S`; the other ten consumers load
  ONLY `CELL`/`OBS`/`REC`.

**OPEN 2026-08-28:** `cell_link_sources` in the script is missing its
`'Vglut2','0824'` row although `D:\Ventral_surface_summary\Vglut2\0824\cell_pooled\
cell_link.mat` exists (2026-08-25), so a rebuild silently drops that curation
(16 grouped cells -> 22 per-ROI). The script is UNTRACKED in git, so the removal
left no history. Restore the row before the next rebuild.

Related: [[project-ventral-selected-cells-260816]], [[project-ventral-summary-archive]],
[[project-cell-link-unapplied-curation]]
