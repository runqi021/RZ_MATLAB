---
name: Never modify existing working scripts in-place
description: When adding a new analysis/feature, always create a NEW standalone .m file first; do not edit the user's existing working script
type: feedback
originSessionId: 0ec053e8-fba2-425a-9aa8-58264c9cddb3
---
When the user asks for a new analysis/feature/visualization to be added to a project, **always create a NEW standalone .m file**. Do not edit fields, parameters, or add new sections inside the user's existing working scripts (e.g., `breath_combined_summary_260329.m`).

**Why:** On 2026-04-13, user explicitly objected: "i hate u for thsi shit. can u not touch my code like this? why cant u just write a new file first." Their existing scripts produce known-good outputs they rely on; modifying the inputPath, struct fields, or appending blocks risks losing their existing figures and forces a revert. The user's existing .m scripts are de facto stable artifacts even when they look like quick analysis files.

**How to apply:**
- Default to creating a new file like `breath_<feature>_NNNNNN.m` in the relevant subfolder (`Breath_summary_timeNphase_260331/`, `IO_population_activity/`, etc.)
- The new script should be fully standalone: own param block, own load logic (read raw SAM/CSV/MAT, OR load already-saved `*.mat` outputs from the original script's `breath_combined_summary/` dir if available — e.g., `fov_map.mat`, `coherence_data.mat`, `roiQC_master_summary.mat`).
- Only modify an existing script if the user explicitly says "modify X" / "add this to X.m" / "edit my script."
- If unsure whether to extend or create new, ask before editing.
- Reverting in-place edits requires re-running the original to regenerate outputs the user may have lost — costly and fragile.

**Saved outputs that downstream scripts can re-use** (so you don't need to rerun the heavy analysis):
- `<inputPath>/breath_combined_summary/fov_map.mat` — `fov_map{i,1}` = sessName, `fov_map{i,2}` = vector of sig ROI IDs (per coherence)
- `<inputPath>/breath_combined_summary/coherence_data.mat` — pooled coherence vectors (all_th_k_*, all_r_k_*, all_is_sig_*)
