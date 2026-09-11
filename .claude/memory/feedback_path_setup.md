---
name: Scripts in subfolders need repo root on MATLAB path
description: Scripts moved into subfolders lose access to repo-root functions (detect_session_fps, helper.dFF_RZ, Chronux, mod/*) — must addpath(repoRoot) and addpath(genpath(mod))
type: feedback
---

When scripts are moved from the repo root into subfolders (e.g., `Breath_summary_timeNphase_260331/`), they lose access to:
- `detect_session_fps.m` (repo root)
- `helper.dFF_RZ` (via `2p_breathing_coherence/+helper/`)
- Chronux functions like `mtspectrumc`, `coherencyc` (`chronux_2_12/`)
- `ghostscript`, `append_pdfs_RZ` (`mod/append_pdfs/`)
- `bluewhitered` and other utilities (`mod/*/`)

**Why:** `fileparts(mfilename('fullpath'))` returns the script's own folder, not the repo root. Relative path lookups break silently. Also `contains(path, repoRoot)` is a bad guard because subfolder paths are substrings of the root path, causing false positives that skip the `addpath`.

**How to apply:**
1. For scripts in subfolders, use `repoRoot = fileparts(fileparts(mfilename('fullpath')))` and call `addpath(repoRoot)` unconditionally (no `contains` guard).
2. The permanent fix is `addpath(genpath('C:\Users\Admin\Desktop\RZ_MATLAB\mod')); savepath;` run once in MATLAB.
3. `dFF_RZ` must always be called as `helper.dFF_RZ(...)` — it lives in the `+helper` package.
