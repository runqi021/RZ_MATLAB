---
name: feedback-matlab-root-shadow-copies
description: Stale copies of a script at the repo ROOT silently shadow the maintained version in a subfolder, because MATLAB's current folder outranks the path
metadata:
  type: feedback
---

MATLAB resolves a function from the **current folder FIRST**, before anything on
the path. The usual pwd here is the repo root `C:\Users\Admin\Desktop\RZ_MATLAB`,
so any leftover copy of a script sitting at the root **wins over** the maintained
one in `analysis_260727\breath_svd\` (or any other subfolder) — even though those
launchers do `addpath(here)` themselves. `addpath` cannot beat the current folder.

Caught 2026-08-28: a root `breath_fixedmetric_to_gui.m` (2026-08-06, pre-full-length)
shadowed the good `analysis_260727\breath_svd\` version and wrote 3000-sample
`breath_pc1.mat` files on a session whose runs were 6000/5095/3000 frames. It read
`S.(METRIC)` (truncated to T = shortest run) instead of `fb_full` + `run_len`.
Moved to `analysis_260727\breath_svd\_superseded_rootcopies\`.

**Why:** the silent failure mode is the dangerous part — the wrong version runs to
completion and produces plausible output. Nothing errors.

**How to apply:**
- Run `which <funcname>` before believing any pipeline output that looks off.
- A stale duplicate at the root is not harmless clutter; move it into
  `_superseded_rootcopies\` (that folder exists for this).
- Fingerprint a stale write by a field the new version adds — here the old copy
  wrote the session-median `fps` into every run instead of per-run `fps_run`.
- Not every root copy is stale: `breath_svd_pc1.m` at the root is NEWER than the
  breath_svd one and is the live edit. Check mtimes before moving anything.

Related: [[feedback-path-setup]], [[project-breath-fixedmetric]], [[feedback-no-inplace-edits]]
