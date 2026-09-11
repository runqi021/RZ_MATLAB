---
name: always-use-helper-prefix-for-helper-package-functions
description: "Functions in +helper/ (at repo root) must be called as helper.funcName(), never bare funcName(); never duplicate the helper as a standalone script"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: be23b187-1863-428f-9dbd-81b57a9328c2
---

`dFF_RZ`, `stackDFF`, `oasis_deconv_and_dff_AR1`, `phase_jackknife` live in `+helper/` at the **repo root** (`C:\Users\Admin\Desktop\RZ_MATLAB\+helper\`, moved here 2026-05-18). They must always be called with the package prefix: `helper.dFF_RZ(F)`, not `dFF_RZ(F)`.

**Why:** `+helper` is a MATLAB namespace package — calls only resolve when the *parent* folder (repo root) is on the path. The user moved it to repo root and deleted the duplicate `dFF_RZ_251202/` copy specifically so there is ONE canonical dFF function. Do not re-introduce a local fallback / bare-name duplicate / "dispatch" wrapper. If `helper.dFF_RZ` fails, **fix the helper**, do not fork it.

**How to apply:**
- Any new script that needs dFF: `addpath(repoRoot)` then call `helper.dFF_RZ(F, 'FPS', fps, 'BaselineWinSec', winSec)` directly.
- Never write a `dFF_RZ_dispatch` / `dFF_RZ_local` / `if exist(...) else ...` fallback wrapper.
- Pre-existing `addpath(fullfile(repoRoot, '2p_breathing_coherence'))` lines are now redundant for helper access but harmless — leave them unless asked.
- When fixing legacy scripts, search for `[^.]dFF_RZ\(` (bare call) and rewrite with `helper.` prefix.
