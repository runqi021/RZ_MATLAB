---
name: project-machine-dklab-rig
description: DESKTOP-VBHRRVC / user dklab — this IS an acquisition rig with the full 2P hardware AND the analysis pipeline on it; how its paths differ from the Admin machine every older memory was written on
metadata: 
  node_type: memory
  type: project
  originSessionId: 5ca8ecc5-4214-46e0-aa0e-606968d5ff57
  modified: 2026-09-11T21:22:44.209Z
---

Set up 2026-09-11 by transferring the MATLAB pipeline + Claude memory from the
**Admin** machine (`C:\Users\Admin\Desktop\RZ_MATLAB`, MATLAB R2021b). Absolute
paths in older memories are that machine's — check before quoting one.

**This machine is a RIG, not just an analysis box.** Verified by device
enumeration, 2026-09-11: PXIe-6361, PXIe-6341, **NI PXIe-7961R** (the resonant
FPGA), **Interface PEX-292144** (the ALPAO DM card, DI+DO), all `Status = OK`.
NI-DAQmx installed (`nicaiu.dll` present in System32 and SysWOW64 — the exact
thing that made SI unlaunchable on the analysis machine), plus
`C:\Program Files\Alpao` and `C:\Program Files\Galil`.
So the auto-acq shakedown and AO stage 0 CAN be run here. Nothing has been.

## Paths here (vs the Admin machine)
| thing | Admin machine | here |
|---|---|---|
| repo | `C:\Users\Admin\Desktop\RZ_MATLAB` | `D:\Data\RZ\RZ-tools\RZ_MATLAB` |
| ScanImage + AO | `D:\RZ_ScanImage_script\SI2018bR1_...` | `C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc` |
| rig scripts | `D:\RZ_ScanImage_script\auto-tile\` | only the repo copy, `scanimage-control\auto-tile\` |
| MATLAB | R2021b | R2024a **and R2017b** (R2017b installed 2025-10-14, almost certainly the one SI2018b runs under) |

**Hardcoded `D:\RZ_ScanImage_script\...` paths in repo scripts are wrong here.**
Known casualty: `ao_motion_260911\dm_motion_calib_260911.m` `Z2C_FILE`.

**The user has said: do not change anything under the ScanImage folder.** Read
only. There is no `D:\RZ_ScanImage_script` here, so `auto-tile\` has no deployed
twin on this machine — the repo copy is the only copy.

## Claude Code setup
- Memory dir `C:\Users\dklab\.claude\projects\D--Data-RZ\memory\`. The folder
  name comes from the **launch directory**, not the repo — launching from inside
  `RZ_MATLAB` makes a different, empty one. **Always start from `D:\Data\RZ`.**
- `.claude/memory/` in the repo is the travelling snapshot: copy in at session
  start, copy back before committing. A copy, not a link.
- The Admin machine's whole `~/.claude` was also dumped at
  `D:\Data\RZ\RZ-tools\.claude\` — read-only archive; sync the repo snapshot.
- CLI is standalone at `D:\Data\RZ\RZ-tools\claude.exe`, autoupdater off —
  see [[claude-code-installed-standalone-on-d]].

## Not installed (2026-09-11)
- **git is not on PATH and not in Program Files** — the pull/commit half of the
  transfer plan cannot run yet.
- No conda, no ffmpeg → the Python legs (Cellpose, CASCADE, OASIS, thermal
  `.ats`, DLC, Basler) are not runnable here.
- C: ~27 GB free of 212 — install to D: (~2.3 TB).
- `E:\260804_shiver_dbh\map` (the stitcher test set with raw TIFFs) IS here.
- Deployment copies hash-checked 2026-09-11: `scanimage-control\{stitch,auto-tile}\`
  byte-identical to `260910-map-registration\` and `commented_originals\`. The old
  slow shadows still exist in `map_stitch_GUI\` and `stitch_working_260320\` —
  see [[feedback_matlab_root_shadow_copies]].

Related: [[project-auto-acq-260910]], [[project-autostitch-ncc-260910]],
[[project-ao-motion-correction-260911]], [[reference-scanimage-rig-control]].
