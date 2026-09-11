---
name: project-auto-acq-260910
description: New map acquisition pair auto_acq_260910.m + acq_pause_queue.m — written 2026-09-10/11, control logic tested, NOT yet run on the rig
metadata:
  type: project
---

`D:\RZ_ScanImage_script\auto_acq_260910.m` + `acq_pause_queue.m` replace
`auto_acq_n_shutterOFF.m` for snake-tile map acquisition. Written 2026-09-10/11.

**Status as of 2026-09-11: the control logic is tested off-rig; every hardware
call is UNEXECUTED.** User knows and will test later. Do not describe it as
validated. Shakedown plan agreed: `numCols=2; numRows=1; doLaser=false`.

What changed from the original:
- 10 s/tile of pure idle (`pause(1)` x10 inside the sub-move loop, stage already
  stopped) -> one `settleSec = 0.5` wait at the target. ~13 min back on 81 tiles.
  The 10 sub-moves are KEPT, unchanged. **`motorFastMotionThreshold` is inert on
  this rig** -- `galil_DMC4040.m:9` sets `TwoStep.Enable = false`, so there is no
  fast/slow velocity switch, and the MDF leaves both velocities empty. (An
  earlier note here claimed the sub-moves existed to stay under that threshold;
  that was wrong.) Their real effect is ten short accel-decel cycles instead of
  one 400 um traverse.
- 11 `motorPosition` queries per tile -> 1. Each one is a live `TP` serial query
  to the Galil; sub-step targets now interpolate from tracked `xPrev`/`yPrev`.
- `keepResonantScannerOn` held true across the run, down while paused, re-armed
  `ResArmSec = 2` before resuming. **The one genuinely new hardware interaction.**
- Live control via root appdata `rz_acq_ctrl` (not a global, no handle either
  side holds) — the two windows can start in either order.

Stop semantics: PAUSE = tile boundary, nothing lost. BRAKE = `hSI.abort()`
mid-grab, CONTINUE **re-runs that tile** (the truncated file would otherwise look
complete to the stitcher with fewer z pages than its header claims). ABORT = same
stop, then ends the run, spins the scanner down, and closes the control window.

**Latent bug found in the ORIGINAL while verifying:** `SI.start()` returns early,
before setting `acqState`, when no channel is selected for Display or Save — so
the wait loop falls through and the script marches the stage over all 81 tiles
acquiring nothing, silently. The new script hard-errors at tile 1 instead.

Gotchas learned:
- Editing `acq_pause_queue.m` while its window is open makes every button inert
  (MATLAB clears the function; the callbacks are closures into the cleared
  workspace). Close and reopen after any edit. This looked like "BRAKE is broken".
- `delete(f)` bypasses `CloseRequestFcn`, orphaning the refresh timer; `refresh`
  now stops its own timer when the figure is gone.
- Off-rig rehearsal driver: `...\scratchpad\acq_rehearsal.m` (scratch, not
  pipeline) — fake timer-driven tile loop, same flag logic, for GUI testing.

See [[reference-scanimage-rig-control]].

## Verification done 2026-09-11 (what is and is not proven)
**Cannot run real ScanImage on this machine**: the launcher dies at `loadlibrary('nicaiu')`
— the NI-DAQmx *runtime driver* is not installed, and SI2018b loads it to read the DAQmx
version BEFORE it checks the `simulated` flag, so simulated mode does not get past it either.
`+dabs/+simulated/Stage.m` exists (`simulated.stage` in MotorRegistry), so a simulated MDF is
otherwise viable if the driver is ever installed. **`checkSystemRequirements()` PASSES on
R2021b** — it only demands win64 and >= R2015a, so the MATLAB-version worry is dead.

Verified by reading the ScanImage source:
- `SI.start()` sets `acqState` SYNCHRONOUSLY before any hardware work, so the wait loop is
  sound — EXCEPT one early `return` before that assignment when no channel is selected for
  Display or Save. The original would then march the stage over all 81 tiles acquiring
  nothing, silently. The new script hard-errors at tile 1 instead.
- `set.keepResonantScannerOn` acts IMMEDIATELY (`resonantScannerActivate(v)`) but only
  `if ~obj.active`, and is **PERSISTED to the class data file** — it IS the GUI checkbox
  `cbKeepScannerOn`, survives restart, and after an abort the script leaves it unticked.
- `parkScanner()` runs `resonantScannerActivate(keepResonantScannerOn)` after EVERY
  acquisition, and `arm()` activates then only ESTIMATES the frequency
  (`zzzEstimateScanFreq`; the measure-after-settle branch is commented out) — **there is no
  settle wait in the grab path**, so a cold scanner means the first frames are taken while it
  is still settling. That is what holding it on avoids.
- `ResArmSec = 2.0` matches ScanImage's own number:
  `resonantScannerWaitSettle(max(2, mdfData.resonantScannerSettleTime))`, commented "should be
  settled after 2 seconds". MDF `resonantScannerSettleTime = 0.5` is the zoom-change case.

**The real script executed end-to-end against mocks** (scratchpad `mock/`, `run_real_script.m`,
`run_real_brake.m`; only the config block changed): snake order and filenames correct,
6 grabs for 6 tiles, **9 `motorPosition` reads vs 66 in the original**, scheduled pause ->
scanner off/on/arm, BRAKE -> 7 grabs with `col01_row00` written twice and exactly 1
`hSI.abort()`, ABORT -> stops at the tile, 1 abort, scanner left OFF, `acqState` idle.

Still unverified, all hardware: real scanner timing, whether an aborted tile's TIFF is
truncated-but-readable, and actual stage landing accuracy vs `posTolUm = 2`.
Shakedown agreed: `numCols = 2; numRows = 1; doLaser = false`.

## Files moved 2026-09-11
`D:\RZ_ScanImage_script\` reorganised to `auto-tile\` (the live commented pair),
`stitch\`, `to-send\` (comment-stripped), `commented_originals\` (backup), `archived\`.
The SI install sits alongside at `SI2018bR1_2018-12-19_4a9264c4fc\`.
