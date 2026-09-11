---
name: Basler dual-camera acquisition GUI + Python backend
description: Working architecture of Bassler/basler_gui.m + basler_dual_acq.py — the live 2P-synchronized behavior-video acquisition system
type: project
originSessionId: 59fed66a-d20f-4e6c-9b4a-706cdbd60cee
---
This is the **primary data-acquisition tool** for behavior video in sync with the 2P scope. It is important and treated carefully by the user.

## Files (`Bassler/` folder)

- `basler_gui.m` — MATLAB `uifigure` front-end
- `basler_dual_acq.py` — Python backend (pypylon + OpenCV)
- Python interpreter path is hardcoded: `C:\Program Files\Python314\python.exe` (not the dlc310 conda env — a plain Python 3.14 with pypylon + opencv-python installed)
- `requirements.txt` — Python deps

## Two-stage architecture

GUI spawns the Python script via `parfeval(backgroundPool, @system, 2, cmd)`. Python does all camera I/O; MATLAB only builds the CLI args, polls for completion, and parses stdout for frame counts.

Python subcommands:
- `list` — enumerate Basler cameras (called at GUI startup to populate SN fields)
- `preview` — single-camera freerun preview window with ROI picker
- `acquire` — session-based dual-camera acquisition (the GUI Record button)
- `acquire-nruns` — untriggered free-run periodic acquisition (dur record + idle × nRuns); driven by `run_basler_n_runs_n_cam.m` (no GUI)

## Trigger modes (three; GUI dropdown = hardware / burst / freerun)

(2026-06 update) GUI now offers **freerun** too — added to the trigger dropdown; `onTriggerChange` grays Burst-N for hardware+freerun. Because free-run has **no trigger gap to end a run**, the bottom panel has a **Free-run dur (s)** field (`freerunDurField`, default 300): when ALL selected cams are freerun, `onRecord` launches `acquire --session <dur>` (no `--single-run`) → one continuous file of that length (or until Stop). When any cam is triggered, it keeps the old `--split-gap G --single-run` arm-and-split behavior (default 4 h session). hardware/burst are Line3 falling-edge as before.

1. **`hardware`** — `TriggerSelector=FrameStart`, `TriggerMode=On`, `TriggerSource=Line3`, `TriggerActivation=FallingEdge`. **One frame per pulse.** FPS field is used only as exposure cap + AVI metadata.

2. **`burst`** — `TriggerSelector=FrameBurstStart`, `TriggerMode=On` on Line3 falling edge, then `AcquisitionBurstFrameCount=N`, `AcquisitionFrameRate=fps`. **N frames at internal fps per pulse, then idle until next pulse.** Designed for "slow 2P × fast behavior" sync: e.g. 2P runs at 5 Hz, each pulse triggers a burst of 30 whisker frames at 150 Hz, keeping the two streams synchronized at the burst boundary.

## Exposure control + fps-mismatch preflight (2026-06-05 update)

Exposure was already adjustable everywhere (GUI `Exp µs` field `controls.exposure`; n-runs `cams(k).exposure`; python `--camN-exposure` → `cam.ExposureTime`). What was added is the **achievable-fps WARNING** (user set 200 Hz, silently got 111–137 Hz — see [[project_orofacial_crop_dff_analysis]]).

- **Why fps varied:** max fps = min(readout∝rows, 1/exposure, link-bandwidth∝pixels). Those recordings were mostly **exposure-limited** (~8–9 ms ⇒ ~111–125 Hz), not 200. Rows are the main *readout* lever; width helps only bandwidth.
- **`basler_dual_acq.py` `open_camera`:** the achievable-fps warning now fires in **ALL trigger modes** (was burst-only). Emits machine-readable `RESULTING_FPS=`, `REQUESTED_FPS=`, `EXPOSURE_US=`, `FPS_MISMATCH=0/1` from `cam.ResultingFrameRate` + an exposure-vs-fps check (`1e6/exposure < fps`).
- **New `check` subcommand:** opens ONE camera with the real roi/exposure/bin/fps/trigger via `open_camera`, prints those tokens, closes — no grab. Preflight only.
- **n-runs `.m`:** after preview, runs `check` per enabled cam, prints achievable fps in the plan; if mismatch → `[Y/n]` flips to **`Proceed anyway? [y/N]`** (default abort).
- **GUI `.m`:** `onRecord` runs `check` per selected cam BEFORE locking/launching; on mismatch shows `uiconfirm` (default Cancel) and aborts if cancelled.
- Parse tokens with `regexp(out,'FPS_MISMATCH=(\d)','tokens','once')` / `'RESULTING_FPS=([\d.]+)'`.

## Critical Line3 config (ace USB GPIO gotcha)

Line3 on Basler ace USB is a **bidirectional GPIO**, not a fixed input. It only receives triggers when `LineMode=Input`. Without this, triggers silently disappear — no error, no frames. `open_camera()` sets:

```python
cam.LineSelector.Value = "Line3"
cam.LineMode.Value = "Input"          # required
cam.LineDebouncerTime.Value = 5.0     # µs, rejects ringing
```

Then clears `TriggerMode=Off` on `AcquisitionStart`, `FrameBurstStart`, and `FrameStart` **before** applying the desired selector — otherwise a prior session's `FrameBurstStart=On` starves a fresh `FrameStart` trigger. Then prints `Line3 input state: HIGH/LOW` + list of active trigger selectors as diagnostics.

## MATLAB anonymous-function closure pitfall (fixed)

The original code set `ValueChangedFcn` for the trigger dropdown on the same line that created the dropdown, when `controls` was still partial (only `.name` existed). Anonymous functions capture variables by value at creation time, so the callback froze a snapshot missing `controls.burstN`/`.fps`. Cam1 silently failed to enter burst. Fix: wire all `ValueChangedFcn`s **after** the `controls` struct is fully populated, at the end of `createCameraPanel`.

## RAM-first lossless pipeline

Frames never touch disk in the grab loop.

```
sensor → pylon transport buffer (MaxNumBuffer=200, ~1.3s @150Hz)
       → Python queue.Queue (unbounded)
       → background writer_thread → FrameWriter → .avi
```

(2026-06 update) `FrameWriter` class replaced the inline `cv2.VideoWriter` in BOTH writer threads. It pipes raw frames to a **multithreaded ffmpeg FFV1** subprocess (`-c:v ffv1 -level 3 -slices 24 -slicecrc 1 -threads 0`), falling back to single-thread `cv2.VideoWriter(FFV1)` only if ffmpeg isn't found. **Why:** single-thread cv2 FFV1 (~25 MB/s) couldn't keep up with 200 Hz×800×400 (~64 MB/s) → the unbounded queue swelled to ~19 GB during a 5-min shoot and `write_q.join()` blocked ~10 min draining it after Stop. Multithreaded ffmpeg (~150-300 MB/s on 8+ cores) keeps up in real time → save tail ≈ 1-2 s regardless of length, RAM stays flat. **Dependency:** needs `ffmpeg.exe` on PATH (or env `BASLER_FFMPEG`); startup prints `[encoder] multithreaded ffmpeg FFV1: <path>` or a WARNING if it fell back. Bonus: killing python (GUI Stop) closes ffmpeg's stdin → it finalizes the file cleanly (cv2 would leave a corrupt AVI).

- `close_run()` calls `write_q.join()` before the AVI is closed → every queued frame is flushed before `timestamps.csv` is written.
- Session end calls `writer.join()` → no frames abandoned.
- After each session, Python prints `WARNING: pylon dropped N frame(s)` if `Statistic_Failed_Buffer_Count > 0` — this is the self-check for "did we actually stay lossless". Non-zero means bump `MaxNumBuffer` further or check USB bandwidth.

## Codec + output layout

- **FFV1** lossless `.avi` (current acquisition, via multithreaded ffmpeg — see RAM pipeline section). Lossless for breathing/whisker tracking; read by cv2/DLC. NOT viewable in Windows Media Player or ImageJ (ImageJ AVI reader has no FFV1) — use **VLC**, or transcode to mp4 (`mp4v`/lossless x264) for scrubbing. (The earlier MJPG-for-WMP plan was superseded by lossless FFV1.)
- Folder structure (per run gets its own subfolder, never pooled):
  ```
  <saveDir>/
    <camName>/
      <camName>_<YYYYMMDD_HHMMSS>_run001/
        <camName>_<YYYYMMDD_HHMMSS>_run001.avi
        timestamps.csv      # frame_idx, camera_timestamp_ns, wall_time_s
  ```
- A new run starts on the first trigger after a `--split-gap` (default 10s) silent gap. Single-run mode (the GUI's default) exits after the first completed recording.

## GUI panel behavior

- Each camera panel has an **Enable** checkbox at top-left. Unchecked → every field in that panel is grayed out and the camera is skipped.
- Whichever camera is enabled is passed to Python as `--cam1-*`. If both are enabled, the second becomes `--cam2-*`. This means single-cam mode works regardless of which physical slot holds the only camera.
- Per-panel **Preview** button runs freerun preview for that cam alone (with `c`/`r`/`q` key handling) and auto-populates the ROI field from what you drew.

## What to never break

- `LineMode=Input` + pre-clear all trigger selectors — removing these brings back the "no triggers detected" bug.
- Callback wiring deferred to end of `createCameraPanel` — removing this brings back the cam1 burst bug.
- `MaxNumBuffer=200` — smaller values drop frames under any brief grab-loop stall.
- Writer-thread `write_q.join()` inside `close_run()` — without it, trailing frames are lost at run boundaries.
