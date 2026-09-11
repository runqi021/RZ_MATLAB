---
name: FFV1 video codec — use Python/OpenCV, not MATLAB
description: Basler AVIs use FFV1 (lossless) codec which MATLAB VideoReader can't decode; route all frame-level work through Python cv2 via the dlc310 env
metadata:
  type: feedback
  originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---

For any frame-level video processing on this project's Basler AVIs, use **Python + OpenCV**, not MATLAB's `VideoReader`.

- Basler acquisition (`basler_dual_acq.py`, `basler_acq.py`, `run_basler_n_runs_n_cam.m`) encodes AVIs with the **FFV1** (lossless) codec — confirmed via `cv2.VideoWriter_fourcc(*"FFV1")` at `basler_dual_acq.py:352, 525` and `basler_acq.py:206`. Mono8 grayscale frames, `isColor=False`, `.avi` container.
- MATLAB's `VideoReader` errors out (`Error using VideoReader/initReader ... FFV1`) on these files.
- Python `cv2.VideoCapture` handles FFV1 natively (ffmpeg backend).
- The existing conda env has cv2 available: `C:\Users\Admin\.conda\envs\dlc310\python.exe`.

**Output codec rule depends on use case:**
- **Acquisition output (Basler → disk):** FFV1 lossless. Best DLC input — no codec artifacts on small features (whisker tips, breathing markers). File size ~3–4× bigger than mp4v but disk is cheap. At ≤100 Hz × small binned ROI the encoder easily keeps up; no drop risk.
- **Annotated/labeled output (e.g. DLC review videos):** `cv2.VideoWriter_fourcc(*"mp4v")` with `.mp4` filename. User saw "blue red line glitching" with MJPG on line drawings; mp4v fixed it.
- **Never use MJPG for new code** — it's lossy AND produces color-edge artifacts on solid markers.

**When writing FFV1 output:** always check `writer.isOpened()` immediately after `cv2.VideoWriter(...)`. If OpenCV lacks ffmpeg/FFV1 support, the writer silently returns an un-opened object and `.write()` calls produce an empty file. Raise a clear error instead. (See `basler_dual_acq.py:352–358, 525–531` for the pattern.)

**How to apply:** Any time the task involves reading frames from the raw Basler AVIs (or writing annotated videos), default to a Python script using OpenCV. Only use MATLAB for CSV-based analysis, ROI overlays on still images, or motion-corrected/transcoded AVIs. For new acquisition or output writers, follow the rule above for codec choice.
