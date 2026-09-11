---
name: reference_ffmpeg_broken_use_imageio
description: The conda ffmpeg/ffprobe on this machine crashes at launch (0xC0000139) — use the imageio_ffmpeg bundled binary instead
metadata: 
  node_type: memory
  type: reference
  originSessionId: 3d6a9b12-9acd-41a0-b81b-4030e0e68644
  modified: 2026-07-23T04:44:19.071Z
---

**The `ffmpeg` on PATH on this machine is BROKEN.**
`C:\Users\Admin\.conda\envs\dlc310\Library\bin\ffmpeg.exe` (and `ffprobe.exe`)
exit immediately with `0xC0000139` / exit 57 — "entry point not found", a DLL
mismatch. Fixing PATH to the env's `Library\bin` does not help. No other conda
env (`flir`, `cellpose-gpu`, `oasis`, `torch113`) ships one.

**Working alternative** — the static build bundled with `imageio_ffmpeg`
(ffmpeg 7.1, has libx264):
```python
import imageio_ffmpeg
exe = imageio_ffmpeg.get_ffmpeg_exe()
# C:\Users\Admin\.conda\envs\dlc310\lib\site-packages\imageio_ffmpeg\binaries\ffmpeg-win-x86_64-v7.1.exe
```
Always *execute* a candidate (`-version`, check returncode) before trusting it —
presence on PATH is not enough here.

**Consequences:**
- [[project_fix_avi_timing]] calls bare `ffmpeg` and will FAIL as-is until it
  gets the same fallback. [[project_make_viewable_video_copies]] already has it.
- OpenCV's *own* bundled ffmpeg is fine — `cv2` decodes FFV1 and writes video
  without touching the broken binary (see [[feedback_ffv1_video_codec]]).
- `cv2` cannot write H.264 here: openh264 DLL is missing/wrong version, so
  `avc1`/`H264`/`X264` fourccs fail to initialise. Working cv2 writers are
  `mp4v` (.mp4), `MJPG` (.avi), `FFV1` (.avi). For real H.264, shell out to the
  imageio_ffmpeg binary.
