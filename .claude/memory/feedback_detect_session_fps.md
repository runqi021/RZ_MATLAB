---
name: Always use detect_session_fps for FPS and zoom
description: Never rely on UI defaults or hardcoded values for FPS/zoomFactor — always call detect_session_fps() to auto-detect from TIF metadata / _meta.mat
type: feedback
---

Always use `detect_session_fps(folderPath)` to get FPS and zoomFactor from TIF metadata or `_meta.mat`. Never rely on UI defaults or hardcoded values — sessions have different zoom levels (4x, 8x, 12x, etc.) and the default will be wrong.

**Why:** The zoom factor directly affects physical-size calculations (um_per_px = 1.7778 / zoomFactor). Using a wrong default (e.g., 4x when the session is 12x) produces completely wrong crop sizes and scale bars. The user had to debug this manually.

**How to apply:** In any GUI or script that needs FPS or pixel size, call `detect_session_fps()` at session load time. It returns `[fps, scan_meta]` where `scan_meta.zoomFactor` and `scan_meta.pixelSize_um` are available. File: `detect_session_fps.m`.
