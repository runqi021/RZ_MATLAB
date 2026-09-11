---
name: reference-fnv-buffer-copy
description: FLIR fnv im.final returns a reused buffer view - must np.array(copy=True) or frames alias
metadata: 
  node_type: memory
  type: reference
  originSessionId: 18f3e613-c0b7-4c78-8b44-95de195de636
---

`fnv.file.ImagerFile.final` (the FLIR FileSDK frame accessor) returns a **view onto an internal buffer that fnv overwrites on the next `get_frame()`**. `np.asarray(im.final, dtype=np.float32)` does NOT copy — so any frame you *store* (append to a list, assign to `prev`) ends up aliasing the same memory and all stored frames become identical to the LAST one read.

Symptoms seen on 260613 thermal:
- Video built by appending frames → a **static "single image"** played repeatedly.
- Frame-to-frame motion energy `|fr - prev|` → **all zeros**.

**Fix:** force a copy whenever you keep a frame across iterations:
`fr = np.array(im.final, dtype=np.float32, copy=True).reshape(H, W)`

Note `thermal_ats_to_mat.py` was already safe because it does `data[i,:] = np.asarray(im.final, ...)` — assigning INTO a preallocated row copies the values immediately. The bug only appears when you retain the view. Scripts fixed: `thermal_inspect_videos.py`, `sync_motion_extract.py`, `sync_combined_video.py`. See [[project_thermal_basler_sync_260613]].
