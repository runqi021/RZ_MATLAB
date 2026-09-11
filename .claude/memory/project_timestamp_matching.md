---
name: Timestamp sources for behavior-calcium pairing
description: How to extract reliable timestamps from Basler AVIs (filename) and ScanImage TIFs (epoch header) for pairing
type: project
---

## Reliable timestamp sources

- **Basler AVI**: start time embedded in filename (`__YYYYMMDD_HHMMSSmmm.avi`). Always reliable -- never modified by downstream tools.
- **ScanImage TIF**: `epoch = [YYYY M D H M S.mmm]` in the `ImageDescription` TIFF tag. Read via `Tiff` class (fast, first IFD only) or `imfinfo` (slow fallback for BigTIFFs).
- **`_meta.mat`**: cached metadata from `detect_session_fps.m`. Has `framesPerSlice`, `numSlices`, `channelSave`, `fps`. Use `framesPerSlice * numSlices` for raw frame count per channel.

## Unreliable sources (do NOT use for pairing)

- AVI file modification time: can be corrupted by DLC (observed 1/12 AVIs with wrong mod time)
- cpSAM `F` matrix size: has frame-dropping already applied (e.g., 2970 instead of 3000)

**Why:** Discovered during implementation of `pair_behavior_to_phys.m` -- AVI mod times failed on real data.
**How to apply:** Always parse AVI timestamps from filename, never from file system metadata. For TIF frame counts, use `_meta.mat` first, `imfinfo` as fallback.
