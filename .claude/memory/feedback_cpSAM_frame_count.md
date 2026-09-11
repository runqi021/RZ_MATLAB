---
name: Do not use cpSAM F matrix for raw frame count
description: cpSAM_output.mat F matrix has TossFrames already removed -- use _meta.mat framesPerSlice instead
type: feedback
---

Never use `size(cpSAM.F, 1)` when you need the raw TIF frame count. The cpSAM pipeline drops `TossFrames` (default 30) before extracting fluorescence, so `F` is shorter than the original recording.

**Why:** Using cpSAM F size for AVI trimming would cut 30 frames short (e.g., 2970 instead of 3000). Discovered during `pair_behavior_to_phys.m` development.
**How to apply:** Use `_meta.mat` fields `framesPerSlice * max(numSlices, 1)` for raw per-channel frame count. Fallback: `numel(imfinfo(tifPath)) / nChannels`.
