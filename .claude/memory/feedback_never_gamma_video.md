---
name: feedback_never_gamma_video
description: Never apply gamma to exported 2P videos — linear intensity only
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 3c68632d-e187-42f7-b413-da2dc5adcb4b
---

NEVER apply a gamma curve to exported 2P video (AVI/movie) frames. Use only a LINEAR
contrast window (percentile clip) to map 16-bit → 8-bit; no `.^gamma`.

**Why:** the user was emphatic ("u never gamma the video, understand?"). Gamma is a nonlinear
intensity distortion that misrepresents the raw calcium signal in a movie; for video the user wants
true-to-data linear intensity (same spirit as the true-quality uncompressed-AVI request — see
[[project_io_population_seq_connectivity]]).

**How to apply:** in any movie writer, scale frames as `(f-lo)/(hi-lo)` clipped to [0,1] then ×255 —
nothing else. Do NOT add a `vidGamma`/`gamma` param to the video path. (Gamma on STATIC display
figures like avg-projection PNGs, e.g. `gamma_disp`, is still fine — this rule is video-only.)
