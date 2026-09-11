---
name: reference-matlab-opengl-export-stubs
description: "MATLAB exportgraphics intermittently fails on this machine's NVIDIA driver, writing 9,223-byte PNG stubs; force opengl('software') in any long figure batch"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 7abf9dd7-af19-40fb-a912-10c9692c4000
  modified: 2026-08-15T05:42:10.598Z
---

MATLAB's **hardware OpenGL** path throws, intermittently, inside `exportgraphics`:

```
com.jogamp.opengl.GLException: Error making context 0x... current
  on Thread AWT-EventQueue-0 ... WindowsWGLContext ... NVIDIA 595.71
  at com.mathworks.hg.util.HGGetframeOutputHelper.generateOutput
```

Measured rate on the 2026-08-14 per-cell batch: **3 casualties in 269 figures** (~1%), non-deterministic and NOT data-dependent — the same cells re-rendered fine immediately afterwards.

**Three distinct failure modes, in increasing order of danger**

1. Caught by `try/catch` → the cell is skipped, and **any previous run's file stays on disk**. A folder listing then looks complete while holding a stale figure. This is how `Vglut2_0810_cell053` kept an 08-13 figure through an 08-14 rerun.
2. A **9,223-byte PNG stub** is written — file exists, looks like output, is blank. Two cells hit this. The byte count appears to be constant for the stub, so it is a usable detector.
3. The whole MATLAB session dies. No crash dump is produced (Java catches it), so the only evidence is the process start time being *later* than the batch's first output.

**Only the raster path is affected.** The PDFs of the same failed figures came out at full size (208-222 KB), because `exportgraphics(..., 'ContentType','vector')` does not touch the driver. A healthy PNG here is 220-410 KB.

**Fix** — at the top of any batch that exports many figures:

```matlab
try, opengl('software'); catch ME, fprintf(2,'%s\n', ME.message); end
```

Verify with `opengl('data')` → `Renderer: GDI Generic`, `Software: 1`. Costs roughly 3x per figure (~30 s vs ~10 s for the per-cell summary) and does not lose cells. This is in `per_cell_summary_260812.m`; do not remove it for speed.

**Always verify a batch by file, never by count.** Check every expected output was written *today* AND is above a plausible size floor. Neither test alone catches all three modes: mode 1 passes the size test, mode 2 passes the existence test.

Related: [[feedback-matlab-offscreen-figure-reuse]] — the other way this project has silently produced blank exports.
