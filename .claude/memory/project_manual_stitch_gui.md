---
name: project_manual_stitch_gui
description: map_stitch_GUI/ — manual grid stitcher + coord appliers; July 2026 rewrite dropped anchor-fit for a fixed step grid
metadata: 
  node_type: memory
  type: project
  originSessionId: d956f206-4229-4c3c-99eb-1cdf98550d6a
  modified: 2026-08-18T19:04:51.106Z
---

`map_stitch_GUI/stitch_manual_gui.m` — manual 2D tile stitcher for overview maps where tile overlap is too small for the phase-correlation auto-stitcher (`NChan_vol_stitch_zUse.m`) to work reliably. At high zoom the overlap strip is too thin and dim for phase correlation, and the auto-stitcher's global LS solve lets a few bad edges corrupt the whole mosaic.

**2026-07-10 rewrite — anchor-fit is GONE.** The earlier version (2026-05-27) placed tiles at stage-µm/`um_per_px` and had an **Anchor → Fit** button that solved a stage-µm→pixel similarity/affine from 2–3 pinned tiles. That is no longer in the file. The GUI is now a plain regular grid:

    x = (colIdx-1)*stepXpx + manX,   y = (rowIdx-1)*stepYpx + manY

- **There is NO Step X / Step Y GUI control** (the line-11 docstring claiming one was stale; corrected 2026-08-18). `stepXpx/stepYpx` come from `S.stepX_um / S.stepY_um` ÷ (`S.umPerPxBase`/`S.zoom`) at the top of the file. They are auto-derived from filename micron spacing **only** when the filenames lack `colXX/rowYY` — with col/row present, `haveCR` wins and your hand-set `stepX_um` is what is used.
- Controls that DO exist: Load, Flip X / Flip Y, tile-click select, arrow-key nudge (Move mode Tile/Row/Col, Move step px, Shift = 10×), Reset tile / Reset ALL, Show borders, Auto B/C, Save.
- Each tile = AVG PROJECTION over its `avgz` Z-pages (single 2D image). Preview is downsampled (`S.dsPreview`); full res only on Save.
- Input layout: `<datasetRoot>/ch<refChID>/avgz/*.tif`. Output → `<datasetRoot>/matlab_stitch/`: `stitched_ch<ch>_avgproj.tif` + `stitch_manual_coords.mat` (x,y,manX,manY,stepXpx,stepYpx,flips,files — auto-offered on reload).

**2026-08-18 fast port — USE THIS ONE.** `260818-map-registration/stitch_manual_gui_fast_260818.m`
(NOT in `map_stitch_GUI/`). Same maths, same output files, rebuilt on a classic figure+OpenGL
instead of uifigure: ~15x faster redraw (the uifigure version spent 915 ms per arrow nudge, only
10.6 ms of it arithmetic). Adds scroll-zoom / right-drag pan that survive a nudge, `dsFast`/`dsIdle`
(blend coarse while nudging, re-blend full-res after `idleSec`), `projMode` avg|max, and
`stitched_ch<ch>_volume.tif` alongside the avg proj.

**It also added `placeMode`, which supersedes the fixed-step grid.** `'stage'` reads
`SI.hMotors.motorPosition` + `SI.hRoiManager.imagingFovUm` from the RAW tile header (looked up
next to `ch<N>/` and one level up) and places every tile at its true micron position -- no step
size, no zoom, no um/px typed in. Fallback chain: header -> filename `_x###_y###` -> col/row
lattice; the status bar says which one it used. `'saved'` reuses curated coords verbatim.

**Flips are MODE-SPECIFIC and stage mode owns the sign.** `recomputePositions` hardcodes the
measured rig convention `bx=(stageX-min)/umPerPx`, `by=(max(stageY)-stageY)/umPerPx`
-- col -> +stage x, row -> **-stage y** ([[reference_stage_axis_convention]]). `setFlipsFor('stage')`
therefore FORCES `flipX=flipY=false`; a `flipY=true` left over from grid curation mirrors the whole
mosaic, measured at 1197 px of Y error. So in stage mode the `S.flipY` in the settings block is inert,
and you must NOT set a flip to "fix" orientation.

**`S.zUse_um` is in MOTOR-Z MICRONS**, looked up against `SI.hStackManager.zs`; order-insensitive
(`[-50 -200]` == `[-200 -50]`, slab inclusive). `S.zUseProj` (raw page indices) is used only when
`zUse_um` is `[]`.

The prep-layer 90-degree mount rotation ([[project_rotated_mount_sessions]], everything from 260721 on)
must NOT be baked into the stitch -- stitch in stage coords, apply anatomy only for display.

**Coord appliers (all read `stitch_manual_coords.mat`, reuse the GUI's linear-blend weights, swap the `ch<N>` token in folder+filename):**
- `stitch_apply_coords_chan.m(root,srcCh,dstCh)` — re-render another channel's avg proj.
- `stitch_apply_coords_gamma.m(root,ch,gamma)` — gamma each Z page **before** the avg proj, using one shared `globalMax` over every page of every tile so seams stay matched. g<1 brightens dim somata. Writes `..._gamma<g>.tif`, leaves the linear mosaic alone.
- `stitch_apply_coords_volume.m(root,chan,zUse)` — apply coords per Z slice → stitched multi-page **volume** (input for X-Z / depth projections).

**Getting `avgz` in the first place:** `tile_deinterleave_N_groupZavg_260131.m` writes both the deinterleaved Z×T stacks (~1.4 GB/tile/channel) and the avgz. When only the stitcher is needed use `tile_avgz_only_260818.m` instead — avgz only, streams pages with `nextDirectory()` instead of `setDirectory()`, `maxFiles` for resumable partial runs. ~73 s per 2.8 GB / 5400-page tile.

Related: [[project_subvolume_render]], [[reference_stage_axis_convention]], [[project_file_reorganization]].
