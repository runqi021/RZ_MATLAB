---
name: project-repr-export-video
description: "repr_export_RZ.m — the calcium + breath video + moving trace composite exporter; N-ROI generic, ylim keyed by ROI, preview-vs-export trap"
metadata: 
  node_type: memory
  type: project
  originSessionId: f145e055-fef7-4364-bc1d-f95f7ce506ec
  modified: 2026-08-16T08:36:51.930Z
---

`repr_export_RZ.m` (repo root) is THE script that merges the calcium movie,
the Basler breathing-camera video and the scrolling traces into one MP4.
`repr_export_RZ_260811.m` is the same file under the dated name — keep them in
sync or delete one. Neither is tracked in git.

## Cell-id variant (2026-08-16)
`analysis_260727\coh_ca_breath\repr_export_cell_260816.m` is a standalone copy
where you **type ids and nothing else** — `fovPath` defaults to `''` and the
RECORDING FOLLOWS FROM THE ID. Identity is the GUI's
(`temporal_phase_cell_gui_260812.m`, same folder): registry
(`event_latency_data.mat` + `cell_merge_overrides_260814`) → session
`cell_link.mat` → plain ROI, registry authoritative wherever it covers the folder.

`idMode='auto'` reads a number as a cell id first, then as a ROI index in
`fovPath`; a registered ROI is relabelled with its cell id, so in Vgat/0730/deep
**`279` and `29` are the same neuron** (cell 279 = ROI 29). Verified against the
live registry: that FOV registers only cells 278/279/280 = ROIs 13/29/33.

**Precedence trap (fixed 2026-08-16):** small ids are usually BOTH a valid ROI in
the typed folder and a cell id somewhere else in the archive — `7` was silently
jumping to `roi3_R_-1000_140_2x_...`. The typed folder now wins, gated on the
folder's real ROI count from cpSAM (`matfile` → `size(m,'F')`), so `279` in a
40-ROI FOV is still recognised as a cell id and moves you.

Clips go to ONE collected `outDir`, named `<recName>_cell279_t35-65s.mp4`.
Downstream keys on the CANONICAL id (cell id where there is one, else ROI index),
so tile label, filename and `YLIM_BY_ID` are all the same number. Ids spread over
different recordings is the one hard error — a video cannot pool.

NOTE the GUI's `get_registry` sets `CACHEKEY = regFile` but compares it to `key`
(path|bytes|datenum), so its cache never hits; the 260816 copy sets `CACHEKEY = key`.

Related but different: `repr_video_gui_RZ.m` is the older interactive sibling
(play/scrub/export). `io_repr_video_RZ.m` is IO calcium with NO behavior cam.
`whisk_breath_video_gui_RZ.m` / `whisk_video_gui_RZ.m` are behavior-only and
produce the `*_breathwhisk` / `*_traceclip` files.

## Layout
Tile stack = 1 breath tile + `numel(roiIdxs)` Ca tiles, left column; scrolling
trace panel on the right. As of 2026-08-11 the tile count follows `roiIdxs` —
it used to be hardwired to exactly 3 Ca ROIs (`for k = 1:3`, `4*TILE_SZ`,
`rowTop(4)`, and the `outFile` sprintf). `roiIdxs = [12]` gives the 2-video
layout.

## Y-limits
- `YLIM_AUTO_MODE = 'robust'` — floor at `prctile(yy,1)`, ceiling at
  `max(movmedian(yy,3))`. The old `min`/`max` let one hot frame flatten every
  real transient. `'minmax'` restores the old behaviour.
- `YLIM_BY_ROI = [roiIdx lo hi; ...]` — keyed by ROI INDEX, not row position.
  The old positional `TRACE_YLIM_MANUAL` silently re-pointed each row at a
  different neuron whenever `roiIdxs` changed.
- Chosen limits are printed every run; copy a line into `YLIM_BY_ROI` to tweak.

## Workflow
`MAKE_LAYOUT_PREVIEW_ONLY = true` renders one frame and stops → eyeball →
adjust ylim → flip to `false` to export. x-window is always manual
(`t_start` / `t_end`).

## TRAP: preview passing does NOT mean export will pass
The preview block `return`s at the end of section 7c. Sections 8+ (outFile
naming, VideoWriter loop, ffmpeg transcode) are never touched by a preview
run, so bugs there only appear on the real export. Check section 8 separately
after any edit.

## Behavior video loading
- AVI discovery tries `Basler_*.avi` → `cam*.avi` → any `*.avi`, excluding
  `_mjpeg.avi` (its own transcode) and `_labeled.avi` (DLC). Old rigs wrote
  `Basler_*`, the current Basler GUI writes `cam1_*`/`cam2_*`.
- FFV1 is transcoded once to `<avi>_mjpeg.avi` because `VideoReader` cannot
  read FFV1; the transcode is cached and reused.
- The transcode stamps the TRUE rate: `ffmpeg -y -r <fps_ca> -i in ...`.
  **`-r` must go BEFORE `-i`** so it reinterprets the input timebase and all
  frames pass through 1:1. On the output side it resamples — with a 300 fps
  header and a 30 fps target that silently drops 9 of every 10 frames.
  Verified on 260806_sst: 6000 frames in → 6000 out, 200.00 s, pixels
  identical to the old cache at frames 0/1500/3330/5999.

## TRAP: the AVI header fps is frequently bogus
260806_sst declares **300 fps** for 6000 frames → a 20 s file, while the
calcium is 199 s. The camera is 2P-frame-triggered, so behavior frame index
= calcium frame index + `toss_beh`, 1:1. Seek by FRAME **in the file's own
timebase**:

    vBeh.CurrentTime = (f0 - 1 + toss_beh) / fps_beh;   % NOT / fps_ca

Dividing by `fps_ca` asked for 111 s in a 20 s file → "Expected CurrentTime
to be a scalar with value <= 20". It only worked before because earlier
sessions happened to declare 30 fps. Same family as the bogus 200 fps
orofacial headers and `camera_fps=50.7`. `selectBreathCropInteractive` and
`estimateBreathClipLimits` are safe — they build sample times from
`v.Duration`/`v.FrameRate`, entirely inside the file's timebase.

## ffmpeg
The script shells out to ffmpeg for the AVI→MP4 transcode. Bare `ffmpeg` on
PATH here is a broken conda build that exits 0xC0000139 before doing anything —
see [[reference-ffmpeg-broken-use-imageio]]. `resolve_ffmpeg()` (local function)
now finds the working binary under
`%USERPROFILE%\.conda\envs\dlc310\...\imageio_ffmpeg\binaries\ffmpeg*.exe`.

Other gotcha: raw Basler AVIs are FFV1 and `VideoReader` cannot open them —
transcode first (`make_viewable_video_copies.py`).
