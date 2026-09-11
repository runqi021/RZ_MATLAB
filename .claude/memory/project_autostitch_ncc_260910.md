---
name: project-autostitch-ncc-260910
description: "NChan_vol_stitch_zUse_260910.m -- why the old auto-stitcher failed, and the bounded-NCC + stage-prior + IRLS design that replaced it"
metadata:
  type: project
---

**`260910-map-registration\NChan_vol_stitch_zUse_260910.m`** — rewrite of the registration
half of `NChan_vol_stitch_zUse.m`. File discovery, both naming conventions, channel pairing,
zUse, blending and the TIFF writer are unchanged.

## Why the old one failed (root cause, measured)
It cut the correlation patch **at the initial guess** and phase-correlated it. A patch of
width P can only express shifts of +/-P/2 before the FFT wraps, so the guess had to be right
to within half the overlap it predicted — the one thing it cannot be.
260909 map: guess 0.8889 um/px -> step 450 px -> patch 62 px -> range +/-31 px; true step
418.75 px, so the needed correction was **-31.2 px = 101% of the range**. Every horizontal
edge aliased or failed `psrMin`, fell back to the wrong guess at weight 0.2, and dragged the
solve onto a lattice 31 px/step too coarse.

## The design that works
1. **Template cut INDEPENDENTLY of the guess** — a strip from the leading edge of tile B,
   width = `TemplateFrac` x expected overlap, searched inside a generous strip of tile A.
   The guess only has to land within `SearchRadiusPx`.
2. **`normxcorr2`, not phase correlation** — linear (no wrap to alias into), normalised, and
   the peak is a correlation coefficient so the threshold means something.
3. **Stage prior on EVERY tile**, not an anchor on one. Fixes the gauge AND stops drift.
4. **IRLS (Huber) global solve** — 36 tiles = 72 unknowns vs 60 edges = 120 equations; the
   redundancy is what lets a bad edge demote itself instead of bending the mosaic.
5. **Loop closure on every 2x2 cycle** — `(A->B)+(B->D)-(A->C)-(C->D)` must be 0. Needs no
   ground truth and catches mutually-inconsistent edges that each look confident alone.

## Traps found the hard way (all cost a full debug cycle)
- **USE THE VESSEL CHANNEL AS REFERENCE.** On 260804 ch3 gives NCC **0.94** and every
  projection variant agrees on the same dx; ch1 gives 0.27-0.49 and the variants disagree.
  I had flipped refChID to 1 and everything collapsed. **ch1 on 260804 has NO visible cells
  at all** — diffuse haze with vessel silhouettes, contrast index 0.02 on all 36 tiles vs
  ch3's 0.02-0.14. (Do not assume ch1 = somata from another session; LOOK.)
- **A high-pass prefilter HURT** (0.941 -> 0.728). Added it to fix a problem that only
  existed because the wrong channel was being correlated. Default 0.
- **PSR is a PHASE-CORRELATION statistic.** An NCC surface is smooth and broad, so a perfect
  match sits only a few SD above its own neighbourhood: `psrMin = 4` rejected 142/144 edges
  that had already passed NCC >= 0.3. Gate on the peak VALUE plus a uniqueness margin
  (peak minus best rival outside `ExclRadiusPx`).
- **Report the LARGEST component, not connectivity to tile 1.** Tile 1 is the acquisition
  start = a CORNER = the likeliest tile to be empty. It said "35 of 36 disconnected" about a
  mosaic whose other 35 were fine.
- **REBUILD THE PRIOR AT THE MEASURED SCALE before solving.** The prior starts at the header
  um/px, which is 8% wrong here; pinning a dark tile to a stretched prior fights every good
  edge. Rebuilding it dropped the post-solve residual **3.99 -> 0.73 px median**. Per axis:
  x and y differ ~6%.
- **Dark tiles**: contrast index `(p99.5-p50)/p50` on the reference projection; below
  `DarkContrast` the tile is pinned to the stage (`DarkPriorBoost`). 11 of 36 on 260804.
  Pinning beats leaving them free — a dark tile still emits one spurious edge, and a single
  unopposed bad edge drags an otherwise unconstrained tile a long way.

## Datasets
- **`E:\260804_shiver_dbh\map`** — THE test set: 6x6=36, **raw TIFFs present**, zoom 1.2,
  FOV 758.55 um, 1.4815 um/px header, 600 um step -> 107 px overlap. avgz 512x512x54.
  Result: 54/60 edges, NCC med 0.93, loop closure 1.27 px, post-solve 0.73 px, clean mosaic.
- `C:\...\260909_ChAT_g8m_Shiverer\map` — 9x9=81, no raw TIFFs, zoom 2, 400 um step.
- **Filename `_x###_y###` vs `SI.hMotors.motorPosition`: max |diff| 0.49 um** over 36 tiles
  (sd 0.30). The commanded position IS the motor position here, so filenames are a fine
  seed; raw headers are only needed if that ever stops being true.

## Measured scale disagreement with the header (open)
260804 edges say **1.6004 um/px (x 1.616, y 1.524)** vs header 1.4815 — **+8%**.
260909 curation says 0.9462 vs 0.8889 — **+6%**. Same direction, different magnitude.
See [[project-structural-stack-260910]] for the curated-map measurement.

Related: [[project-manual-stitch-gui]], [[reference-stage-axis-convention]].

## Folder renamed 2026-09-10: `260818-map-registration` -> `260910-map-registration`
It now holds ONLY the three new scripts; the old `NChan_vol_stitch_zUse.m`,
`stitch_manual_gui_fast_260818.m` and `tile_deinterleave_N_groupZavg_260131.m` were
deleted from it. Nothing breaks — every script resolves its own location with
`mfilename('fullpath')` — but a few COMMENTS still name `stitch_manual_gui_fast_260818`.

**THE SLOW ORIGINALS STILL EXIST ELSEWHERE** and will be found first if the path order
favours them: `map_stitch_GUI\tile_deinterleave_N_groupZavg_260131.m`,
`stitch_working_260320\tile_deinterleave_N_groupZavg_260131.m`,
`map_stitch_GUI\NChan_vol_stitch_zUse.m`, `stitch_working_260320\NChan_vol_stitch_zUse.m`,
`map_stitch_GUI\stitch_manual_gui.m`. See [[feedback_matlab_root_shadow_copies]] —
MATLAB's current folder outranks addpath, so a stale copy can silently run instead.

## tile_deinterleave_avgz_260910.m — the speed work (2026-09-10)
**TWO WRONG DIAGNOSES BEFORE THE RIGHT ONE. Measure, ALTERNATING, before claiming.**
- `setDirectory` per page is NOT the bottleneck: 29.8 vs 27.7 s for nextDirectory,
  **1.08x**, stable across three alternating rounds. libtiff caches the directory
  offsets it has visited and the pages are read ASCENDING, so the rewalk never happens.
  The quadratic blowup is real only for random/backward access — which is what a first
  probe accidentally sampled, producing a bogus "27x".
- Write-cache saturation is NOT it either: 10 x 1 GB to E: ran 130 MB/s dead flat.
- **It is MATLAB's `Tiff` WRITER**, 10-16x slower than raw bytes. One 2700-page 1.32 GB
  file: raw fwrite 10.5 s (129 MB/s) | Tiff 1 strip/page 103.8 s | Tiff 64-row strips
  169.2 s (the ORIGINAL setting) | Tiff BigTIFF 134.7 s | hand-laid TIFF 22.6 s.
  **BigTIFF is SLOWER** — only reach for it past the 4 GB classic-TIFF limit.

Fixes, both on the write path: `RowsPerStrip = ImageLength` (1.63x) and a direct TIFF
writer — pixels streamed from byte 8, all IFDs appended at the end (4.6x).

**THE DIRECT WRITER MUST TRANSPOSE.** MATLAB is column-major, TIFF stores raster rows;
writing `P` instead of `P.'` gives a file with perfectly valid metadata and every page
silently transposed — invisible on a square tile. Caught only by verifying on
non-square 64x48 pages with different content per page.

**Measured, 260804, one tile (2 ch x 54 z x 50 t = 5400 pages):** original ~190 s ->
53 s with deint -> **35 s avgz-only**, output byte-identical to the production avgz
(max|new-old| = 0). `writeDeint = false` is the default now: the deint file is a
regenerable cache (the raw TIFF is still on disk) and NOTHING downstream reads it —
the stitcher and both GUIs open `ch<N>/avgz` only. 57 MB/tile vs 2.5 GB.
Batch: 36 tiles ~21 min / 2 GB; 81 tiles ~47 min / 4.6 GB.

Untested: the volume-interleaved branch (falls back to the indexed path with a warning),
single-channel datasets, and the other 35 tiles.
Headroom left: avgz-only reads at 77 MB/s against the drive's 129 — the read-side twin
of the same per-page TIFF overhead, worth ~1.5x if ever wanted.

## N-channel generalisation of the deinterleaver (2026-09-11)
`ch1_id`/`ch2_id` are GONE. One setting now: `chUse = []` = every channel the file saved.
`SI.hChannels.channelSave` already carries everything — its LENGTH is the page cycle and its
ORDER is each channel's position in it — so nothing is hard-coded to two. Core rewritten over
`numel(cIdx)`: one accumulator plane per channel (`acc(:,:,slot)`), one writer per channel,
and a `slotOf` map so unwanted channels' pages are stepped over unread. The output-folder
filter matches any `ch<N>` dir by regex, since the count is unknown until each header is read.

Tested with SYNTHETIC ScanImage stacks (E: was unmounted) where every page is a constant
`ch*100 + z*10 + t`, so a mis-deinterleaved page shows up as a wrong number: 3-channel auto,
1-channel (`channelSave = [2]`), and 3-saved-with-`chUse = [1 3]` — all exact, non-square
8x6 frames un-transposed.

**Latent bug fixed on the way:** `indexed_path` (the volume-interleaved branch) was called
with FILENAMES in parameters that expected open `Tiff` objects, and one argument short of its
own signature. It would have errored on the first volume-interleaved tile with `writeDeint`
on. Dormant because `writeDeint` defaults false. It opens its own handles now.

## to-send copies are COMMENT-STRIPPED (2026-09-11)
`D:\RZ_ScanImage_script\to-send\{auto-tile,stitch}\` holds comment-free copies for sharing.
Stripped with a MATLAB-aware tokenizer (scratchpad `strip_comments.py`) — a naive `%` strip
destroys `sprintf('%s_col%02d_row%02d_x%.0f_y%.0f')`, the line that builds every tile
filename the stitcher parses. Verified by `mtree` node counts identical to the commented
originals (2653, 1616, 1481, 8564, 8234), i.e. only comments went.
**Edit the COMMENTED original, then re-strip** — never edit the stripped copy.
