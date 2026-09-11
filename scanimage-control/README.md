# scanimage-control

Mirror of `D:\RZ_ScanImage_script` on the imaging computer — the rig-side scripts,
version controlled. The ScanImage 2018b install itself is NOT here (third-party,
~850 .m files); it lives beside these on the rig at
`D:\RZ_ScanImage_script\SI2018bR1_2018-12-19_4a9264c4fc\`.

## Layout

| folder | what it is |
|---|---|
| `auto-tile/` | **the live pair.** `auto_acq_260910.m` (snake-tile map acquisition) + `acq_pause_queue.m` (live pause/brake/abort window). Run these on the rig. |
| `stitch/` | stitching scripts as deployed on the rig |
| `to-send/` | comment-stripped copies for sharing. Generated, never edited by hand. |
| `commented_originals/` | backup of the commented `auto-tile` pair |
| `archived/` | superseded: the original `auto_acq_n_shutterOFF.m`, earlier tiling attempts, laser-control vendor code |

## Which copy is authoritative

The stitching scripts also exist at `../260910-map-registration/`, which is where
they are **edited**. `stitch/` and `to-send/stitch/` are deployment copies.

> Edit the commented original, then re-copy and re-strip. Never edit a stripped
> copy or a deployment copy — `scanimage-control/stitch/tile_deinterleave_avgz_260910.m`
> had already silently fallen a version behind (still the old two-channel
> `ch1_id`/`ch2_id` form) before this was written down.

MATLAB's current folder outranks `addpath`, so a stale duplicate can run instead
of the one you meant. That is the failure mode these copies invite.

## Status

`auto_acq_260910.m` and `acq_pause_queue.m` have been tested off-rig against mocks
(snake order, filenames, pause/brake/abort, retry) but **no hardware call has ever
been executed**. Shakedown before trusting a real map:
`numCols = 2; numRows = 1; doLaser = false;`
