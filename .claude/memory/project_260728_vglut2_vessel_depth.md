---
name: project_260728_vglut2_vessel_depth
description: "260728_vglut2_soma-g8s: behaviour/phys pairing (fileparts folder-name bug), and the vessel depth series where filename z is DEPTH but ScanImage zs is MOTOR position"
metadata: 
  node_type: memory
  type: project
  originSessionId: b2f03943-7192-4cdf-a370-82cf33218d5f
  modified: 2026-07-29T18:37:08.603Z
---

Dataset acquired 2026-07-28, both halves on **D:** — `D:\260728_vglut2_soma-g8s\{phys,vessel}` and `D:\260728_vglut2_soma-g8s_falling\cam1`. 19 phys FOVs, 19 behaviour AVIs, paired 2026-07-29 (all dt < 1.5 s, one AVI trimmed 3009->3007).

**`pair_behavior_to_phys.m` has a latent folder-name bug.** Step 1 tests `strcmp(tifBase, parentName)` with `[~, parentName] = fileparts(allTifs(i).folder)`. On a FOLDER path `fileparts` reads the last dot as an extension, so `fov1_1.7x_y1650x1320_z30_..._00001` returns `fov1_1` and the FOV is silently dropped — no warning, it just reports fewer FOVs. Any zoom or power with a decimal point hits it (1.7x, 2.1x, 2.5x, 15.5lp, 17.5lp): **11 of 19 FOVs on 260728**. Fix is to take the last path component without `fileparts`:
```matlab
parts = split(string(p), filesep);  parts(strlength(parts)==0) = [];  nm = char(parts(end));
```
Same bug at the display line `[~, fovName] = fileparts(fov(fi).folder)`. Also in that script: `writelines` (the trim path) **does not exist in R2021b** — added R2022a — so trimming crashes; and the ffmpeg auto-detect finds the broken conda binary first (see [[reference_ffmpeg_broken_use_imageio]]), so FFV1 trims fail silently before that. Working binary: `C:\Users\Admin\.conda\envs\dlc310\lib\site-packages\imageio_ffmpeg\binaries\ffmpeg-win-x86_64-v7.1.exe` (the ffmpeg search now globs `envs\*\lib\site-packages\imageio_ffmpeg\binaries\ffmpeg*.exe` first, so the version can change). **All four fixed in the repo script 2026-07-29 and re-verified: 19/19 FOVs found, was 8.**

**`_meta.mat framesPerSlice` is the CONFIGURED frame count, not what was written — and Runqi stops acquisitions by hand routinely, so this is a recurring case, not a one-off.** Measured across the 19 FOVs of 260728: 18 match the raw TIF page count exactly, and the one that does not is the manually-stopped recording — `fov1_2.1x_..._z30_3000f_13lp_00001` has `framesPerSlice = 6000` in `_meta.mat` but only **3007** pages on disk (folder named "3000f"; ScanImage was set to 6000 and stopped early).

**Rule Runqi confirmed 2026-07-29: when the header and the file disagree, trust the actual data.** `pair_behavior_to_phys.m` used to prefer `_meta.mat`, which on a stopped run compares the AVI against 6000 instead of 3007 and skips the trim in exactly the case that needs it. Now fixed: a new `count_tif_pages()` local walks the TIFF directory chain and the raw file is authoritative; `_meta.mat` is only a fallback when the raw TIF is unreadable, and a disagreement prints `stopped early -- N frames on disk vs M configured; using disk`. The walk costs **1.0 s** on a 3 GB / 6000-page BigTIFF versus 4.3 s for `imfinfo`, because it never builds per-page metadata structs — so ~20 s for a 19-FOV session, and there is no longer any reason to trust the header.

This also removes the earlier "run pairing before dffQC" ordering constraint — the script is now correct in either order. Verified end-to-end on 260728: AVI and timestamps.csv both trimmed 3009 → 3007 to match the TIF.

## Layout independence (2026-07-29)

Runqi's behaviour layout varies by what he did that day — MATLAB GUI or Pylon Viewer, `_falling/cam1/<run>/*.avi` or a flat pile of AVIs, sometimes a timesheet .xlsx alongside — and he wants the script to just cope. Surveying every AVI on D: and E: found **exactly two filename patterns**, and they can be MIXED IN ONE SESSION (260721_Sert has one Pylon file among GUI ones), so format must be decided per file:
- `<cam>_YYYYMMDD_HHMMSS_runNNN.avi` — MATLAB GUI, has a sibling `timestamps.csv`
- `<model>__<serial>__YYYYMMDD_HHMMSSmmm.avi` — Pylon Viewer, 3-digit ms, no run suffix

`pair_behavior_to_phys.m` now handles this:
- **`behavDir = ""` auto-detects** — looks for a sibling of the session folder whose name starts with the session name AND actually contains AVIs; falls back to the session folder itself. Nesting inside is irrelevant (the search is recursive) and non-video files are ignored because only `*.avi` is globbed.
- **Per-camera matching.** Each camera records the SAME FOV, so a single shared used-mask made cam1 claim a FOV and left camera 2's video unmatched. Now one used-mask per camera, keyed on the label captured before the timestamp — which works for any `--cam2-name` (`basler_dual_acq.py` takes it freely, so it is NOT necessarily "cam2"). No dual-cam data exists yet, but the backend supports it.
- **AVIs already inside physDir are ignored** as previous-run copies, so behavDir can safely point at the session root — or at physDir, which is the 260721_Sert layout where there is no separate behaviour folder.
- If every video found is already paired it prints `Nothing to do` and returns instead of asserting.

Verified on both layouts: 260728 (sibling `_falling`, auto-detected, 19/19, camera `cam1`) and 260721_Sert (no sibling, 27 FOVs across nested `baseline/`,`roomAir/`,`reoxygenated/`, all 27 recognised as already paired).

## Vessel depth series — two conventions that disagree

`vessel/` holds a depth series acquired as **overlapping stacks at different fixed laser powers**, so no single stack spans the full depth at one power. roi1 has 4 rounds (11/21/35/48 % lp); roi2 has only the first (unfinished).

**Filename z is DEPTH, positive downward, with the sign dropped. ScanImage `zs` is MOTOR position, negative downward.** They are not the same number:

| file | filename z | `SI.hStackManager.zs` |
|---|---|---|
| `roi1_z-10-100_lp11` | −10..100 | **+10 → −100** |
| `roi1_z40-160_lp21` | 40..160 | **−40 → −160** |
| `roi1_z90-210_lp35` | 90..210 | **−90 → −210** |
| `roi1_z140-265_lp48` | 140..265 | **−140 → −265** |

Always use `zs`. All roi1 stacks: 512×512, zoom 2, framesPerSlice 30, 1 µm step, 1 channel.

**Power is CONSTANT within each stack** — `pzAdjust = false` and `stackStartPower == stackEndPower`. But `SI.hBeams.lengthConstants = -374.237` is still present and stale. `260114_shi_temp/F_attenuation.m` applies `P0*exp((z-z0)/LC)` **unconditionally**, which on these files invents a ~35% power rise by the bottom of a flat-power stack. Any script reusing that code must branch on `pzAdjust` first.

## What the roi1 depth series actually shows (2026-07-29)

**The top of every round is SATURATED and the tail is invalid there** — 37/111, 45/121, 39/121, 16/126 slices for 11/21/35/48 % lp, clipping at the int16 ceiling 32767. It does not look like a failure: the tail goes *flat* at ~3.2e4, which reads as a plateau. Saturation must be counted on the RAW frames, not the 30-frame average (a pixel clipped in 10 of 30 frames averages well below the ceiling and looks clean).

**ScanImage % is a Pockels setpoint, not power, and it matters more than the saturation.** Correcting by (%)² leaves the three overlapping pairs disagreeing by 2.50 / 1.67 / 1.05 — a 2.4-fold spread, so no single power law fits. Re-doing the same ratios with the measured %→mW table from `laser_power_vs_z_260530.m` collapses them to **0.83 / 0.75 / 0.75** — near-constant. The mW calibration is clearly the right denominator in shape; the residual ~0.78 per step is a separate uniform effect (bleaching between rounds and/or sub-quadratic response at these powers — the implied exponent still drifts 1.84→1.68→1.41, so it is not a clean power law either). Verify the calibration is current for this rig before quoting absolutes.

`vessel_mode_tail_vs_z_260728.m` (repo root, written 2026-07-29) computes mode (first histogram peak from the left) and tail (mean of brightest 1000 px above the 0.5th pct) per z, one colour per power round, on a shared reversed-x `zs` axis, plus P²-corrected panels and a printed overlap-agreement check between adjacent rounds. Metric definitions copied from `F_attenuation.m`; see [[project_dffQC_pipeline]] for the phys side.

## The mode metric is BROKEN — measured 2026-07-29, do not reuse it

The earlier claim "mode is the robust metric, the blue 11 % spike is saturation" is **wrong**. `vessel_hist_raw_260729.m` (repo root) dumps the full pixel-value distribution of every slice on ONE fixed 8-ADU axis (cached to `vessel/vessel_hist_raw_roi1.mat`, delete to recompute), and it shows three separate things:

- **The distribution has a hard left edge at ~150 ADU in all four rounds** — the digitiser offset. It does not move with depth.
- **`findpeaks` cannot return an edge maximum.** In the 11 % round the frame-averaged histogram *decreases monotonically* from that edge — there is no interior background peak at all — so "leftmost peak" skips it and returns the next interior bump, the bright-tissue peak at ~2.7e4. Measured saturation at those slices is only **1–5 %**, so saturation is not the cause. 20 of 111 slices, ratio up to **172×**.
- **Plain argmax fails in the opposite regime.** At the top of the 21 / 35 % rounds up to **9.6 %** of pixels sit at the ADC ceiling, so the most-common bin *is* 32764. Neither estimator is safe on its own; both need the ceiling bin excluded and a monotone-decreasing histogram handled.

**The per-slice adaptive binning is the mechanism.** The script re-derives edges as `prctile(I,[0.5 99.5])` / 512 bins, so bin width runs **~60 ADU at the surface down to ~2.8 ADU at depth**. At 60 ADU the 30-frame-averaged noise peak is 1–2 bins wide and lands at the left edge of the window.

**Worse, "mode" is not the same quantity across rounds.** At 11 % it is pinned at 155–175 ADU, dead flat with depth — it is the offset, there is no background above it. At 21/35/48 % it is a genuine interior peak (out-of-focus scattered background) that decays smoothly. So mode/P² can never line up across all four rounds.

**The median is the estimator that behaves** — smooth, monotone, no discontinuity in any round (summary panel). Whatever replaces the mode: subtract the ~150 ADU offset *before* any P² division (it is additive and power-independent, so raw-ADU/P² is wrong at depth where signal approaches the offset).
