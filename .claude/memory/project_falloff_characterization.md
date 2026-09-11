---
name: project_falloff_characterization
description: "Tissue optical penetration from vessel depth series: C:\\fall-off\\ data root, falloff_260804.m, laser_power_calibration.m (date-aware LUT), and the finding that attenuation length is NOT constant with depth"
metadata: 
  node_type: memory
  type: project
  originSessionId: f3c40add-59ae-4001-987d-927c5eb434ad
  modified: 2026-08-05T04:53:32.361Z
---

Systematic characterization of tissue optical penetration from **vessel depth series** — one ROI, imaged repeatedly down the same column, either as one power-ramped stack or several overlapping fixed-power stacks. Started 2026-08-04.

**Data root is `C:\fall-off\<session>\`** (e.g. `C:\fall-off\260728_vglut2_vessel`), not the session data drive. Runqi copies the vessel stacks there. One ROI per series, grouped by the `roiN` filename prefix.

## Scripts (repo root)

- **`falloff_260804.m`** — the falloff plot. Standalone, `falloff_260804(dataDir, roiTag)`. Caches per-slice metrics to `falloff_<roi>_metrics.mat` (delete to re-read; the full read is ~2.5 min for 7.6 GB). Four panels: raw tail → power-normalised → stitched + global exp fit → **local ℓ(z)**.
- **`laser_power_calibration.m`** — the ONLY place the %→mW tables live. `laser_power_calibration(pct, acqDate)`, or pass a table name (`'pre_260723'` / `'post_260723'`) to force one explicitly.
- **`laser_power_calibration_compare_260804.m`** — analysis of the two LUTs on their own (no imaging data).

## The power LUT changed on 2026-07-23

Two measured tables, and picking the wrong one is silent and **squared** downstream. See [[reference_laser_power_lut]] for the numbers. `laser_power_vs_z_260530.m` hard-codes the old one.

## Metric decisions (differ from `260114_shi_temp/F_attenuation.m`)

- **tail = brightest `TAIL_FRAC` = 0.4 % of pixels**, a fraction not the old hard-coded 1000 px, so a zoom or FOV-size change stays comparable.
- **Mode was dropped entirely** — it is broken, see [[project_260728_vglut2_vessel_depth]].
- **Offset subtracted before dividing by power.** Additive detector floor, measured **120.3 ADU** on 260728 (estimated from the deepest slices of the lowest-power round; a shutter-closed acquisition would make it a measurement — `DARK_ADU` param).
- **Saturated slices excluded from every fit**, counted on RAW frames, drawn dotted.
- **Per-round stitch scales printed**, measured from the depth overlap. They are the honesty check on the power correction.

## Result on 260728 roi1 (old LUT, by explicit choice)

Stitch scales **1.00 / 1.20 / 1.60 / 2.13** — each successive round needs a 20–33 % boost, i.e. gives *less* than mW² predicts. Global fit ℓ = **39.8 µm**, R² = 0.937 over 30–275 µm.

**But ℓ is not constant — it rises monotonically from ~22 µm at 75 µm depth to ~95 µm at 250 µm.** The single exponential describes neither end. This is *not* a stitching artifact: a per-round scale is a vertical offset in log space, so it cannot tilt the local slope *within* a round, and the rise continues right through the 48 % round which alone spans 150–275 µm.

**Most likely cause is the metric, not the tissue.** The brightest 0.4 % of pixels are in-focus vessel wall near the surface but become out-of-focus scattered haze at depth, and haze decays more slowly than ballistic signal. So the tail silently changes what it measures with depth. Testing that is the open question — it needs a metric that stays in-focus (e.g. contrast or a vessel-mask-restricted statistic) rather than "brightest N".
