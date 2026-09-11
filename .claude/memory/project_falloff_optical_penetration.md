---
name: project_falloff_optical_penetration
description: "Optical-penetration project: goal (wildtype vs shiverer 2P depth limit), C:\\fall-off data layout, acquisition convention, and the standard falloff + projection method in falloff-analysis-260805/"
metadata: 
  node_type: memory
  type: project
  originSessionId: f3c40add-59ae-4001-987d-927c5eb434ad
  modified: 2026-08-17T07:21:55.503Z
---

# Goal

Systematically characterize **how deep two-photon imaging reaches, and why**, comparing **wildtype vs shiverer** mice (shiverer = myelin-deficient, so the prediction is less scattering / longer penetration). Other treatment groups may be added later — keep the pipeline group-agnostic; group is a property of the SESSION FOLDER, not hard-coded.

Measured on **vessel-labelled stacks** (dye-filled vasculature), one ROI = one column imaged repeatedly down.

Two standard deliverables per ROI:
1. **falloff summary** — 6-panel intensity-vs-depth figure with the attenuation length
2. **projections summary** — depth-corrected max projections (x-y bins + sagittal + coronal)

# Data layout

Root is **`C:\fall-off\<session>\`** (NOT the acquisition drive; Runqi copies stacks here).

**TWO acquisition styles exist and the scripts handle both:**

| style | sessions | shape |
|---|---|---|
| **A: overlapping fixed-power stacks** | `260728_vglut2_vessel`, `260804_shiver_dbh_vessel`, `260806_sst_vessel` | one ROI = 4 stacks (260806: 3), each deeper + higher power. 512×512, zoom 2, 30 f/slice, **1 µm step** |
| **B: one full ramped stack per FOV** | `260114_shiver_vessel`, `251104_wt_fitc` | `pzAdjust = true`, power ramps exponentially (10→40 %, 7→26 %) via `lengthConstants`. zoom 1.2, **10 µm step**, 50–100 f/slice |

FOV grouping is by a tag passed to the script. Prefix works for `roi1`/`tile15`; when the distinguishing part is in the MIDDLE of the name (251104: `..._col00_row02_...`) pass a **glob** `'*col00_row02*'` — `match_stacks()` handles both and sanitises the tag for output filenames.

## Two-channel sessions (from 260806)

`260806_sst_vessel` is the first with `channelSave = [1;3]`: **channel 1 = the SST cell label, channel 3 = the vessel dye**. Pages are **interleaved per frame** — `page = ((slice-1)*fpsl + frame-1)*nCh + chIdx` — so reading the wrong channel does not error, it silently measures the other label. All three readers (`falloff_260804`, `falloff_projections_260805`, `spatial_fourier_260806`) now take a trailing `channel` arg = the **SI channel number**, and the `SESSIONS` table in `falloff_run_all_260805` has a 7th column for it (`[]` on single-channel sessions, which behave exactly as before). Caches store `CHANNEL` and re-read on a mismatch — a cache from another channel is a different label, not a variant.

## The z convention — the thing that keeps biting

There is **no single rule**; `DEPTH_MODE` selects it, and getting it wrong costs a full re-read.

- **`'filename'`** (style A) — depth = filename z token + uniform `stackZStepSize`. The filename z IS depth below pia.
- **`'motor'`** (style B) — depth = `Z_SURFACE − (stackZStartPos + k·stackZStepSize)`. Here the filename z is MOTOR z (`z50_-700`, step `-10`), not depth.

Never use per-slice `SI.hStackManager.zs`: ScanImage z is inverted, the returned steps are non-uniform (0.995 µm), and **the motor zero moves between stacks** — on 260804 the implied surface was 75/100/100/125 µm across four stacks of one ROI.

**Style C (added 2026-08-17): one ramped stack of a SOMA label.** `C:\260730_shiver_vgat_g8m` (shiverer, vgat-G8m) — a single 12.4 GB stack, `roi1_stack_00001.tif`, no z token in the name, 471 slices x 50 frames, 1 µm step, motor z **+20 → −450**, zoom 1.2 (1.4815 µm/px), one saved channel, `pzAdjust` 10 → 50 % with `lengthConstants = −292.027`, `zPowerReference = +20`. Script: **`falloff_gcamp_260730.m`** — standalone, reads the stack ONCE and emits both the falloff and the projections figure. Walk pages with `nextDirectory`, not `setDirectory(t,k)`: the latter is quadratic over 23 550 IFDs. Read = 212 s.

**The raw-(top−bot)-peak pia rule BREAKS on a ramped stack**: power climbs 5x with depth, so raw signal rises by construction and the peak landed at motor **−256**, 250 µm into the tissue. Search on **(top−bot)/mW²** instead. On a soma label even that peak is biased deep (sparse L1 → peak at motor −13), so the default marker is the **half-rise** of the normalised curve, the blurred tissue/water edge where the beam enters: motor **+7.3** here, against the +0 the range implies. `SURF_MODE = 'halfrise' | 'peak'` selects; both are always printed. The choice barely touches ℓ (a shift changes the intercept, not the slope): 78.0 µm at half-rise.

**`Z_SURFACE` (the pia) must be found from the data, not assumed.** `zPowerReference` / `stackZStartPos` is where the ramp is referenced, often set ABOVE the pia. Find it as the depth where **(top−bot) peaks** in the shallowest stack:
- 260728, 260114 → pia at motor **0** (260114 stack starts 50 µm above it)
- 251104 → pia at motor **−20** (peak at −20 in 3 of 4 FOVs, −30 in the fourth)

Stacks routinely start above the pia (negative depth). Those points are drawn faded and excluded from fits.

**Never set a fixed `Z_TOP` truncation.** `Z_TOP = 10`, written for 260728, silently deleted 115 of 151 slices on 260804. It is `Inf` now.

## Overlap ownership

Consecutive stacks overlap. Each overlap is split at its **MIDPOINT** — shallower stack keeps the top half, deeper takes the bottom half — so every depth belongs to exactly one stack. Same rule in both scripts. Per-stack fits still use all of their own slices; the split only governs the POOLED fit. Before this, overlap depths entered the pooled fit twice with values ~25 % apart.

# Standard method

Scripts: **`falloff-analysis-260805/`** (moved there 2026-08-05; each addpaths its own dir + the repo root).

| file | role |
|---|---|
| **`falloff_run_all_260805.m`** | **THE BATCH SCRIPT — start here.** `falloff_run_all_260805` runs everything; `(root,'falloff'\|'proj'\|'scatter')` runs one stage. Holds the SESSIONS table (group, LUT, zSurface, depthMode, FOV tags). Adding a session = adding a row. |
| `falloff_260804.m` | ONE FOV, falloff figure. `(dataDir, roiTag, zSurface, powerLut, fitDepth, depthMode)` |
| `falloff_projections_260805.m` | ONE FOV, projections figure. `(dataDir, roiTag, powerLut, ell, zSurface, depthMode, mount)` |
| `falloff_summary_scatter_260805.m` | cross-session ℓ scatter; group table at the top. Since 2026-08-18 it reads BOTH writers (`Lp`/`r2p` from `falloff_260804`, `ell`/`r2` from `falloff_gcamp_260730`) and accepts an ABSOLUTE session path, so 260730 (which lives outside `C:all-off`) is a normal row. Axis label is now `λ`, not an italic ell that read as I/1. |
| `laser_power_calibration.m` | %→mW, date-aware. See [[reference_laser_power_lut]] |
| `laser_power_calibration_compare_260804.m` | LUT analysis only, no imaging data |
| `vessel_hist_raw_260729.m` | diagnostic: full pixel distribution per z per stack |
| `vessel_mode_tail_vs_z_260728.m` | SUPERSEDED first attempt; mode metric is broken |

Both main scripts cache to the data folder (`falloff_<tag>_metrics.mat`, `falloff_proj_<tag>_vol.mat`) — delete to force a re-read (~2 min per 8 GB). The metrics cache also holds a **per-slice fixed-bin histogram**, so any distribution depth can be replotted for free.

## Metrics

- **top 5 %** = mean of the brightest 5 % of pixels of the frame-averaged slice; **bot 5 %** = mean of the dimmest 5 %. Fractions, not pixel counts. Widening 1 %→5 % barely moved ℓ but noticeably tightened R² (0.91→0.99 on the worst stacks).
- **top − bot** removes the additive floor.
- Normalise by **mW²**, mW from the date-aware LUT — never by percent.
- **ℓ** = decay constant of (top−bot)/mW² vs depth, per stack and pooled, log-linear OLS. **Literature EAL = 2ℓ.**
- **FIT_DEPTH is group-dependent**: **10–150 µm for wildtype, 10–500 µm for shiverer** — set by Runqi, passed per call. WT was 10–200 until 2026-08-06; shortening it to 150 dropped every WT ℓ by ~4 µm and **collapsed the WT spread from ±1.5 to ±0.4 µm**, because past ~150 µm the WT signal is into the floor and each FOV runs out at a different depth. It also moves the projection gain clamp to 150 µm, so the deep half of a WT side view renders darker than it used to.
- **Attenuation = F(0)/F(z)** with ONE F(0), literally the value at depth 0, shared by every stack. Raw on a LINEAR axis (matches published panels; segments reset at each power change); power-normalised on log with the pooled fit.

## Falloff figure layout (4 rows, settled 2026-08-05)

1. **F distribution** at 50, 150, 250 … µm (or 25, 75, 125 … if the stack never reaches 250). Lines, not dots; empty bins are NaN so the line breaks. x = −50…32767 ADU.
2. raw top 5 % | raw bot 5 %
3. top−bot | (top−bot)/mW² + per-stack and pooled fits
4. attenuation raw (linear y) | attenuation power-normalised (log y)

Rows 2–4 x-limits fixed **−50…600 µm with ticks at 0:100:600** (the −50 margin only keeps the faded above-pia points visible). Two **jet** colourbars on FIXED ranges so any two sessions can be laid side by side: **depth 0–600 µm** (reversed, 0 at top) and **power 10–350 mW**. The power legend keeps the discrete mW values; a ramped stack is coloured by its median mW and labelled with the range.

**Everything is on a fixed 0–600 µm axis at 100 µm ticks — falloff and projections both.** With 1×1×1 µm voxels x, y and z are then the same scale in every figure, which is what makes sessions comparable at a glance. A shallower session gets empty "no data" panels; a deeper one is cropped at 600.

## Projections figure

1 µm **isotropic**: lateral (1.7778/zoom) resampled **bicubic**; if the axial step is coarser (10 µm on style B) interpolate along z with **pchip**. Reslicing to y-z / z-x is then a pure permute. **Grayscale, not jet.**
**MAX** projections: x-y per 100 µm bin (six FIXED bins, 0–600), then **x-z coronal view** and **y-z sagittal view, in that order**. Side views carry a real-micron axis, depth 0–600 at 100 µm ticks, and are drawn with **`daspect [1 1 1]`** — the voxels are isotropic so the display must be too, or the vessels are sheared. Display **0.5–99.5 %** per panel.

**WHICH SIDE VIEW IS WHICH DEPENDS ON THE MOUNT.** Image row = stage y, image column = stage x always (scan mirrors are fixed to the rig); the mount only sets what stage x/y MEAN. Coronal = ML vs depth = max along the AP axis; sagittal = AP vs depth = max along ML. So the two views **swap** between old and rotated sessions:

| mount | sessions | AP axis | ML axis |
|---|---|---|---|
| `old` | 251104, 260114 | stage y (rows) | stage x (cols) |
| `rot90` | 260728, 260804 | stage x (cols) | stage y (rows) |

`MOUNT = 'auto'` derives this from the **YYMMDD prefix of the session folder** against the 2026-07-21 boundary and prints what it chose. The rule is temporal — never inferred from where the coordinates land. See [[project_rotated_mount_sessions]].

**Direction, not just name.** On a rot90 session the ML axis is the image ROW, and the image row runs along NEGATIVE stage y ([[reference_stage_axis_convention]]), so the coronal view goes **lateral → medial** left to right. Runqi caught this on 2026-08-06: it was labelled "medio-lateral", i.e. backwards. The rot90 coronal xlabel now reads `lateral \rightarrow medial`. The **old** mount was left as `x, medio-lateral` — its ML direction has not been measured, so it is not to be relabelled on the assumption that it mirrors rot90. The AP label is right on both. Still inconsistent and not worth touching unless asked: the panel titles ("x-z coronal") and the AP xlabel keep stage-letter prefixes that are swapped under rot90, where ML is stage y and AP is stage x.

**260730 axis directions, stated by Runqi (2026-08-17), not derived:** in the image **+y (row) = more LATERAL** and **+x (col) = CAUDAL, −x = rostral**. So rows = ML (medial → lateral) and cols = AP (rostral → caudal); coronal = max along cols, sagittal = max along rows. Note the ML direction is the OPPOSITE of the `rot90` label ("lateral → medial") derived for 260728/260804 — do not carry either direction across sessions without asking.

Three things the correction MUST do (each was a real bug):
1. **Per-slice power.** A ramped stack normalised by `stackStartPower` alone leaves deep planes ~100–150× too bright — this is what made the superficial layer look dark and the deep end white.
2. **Subtract the per-plane floor BEFORE `exp(z/ℓ)`**, because ℓ was measured on (top−bot). The floor also grows with depth, so amplifying it brightens the deep end for a non-vessel reason.
3. **Clamp the gain beyond `FIT_DEPTH(2)`.** Past the fitted range the signal stops decaying (sometimes rises), so continuing the exponential just explodes background.

The cached volume is RAW ADU + per-plane setpoint, so a re-measured LUT costs a replot, not a re-read.

# Results (2026-08-06, WT fit shortened to 10–150)

| session | group | n | ℓ (µm) | fit |
|---|---|---|---|---|
| 251104_wt_fitc | wildtype | 2 | 24.5, 24.7 | 10–150 |
| 260728_vglut2_vessel | wildtype | 1 | 24.6 | 10–150 |
| 260806_sst_vessel | wildtype | 1 | 25.5 | 10–150 |
| 260114_shiver_vessel | shiverer | 3 | 78.4, 80.7, 86.7 | 10–500 |
| 260804_shiver_dbh_vessel | shiverer | 2 | 61.4, 63.0 | 10–500 |
| 260730_shiver_vgat_g8m (GCaMP soma, not vessel) | shiverer | 1 | 78.0 (R² 0.996) | 10–500 |

**wildtype ℓ = 24.8 ± 0.4 µm (n=4, EAL 50)** vs **shiverer ℓ = 74.7 ± 10.2 µm (n=6, EAL 149)** — as of 2026-08-18 the shiverer pool INCLUDES 260730 (the soma label); it was 74.0 ± 11.2 at n=5 without it — ~3× longer in shiverer, the predicted direction for less myelin. Figure: `C:\fall-off\falloff_summary_scatter.png`, made by `falloff_summary_scatter_260805.m` (group table at the top of that script). (At the old 10–200 WT window these were 29.5/31.3, 28.0, 28.3 → 29.3 ± 1.5.)

260806 was entered as **wildtype** (SST-cre is a normal-myelin line; the shiverer folders say so in the name) — confirm with Runqi if that is ever in doubt, since group also sets the fit window. Its 15 % stack clips at the ADC ceiling over roughly −10 → +25 µm (5 % of pixels at the shallowest), but ℓ is insensitive: at the old 10–200 window, fitting from 10/20/25/30/40 µm gave 28.3/28.4/28.7/28.9/29.5 µm. No exclusion applied, per [[feedback_falloff_analysis_style]].

Caveats: n is small; the two shiverer sessions do NOT overlap each other (78–87 vs 61–63) and differ in acquisition style, LUT and fit window; **ℓ still climbs with depth/power within a session** in style-A data, which no metric variant removed.

Deleted 2026-08-05 as useless edge FOVs (Runqi has copies): 251104 col00 + col05, and 260114 tiles 15/16/21/22/28 — raw TIFFs included.

# Hard-won preferences (see [[feedback_falloff_analysis_style]])

Do **not** exclude or fade slices for saturation, do **not** stitch/rescale stacks to make them agree, do **not** add analyses that were not asked for.
