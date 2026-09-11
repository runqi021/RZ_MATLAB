# 2026-07-29 — Sert + vglut2 into the archive, population histogram, 8-panel per cell

Continues `PROGRESS_REPORT.md`. Four things were asked for; three are complete, one is
blocked on a step only you can do.

---

## 1. Archive — DONE (copy, sources untouched)

`archive_into_summary_260729.m`. Two modes: `'template'` writes an editable
`site_mapping_260729.csv`, `'copy'` does exactly what that csv says. **13 recordings
copied**, 18 correctly skipped.

```
Ventral_surface_summary\
  Sert\0721\cell1 (3)  cell3 (2)  cell4 (1)  cell5 (1)  cell6 (1)  cell7 (2)
  Vglut2\0728\cell1 (1)  cell2 (2)
```

Site convention read off the existing archive: in `Vglut2\0224` the pFN recordings are
filed `cell1..cell5`, one imaging site each, and `IO` is its own label. So the new pFN
data goes in as `cellN`, not a `pFN` folder.

**Two judgement calls you should check:**

- **`ro7` and `roi7` both went to `cell7`.** Their names differ only by the `i`, but
  the coordinates have opposite x sign — `x1000y-1000` vs `x-1000-1000`. If those are
  two different sites, split them; the fix is one edit to the csv and a re-run.
- **`ro8_2.3x_...` was skipped** — no `ca_spike_data.mat`, so no active cells. Run the
  detector on it and re-run the template if you want it in.

Skipped and why: MAO ×3 (part of IO, not pFN, as you said), 12 vglut2 recordings with
no spike detection, 2 with 0 active ROIs, and the excluded first Sert FOV (§2).

To undo: delete `Ventral_surface_summary\Sert\0721` and `Vglut2\0728`. Nothing was
moved or overwritten — the copy skips any folder that already exists.

## 2. The first Sert FOV is excluded everywhere — CONFIRMED

`ro1_2x_x1018y1130_18lp_3000f_00001`. You said it used a different zero reference;
that is confirmed three independent ways:

- acquired **15:24:58**, ~20 min before the next recording (15:45)
- the **only** recording in the session with no Basler `.avi` — the camera was not
  running yet
- its `motorPosition` y (1018.3) uniquely matches its own folder-name x (1018), while
  the later `ro1_*` recordings at the nominally same site do not

Now in `coh_cfg_260727.cfg.excludeRecordings`, honoured by the population histogram,
the 8-panel figure and the archive template. It was already in the matcher's
`tossFolders`; this extends it to the analysis and the archive.

**Why this matters beyond tidiness:** the cartoon map applies ONE stage→canvas affine
to every animal, so the whole spatial map assumes every session zeroed on the same
anatomical landmark. A recording that did not cannot be placed on it.

## 3. Population histogram — DONE

`population_hist_260729.m` → `Ventral_surface_summary\analysis_260729\`.
Fixed **3 s window (±1.5 s)**, same for every dataset — deliberately not IBI-scaled,
because the median IBI is 2.15 s in Sert and 0.63 s in the older ventral recordings, so
a scaled window would put different amounts of physiology on the same axis.

**Active cell = `nnz(spike_train>0) > 5`, matching your existing ventral criterion**
(`spike_trigger_dFF.m:34`, `temporal_phase_perROI.m:56`). Deliberately looser than the
≥20 events the per-cell PETH test needs: pooling means a 6-event cell contributes 6
real events, whereas a per-cell shuffle test on 6 events has no power.

| genotype | active ROI | rec | events | IBI | observed peak vs 95% band |
|---|---|---|---|---|---|
| **Sert** | 57 | 10 | 2761 | 2.15 s | **190 vs 96 — far above** |
| **Vglut2** | 138 | 13 | 3157 | 0.63 s | **229 vs 192 — above, periodic** |
| ChAT | 22 | 7 | 1348 | 1.33 s | 136 vs 148 — inside |
| Sst | 4 | 3 | 222 | 1.03 s | 15 vs 18 — inside |
| Vgat | 4 | 3 | 201 | 2.50 s | 6.3 vs 8.1 — inside |
| ALL | 225 | 36 | 7689 | 1.58 s | 450 vs 433 — marginal |

Sert is the clearest: a single sharp peak ~0.3 s after onset (~0.15 s after peak),
well clear of the band. Vglut2 shows **repeating peaks at ±1 and ±2 IBI** — the fast
breath rate (0.63 s) fits several cycles in a 3 s window, so its breath-locked
rhythmicity is visible directly. Sst and Vgat have only 4 active ROIs each and are
not interpretable either way.

**Cosmetic caveat:** the bars and the band both carry a sawtooth ripple, from ~30 fps
frames in 50 ms bins. It is present identically in observed and shuffle so inference
is unaffected, but for a figure you would want exposure normalisation as the per-cell
scripts use.

## 4. 8-panel per-cell figure — DONE

`percell_8panel_260729.m` → `breath_time\percell_8panel\`, **48 active cells** on Sert
baseline.

```
row 1   onset dF/F avg | peak dF/F avg | onset event hist | peak event hist
row 2   onset dF/F map | peak dF/F map | power spectrum   | spike-triggered dF/F
```

All triggered panels share the same fixed 3 s window. Grey band = circular-shift
shuffle 95% throughout, so "is this real" is readable without a p-value.

Deliberate choices worth knowing:

- **The dF/F heatmaps are chronological and use every accepted breath.** Not sorted by
  event latency, not restricted to event-containing cycles — either of those
  manufactures a diagonal band out of noise, which is what made the older heatmaps
  misleading.
- **The power spectrum is per-recording then averaged** on a common 0.1–14 Hz grid,
  because fps differs between recordings. The breath rate is marked: a peak there is
  breath-locked, a peak elsewhere would be the cell's own rhythm.
- **The spike-triggered dF/F is partly circular** — the events were detected from that
  same dF/F — so read it as the event kernel / a detection sanity check, not evidence.

Example, cell 159: dF/F and events both clear the band, spectrum peaks exactly at the
0.39 Hz breath rate with nothing above 2 Hz, clean event kernel.

## 5. Cartoon map — DONE, Sert is on it

`Ventral_surface_ROI_on_cartoon_260729.m` → `Desktop\ROI_on_cartoon.{png,pdf}` +
`ROI_on_cartoon_data.mat`. **37 FOVs, 1263 ROIs, 32 sites drawn; Sert = 402 ROIs in
orange.** Three things had to be solved.

### 5a. The missing calibration was recovered, not re-clicked

`ROI_on_vessel_data.mat` — the stage→anatomy dictionary (`ax`/`ay` affine to
stitched-canvas px, `bx`/`by` = hand-clicked basilar branch = anatomical origin,
`theta_rad` = vessel-axis rotation, `um_per_px_stitch`) — was missing, and its producer
`Ventral_surface_ROI_on_vessel_260528.m` is interactive.

But a **prior successful run had saved all six variables into its own output**,
`Desktop\ROI_on_cartoon_data.mat`. So `recover_vessel_calib_260729.m` reads them back and
writes the file. This is recovery, not re-derivation: the new map lands in the *same*
anatomical frame as the old one, which is the whole point if the added Sert cells are to
be comparable. Recovered: `ax=[1.03797 1629.81]`, `ay=[-1.06443 3581.34]`,
`bx,by=1634.53,3633.19`, `theta_rad=0.083686` (4.79°), `0.8889 µm/px`. It refuses to
overwrite an existing file, so re-click only if the stitch itself changes.

### 5b. Two sessions were mounted rotated 90° — Sert/0721 AND Vglut2/0728

The map used to assume every animal was mounted the same way, i.e. stage x = lateral and
+stage y = rostral. For these two sessions that is wrong:

```
-stage x = anterior (rostral)      stage y = lateral (+y = left)
```

Encoded as `AXMAP_ROT90 = [0 1; -1 0]` sending `[stage x; stage y]` → `[lateral; rostral]`,
keyed on **session (scan dir + date), not genotype** — `Vglut2\0224` and `Vglut2\1124` are
*not* rotated, so a genotype-level flag would have wrongly moved the older pFN data.

Three consequences, each easy to miss:

- The **lateral sign does not matter** — everything is folded into the left hemisphere with
  `-abs()`. The real content is "rostral = −stage x".
- The **FOV footprint transposes.** Scan mirrors are fixed to the rig, so a FOV column is
  always stage x and a row always stage y; the mount enters only at stage→anatomy. Sert's
  `ro7` is 320×512 px: it spanned 569 µm of stage x, so it now spans 569 µm *rostrally*
  and 356 µm laterally. `rot90(...,-1)` carries texture and extent across together. Sert
  has genuinely non-square FOVs (320×512, 512×380), so this is not academic.
- **Independent corroboration, from your own filenames**: the two rotated sessions write
  the coordinate pair y-first (`y1900x1320`) and their `motorPosition` satisfies
  `motorX = -nameY`, `motorY = nameX`. Filenames are otherwise unreliable in these
  sessions (one Sert name has a flipped sign) — `motorPosition` is the only authority.

Before the fix Sert was two disconnected blobs, one at y≈+500 and one at y≈−1100. After,
it is a single contiguous rostro-caudal strip, 400–1200 µm lateral, y ≈ 0 → +1400.

### 5c. One patch per SITE, so a re-imaged site is not drawn twice

Several sites were re-imaged at a second Z, and `ro7`/`roi7` were recorded twice under
near-identical names. Each stamped its own 0.85-alpha patch at the same place, and piled
up they read as one dark smear. FOV centres within `mergeSiteUm = 100` (same session +
group) now collapse to one patch, drawn from the recording with the **most ROIs**.

**Every ROI is still in `R`** — two Z planes are two different sets of cells; only the
*rendering* was duplicated, so this is a rendering fix and culls no data.

Five merges, and one of them found a mistake:

| dropped | kept | apart |
|---|---|---|
| Sert `ro1_…z10` | `ro1_…z15` | 42 µm |
| Sert `roi7_1.8x` | `ro7_1.6x` | 10 µm |
| Vglut2/0728 `roi2_…z-20` | `roi2_…z-15` | 28 µm |
| Vglut2/1124 `roi3_-1000-300-40` | `…-60` | 0 µm |
| Sert `ro3_4x_x700y-1233` | `ro1_…z15` | **73 µm** |

- **`cell7` was filed correctly** — its two recordings are 10 µm apart, so despite the
  opposite-sign x in their *names* they are one site. That answers the open question from
  §1.
- **`ro3_4x_x700y-1233` is in the wrong archive folder.** By stage position it is a 4×
  zoom-in on the **cell1** site, while its supposed sibling `ro3_2.5x` is 312 µm away.
  Harmless for the analysis (cell identity comes from the morphology matcher, not the
  site folder) but worth correcting for tidiness: edit `site` to `cell1` in
  `site_mapping_260729.csv` and re-run the copy.

### 5d. The coherence foundation was re-run so Sert gets significance — DONE

`Ventral_surface_coherence_polar_svd_260729.m` (new copy of the `_260606` script) →
`Ventral_surface_summary\coherence_polar_svd_260729\`. Same analysis, same parameters,
same archive — only `scan_dirs`/`groups` extended to Sert, plus an orange tint and the
`coh_cfg_260727` exclusion list honoured. **The 260606 output is untouched**, so the two
stay comparable.

| run | ROIs | significant |
|---|---|---|
| 260606 | 297 | 10 — 9 Vglut2, 1 ChAT |
| **260729** | **368** | **14 — 9 Vglut2, 4 Sert, 1 ChAT** |

Every archive recording has both PC1 breath files (`breath_peak_pc1.mat`,
`breath_insp_start_pc1.mat`) except `ChAT\0522\7N`, which has none and is skipped as
before. `confC = 0.827` (α=0.001, TW=4, K=7) — unchanged, so the bar did not move.

**The 4 significant Sert cells all sit in the pre-peak half of inspiration**, phase +1.68
to +2.50 rad (peak = π ≈ 3.14), r = 0.83–0.93. That is the same phase range as the 9
Vglut2 pFN cells (+1.31 to +2.18), which is why both draw green on the map. All 4 are at
the `cell1` site.

The cartoon map's `dataMat` is repointed at the 260729 output. Also fixed on the way: this
script's `sig_rois.csv` split labels on `#`, which never occurs in a
`group/date/recName/rid` label, so its `fov` and `roi` columns both contained the whole
label. Now splits on `/` and adds a `date` column.

### 5e. The legend reports CELLS, not masks

A Cellpose mask is not a cell — the same neuron re-imaged at a second Z or zoom yields one
mask per recording. The map now resolves identity through the matcher's `cell_link.mat`,
keyed on the portable `(recording name, maskL label)` pair, registered per session in
`cell_link_sources`. Two statuses it carries are honoured:

- `tossed` (`cell_id` NaN) — rejected during your curation, **not** a cell. Dropped from
  the map entirely (13 masks across 6 Sert recordings).
- `ungrouped` — a real cell seen in only one recording. Counts as one cell.

Sessions with no matcher run fall back to one cell per mask, which is what
`cell_link_260727` itself does without curation, so the count degrades gracefully instead
of silently reverting to mask counts.

**Sert: 402 masks → 13 tossed → 389 → 296 cells.** Title now reads
*1157 cells (1250 masks) / 37 FOVs, 32 sites*.

| group | masks | cells | |
|---|---|---|---|
| IO | 679 | 679 | no matcher run |
| ChAT | 39 | 39 | no matcher run |
| Vglut2 | 118 | 118 | 14 masks matched (the 0728 recordings; identity grouping, as you asked) |
| Vgat | 13 | 13 | no matcher run |
| Sst | 12 | 12 | no matcher run |
| **Sert** | **389** | **296** | all matched + hand-curated |

Carrying that further down the funnel, for Sert: **296 cells → 48 active (>5 events) → 34
testable (≥20 events) → 18 breath-modulated**. The 34 is an independent check — it is the
same number the PETH reported.

### 5f. Two smaller things

- **Outlines survive patch merging.** The site dedup suppresses only the tinted blob; the
  white significant-ROI outlines still draw for every recording, each in its own extent.
  Without this, 3 of the 4 significant Sert ROIs would have lost their outline because
  their recording's patch was the merged-away one. Masks from different recordings cannot
  be overlaid on a single patch anyway — different zoom, pixel size and Z.
- The legend moved inside the axes (`southeast`); `eastoutside` was covering the
  colorbar's `preferred phase (rad)` label.

**One thing left as-is:** the non-significant grey dots on the mirrored right hemisphere
have never actually rendered — `scatter(..., 'MarkerEdgeColor','none')` without `'filled'`
draws nothing at all. That is pre-existing in the 260530 original, and showing only the 14
significant ROIs is arguably the better figure, so I did not change it. Say the word if
you want all 1236 drawn.

---

## 6. Coherence vs Rayleigh on the linear phase — DONE

`Ventral_surface_polar_coh_vs_rayleigh_260729.m` →
`Ventral_surface_summary\polar_coh_vs_rayleigh_260729\`. Two polar panels, **same 368
ROIs, same phase, scored two ways** — the phase-domain counterpart to the absolute-time
PETH.

**The phase is the point.** Both panels read `piecewise_phase_local`: feet at 0/2π, peaks
at π, **linear ramps in time between consecutive landmarks**. Interpolating peak↔trough
warps every cycle onto its own [0, 2π] regardless of duration, which is what makes phase
usable when breathing is irregular — a long breath and a short breath land in the same
phase bin. That is the middle ground between assuming a fixed breath frequency and
dropping phase entirely.

- **A — coherence**: `coherencyc(cos φ, spikes)`, r = band-mean |C|, dashed confC. A
  waveform measure; needs a linear relationship at breath frequency.
- **B — Rayleigh**: the phase *at each Ca event*, no spectral step, no band choice.
  r = **log Z**.

### Two things had to be got right

**Occupancy is ~7× uneven, so correcting it is not optional.** Measured on Vglut2/1124:
27–47 frames per 10° bin through the fast inspiratory half, **281–304 through the slow
expiratory half**. A cell firing uniformly in *time* would look strongly
expiration-locked. So the resultant is built from the events/frames rate profile. The
angle stays in raw linear-phase units, so 0 = onset and π = peak in every recording and
angles remain poolable and comparable with panel A. (Remapping phase by its own ECDF also
debiases, but then the peak lands at a different angle per recording and the pooled plot
stops meaning one thing.)

**A weighted R̄ with raw n breaks the analytic Rayleigh Z — I hit this and fixed it.**
`Z = n·R̄²` assumes n independent, *equally weighted* angles. With 7× weights, a couple of
events in a sparse inspiratory bin get huge leverage, and the first run produced
**p = 5e-79 where the circular-shift shuffle said p = 0.05.** The fix is Kish's effective
sample size, `nEff = (Σw)²/Σw²`, which reduces to n under uniform occupancy:

| | before | after |
|---|---|---|
| median(log₁₀ p analytic − shuffle) | wildly negative | **−0.13 decades** |
| Spearman ρ, analytic vs shuffle | — | **0.957** |
| median nEff / n | — | 0.66 |

Worth ruling out explicitly: `spike_train` is binary with single-frame events (run length
always 1), so `nnz` is a true event count — event inflation was *not* the cause.

The analytic threshold is still mildly liberal (83 vs 52 at α=0.05), so the dashed circle
is labelled as a guide and the shuffle-confirmed count is printed next to it. **The
authority is `ray_p_shuffle` in `polar_scores.csv`**, which has all 368 ROIs with both
scores, both p-values, `n_eff` and a BH flag.

Why log Z earns a single dashed circle: p ≈ exp(−Z) so Z_crit = −ln α is *n-independent*,
unlike R̄ whose critical value moves with event count. And Z ≤ nEff means an ROI needs ≥4
events to reach α=0.05 at all, so the statistic is self-limiting and no extra event gate
is imposed.

### Results

**Rayleigh is far more sensitive than coherence**: 65 Rayleigh-only vs **1**
coherence-only, 13 significant in both, 289 neither. That ordering makes sense — coherence
is band-limited and needs a waveform relationship, Rayleigh only needs phase
concentration.

**The two agree on phase**, which is the real reassurance: circular-mean offset
**−0.33 rad (sd 0.56)** over the 13 significant in both, −0.23 rad over all 368. So the
phase estimate is robust to the choice of score.

Shuffle-confirmed circular-mean phase per group:

| group | n | phase | concentration |
|---|---|---|---|
| **Vglut2** | 8 | **+1.66 rad (95°)** — mid-inspiration | 0.96 |
| **Sert** | 26 | **−3.10 rad (178°)** — at the inspiratory peak | 0.72 |
| ChAT | 2 | +3.10 rad (177°) | 0.74 |
| IO | 16 | −1.88 rad (−108°) — expiratory, loose | 0.56 |

**Sert lags Vglut2 pFN by ~85° of the breath cycle** — Vglut2 fires on the inspiratory
rise, Sert at the peak. Note the Sert n of 26 here is far above the 4 that cleared the
coherence threshold in §5d, which is the sensitivity difference above, not a
contradiction.

### Panel B marks: circular median + IQR, not the resultant vector

Each dot is the occupancy-**weighted circular median** of that ROI's event phases and the
arc is its weighted circular **IQR** — a direct summary of the per-ROI phase histogram,
which is the object you actually read. Median by grid search (no closed form exists on
the circle); quartiles from the weighted percentile of `angle(exp(1i*(a−med)))`. Both use
the continuous event phase, not the 10° bin centre.

Arcs are drawn **only for the 47 shuffle-confirmed ROIs**. Drawing all 368, or even the 78
above the liberal analytic circle, is mostly near-uniform cells whose IQR approaches the
full circle; they pile into concentric rings that hide the tuned cells. Every ROI is still
plotted as a dot.

| group | n arcs | circular median | median IQR width |
|---|---|---|---|
| Vglut2 | 8 | +97° | 61° |
| Sert | 25 | +175° | 75° |
| ChAT | 2 | +171° | 50° |
| IO | 12 | −104° | 68° |

Weighted vs unweighted median differs by only ~18° in practice, so the occupancy
correction sharpens the statistic without distorting the reported angle.

### Pooled per CELL, and the coherence panel hidden

**The plot is now one panel: log Z on the linear phase.** `showCoherencePanel = false`.
Coherence is still *computed*, so `coh_r` / `coh_th` survive in the CSV and the comparison
stays available — set the flag true to bring panel A back.

**And it now pools per cell, which it did not before.** A cell re-imaged at a second Z
appears once per recording, and scoring each appearance separately double-counted it *and*
made significance depend on which recording you looked at:

```
Sert cell 57   logZ [0.98, 2.24, 2.83]   p_shuf [0.166, 0.006, 0.020]
Sert cell 17   logZ [1.36, 0.18]         p_shuf [0.036, 0.214]
```

Phase is exactly the coordinate where pooling is legitimate — 0 = onset and π = peak mean
the same landmark in every recording, which absolute time does not. So event phases are
concatenated across a cell's recordings, **each keeping its own recording's occupancy
weights** (cycle shapes differ), and `n_eff` accumulates so a cell seen three times is
scored with the power of all three. The circular-shift null shifts each recording's train
independently and re-pools. Identity comes from the same `cell_link_sources` registry as
the cartoon map; `tossed` rows are dropped.

**367 ROI-observations → 358 cells.** Only Sert actually merges (**67 → 58**, 7 cells span
more than one recording) — the other sessions have no matcher run, so pooling is a no-op
there. The conflicts resolve into one number each:

| cell | rec | events | n_eff | logZ | p_shuffle | before |
|---|---|---|---|---|---|---|
| c57 | 3 | 298 | 131 | 2.82 | **0.004** | 0.166 / 0.006 / 0.020 |
| c159 | 3 | 90 | 38 | 2.33 | **0.022** | 0.08 / 0.22 / 0.07 — n.s. in *every* one |
| c17 | 2 | 84 | 46 | 1.38 | **0.018** | 0.036 / 0.214 |
| c13 | 2 | 47 | 27 | −1.94 | 0.840 | consistently n.s. |

c159 is the case worth noting: not significant in any single recording, significant pooled.
That is the power gain, not a threshold artefact.

**How the collapse actually works** — individual EVENTS are concatenated, not
per-recording scores. For each of a cell's observations, append to one growing list: the
phase φ at each event (from that recording's own φ(t)) and its weight 1/occ(bin) (from
that recording's own occupancy histogram). Then a single weighted resultant over the whole
list. Worked example, Sert cell 57:

| recording | events | Σw | Σw² | n_eff | median |
|---|---|---|---|---|---|
| ro1 z10 | 119 | 0.703 | 0.01058 | 46.7 | −178° |
| ro1 z15 | 125 | 0.534 | 0.00456 | 62.4 | −123° |
| ro3 4x | 54 | 0.305 | 0.00301 | 30.9 | −178° |
| **pooled** | **298** | **1.541** | **0.01815** | **130.9** | **−171°** |

`R̄ = 0.3576`, `Z = 130.9 × 0.3576² = 16.74`, `logZ = 2.82`. n_eff is close to additive
(46.7+62.4+30.9 = 140 vs 130.9 pooled), which is where the power comes from.

Because `w = 1/occ`, a longer recording has *smaller* weights but proportionally *more*
events, so Σw is set mainly by firing rate while n_eff grows with duration:
**direction is rate-weighted, power is duration-weighted.** That is the behaviour you want.

**The assumption checked.** Pooling presumes a cell's phase preference is stable across its
recordings; if it drifted, concatenating would blur a real preference into a weak one.
Across all 7 multi-recording cells the per-recording medians agree to within 42° — well
inside the typical 61–77° IQR — except c13 at 80°, and c13 is not modulated at all
(logZ −1.94, p 0.84), so its per-recording medians are noise directions. Nothing is being
blurred.

### The GCaMP lead correction was REMOVED

`ca_lag_sec` was `0.1` — every event shifted 3 frames earlier before scoring, on the
argument that detection runs on raw dF/F (`spike_params.detect_on = 'Raw dFF'`,
`snapWin 0.2 s`), so an event time lags the underlying spike. The argument is sound; a
**fixed time shift** is not, because these recordings span a 6× range of breath rate while
phase is normalised per cycle. The same 100 ms is worth:

| group | median cycle | ° if event in INSP | ° if in EXP |
|---|---|---|---|
| Vglut2 | 0.47 s | **103°** | 62° |
| Sst | 1.03 s | 84° | 25° |
| ChAT | 1.33 s | 63° | 17° |
| IO | 1.58 s | 57° | 15° |
| Sert | 2.15 s | 65° | 10° |

Range across all 36 recordings: **120° in the fastest breather, 7° in the slowest**. So it
was not a common correction — it rotated the groups differentially (Vglut2 +88°, ChAT +51°,
Sert +40°, IO +1°) and reshaped the between-group ordering.

It was caught on a single cell. `Vglut2/0224/left_pFN_roi1` ROI 3 plainly fires *after* the
breath peak in the raw trace, and the older per-ROI figure reads 180°/215°, but the polar
plot put it at 125°. Measured directly:

| | coherence phase | coherence r | Rayleigh median |
|---|---|---|---|
| lag = 0 | **207°** | 0.900 | **220°** |
| lag = 100 ms | 125° | 0.900 | 120° |

`r` is identical either way, as a pure time shift must leave coherence magnitude alone.
Removing the shift also made that cell far stronger — `logZ 1.64 → 2.81`,
`p_shuffle 0.327 → 0.002`.

**Phases are now those of the detected dF/F event** — what is actually measured, with no
model-dependent correction. `ca_lag_sec = 0` restores nothing; set it back to 0.1 to
recover the old behaviour, but note it is not a rigid rotation of the plot.

Shuffle-confirmed **cells** per group (no lead correction):

| group | confirmed / total | median phase | median IQR |
|---|---|---|---|
| Vglut2 | 8 / 17 | **+184°** | 51° |
| Sert | 19 / 58 | **+209°** | 71° |
| ChAT | 1 / 3 | +222° | 60° |
| IO | 7 / 272 | +260° | 59° |

Everything now sits at or just after the inspiratory peak (180°), in the order
**Vglut2 → Sert → ChAT → IO**, which is both cleaner and consistent with what the raw
traces show. Total shuffle-confirmed went 45 → 36; some cells gain, some lose.

For reference, the same table with the old 100 ms shift was Vglut2 +97°, Sert +168°,
ChAT +171°, IO +261°.

Outputs split in two: **`polar_scores_percell.csv`** (what the panel plots) and
`polar_scores_perobs.csv` (unpooled coherence, for reference).

Dashed circles are **significance thresholds** — α = 0.05 at logZ 1.10 and α = 0.001 at
logZ 1.93 — now labelled on the figure with a white background box, since they cross the
data cloud.

**One trap worth recording:** panel A's arc is a jackknife *confidence interval* on the
phase (median 21° wide) while panel B's arc is the *dispersion* of the data (median 61°).
They differ ~3× and must not be compared as if they were the same thing. Hiding panel A
removes the hazard.

### Worked case: the ramping Vglut2 cell

`Vglut2/1124/roi5_1400-1230-0_x4.4_15lp_6000f_00001` **ROI 5** — the largest ramp in the
whole Vglut2 set (`ramp = +0.371 dF/F per s`, peak amplitude 0.463). It **is** in the polar
plot and is one of the strongest points: `coh r = 0.980` (confC 0.83), `logZ = 3.70`
(crit 1.10), `p_shuffle = 0.050`, 96 events, median **+72°**, IQR +45…+106°.

It is *not* pre-inspiratory: dF/F leaves baseline at **−0.23 s** (just before onset) but
half-rises at **+0.10 s** and peaks at **+0.27 s**; events sit a median **+0.30 s** after
onset. A ramp that begins slightly before onset, inspiratory-locked rather than pre-I.

This cell is also the cleanest demonstration of why occupancy correction is mandatory. Its
cycle is 2.83 s and inspiration occupies only **11%** of it (669 of 5902 frames), yet **65
of 96 events** fall in that 11% — an event rate of 0.097/frame in inspiration against
0.006/frame in expiration, **16×**. On raw counts alone the 65/31 split looks unremarkable;
normalised by occupancy it is a strong inspiratory lock.

---

## 7. What the breath waveform looks like in phase coordinates

`breath_amplitude_by_phase_260730.m` → `polar_coh_vs_rayleigh_260729\breath_amplitude_by_phase.*`.
Per-recording SVD PC1, **z-scored within each recording** (shape on a common vertical
scale, deliberately not the acquisition-comparable fixed-metric version), binned into 72
phase bins of 5°. 36 recordings.

**All 36 collapse onto essentially the same curve** — that is the main result, and it is a
strong validation of the peak/trough GUIs: if a recording's landmarks were wrong, its
waveform would not sit on top of the others.

**The profile is built by PER-CYCLE INTERPOLATION, not by binning frames.** This matters
and the first version got it wrong. Inspiration is ~21% of a ~1.5 s cycle, so only ~8.5
frames span 180° — **~21° between consecutive frames against 5° bins**. Most inspiratory
bins are therefore empty in any given recording, and *which* ones are populated depends on
that recording's frames-per-inspiration. The grand mean at a sparse bin is then an average
over a small, biased subset of recordings, which produced a spurious **"ditch" at 172°**:
that bin had data from **5 of 36 recordings** while its neighbours had 15 and 36 — and
those same 5 read 2.89 at 160° where all 29 read 3.26, so they were systematically low.

Interpolating each *cycle* onto the common grid instead gives every bin a value from every
cycle: **72 of 72 bins now have all 36 recordings**, 3724 cycles total (median 78 per
recording). The ditch is gone and the waveform peaks at 178° and troughs at 2° — i.e.
exactly on the landmarks, which is the sanity check that the fix is right (binning had put
the peak at 168°).

Second bug behind the same fix: the within-cycle phase must come from **unwrapped**
`phi − phi(cycle start)`, not `mod(phi,2π)`. With the wrapped version a cycle's last frame
folds back to ~0 and its first can land at ~2π−ε from floating-point `mod` of `2kπ`, either
of which breaks monotonicity and makes the cycle get discarded — only ~11 of ~80 cycles per
recording survived.

Grand mean, z units:

| | |
|---|---|
| expiratory plateau | −0.47 |
| peak | +3.19, at **178°** |
| rises through 10% / 50% | 68° / 112° |
| falls back through 50% / 10% | 192° / 212° |
| rise occupies | 110° of phase |
| fall occupies | **35°** of phase |

**Inspiration averages 21% of the cycle** (range 10–39%) yet is given 180° of the phase
axis, so the warping stretches inspiration ~2.4× and compresses expiration ~1.6×. That is
why the rise looks gradual and the fall abrupt. Occupancy per 5° bin: 0.57% in the
inspiratory half vs 2.21% in the expiratory half, **3.9×** — the imbalance that forces the
occupancy correction in the Rayleigh.

**A limitation this exposes.** From ~220° to 360° the waveform is a flat plateau at
baseline: roughly 40% of the phase axis carries no breath signal at all. Phase is still
well *defined* there (it is a linear position in time between the peak and the next
onset), but it is not anchored to any feature of the waveform. So a cell whose median
phase lands at, say, 260° is firing *somewhere in the expiratory pause*, and the precise
number should not be over-read. This matters most for IO, whose median sits at −99°
(= 261°), squarely in the plateau.

**Why some recordings look jagged in the occupancy panel — it is breath rate, not
genotype.** Jaggedness in the expiratory plateau tracks "bins per frame" at
**Spearman ρ = +0.949**, essentially deterministically:

| | cycle | frames/cycle | bins between frames | jaggedness |
|---|---|---|---|---|
| Vglut2 0224, IO 0124 | 0.43–0.47 s | 13–14 | **5.1–5.5** | 0.68–0.88 |
| Vglut2 0728 | 0.60–0.67 s | 18–20 | 3.6–4.0 | 0.52–0.57 |
| Sert / Vgat | 1.9–2.8 s | 56–91 | 0.8–1.3 | 0.05–0.15 |

At a fixed 30 fps a 0.47 s breath gives only 14 camera frames per whole cycle spread over
72 bins, so ~5 bins sit between consecutive frames: most bins can never be hit and the
ones that are take all the counts. It reads as "Vglut2" only by coincidence —
`Vglut2/1124` (2.83 s cycle) is among the smoothest and `IO/0124` at 14 frames/cycle is
just as rough.

**This does NOT affect the Rayleigh**, which bins at 36 not 72. At 36 bins the fastest
breather has *zero* empty bins (min 4 frames, max 303) where 72 bins had 4 empties, and
`n_eff/n` is actually higher for fast breathers (median 0.72) than slow ones (0.52),
because their inspiration is a larger fraction of the cycle (35–39% vs 10–16%) so the
occupancy weights are more even. The aliasing is confined to the display histogram.

**Where each type fires on the waveform** (panel 4, shuffle-confirmed cells):

```
Vglut2   65-130 deg    on the inspiratory RISE, before the peak
Sert    110-230 deg    spanning the PEAK
IO      230-270 deg    after the fall, in the expiratory plateau
```

A clean sequence: Vglut2 → Sert → IO.

---

## Bugs fixed on the way

| where | what |
|---|---|
| `archive_into_summary` | `readtable` parses `"0721"` as the number **721**; `char(721)` would have silently created a garbage folder name |
| `archive_into_summary` | a blank csv cell reads back as `<missing>`, not `""`, so the skip test never fired and `char(<missing>)` errored |
| `percell_8panel` | `o.spikes` is stored **sparse**; `double()` preserves that and `fft` rejects it |
| `percell_8panel` | the spike-triggered window was in **frames**, so trials from 30 / 42 / 47 fps recordings could not be stacked. Now a fixed 25 ms time grid |
| `percell_8panel` | the event-histogram panels inherited the spectrum's 0–14 axis. Triggered-panel limits are now enforced after all plotting |
| `Ventral_..._cartoon_260729` | a Python patch turned `\260` into an octal escape (`°`) inside a comment |
