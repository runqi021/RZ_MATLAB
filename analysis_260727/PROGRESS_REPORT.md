# analysis_260727 — progress report

**2026-07-27** · Runqi Zhang / Claude

Breath × calcium analysis rebuilt around cell identity, re-founded on an
absolute-time statistic, and retargeted to the `baseline` condition folder.

---

## 1. What exists now

```
RZ_MATLAB\analysis_260727\
    run_analysis_260727.m         driver: whole chain, any dataset
    cell_metadata_260727.m        master metadata join (ROI-level + cell-level)
    cell_pair_morph\              cross-FOV cell matching + curation PORT
    coh_ca_breath\                coherence (secondary) + the cell-identity join
    breath_time\                  PRIMARY: inspiration-triggered PETH
    phase_rayleigh\               RETIRED — kept and runnable, not in the chain
```

```matlab
run_analysis_260727()                                    % config default
run_analysis_260727('C:\...\Ventral_surface_summary','Ventral')
```

Output always lands in `<dataset>\analysis_260727\`. Cell matching is optional —
with no `roi_match_results.mat` the chain falls back to identity grouping.

**Currently pointed at** `D:\260721_Sert_soma_G8s\phys\baseline`.

---

## 2. The primary statistic

For the observed PETH `r(t)`, with `rs` = smoothed and `rbar = mean(r)`:

```
T_exc = max over the search window of [ rs(t) - rbar ]
T_sup = max over the search window of [ rbar - rs(t) ]
M_exc = ( T_exc,obs - mean(T_exc,shuffle) ) / std(T_exc,shuffle)
p_exc = ( 1 + #{ T_exc,s >= T_exc,obs } ) / ( 1 + nShuffle )
```

Every shuffle recomputes the same quantity **including its own max**, so the search
over candidate latencies is paid for by the null. There is no designated
respiratory baseline; each cell is compared only against shuffled versions of
itself. Excitation and suppression are separate throughout — a signed summary
exists but both originals are kept.

Significance is called on the **empirical p** (rank-based, valid at any n).
`M` is an effect size and is unreliable below ~40 events, where the null of a
max-statistic is discrete and right-skewed.

## 2b. What "zero" is, and the two alignments

The upstream field names are a trap and are worth stating once:

| landmark | file → field | what it is |
|---|---|---|
| `foot_idx` | `breath_insp_start_pc1.mat` → `insp_start_idx` | **inspiration onset** — the foot of the breath PC1, where the rise begins |
| `peak_idx` | `breath_peak_pc1.mat` → **`insp_onset_idx`** | the **inspiratory peak**. The field says "onset" but it is the peak |

They are **267 ms apart** (CV 0.15, SD ≈ 40 ms). `TRIGGER = 'onset' | 'peak'` selects
which is zero; output goes to `breath_time\<TRIGGER>\` so both coexist.

**This is a result, not a robustness check.** Because the onset→peak interval varies,
a cell locked to the peak is smeared by that jitter when aligned on the onset, and
vice versa. Whichever alignment gives the tighter latency distribution is the
landmark the cell is actually locked to.

**Verdict on baseline: 24 of 34 cells are tighter PEAK-aligned, 4 onset-aligned,
6 tied** (among significant: 16 peak / 3 onset / 4 tied). Median latency MAD 264 ms
peak vs 320 ms onset. `M_exc` barely discriminates (3.25 vs 3.33) because it mixes
dispersion with firing rate and with how much evidence the cell provided —
**MAD is the discriminator**, not `M_exc`.

Latency MAD is quantised, so exact ties are common and are now reported as `tied`
rather than silently assigned to whichever side is compared first (that bug briefly
made the count read 28 vs 6).

> **This verdict reversed when latency assignment was fixed.** With events assigned
> to the *preceding* trigger the answer came out 25-of-34 onset-locked; with
> **nearest**-trigger assignment it is 28-of-34 peak-locked. Preceding-assignment
> systematically penalised the peak alignment: an event during inspiration (between
> foot and peak) has no preceding *peak* in its own breath, so it was charged to the
> previous peak at a latency of nearly a whole cycle, inflating peak-MAD for
> exactly the cells that fire around inspiration. Nearest-assignment removes that
> artifact. The earlier onset-locked conclusion was an artifact of the assignment
> rule, not a property of the cells.

## 2c. Latency is signed, and measured to the NEAREST trigger

Each event is assigned to its **nearest** accepted trigger, giving a **signed**
latency, rather than to the preceding one. An event 80 ms *before* an onset is
−80 ms — pre-inspiratory — where preceding-assignment would call it +2.05 s of the
previous breath. That is both wrong and destroys the measured precision of any
pre-inspiratory cell, which in a breathing circuit are the cells of interest.

Each event still belongs to exactly one breath, so nothing is double counted.
Signed latency is bounded by roughly ±half the local inter-trigger interval — a
larger distance implies a nearer trigger — so the precision search runs over
[−cycle/2, +cycle/2], and a preferred latency landing in the second half of the
cycle is wrapped to the equivalent negative latency to the *next* breath
(`pref_latency_nearest_ms`).

Consistency: median onset latency 442 ms − median peak latency 142 ms = 300 ms
against a median onset→peak interval of 267 ms. Within one 50 ms bin, and these are
medians of quantised latencies across cells rather than a per-cell subtraction, so
the agreement is as close as the resolution allows — both runs are describing the
same event in two reference frames.

`breath_time_trigger_compare_260727.m` does this comparison and writes
`trigger_compare.csv/.png/.pdf`. `breath_time_peth_percell_compare_260727.m` puts
both alignments in ONE figure per cell — row 1 onset, row 2 peak, matched columns,
with the per-cell verdict in the title — so comparing them needs no flipping
between folders.

## 2d. The modulation-score tail plot

`breath_time_modulation_tail_260727.m`, one panel per alignment, each showing TWO
survival functions of `M_exc`:

- **grey = the shuffle null**, `P(M ≥ x | H0)`, pooled from every tested cell's own
  null after z-scoring each against itself — 46,240 values;
- **coloured = the observed cells**, fraction with `M_exc ≥ x`.

The dashed α lines at 0.05 / 0.01 / 0.005 are read off the **grey null** curve;
where it crosses each is the critical score, dropped as a vertical guide.

Getting this right matters. On the *observed* survival a horizontal line at 0.05
would mean "the top 5% of the cells I happened to record" — a sample quantile, not
a threshold — and with N = 34 the observed curve floors at 1/34 = 0.029, so 0.01
and 0.005 could never be reached at all. On the null curve they are genuine tail
probabilities.

**Empirical critical scores (baseline, onset):** α = 0.05 → 1.77, 0.01 → 3.05,
0.005 → 4.01, with 23 / 18 / 15 of 34 cells above them. Compare the normal-theory
values 1.645 / 2.33 / 2.58: the null is strongly right-skewed, so assuming
normality would have been badly anti-conservative in the far tail (2.58 vs 4.01).
This is the same skew that makes `M_exc` unreliable as an effect size at low n.

These α guides are population-level. Per-cell significance remains `p_exc`/`q_exc`
in the CSV, from that cell's own null with BH across cells.

## 2e. Per-recording testing, and the simple per-cell figure

`breath_time_percell_hist_260727.m`. Everything above pools a cell's recordings
before testing. This tests **each recording separately** and draws one figure per
cell:

```
rows    = the recordings that cell appears in (+ a POOLED row when >1)
left    = histogram of Ca events vs time from inspiration ONSET
right   = histogram of Ca events vs time from the inspiratory PEAK
```

Just event counts per time bin and a line at zero — no shuffle envelopes, no null
distributions. Each title carries that recording's n, `M_z` and shuffle p, so you
can see whether the modulation **reproduces across recordings** rather than only
the pooled answer. Example, cell 87 in 3 recordings: onset `M_z` = 4.9 / 4.2 / 4.9,
all p ≤ 0.007.

102 (recording × trigger) tests on baseline, 52 significant at BH q<0.05. Written
to `per_recording_tests.csv`. 44 cells get a figure (more than the 34 pooled-tested
cells, because the per-recording floor is 10 events rather than 20).

Two things to read correctly on these panels: the bars are **raw counts** (what
"histogram" means, and easy to sanity-check) while the test underneath uses the
exposure-normalised rate; and the rise near the right edge is the **next breath**
entering the window, which extends to 1.25 cycles.

`breath_time_peth_percell_compare_260727.m` (2×5 panels, shuffle envelopes, null
distributions) is kept for auditing the statistics, not for looking at data.

## 3. Result — baseline (11 recordings, 447 ROIs, 338 cells)

Onset-aligned (the tighter alignment for most cells):
**18 of 34 tested cells significantly excited**, `M_exc` median 5.16 (max 12.33),
preferred latency median **417 ms** after inspiration onset — i.e. ~150 ms after the
inspiratory peak, in early expiration. Precision significant in 23; suppression
testable in only 5 cells and significant in 2.

Recruitment median 0.238, precision 0.345, latency MAD 233 ms (over significant
cells). Peak-aligned gives 19 significant at 142 ms, the same event re-referenced.

| method | significant |
|---|---|
| per-ROI coherence | 4 / 76 |
| **time PETH (primary)** | **18 / 34** |

Baseline breathing is markedly more regular than the mixed-condition set was —
cycle CV **0.19** vs 0.33, inspiration CV 0.15 vs 0.21 — which is what you would
expect once conditions are no longer pooled.

**`corr(M_exc, precision_z) = 0.93`.** The two are effectively one axis: both are
shuffle-standardised and both scale with event count. So the modulation-vs-precision
scatter is largely a single dimension, and **recruitment vs 1/latency-MAD is the
plot that carries independent biology** — those are raw quantities with no shuffle
standardisation, and they do dissociate (cell 209: 33 events, recruitment 0.18,
MAD 67 ms — fires rarely, always at the same moment; cell 177: 136 events,
recruitment 0.33, MAD 600 ms — broadly modulated, poorly timed).

## 4. The data reorganisation, and the curation port

`phys\` was split into `baseline\`, `roomAir\`, `reoxygenated\`. The hand curation
was made on the old 26-recording set, so it had to move with the data.

**It ported exactly — 447 of 447 observations, every centroid bit-identical.**
`roi_curation_port_260727.m` does this generally, keyed on the invariant that a
curation decision is really

```
(recording folder NAME, ROI index)  ->  which cell
```

Neither term depends on where the folder sits in the tree, so a curation survives
**any** amount of moving recordings between parent folders. Re-curating after a
reorganisation is never necessary.

What it will *not* survive is re-running the segmentation, which renumbers ROIs —
so the port verifies that every matched ROI has an identical centroid on both
sides and **refuses** rather than silently attaching decisions to wrong neurons.

Surviving: 186 grouped, 245 ungrouped, 16 tossed → **76 multi-recording cells**
(62 pairs, 12 triples, 1 quad, 1 quint). 17 groups shrank to one member because
their other members are in `roomAir\`/`reoxygenated\`; those stay marked settled.
The original `phys\roi_match_out_260724` was never written to.

## 4b. Cell ids are inherited, not recomputed

`cell_id` was **positional** — the rank of a group within whatever subset was being
analysed — so it renumbered on every subset. Baseline cell 87 was full-set cell
159, and across the 447 shared observations the two numberings agreed **2 times
(0.4%)**. Figures in `baseline\...\percell_hist\` and `phys\...\cell_activity_summary\`
with the same filename were different neurons.

Fixed at the source. `roi_curation_port_260727.m` now writes `curated.srcCellId`,
the SOURCE cell id for every ported observation, and `cell_link_260727.m` and
`cell_activity_summary_260727.m` use it when present. **One neuron now carries one
number in the full set and in every condition folder**, so `baseline`, `roomAir`
and `reoxygenated` will also be mutually comparable.

Ids are consequently **sparse** — a gap is a cell whose members all live in another
condition folder — so every loop over cells must skip empties. `n_cells` counts
non-empty; `max_cell_id` is reported separately.

Verified: cell 159 is the same neuron in every figure folder —

```
baseline\cell_activity_summary\cell_159_03x.png
baseline\breath_time\percell_hist\cell_159.png
baseline\breath_time\onset\percell\cell_159_Mexc012.3.png
baseline\breath_time\peak\percell\cell_159_Mexc009.0.png
phys\cell_activity_summary\cell_159_03x.png          <- full set, same cell
```

Two supporting pieces:

- **`cell_key`** in every table — the anchor member's `<recording>#<roi_index>`,
  e.g. `ro1_2x_x1018y1130_z10_20lp_6000f_00001#25`. Content-derived, so it is
  meaningful even for an analysis that was never ported.
- **`cell_id_crosswalk_260727.m`** for two analyses where one was NOT ported. It
  matches cells by shared members and labels each pair exact / subset / split.
  Baseline vs full set: 287 exact, 51 subset, 0 split, 0 unmatched.

## 5. The metadata tables

Results used to be scattered across four CSVs with no single source of truth.
`cell_metadata_260727.m` joins them once:

- **`roi_metadata.csv`** — one row per ROI observation. Identity (`cell_id`,
  `cell_size`, `status`), **both coordinate frames** (`x/y/z_um` stage, comparable
  across recordings; `cx/cy_px` within its own FOV, for drawing on that image),
  its own activity, its cell's scores repeated on every member row so a spatial
  map needs no join, and its per-ROI coherence.
- **`cell_metadata.csv`** — one row per cell, mean coordinates, all scores once.

Use the ROI table for anatomy, the cell table for statistics: an ROI is not an
independent unit, a cell is.

---

## 5b. The chain builds its own prerequisites, and says what is really missing

Two additions after `260728_vglut2` failed with `cell_pool.mat not found`:

**`ensure_pool_260727.m`** — every script that needs the pooled container now calls
this instead of asserting the file exists. If `cell_pool.mat` is absent it runs
`cell_link_260727` then `cell_pool_260727` itself, and `cell_link` already falls
back to identity grouping (each ROI its own cell) when a dataset has no cross-FOV
matching. So **no curation step is needed for a dataset with few active cells**, and
any downstream script can be run first — the chain assembles behind it.

**`analysis_preflight_260727.m`** — run first on any new dataset. Per recording it
checks meta / dFF / cpSAM / ca_spike_data / breath_pc1 / breath_peak / breath_onset
/ cam .avi, counts how many are READY, and names the upstream step to run. It is now
the driver's step 0, and `ensure_pool` calls it before erroring.

**Why this was needed.** `cell_pool.mat` is DERIVED and would have built itself; the
error named a symptom. The real gap on 260728_vglut2 was that **no recording had a
breathing trace** — all 19 Basler `cam1*.avi` files were present but had never been
through `breath_svd_pc1` and the peak/onset GUIs, so there were no triggers to align
calcium to. Auto-building the pool would only have walked into the next failure.

```
260728_vglut2   meta 19/19  dFF 19/19  cpSAM 19/19  avi 19/19
                ca_spikes  6/19   breath triggers 0/19   -> READY 0 of 19
baseline        ca_spikes 11/12   breath triggers 12/12  -> READY 11 of 12
```

## 6. Bugs found and fixed

**In pre-existing code** (fixed in the `analysis_260727` copies; originals untouched):

| # | where | what |
|---|---|---|
| 1 | `roi_review4_ungrouped` | its save built cells from `grpOf>0` only, so **all 274 ungrouped ROIs vanished** from the cell list and CSV the moment round 4 saved |
| 2 | `roi_review4_ungrouped` | a phase-2 "skip" left `grpOf==0`, indistinguishable from never-reviewed, so reviewed-and-kept ROIs reappeared every run |
| 3 | `coherence_polar_general` | `sig_rois.csv` split its label on `#`, a character the label never contains |
| 4 | `dFF_heatmaps` | `doSave` hard-coded `false`, and its save block referenced three commented-out figures, so it crashed the first time saving was attempted |

**Introduced and fixed during development** (mine):

| # | what |
|---|---|
| 5 | `repoRoot` needed a third `fileparts` — copies sit two levels below the repo root |
| 6 | `struct('rec',rec,...)` with a struct-array value returns a struct *array*, not a scalar struct |
| 7 | several `table()` row/column orientation mismatches |
| 8 | RGBA 4-element `'Color'` in `polarplot`/`plot` |
| 9 | `mod_depth = (max−min)/(max+min)` saturates at 1 whenever any bin is empty — nearly every cell |
| 10 | PETH normalised by **triggers** not frame exposure → a sawtooth from 33.3 ms frames in 25 ms bins. Pure aliasing |
| 11 | PETH window hard-coded −0.5…2.0 s: 1.1 cycles for Sert, **4.4 cycles** for Ventral, where the "baseline" then held the previous inspiration |
| 12 | shuffle count fixed at 1000 → smallest possible p 9.99e−4, but BH rank 1 with 296 tests needs 1.7e−4. **No cell could be called at all** |
| 13 | median baseline is **0** when most bins are empty → peak/baseline reached 3.7e15 |
| 14 | inclusion by rate alone admitted 4-event cells on short recordings (median 7 on Ventral) |
| 15 | **zero-event observations were dropped** — on baseline that is **315 of 389 observations**, so every event rate and every recruitment denominator was inflated by however many recordings a cell was silent in |
| 16 | circular FFT correlogram allowed wraparound between recording end and start. Fixed by putting **only boundary-safe triggers** in the trigger vector, after which the circular correlation *is* the linear one |
| 17 | adaptive bin width made latency incomparable between datasets and quantised it so coarsely the IQR came out one bin wide. Bin width is now fixed, smoothing specified in seconds |
| 18 | `colorbar`'s `Label` getter fires a listener that errors on a scatter with CData in R2021b — label moved into the title |
| 19 | `cell_metadata` still read `breath_time\breath_time_peth_data.mat` after the trigger split, so it silently joined the stale pre-split run. It now takes a `TRIGGER` and asserts the file exists |
| 20 | **per-cell figures came out blank after the 6th cell.** `hf = hOff; clf(hf)` clears an off-screen figure but does **not make it current**, so `subplot()` kept drawing into the last *visible* figure — 28 of 34 exports were an identical 4778-byte blank. Fixed with `set(0,'CurrentFigure',hf)`. The same pattern is still present in `phase_rayleigh_percell_260727.m` (retired, unfixed) |
| 21 | latency measured to the **preceding** trigger rather than the nearest, which mis-assigns pre-trigger events by nearly a whole cycle and **reversed the onset-vs-peak verdict** (see §2b) |

**A correction to something reported mid-session:** "~10% of cycles are junk
low-amplitude detections" used a *global* amplitude median across recordings,
conflating between-recording brightness with within-recording pauses. The correct
per-recording test gives **6 of 2690 (0.2%)** on the mixed set and **0 of 866** on
baseline. Long pauses are real (11.71 s, 6.70 s, 5.17 s) but rare, and detection is
**not** meaningfully splitting them.

## 7. Why phase and coherence are secondary

Cycle-interpolated phase is robust to rate drift, but only **15% of each cycle**
falls in `[0, π)` — half the phase axis — with a 24-bin occupancy max/min ratio
averaging **12×**. An unmodulated cell then produces a large, significant *raw*
Rayleigh vector pointing into expiration, and the ECDF correction needed to fix
that shifted preferred phase by a median of 40°. Under irregular breathing or a
pause, normalised phase also erases the real latency and the real pause duration.

Coherence needs breathing to occupy a stable frequency band and weakens when the
rate drifts, even for a cell locked to every breath.

Absolute time has neither problem: every trigger contributes to every bin, so the
null is flat with no correction. All three were run and agree; the simplest is
reported.

## 8. Known gaps

- `ro1_2x_x1018y1130_18lp_3000f_00001` is in `tossFolders` as a re-zeroed FOV, so
  its 44 ROIs get no cell identity and drop out of the pool. Its calcium and breath
  data are fine — **a hybrid that admits unmatched recordings as identity singletons
  would recover them.** Not built.
- `ro8_2.3x_x-500-1500_z0_27lp_6000f_00001` has no `ca_spike_data.mat`.
- **The PETH does not lead-shift for GCaMP delay; coherence shifts by 100 ms.** And
  `spike_trigger_dFF` reported time-to-peak 0 ms with a 163 ms rise, suggesting the
  detector marks the transient peak rather than its onset — so 417 ms plausibly
  includes ~150 ms of indicator delay. Reported raw; the two pipelines should be
  reconciled.
- Suppression is testable in only 5 of 34 cells. `T_sup ≤ rbar` by construction, so
  at low rates nearly every shuffle hits that ceiling; the gate reports
  `sup_tested = false` rather than a meaningless column.
- Ventral_surface_summary runs but has no cell matching, hence no pooling.

## 9. Next

1. Run the matcher on `roomAir\` and `reoxygenated\`, port the same curation into
   each, and compare conditions — the port makes this cheap.
2. Recover the re-zeroed recording via the identity-singleton hybrid.
3. Reconcile the GCaMP delay between pipelines.
