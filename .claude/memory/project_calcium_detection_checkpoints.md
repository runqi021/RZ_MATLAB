---
name: project_calcium_detection_checkpoints
description: "Checkpoint log for the calcium dF/F event detection effort — what is SETTLED, what is OPEN, and which avenues are measured-closed; read this first before resuming"
metadata: 
  node_type: memory
  type: project
  originSessionId: a1713ea0-2ee8-4943-8d90-ba971a8be969
  modified: 2026-07-24T01:51:26.292Z
---

Running checkpoint log for pushing calcium event detection toward the physical limit. Agent: `.claude/agents/calcium-dff-detector.md`. Companions: [[project_ca_spike_detector_260722]], [[project_cascade_gcamp8s_detection]], [[reference_oasis_on_this_machine]].

**Read the SETTLED and CLOSED lists before proposing anything — several of these cost a full failed cycle to establish.**

---

## CHECKPOINT 1 (2026-07-22) — the kernel needs no labels
Measured the GCaMP8s decay from isolated ≥4σ hand-labelled events, then reproduced it from the **autocorrelation alone**:

| group | supervised τ | label-free τ |
|---|---|---|
| ChAT | 0.49 s | 0.41 s |
| Sst | 0.64 s | 0.69 s |
| Vgat | 0.48 s | 0.48 s |
| Vglut2 | 0.91 s | 0.63 s |

**Significance:** the whole pipeline can run unsupervised. Also: per-session τ varies as much as between cell types (Vglut2 0.36→1.20), so estimate per session, not per group. τ must be fit on autocovariance lags ≥1 — lag 0 is the only lag white noise touches.

## CHECKPOINT 2 (2026-07-22) — motion is NOT the artifact source
Motion-speed AUC at labelled events: **median 0.500** over 25 sessions. `sx`/`sy` from NoRMCorre carry no information about which events are real.

**Significance:** killed the entire premise that a motion-shift veto would work. Artifact evidence must be spatial (core/ring, ∇F projection) from the movie. Only `Vglut2/1124/cell1` (0.65) and `Vglut2/0224/cell1` (0.62) show real motion-locking.

## CHECKPOINT 3 (2026-07-22) — labels are a reference, not truth
`ChAT/0521/cell1` and `ChAT/0523/cell1` are largely noise, established three independent ways: sign-flipped FDR ≈ 1.0; zero usable events in the kernel fit; CASCADE finds 4 events where 219 were marked.

**Significance:** when a detector disagrees with the labels, decide which is wrong before assuming it is the detector. 12.8% of all labels sit below 3σ.

## CHECKPOINT 4 (2026-07-22) — best supervised-agreement detector so far
OASIS AR(1), per-session g, threshold 1.5×sn, lag +4 frames, tolerance ±6 frames: **per-ROI median F1 0.96**, precision 0.91 / recall 0.76 over 100 sampled curated ROIs.

**Caveat:** tuned ON these labels, so in-sample. Fails on dense/rhythmic ROIs (~35% of ROIs fire >20 events/min).

## CHECKPOINT 5 (2026-07-22) — CASCADE runs, out of sample
Deep net trained on real ephys ground truth; needs no labels. `GC8s_EXC_30Hz_smoothing25ms_high_noise` for all groups. Agreement with hand labels: per-ROI median F1 0.55 (out of sample, and it counts APs where labels count events). Detections are sound — median spike-triggered rise **7.5σ**, only 4.7% below 1σ.

**Significance:** an independent, label-free detector exists and its disagreements are diagnostic.

## CHECKPOINT 6 (2026-07-22) — within-mask reweighting is CLOSED
Extracted per-ROI pixel×time matrices from the registered movie and tested SNR-optimal weighting.

- first spatial PC: **0.14–0.78×** (WORSE — PC1 of a noise-dominated ROI is the noise mode; clipping concentrates weight and cuts effective averaging)
- matched filter `w ∝ cov(y_i,m)/σ_i²`: **0.99–1.05×** (no gain)

**Significance:** cpSAM masks are already tight, so there is no dim-edge dilution to recover. Predicted 1.5–3×, measured none. **Do not redo this.** Also established: `cellpose-gpu`'s BLAS is broken outright (matmul included), so extraction must run in `dlc310`.

## CHECKPOINT 7 (2026-07-23) — FOV low-rank denoising WORKS, and only on the weak ROIs
`pmd_denoise_traces.py`, patch 32 / stride 16, keep components with lag-1 temporal autocorrelation ≥ 0.25.

**Verified the separation is clean:** in a patch containing an active cell, component 0 has lag-1 acorr **0.966** and correlates **+0.996** with the ROI trace; components 1-9 are white (acorr ~0) and uncorrelated. Rank 1 is the right answer per active patch.

**Gain is confined to low-SNR ROIs**, which is exactly where detection was failing:

| raw tSNR | n | median gain |
|---|---|---|
| 0-2 | 1 | 1.81x |
| **2-4** | **32** | **1.24x** (2.6 -> 3.2) |
| 4-8 | 7 | 0.95x |
| 8-16 | 7 | 0.96x |
| 16+ | 1 | 1.00x |

**Significance:** first measured improvement on the failing tail — it lifts events from below 3 sigma to above it. Above tSNR ~4 there is no gain (slightly negative; the reconstruction sheds a little real signal), so **apply it selectively to low-SNR ROIs, not globally**.

**Two implementation traps, both of which silently produced "the method does not work":**
- **Randomized SVD MUST use power iterations** (q=2). Shot noise flattens the singular spectrum, so a single-pass range finder returns noise directions, every component tests as white, all get rejected, and the movie collapses to the per-pixel mean (mean kept rank 0.1/20, tSNR -> 0).
- A full economy SVD per patch is ~30 min/session over ~700 patches; the randomized version is ~100x cheaper (~2-5 min/session).

## CHECKPOINT 8 (2026-07-23) — split-half harness built; and a normalisation bug that invalidates earlier matched-filter results
`split_half_validate.py`. Splits each ROI's pixels into interleaved (checkerboard) halves: same signal, INDEPENDENT shot noise. Gives reliability `corr(A,B)` (Spearman-Brown corrected) = a label-free detectability index, a consensus reference for recall, and an objective for parameter search — on all ~2100 ROIs including the ~1500 never curated.

**Reliability is high and uniform:** Sert median 0.78-0.90, ~100% of ROIs above the 0.10 gate. Ventral similar. So these traces do carry strong shared signal.

**THE BUG, and it reaches backwards:** normalising a matched-filter output by `robust_sn` (frame-to-frame differences) is WRONG. That estimator assumes WHITE noise; matched filtering makes the noise strongly correlated, so it underestimates sigma badly and a nominal "3 sigma" lands deep inside the noise. Symptom: ~85 detections/ROI-min and sign-flip FDR ~0.85 **even on Ventral, where detection provably works** (that cross-dataset check is what exposed it — the fault was the statistic, not the data). Correct scale: for a unit-norm kernel, filtered white-noise sigma equals the INPUT trace's sigma. After the fix, rate fell to a plausible 16-27/min.

**This means the matched-filter arm of CHECKPOINT 5's benchmark was unfairly penalised** — `mf`/`mfs` were z-scored by MAD of their own filtered output. Their F1 (0.53/0.52 vs raw 0.45, OASIS 0.66) is a LOWER BOUND and the comparison should be redone before concluding anything about matched filtering.

**Still unresolved:** sign-flip FDR remains 0.72-0.81 at the chosen threshold, and the objective `recall*(1-FDR)` keeps selecting low thresholds. Either the objective needs to be a proper F1 (precision = 1-FDR), or the zero-mean kernel's negative lobes are generating spurious negative-going detections and inflating the null. Five iterations on this objective without convergence — **stop and think before patching again.**

Also note: a lag-1-autocorrelation gate and a zero-mean kernel are both fine ideas that were each tried and each failed to fix the FDR, because the real fault was the normalisation. Diagnose the scale before redesigning the statistic.

## CHECKPOINT 9 (2026-07-23) — SHIPPED: `ca_detect_ship.py`, both datasets, plus the tau gate
Ran end to end on **2009 ROIs / 49 sessions**. Everything estimated, nothing tuned: per-ROI tau (autocovariance lags >= 1), g, sn (diff-based on RAW dF/F only), OASIS at 1.5*sn, +4 frame lag, events segmented by trace-return.

| dataset | ROIs | usable | events |
|---|---|---|---|
| Sert (D:\260721_Sert_soma_G8s\phys) | 1111 | 683 (61%) | 7610 |
| Ventral_surface_summary | 898 | 653 (73%) | 6837 |

Outputs per session: `ca_events_auto.mat` (v7, MATLAB-ready) + `.npz`; per root `_detect_260723/summary.csv` + gallery.

**THE KEY FINDING — tau is itself a detectability gate, and it is COMPLEMENTARY to split-half reliability.** GCaMP8s somatic decay is physically bounded, so an ROI whose autocorrelation returns tau outside ~0.3-1.5 s has no calcium kinetics in it. Before this gate, noise-only ROIs produced the MOST detections, because a bad tau also corrupts the deconvolution: tau 0.13 -> g 0.77 makes OASIS fire on every noise bump (`ro5 ROI29`: 96 spurious events). Adding the gate cut one session from 1013 to 208 events and the gallery went from "covered in detections on noise" to clean onset marking.

**Split-half reliability did NOT catch these** — `ro5 ROI29` scored rel 0.91 while being pure noise, because both pixel halves see the same shared non-neural fluctuation. Reliability certifies against SHOT NOISE ONLY. Use both gates; neither alone is sufficient.

**Also fixed:** missing reliability now means UNKNOWN, not unusable. Defaulting it to fail silently excluded 807 of 898 Ventral ROIs on the first run purely because `split_half` had not been computed there yet.

**Remaining visible weaknesses in the gallery:** ROIs near the tau gate edge (~0.30-0.35) still look messy; the dense/rhythmic regime is still not well served; and at 1.5*sn some obvious large transients are missed (conservative).

## CHECKPOINT 10 (2026-07-23) — measured quality of the shipped detector (100-ROI pooled gallery)
`ca_detect_gallery100.py` — 100 random usable ROIs drawn across BOTH datasets, out to `_detect_260723/gallery100/` (10 pages PNG + PDF).

**Quality, label-free:** median event amplitude **5.1 sigma** (pct 10/25/50/75/90 = 2.7/3.7/5.1/8.0/10.9); 63% of ROIs have median rise >= 4 sigma; only 2% below 2 sigma. Mark-to-peak offset median -1 frame, so timing is right.

**THE OPEN PROBLEM, now quantified:** recall against *obvious* (>=5 sigma, label-free) transients is only **0.67 pooled**, per-ROI median 0.67, and in **31% of ROIs more than half the obvious transients are missed**. The detector is too conservative at 1.5*sn. This is the single biggest remaining defect and it is measurable without any labels, so it can serve directly as the objective for tuning.

**METRIC TRAP (cost a wrong conclusion):** measuring event amplitude with a baseline window at `x[t-3:t+1]` reports **1.1 sigma** where the true value is **5.1 sigma**. Events are timed at the transient PEAK (detector adds +4 frames to the OASIS onset), so that window sits ON the rising phase. Take the baseline at `t-14..t-8` and search the peak in `t-6..t+12` — valid whether the mark lands on onset or peak. The same formula gave 7.5 sigma for CASCADE only because CASCADE marks onsets.

## NEXT (identified, not yet done) — refine the ROI mask to the signal footprint
The spatial footprint `Vt[keep]` of the signal component is an SNR-optimal mask, and it correlates +0.996 with the true trace. cpSAM outlines are drawn on an average projection and need not match where the signal actually is. This is likely worth MORE than the denoising, is the same computation, and directly answers "slightly adjust ROI to get better signal". Note this is NOT the same as within-mask reweighting (CHECKPOINT 6, closed) — that redistributed weight inside a fixed outline; this changes which pixels belong.

---

## OPEN — the remaining avenue
**FOV-level low-rank denoising (PMD, Buchanan et al. 2019).** Within a patch, calcium is low-rank and temporally autocorrelated; shot noise is full-rank and white. Keep only components whose temporal singular vector is autocorrelated. `pmd_denoise_traces.py` implements this without needing `trefide` to compile. **This is the only untested axis that attacks the low-amplitude tail** where the user's complaint (simultaneous misses and false positives) actually lives.

Also untested: neuropil subtraction, and whether the sliding-median dF/F baseline injects avoidable noise.

## Where the limit actually is
Labelled events: median SNR **7.2σ**. Detections: median rise **7.5σ**. The failures are the **12.8% of events below 3σ** and the ROIs dominated by them. A tail problem, not a wall. If PMD also returns ~1×, the honest deliverable becomes a **per-ROI detectability gate** — the one-sidedness ratio already flags 24% of 868 ROIs as having no calcium-like asymmetry at all — rather than a detector pretending to work everywhere.

## Datasets
- `Ventral_surface_summary` — 26 sessions, 898 ROIs, 362 curated.
- `D:\260721_Sert_soma_G8s\phys` — 23 sessions with dF/F + cpSAM, 11 curated, some 18,000 frames. Added 2026-07-22, not yet characterised.
