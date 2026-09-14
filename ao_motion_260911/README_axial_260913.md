# Axial (z) motion measurement and correction for the DM — 2026-09-13

## STATUS: the closed loop is VALIDATED on real optics (2026-09-13, 20:09)

Static prep, disturbance injected with the DM, reference −12…+12 µm, 60 Hz loop,
gain 0.7, extrapolation 1.0. **Estimator gain 0.932 ± 0.038.**

| f | sine residual | sharp broadband | sharp fundamental |
|---|---|---|---|
| 0.5 Hz | 8% | — | — |
| 1 Hz | 16% | **31%** | 14% |
| 2 Hz | 32% | **57%** | 40% |
| 3 Hz | 44% | **76%** | 60% |
| 8 Hz | 115% — amplifies | — | — |

**At 1 Hz the loop removes 69% of realistic breathing-shaped motion; at 2 Hz,
43%.** Broadband is consistently worse than the fundamental because brain motion
is not sinusoidal — the loop cleans the breath rate and leaves the sharp
inspiratory transient. A sine-only test claims 16% at 1 Hz and is wrong by 2×.

**Not yet tried, in order of expected payoff:**
1. `EXTRAP_GAIN = 0` — results are still ~1.4× worse than the plain-delay
   prediction at *every* frequency, so the derivative term is costing, not
   helping. Free test.
2. `LOOP_GAIN` 0.7 → 1.0, watching for ringing near 5 Hz.
3. **Frame rate is the fundamental lever.** Residual ≈ `2·sin(πfτ)`, so halving
   τ halves it. 120 Hz (~124 lines) would take 3 Hz from ~44% to ~20% and move
   amplification from 8 Hz to ~16 Hz, above the breathing harmonics. Costs FOV
   height.

**Never run on a live animal.** The mouse died before the measurement session,
so the excursion per breath is still unmeasured and `dm_axial_motion_acq_260913`
and `dm_axial_analyze_260913` have never seen real breathing data.

---

Answers one question: **how far does the brain move in z per breath?** That is the
remaining go/no-go for DM-based real-time axial correction. Everything else about
the DM was settled on 2026-09-11.

## Files

| file | runs where | what it does |
|---|---|---|
| `dm_axial_motion_acq_260913.m` | **on the rig** | acquires four blocks, logs raw frames |
| `dm_axial_estimator_260913.m` | offline | the algorithm: frame → dz, dx, dy, gain |
| `dm_axial_analyze_260913.m` | offline | driver: reference → dz(t) → dz(phase) → verdict |
| `test_axial_estimator_260913.m` | offline, no data | synthetic ground-truth tests |

**Measured performance** (`test_axial_estimator_260913`, run 2026-09-13, R2024a,
synthetic volume of 60 blobs, 2.5 px lateral / 3 µm axial PSF):

| test | result |
|---|---|
| pure axial, noiseless | slope **1.0002**, bias 0.002 µm, residual **0.006 µm** |
| sign convention | +3 µm renders as +3.01 µm — correction is `-dz` |
| 6 px lateral, zero true dz | leakage **< 0.005 µm** |
| ±30% brightness, zero true dz | leakage **< 0.005 µm** |
| Poisson noise, low photon count | bias +0.023 µm, scatter **0.072 µm** |
| `phase_shift` sign + off-by-one | correct |

For scale: the bench DM demo reproduced a stage stack to 0.73 ± 0.13 µm RMS. The
estimator is not the limiting term.

**Two tests that changed the design:**

| test 8 — activity leaking into dz | plain | robust refit |
|---|---|---|
| 5% of cells firing 3× | 0.006 µm | 0.001 µm |
| 15% firing | 0.219 µm | 0.188 µm |
| 40% firing | 0.382 µm | **0.416 µm** |

The robust refit does **not** solve this — it is nearly free and helps at low
activity, but a transient covers a whole cell, so at 15–40% active the outliers
are not a trimmable minority, and at 40% it is worse than plain. **Register on a
structural channel.** On a functional channel, breath-locked *activity* is
indistinguishable from breath-locked *motion* — which is the exact thing being
measured.

| test 9 — vessel geometry | axial SNR | recovered slope |
|---|---|---|
| penetrating vessels (along z) | **0.000** | **0.000** |
| in-plane vessels | 4.111 | 1.073 |

A vessel running down the z axis looks identical at every plane, so it carries
zero axial information however bright it is. **Pick a FOV with vessels running
across the plane, or branch points.** There is no recovering from this in
analysis, so the acq script probes it before the run: two frames 5 µm apart must
actually look different (`r > 0.97` triggers a warning and a 10 s window to
Ctrl-C and move).

## How many scripts do I run, and when?

**Two. One on the rig, one afterwards.** You do **not** press Grab in ScanImage —
the acquisition script calls `startGrab` itself, four blocks' worth, and drives
the stage and the DM between them. Pressing Grab by hand would fight it.

By hand, before running anything:

1. Start ScanImage normally (so `hSI` and `hAOROIctrl` exist in the base workspace).
2. **Load and activate the system aberration correction**, so
   `DMcommand_sys_aber_cor_930` is in the base workspace. The script reads the
   baseline from that variable and errors if it is missing.
3. Select channels for **Save** — both the vessel channel and the GCaMP channel.
   Select at least one for **Display**.
4. Find your FOV, set focus and power. **Choose the field on vessel geometry**
   (see test 9 below), not on how pretty the cells are.
5. Then, at the MATLAB prompt:

```matlab
dm_axial_motion_acq_260913            % ~6 min, drives everything itself
```

It prints the session folder when it finishes. Then, any time afterwards, on any
machine:

```matlab
R = dm_axial_analyze_260913('D:\Data\RZ\ao_axial_260913\axial_<stamp>');
```

`test_axial_estimator_260913` needs no rig and no data; run it whenever you want
to re-check the algorithm.

**On the first session, run block C alone** (set `DO_STAGE_STACK`, `DO_DM_STACK`,
`DO_INJECTION` to `false`) and analyse it, just to prove the TIFF reader works on
a real logged file. Then do the full run. That costs two minutes and protects the
animal.

## Before you press go

- **The system aberration correction must already be loaded and active**, i.e.
  `DMcommand_sys_aber_cor_930` exists in the base workspace. The baseline is read
  from that variable. There is no read-back path on this hardware — `DMctrl`
  caches nothing and `AOdata`'s getters are broken in this install — so if the
  variable is missing the script errors rather than guessing. A fallback to
  `zeros(1,97)` would silently measure everything about a **flat mirror** and
  then flatten the DM on cleanup.
- **A channel must be selected for Display or Save in ScanImage.** Otherwise
  `startGrab` returns with `acqState` still `'idle'`, and the obvious wait loop
  falls straight through and "acquires" nothing. `auto_acq_260910` hit this on
  the rig; the acq script now errors on it.
- Check `CHANNEL` at the top matches the structural label you want to register on.
- `OUT_ROOT` defaults to `D:\Data\RZ\ao_axial_260913`.

## The four blocks

| | block | stage | DM | why |
|---|---|---|---|---|
| A | stage z-stack | **moves** | baseline | the metric ruler — encoder-backed |
| B | DM z-stack | parked | **steps** | the registration **reference** |
| C | time series, 2 min | parked | baseline | **the measurement** |
| D | injected staircase | parked | **known steps** | the positive control |

**Why the reference is taken with the DM and not the stage.** The 2026-09-11 demo
found peak correlation rising monotonically with defocus (0.345 → 0.526): DM
defocus is not optically identical to stage translation, it carries some
spherical aberration, so the PSF differs slightly with coefficient. Harmless for
holding a plane, but it is exactly the mismatch a registration reference must not
have. Block B shares its optical path with block C. Block A only puts a micron
label on the coefficient axis — and doing so **re-measures the gain in tissue**,
which was the outstanding next step from 2026-09-11 (all four bench FOVs were a
static test sample).

**Why block D is not optional.** Without it, a null result in the breathing
analysis is uninterpretable — you cannot tell "the brain does not move axially"
from "the estimator is not sensitive enough."

## The estimator

Not argmax-over-slices. That quantises to the plane step, the correlation-vs-z
peak is broad and shallow (median r was 0.35–0.53 on the bench), and it costs K
correlations per frame. It is kept, but only as a **coarse tier** for lock-on.

The workhorse linearises about the nearest reference plane:

```
I ≈ (1+g)·R₀ + dz·∂R/∂z + dx·∂R/∂x + dy·∂R/∂y + b
```

Five unknowns, one weighted least-squares solve. The weights come from the
*reference* (expected photon noise ~ expected intensity), not the frame, so the
whole normal-equation inverse is precomputed per plane and running a frame is one
5×Npix matrix-vector product.

Three terms that look droppable and are not:

- **dx, dy** — nothing corrects them in real time, but lateral motion projects
  onto `∂R/∂z` wherever structure is oblique. Estimating them is how you keep
  them *out* of dz.
- **g** — defocus dims and blurs, so any laser/PMT/bleaching drift reads as
  defocus. Biggest single confound; costs one column.
- **b** — PMT baseline. SI digitises signed with dark near 0.

**Strips, not frames.** A resonant frame takes ~33 ms — about a tenth of a breath
cycle — so a frame is not a snapshot, it is z-*sheared* top to bottom. Solving per
strip (8 × 64 lines) removes the bias and gives ~240 Hz instead of 30.

**Sign convention.** `dz` is the displacement *of the sample*; the correction to
send is `-dz`. The lateral helper `phase_shift` returns the **correction**, not
the displacement, and is negated where used. Backwards in a closed loop, the loop
drives the motion instead of cancelling it and runs away.

## Motion in the reference stack — measured, and smaller than expected

Block B is acquired while the animal breathes, so each plane is the true structure
convolved with the breathing displacement. The expectation was that this blurs the
axial profile, shrinks `∂R/∂z`, and biases every dz estimate toward zero.

**It does not.** From `test_axial_estimator_260913`, tests 6 and 6b:

| reference | recovered slope |
|---|---|
| sharp (gated) | 0.999 |
| blurred over 2 µm | 0.999 |
| blurred over 4 µm | 1.004 |
| 40 frames/plane, noisy | 1.002 |
| 10 frames/plane | 0.989 |
| 2 frames/plane | 0.995 |

Smoothing is a deterministic linear operation, so least squares with a smoothed
predictor stays *consistent* — regression dilution needs **noise** in the
predictor, not blur. And reference noise washes out too, because the fit averages
over thousands of pixels. Blur costs a constant offset and some sensitivity, not a
gain error.

So **gating is defensive, not load-bearing**, and the acquisition does not depend
on it working. Set `opts.gate = false` and the numbers should barely move; if they
move a lot, something else is wrong and that is worth knowing. Raw frames are kept
regardless (`logAverageFactor = 1`) because the gate is then retunable offline
without re-imaging, and because the per-frame lateral shifts are where breath phase
comes from. That costs disk and nothing else.

The gate uses the **mode** of the lateral displacement, not the median: a breathing
trace dwells at end-expiration, so the density peak *is* end-expiration while the
median sits partway up the excursion.

*Caveat on those numbers:* the synthetic volume is isolated Gaussian blobs at good
SNR. Dense neuropil with a large background pedestal is a harder case, and block D
is what checks it on real data.

## Breath phase comes out of the image

No camera. This machine has no conda and no ffmpeg, so the Basler and thermal legs
don't run here — and it turns out not to matter. The estimator already solves dx,
dy per strip at ~240 Hz, and lateral breathing motion is large and obvious. Phase
from lateral, displacement from axial: two independent quantities from the same
frames, nothing to synchronise.

## The confound the whole session exists to rule out

**A breath-locked signal is not necessarily breath-locked motion.** This prep has
already produced one that wasn't: deep Vgat, top SVD mode at r = 0.96 with chest
breathing, *survived* motion correction, 4× stronger in background than in somata
→ optical, not displacement. Anything of that shape projects onto `∂R/∂z` and gets
reported as axial motion. The gain term only absorbs *global* brightness; a
spatially structured intensity change walks straight through it.

Section **4b** is the discriminator. A real axial displacement must change the
image by `dz·∂R/∂z` — a pattern we have *measured*, not assumed. So the observed
breath-extreme frame difference is projected onto `∂R/∂z`:

- variance explained > 50% → consistent with real axial motion
- 20–50% → mixed; the p2p is an **upper bound**
- < 20% → the optical-artifact signature; **do not quote the p2p as a displacement**

This cuts both ways: if it comes back high, it is also direct evidence about the
breath-locked artifact question, because axial drift of structures through the
plane is the standing competing explanation for it.

## Where this goes next

The measured `dz(phase)` is the input to a **feed-forward** corrector, which is
what the frame rate forces. At 30 Hz, one frame of delay against 3 Hz motion is
36° of phase lag — a purely reactive loop chases the motion a third of a cycle
late and can make things worse. Breathing is phase-locked and predictable, so:
fit `dz(φ)`, clock it out with `sendpatterndata` + `triggeron` (hardware, zero
latency), and keep feedback only for the residual and for drift.

Split by frequency: thermal drift wanders ~1 µm over tens of minutes (that's what
the scattered +1.58/+0.38/+0.72/+0.02 stack-to-stack offsets were). Send DC to the
stage, keep the DM for the 2–5 Hz AC. Protects the stroke budget, which is the
scarce resource.

**Re-centre the operating point before relying on it.** The envelope is −14 to
+24 µm, strongly asymmetric with the negative side tight. A static defocus bias on
the baseline plus a compensating stage park buys roughly ±19 both ways — worth
doing before discovering breathing pushes the wrong direction.

## Traps encoded in these scripts

- `startGrab` can return with `acqState == 'idle'` → guarded, errors.
- Pages interleave channels **per frame**: `page = (frame-1)·nCh + chIdx`. A
  per-frame stride still returns correctly-sized images and never errors — it
  silently averages the *other* label. Parsed from `channelSave`, asserted.
- Walk pages with `nextDirectory`, never `setDirectory(t,k)` in a loop —
  quadratic in page count.
- Frame rate from TIFF `frameTimestamps_sec` first, SI header second, **never**
  hardcoded.
- z is 0.3125 µm/count, so 1 µm is 3.2 counts (−6.25%) while **1.25 µm is exactly
  4 counts**. Block A uses exact steps.
- Stage z-stack is **unidirectional** with a run-up: backlash is 0.205 µm on
  reversal, 16% of a 1.25 µm step, and monotonic stacks never pay it.
- `hann`, `hilbert`, `prctile`, `imgaussfilt` all live in toolboxes that are not
  guaranteed on the rig MATLAB. All four are written out longhand.
- Phase correlation puts the zero-shift peak at index (1,1), not the centre.
  Off-by-one here is invisible in a slope and wrong in every displacement.
- Circular mean for frame-level phase — a plain mean of angles is wrong across
  the ±π wrap.
- `opengl('software')` before `exportgraphics`, and the PNG is verified by **size**
  (an NVIDIA OpenGL stub is a ~9 kB blank), never by the file merely existing.
- Every DM command is checked against the ±1 limit *before* sending; a clipping
  step is **skipped, never clamped** — a clamped mirror is not producing the
  Zernike you asked for. Injection steps are all pre-checked before the
  acquisition starts, not mid-run.
- `onCleanup` restores DM baseline, stage position, and every SI setting touched,
  on error and on Ctrl-C.

## Open / untested

- **The estimator is tested; the acquisition is not.** `test_axial_estimator_260913`
  passes on synthetic data, and all three files pass `checkcode`. But no part of
  `dm_axial_motion_acq_260913.m` has touched hardware, and the SI TIFF reader in
  the analysis has never seen a real logged file — the `channelSave` parsing and
  `frameTimestamps_sec` extraction are the most likely things to need a fix on
  first contact. Run block C alone first (set the other three flags false) and
  check the reader works before spending an animal on a full session.
- `hSI.hMotors.backlashCompensation` still unchecked (one line, no motion). If
  non-zero, the 0.205 µm figure is a residual after SI's own compensation.
- The bench slope scatter (sd 0.052 over four FOVs, 1838 an outlier at 0.87) is
  still unexplained. The bleaching/regression-dilution story was *ruled out* by
  `corr(slope, median r) = −0.75`, which is the wrong way round. Do not repeat it
  as established.
- Whether the DM's spherical-aberration drift with coefficient matters over the
  small range actually used here. Block B vs block A in tissue will show it.
