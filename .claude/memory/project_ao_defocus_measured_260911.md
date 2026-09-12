---
name: project-ao-defocus-measured-260911
description: "MEASURED on the rig 2026-09-11: the ALPAO DM moves focus 6.76 um per um of Zernike defocus, linear over -14 to +24 um, and reproduced a stage z-stack with 0.65 um RMS without moving the stage. Coefficients are PtV, not RMS."
metadata: 
  node_type: memory
  type: project
  originSessionId: c2f9c05c-e8f0-4298-9f58-57d55a9092fe
  modified: 2026-09-12T01:37:09.238Z
---

Stage 0 of [[project-ao-motion-correction-260911]] is DONE and the DM demonstrably
controls focus. Scripts in `matlab-wroking\ao-motion\` (repo `ao_motion_260911\`):
`dm_motion_calib_260911.m` (gain) and `dm_zstack_demo_260911.m` (proof).

## THE NUMBERS
- **Gain 6.656 um of focus per um of Zernike defocus coefficient** (Z2C row 3),
  from `dm_motion_calib_260911_1824.mat`. The stack demo then measured a slope of
  0.9845 against it, so the **refined value is 6.76 um/um**.
- **Linear over the whole tested range.** Residuals about a straight line show
  scatter, no curvature, at both +-2 coefficient and across the 39 um demo.
- **Lateral cross-talk < 1 px** for coefficient in [-1, +1]. Defocus is
  essentially pure axial -- it will not fight the existing NoRMCorre lateral
  correction. Breaks down at c = -2 (30 px of dx); c = +2 is still fine, so the
  breakdown is asymmetric and probably image degradation, not beam walk.
- **Usable coefficient -2.10 to +3.61 = -14 to +24 um of focus**, measured by
  what the CMD_LIMIT check actually accepted on top of the baseline
  (max|cmd| 0.673). Strongly ASYMMETRIC; the negative side is the tight one.

## *** COEFFICIENTS ARE PtV MICRONS, NOT RMS *** (the long-open question)
[[project-ao-motion-correction-260911]] carried both readings, a factor of ~4 on
every range estimate: RMS predicted 30-55 um of axial range, PtV predicted 8-14.
**Measured ~10.4 um (and -14/+24 asymmetric).** That is PtV. Question closed by
measurement.

## THE DEMONSTRATION (`dm_zstack_demo_260911_1830.mat`, 41 planes, -15 to +25 um)
Same field imaged twice: once by stepping the STAGE, once with the stage PARKED
and only defocus commanded. Every DM plane was then correlated against the whole
stage stack to find which stage z it really matched.

```
matched_z = 0.9845 * nominal_z + 1.575
residual RMS 0.648 um over 39 planes
```
**Sub-micron residual, below the 1 um step size.** The DM puts the focus where it
is asked to, across 39 um, with no stage motion.

Two things to know about that run:
- **Offset +1.58 um**, DM stack sitting above the stage stack. The two stacks ran
  sequentially ~10 min apart, so thermal drift is the leading candidate (too big
  for the 0.205 um stage backlash). UNRESOLVED — a second FOV was being acquired
  to test whether it reproduces. If it does, it is a real defocus zero-point
  offset; if it changes, it was drift.
- **Peak correlation rises MONOTONICALLY with defocus, 0.345 -> 0.526.** Not
  noise. Most likely because DM defocus is not optically identical to stage
  translation: the stage translates the sample through an unchanged focus, while
  the DM reshapes the wavefront and picks up some spherical aberration, so the
  PSF differs slightly with coefficient. **Harmless for motion correction** (what
  matters is imaging the right PLANE, and 0.65 um RMS says it does) but it would
  matter for acquiring z-stacks with the DM.

## What this means for the project
The DM has ~-14/+24 um of axial authority at kHz bandwidth, linear, with no new
optics and no interference with lateral MC. Breathing displacement is expected to
be a few um, so this is headroom rather than a squeeze. **The remaining go/no-go
is measuring how far the sample actually moves axially per breath**, which needs
an animal and is the same measurement shape (z-stack + time series registered
into it).

Why z and not x/y: lateral motion is already recoverable with NoRMCorre; axial is
not, because you imaged different tissue. Worse, it FABRICATES dF/F as structures
drift through the plane — which is a competing explanation for
[[project-breath-locked-optical-artifact]] (r=0.96 with breathing, SURVIVED motion
correction, stronger in background than somata).

Related: [[project-ao-motion-correction-260911]],
[[reference-galil-position-registers]], [[feedback-measure-raw-not-converted]].
