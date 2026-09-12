---
name: project-ao-defocus-measured-260911
description: "MEASURED on the rig 2026-09-11: the ALPAO DM moves focus 6.76 um per um of Zernike defocus, linear over -14 to +24 um, and reproduced a stage z-stack with 0.65 um RMS without moving the stage. Coefficients are PtV, not RMS."
metadata: 
  node_type: memory
  type: project
  originSessionId: c2f9c05c-e8f0-4298-9f58-57d55a9092fe
  modified: 2026-09-12T01:51:43.964Z
---

Stage 0 of [[project-ao-motion-correction-260911]] is DONE and the DM demonstrably
controls focus. Scripts in `matlab-wroking\ao-motion\` (repo `ao_motion_260911\`):
`dm_motion_calib_260911.m` (gain) and `dm_zstack_demo_260911.m` (proof).

## THE NUMBERS
- **Gain ~7.0 um of focus per um of Zernike defocus coefficient** (Z2C row 3).
  The calibration gave 6.656; four independent stack demos then measured slopes
  of 0.9845 / 0.9688 / 0.8685 / 0.9580 against it (mean 0.945, sd 0.052), i.e. an
  implied gain of 6.76 / 6.87 / 7.66 / 6.95, **mean 7.04**. Quote **~7 um/um
  +-5%**, not three digits. Every run gives slope < 1, so 6.656 is consistently
  an UNDERESTIMATE.
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

## REPRODUCED OVER FOUR FOVs (1830, 1835, 1838, 1842)
| run | slope | offset um | resid RMS | med r |
|---|---|---|---|---|
| 1830 | 0.9845 | +1.58 | 0.648 | 0.481 |
| 1835 | 0.9688 | +0.38 | 0.636 | 0.476 |
| 1838 | 0.8685 | +0.72 | 0.912 | 0.514 |
| 1842 | 0.9580 | +0.02 | 0.709 | 0.451 |

**Residual RMS 0.73 +- 0.13 um across four independent fields** — that is the
robust result, sub-micron every time.

- **The OFFSET IS DRIFT, settled.** +1.58 / +0.38 / +0.72 / +0.02, scattered
  about zero with no repeated value. A real defocus zero-point offset would have
  reproduced. It is thermal drift between the two sequentially-acquired stacks.
- **Slope scatter (sd 0.052) is NOT explained.** 1838 is the outlier at 0.87.
  The obvious story — regression dilution from bleaching, where a noisier match
  pulls estimates toward the stack centre and shrinks the slope — is ruled out by
  the data: `corr(slope, median r) = -0.75`, i.e. 1838 has the HIGHEST
  correlation and the LOWEST slope, the wrong way round. With n=4 that is weak
  evidence, but do not repeat the bleaching explanation as if it were
  established. Something field-specific, possibly how structure is distributed
  in z. OPEN.
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
