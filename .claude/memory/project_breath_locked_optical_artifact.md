---
name: project-breath-locked-optical-artifact
description: "Deep Vgat FOV — dominant SVD mode tracks chest breathing at r=0.96, survives MC, and is 4x stronger in background than somata = optical artifact"
metadata: 
  node_type: memory
  type: project
  originSessionId: f145e055-fef7-4364-bc1d-f95f7ce506ec
  modified: 2026-08-12T19:12:47.860Z
---

Found 2026-08-12 on `D:\Ventral_surface_summary\Vgat\0730\deep\roi1_3x_x1350y850_z265_3000f_31lp_00001`
(z265, shiverer) with `SVD\svd_motion_artifact_260812.m`.

## The finding
The LARGEST pixel-SVD mode of the raw movie tracks the **chest-wall breathing
trace** (`breath_pc1.mat`, behaviour camera) at **r = 0.914 at lag 0, 0.96 at
lag −1 frame**. Overlaying them, every breath produces a matching excursion.

It **survives motion correction**: r 0.96 → 0.88, variance 0.41% → 0.32% after
rigid + non-rigid in-plane MC.

## Why in-plane MC cannot fix it
- in-plane translation R² (regress U_k on [Gx Gy] of the mean image) = **0.10**
- defocus / z R² (regress on [∇²I, I]) = **0.02**

So it is neither a rigid shift nor a simple defocus or gain change — which is
exactly why frame-wise in-plane registration leaves it untouched. Do NOT read a
low gradient R² as "clean"; it only rules out in-plane translation.

## The decisive test: background vs somata
GCaMP is only in cells, so background pixels cannot carry a calcium signal.
Using the cpSAM `maskL` (cropped by the SVD's `pad`) and taking background =
non-soma AND below the 40th percentile of non-soma intensity:

    mode-1 modulation as dF/F:  soma 0.0100   background 0.0421
    background carries 420% of the soma modulation

Fractional modulation is **4.2x larger in background than in somata**. Genuine
calcium cannot do that. The mode is optical/mechanical, not neural.

Because it is fractional, **the artifact hits DIM ROIs hardest** — a small
absolute intensity swing is a large dF/F on a faint cell.

## Replication — 2nd deep FOV, same session
`Vgat\0730\deep\roi1_2.4x_x1300y900_z250_3000f_30lp_00001` (z250):

| | raw | mcmc |
|---|---|---|
| mode 1 breath \|r\| | 0.984 @ lag 0 | 0.945 @ lag 0 |
| mode 1 bg/soma | 6.27 | 1.72 |
| shiftR² / defocR² | 0.08 / 0.11 | 0.04 / 0.02 |

Same picture as z265, so this is a property of the PREPARATION, not one
recording. Two refinements from having both FOVs:

- **raw is robust**: mode 1 is background-dominated in raw in both FOVs
  (4.20 and 6.27) — it is an artifact there, not in doubt.
- **mcmc is variable**: 0.76 at z265 but 1.72 at z250. So MC's effect on the
  SPATIAL signature is inconsistent; do not conclude from a single FOV that MC
  made the mode cell-like. The TEMPORAL signature is untouched by MC in both
  (r stays ~0.95, lag stays 0).
- **lag separates the modes**: mode 1 peaks at lag 0 (simultaneous = mechanical);
  modes 2-3 are soma-dominated (bg/soma 0.10-0.59) and lag 170-400 ms, which is
  what a GCaMP transient following respiratory-locked firing would do.

ROI damage is worse here: **31 of 46 ROIs** have |r| > 0.3 with |MC shift|
(max 0.83), vs 17/35 at z265.

## THE KEY RESULT — it is DEPTH-dependent, and MC fixes the surface
Surface control `Vgat\0730\cell1\roi1_3.2x_x1100y1050_z6_6000f_12lp_00001`
(z6, SAME animal + session as the deep FOVs):

| | surface z6 | deep z250 | deep z265 |
|---|---|---|---|
| MC shift range | ±1 px | ±3 px | ±3/4 px |
| raw mode1 breath \|r\| | 0.952 @ 0 | 0.984 @ 0 | 0.960 @ −1 |
| raw mode1 bg/soma | 2.48 | 6.27 | 4.20 |
| **mcmc mode1 breath \|r\|** | **0.028** | **0.945** | **0.879** |
| mcmc mode1 bg/soma | 0.18 | 1.72 | 0.76 |
| mcmc most-breath mode | #6, r 0.63 | #1, r 0.945 | #1, r 0.879 |
| ROIs \|r\|>0.3 vs shift | **0 / 8** | 31 / 46 | 17 / 35 |

**The artifact is in the RAW movie at every depth. In-plane MC removes it
completely at the surface (r 0.95 → 0.03) and fails at depth (stays 0.88-0.95).**
At z6 the breathing motion is small and essentially in-plane (±1 px), which
rigid/non-rigid registration handles; at z250-265 it is large and has an
out-of-plane component in-plane MC cannot touch.

At the surface, mcmc mode 1 becomes 6.85 Hz (cardiac band) and breathing is
gone from the leading modes entirely; the raw mode 1 was 0.37 Hz = the breath
rate.

**Consequence: the Ventral_surface archive that feeds the Rayleigh / polar /
stack analyses is NOT explained by this artifact** — those are surface sites.
The DEEP recordings are the contaminated ones. Check any deep FOV before
trusting its breath phase statistics; project out mode 1 if needed.

Also: this surface FOV holds the Vgat cell `.../1` with logZ 1.80 (2nd-highest
Vgat). Its recording is clean after MC (0/8 ROIs motion-correlated), so that
cell's phase locking is NOT an artifact.

## Why this matters
17 of 35 ROIs in this FOV have |r| > 0.3 between their dF/F and |MC shift|.
A breathing-locked artifact that survives MC is precisely what fakes
breath-locked calcium, i.e. it inflates the Rayleigh phase-locking that
[[project-stack-dff-logz]] and the polar figures select on. Worth running this
check on any FOV before trusting its logZ.

Contrast with the vglut2 surface FOV in [[project-pixel-svd-movie-260811]],
where the top modes were 10 Hz cardiac edge dipoles instead.
