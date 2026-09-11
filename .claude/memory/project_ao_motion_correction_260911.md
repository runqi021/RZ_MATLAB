---
name: project-ao-motion-correction-260911
description: "Real-time motion correction using the existing ALPAO DM instead of new scanners — plan, measured stroke budget, and the calibration script (stage 0 of 5, nothing run on hardware yet)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0c9f4688-9d18-4f93-9896-d7dc4ce25892
  modified: 2026-09-11T19:39:54.229Z
---

Goal: real-time (closed-loop) motion correction on the 2P rig, prompted by Kaspar
Podgorski's SLAP2. **Status 2026-09-11: PLANNING + stage 0 written, NOTHING has run on
hardware.** Do not describe any of this as working.

## The key finding: ScanImage 2018b ALREADY ships the framework
No porting needed, no 2026 upgrade needed. Present in the install:
`+components/MotionManager.m` (104 KB), `+motionEstimators/{Gpu,Marius,Parallel,Simple}
MotionEstimator.m` + `+suite2P/`, `+motionCorrectors/{Simple,Marius}MotionCorrector.m`,
`+interfaces/IMotion{Estimator,Corrector,EstimatorResult}.m`. `hSI.hMotionManager` is in the
component tree. `GpuMotionEstimator` takes a reference VOLUME and reads Z off which slice
correlates best (`dz = obj.zs(zidx) - obj.zs(zIdxMaxCorr)`), exactly the SLAP2 shape.

## THE BLOCKER IS ACTUATORS, NOT SOFTWARE
`MotionManager` applies corrections through three paths:
`getMotorCorrection()` (Galil, slow), `getGalvoVoltCorrectionXY()` (needs an X **and** Y
galvo), `getFastZVoltCorrection()` (needs FastZ).
From this rig's MDF: **`galvoAOChanIDX = []`** (NO X galvo — resonant-X + Y galvo only, "RG"
not "RGG") and **`actuators(1).controllerType = ''`** (NO FastZ). So two of three axes have
no fast actuator and fall back to the stage. That is why both shipped correctors default to
`runningAverageLength_s = 3`, `correctionInterval_s = 10` — they are DRIFT correctors.
Slow-drift correction works today with zero new hardware. Breathing-rate does not.

Adding an X galvo = resonant + 2 galvos (not 3). Cheap electrically (one AO channel on
PXI1Slot4; the MDF's `xGalvoAngularRange`/`galvoVoltsPerOpticalDegreeX`/`galvoParkDegreesX`
are already filled in, only the chan ID is empty). Expensive optically: a relay to conjugate
the resonant mirror onto it, **which moves where the DM and Shack-Hartmann sit optically and
would likely force a Z2C recalibration**. That is the real cost, not the galvo.

## The AO route (chosen: "lets try the ao first")
Use the ALPAO **DM97 BAX331** already in the beam and already conjugate to the back aperture.
Tip/tilt -> lateral, defocus -> axial, kHz bandwidth, no new optics, no risk to the AO arm.
For 2P specifically defocus is free: no descanned pinhole, so shifting the excitation focus
needs no detection-side change.

**It is complementary BY CONSTRUCTION.** The system aberration correction
(`water_lens_System_aberration_correction_8mm\SystemCorrection_20260117_930_Iter30_FITC_Z60_
P3_20K_water.mat`, 8 mm water lens, 930 nm, FITC, 2026-01-17) has `Zstart = 4` — the
optimisation *starts at Zernike mode 4* — and `zernike_current(1:3) = [0 0 0]` exactly.
Tip/tilt/defocus are not aberrations, so they are deliberately excluded. Commands simply add:
`cmd = DMcommand_sys_aber_cor_930 + [dTiltY dTipX dDefocus] * z2c(1:3,:)`.
(A newer 10x variant exists: `..._8mm_10x\...20260326_930_...`. Check which is loaded.)

**What they DO share is stroke.** The system correction already uses `max|cmd| = 0.673` of
the +-1 range. Remaining headroom, computed per actuator (asymmetric — the baseline is):
| mode | on top of the system correction | on a flat mirror |
|---|---|---|
| tiltY   | -3.63 to +3.59 um | +-7.98 |
| tipX    | -3.37 to +4.58 um | +-7.65 |
| defocus | -2.48 to +4.22 um | +-5.63 |

## Z2C mode identity — MEASURED, not assumed
`BAX331-Z2C.mat` is [96 modes x 97 actuators]; coefficients are in **microns**
(`Vmatrix5micro = 5*eye(30)` in `DM_SHWS_calibration_Ctr.m`). Identified each row by laying
its 97 commands back onto the DM97 grid (11x11 minus corners = 97) and projecting onto
x-ramp / y-ramp / r^2:
**row 1 = tilt Y (0.988), row 2 = tip X (0.986), row 3 = defocus (0.959).**
Rows 4-5 astigmatism, 6-7 coma (they project partly onto tilt, as coma should), row 10
spherical. **RMS-vs-PtV is UNRESOLVED** — a factor of ~4 in every range number. Optimistic
(RMS): ~+-9 um lateral, 30-55 um axial. Pessimistic (PtV): ~+-2.3 um lateral, 8-14 um axial,
which could be MARGINAL for breathing. The calibration settles it.

## Latency / why feed-forward may be needed
At 30 Hz one frame of delay against 3 Hz motion is **36 deg of phase lag**. SLAP2 ran at
200 Hz to make this negligible. A purely reactive corrector at this frame rate chases the
motion a third of a cycle late and could make things worse. But breathing is phase-locked and
predictable, and the breath-phase modelling is already done — so the right corrector here may
be **feed-forward** (predict from breath phase, apply ahead of the frame) with feedback only
on the residual. `DMctrl.sendpatterndata` + `triggeron` clocks out a PRELOADED sequence in
hardware (zero latency) and `senddata` sends an arbitrary vector immediately — the hybrid
uses both. A preloaded pattern ALONE cannot track cycle-to-cycle variation (sniff vs basal).
DM is on a **PEX-292144 PCIe** card: a send is tens of us, settling << 1 ms. The DM is never
the bottleneck; the frame rate and the estimator are.

## The pipeline
```
1 park, take ONE dense z-stack with the Z motor   -> axial ruler AND registration reference
2 step the DM, register each frame into it        -> gain matrix G      [script written]
3 custom IMotionCorrector: d -> G^-1 -> baseline + motion -> senddata   [does not exist]
4 close the loop on a STILL prep, measure residual
5 only then a breathing animal
```
Requirement: a sample with correlatable structure that HOLDS STILL, and one z-stack
(~+-25 um at 1 um). No calibration grid, no brain map, no stage moves beyond the stack.

**Work in PIXELS laterally, um axially.** The corrector needs "image shift -> coefficient",
so if the gain is measured in px/coefficient the units cancel and `objectiveResolution`
never enters — which also sidesteps the measured ~6% x/y scale disagreement (the resonant
and galvo axes are not calibrated together; see [[project-autostitch-ncc-260910]]).
If physical um are ever wanted, **the stage is the ruler** (move a known distance, measure
the pixel shift) — never `objectiveResolution = 38.9`, which is one number for both axes.

## Stage 0 script (written, NOT run)
`RZ_MATLAB\ao_motion_260911\dm_motion_calib_260911.m`. Sends `baseline + c*Z2C(row,:)` for
tilt/tip/defocus, grabs, measures shift, fits a 3x3 gain matrix (off-diagonals = cross-talk;
the corrector must invert the WHOLE matrix, not use three scalars). Baseline restored by
`onCleanup` even on error/Ctrl-C. A step that would clip is **skipped, not clamped** — a
clamped mirror is not producing the Zernike you asked for and would silently bend the fit.
Still needs rewriting to register into the stack rather than measure lateral and axial
separately.

Analysis helpers tested against synthetic data with known answers. **Two findings:**
- `phase_shift` had an **off-by-one**: MATLAB indexes the zero-shift peak at (1,1), so every
  measurement was 1 px out. Fixed. Invisible in a slope, wrong in every displacement.
- **SIGN:** `phase_shift(ref,im)` returns the shift mapping `im` BACK onto `ref` — the
  CORRECTION, not the displacement. Backwards in a corrector and the loop drives the motion
  instead of cancelling it and runs away.
- Correlation-based z matching is **degenerate about focus** if the structure is a symmetric
  blur (a bead or slide could be). Real tissue is not — cells enter and leave.

Related: [[project-auto-acq-260910]], [[reference-scanimage-rig-control]],
[[project-breath-svd-pc1]].
