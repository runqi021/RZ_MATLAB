---
name: project-stage-backlash-260911
description: "MEASURED on the rig 2026-09-11: z has a 0.1 um encoder and ~0.205 um backlash on reversal, flat across 5-100 um steps. x and y have NO encoder readout at all, so their backlash is unmeasured."
metadata: 
  node_type: memory
  type: project
  originSessionId: c2f9c05c-e8f0-4298-9f58-57d55a9092fe
  modified: 2026-09-12T00:36:33.484Z
---

First real characterisation of the Galil stage, 2026-09-11. Scripts:
`matlab-wroking\stage-cal\stage_backlash_260911.m` (acquire, 921 moves, no
acquisition) + `stage_backlash_analyze_260911.m` (offline, re-derives everything
from raw counts). Repo: `stage_cal_260911\`. Data:
`stage_backlash_260911_1701.mat`.

## RESULT: z
- **Encoder resolution 0.09987 um/count (~0.1 um)**, fitted from `tpCnt` vs
  commanded target over 258 points. **This is NOT the stepper quantum** —
  `positionDeviceUnits` gives 0.3125 um/count for z, 3.13x coarser. Using the
  stepper quantum on TP inflates every um figure by that factor.
- **Backlash = 2.05 +- 0.22 counts = 0.205 +- 0.022 um** over steps 5, 10, 20,
  50, 100 um. Five independent step sizes agreeing to +-0.22 counts is the real
  evidence; within a single condition the residual scatter is 2-5 counts.
- Below 5 um the measured backlash falls off (2 um and 1 um both -1.13, 0.5 um
  -0.73, 0.1 um 0.00) — the reversal is not cleanly resolved at a few counts.
  Quote the plateau, not the mean over all steps.
- **TD control = 0.000 counts everywhere.** TD is the commanded step count and
  up/down visit the same targets, so a non-zero value would have meant the
  command itself was direction-dependent. It is clean, so the TP effect is
  mechanical.

## The step SIZE is not direction-dependent (z) — only the reversal is
Separate question from backlash, and separately answered from the same data
(`incr_for` in the analyse script: consecutive TP increments within each
direction run, first increment after a reversal dropped because it carries the
backlash). Measured |step| up minus down, in um, by commanded step:
1:+0.044  2:+0.044  5:+0.144  10:-0.444  20:+0.089  50:-0.333  100:+1.099.
**The sign flips** and every value is a few encoder counts, inside the 2-5 count
per-condition scatter — i.e. noise. Compare the backlash over the same steps:
-2.00 -2.40 -2.00 -2.07 -1.78 counts, one sign, one magnitude. That is the
signal.

**So the direction effect does NOT accumulate along a travel.** In a snake
raster backlash shifts a whole row by a constant; an asymmetric step size would
have sheared it progressively. It is the former.

Also visible in that table, and NOT a direction effect: small steps undershoot
in BOTH directions — 0.5 um commanded measures 0.37, 1 um measures ~0.84, 2 um
~1.74, converging to within ~1% by 20 um. That is the 0.3125 um stepper quantum
(see the quantisation note in [[reference-galil-position-registers]]).

**Practical consequence: keep z-stacks UNIDIRECTIONAL.** 0.205 um is 20% of a
1 um z step and 4% of a 5 um one, but only on a direction reversal. Monotonic
stacks never pay it.

## RESULT: x and y — NOT MEASURED, and not measurable this way
`TP` never changes on x or y: the scale fit returns -0.0000 counts/um and every
backlash row is exactly 0.00. **There is no encoder readout on those axes.**
The zeros are absence of measurement, not absence of backlash — exactly the
failure mode [[reference-galil-position-registers]] warns about. Measuring x/y
backlash needs an external ruler: the 2P image (phase-correlate a frame before
and after a reversal) or a dial indicator. Re-running this script will never
produce an x/y number.

## OPEN: backlashCompensation
`Motors.moveStartRelative` (`Motors.m:725-746`) has **built-in backlash
compensation** — if `hSI.hMotors.backlashCompensation` is non-zero it overshoots
and comes back on every direction reversal. Both this script and
`auto_acq_260910` use that call. **Not yet checked.** If it is non-zero the
0.205 um is the RESIDUAL after SI's own compensation, not the raw stage.

## Two unit traps this cost
1. `positionDeviceUnits` is in **metres** — see
   [[reference-galil-position-registers]]. Cost one completely bogus probe.
2. The **encoder quantum != the stepper quantum** (0.1 vs 0.3125 um on z). The
   acquiring script's um columns are all wrong by 3.13x; its raw `tdCnt`/
   `tpCnt`/`rpCnt` are correct, which is the only reason no re-run was needed.
**Always record raw counts and convert at analysis time.**

Also flawed in the first analysis: reporting `std` of `tpCnt` within a condition
as "repeatability". It came out as exactly targetSpread x counts/um — it was
measuring how far apart the targets were. Detrend first (`pooled_resid`).

## The hardware (identified 2026-09-11)
Motors are **API Gettys / Rapidsyn 23D-6102BN**, a bare NEMA 23 stepper —
5 VDC, 1 A, 1/4" shaft, **no integral encoder in any listing**. The SAME model
appears on NEAT X-Y stages (e.g. XYR-8080) and on NEAT Z-stage blocks, so x/y
and z almost certainly run identical motors. **Therefore z's 0.1 um readout is
not in the motor** — it is a separate encoder bolted to the z stage.

The stage is **NEAT** (New England Affiliated Technologies, Lawrence MA;
acquired by Kollmorgen 1999, now Dover Motion). Dover's FAQ confirms
**optional rotary or linear encoders** on these stages, and encoder-equipped
NEAT XY stages exist on the used market. So x/y encoders were an ORDER OPTION —
this rig may or may not have them fitted.

**Fastest way to settle it, better than any script:** look at how z's encoder is
physically implemented (a housing with its own cable on the motor's rear =
rotary; a scale strip + readhead along the stage body = linear), then look for
the same thing on x and y. Cable count is the giveaway — a bare stepper has ONE
cable, an encoder adds a second connector. If x/y have the hardware, this is a
cabling / Galil-configuration problem, not a missing part.
`stage-cal\galil_encoder_probe_260911.m` (read-only: ID, TP, TD, RP, MG of
internals, all four axes A-D since the driver only ever queries A-C) checks the
electrical side; its hand-push test — move the stage by hand, see if TP changes
with the motor idle — is the decisive electrical check. **Not yet run.**

Step angle is UNCONFIRMED. The 0.78125 um/count derivation assumes the NEMA 23
standard 1.8 deg/step (200 steps/rev); no datasheet found states it.

## Still open after all this
- `hSI.hMotors.backlashCompensation` — never checked. One line, no motion.
- whether x/y have encoders fitted (see above).
- **x/y step error and backlash are entirely unmeasured.** x is the axis the
  snake raster in `auto_acq_260910` actually reverses on, over ~3.6 mm. Only the
  image can measure it — phase-correlate a frame before and after a reversal.
- position dependence / leadscrew error: untested. z's travel limit (+-300 um)
  is far too short to see a 1-2 mm leadscrew period, and x/y have no ruler.
- The stitcher's +6-8% scale disagreement is **NOT** explained by this work. The
  MDF's 0.7815-vs-0.78125 rounding is only 0.032%. Still open.

Related: [[reference-galil-position-registers]],
[[feedback-measure-raw-not-converted]], [[project-auto-acq-260910]],
[[project-autostitch-ncc-260910]], [[project-ao-motion-correction-260911]].
