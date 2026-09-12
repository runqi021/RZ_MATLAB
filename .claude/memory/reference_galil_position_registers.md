---
name: reference-galil-position-registers
description: "hSI.hMotors.motorPosition on this rig is Galil TD - an OPEN-LOOP step-pulse count, not a measured position. It cannot detect backlash, lost steps or leadscrew error."
metadata: 
  node_type: memory
  type: reference
  originSessionId: c2f9c05c-e8f0-4298-9f58-57d55a9092fe
  modified: 2026-09-12T00:10:25.581Z
---

**`hSI.hMotors.motorPosition` is not a measurement.** Read path, traced in the
SI2018b source on the rig 2026-09-11:

```
hMotors.motorPosition -> hMotor(1).positionRelative -> hLSC.positionAbsoluteRaw
   -> getStepperCount()  ->  Galil "TD"
```
`+dabs\+galil\DMC4040.m:132-134` (the getter) and `:164-165` (the TD command).

For a **stepper** axis Galil `TD` is the number of step pulses the controller
GENERATED. It increments whether or not the stage actually moved. So it echoes
the command: ask for 400 um and it will report ~400 um, in both directions,
every time, even with gross backlash or a stalled motor. **Any "is the stage
accurate?" test built on `motorPosition` alone returns a table of zeros and
proves nothing.**

This also reframes an older note: "the commanded position IS the motor position
here, so filenames are a fine seed" ([[project-autostitch-ncc-260910]]) is true
but circular — the tile filename is built FROM motorPosition, so agreement
between them was never evidence about the stage.

## The two registers ScanImage does NOT use
Both are public methods on `hSI.hMotors.hMotor(1).hLSC`:
| method | Galil cmd | what it is |
|---|---|---|
| `getEncoderPosition()` | `TP` | encoder register — REAL position, **if encoders are wired** |
| `getReferencePosition()` | `RP` | commanded reference |
| `getStepperCount()` | `TD` | step pulses generated (what SI uses) |

If `TP` is live, **`TD - TP` is the direction-dependent deviation** directly.
If it is not, no register on this controller can measure backlash and the only
rulers left are the 2P image (phase-correlate a frame before/after a move) or an
external dial indicator. **Probe before trusting either branch.**

## Units — THE MDF LITERAL AND THE LIVE PROPERTY DIFFER BY 1e6
MDF `MDF_PXI_resonant_scanner_20210322.m:61` reads
`motors(1).positionDeviceUnits = [0.7815 -0.7815 0.3125]` with the comment
"in meters", but **`hLSC.positionDeviceUnits` returns `[7.815e-07 -7.815e-07
3.125e-07]`** — SI stores that literal as METRES. So:

```matlab
devUm = hLSC.positionDeviceUnits * 1e6;   % 0.7815, -0.7815, 0.3125 um/count
```

**Using the property directly as um/count scales every register reading by
1e-6.** It cost a completely bogus probe run on 2026-09-11: a real 25.6-count
move printed as `+0.000` and the script concluded "the encoder is dead" when
nothing had been shown either way. Sanity-check any device-unit value against a
plausible 0.01-100 um/count range, and **report raw COUNTS as the primary
number** — a count is unambiguous, a converted micron is not.

**y is NEGATIVE** (matches [[reference-stage-axis-convention]]) and **z is finer
than x/y** (0.3125 vs 0.7815 um/count).

## Quantisation follows from the quantum, and it is NOT backlash
Only whole counts can be commanded, so a step that is not an integer number of
them is rounded — a fixed, **direction-INDEPENDENT** error that cancels in
`mean(err down) - mean(err up)` but dominates the individual columns.
Computed for the candidate steps: on **z** (0.3125) 0.5 um is +25%, 1 and 2 um
are -6.25%, while **5 um = exactly 16 counts** — which is the mechanical reason
5 um is the best small z step. Exact alternatives below it are 1.25 (4 counts)
and 2.5 (8 counts), NOT 1 or 2. On **x/y** (0.7815) anything under ~10 um is
quantisation-dominated: 0.5 um is +56%, 1 um is -22%, 5 um -6.2%, while
50/100/200/500 are all +0.03%. 0.1 um is below one count on every axis. `invertDim = '+++'`, `moveCompleteDelay = 0`.
`resolutionRaw` returns NaN — do not read resolution off the LSC.
`hMotors.moveStartRelative` takes coordinates in the SAME frame as
`motorPosition` despite the name; that is the call `auto_acq_260910` uses and it
is proven on the rig.

## The script
`matlab-wroking\stage-cal\stage_backlash_260911.m` (repo: `stage_cal_260911\`).
`MODE='probe'` does ONE small move and reports whether TP tracks it — run that
first. `MODE='full'` sweeps step sizes per axis, approaching each target from
both directions, and reports backlash = mean(err down) - mean(err up).
Safety: `DRY_RUN=true` default, `onCleanup` returns the stage to where the user
left it, and every target is checked against a max-distance-from-centre budget
(user-stated 2026-09-11: **x/y +-3000 um about a centre at start+[500 500],
z +-300 um** — widened from +-1500 late that day; worst planned case uses 67% of
it) BEFORE anything moves — refused, never clipped. No acquisition: laser, shutter, scanner and DM are untouched.
**Written 2026-09-11, NOT YET RUN.**

Caveat this puts on other work: the AO plan's "the stage is the ruler you
already trust" ([[project-ao-motion-correction-260911]]) assumes commanded stage
moves land where asked. That is exactly what is now in question — settle the
probe before using stage moves to calibrate anything.
