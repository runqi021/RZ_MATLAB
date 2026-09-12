---
name: feedback-measure-raw-not-converted
description: "Method rules learned the hard way on the 2026-09-11 stage characterisation: record RAW device units, never trust a readback that echoes the command, always carry a null control, and design metrics as differences."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c2f9c05c-e8f0-4298-9f58-57d55a9092fe
  modified: 2026-09-12T00:36:27.004Z
---

Six rules, every one of them paid for during the 2026-09-11 stage-backlash work
([[project-stage-backlash-260911]]). They apply to any instrument measurement,
not just the stage.

**1. Record RAW device units. Convert only at analysis time.**
Two separate unit bugs hit the same script. `positionDeviceUnits` turned out to
be in METRES (a real 25-count move printed as `+0.000` and produced a confident
"the encoder is dead" verdict that was pure fiction). Then the encoder quantum
turned out not to be the stepper quantum — 0.1 vs 0.3125 um/count — so every um
figure was 3.13x too large. **Neither needed a re-run, purely because the raw
counts were in the .mat.** Save `tdCnt`/`tpCnt`-style raw fields always.

**2. A readback that always agrees with the command is not a measurement.**
`motorPosition` is Galil TD, the count of pulses emitted. It matched the
expected counts/um to four decimals on all three axes — which looked like
excellent accuracy and was actually just the controller repeating itself. Ask of
any sensor: *could this value ever disagree with what I asked for?* If not, it
measures nothing. See [[reference-galil-position-registers]].

**3. Carry a null control that MUST come out zero.**
The analysis computed backlash on TP (the encoder) **and** on TD (the command).
TD came out 0.000 everywhere, as it must, since both directions visit identical
targets. That single column is what rules out "the command itself was
asymmetric" and makes the TP result mechanical rather than an artifact.

**4. Design the metric as a DIFFERENCE where possible.**
Backlash = mean(down) - mean(up) over the same target set. This survived BOTH
unit bugs and a huge unknown frame offset (+8874 um between the encoder origin
and SI's soft zero), because offsets and scale factors cancel in a difference.
The absolute `err UP` / `err DOWN` columns were garbage the whole time; the
difference of them was correct throughout.

**5. Do not average across regimes.**
The first backlash summary averaged all step sizes together, including ones
below the device resolution, and reported "varies with step size". Quoting the
plateau (steps >= 5 um) instead gave -2.05 +- 0.22 counts, flat. Sub-resolution
points are not weak data, they are a different measurement.

**6. Detrend before quoting scatter.**
`std` of the readings within a condition was reported as "repeatability". It
came out as exactly targetSpread x counts/um — it was measuring how far apart
the targets were, and said nothing about the stage. Fit the line, quote the RMS
residual about it.

**And one on identifiability:** from (commanded, readback) alone, exactly ONE
quantity is unrecoverable — the single global scale factor. Everything else
(direction asymmetry, step-size dependence, position dependence, repeatability)
is measurable. Do not dismiss a dataset as circular when only one number is.
Here the fitted encoder quantum landed on 0.09987 um/count, within 0.13% of a
round 0.1 — itself evidence that the scale is right, since a real scale error
would have produced an odd number.
