---
name: reference_laser_power_lut
description: "The two measured ScanImage %-to-mW calibrations (pre/post 2026-07-23), why mW not %, and how to interpolate"
metadata: 
  node_type: memory
  type: reference
  originSessionId: f3c40add-59ae-4001-987d-927c5eb434ad
  modified: 2026-08-06T06:20:08.430Z
---

ScanImage `%` is a **Pockels setpoint, not power**. Held in `laser_power_calibration.m` (repo root) — the only place these tables should live.

**Two tables, boundary 2026-07-23** (set inclusive: 7/23 itself uses the new one — unconfirmed which side of that day the re-measurement fell).

- `pre_260723` — 0.10–90 %, 2.3–392 mW
- `post_260723` — **re-measured 2026-08-05, 0–85 %, 0–371 mW.** Starts at a TRUE zero (0 mW at 0–1 %), so unlike the old table it is usable at the low end. Peaks at 371 mW near 81 % and rolls over by 85 %. This SUPERSEDED an earlier 1–50 % post-260723 table; anything computed before 2026-08-05 with `post_260723` used the old numbers.

**The rig barely lost power at the top; nearly all the change is at the bottom.** 48 %: 262→248 mW (0.95). 11 %: 19.6→15.4 mW (0.79, and **0.62 in signal**). So the LUT choice matters most for *low-power* rounds.

**The old table's low end is fake** — it reads 2.3 mW at both 0.1 % and 1 %, a meter floor / Pockels leakage. Unusable below ~3 %. The new one resolves 0.6 mW at 1 %.

**mW is nowhere near linear in %.** Local exponent m in mW ∝ %^m: ≈**1.9–2.0 around 10–20 %**, falling to ~1.3 at 50 %. Since signal ∝ mW², signal goes as ~**%^3.8** there. Using % as the denominator of a P² correction is not a scale error — it is the wrong exponent. Concretely for the 260728 rounds, the correct divisor is 13.9 / 6.6 / 2.8 (new LUT) where naive %² gives 3.6 / 2.8 / 1.9.

**The old curve ROLLS OVER at 75 % (392 mW peak)** — above that more setpoint loses power. The new table stops at 50 %, so nothing above 50 % on post-7/23 data can be quoted in mW at all.

## Interpolation: `pchip`

Measured. Above ~10 % the choice is worth <0.5 %; at 11 % linear vs pchip differ 2.1 % (4 % in signal); the worst disagreement is **12.5 % at 1.7 %**, where the curve is steeply convex and a linear chord cuts across it. Spline can overshoot where the point spacing changes (up to 2.7 % near 2 %); pchip is shape-preserving and cannot make power fall as setpoint rises.

**Table density**: leave-one-out on the 18 interior points of the new table gives median error 0.36 %, worst 8.2 % (at 3 %). Above ~5 % the spacing is the limit, not the interpolant.

**Outside the table range is extrapolation, not interpolation** — pchip extrapolates badly on a curve that physically rolls over. `laser_power_calibration.m` clamps and warns instead.
