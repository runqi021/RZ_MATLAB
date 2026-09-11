---
name: reference-scanimage-rig-control
description: "D:\RZ_ScanImage_script layout: how ScanImage launches, the Motors/stage API, the adaptive-optics add-on, and auto_acq_n_shutterOFF.m (the script that acquired the maps)"
metadata:
  type: reference
---

Local copy of the imaging computer's ScanImage control, at **`D:\RZ_ScanImage_script`**.

## Launch
Just **`scanimage`** (`SI2018bR1_2018-12-19_4a9264c4fc\scanimage.m`). Assigns into base:
`hSI` (scanimage.SI, the model), `hSICtl` (scanimage.SIController, the GUI layer),
`hAOROIctrl` (AOcontrol.AOROIctrl). Config: `MDF_PXI_resonant_scanner_20210322.m` +
`PantongYao_resonant_scanner.cfg/.usr` — that cfg path appears in every TIFF header.

## Stage  (`hSI.hMotors` = scanimage.components.Motors; GUI `guis\motorControlsV5.m`)
Galil DMC-4040 on COM1, XYZ.
- **`positionDeviceUnits = [0.7815 -0.7815 0.3125]` — the Y term is NEGATIVE.** This is the
  hardware origin of the "image ROW runs along -stage y" convention the stitchers use
  ([[reference_stage_axis_convention]]). It is in the MDF, not a display choice.
- **`moveStartRelative(pos)` takes ABSOLUTE ScanImage coordinates despite its name.** It
  reads current position, overwrites non-NaN entries with `pos`, converts to motor space;
  "relative" refers to the downstream motor-space command. NaN = leave that axis alone.
  `hMotors.motorPosition = [x y z]` calls the same thing and waits.

## Adaptive optics (`+AOcontrol`, an add-on, NOT stock ScanImage)
`DMctrl` ALPAO deformable mirror serial **BAX331** via the `asdkDM` mex; `SHctrl`
Shack-Hartmann on an Andor EMCCD (cooler, EM gain); `AOdata` + `calculateWF` reconstruct
the wavefront; `AOROIctrl` is the controller -- per-ROI wavefronts, separate SYSTEM vs
SAMPLE correction files, a Z2C (Zernike->command) matrix, autoROI, and AOAqstart/AOAqend
hooked into the grab. `Aqstart`/`Aqstop`/`Framedone` are SI user-function hooks.
Top level holds sensorless-AO-by-descent and DM/SHWS calibration scripts; the measured
corrections live in `calibration-*`, `*System_aberration_correction*`, `WFS_REF`, `mshwfs`.

## `auto_acq_n_shutterOFF.m` -- THE script that acquired the maps
Chameleon on COM4 -> on, shutter open, 930 nm, GDD 12500, blocks until ready. Then a snake
over numCols x numRows at xStep/yStep um (even rows L->R, odd R->L). Origin = whatever
`motorPosition` reads at start, so the corner is framed by hand first.
- **The filename `_x###_y###` is the ACTUAL motor position read AFTER the move**
  (`pos = hMotors.motorPosition`), printed `%.0f`. NOT the commanded position. That is why
  it matches `SI.hMotors.motorPosition` in the header to within 0.49 um -- pure rounding.
- Each move is split into **10 equal ABSOLUTE sub-steps**, each `moveWaitForFinish` +
  `pause(1)` -> **>=10 s per tile of pause alone** (~13.5 min over 81 tiles).
- `z0` is captured once and resent on every move: Z is pinned to the starting focus for the
  whole grid, no drift compensation over a 2+ hour run.
- **The trailing shutter-off line is commented out DELIBERATELY** (RZ, 2026-09-10) -- do not
  "fix" it. Despite the filename the script leaves the laser on and the shutter open at the
  end. For the record the commented form would also not work as written
  (`while strcmp(acqState,'idle')` is true the moment the grid ends and nothing inside
  changes it), so if it is ever wanted it needs to be a plain sequence, not a loop.
- 260909 ChAT map = this script with numCols/numRows 9, step 400, prefix
  '260909_ChAT_g8m_shiverer'.

Other folders: `AutoMotor\` earlier tiling attempts; `auto_RZ\` laser-GUI and dye-imaging
experiments; `LaserControl-master\` Chameleon control; a 256 KB `SIController.m` copy sits
at the top level beside the install.

See [[project_autostitch_ncc_260910]], [[reference_stage_axis_convention]].
