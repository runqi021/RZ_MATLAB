---
name: Data folder structure and naming conventions
description: How calcium imaging (phys) and behavior video (falling_edge) folders are organized, named, and related
type: project
---

## Two-folder acquisition layout

Each experiment session produces two top-level folders on the Desktop:

1. **Calcium folder**: `YYMMDD_<line>_<indicator>/phys/`
   - Contains FOV subfolders, each named to match its raw TIF
   - Sometimes has intermediate category subfolders (e.g., `phys/IO/`, `phys/maybe_breathing/`)
   - FOV folder contents: raw `.tif`, `_meta.mat`, `_preproc.tif`, `_MC.tif`, `_cpSAM_output.mat`, `_dFF.mat`, etc.

2. **Behavior folder**: `YYMMDD_<line>_<indicator>_falling_edge/`
   - Root level: original Basler AVI files (one per recording)
   - AVI naming: `Basler_acA1300-200um__<serial>__YYYYMMDD_HHMMSSmmm.avi`
   - DLC project subfolder: `<project>-<scorer>-<date>/videos/` contains DLC copies of AVIs + output CSVs/H5/pickle

## Pairing behavior to calcium

- Both acquisitions are hardware-triggered simultaneously (same start time)
- Always more calcium recordings than behavior (some calcium FOVs have no behavior)
- Sometimes behavior recording runs longer (forgot to stop) -- trim to match calcium frame count
- Clock offset between the two computers: typically 3-18 sec, up to ~80 sec for first recording of session

**Why:** This is foundational context for `pair_behavior_to_phys.m` and any future cross-modal analysis scripts.
**How to apply:** When writing scripts that need both calcium and behavior data, expect this folder layout. Use `pair_behavior_to_phys.m` to automate the pairing.
