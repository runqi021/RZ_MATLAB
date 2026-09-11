---
name: project_subvolume_render
description: "map_subvolume_render_260530.m — oblique-MIP \"volume render\" of a sub-volume cropped from map avgz tiles"
metadata: 
  node_type: memory
  type: project
  originSessionId: 939960db-c16e-4613-8c7e-312ac4505fd5
---

`map_stitch_GUI/map_subvolume_render_260530.m` — renders a 3D oblique-MIP "volume render" (vasculature look, red wireframe box + "X×Y×Z µm³" label) of a sub-volume cropped around a stage coordinate. Built to mimic a paper figure. Related to [[project_manual_stitch_gui]].

**Pipeline:** assembles the box directly from per-tile `chN/avgz` multi-page stacks using stage coords in the filenames (NO prior stitch); overlapping tiles MAX-blended → oblique orthographic MIP by ray-casting in real microns (interp3); box corners projected through the same transform so the wireframe always aligns.

**Key conventions:**
- `parse_xy` reads BOTH filename schemes: explicit `_x-1500_y-499` (live_251104 dataset) AND `tileNN_<X>_<Y>` (260114 dataset).
- `S.bothSides=true` renders both `+|X|` and `-|X|` (same Y/box), saving each PNG.
- `S.medialLeft=true` orients so MEDIAL (midline, X→0) is on image-LEFT, LATERAL on right; the -X hemisphere is mirrored via M=diag([-1 1 1]) on the screen basis. +X is already correct.
- `S.interp` ('cubic' default) used for BOTH tile-assembly interp2 and ray-cast interp3. makima = good alt; spline rings; Frangi (`S.vesselness='fibermetric'`, 3D, resamples Z to isotropic first) was tried and REJECTED by user as "looks fake".
- Voxel size: both datasets are 1.4815 µm/px metadata (user calibrated **1.4648**), 10 µm Z. 260114 = 76 slices (750µm); live_251104 = 23 slices (220µm). Auto-detect via `read_si_meta` if umPerPx/dz set to [].
- `S.xSign=+1, S.ySign=-1` (stage-Y opposite image rows, matches GUI flipY).

**Caveat:** no seam refinement/intensity matching — pure stage-coordinate placement + MAX blend. Fine because stage coords match motor metadata to sub-micron.
