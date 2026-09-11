---
name: project-rotated-mount-sessions
description: "Sert/0721 and Vglut2/0728 were mounted rotated 90 deg, so their stage x/y carry swapped anatomical meaning — any spatial analysis must remap them"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-07-31T18:53:58.349Z
---

**THE RULE IS TEMPORAL: everything acquired from 260721 onward is on the new mount.**
Confirmed by Runqi 2026-07-31 when `Vgat\0730` was added. Assume any NEW session is
rotated and add it to `rotated_sessions`; do not try to infer it.

**Do NOT infer the mount from where the stage coordinates land.** Tried on Vgat/0730 —
its footprint (x −1008, y +1143) sits squarely inside the old-mount cloud (ChAT −1018/139,
Sst −1041/731, Vglut2/0224 −1634/1515) and nowhere near the other rotated sessions, which
made it look un-rotated. It is rotated. The mount is a fact about the rig on the day and
coordinate overlap says nothing about it.

Three sessions were mounted **rotated 90°** relative to the older ventral-surface
sessions, so their stage axes do NOT mean what the rest of the dataset assumes:

| session | mapping |
|---|---|
| `Sert\0721`, `Vglut2\0728`, `Vgat\0730` | **−stage x = anterior (rostral)**, **stage y = lateral (+y = left)** |
| everything else | stage x = lateral, +stage y = rostral |

As a 2×2 sending `[stage x; stage y]` to `[lateral; rostral]`: `[0 1; -1 0]` for the
rotated sessions, identity otherwise. Applied in
`analysis_260727\Ventral_surface_ROI_on_cartoon_260729.m` as `AXMAP_ROT90` /
`rotated_sessions`, keyed on **session (scan dir + date), not genotype** — `Vglut2\0224`
and `Vglut2\1124` are NOT rotated, so a genotype-level flag would wrongly move the older
pFN data.

**Three things follow, and all three are easy to get wrong:**

1. The **lateral sign is irrelevant** to the map — every FOV is folded into the left
   hemisphere with `-abs()`. The only consequential content is "rostral = −stage x".
2. The **FOV footprint transposes**. Scan mirrors are fixed to the RIG, so a FOV column
   is always stage x and a row always stage y; the mount enters only at stage→anatomy.
   A 320×512 px FOV spanning 569 µm of stage x therefore spans 569 µm *rostrally* and
   356 µm laterally. `rot90(...,-1)` carries texture and extent across together. Sert
   has non-square FOVs (320×512, 512×380), so this is not academic.
3. **Corroborated by the experimenter's own filenames**: rotated sessions write the pair
   y-first (`y1900x1320`) and satisfy `motorX = -nameY`, `motorY = nameX`. Filenames are
   otherwise unreliable in these sessions (one Sert name has a flipped sign) —
   **`motorPosition` in `_meta.mat` is the only authority**.

Symptom if unfixed: Sert appears as two disconnected blobs (one at y≈+500, one at
y≈−1100). Corrected, it is a single contiguous rostro-caudal strip 400–1200 µm lateral,
y ≈ 0 → +1400, i.e. the parapyramidal/raphe region.

See [[project_ventral_summary_archive]], [[project_analysis_260727_cell_pooled_coherence]].
