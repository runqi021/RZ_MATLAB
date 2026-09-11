---
name: feedback-roi-curation-port
description: "a hand curation is keyed on (folder name, ROI index) and therefore survives ANY directory reorganisation — never re-curate after moving recordings"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-07-29T05:53:53.140Z
---

A cross-FOV curation decision is really a statement about

    (recording folder NAME, ROI index)  ->  which cell

The ROI index comes from that recording's own cpSAM `maskL`, and the folder name
travels with the folder. Neither depends on where the folder sits in the tree.
**So a curation survives any amount of moving recordings between parent folders.
Re-curating after a reorganisation is never necessary.**

`analysis_260727\cell_pair_morph\roi_curation_port_260727.m` does it:
re-run the matcher on the new folder, then port. On the 2026-07-27 split of
`phys\` into `baseline\`/`roomAir\`/`reoxygenated\` it moved 447 of 447
observations and preserved 76 multi-recording cells.

**The integrity check is the important part.** What a curation does NOT survive is
re-running the SEGMENTATION, which renumbers ROIs. So the port verifies that every
matched ROI has a bit-identical centroid on both sides and REFUSES if not — porting
onto renumbered ROIs would silently attach decisions to the wrong neurons and look
completely plausible.

**Why:** Runqi worked this out unprompted ("its foldername roi number into which
group, so if i havent change the folder name and just moving them around into other
subfolders than we should always be able to port them") — it is the right mental
model and should be the default assumption.

**How to apply:** when recordings get reorganised, never propose re-curating. Run
the matcher on the new layout, port, and check the centroid assertion passed.
Groups that lose members to other condition folders keep the surviving members and
stay marked settled rather than reverting to the review queue.

**CELL IDs MUST BE INHERITED, NOT RECOMPUTED.** `cell_id` is positional — the rank
of a group within whatever subset is being analysed — so it silently renumbers on
every subset. Baseline cell 87 was full-set cell 159, and across 447 shared
observations the two numberings agreed only twice (0.4%). The port now writes
`curated.srcCellId` (the SOURCE cell id per observation) and `cell_link_260727` /
`cell_activity_summary_260727` use it when present. Ids are then SPARSE — gaps are
cells with no member in this condition folder — so every loop over cells must skip
empties. Result: one neuron carries one number in the full set and in every
condition folder, and figure filenames line up across all of them.
`cell_id_crosswalk_260727.m` maps ids between two analyses when one was NOT ported.

See [[project_analysis_260727_cell_pooled_coherence]], [[project_roi_crossfov_matching]].
