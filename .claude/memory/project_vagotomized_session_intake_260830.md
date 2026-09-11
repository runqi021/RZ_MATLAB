---
name: project-vagotomized-session-intake-260830
description: "The per-session intake drill for vagotomized Sert/ChAT sessions (stage -> pair -> cell_link -> append registry -> per-cell figures), plus the site-grouping trap and the recurring unlinked-append bug"
metadata: 
  node_type: memory
  type: project
  originSessionId: 41ccab4b-4965-4c62-809e-6991e34f1bb8
  modified: 2026-08-30T09:16:27.040Z
---

The repeatable drill for taking a freshly acquired session (e.g.
`C:\Users\Admin\Desktop\260829_Sert-soma-g8s_vagotomized`) all the way to
per-cell figures. Established over Vglut2/0824, ChAT/0826, Sert/0828.

**1. Stage into the archive** `D:\Ventral_surface_summary\<Genotype>\<MMDD>\<site>\<recording>`.
Only 10 light files per recording — **no big TIFFs** (Sert/0828 = 33 MB for 3 recordings):
`breath_pc1.mat`, `breath_peak_pc1.mat`, `breath_insp_start_pc1.mat`,
`ca_spike_data.mat`, `*_ch1_dFF.mat`, `*_ch1_meta.mat`,
`*_ch1_preproc_MC_MC_AVG_for_CP.tif`, `*_..._AVG_ROIlabel.tif`,
`*_..._AVG_ROImask.tif`, `*_..._cpSAM_output.mat`.
Site folders **must be named `cell*`** — `make_cell_link_260824.m:91` enumerates them
with `dir(fullfile(archiveSession,'cell*'))`. The site folder does NOT appear in the
registry cell key (keys are `roi:<Genotype>/<MMDD>/<recording>/<roi>`), so site naming
is organisational only — it decides what gets matched against what.

**2. Pair** — `roi_pair_morph_match_260727.m` (+ `cell_cfg_260727.m`,
`roi_review3_merge_260727.m`, `roi_review4_ungrouped_260727.m`) per multi-recording site
→ `roi_match_out_<date>\roi_match_curated.mat`.

**3. `make_cell_link_260824.m`** — edit `archiveSession` / `curatedMat` / `genotype` /
`dateStr` → `<session>\cell_pooled\cell_link.mat`. Tossed ROIs must stay in with
`cell_id = NaN`; dropping them puts the curation's rejects back into every analysis.

**4. `append_0824_to_registry_260825.m`** — reusable, edit its USER block
(`sessionRel`, `scanName`, `recDate`, `linkFile`). Appends AFTER the last existing cell so
no existing figure filename goes stale; backs the registry up first and refuses to
double-append. Writes `cell_id_map_<MMDD>.csv`.

**5. `analysis_260806\per_cell_summary_260812.m`** with `selMode='active'` →
`Sert_<MMDD>_cell<NNN>.{png,pdf,_avgproj.png}` in
`D:\Ventral_surface_summary\per-cell-summary_active_260812\`.

## Registry state 2026-08-30
330 cells / 2274 obs / 98 recordings. Sert/0828 = cells **301–330**, so the next session
appends from **331**. Key prefixes: `Sert/ Sst/0 Vgat/ Vglut` = 104 LINKED cells (key is
the cell_link group); `roi:` = 226 per-ROI cells.

## BUG — Sert/0828 was appended UNLINKED (open)
All 30 Sert/0828 cells have `roi:` keys and **exactly 1 observation each**. Both
`roi_match_out_260724\` and `roi_match_out_260828\` contain a `roi_match_curated.mat`, but
`D:\Ventral_surface_summary\Sert\0828\cell_pooled\cell_link.mat` **does not exist**, so
step 4 ran with `linkFile=''`. Its three recordings are one site at three z-planes, so its
cells are almost certainly split one-per-recording right now. This is the THIRD occurrence
of the same failure — see [[project-cell-link-unapplied-curation]] for Vglut2/0224 and
1124. **Step 3 is the step that gets skipped; check for `cell_pooled\cell_link.mat`
before every append, and check `numel(CELL(i).obs) > 1` for at least some cells after.**

## Site grouping — measured, then confirmed by the matcher
Cluster recordings by `motorPosition`, never by the roi number in the folder name
(the standing rule, see [[project-ventral-summary-archive]]). Single-linkage at 150 µm
gives 11 sites for 260829's 24 recordings; two folder names lie:
- **`roi12_8x_z0_15lp_6000f` is at (−828, −1139)**, 28 µm from `roi9_6x` — the other three
  roi12 recordings sit at (−620, +910), over 2 mm away. It belongs to the roi9 site.
- **roi7 and roi14 are ONE site** — all six recordings within ~100 µm of (−1090, −830).
- `roi4` splits across TWO sites (−1083,+771) and (−1241,+988), 265 µm apart, so site
  folders cannot simply be named after roi numbers without colliding.

Multi-recording sites (the only ones needing a pair pass): {roi9×2 + roi12_8x_z0},
{roi12×3}, {roi1, roi2, roi15×2}, {roi7×2, roi14×4}, {roi6, roi11} = 16 recordings.
**Over-grouping is the safe direction** — the matcher plus the review GUI reject false
pairs, whereas splitting a site silently loses merges.

## 260829 — COMPLETED 2026-08-30, and the scripts that did it
`analysis_260829\` holds the three dated scripts, each a copy of the 0824/0825 originals
with the USER block retargeted (originals untouched):
`stage_0829_to_archive_260830.m`, `make_cell_link_260829.m`,
`append_0829_to_registry_260830.m`.

Numbers: 24 recordings, 515 ROIs, `ifSpike` keeps 89. Staged to
`D:\Ventral_surface_summary\Sert\0829\cell01..cell11`, 110 MB. cell_link = **296 cells,
108 spanning >1 recording**, 472 obs assigned / 43 tossed. Registry 330 → **366 cells,
ids 331–366** (36 ACTIVE by the append's `activeMinEv = 5`), 2746 obs, 122 recordings.
All 36 carry LINKED `Sert/0829#cNN` keys, 18 of 36 seen in >1 recording, max 6 obs —
verify exactly this after every append. Backup:
`event_latency_data_backup_pre0829_Sert-0829.mat`.

"Curated" in this workflow means `ifSpike`, the calcium_spike_gui keep/reject flag —
NOT ROI matching. Ask which one is meant.

## TRAP — `make_cell_link` dies on any recording folder with a DOT in its name
`roi_match_results.fov_name` is `fileparts`-truncated, so `roi12_5x_z-5_13.5lp_3000f_00001`
becomes `roi12_5x_z-5_13` and `roi14_4.2x_z-5_15lp_3000f_00001` becomes `roi14_4`. Those
stems match no folder, and `make_cell_link_260824.m`'s archive lookup errors with
"N recording(s) are in the curation but not under ...". Five of 260829's 24 names carry a
dot. **Fix, already in `make_cell_link_260829.m`:** rebuild the names from
`match.fov_folder` in `roi_match_results.mat` (which keeps full paths), re-appending the
`fileparts` extension, and assert every truncated name is a `startsWith` prefix of its
replacement before substituting. Never join on `fov_name` — it can also COLLIDE
(see [[project-cell-link-unapplied-curation]]).

## 6. Population select — groups are GENOTYPE-wide, which is wrong post-vagotomy
`popsel_precompute_260816.m` picks cells with `strcmp(p{1}, GROUP)` — the genotype
only — so `GROUP='Sert'` pools Sert/0721 (INTACT, 48) with Sert/0828 (30) and
Sert/0829 (36) = 114 cells ACROSS the manipulation. `popsel_decisions_Sert.csv`
(16 include / 10 exclude / 22 undecided) is the 48-cell INTACT curation; never mix it.

`analysis_260829\popsel_precompute_sessions_260831.m` selects on genotype AND a list
of session dates, writing `popsel_cache_<name>.mat` under the same naming, so
`popsel_gui_260816.m` opens it by setting its `GROUP` — no GUI change needed.
Built 2026-08-31: **`Sert_vagotomized` = Sert/0828 + 0829, 66 cells**, all with
rendered figures. `popsel_gui_260816.m`'s `GROUP` now points at it, and it writes its own
`popsel_decisions_Sert_vagotomized.csv`.

**0828 is in by RZ decision (2026-08-31), with a known caveat.** Its 30 registry cells are
11+12+7 ROIs from three z-planes (z-5/-10/-15) of ONE site with zero cross-recording
links, because its curation never became a cell_link. Planes 5 um apart repeat most
neurons, so a neuron can enter the population mean up to 3x — RZ judged this acceptable.
0829 by contrast is properly linked (18 of 36 cells span >1 recording), so **the two
sessions are not on the same footing inside this group**; say so in any figure caption.
To level them later: transcribe 0828's curated merges into `cell_merge_overrides_260814.m`
keyed on `(recName, roi)`. It is applied at render time by `apply_cell_merges_260814`
(union-find) and renumbers ONLY merged cells — a merged cell takes `min(origCells)` — so
0829's ids 331-366 and their figures stay valid. Then rebuild the cache with
`overwrite = true`.

See [[feedback-cell-ids-not-stable]], [[project-ventral-selected-cells-260816]],
[[project-roi-crossfov-matching]], [[project-analysis-260727-cell-pooled-coherence]].
