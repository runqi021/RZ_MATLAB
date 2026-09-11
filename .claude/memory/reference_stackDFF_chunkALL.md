---
name: stackDFF_chunkALL_260413.m — pooled stacked dFF chunks
description: Standalone script in Breath_summary_timeNphase_260331/ that stacks sig-ROI dFF chunks across multiple inputPaths with xcorr alignment + parallelogram shear plot
type: reference
originSessionId: 0ec053e8-fba2-425a-9aa8-58264c9cddb3
---
# `stackDFF_chunkALL_260413.m`

Path: `Breath_summary_timeNphase_260331/stackDFF_chunkALL_260413.m`

Built 2026-04-13 as a **new standalone script** (not as edits to `breath_combined_summary_260329.m`).

## What it does
Pools the significant ROIs across one or more experiment folders, splits each ROI's full dFF into fixed-duration chunks, optionally aligns chunks to a common breath phase via iterative xcorr template-matching, and renders a sheared "parallelogram" stacked-trace figure.

## Inputs (already produced upstream)
- `<inputPath>/breath_combined_summary/fov_map.mat` — `fov_map{i,1}` = sessName, `fov_map{i,2}` = vector of sig ROI IDs (per coherence). Produced by `breath_combined_summary_260329.m`.
- Per-session: `*cpSAM_output.mat` (re-runs `helper.dFF_RZ` to get dFF), `ca_spike_data.mat` (uses `roi_spikes(roi_id).spike_idx`).

## Output
- `<inputPaths(1)>/breath_combined_summary/stacked_dFF_chunks.pdf`

## Key params (all at top of script)
| Param | Default | Meaning |
|---|---|---|
| `inputPaths` | string array | Pool sig ROIs from all listed folders; FOV index is **global** (counts across inputPaths in the listed order). |
| `chunk_sec` | 32 | Chunk duration (s). |
| `chunk_min_n` | 3 | Skip FOV if `floor(T/L) < chunk_min_n` (the "mod >2" rule). |
| `align_xcorr` | true | Iterative xcorr alignment to rolling template. Uses `circshift` (wrap-around — appropriate for periodic breathing). |
| `max_shift_sec` | 5 | Max allowable shift per chunk (s). |
| `n_align_passes` | 3 | Template→align→retemplate iterations. |
| `dFF_scale` | 0.3 | Scale-bar height in ΔF/F (label + bar length). |
| `spacingFrac` | 0.15 | Row spacing fraction. Lower → more invasion of neighboring rows. |
| `shear_from_vert_deg` | 30 | Parallelogram lean. **Geometrically too steep for >~30 rows** — auto-clamps; reduce for many rows. |
| `ax_w_cm` | 36 | Axes width (cm). Wider axes can fit more shear without clamping. |
| `ax_h_cm_per_row` | 0.35 | Axes height per row (cm). |
| `exclude_chunks` | [] (Nx3) | Each row `[global_FOV_idx, ROI_id, chunk_idx]` to drop. |

## Key implementation notes
- **No-spike filter granularity = chunk** (drops only the chunk, not the whole ROI).
- **Global FOV index** counts across `inputPaths` in the listed order (so adding a folder before another shifts all downstream FOV indices — keep the exclusion list in sync).
- **Chunks resampled to common grid** at `fps_img` (default 30) so different per-session fps still align.
- **xcorr alignment** uses zero-mean signals, max-correlation lag clamped to `±max_shift_sec`. Iterates `n_align_passes` (template = mean of currently-aligned stack).
- **circshift** wraps the cycle continuation from the other end (no NaN edges, no flat-line artifacts that NaN-padding caused).
- **Parallelogram shear** is computed in pixel/cm units after fixing axes size (`ax_w_cm × ax_h_cm`); per-row x-shift = `gap_cm × tan(shear_from_vert_deg)`. If the resulting total shear exceeds `ax_w_cm`, the script clamps to 80 % of axes width and prints a warning.
- **Custom plotting** (does NOT call `helper.stackDFF`) so the per-row x-offset can be applied directly. Y-tick labels are at the offset positions of selected rows (every 5th + last).

## Tested folders
- `D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing` (7 FOVs)
- `D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing` (16 FOVs)

Combined run: 23 FOVs → 11 sig ROIs → 45 chunks at 32 s.

## Don't touch
This script does **not** modify `breath_combined_summary_260329.m`. Per [no-inplace-edits feedback](feedback_no_inplace_edits.md), all new analyses go into separate files like this one.
