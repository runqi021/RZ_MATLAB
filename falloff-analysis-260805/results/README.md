# falloff analysis — archived outputs (copied 2026-08-18)

Copies of every FIGURE and every FITTED DATAPOINT of the optical-penetration
project, organised by experimental folder, so the analysis survives the raw
TIFFs being moved off this machine.

Source folders at the time of copying:

| archive folder | copied from |
|---|---|
| `251104_wt_fitc`           | `C:\fall-off\251104_wt_fitc` |
| `260728_vglut2_vessel`     | `C:\fall-off\260728_vglut2_vessel` |
| `260114_shiver_vessel`     | `C:\fall-off\260114_shiver_vessel` |
| `260804_shiver_dbh_vessel` | `C:\fall-off\260804_shiver_dbh_vessel` |
| `260806_sst_vessel`        | `C:\fall-off\260806_sst_vessel` |
| `260730_shiver_vgat_g8m`   | `C:\260730_shiver_vgat_g8m` |
| `260810_vglut2-soma-g8s\roi1_stack_z-5-70_30f_00001` | same, nested as in the original |

Root: `falloff_summary_scatter.png/.pdf/.mat` — the cross-session figure.

## What is here

- `falloff_<tag>.mat`        the FIT: `Lp` (lambda), `r2p`, `FIT_DEPTH`, `POWER_LUT`.
                             This is the datapoint the scatter plots.
- `falloff_<tag>.png/.pdf`   4-row intensity-vs-depth figure
- `falloff_proj_<tag>.*`     depth-corrected MAX projections
- `falloff_gcamp_*`          260730 only (the standalone ramped-GCaMP script)
- `spfft_<tag>.*`            spatial Fourier figures
- `vessel_hist_raw_*`        260728 only
- `*_metrics.mat`            per-slice top/bot + histograms. KEPT ON PURPOSE:
                             `falloff_260804.m` loads this instead of reading the
                             TIFFs, so lambda can be REFITTED over a different
                             depth window after the raw data is gone.
- `*_cache.mat`              same idea for the spatial-Fourier script.

## What is NOT here

`*_vol.mat` — the 1x1x1 um isotropic volumes, ~2.6 GB in total. Rebuildable only
from the raw TIFFs. Consequence: once the raw data moves, the PROJECTION FIGURES
here are final — they can be re-styled only by re-reading the raw stacks. The
falloff fits are not affected, because `*_metrics.mat` is kept.

## Re-running the scatter from this archive

`falloff_summary_scatter.mat` holds `T` (every FOV: session, group, label type,
lambda, R^2, fit range) and `SESSIONS`, so the numbers stand alone.

To re-run `falloff_summary_scatter_260805('<this folder>')` after the raw folders
move, note that two rows of its `SESSIONS` table are ABSOLUTE paths
(`C:\260730_...`, `C:\260810_...`) and therefore ignore the ROOT argument. Point
them at this archive, or make them plain folder names, before re-running.

## Values as archived

| session | FOV | lambda (um) | R^2 | fit range |
|---|---|---|---|---|
| 251104_wt_fitc            | col01_row02 | 24.7 | 0.984 | 10-150 |
| 251104_wt_fitc            | col04_row02 | 24.5 | 0.965 | 10-150 |
| 260728_vglut2_vessel      | roi1        | 24.6 | 0.989 | 10-150 |
| 260806_sst_vessel         | roi1        | 25.5 | 0.994 | 10-150 |
| 260810_vglut2-soma-g8s    | roi1 (GFP)  | 20.9 | 0.997 | 10-65  |
| 260114_shiver_vessel      | tile20      | 80.7 | 0.975 | 10-500 |
| 260114_shiver_vessel      | tile23      | 78.4 | 0.988 | 10-500 |
| 260114_shiver_vessel      | tile27      | 86.7 | 0.965 | 10-500 |
| 260804_shiver_dbh_vessel  | roi1        | 63.0 | 0.946 | 10-500 |
| 260804_shiver_dbh_vessel  | roi2        | 61.4 | 0.966 | 10-500 |
| 260730_shiver_vgat_g8m    | roi1 (GFP)  | 78.0 | 0.996 | 10-500 |

wildtype  n=5  lambda = 24.0 +- 1.8 um   (EAL = 2*lambda = 48.1 um)
shiverer  n=6  lambda = 74.7 +- 10.2 um  (EAL = 2*lambda = 149.4 um)
