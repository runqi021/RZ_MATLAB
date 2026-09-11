# breath_svd — getting a breathing trace and its events out of the Basler video

Everything upstream of `analysis_260727`. This folder turns a behaviour `.avi`
into the three files the rest of the chain consumes, and nothing here knows
about calcium.

```
<recording folder>\
    breath_pc1.mat             the breathing trace          <- produced here
    breath_peak_pc1.mat        inspiratory peaks   TRIGGER  <- produced here
    breath_insp_start_pc1.mat  inspiration onsets  TRIGGER  <- produced here
```

`analysis_preflight_260727.m` checks for exactly these three. A recording joins
the breath × Ca analysis only if it has both triggers plus `ca_spike_data`.

---

## Which script makes the per-video breath trace?

**`breath_svd_pc1.m`.** That is the current, per-video route, and every dataset
analysed so far was built with it. It is what `analysis_preflight_260727.m`
names when a recording is missing its trace.

```matlab
breath_svd_pc1          % drag an ROI per video -> breath_pc1.mat + .png
breathing_peak_gui_pc1  % Browse Master Folder -> breath_peak_pc1.mat
breathing_trough_gui_pc1%                      -> breath_insp_start_pc1.mat
```

fps is never typed in or taken from the AVI header — the breath cam is
2P-frame-triggered, so `detect_session_fps` reads the imaging rate out of the
ScanImage metadata. The header fps on this rig is wrong (300 in a 30 Hz file on
260728), so anything that trusts it is wrong by 10×.

## Then what is `breath_fixedmetric_*`?

A **second, newer front end** (built 2026-07-24) that writes the *same*
`breath_pc1.mat`, so the peak/trough GUIs and everything downstream are
unchanged. It exists for one reason: **per-video PC1 is not comparable between
videos.** Each video gets its own spatial basis, its own arbitrary sign and its
own arbitrary scale, so "breath amplitude" means something different in every
file. Z-scoring does not fix that — it fixes the scale and destroys the
amplitude.

| | `breath_svd_pc1` | `breath_fixedmetric_*` |
|---|---|---|
| unit of work | ONE video | ONE session of many runs |
| ROI | one per video | **one per session** |
| basis | re-fit per video | frozen once on the baseline runs |
| output | PC1 (arbitrary units) | `disp` (PIXELS), `fb`, `pv` |
| amplitude comparable across runs | no | yes |
| right when | videos differ, or you only need timing | many runs of the SAME view |

`pv` inside the fixed-metric output **is** per-video PC1 — kept only so the
other two can be benchmarked against it.

Use `breath_svd_pc1` when each recording is a different view (e.g. one video per
imaging FOV) or you only need breath timing. Use `breath_fixedmetric_*` when you
are comparing breath *depth* across runs, animals or hours.

### Running the fixed-metric route

```matlab
breath_fixedmetric_gui('D:\<session>\cam1')   % ROI -> extract -> analyze -> bridge
```

One call does all four steps. It wants a folder of **run subfolders, each with
one `.avi` and a plain `timestamps.csv`** — i.e. the acquisition folder as the
Basler GUI writes it, *not* a `phys\` folder after pairing (the pairing step
renames `timestamps.csv` to `<avi>_timestamps.csv`, which the extractor will not
find; it then silently falls back to `FPS_FALLBACK = 60`).

`STRIDE` is tuned for the 60 fps chest rig. On a 30 fps 2P-triggered camera set
`STRIDE = 1`, or the Nyquist assertion aborts — correctly, since 15 Hz cannot
carry an 8 Hz band.

Detect on `fb`, read amplitude off `disp`. The GUIs save event *indices*, and
`breath_fixedmetric_analyze` applies those same indices to every metric, so the
choice of detection trace costs nothing in amplitude.

`breath_fixedmetric_dispaxis.m` is the diagnostic for the measurement itself:
the `(dx, dy)` cloud with the frozen motion axis through it, on-axis vs
orthogonal spectra, a sub-pixel magnitude histogram, and per-run axis drift. It
is a pure re-plot of `breath_fixedmetric.mat` and decodes no video.

## Two things that look like results and are not

- **A flat tail is saturation, not a plateau.** Wherever the sensor clips, the
  brightest-pixel metrics stop measuring the sample and report the ADC ceiling.
- **A per-run axis angle only means something if that run has an axis.** When
  the displacement cloud is isotropic the top singular vector is a random
  direction; a low `|cos|` there is not evidence the axis moved.

## Folders

| | |
|---|---|
| `legacy_dlc\` | `breathing_peak_gui.m` / `breathing_trough_gui.m` — the DLC-dot predecessors of the `_pc1` GUIs. Still correct for DLC-tracked sessions; the `_pc1` versions are the same workflow reading PC1 instead of dots. |
| `diagnostics\` | Method-development scripts from June 2026 (onset-method comparison, PC1-vs-DLC, doublet and diff-image checks). Nothing calls them; kept because they document *why* the current choices were made. |

## Paths

These files sit **two levels** below the repo root, so `repoRoot` needs two
`fileparts` — one resolves to `analysis_260727` and `detect_session_fps`,
Chronux and `+helper` are not found. Same correction the rest of
`analysis_260727` carries. `orofacial_crop_extract.py` deliberately stayed at
the repo root: it is the shared cropper, used by the orofacial analysis too.
