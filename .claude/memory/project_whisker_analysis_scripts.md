---
name: Whisker analysis + label-video script suite
description: Per-whisker analysis (traces, PSDs, averages) and cross-marker label-video pipeline — what each script does and the conventions they share
type: project
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
## Analysis scripts (MATLAB, `DLC_tracking/`)

- **`whisker_LR_batch.m`** — older. vL2 vs vR1 comparison across runs, 5×6 grid, Y demeaned, lin-lin PSD TW=6. Filenames `whisker_<N>runs_{traces,spectra}.pdf` in the data folder.

- **`whisker_L_3w_batch.m`** — current primary. Whisker-centric: each whisker processed independently with its own `run_idx_<vW>` list at the top. User picks subset via `whiskers = {'vL1','vL3','vR1','vR2'}`. Grid auto-sizes per whisker (cols = min(6, nPanels)). Outputs per whisker: `traces_XY.pdf` (raw xy trajectory), `spectra.pdf` (per-run PSD_X+PSD_Y grid), `spectra_avg.pdf` (mean±SEM across runs). Log-y average is in the source but commented out. xlim/ylim are shared across all whiskers for comparability. Bandpass on x/y is optional (`bpBand = []` default); use MATLAB's built-in `bandpass()` (FIR) rather than `designfilt` IIR to avoid ringing when upper cutoff is near Nyquist.

## Label-video scripts (Python, `DLC_tracking/`, run via dlc310 conda env)

- **`whisker_label_video_test.py`** — single-video overlay. Draws + (plus) crosses at vL1 (blue BGR `(217,115,0)`) and vR1 (red BGR `(25,25,217)`), rotates canvas 150° CCW with black padding, crops, writes `.mp4` with `mp4v` codec.

- **`whisker_select_crop.py`** — interactive `cv2.selectROI` helper. Opens a middle frame rotated the same way, prints `crop_xywh = [x,y,w,h]` to paste into the other scripts. Must be run in user's terminal session (`! <cmd>`) so the ROI window shows up interactively.

- **`whisker_label_video_batch.py`** — batch version: loops over all 30 runs with the same rotation/crop/colors. ~1 minute per video on the A1000 laptop (≈30 min total).

Shared conventions:
- Colors match the analysis plots (vL1 blue = `[0 0.45 0.85]` RGB, vR1 red = `[0.85 0.10 0.10]` RGB).
- Rotation: `cv2.getRotationMatrix2D` with an expanded canvas so nothing clips, black borders.
- Output saved into subfolder `representative_whisker_label_video/` inside the data folder.
- MJPG codec caused visible line-glitching artifacts on the drawn markers → switched to `mp4v` (.mp4) which is stable.

**Why:** User iterated through many variants (vL2 vs vR1, then per-side multi-whisker, then per-whisker independent run selection). Current state is the whisker-centric pipeline because different whiskers often have different usable subsets of runs. Label videos use only two whiskers (vL1, vR1) as visual representatives.

**How to apply:** For a new analysis, edit the per-whisker `run_idx_<vW>` variables at the top of `whisker_L_3w_batch.m`, and set `whiskers` to the subset you want. For label videos on a new dataset, regenerate `crop_xywh` via `whisker_select_crop.py` once, paste into both label scripts, run.
