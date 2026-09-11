---
name: Chronux multitaper PSD conventions for whisker analysis
description: Standard TW, tapers, fpass, and plotting scale used for all whisker spectra in this project; how they relate to fs=90 Hz and the whisking band
type: project
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
All whisker spectral analysis in this project uses **Chronux `mtspectrumc`** with the following conventions (set at the top of `whisker_L_3w_batch.m`, `whisker_LR_batch.m`, and similar):

```matlab
params_mt.Fs     = fps;              % 90 Hz (Basler camera rate)
params_mt.tapers = [TW, 2*TW - 1];   % TW=6 → 11 tapers (good bias/variance at ~60 s recordings)
params_mt.pad    = 0;                % no FFT zero-padding
params_mt.fpass  = [5 44];           % display 5-44 Hz (must stay <= Nyquist = 45 Hz)
params_mt.err    = 0;                % no error bars from Chronux (we compute SEM across runs ourselves)
```

- **TW = 6** is the standard. Earlier scripts used TW=8; TW=6 gives narrower peaks at this recording length.
- **`fpass(2) = 44`** (not 45) — keep strictly below Nyquist to avoid edge artifacts in the returned spectrum.
- **Plot scale: linear-linear** (not log). User's preference — log-y version is commented out in `whisker_L_3w_batch.m`.
- **Per-axis PSD + sum**: `mtspectrumc` is applied to demeaned X and demeaned Y separately; per-whisker "summed PSD" = PSD_X + PSD_Y. This captures total motion power without rectification. Never rectify before FFT (see `project_whisker_signal_choices.md`).
- **Average across runs**: mean ± SEM across runs computed in MATLAB after the per-run PSDs are stacked into `[nFreq × nRuns × nWhiskers]`. SEM = std / sqrt(nRuns).

**Why:** TW=6 with 11 tapers + 5400 samples at fs=90 gives good frequency resolution for the whisking band (5-20 Hz) without over-smoothing. Linear y lets the whisking peak dominate visually, which is what the user wants.

**How to apply:** For any new whisker spectrum figure, reuse this exact `params_mt` block. Change `TW` only if recording length or expected peak width changes. Always plot lin-lin unless user explicitly asks for log.
