# Agent Memory: RZ_MATLAB Calcium Imaging Pipeline

## GUI Conventions (from breathing_peak_gui.m)
- All state in struct `S`, factory defaults in struct `DEFAULTS`
- `nCtrlRows` must match actual `r=r+1` count (enforced by assert)
- `Interpreter='none'` on all axis titles/labels
- `disableDefaultInteractivity` / `enableDefaultInteractivity` for click modes
- `strrep(sessName,'_',' ')` in figure titles
- In-place XData/YData updates via pre-allocated scatter/line handles
- `setpref`/`getpref` for persistent params between sessions
- Session scanning: `dir(fullfile(master,'**','pattern'))` then unique folders
- Save pattern: per-session .mat + master .mat in master folder

## Key File Structures
- `*_dFF.mat` contains: `dFF` [TxN], `dFFout` (struct with .t_dff, .spikes_oasis), `F_raw`, `F`, `params`
- `*_cpSAM_output.mat` contains: `F` [TxN], `t`, `maskL` [HxW uint16 label image], `avgPath`, `maskPath`, `opt`
- `maskL`: integer labels, 0=background, regionprops for centroids, bwboundaries for contours

## Created Files
- `calcium_spike_gui.m` — calcium spike detection GUI with per-ROI inspection (2026-03-01, updated 2026-03-02)
  - Paired files: `ca_spike_data.mat` (per session), `ca_spike_master.mat` (master folder)
  - Fields: `roi_spikes(n).spike_idx/.spike_t/.spike_amp/.spike_train/.n_spikes`, `spike_params`
  - v2 changes: removed video panel/TIF loading, removed contrast sliders, added Raw F panel (axRawF),
    avg projection loaded from SAM file (meanImg/mean_img/Ims/meanIm/avgImg fields), auto-contrast from ROI mask pixels,
    yellow ROI boundary, 3-row rightGL layout (0.22x/0.42x/0.36x), scanForDFFSessions no longer returns .tifFile
