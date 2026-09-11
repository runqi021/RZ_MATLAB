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
- `*_cpSAM_output.mat` contains: `F` [TxN], `t` [Tx1], `maskL` [HxW uint16 label image], `avgPath`/`ovPath1`/`ovPath2` (string PATHS, can be empty ""), `maskPath` (char PATH), `opt` (struct). NO centroids/areas stored — derive from maskL.
- `maskL`: integer labels, 0=background, regionprops(maskL,'Centroid','Area') for centroids ([x y] px) + area (px); bwperim/bwboundaries for contours. N labels == size(F,2). Verified maskL=[512x256] HxW.
- Numeric mean/avg-projection image is NOT in the .mat — it is the 16-bit `*_AVG_for_CP.tif` (same HxW as maskL), found via dir('*_AVG_for_CP.tif'). Downstream reads it with imread (e.g. breathing_population_manifold.m, calcium_spike_gui.m).
- `*_meta.mat` (saved by Batch_dffQC_260325.m + detect_session_fps.m) fields: fps, scanFrameRate_raw, zoomFactor, pixelSize_um (=1.7778/zoom um/px), laserPower_pct, channelSave, channelInputRanges, motorPosition, framesPerSlice, numSlices, source_tif. source='tif_imfinfo'|'qc' OR has source_tif => trusted.
  - VERIFIED: `motorPosition` = 3-elem [X Y Z] stage microns (e.g. [-987.034 290.718 -40.3125]). This IS the absolute FOV position. Processed TIFFs carry NO metadata — _meta.mat / raw TIFF are the only source of truth.
- Recording FOLDER NAMES often encode stage coords + zoom but INCONSISTENTLY: newer `roiN_4x_x-1120y600z-5_...`, older `-1050-760-z5_5x_...` or `600-1000_z-5_...`, and some have NO coords (`ROI3_5x_3000f_13lp`). Trust _meta.mat motorPosition, not folder name.

## Breath-tuning "hit" labels (breath_combined_summary_260329.m)
- Per-ROI significance: `is_sig = r_k > confC_val` where r_k = Chronux coherence magnitude in breath band (spike-train or waveform per `coherence_use_spike`), confC_val from coherencyc (alpha_sig=0.01). th_k = phase angle (rad).
- `fov_map.mat` (in <inputPath>/breath_combined_summary/): cell [nFOV x 2], col1=sessName (char/string, matches recording folder name), col2=sig ROI ID vector. This is the closest thing to a registration table (per-session, NOT cross-session).
- `coherence_data.mat` (same folder, newer runs only): struct `coherence_data` with all_th_k_wave/spk, all_r_k_*, all_rlo/rhi/dphi_*, all_is_sig_wave/spk, all_roi_sess, all_roi_id, confC_mean, labels_all, inputPath. Pooled across sessions in iteration order.
- `master_cosine_fit_phase.m` saves master_cosine_fit_phase_{sig,all}.mat: mu_fit,a_fit,b_fit,nSpk,all_sig_vec,all_exp_vec.

## Cell registration (no existing code)
- NO register/IoU/cross-FOV matching code exists. fov_map.mat is per-session sig-ROI list only. A real cross-session/cross-mouse cell-identity table must be built new. Inputs available: motorPosition [XYZ um] + pixelSize_um (-> absolute centroid coords), maskL (shape/area/centroid), *_AVG_for_CP.tif (morphology), is_sig + th_k (hit + phase for coloring).

## Created Files
- `calcium_spike_gui.m` — calcium spike detection GUI with per-ROI inspection (2026-03-01, updated 2026-03-02)
  - Paired files: `ca_spike_data.mat` (per session), `ca_spike_master.mat` (master folder)
  - Fields: `roi_spikes(n).spike_idx/.spike_t/.spike_amp/.spike_train/.n_spikes`, `spike_params`
  - v2 changes: removed video panel/TIF loading, removed contrast sliders, added Raw F panel (axRawF),
    avg projection loaded from SAM file (meanImg/mean_img/Ims/meanIm/avgImg fields), auto-contrast from ROI mask pixels,
    yellow ROI boundary, 3-row rightGL layout (0.22x/0.42x/0.36x), scanForDFFSessions no longer returns .tifFile
