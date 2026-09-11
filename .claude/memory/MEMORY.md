# Project Memory: RZ_MATLAB Calcium Imaging Pipeline

## Separate Project: CMI-Flu 2026 Challenge (Python, different repo)
- [CMI-Flu 2026 Prediction Challenge](project_cmi_flu_2026_challenge.md) — private repo runqi021/cmi-flu-2026 at D:\RUNQI\cmi-flu-2026; Part I closes 2026-10-06
- [CMI-Flu validation rules](reference_cmi_flu_validation_rules.md) — CMI-PB precedent: baselines beat elaborate models; LOCO + baseline ladder mandatory
- [CMI-Flu strategy + immune models](project_cmi_flu_strategy_and_immune_models.md) — docs/cmi_flu_strategy.md + reproduced models; validation>architecture
- [Kaggle MCP auto-submission](reference_kaggle_mcp_submission.md) — `.mcp.json` loads ONLY from D:\RUNQI\cmi-flu-2026; two-call flow, `sandbox:true` for free format checks
- Prior-attempt understanding lives in the repo at `research/understanding_ledger.md` — keep it updated

## Key Files
- `Breathing_triggered_avg_260226.m` (+ `Batch_` wrapper) — main single-file breathing analysis
- `2p_breathing_coherence/+helper/dFF_RZ.m` — core dF/F computation
- `CLAUDE.md` — full project architecture reference
- Peak-editor convention: `insp_onsets_final` (auto + added − deleted) overrides `insp_onsets` downstream; `Sb` = event train rebuilt from it

## User Preferences
- MATLAB only, no build system, scripts run directly
- Parameters hardcoded at top of scripts as user-editable vars
- Output files saved into `folderPath` alongside input data

## dffQC Pipeline
- [dffQC naming and architecture](project_dffQC_pipeline.md) — single-Z pipeline: rigid + non-rigid MC, cpSAM, dFF with QC plots

## Analysis Architecture
- [Breath-triggered heatmap 260806](project_breath_trig_heatmap_260806.md) — analysis_260806/: single-FOV + pooled-session + cross-genotype; NO sorting, windows in breath cycles, W_hz multitaper, nothing normalised
- [breathing_is_motion module](project_breathing_is_motion.md) — DLC-free breathing from MC |shift|; alignment & Chronux gotchas
- [Phase & time analysis plan](project_phase_time_analysis.md) — per-ROI / per-FOV / per-experiment, phase and time domains
- [Whisker L/R DLC layout](project_whisker_LR_layout.md) — L (Desktop, vL1/vL3) + R (data-folder, vR1/vR2); scorer name Apr8 = project creation, not training
- [Whisker analysis + label-video scripts](project_whisker_analysis_scripts.md) — whisker_L_3w_batch.m + whisker_label_video_{test,batch,select_crop}.py
- [Whisker signal choices](project_whisker_signal_choices.md) — rejected PC1/z-score/rectified; accepted X/Y demeaned + sum-of-PSDs
- [Chronux multitaper PSD conventions](project_chronux_multitaper_conventions.md) — TW=6, tapers=[TW, 2TW-1], fpass=[5 44], lin-lin, per-axis + sum
- [Coherence TW design theory](reference_coherence_TW_design.md) — reference only: 2W<FWHM rule, TW trade-offs at T≈100s
- [Large-scale significance options](reference_large_scale_significance.md) — reference only: Bonferroni/FDR/shuffle/population as N scales 30→20k
- [Rayleigh on linear breath phase](reference_rayleigh_on_linear_phase.md) — occupancy ~7x uneven so correction mandatory; weighted Rbar needs Kish n_eff not raw n; Sert 178°, Vglut2 95°
- [Cosine phase-tuning fit](reference_cosine_phase_tuning.md) — a + b·cos(phase−μ); `master_cosine_fit_phase.m`; 24 bins, nDrop=30
- [repr_export_RZ.m composite video](project_repr_export_video.md) — calcium + breath video + moving trace → MP4; N-ROI generic; ylim keyed by ROI; preview returns before the export section
- [Representative traces figure](project_representative_traces_260527.md) — hand-picked SST/Vgat ROIs, breath+dFF overlay + gamma avg proj; single-MC pairing + front-drop align
- [Rotated-mount sessions](project_rotated_mount_sessions.md) — Sert/0721 + Vglut2/0728 are 90° off; keyed per SESSION; motorPosition is the only authority
- [ChAT ROI quick view](project_chat_roi_quickview.md) — 5-panel breath/Ca view; heatmap sort = post-insp mean dF/F (NOT dt-nearest)
- [**Breath-locked optical artifact**](project_breath_locked_optical_artifact.md) — deep Vgat: top SVD mode tracks chest breathing r=0.96, survives MC, 4x stronger in background than somata → optical, not neural; fakes breath phase-locking
- [Pixel-SVD of the 2P movie](project_pixel_svd_movie_260811.md) — `SVD\Run_SVD_recon_260811.m`; Gram-matrix SVD; top modes are 10 Hz cardiac that MC does NOT remove
- [svd_breath_motion file inventory](project_svd_breath_motion_files.md) — 13 scripts after the 2026-06-06 rename/de-dup; old→new names
- [General breath×Ca coherence pipeline](project_coh_ca_breath_general.md) — coh_ca_breath_260722/ = FLAT single-genotype version; _260606 kept for the nested layout
- [Temporal-trace trajectory (PCA + Isomap)](project_traj_temporal_trace.md) — delay-embed the 1-D trace into an attractor loop; NOT neuron state-space
- [Ventral temporal-phase + batch](project_ventral_temporal_phase.md) — coherence dot must replicate polar stats (TW=4, f_breath clamp, Vglut2/1124 +1-frame fix)
- [Orofacial crop dF/F analysis](project_orofacial_crop_dff_analysis.md) — per-pixel dF/F not motion-energy; 200fps header WRONG, real ~111.3fps
- [Fix AVI timing](project_fix_avi_timing.md) — fix_avi_timing.py: lossless re-time from timestamps.csv → *_fixedfps.mkv
- [Viewable video copies](project_make_viewable_video_copies.md) — FFV1 → H.264 mp4 at TRUE fps; works with no timestamps.csv via _meta.mat scanFrameRate_raw
- [Acquisition-comparable breathing metric](project_breath_fixedmetric.md) — rigid displacement in px replaces per-video PC1; z-score does NOT fix comparability; baseline drifts 0.48x/75min
- [Breath SVD PC1 (DLC-free)](project_breath_svd_pc1.md) — **THE current per-video breath trace**; in `analysis_260727\breath_svd\`; fps from detect_session_fps
- [Breath sign flips but analyses are immune](project_breath_sign_is_analysis_immune.md) — polarity differs 260806 vs 260807; all downstream uses onset EVENTS not the signed trace; SST uses fixedmetric (`src_metric='fb'`) not PC1
- [Thermal .ats breathing pipeline](project_thermal_ats_breathing.md) — FLIR .ats → fnv (env `flir`) → thermal_ats_to_mat.py + thermal_breath_svd.m; Python→.mat→MATLAB split
- [Thermal+Basler sync 260613](project_thermal_basler_sync_260613.md) — 400Hz dual cam; clocks not comparable→align by shared signal; thermal=breathing, Basler=whisker
- [Thermal nostril-ROI breathing (planned)](project_thermal_nostril_roi_breathing.md) — UPCOMING: DLC nostril track → ROI → breathing in °C
- [260615 whisk + breath coordination](project_260615_whisk_breath_coordination.md) — DLC repointed C:\DLC→D:; vL1/vL2=base/tip; atan2d angle; thermal flip in labels
- [IO population sequence + connectivity](project_io_population_seq_connectivity.md) — sequence-first, PCA-after; ca_spike_data is a per-ROI struct array
- [Bilateral whisk sync + breath](project_bilateral_whisk_breath_pipeline.md) — PEAK events (not onset); cameras frame-synced (NO lag shift); MINDIST 0.02 too loose
- [Whisk-breath coordination goal](project_whisk_breath_coordination_goal.md) — L/R synchrony + phase resetting + sniff-vs-basal; intervening-whisk coherence INFEASIBLE
- [Whisk-breath pooled-xcorr method](project_whisk_breath_xcorr_epoch_merge.md) — pooled masked exact-Pearson xcorr + full-trace shuffle null; MERGE_GAP=0
- [Phase resetting method (Moore S4/S5)](reference_phase_resetting_method.md) — φ_reset=2π·t_reset/T0, Δφ=2π(T1−T0)/T1; PRC + period-shift; 30ms drive correction
- [Whisk-breath coherence params](reference_whisk_breath_coherence_params.md) — Moore 1s basal/500ms sniff TW=1; Deschênes 2s TW=2

## Spike Detection
- [**Detection checkpoints — READ FIRST**](project_calcium_detection_checkpoints.md) — what is settled, measured-CLOSED, still open; agent `calcium-dff-detector`
- [CASCADE GCaMP8s auto-detection](project_cascade_gcamp8s_detection.md) — THE label-free detector: CascadeTorch in dlc310; model per cell type; return-order/empty-array traps
- [Calibrated Ca spike detector](project_ca_spike_detector_260722.md) — kernel fit, OASIS vs raw vs matched filter, QC gallery; MC sx/sy do NOT predict spikes → veto must be spatial
- [OASIS on this machine](reference_oasis_on_this_machine.md) — pyenv impossible in R2021b; only `oasisAR1` works, via `ca_oasis_run.py` subprocess

## Optical penetration / falloff  (wildtype vs shiverer)
- [**Optical penetration — READ FIRST**](project_falloff_optical_penetration.md) — `C:\fall-off\<session>`, two acquisition styles, `DEPTH_MODE` rule, 2-channel sessions, ML label direction. **WT ℓ = 24.8 ± 0.4 µm (fit 10–150) vs shiverer 74.0 ± 11.2 µm (10–500).** Never trust motor zs; pia = where (top−bot) peaks
- [**WT vs shiverer vessel figures + boxed falloff**](project_falloff_render_and_box_260901.md) — cube/ray-cast renders, above-pia extension without a 13GB re-read, saturation in MIPs, per-plane norm, real z0 NOT applied, scatter was missing the E-drive session (now n=5)
- [Structural 3-D stack pipeline](project_structural_stack_260910.md) — `structural_stack_260910/`: SI header → group avg per z → 1 µm isotropic → 100 µm-bin MIPs + coronal/sagittal + true 3-D ray-cast, one figure per channel
- [Falloff analysis style](feedback_falloff_analysis_style.md) — no saturation exclusion, no stitching, no unrequested additions, check metadata before a 25 GB re-run
- [Laser power LUT (%→mW)](reference_laser_power_lut.md) — two tables, boundary 2026-07-23; mW ∝ %^1.9 so signal ∝ %^3.8 (using % is the WRONG EXPONENT); pchip
- [Early falloff findings](project_falloff_characterization.md) — superseded: why the mode metric is broken, first 260728 numbers

## Cross-FOV ROI matching (cell registration)
- [**Stage axis convention — MEASURE IT**](reference_stage_axis_convention.md) — image ROW maps to NEGATIVE stage y; wrong sign halved Vgat/0730 matches; `roi_stitch_check_axes_260801.m`
- [Stitch-first ROI grouping](project_roi_stitch_coordinate_first_matching.md) — `roi_stitch_260801/`: position is the metric, appearance only reports; phase-correlation auto-register
- [ROI cross-FOV matching pipeline](project_roi_crossfov_matching.md) — roi_pair_morph_match + roi_match_review_260724; coord=prefilter, morphology=metric; **nccMergeThr=0.40 made one 236 µm false merge — gate on offset-corrected position residual too**
- [**Curation exists ≠ curation applied**](project_cell_link_unapplied_curation.md) — Vglut2/0224 + 1124 were curated but cell_link never ran; `cell_merge_overrides_260814.m`; matcher `fov_name` is fileparts-truncated and collides
- [analysis_260727 breath×Ca pipeline](project_analysis_260727_cell_pooled_coherence.md) — PRIMARY = absolute-time inspiration-triggered PETH + circular-shift null; onset vs peak trigger
- [Ventral_surface_summary archive](project_ventral_summary_archive.md) — Genotype/MMDD/site layout; ACTIVE = >5 events; self-contained since 2026-07-31; 13 blocks / 52 recordings; site guess from filename is routinely wrong, use motorPosition
- [**Session intake drill (stage->pair->link->append->figs)**](project_vagotomized_session_intake_260830.md) — the repeatable per-session chain; Sert/0828 went in UNLINKED; make_cell_link dies on dotted folder names
- [**Selected-cell chain 260816**](project_ventral_selected_cells_260816.md) — archive moved to `D:\Ventral_surface_summary`; per-cell figs -> popsel GUI -> polar_selected/ventral_map (147 cells)
- [**Cell IDs are NOT stable across rebuilds**](feedback_cell_ids_not_stable.md) — positional via dir() order; 256/278 moved on 2026-08-28; translate by cell_key. Vglut2/0824 cell_link row missing from the script
- **Work log: `analysis_260727\PROGRESS_260731.md`** — current state, judgement calls to check, Next list
- [Curation ports across folder moves](feedback_roi_curation_port.md) — (folder name, ROI index)→cell survives any reshuffle; port + verify centroids, never re-curate
- [roi_review4 save bugs](feedback_roi_review4_save_bugs.md) — dropped ungrouped ROIs on save + forgot skips; fixed only in the analysis_260727 copy
- [MATLAB figure-reuse trap](feedback_matlab_offscreen_figure_reuse.md) — clf(hf) does NOT make hf current; batch exports come out blank

## Stitching
- [**AO real-time motion correction**](project_ao_motion_correction_260911.md) — SLAP2-style closed loop using the ALPAO DM, NOT new galvos; SI2018b already ships MotionManager; rig has NO X galvo and NO FastZ; stroke budget + Z2C rows 1/2/3 = tiltY/tipX/defocus; **nothing run on hardware**
- [**Map acquisition pair 260910**](project_auto_acq_260910.md) — auto_acq_260910.m + acq_pause_queue.m; 10s/tile idle removed, live pause/brake/abort; **UNTESTED on the rig as of 2026-09-11**
- [ScanImage rig control + auto-acq](reference_scanimage_rig_control.md) — D:\RZ_ScanImage_script; `scanimage` gives hSI/hSICtl/hAOROIctrl; moveStartRelative is ABSOLUTE; filename x/y is the real motor position; AO = deformable mirror + Shack-Hartmann
- [Auto-stitch: bounded NCC + stage prior](project_autostitch_ncc_260910.md) — why the old phase-corr stitcher aliased; use the VESSEL channel; rebuild the prior at the measured scale
- [Manual stitch GUI](project_manual_stitch_gui.md) — coordinate-based 2D stitcher; anchor-fit pins 2-3 tiles → solves stage→px transform → snaps all
- [Sub-volume oblique render](project_subvolume_render.md) — 3D oblique-MIP of a box cropped from avgz tiles; bothSides +X/-X, medial-left; Frangi rejected
- [Diagonal-wave stack of logZ>=2 cells](project_stack_dff_logz.md) — `analysis_260806\stack_dff_logz_260812.m`; IO is a SITE not a genotype; colours from the polar file; **a MATLAB figure canvas is capped at the screen — a taller axes is silently cut and print() stretches it**
- [stackDFF_chunkALL_260413.m](reference_stackDFF_chunkALL.md) — pools sig ROIs across inputPaths, xcorr-aligns chunks, sheared-parallelogram stacked plot

## Experiments
- [260521_ChAT_soma_G8s DLC](project_260521_chat_soma_dlc.md) — 3-dot breathing project trained 2026-05-26; test RMSE 1.39px; 10 videos @98.6fps
- [ChAT_analysis breath alignment](project_chat_analysis_breath_alignment.md) — breath cam 2P-triggered frame/frame; stored camera_fps=50.7 is bogus; toss-30-then-truncate
- [Breath rate across sessions](project_breathing_freq_across_sessions.md) — per-session 1/ITI from DLC; fps=50.68 from timestamps NOT 90
- [260728 vglut2 + vessel depth series](project_260728_vglut2_vessel_depth.md) — `fileparts` on folder paths silently drops any FOV with a dot; vessel filename z = DEPTH vs SI `zs` = MOTOR

## Pipeline Setup & Data Organization
- [Data folder structure](project_data_folder_structure.md) — calcium (phys/) vs behavior (falling_edge/) two-folder layout
- [Timestamp matching](project_timestamp_matching.md) — AVI filename timestamp + ScanImage epoch; AVI mod times unreliable
- [pair_behavior_to_phys.m](reference_pair_behavior_to_phys.md) — pairs behavior AVIs + DLC CSVs to calcium FOVs, copies, trims
- [Basler dual-camera GUI + backend](project_basler_acquisition_gui.md) — basler_gui.m + basler_dual_acq.py; Line3 falling-edge; **FFV1 lossless** .avi + timestamps.csv per run
- [File reorganization (2026-04)](project_file_reorganization.md) — breathing → Breath_summary_timeNphase_260331/, IO → IO_population_activity/

## Critical Rules
- [TIFF metadata only in raw files](feedback_tiff_metadata.md) — processed TIFFs have NO metadata; use _meta.mat or raw TIFF
- [Always use detect_session_fps](feedback_detect_session_fps.md) — never hardcode/default FPS or zoom
- [Subfolder scripts need repo root on path](feedback_path_setup.md) — addpath(repoRoot) unconditionally; no contains() guard
- [**Root shadow copies beat the path**](feedback_matlab_root_shadow_copies.md) — MATLAB's current folder outranks addpath; a stale root .m silently runs instead of the subfolder one
- [Always use helper. prefix](feedback_helper_package_prefix.md) — helper.dFF_RZ(), never bare dFF_RZ()
- [cpSAM F != raw frame count](feedback_cpSAM_frame_count.md) — TossFrames already removed; use _meta.mat framesPerSlice
- [Keep rerun logic simple](feedback_rerun_logic.md) — rerun = delete + redo, no smart checks
- [DLC train_network: use defaults](feedback_dlc_pytorch_train_defaults.md) — don't pass maxiters/saveiters/displayiters
- [Always save scripts to disk](feedback_script_persistence.md) — never deliver unsaved editor buffers
- [Never modify existing working scripts in-place](feedback_no_inplace_edits.md) — new analyses get a NEW standalone .m file
- [FFV1 codec → Python/OpenCV](feedback_ffv1_video_codec.md) — cv2 not VideoReader; for output use mp4v not MJPG
- [conda ffmpeg is BROKEN](reference_ffmpeg_broken_use_imageio.md) — use imageio_ffmpeg.get_ffmpeg_exe(); cv2 can't write H.264 here
- [FLIR fnv buffer needs copy](reference_fnv_buffer_copy.md) — im.final is a reused VIEW; np.array(copy=True) or all frames alias the last
- [DLC environment + GPU](reference_dlc_environment.md) — dlc310 conda env; RTX A1000 6GB; autocast doesn't help
- [DLC output gotchas](feedback_dlc_outputs_gotchas.md) — clear old model folder before retrain; use `dorsal_R` not `_R`; snapshot suffix marks vintage
- [Don't pester about git tracking](feedback_git_tracking_priority.md) — what matters is the code on disk being good
- [Never gamma exported videos](feedback_never_gamma_video.md) — 2P video frames get LINEAR contrast only; gamma OK on static figures
- [**MATLAB OpenGL export stubs**](reference_matlab_opengl_export_stubs.md) — NVIDIA path fails ~1% of exportgraphics PNGs as 9,223-byte blanks or leaves a STALE file; force `opengl('software')`; verify batches by timestamp AND size, never by count
- [Vectorize over loops; make it fast](feedback_vectorize_over_loops.md) — tic/toc the bottleneck, batch as chunked matrix ops
