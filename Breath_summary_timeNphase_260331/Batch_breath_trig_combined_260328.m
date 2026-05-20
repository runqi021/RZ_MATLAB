%% Batch_breath_trig_combined_260328.m
%  Batch breath-triggered combined analysis.
%
%  Auto-detect mode:
%    - inputPath contains ca_spike_data.mat → sub mode (single folder)
%    - Otherwise → master mode (discover all subfolders with ca_spike_data.mat)
%
%  Per subfolder:
%    Page 0: FOV avg projection with all ROI outlines, 50um scale bar
%    Starter pages: per-ROI avg crop (20%) + breath+dF/F full trace (80%), max 5/page
%
%  Per spiking ROI:  4×2 figure (8 panels)
%    Row 1: breath overlay + avg dF/F in time (col1) | breath overlay + avg dF/F in phase (col2)
%    Row 2: sorted heatmap+Ca raster in time (col1) | sorted heatmap+Ca raster in phase (col2)
%    Row 3: spike count histogram (col1) | spike probability histogram [0,4pi] + cosine fit (col2)
%    Row 4: PSD+coherence inset (col1) | polar (col2)
%
%  Per subfolder summary:  3 rows × 1 col
%    1. Overlaid breathing cycles in phase [-2π, 2π]
%    2. Overlaid cosine fits (all spiking ROIs)
%    3. Overlaid polar plot (all spiking ROIs)
%
%  PDF merging:
%    Per-FOV: p0 overview + starter pages + ROI PDFs + summary → merged PDF
%    Master:  all per-FOV merged PDFs → giant summary
%
%  Dependencies:
%    - helper.dFF_RZ
%    - Chronux (mtspectrumc, coherencyc)
%    - ghostscript on MATLAB path (used by append_pdfs_RZ)
clear all; clc; close all;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(fullfile(repoRoot, 'mod', 'append_pdfs'));
% ======================================================================

%% ============ SHARED PARAMETERS (same in Batch_breath_trig_combined) =========
%D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing
%D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing
%C:\Users\Admin\Desktop\260330_sst_soma_g8s\phys\maybe_breathing
%D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed
%D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing

inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx         = [];             % session indices to skip (master mode)
fps_img         = 30;             % imaging frame rate (Hz)
nDrop           = 30;             % frames to drop from start
fmin            = 0.01;           % Hz (PSD lower bound)
fmax            = 14;             % Hz (PSD upper bound; [] -> fs/2)
TW              = 5;              % multitaper time-bandwidth (K = 2*TW-1 tapers)
f_breath_search = [1.5 2.5];     % Hz (search band for breathing peak; FWHM narrows for coherence)
fwhm_factor     = 1.5;           % coherence band = fwhm_factor x FWHM (centered on peak)
BaselineWinSec  = 20;            % dF/F sliding median baseline window (seconds)
prcLim          = [0.1 99.9];    % percentile clamp for dF/F heatmap
nPhaseBins      = 200;            % bins for phase-normalized breath cycles
coherence_use_spike = true;      % true = coherence on spike train, false = on dF/F waveform
win_sec         = 0.5;            % half-window for peak-trig avg (seconds around breath peak)
alpha_sig       = 0.001;

%% ============ SCRIPT-SPECIFIC PARAMETERS ================================
PixelSizeBase  = 1.7778;        % um/px at 1x zoom (ScanImage default)
crop_um        = 20;            % ROI crop box size (um)
gamma_val      = 0.5;           % gamma correction
clip_lo        = 0.5;           % brightness clip percentile (low)
clip_hi        = 99.5;          % brightness clip percentile (high)
signedSort     = true;          % false = sort by |dt|, true = signed dt
mergePerFOV    = true;          % create per-FOV merged PDF
mergeAllFOV    = true;          % create giant merged PDF across all FOVs
% ====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% 1 ── Auto-detect mode & discover sessions
if exist(fullfile(inputPath, 'ca_spike_data.mat'), 'file') == 2
    % Sub mode: single folder
    fprintf('Sub mode: processing single folder.\n');
    allMat = dir(fullfile(inputPath, 'ca_spike_data.mat'));
    isMasterMode = false;
else
    % Master mode: discover all subfolders
    allMat = dir(fullfile(inputPath, '**', 'ca_spike_data.mat'));
    isMasterMode = true;

    % Apply skipIdx
    if ~isempty(skipIdx)
        skipIdx = skipIdx(:)';
        keep = true(numel(allMat), 1);
        valid_skip = skipIdx(skipIdx >= 1 & skipIdx <= numel(allMat));
        keep(valid_skip) = false;
        allMat = allMat(keep);
    end
end

nSess = numel(allMat);
fprintf('Found %d session(s) with ca_spike_data.mat\n', nSess);
assert(nSess > 0, 'No ca_spike_data.mat found under %s', inputPath);

nProcessed = 0;
nSkipped   = 0;
nFailed    = 0;
nValidFov  = 0;
failedSessions = {};
allFovMergedPdfs = {};

%% 2 ── Main session loop
for kk_sess = 1:nSess
    folderPath = allMat(kk_sess).folder;
    [~, sessName] = fileparts(folderPath);

    fprintf('\n========== [%d/%d] %s ==========\n', kk_sess, nSess, sessName);
    close all;

    sessRoiPdfs   = {};
    sessSummPdf   = '';

    try
        %% ---- Auto-find required files ----
        csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));

        if isempty(csv_hits)
            warning('No DLC CSV found in %s, skipping.', sessName);
            nSkipped = nSkipped + 1; continue;
        end
        if isempty(sam_hits)
            warning('No SAM output found in %s, skipping.', sessName);
            nSkipped = nSkipped + 1; continue;
        end
        if isempty(bp_file)
            warning('No breath_peak_data.mat found in %s, skipping.', sessName);
            nSkipped = nSkipped + 1; continue;
        end

        dlc_csv = fullfile(folderPath, csv_hits(1).name);
        SAM_path = fullfile(folderPath, sam_hits(1).name);

        fprintf('DLC CSV : %s\n', csv_hits(1).name);
        fprintf('SAM     : %s\n', sam_hits(1).name);

        %% ---- Auto-detect fps, load SAM -> dFF ----
        fs = detect_session_fps(folderPath, fps_img);

        SAMload = load(SAM_path);
        F_raw = SAMload.F;
        %F_raw(1:min(nDrop, size(F_raw,1)), :) = [];  % batch QC already saved with tossFrames=0
        dFFout = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF = double(dFFout.dFF);     % [T × N]
        maskL = SAMload.maskL;

        [T, N_roi] = size(dFF);
        if isempty(fmax), fmax_use = fs / 2; else, fmax_use = fmax; end

        fprintf('dFF: [%d x %d] (T x N_roi)\n', T, N_roi);

        %% ---- Read zoom factor from source TIFF → um_per_px ----
        tif_hits = dir(fullfile(folderPath, '*.tif'));
        um_per_px = PixelSizeBase;   % fallback: 1x zoom

        for tt = 1:numel(tif_hits)
            try
                info_t = imfinfo(fullfile(folderPath, tif_hits(tt).name));
                meta_str = "";
                if isfield(info_t(1),'Software') && ~isempty(info_t(1).Software)
                    meta_str = meta_str + string(info_t(1).Software) + newline;
                end
                if isfield(info_t(1),'ImageDescription') && ~isempty(info_t(1).ImageDescription)
                    meta_str = meta_str + string(info_t(1).ImageDescription) + newline;
                end
                tok = regexp(meta_str, ...
                    "SI\.hRoiManager\.scanZoomFactor\s*=\s*([-+]?\d+\.?\d*)", ...
                    "tokens", "once");
                if ~isempty(tok)
                    zf = str2double(tok{1});
                    if isfinite(zf) && zf > 0
                        um_per_px = PixelSizeBase / zf;
                        fprintf('Zoom = %.2fx  ->  um_per_px = %.4f\n', zf, um_per_px);
                        break;
                    end
                end
            catch
            end
        end

        %% ---- Load avg projection ----
        avg_hits = dir(fullfile(folderPath, '*_AVG_for_CP.tif'));
        if ~isempty(avg_hits)
            avgImg = double(imread(fullfile(folderPath, avg_hits(1).name)));
        else
            avgImg = [];
        end

        %% ---- Load calcium spikes ----
        caLoad      = load(fullfile(folderPath, 'ca_spike_data.mat'));
        ifSpk       = caLoad.ifSpike;
        roiSpk_id   = find(ifSpk == 1);
        caData      = caLoad.roi_spikes;
        nSpkROI     = numel(roiSpk_id);

        if nSpkROI == 0
            fprintf('No spiking ROIs, skipping session.\n');
            nSkipped = nSkipped + 1; continue;
        end
        fprintf('Spiking ROIs: %d of %d\n', nSpkROI, N_roi);

        %% ---- Load breath peaks ----
        breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));
        insp_onsets = sort(breath_peak.insp_onset_idx(:));
        fp          = breath_peak.findpeak_params;
        fs_breath = fs;  % breathing camera is frame-triggered 1:1 by calcium

        % Adjust for dropped frames
        insp_onsets(insp_onsets <= nDrop) = [];
        insp_onsets = insp_onsets - nDrop;

        if numel(insp_onsets) < 2
            warning('Fewer than 2 breath peaks after drop, skipping %s.', sessName);
            nSkipped = nSkipped + 1; continue;
        end

        %% ---- Reconstruct breathing signal from DLC ----
        dlc_raw  = readmatrix(dlc_csv, 'NumHeaderLines', 3);
        data_dlc = dlc_raw;
        data_dlc(1:min(nDrop, size(data_dlc,1)), :) = [];

        dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
        dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
        nActiveDots = numel(fp.dot_selection);
        traces   = NaN(size(data_dlc,1), nActiveDots);

        for d = 1:nActiveDots
            di  = dot_idx.(fp.dot_selection{d});
            xc  = data_dlc(:, dot_cols(di,1));
            yc  = data_dlc(:, dot_cols(di,2));
            pc  = data_dlc(:, dot_cols(di,3));
            bad = pc < fp.likelihood_thr;
            switch fp.coord_types{di}
                case 'x',         sig =  xc;
                case 'y',         sig =  yc;
                case '-x',        sig = -xc;
                case '-y',        sig = -yc;
                case 'magnitude', sig = sqrt(xc.^2 + yc.^2);
                otherwise,        sig =  xc;
            end
            sig(bad) = NaN;
            if sum(~isnan(sig)) >= 2
                sig = fillmissing(sig, 'linear', 'EndValues', 'nearest');
            end
            traces(:,d) = sig;
        end

        switch fp.combine_method
            case 'sum',  breath_dlc = sum(traces,  2, 'omitnan');
            otherwise,   breath_dlc = mean(traces, 2, 'omitnan');
        end
        if fp.inverted, breath_dlc = -breath_dlc; end
        breath_dlc = detrend(breath_dlc);
        breath_dlc = (breath_dlc - mean(breath_dlc)) / std(breath_dlc);

        nB = numel(breath_dlc);
        if nB >= T, breath_dlc = breath_dlc(1:T);
        else,       breath_dlc(end+1:T) = 0; end

        %% ---- Apply trim from calcium spike GUI (if set) ----
        if isfield(caLoad, 'trim_info') && isstruct(caLoad.trim_info) ...
                && isfield(caLoad.trim_info, 'trimStart_sec')
            ti = caLoad.trim_info;
            t_full = (0:T-1)' / fs;
            keepMask = true(T, 1);
            if ti.trimStart_sec > 0
                keepMask = keepMask & t_full >= ti.trimStart_sec;
            end
            if isfinite(ti.trimEnd_sec)
                keepMask = keepMask & t_full <= ti.trimEnd_sec;
            end
            trimIdx = find(keepMask);
            if numel(trimIdx) < T
                % Remap breath peak indices
                oldToNew = zeros(T, 1);
                oldToNew(trimIdx) = 1:numel(trimIdx);
                insp_onsets = insp_onsets(insp_onsets >= trimIdx(1) & insp_onsets <= trimIdx(end));
                insp_onsets = oldToNew(insp_onsets);

                dFF        = dFF(trimIdx, :);
                breath_dlc = breath_dlc(trimIdx);
                T          = numel(trimIdx);
                fprintf('  GUI trim applied: %d frames kept (%.1f – %.1f s)\n', ...
                    T, ti.trimStart_sec, ti.trimEnd_sec);
            end
        end

        %% ---- Breathing PSD + peak + FWHM (Chronux) ----
        params_breath.Fs     = fs_breath;
        params_breath.tapers = [TW, 2*TW - 1];
        params_breath.pad    = 0;
        params_breath.fpass  = [fmin, min(fmax_use, fs_breath/2)];
        params_breath.err    = [2, 0.05];
        [Sb, fbC, Serr_b] = mtspectrumc(breath_dlc, params_breath);
        fbC = fbC(:); Sb = Sb(:);

        mask_bp      = fbC >= f_breath_search(1) & fbC <= f_breath_search(2);
        [~, rel_idx] = max(Sb(mask_bp));
        idx_bp       = find(mask_bp, 1) + rel_idx - 1;
        f_breath_peak = fbC(idx_bp);
        half_max      = Sb(idx_bp) / 2;

        idx_lo = idx_bp;
        while idx_lo > 1          && Sb(idx_lo) > half_max, idx_lo = idx_lo - 1; end
        idx_hi = idx_bp;
        while idx_hi < numel(fbC) && Sb(idx_hi) > half_max, idx_hi = idx_hi + 1; end
        f_fwhm = [fbC(idx_lo), fbC(idx_hi)];
        f_fwhm = [max(f_fwhm(1), f_breath_search(1)), min(f_fwhm(2), f_breath_search(2))];
        freq_res = fbC(2) - fbC(1);
        min_bw   = 4 * freq_res;
        if diff(f_fwhm) < min_bw
            f_fwhm = [f_breath_peak - min_bw/2, f_breath_peak + min_bw/2];
            f_fwhm = [max(f_fwhm(1), f_breath_search(1)), min(f_fwhm(2), f_breath_search(2))];
        end
        bw_scaled = diff(f_fwhm) * fwhm_factor;
        f_coh_band = [f_breath_peak - bw_scaled/2, f_breath_peak + bw_scaled/2];
        f_coh_band = [max(f_coh_band(1), f_breath_search(1)), min(f_coh_band(2), f_breath_search(2))];
        fprintf('Breathing peak: %.3f Hz  FWHM: [%.3f, %.3f]  coh band (x%.1f): [%.3f, %.3f] Hz\n', ...
            f_breath_peak, f_fwhm(1), f_fwhm(2), fwhm_factor, f_coh_band(1), f_coh_band(2));

        %% ---- Inspiration timing & windows ----
        b_frames = sort(insp_onsets(:));
        Delta_k  = diff(b_frames) / fs;
        if ~isempty(win_sec)
            win_sec = mean(Delta_k) * 2.5;
        end
        win      = round(win_sec * fs);
        t_w      = (-win:win)' / fs;
        nW       = 2*win + 1;

        edge_margin = max(win, nDrop);
        insp_valid = insp_onsets(insp_onsets > edge_margin & insp_onsets <= T - edge_margin);
        nInsp      = numel(insp_valid);
        fprintf('Breath-trig window: +/-%.3f s  (%d valid onsets)\n', win_sec, nInsp);

        if nInsp < 1
            warning('No valid insp onsets, skipping %s.', sessName);
            nSkipped = nSkipped + 1; continue;
        end

        nValidFov = nValidFov + 1;
        fov_id    = nValidFov;

        %% ---- Pre-compute breath-triggered dFF for all ROIs ----
        dff_trig = cell(N_roi, 1);
        for i = 1:N_roi
            seg = zeros(nInsp, nW);
            for e = 1:nInsp
                c = insp_valid(e);
                seg(e,:) = dFF(c-win : c+win, i);
            end
            dff_trig{i} = seg;
        end

        %% ---- Breath-triggered autocorrelogram ----
        breath_trig_seg = zeros(nInsp, nW);
        for e = 1:nInsp
            c = insp_valid(e);
            breath_trig_seg(e,:) = breath_dlc(c-win : c+win);
        end
        breath_autocorr = mean(breath_trig_seg, 1);
        if max(abs(breath_autocorr)) > 0
            breath_autocorr = breath_autocorr / max(abs(breath_autocorr));
        end

        % xcorr-based autocorrelogram (for row 1 overlay)
        max_lag_frames = round(win_sec * fs);
        [breath_acf, breath_acf_lags] = xcorr(breath_dlc, max_lag_frames, 'coeff');
        breath_acf_t = breath_acf_lags / fs;

        % one-cycle mapping
        one_cycle  = mean(Delta_k);
        half_cycle = one_cycle / 2;
        time_to_phase = @(t) t * 2*pi / one_cycle;
        acf_mask = abs(breath_acf_t) <= half_cycle;
        acf_1cyc = breath_acf(acf_mask);

        % PETH bins
        bin_dt     = 1 / fs;
        peth_edges = -win_sec : bin_dt : win_sec;

        %% ---- Chronux params for per-ROI PSD + coherence ----
        params_mt.Fs     = fs;
        params_mt.tapers = [TW, 2*TW - 1];
        params_mt.pad    = 0;
        params_mt.fpass  = [fmin, fmax_use];
        params_mt.err    = [2, alpha_sig];

        params_coh.Fs     = fs;
        params_coh.tapers = [TW, 2*TW - 1];
        params_coh.pad    = 0;
        params_coh.fpass  = f_coh_band;
        params_coh.err    = [2, alpha_sig];

        % Resample breath to calcium fs if needed
        if abs(fs_breath - fs) > 0.1
            breath_rs = resample(breath_dlc, round(fs), round(fs_breath));
        else
            breath_rs = breath_dlc;
        end
        T_coh     = min(T, numel(breath_rs));
        breath_rs = breath_rs(1:T_coh);

        % confC: compute once
        sig0 = dFF(1:T_coh, roiSpk_id(1)) - mean(dFF(1:T_coh, roiSpk_id(1)));
        [~, ~, ~, ~, ~, ~, ~, confC_val] = coherencyc(breath_rs, sig0, params_coh);

        cmap_roi = lines(nSpkROI);

        %% ---- Figure layout constants (cm) ----
        sq       = 3.2;       % square side
        margin_l = 1.8;
        margin_r = 0.6;
        margin_b = 1.2;
        margin_t = 1.8;
        vgap     = 1.2;
        hgap     = 2.2;

        fig_w = margin_l + 2*sq + hgap + margin_r;
        fig_h = margin_b + 4*sq + 3*vgap + margin_t;

        % y positions (bottom-up)
        row4_y = margin_b;
        row3_y = row4_y + sq + vgap;
        row2_y = row3_y + sq + vgap;
        row1_y = row2_y + sq + vgap;

        % x positions
        col1_x = margin_l;
        col2_x = col1_x + sq + hgap;

        %% ---- Breathing waveform for polar overlay ----
        nPhaseBins_polar = 200;
        phase_1cyc_p = linspace(0, 2*pi, nPhaseBins_polar+1);
        phase_1cyc_p = phase_1cyc_p(1:end-1);
        nCyc_p = numel(b_frames) - 1;
        cycle_mat_p = nan(nCyc_p, nPhaseBins_polar);
        for c = 1:nCyc_p
            i1 = b_frames(c); i2 = b_frames(c+1);
            if i1 >= 1 && i2 <= T
                chunk = breath_dlc(i1:i2);
                cycle_mat_p(c,:) = interp1(linspace(0,2*pi,numel(chunk)), chunk, phase_1cyc_p, 'linear');
            end
        end
        mean_cycle_p = mean(cycle_mat_p, 1, 'omitnan');
        %breath_norm_polar = (mean_cycle_p - min(mean_cycle_p)) / (max(mean_cycle_p) - min(mean_cycle_p) + eps);

        %% ---- Preallocate summary collectors ----
        % For cosine fits
        summ_th_fit  = cell(nSpkROI, 1);
        summ_fit_ext = cell(nSpkROI, 1);
        summ_col     = nan(nSpkROI, 3);
        summ_valid   = false(nSpkROI, 1);

        % For polar
        summ_th_k = nan(nSpkROI, 1);
        summ_r_k  = nan(nSpkROI, 1);
        summ_rlo  = nan(nSpkROI, 1);
        summ_rhi  = nan(nSpkROI, 1);
        summ_dphi = nan(nSpkROI, 1);

        %% ---- Loop over spiking ROIs ----
        for kk_roi = 1:nSpkROI
            roi_id = roiSpk_id(kk_roi);
            col    = cmap_roi(kk_roi,:);
            sig_c  = dFF(:, roi_id) - mean(dFF(:, roi_id));
            ca_ev  = caData(roi_id).spike_idx(:);

            if isempty(ca_ev)
                fprintf('  ROI#%02d : no Ca events, skipping.\n', roi_id);
                continue;
            end

            % --- PSD ---
            [Smt, fmt, Serr_mt] = mtspectrumc(sig_c, params_mt);
            fmt = fmt(:); Smt = Smt(:);

            % --- Waveform coherence (dF/F vs breathing) ---
            [~, C_coh, phi_coh, ~, ~, ~, f_coh, ~, phistd_coh, Cerr_coh] = ...
                coherencyc(breath_rs, sig_c(1:T_coh), params_coh);
            f_coh = f_coh(:); C_coh = C_coh(:); phi_coh = phi_coh(:);

            mask_band = f_coh >= f_coh_band(1) & f_coh <= f_coh_band(2);
            if ~any(mask_band), mask_band = true(size(f_coh)); end
            r_k  = mean(C_coh(mask_band));
            th_k = angle(mean(exp(1i * (-phi_coh(mask_band)))));
            dphi = mean(phistd_coh(mask_band));
            rlo  = max(0,   mean(Cerr_coh(1, mask_band)));
            rhi  = min(1,   mean(Cerr_coh(2, mask_band)));

            % --- Spike-train coherence (binary spike train vs breathing) ---
            spk_train = double(caData(roi_id).spike_train(:));
            spk_train = spk_train(1:min(T_coh, numel(spk_train)));
            if numel(spk_train) < T_coh
                spk_train(end+1:T_coh) = 0;
            end
            spk_train_zm = spk_train - mean(spk_train);

            [~, C_spk, phi_spk, ~, ~, ~, ~, ~, phistd_spk, Cerr_spk] = ...
                coherencyc(breath_rs, spk_train_zm, params_coh);
            C_spk = C_spk(:); phi_spk = phi_spk(:);

            r_k_s  = mean(C_spk(mask_band));
            th_k_s = angle(mean(exp(1i * (-phi_spk(mask_band)))));
            dphi_s = mean(phistd_spk(mask_band));
            rlo_s  = max(0, mean(Cerr_spk(1, mask_band)));
            rhi_s  = min(1, mean(Cerr_spk(2, mask_band)));

            % --- Breath-trig dF/F snippets ---
            seg_dff  = dff_trig{roi_id};
            dff_mean = mean(seg_dff, 1);

            % --- dF/F snippets around each breath onset (for heatmap) ---
            nBreath    = numel(b_frames);
            win_frames = round(win_sec * fs);
            dff_roi    = dFF(:, roi_id);
            snippets   = nan(nBreath, 2*win_frames + 1);
            for b = 1:nBreath
                bp = b_frames(b);
                i1 = bp - win_frames;
                i2 = bp + win_frames;
                if i1 >= 1 && i2 <= T
                    snippets(b, :) = dff_roi(i1:i2);
                end
            end
            valid_snip = ~all(isnan(snippets), 2);
            nValid     = sum(valid_snip);
            t_snip     = (-win_frames:win_frames) / fs;
            mean_dff_snip = mean(snippets(valid_snip, :), 1, 'omitnan');

            valid_idx      = find(valid_snip);
            b_frames_valid = b_frames(valid_idx);
            win5           = round(win_sec * fs);
            ca_rel_valid   = cell(nValid, 1);
            for b = 1:nValid
                bp = b_frames_valid(b);
                nearby_ca = ca_ev(ca_ev >= (bp - win5) & ca_ev <= (bp + win5));
                ca_rel_valid{b} = (nearby_ca - bp) / fs;
            end

            % --- Ca events in breath intervals (for dt / phase) ---
            [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
            keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
            ca_in   = ca_ev(keep_ca);
            bin_in  = bin_idx(keep_ca);
            nEv_in  = numel(ca_in);

            % dt_last / dt_next per spike (for phase)
            dt_last_spk = nan(numel(ca_ev), 1);
            dt_next_spk = nan(numel(ca_ev), 1);
            for e = 1:numel(ca_ev)
                c   = ca_ev(e);
                bef = insp_onsets(insp_onsets <= c);
                aft = insp_onsets(insp_onsets >  c);
                if ~isempty(bef), dt_last_spk(e) = (c - bef(end)) / fs; end
                if ~isempty(aft), dt_next_spk(e) = (aft(1) - c)   / fs; end
            end
            keep3 = ~isnan(dt_last_spk) & ~isnan(dt_next_spk);
            xl3   = dt_last_spk(keep3);
            xn3   = dt_next_spk(keep3);

            % --- PETH ---
            all_ca_rel  = cell2mat(ca_rel_valid);
            peth_counts = histcounts(all_ca_rel, peth_edges);
            peth_rate   = peth_counts / (nValid * bin_dt);
            peth_centers = peth_edges(1:end-1) + bin_dt/2;

            % --- Sort by dt_nearest_ca ---
            dt_nearest_ca = nan(nValid, 1);
            for b = 1:nValid
                bp = b_frames_valid(b);
                before = ca_ev(ca_ev <  bp);
                after  = ca_ev(ca_ev >= bp);
                dt_bef = Inf;  dt_aft = Inf;
                if ~isempty(before), dt_bef = (bp - before(end)) / fs; end
                if ~isempty(after),  dt_aft = (after(1) - bp)    / fs; end
                if signedSort
                    if dt_bef <= dt_aft
                        dt_nearest_ca(b) = -dt_bef;
                    else
                        dt_nearest_ca(b) =  dt_aft;
                    end
                else
                    dt_nearest_ca(b) = min(dt_bef, dt_aft);
                end
            end
            dt_near_sort = dt_nearest_ca;
            dt_near_sort(isnan(dt_near_sort)) = Inf;
            [~, sort_near_b] = sort(dt_near_sort, 'ascend');

            % sorted snippet matrix
            snip_valid  = snippets(valid_idx, :);
            snip_sorted = snip_valid(sort_near_b, :);
            cHigh = prctile(snip_sorted(:), prcLim(2));

            if nValid < 1 || isempty(mean_dff_snip) || isempty(all_ca_rel) || isempty(xl3)
                fprintf('  ROI#%02d : insufficient data for figure, skipping.\n', roi_id);
                continue;
            end

            % --- Phase histogram + cosine fit data ---
            phi_spikes = 2*pi * xl3 ./ (xl3 + xn3);
            nBins    = 20;
            edges_ph = linspace(0, 2*pi, nBins+1);
            ctrs     = (edges_ph(1:end-1) + edges_ph(2:end)) / 2;
            counts   = histcounts(phi_spikes, edges_ph);
            nCycles  = numel(b_frames) - 1;
            counts_pc = counts / max(nCycles, 1);   % spikes per cycle per bin
            ctrs_ext   = [ctrs, ctrs + 2*pi];
            cnts_ext   = [counts_pc, counts_pc];
            r_vec = mean(exp(1i * phi_spikes));
            phi0  = angle(r_vec);
            Rlen  = abs(r_vec);
            B_cos = numel(phi_spikes) / (nBins * max(nCycles, 1));
            A_cos = 2 * B_cos * Rlen;
            th_fit  = linspace(0, 4*pi, 300);
            fit_ext = B_cos + A_cos * cos(th_fit - phi0);

            % Store for summary
            summ_th_fit{kk_roi}  = th_fit;
            summ_fit_ext{kk_roi} = fit_ext;
            summ_col(kk_roi,:)   = col;
            summ_valid(kk_roi)   = true;
            if coherence_use_spike
                summ_th_k(kk_roi)  = th_k_s;
                summ_r_k(kk_roi)   = r_k_s;
                summ_rlo(kk_roi)   = rlo_s;
                summ_rhi(kk_roi)   = rhi_s;
                summ_dphi(kk_roi)  = dphi_s;
            else
                summ_th_k(kk_roi)  = th_k;
                summ_r_k(kk_roi)   = r_k;
                summ_rlo(kk_roi)   = rlo;
                summ_rhi(kk_roi)   = rhi;
                summ_dphi(kk_roi)  = dphi;
            end

            % Per-ROI polar overlay uses the same source as summary
            if coherence_use_spike
                th_k_plot = th_k_s;  r_k_plot = r_k_s;
                rlo_plot  = rlo_s;   rhi_plot = rhi_s;   dphi_plot = dphi_s;
            else
                th_k_plot = th_k;    r_k_plot = r_k;
                rlo_plot  = rlo;     rhi_plot = rhi;     dphi_plot = dphi;
            end

            fprintf('  ROI#%02d : %d Ca events, %d valid breaths\n', ...
                roi_id, nEv_in, nValid);

            % ================================================================
            % FIGURE  4 rows x 2 columns
            % ================================================================
            fig = figure('Color', 'w', 'Visible', 'off', ...
                'Name', sprintf('FOV#%02d | %s | ROI#%02d', fov_id, sessName, roi_id), ...
                'Units', 'centimeters', 'Position', [1 0.5 fig_w fig_h], ...
                'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
                'PaperPosition', [0 0 fig_w fig_h]);
            set(fig, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 8);

            pageLabel = sprintf('FOV#%02d | %s | ROI#%02d', fov_id, sessName, roi_id);
            annotation(fig, 'textbox', [0.02 0.975 0.96 0.02], ...
                'String', pageLabel, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'Interpreter', 'none', ...
                'FontName', 'Arial', ...
                'FontSize', 8, ...
                'Color', [0 0 0]);

            % ---- (4,1) PSD + coherence inset ----
            ax_22 = axes(fig, 'Units','centimeters', ...
                'Position', [col1_x row4_y sq sq]);
            hold(ax_22, 'on');
            yLo = Serr_mt(1,:)'; yHi = Serr_mt(2,:)';
            patch(ax_22, [fmt; flipud(fmt)], [yLo; flipud(yHi)], col, ...
                'EdgeColor','none', 'FaceAlpha',0.2);
            plot(ax_22, fmt, Smt, 'Color', col, 'LineWidth', 1.5);
            yLo_b = Serr_b(1,:)'; yHi_b = Serr_b(2,:)';
            patch(ax_22, [fbC; flipud(fbC)], [yLo_b; flipud(yHi_b)], [0.3 0.3 0.3], ...
                'EdgeColor','none', 'FaceAlpha',0.15);
            plot(ax_22, fbC, Sb, 'k', 'LineWidth', 2);
            ax_22.XScale = 'log'; ax_22.YScale = 'log';
            psd_xmin = min(fmt(1), fbC(1));
            ax_22.XLim = [0.01, 14];
            all_ticks  = [0.01 0.1 1 10];
            all_labels = {'10^{-2}','10^{-1}','10^0','10^1'};
            keep_t     = all_ticks >= psd_xmin & all_ticks <= fmax_use;
            ax_22.XTick      = all_ticks(keep_t);
            ax_22.XTickLabel = all_labels(keep_t);
            grid(ax_22, 'on');
            xlabel(ax_22, 'Freq (Hz)'); ylabel(ax_22, 'PSD');
            hold(ax_22, 'off'); box(ax_22, 'on');
            pbaspect(ax_22, [1 1 1]);

            % ---- (1,1) Breath overlay + avg dF/F (time domain) ----
            % Normalize each breath snippet to [0,1]
            breath_snip_norm = breath_trig_seg;
            for e_ = 1:size(breath_snip_norm, 1)
                mn_ = min(breath_snip_norm(e_,:));
                mx_ = max(breath_snip_norm(e_,:));
                if mx_ > mn_
                    breath_snip_norm(e_,:) = (breath_snip_norm(e_,:) - mn_) / (mx_ - mn_);
                end
            end
            dff_roi_snips = dff_trig{roi_id};
            mean_dff_t = mean(dff_roi_snips, 1, 'omitnan');

            ax_21 = axes(fig, 'Units','centimeters', ...
                'Position', [col1_x row1_y sq sq]);
            hold(ax_21, 'on');
            plot(ax_21, t_w, breath_snip_norm', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
            xline(ax_21, 0, 'k--', 'LineWidth', 1);
            yyaxis(ax_21, 'right');
            plot(ax_21, t_w, mean_dff_t, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
            ylabel(ax_21, 'dF/F');
            ax_21.YAxis(2).Color = [0.2 0.7 0.2];
            yyaxis(ax_21, 'left');
            ylabel(ax_21, 'Breath (norm)');
            ax_21.YAxis(1).Color = [0 0 0];
            hold(ax_21, 'off');
            xlim(ax_21, [-win_sec, win_sec]);
            title(ax_21, 'Time domain');
            box(ax_21, 'on'); pbaspect(ax_21, [1 1 1]);

            % ---- (1,2) Breath overlay + avg dF/F (phase domain) ----
            % Phase-normalize dF/F for this ROI
            dff_roi_vec = dFF(:, roi_id);
            nCyc_r = numel(insp_valid) - 1;
            dff_phase_roi = nan(nCyc_r, nPhaseBins_polar);
            for c_ = 1:nCyc_r
                i1_ = insp_valid(c_); i2_ = insp_valid(c_+1);
                if i1_ >= 1 && i2_ <= T
                    chunk_ = dff_roi_vec(i1_:i2_);
                    dff_phase_roi(c_,:) = interp1(linspace(0,2*pi,numel(chunk_)), chunk_, phase_1cyc_p, 'linear');
                end
            end
            mean_dff_ph = mean(dff_phase_roi, 1, 'omitnan');

            % Normalize breath cycles to [0,1] for display
            cycle_disp = cycle_mat_p;
            for c_ = 1:size(cycle_disp, 1)
                mn_ = min(cycle_disp(c_,:)); mx_ = max(cycle_disp(c_,:));
                if mx_ > mn_, cycle_disp(c_,:) = (cycle_disp(c_,:) - mn_) / (mx_ - mn_); end
            end
            phase_2cyc = [phase_1cyc_p, phase_1cyc_p + 2*pi];
            cycle_disp_2x = [cycle_disp, cycle_disp];
            mean_dff_ph_2x = [mean_dff_ph, mean_dff_ph];

            ax_22b = axes(fig, 'Units','centimeters', ...
                'Position', [col2_x row1_y sq sq]);
            hold(ax_22b, 'on');
            plot(ax_22b, phase_2cyc, cycle_disp_2x', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
            xline(ax_22b, 2*pi, 'k--', 'LineWidth', 1);
            yyaxis(ax_22b, 'right');
            plot(ax_22b, phase_2cyc, mean_dff_ph_2x, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
            ylabel(ax_22b, 'dF/F');
            ax_22b.YAxis(2).Color = [0.2 0.7 0.2];
            yyaxis(ax_22b, 'left');
            ylabel(ax_22b, 'Breath (norm)');
            ax_22b.YAxis(1).Color = [0 0 0];
            hold(ax_22b, 'off');
            xlim(ax_22b, [0, 4*pi]);
            set(ax_22b, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
                'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
            title(ax_22b, 'Phase domain');
            box(ax_22b, 'on'); pbaspect(ax_22b, [1 1 1]);

            % ---- (2,1) Sorted dF/F heatmap + Ca raster overlay (time domain) ----
            ax_31 = axes(fig, 'Units','centimeters', ...
                'Position', [col1_x row2_y sq sq]);
            imagesc(ax_31, t_snip, 1:nValid, snip_sorted);
            set(ax_31, 'YDir', 'reverse');
            colormap(ax_31, flipud(gray(256)));
            caxis(ax_31, [0 cHigh]);
            hold(ax_31, 'on');
            xline(ax_31, 0, 'k--', 'LineWidth', 1);
            % overlay Ca spike raster
            for row = 1:nValid
                b_s     = sort_near_b(row);
                t_ticks = ca_rel_valid{b_s};
                for t_t = 1:numel(t_ticks)
                    plot(ax_31, [t_ticks(t_t) t_ticks(t_t)], [row-0.4 row+0.4], ...
                        'Color', col, 'LineWidth', 1.5);
                end
            end
            hold(ax_31, 'off');
            xlim(ax_31, [-win_sec, win_sec]);
            ylim(ax_31, [0.5, nValid + 0.5]);
            set(ax_31, 'XTickLabel', []);
            ylabel(ax_31, 'Breath # (sorted)');
            % thin colorbar manually positioned
            cb = colorbar(ax_31);
            ylabel(cb, 'dF/F');
            set(ax_31, 'Units','centimeters', 'Position', [col1_x row2_y sq sq]);
            cb.Units = 'centimeters';
            cb.Position = [col1_x + sq + 0.1, row2_y, 0.15, sq];

            % ---- (2,2) Sorted dF/F heatmap + Ca raster (phase domain) ----
            % Sort cycles by nearest spike phase
            nCyc_r_ = size(dff_phase_roi, 1);
            valid_ph_rows = ~all(isnan(dff_phase_roi), 2);
            valid_ph_idx = find(valid_ph_rows);

            % For each valid cycle, find spikes and their phases
            spk_ph_per_cyc = cell(numel(valid_ph_idx), 1);
            min_spk_ph = nan(numel(valid_ph_idx), 1);
            for cc_ = 1:numel(valid_ph_idx)
                ci_ = valid_ph_idx(cc_);
                if ci_ <= numel(insp_valid)-1
                    bp1_ = insp_valid(ci_); bp2_ = insp_valid(ci_+1);
                    spk_ = ca_ev(ca_ev >= bp1_ & ca_ev < bp2_);
                    if ~isempty(spk_)
                        spk_ph_per_cyc{cc_} = 2*pi * (spk_ - bp1_ + 0.5) / (bp2_ - bp1_);
                        min_spk_ph(cc_) = min(spk_ph_per_cyc{cc_});
                    end
                end
            end

            % Sort by min spike phase (cycles without spikes go to end)
            min_spk_ph(isnan(min_spk_ph)) = Inf;
            [~, sort_ph_ord] = sort(min_spk_ph, 'ascend');

            % Show ALL valid cycles (spiking sorted to top, non-spiking at bottom)
            dff_ph_sorted = dff_phase_roi(valid_ph_idx(sort_ph_ord), :);
            dff_ph_sorted_2x = [dff_ph_sorted, dff_ph_sorted];

            ax_32b = axes(fig, 'Units','centimeters', ...
                'Position', [col2_x row2_y sq sq]);
            if ~isempty(dff_ph_sorted_2x)
                imagesc(ax_32b, phase_2cyc, 1:size(dff_ph_sorted_2x,1), dff_ph_sorted_2x);
                set(ax_32b, 'YDir','reverse');
                colormap(ax_32b, flipud(gray(256)));
                cH_ph = prctile(dff_ph_sorted_2x(:), 99.9);
                if cH_ph > 0, caxis(ax_32b, [0 cH_ph]); end
                hold(ax_32b, 'on');
                % Overlay spike phase ticks
                for rr_ = 1:numel(sort_ph_ord)
                    spk_p = spk_ph_per_cyc{sort_ph_ord(rr_)};
                    if ~isempty(spk_p)
                        for sp_ = 1:numel(spk_p)
                            plot(ax_32b, [spk_p(sp_) spk_p(sp_)], [rr_-0.4 rr_+0.4], ...
                                'Color', col, 'LineWidth', 2);
                            plot(ax_32b, [spk_p(sp_)+2*pi spk_p(sp_)+2*pi], [rr_-0.4 rr_+0.4], ...
                                'Color', col, 'LineWidth', 2);
                        end
                    end
                end
                xline(ax_32b, 2*pi, 'k--', 'LineWidth', 0.8);
                hold(ax_32b, 'off');
            end
            xlim(ax_32b, [0, 4*pi]);
            if ~isempty(dff_ph_sorted_2x)
                ylim(ax_32b, [0.5, size(dff_ph_sorted_2x,1)+0.5]);
            end
            set(ax_32b, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
                'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
            ylabel(ax_32b, 'Cycle (sorted)');
            set(ax_32b, 'XTickLabel', []);
            cb_ph = colorbar(ax_32b);
            ylabel(cb_ph, 'dF/F');
            set(ax_32b, 'Units','centimeters', 'Position', [col2_x row2_y sq sq]);
            cb_ph.Units = 'centimeters';
            cb_ph.Position = [col2_x + sq + 0.1, row2_y, 0.15, sq];

            % ---- (3,2) Spike probability histogram [0, 4pi] ----
            ax_32 = axes(fig, 'Units','centimeters', ...
                'Position', [col2_x row3_y sq sq]);
            hold(ax_32, 'on');
            nBinsProb = 12;
            edges_prob = linspace(0, 2*pi, nBinsProb + 1);
            ctrs_prob  = (edges_prob(1:end-1) + edges_prob(2:end)) / 2;
            cnt_prob   = histcounts(phi_spikes, edges_prob);
            prob_prob  = cnt_prob / max(sum(cnt_prob), 1);

            % Cosine fit
            mu_prob = angle(mean(exp(1i * phi_spikes)));
            a_prob  = mean(prob_prob);
            b_prob  = 2 * mean(prob_prob .* cos(ctrs_prob - mu_prob));

            % Duplicate to [0, 4pi]
            ctrs_prob_dup = [ctrs_prob, ctrs_prob + 2*pi];
            prob_prob_dup = [prob_prob, prob_prob];

            bar(ax_32, ctrs_prob_dup, prob_prob_dup, 1, ...
                'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            x_fit_prob = linspace(0, 4*pi, 500);
            y_fit_prob = a_prob + b_prob * cos(x_fit_prob - mu_prob);
            plot(ax_32, x_fit_prob, y_fit_prob, 'r-', 'LineWidth', 2);
            xline(ax_32, 2*pi, 'k--', 'LineWidth', 1);
            xlim(ax_32, [0, 4*pi]);
            set(ax_32, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
                'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
            xlabel(ax_32, 'Phase (rad)');
            ylabel(ax_32, 'P(spike)');
            hold(ax_32, 'off'); box(ax_32, 'on');
            pbaspect(ax_32, [1 1 1]);

            % ---- (3,1) Spike count histogram ----
            ax_41 = axes(fig, 'Units','centimeters', ...
                'Position', [col1_x row3_y sq sq]);
            hold(ax_41, 'on');
            dt_bin_sc = t_snip(2) - t_snip(1);
            sc_edges = [t_snip(:) - dt_bin_sc/2; t_snip(end) + dt_bin_sc/2];
            sc_counts = histcounts(all_ca_rel, sc_edges);
            bar(ax_41, t_snip, sc_counts, 1, ...
                'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            xline(ax_41, 0, 'k--', 'LineWidth', 1);
            xlim(ax_41, [-win_sec, win_sec]);
            xlabel(ax_41, 'Time from breath peak (s)');
            ylabel(ax_41, 'Spike count');
            hold(ax_41, 'off'); box(ax_41, 'on');
            pbaspect(ax_41, [1 1 1]);

            % ---- (4,2) Polar plot ----
            fig_pos = fig.Position;
            polar_pos_norm = [col2_x/fig_pos(3), row4_y/fig_pos(4), ...
                              sq/fig_pos(3),     sq/fig_pos(4)];
            ax_42 = polaraxes(fig, 'Position', polar_pos_norm);
            hold(ax_42, 'on');
            th_circ = linspace(0, 2*pi, 360);
            polarplot(ax_42, th_circ, repmat(confC_val, 1, 360), 'k--', 'LineWidth', 1.2);
            %polarplot(ax_42, [phase_1cyc_p, 2*pi], [breath_norm_polar, breath_norm_polar(1)], 'k-', 'LineWidth', 0.8);
            polarplot(ax_42, [th_k_plot th_k_plot], [rlo_plot rhi_plot], '-', 'Color', col, 'LineWidth', 2);
            phi_arc = linspace(th_k_plot - dphi_plot, th_k_plot + dphi_plot, 40);
            polarplot(ax_42, phi_arc, r_k_plot * ones(1,40), '-', 'Color', col, 'LineWidth', 1.5);
            ph = polarplot(ax_42, th_k_plot, r_k_plot, 'o', ...
                'MarkerFaceColor', col, 'MarkerEdgeColor','k', 'MarkerSize', 3);
            ph.UserData = struct('roi', roi_id, 'phistd', dphi_plot, 'rlo', rlo_plot, 'rhi', rhi_plot);
            ax_42.RLim              = [0, 1];
            ax_42.ThetaZeroLocation = 'top';
            ax_42.ThetaDir          = 'clockwise';
            ax_42.RAxisLocation     = 180;
            ax_42.FontSize          = 7;
            hold(ax_42, 'off');

            dcm = datacursormode(fig);
            dcm.UpdateFcn = @(~, info) polarTip(info);

            drawnow;

            % ---- Coherence inset (bottom-left of PSD panel) ----
            rp   = ax_22.Position;
            inW  = 0.38 * rp(3);
            inH  = 0.38 * rp(4);
            inX  = rp(1) + 0.12 * rp(3);
            inY  = rp(2) + 0.12 * rp(4);
            ax_ins = axes(fig, 'Units','centimeters', 'Position',[inX inY inW inH]);
            hold(ax_ins, 'on');
            plot(ax_ins, f_coh, C_coh, '--', 'Color', col, 'LineWidth', 1);
            plot(ax_ins, f_coh, C_spk, '-', 'Color', col, 'LineWidth', 1.2);
            yline(ax_ins, confC_val, 'k', 'LineWidth', 0.8);
            ax_ins.XLim  = [f_coh_band(1), f_coh_band(2)];
            ax_ins.YLim  = [0, 1];
            ax_ins.XTick = [f_coh_band(1), f_coh_band(2)];
            ax_ins.XTickLabel = {sprintf('%.2f',f_coh_band(1)), sprintf('%.2f',f_coh_band(2))};
            ax_ins.YTick = [0, 1];
            ylabel(ax_ins, '|C|');
            ax_ins.XScale  = 'linear';
            ax_ins.FontSize = 6;
            ax_ins.Color = [1 1 1 0.85];
            grid(ax_ins, 'on');
            hold(ax_ins, 'off'); box(ax_ins, 'on');
            axis(ax_ins, 'square');

            % ---- Save per-ROI figure ----
            outDir = fullfile(folderPath, 'breath_trig_combined');
            if ~exist(outDir, 'dir'), mkdir(outDir); end

            outFig = fullfile(outDir, sprintf('breath_trig_combined_ROI%02d.fig', roi_id));
            savefig(fig, outFig);

            outPdf = fullfile(outDir, sprintf('breath_trig_combined_ROI%02d.pdf', roi_id));
            exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','none');

            fprintf('    Saved: %s\n', outPdf);
            sessRoiPdfs{end+1} = outPdf; %#ok<SAGROW>

            close(fig);
        end  % end ROI loop

        %% ---- Per-subfolder summary (3 rows × 1 col, squared panels) ----
        summ_use = find(summ_valid);
        if isempty(summ_use)
            fprintf('  No valid spiking ROIs for summary, skipping.\n');
            nSkipped = nSkipped + 1;
            continue;
        end

        % --- Phase-normalize each breathing cycle ---
        phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
        phase_1cyc = phase_1cyc(1:end-1);   % [0, 2pi)
        nPeaks_  = numel(b_frames);
        nCyc_    = nPeaks_ - 1;
        cycle_mat = nan(nCyc_, nPhaseBins);
        for c = 1:nCyc_
            i1 = b_frames(c);
            i2 = b_frames(c + 1);
            if i1 >= 1 && i2 <= T
                chunk = breath_dlc(i1:i2);
                orig_phase = linspace(0, 2*pi, numel(chunk));
                cycle_mat(c, :) = interp1(orig_phase, chunk, phase_1cyc, 'linear');
            end
        end
        valid_cyc = ~all(isnan(cycle_mat), 2);
        cycle_mat = cycle_mat(valid_cyc, :);
        nValidCyc = size(cycle_mat, 1);
        mean_cycle = mean(cycle_mat, 1, 'omitnan');

        % duplicate to [-2pi, 2pi]
        phase_full = [phase_1cyc - 2*pi, phase_1cyc];
        cycle_dup  = [cycle_mat, cycle_mat];
        mean_dup   = [mean_cycle, mean_cycle];

        % --- Page 0: FOV overview — avg projection with all ROI outlines ---
        outDir_st = fullfile(folderPath, 'breath_trig_combined');
        if ~exist(outDir_st, 'dir'), mkdir(outDir_st); end

        if ~isempty(avgImg)
            p0_w = 18; p0_h = 18;
            fig_p0 = figure('Color','w', 'Visible','off', ...
                'Name', sprintf('FOV#%02d | %s | overview', fov_id, sessName), ...
                'Units','centimeters', 'Position', [1 1 p0_w p0_h], ...
                'PaperUnits','centimeters', 'PaperSize', [p0_w p0_h], ...
                'PaperPosition', [0 0 p0_w p0_h]);

            ax_p0 = axes(fig_p0, 'Units','normalized', 'Position', [0.05 0.05 0.9 0.88]);

            % Gamma-corrected, clipped avg projection
            lo_p0 = prctile(avgImg(:), 0.5);
            hi_p0 = prctile(avgImg(:), 99.5);
            if hi_p0 <= lo_p0, hi_p0 = lo_p0 + eps; end
            fovImg = (avgImg - lo_p0) / (hi_p0 - lo_p0);
            fovImg = max(0, min(1, fovImg));
            fovImg = fovImg .^ 0.6;

            % RGB with ROI outlines (yellow)
            fovRGB = repmat(fovImg, [1 1 3]);
            for ri = 1:N_roi
                perim_ri = bwperim(maskL == ri);
                for ch = 1:3
                    plane = fovRGB(:,:,ch);
                    plane(perim_ri) = [1 1 0] * (ch == [1;2;3]);
                    fovRGB(:,:,ch) = plane;
                end
            end

            imshow(fovRGB, 'Parent', ax_p0);
            hold(ax_p0, 'on');

            % ROI labels
            for ri = 1:N_roi
                [ry_, rx_] = find(maskL == ri);
                if ~isempty(ry_)
                    text(ax_p0, mean(rx_), mean(ry_), sprintf('%d', ri), ...
                        'Color', [1 1 0], 'FontSize', 7, ...
                        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            end

            % 50 um scale bar
            sb_px_50 = round(50 / um_per_px);
            [imgH_p0, imgW_p0, ~] = size(fovRGB);
            sb_y0 = imgH_p0 - 10;
            plot(ax_p0, [10 10+sb_px_50], [sb_y0 sb_y0], 'w-', 'LineWidth', 3);
            text(ax_p0, 10+sb_px_50/2, sb_y0-12, '50 \mum', 'Color', 'w', ...
                'FontSize', 9, 'HorizontalAlignment', 'center');
            hold(ax_p0, 'off');

            title(ax_p0, sprintf('FOV#%02d  %s  (%d ROIs)', fov_id, strrep(sessName,'_',' '), N_roi), ...
                'Interpreter', 'none', 'FontSize', 10);

            p0Pdf = fullfile(outDir_st, sprintf('FOV%02d_%s_p0_overview.pdf', fov_id, sessName));
            exportgraphics(fig_p0, p0Pdf, 'ContentType','vector', 'BackgroundColor','none');
            fprintf('  FOV overview page saved: %s\n', p0Pdf);
            close(fig_p0);
        end

        % --- Starter pages: crop | breath+dFF trace (70%) | spike/non-spike (square) ---
        sessStarterPdfs = {};
        t_full_st = (0:T-1)' / fs;
        maxPerPage = 5;
        nPages = ceil(nSpkROI / maxPerPage);

        for iPage = 1:nPages
            roi_start = (iPage-1)*maxPerPage + 1;
            roi_end   = min(iPage*maxPerPage, nSpkROI);
            nThisPage = roi_end - roi_start + 1;

            st_w = 28; st_h = 4 * nThisPage + 2;
            st_ml = 0.8; st_mr = 0.5; st_mb = 1.2; st_mt = 1.5; st_vgap = 0.3;
            st_hgap = 0.3;
            usable_w = st_w - st_ml - st_mr;
            row_h = (st_h - st_mb - st_mt - (nThisPage-1)*st_vgap) / nThisPage;

            trace_w_cm = usable_w * 0.70;
            spk_w_cm   = row_h;   % square
            crop_w_cm  = usable_w - trace_w_cm - spk_w_cm - 2*st_hgap;
            crop_x     = st_ml;
            trace_x    = crop_x + crop_w_cm + st_hgap;
            spk_x      = trace_x + trace_w_cm + st_hgap;

            fig_st = figure('Color','w', 'Visible','off', ...
                'Name', sprintf('FOV#%02d | %s | starter p%d', fov_id, sessName, iPage), ...
                'Units','centimeters', 'Position', [1 1 st_w st_h], ...
                'PaperUnits','centimeters', 'PaperSize', [st_w st_h], ...
                'PaperPosition', [0 0 st_w st_h]);

            if nPages > 1
                annotation(fig_st, 'textbox', [0.02 0.97 0.96 0.03], ...
                    'String', sprintf('%s  breath + dF/F  (page %d/%d)', strrep(sessName,'_',' '), iPage, nPages), ...
                    'EdgeColor','none', 'FontSize',10, 'Interpreter','none');
            else
                annotation(fig_st, 'textbox', [0.02 0.97 0.96 0.03], ...
                    'String', sprintf('%s  breath + dF/F', strrep(sessName,'_',' ')), ...
                    'EdgeColor','none', 'FontSize',10, 'Interpreter','none');
            end

            for jj_st = roi_start:roi_end
                rid = roiSpk_id(jj_st);
                row_idx = jj_st - roi_start;  % 0-based
                row_y = st_h - st_mt - (row_idx+1)*row_h - row_idx*st_vgap;

                % Col 1: avg projection crop for this ROI
                ax_crop = axes(fig_st, 'Units','centimeters', ...
                    'Position', [crop_x row_y crop_w_cm row_h]);
                if ~isempty(avgImg)
                    [ry_, rx_] = find(maskL == rid);
                    if ~isempty(ry_)
                        cy_ = round(mean(ry_)); cx_ = round(mean(rx_));
                        crop_px_ = round(crop_um / um_per_px / 2);
                        [imgH_, imgW_] = size(avgImg);
                        r1_ = max(1, cy_-crop_px_); r2_ = min(imgH_, cy_+crop_px_);
                        c1_ = max(1, cx_-crop_px_); c2_ = min(imgW_, cx_+crop_px_);
                        rp_ = avgImg(r1_:r2_, c1_:c2_);
                        lo_ = prctile(rp_(:), clip_lo); hi_ = prctile(rp_(:), clip_hi);
                        if hi_ <= lo_, hi_ = lo_ + eps; end
                        rp_ = max(0, min(1, (rp_ - lo_) / (hi_ - lo_)));
                        rp_ = rp_ .^ gamma_val;
                        mc_ = maskL(r1_:r2_, c1_:c2_);
                        perim_ = bwperim(mc_ == rid);
                        rpRGB_ = repmat(rp_, [1 1 3]);
                        for ch = 1:3
                            pl_ = rpRGB_(:,:,ch); pl_(perim_) = [1 1 0]*[ch==1; ch==2; ch==3];
                            rpRGB_(:,:,ch) = pl_;
                        end
                        imshow(rpRGB_, 'Parent', ax_crop);
                        hold(ax_crop, 'on');
                        sb5 = round(5 / um_per_px);
                        [crH_, crW_, ~] = size(rpRGB_);
                        sb_y5 = crH_ - 4;
                        plot(ax_crop, [4 4+sb5], [sb_y5 sb_y5], 'w-', 'LineWidth', 2);
                        text(ax_crop, 4+sb5/2, sb_y5-6, '5 \mum', 'Color','w', ...
                            'FontSize', 6, 'HorizontalAlignment','center');
                        hold(ax_crop, 'off');
                    else
                        text(ax_crop, 0.5, 0.5, '?', 'HorizontalAlignment','center');
                        axis(ax_crop, 'off');
                    end
                else
                    axis(ax_crop, 'off');
                end
                title(ax_crop, sprintf('ROI#%02d', rid), 'FontSize', 7, 'Interpreter','none');

                % Col 2: breath + dF/F full trace (yyaxis)
                dff_roi_st = dFF(:, rid);
                ca_ev_st   = caData(rid).spike_idx(:);
                ca_ev_st   = ca_ev_st(ca_ev_st >= 1 & ca_ev_st <= T);

                ax_tr = axes(fig_st, 'Units','centimeters', ...
                    'Position', [trace_x row_y trace_w_cm row_h]);

                % Left axis: dF/F + spike dots
                yyaxis(ax_tr, 'left');
                hold(ax_tr, 'on');
                plot(ax_tr, t_full_st, dff_roi_st, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 0.4);
                if ~isempty(ca_ev_st)
                    plot(ax_tr, t_full_st(ca_ev_st), dff_roi_st(ca_ev_st), ...
                        'r.', 'MarkerSize', 4);
                end
                hold(ax_tr, 'off');
                dff_mn = min(dff_roi_st); dff_mx = max(dff_roi_st);
                dff_range = dff_mx - dff_mn + eps;
                dff_pad = dff_range * 0.5;   % 50% padding each side → dFF in middle 50%
                ylim(ax_tr, [dff_mn - dff_pad, dff_mx + dff_pad]);
                ylabel(ax_tr, 'dF/F', 'FontSize', 5);
                ax_tr.YColor = [0.2 0.7 0.2];

                % Right axis: breathing (no labels)
                yyaxis(ax_tr, 'right');
                plot(ax_tr, t_full_st, breath_dlc, 'k-', 'LineWidth', 0.3);
                set(ax_tr, 'YTick', [], 'YColor', 'none');

                xlim(ax_tr, [0, t_full_st(end)]);
                if jj_st == roi_end
                    xlabel(ax_tr, 'Time (s)');
                else
                    set(ax_tr, 'XTickLabel', []);
                end
                box(ax_tr, 'on');
                set(ax_tr, 'FontSize', 5);

                % Col 3: spike vs non-spike dF/F overlay (square)
                nCyc_st = numel(insp_valid) - 1;
                if nCyc_st >= 1
                    cyc_traces_st = cell(nCyc_st, 1);
                    cyc_times_st  = cell(nCyc_st, 1);
                    has_spike_st  = false(nCyc_st, 1);
                    for c_st = 1:nCyc_st
                        i1_ = insp_valid(c_st); i2_ = insp_valid(c_st + 1);
                        cyc_traces_st{c_st} = dff_roi_st(i1_:i2_);
                        cyc_times_st{c_st}  = (0:(i2_ - i1_))' / fs;
                        has_spike_st(c_st)  = any(ca_ev_st >= i1_ & ca_ev_st < i2_);
                    end

                    ax_spk = axes(fig_st, 'Units','centimeters', ...
                        'Position', [spk_x row_y spk_w_cm row_h]);
                    hold(ax_spk, 'on');
                    for c_st = find(~has_spike_st)'
                        plot(ax_spk, cyc_times_st{c_st}, cyc_traces_st{c_st}, ...
                            'Color', [0 0 0 0.3], 'LineWidth', 0.5);
                    end
                    for c_st = find(has_spike_st)'
                        plot(ax_spk, cyc_times_st{c_st}, cyc_traces_st{c_st}, ...
                            'Color', [1 0 0 0.3], 'LineWidth', 0.5);
                    end
                    axis(ax_spk, 'square');
                    hold(ax_spk, 'off');
                    if jj_st == roi_end
                        xlabel(ax_spk, 'Time (s)');
                    else
                        set(ax_spk, 'XTickLabel', []);
                    end
                    ylabel(ax_spk, 'dF/F');
                    box(ax_spk, 'on');
                    set(ax_spk, 'FontSize', 5);
                end
            end
            drawnow;

            starterPdf = fullfile(outDir_st, sprintf('FOV%02d_%s_starter_p%d.pdf', fov_id, sessName, iPage));
            exportgraphics(fig_st, starterPdf, 'ContentType','vector', 'BackgroundColor','none');
            fprintf('  Starter page %d saved: %s\n', iPage, starterPdf);
            sessStarterPdfs{end+1} = starterPdf; %#ok<SAGROW>
            close(fig_st);
        end

        % --- Summary figure: 3 rows x 1 col ---
        summ_sq = 5.0;   % slightly larger for summary panels
        summ_ml = 2.0;
        summ_mr = 0.8;
        summ_mb = 1.4;
        summ_mt = 1.8;
        summ_vg = 1.6;

        summ_fig_w = summ_ml + summ_sq + summ_mr;
        summ_fig_h = summ_mb + 3*summ_sq + 2*summ_vg + summ_mt;

        summ_row3_y = summ_mb;
        summ_row2_y = summ_row3_y + summ_sq + summ_vg;
        summ_row1_y = summ_row2_y + summ_sq + summ_vg;

        figS = figure('Color','w', 'Visible','off', ...
            'Name', sprintf('FOV#%02d | %s | Summary', fov_id, sessName), ...
            'Units','centimeters', 'Position', [1 0.5 summ_fig_w summ_fig_h], ...
            'PaperUnits','centimeters', 'PaperSize', [summ_fig_w summ_fig_h], ...
            'PaperPosition', [0 0 summ_fig_w summ_fig_h]);
        set(figS, 'DefaultAxesFontSize', 8, 'DefaultTextFontSize', 9);

        summLabel = sprintf('FOV#%02d | %s | Summary', fov_id, sessName);
        annotation(figS, 'textbox', [0.02 0.975 0.96 0.02], ...
            'String', summLabel, ...
            'EdgeColor','none', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'Interpreter','none', ...
            'FontName','Arial', 'FontSize', 9, 'Color', [0 0 0]);

        % ---- Summary Row 1: Overlaid breathing cycles ----
        axS1 = axes(figS, 'Units','centimeters', ...
            'Position', [summ_ml summ_row1_y summ_sq summ_sq]);
        hold(axS1, 'on');
        plot(axS1, phase_full, cycle_dup', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.5);
        plot(axS1, phase_full, mean_dup, 'k', 'LineWidth', 2);
        xline(axS1, 0, 'r--', 'LineWidth', 1);
        xline(axS1, -2*pi, 'r--', 'LineWidth', 1);
        xline(axS1, 2*pi, 'r--', 'LineWidth', 1);
        hold(axS1, 'off');
        xlim(axS1, [-2*pi, 2*pi]);
        set(axS1, 'XTick', [-2*pi -pi 0 pi 2*pi], ...
            'XTickLabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
        xlabel(axS1, 'Phase (rad)');
        ylabel(axS1, 'Breathing (z-scored)');
        title(axS1, sprintf('%d cycles overlaid', nValidCyc), 'Interpreter','none');
        box(axS1, 'on');
        pbaspect(axS1, [1 1 1]);

        % ---- Summary Row 2: Overlaid cosine fits ----
        axS2 = axes(figS, 'Units','centimeters', ...
            'Position', [summ_ml summ_row2_y summ_sq summ_sq]);
        hold(axS2, 'on');
        for jj = 1:numel(summ_use)
            kk_r = summ_use(jj);
            plot(axS2, summ_th_fit{kk_r}, summ_fit_ext{kk_r}, ...
                'Color', summ_col(kk_r,:), 'LineWidth', 1.5);
        end
        hold(axS2, 'off');
        xlim(axS2, [0, 4*pi]);
        set(axS2, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
            'XTickLabel', {'0', '\pi', '2\pi', '3\pi', '4\pi'});
        xlabel(axS2, 'Phase (rad)');
        ylabel(axS2, 'Spikes/cycle');
        title(axS2, 'Cosine fits (all spiking ROIs)', 'Interpreter','none');
        box(axS2, 'on');
        pbaspect(axS2, [1 1 1]);

        % ---- Summary Row 3: Overlaid polar ----
        figS_pos = figS.Position;
        polar_summ_norm = [summ_ml/figS_pos(3), summ_row3_y/figS_pos(4), ...
                           summ_sq/figS_pos(3), summ_sq/figS_pos(4)];
        axS3 = polaraxes(figS, 'Position', polar_summ_norm);
        hold(axS3, 'on');
        % Confidence circle
        th_circ_s = linspace(0, 2*pi, 360);
        polarplot(axS3, th_circ_s, repmat(confC_val, 1, 360), 'k--', 'LineWidth', 1.2);
        % Breathing waveform
        %polarplot(axS3, [phase_1cyc_p, 2*pi], [breath_norm_polar, breath_norm_polar(1)], 'k-', 'LineWidth', 0.8);
        % Each spiking ROI
        for jj = 1:numel(summ_use)
            kk_r = summ_use(jj);
            c_r  = summ_col(kk_r,:);
            % Radial error bar
            polarplot(axS3, [summ_th_k(kk_r) summ_th_k(kk_r)], ...
                [summ_rlo(kk_r) summ_rhi(kk_r)], '-', 'Color', c_r, 'LineWidth', 2);
            % Phase arc
            phi_arc_s = linspace(summ_th_k(kk_r) - summ_dphi(kk_r), ...
                                 summ_th_k(kk_r) + summ_dphi(kk_r), 40);
            polarplot(axS3, phi_arc_s, summ_r_k(kk_r)*ones(1,40), ...
                '-', 'Color', c_r, 'LineWidth', 1.5);
            % Dot
            polarplot(axS3, summ_th_k(kk_r), summ_r_k(kk_r), 'o', ...
                'MarkerFaceColor', c_r, 'MarkerEdgeColor','k', 'MarkerSize', 3);
        end
        axS3.RLim              = [0, 1];
        axS3.ThetaZeroLocation = 'top';
        axS3.ThetaDir          = 'clockwise';
        axS3.RAxisLocation     = 180;
        axS3.FontSize          = 8;
        hold(axS3, 'off');

        drawnow;

        % ---- Save summary ----
        outDir = fullfile(folderPath, 'breath_trig_combined');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        summFigPath = fullfile(outDir, ...
            sprintf('FOV%02d_%s_summary.fig', fov_id, sessName));
        savefig(figS, summFigPath);

        summPdfPath = fullfile(outDir, ...
            sprintf('FOV%02d_%s_summary.pdf', fov_id, sessName));
        exportgraphics(figS, summPdfPath, 'ContentType','vector', 'BackgroundColor','none');
        fprintf('  Summary saved: %s\n', summPdfPath);

        sessSummPdf = summPdfPath;
        close(figS);

        %% ---- Per-FOV merged PDF ----
        if mergePerFOV && ~isempty(sessRoiPdfs)
            % Sort ROI PDFs numerically
            roiNums = nan(size(sessRoiPdfs));
            for ii = 1:numel(sessRoiPdfs)
                tok = regexp(sessRoiPdfs{ii}, 'ROI(\d+)\.pdf$', 'tokens', 'once');
                if ~isempty(tok)
                    roiNums(ii) = str2double(tok{1});
                end
            end
            [~, ord] = sort(roiNums);
            sessRoiPdfs = sessRoiPdfs(ord);

            % Assemble PDF list: p0 + starter pages + ROI pages + summary
            pdfList = {};
            if exist('p0Pdf','var') && ~isempty(p0Pdf) && isfile(p0Pdf)
                pdfList{end+1} = p0Pdf;
            end
            pdfList = [pdfList, sessStarterPdfs, sessRoiPdfs];
            if ~isempty(sessSummPdf) && exist(sessSummPdf, 'file') == 2
                pdfList{end+1} = sessSummPdf;
            end

            fovMergedPdf = fullfile(outDir, ...
                sprintf('FOV%02d_%s_breath_trig_combined_ALL.pdf', fov_id, sessName));

            if exist(fovMergedPdf, 'file') == 2
                delete(fovMergedPdf);
            end

            append_pdfs_RZ(fovMergedPdf, pdfList{:});
            fprintf('  Merged FOV PDF: %s\n', fovMergedPdf);

            allFovMergedPdfs{end+1} = fovMergedPdf; %#ok<SAGROW>
        end

        nProcessed = nProcessed + 1;

    catch ME
        nFailed = nFailed + 1;
        failedSessions{end+1} = sessName; %#ok<SAGROW>
        fprintf('ERROR in %s:\n  %s\n', sessName, ME.message);
        for me_k = 1:numel(ME.stack)
            fprintf('  in %s (line %d)\n', ME.stack(me_k).name, ME.stack(me_k).line);
        end
    end
end  % end session loop

%% 4 ── Giant summary across all FOVs
if mergeAllFOV && ~isempty(allFovMergedPdfs)
    giantPdf = fullfile(inputPath, 'breath_trig_combined_GIANT_SUMMARY.pdf');

    if exist(giantPdf, 'file') == 2
        delete(giantPdf);
    end

    append_pdfs_RZ(giantPdf, allFovMergedPdfs{:});
    fprintf('\nGiant all-FOV summary saved:\n  %s\n', giantPdf);
end

%% 5 ── Final report
fprintf('\n====== BATCH COMPLETE ======\n');
fprintf('Processed: %d | Skipped: %d | Failed: %d\n', nProcessed, nSkipped, nFailed);

if ~isempty(allFovMergedPdfs)
    fprintf('Per-FOV merged PDFs created: %d\n', numel(allFovMergedPdfs));
end

if nFailed > 0
    fprintf('Failed sessions:\n');
    for f = 1:numel(failedSessions)
        fprintf('  - %s\n', failedSessions{f});
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function txt = polarTip(info)
% POLARTIP  Datacursor callback for polar plots with ROI metadata.
    src = info.Target;
    ud  = [];
    if isprop(src,'UserData'), ud = src.UserData; end
    th = info.Position(1);
    r  = info.Position(2);
    if isstruct(ud) && isfield(ud,'roi')
        txt = sprintf('ROI %d\nPhase = %.3g rad (%.1f%c) %c %.3g rad\nCoh = %.3g  [%.3g, %.3g]', ...
            ud.roi, th, rad2deg(th), char(176), char(177), ud.phistd, r, ud.rlo, ud.rhi);
    else
        txt = sprintf('Phase = %.3g rad (%.1f%c)\nCoh = %.3g', th, rad2deg(th), char(176), r);
    end
end

function append_pdfs_RZ(output, varargin)
% APPEND_PDFS_RZ  Concatenate multiple PDFs into one, silently.
% Includes a dummy last-page hack to avoid Ghostscript dropping
% the final input PDF in some installations.

    if nargin < 2
        error('append_pdfs_RZ:NotEnoughInputs', ...
              'Usage: append_pdfs_RZ(output, input1, input2, ...)');
    end

    % Remove empties / non-existing files
    keep = false(size(varargin));
    for i = 1:numel(varargin)
        keep(i) = ischar(varargin{i}) || isstring(varargin{i});
        if keep(i)
            keep(i) = exist(varargin{i}, 'file') == 2;
        end
    end
    inputFiles = varargin(keep);

    if isempty(inputFiles)
        error('append_pdfs_RZ:NoValidInputs', 'No valid input PDFs found.');
    end

    % Determine append vs create
    if exist(output, 'file') == 2
        tmpOutput  = [tempname '.pdf'];
        inputFiles = [{output} inputFiles];
    else
        tmpOutput  = output;
    end

    % ----- Create dummy page to protect the last real file -----
    dummyPDF = [tempname '_dummy.pdf'];
    hf = figure('Visible','off');
    axis off;
    set(hf, 'PaperPosition', [0 0 8.5 11]);
    print(hf, dummyPDF, '-dpdf', '-fillpage');
    close(hf);

    inputFiles = [inputFiles {dummyPDF}];
    % -----------------------------------------------------------

    % Ghostscript command file
    cmdfile = [tempname '.txt'];
    fh = fopen(cmdfile, 'w');

    fprintf(fh, ...
        '-q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="%s" -f', ...
        tmpOutput);

    for k = 1:numel(inputFiles)
        fprintf(fh, ' "%s"', inputFiles{k});
    end
    fclose(fh);

    % ---- Silent Ghostscript call ----
    try
        evalc('ghostscript(''@"'' + string(cmdfile) + ''"'' )');
    catch ME
        if exist(cmdfile,'file'), delete(cmdfile); end
        if exist(dummyPDF,'file'), delete(dummyPDF); end
        error('append_pdfs_RZ:GhostscriptFailed', ...
            'Ghostscript PDF append failed: %s', ME.message);
    end
    % ---------------------------------

    % Cleanup
    if exist(cmdfile,'file'), delete(cmdfile); end
    if exist(dummyPDF,'file'), delete(dummyPDF); end

    % Move temp into final if needed
    if ~strcmp(tmpOutput, output)
        movefile(tmpOutput, output, 'f');
    end
end
