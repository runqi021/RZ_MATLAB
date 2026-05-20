% breath_coherence_dFF_vs_spike_260410.m
%  Side-by-side comparison of breathing coherence computed two ways:
%    (1) dF/F waveform vs breathing  (coherencyc  -- continuous vs continuous)
%    (2) Calcium spike train vs breathing (coherencycpb -- continuous vs binned point process)
%
%  Produces a 3-panel figure:
%    Left polar  -- dF/F coherence  (r = magnitude, theta = phase)
%    Right polar -- spike train coherence (same format)
%    Bottom      -- scatter of dF/F coherence vs spike coherence (per ROI)
%
%  Session discovery & data loading follow the same pattern as
%  breath_combined_summary_260329.m.
%
%  Dependencies:
%    - helper.dFF_RZ
%    - Chronux (mtspectrumc, coherencyc, coherencycpb)
%    - detect_session_fps.m

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
% ======================================================================

%% ========================= USER PARAMETERS ===========================
inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx         = [];             % session indices to skip (master mode)
fps_img         = 30;             % imaging frame rate (Hz) -- used as fallback only
nDrop           = 30;             % frames dropped from start (by cpSAM)
TW              = 5;              % multitaper time-bandwidth product (K = 2*TW-1 tapers)
f_breath_search = [1.5 2];       % Hz search band for breathing peak; FWHM narrows for coherence
fwhm_factor     = 1;              % coherence band = fwhm_factor x FWHM (centered on peak)
BaselineWinSec  = 20;             % dF/F sliding median baseline window (seconds)
fmin            = 0.01;           % Hz (PSD lower bound)
fmax            = 14;             % Hz (PSD upper bound; [] -> fs/2)
doSave          = true;           % true = save PDF + .fig, false = just plot
% ======================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% 1 -- Auto-detect mode & discover sessions
if exist(fullfile(inputPath, 'ca_spike_data.mat'), 'file') == 2
    fprintf('Sub mode: processing single folder.\n');
    allMat = dir(fullfile(inputPath, 'ca_spike_data.mat'));
    isMasterMode = false;
else
    allMat = dir(fullfile(inputPath, '**', 'ca_spike_data.mat'));
    isMasterMode = true;

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

%% 2 -- Preallocate master collectors
% dF/F coherence collectors
all_th_dff     = [];
all_r_dff      = [];
all_rlo_dff    = [];
all_rhi_dff    = [];
all_sig_dff    = [];    % logical: significant?

% Spike train coherence collectors
all_th_spk     = [];
all_r_spk      = [];
all_rlo_spk    = [];
all_rhi_spk    = [];
all_sig_spk    = [];    % logical: significant?
all_zerosp     = [];    % logical: ROI had zero spikes in coherence window?

% Shared metadata
all_roi_sess   = [];    % session index per ROI
all_roi_id     = [];    % ROI index within session

% Per-session confidence thresholds
confC_dff_vals = nan(nSess, 1);
confC_spk_vals = nan(nSess, 1);

% Session metadata
sessNames      = cell(nSess, 1);
validSess      = false(nSess, 1);

%% 3 -- Main session loop
for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, sessName] = fileparts(folderPath);
    sessNames{kk} = sessName;

    fprintf('\n========== [%d/%d] %s ==========\n', kk, nSess, sessName);

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
            continue;
        end
        if isempty(sam_hits)
            warning('No SAM output found in %s, skipping.', sessName);
            continue;
        end
        if isempty(bp_file)
            warning('No breath_peak_data.mat found in %s, skipping.', sessName);
            continue;
        end

        dlc_csv  = fullfile(folderPath, csv_hits(1).name);
        SAM_path = fullfile(folderPath, sam_hits(1).name);

        fprintf('DLC CSV : %s\n', csv_hits(1).name);
        fprintf('SAM     : %s\n', sam_hits(1).name);

        %% ---- Auto-detect fps, load SAM -> dFF ----
        fs = detect_session_fps(folderPath, fps_img);
        fs_breath = fs;  % breathing camera is frame-triggered 1:1

        SAMload = load(SAM_path);
        F_raw   = SAMload.F;
        dFFout  = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF     = double(dFFout.dFF);     % [T x N]
        [T, N_roi] = size(dFF);
        if isempty(fmax), fmax_use = fs / 2; else, fmax_use = fmax; end

        fprintf('dFF: [%d x %d] (T x N_roi)\n', T, N_roi);

        %% ---- Load calcium spikes ----
        caLoad    = load(fullfile(folderPath, 'ca_spike_data.mat'));
        ifSpk     = caLoad.ifSpike;
        roiSpk_id = find(ifSpk == 1);
        caData    = caLoad.roi_spikes;
        nSpkROI   = numel(roiSpk_id);

        if nSpkROI == 0
            fprintf('No spiking ROIs, skipping session.\n');
            continue;
        end
        fprintf('Spiking ROIs: %d of %d\n', nSpkROI, N_roi);

        %% ---- Load breath peaks ----
        breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));
        insp_onsets = sort(breath_peak.insp_onset_idx(:));
        fp          = breath_peak.findpeak_params;

        % Adjust for dropped frames
        insp_onsets(insp_onsets <= nDrop) = [];
        insp_onsets = insp_onsets - nDrop;

        if numel(insp_onsets) < 2
            warning('Fewer than 2 breath peaks after drop, skipping %s.', sessName);
            continue;
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
        breath_dlc = breath_dlc - mean(breath_dlc);

        nB = numel(breath_dlc);
        if nB >= T, breath_dlc = breath_dlc(1:T);
        else,       breath_dlc(end+1:T) = 0; end

        %% ---- Breathing PSD + FWHM (Chronux) ----
        params_breath.Fs     = fs_breath;
        params_breath.tapers = [TW, 2*TW - 1];
        params_breath.pad    = 0;
        params_breath.fpass  = [fmin, min(fmax_use, fs_breath/2)];
        params_breath.err    = [2, 0.01];
        [Sb, fbC, ~] = mtspectrumc(breath_dlc, params_breath);
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
        bw_scaled  = diff(f_fwhm) * fwhm_factor;
        f_coh_band = [f_breath_peak - bw_scaled/2, f_breath_peak + bw_scaled/2];
        f_coh_band = [max(f_coh_band(1), f_breath_search(1)), min(f_coh_band(2), f_breath_search(2))];

        fprintf('Breathing peak: %.3f Hz  FWHM: [%.3f, %.3f]  coh band (x%.1f): [%.3f, %.3f] Hz\n', ...
            f_breath_peak, f_fwhm(1), f_fwhm(2), fwhm_factor, f_coh_band(1), f_coh_band(2));

        %% ---- Coherence parameters ----
        params_coh.Fs     = fs;
        params_coh.tapers = [TW, 2*TW - 1];
        params_coh.pad    = 0;
        params_coh.fpass  = f_coh_band;
        params_coh.err    = [2, 0.01];

        % Resample breathing to match imaging fs if needed
        if abs(fs_breath - fs) > 0.1
            breath_rs = resample(breath_dlc, round(fs), round(fs_breath));
        else
            breath_rs = breath_dlc;
        end
        T_coh     = min(T, numel(breath_rs));
        breath_rs = breath_rs(1:T_coh);

        %% ---- Compute confC once per session (from first spiking ROI) ----
        sig0 = dFF(1:T_coh, roiSpk_id(1)) - mean(dFF(1:T_coh, roiSpk_id(1)));
        [~, ~, ~, ~, ~, ~, ~, confC_dff] = coherencyc(breath_rs, sig0, params_coh);
        confC_dff_vals(kk) = confC_dff;

        % For spike coherence confC: compute from first ROI with spikes
        % NOTE: coherencycpb only computes confC when nargout >= 9 (unlike
        % coherencyc which does it at nargout >= 8), so request 9 outputs.
        spk0 = caData(roiSpk_id(1)).spike_train;
        spk0 = spk0(1:min(numel(spk0), T_coh));
        if numel(spk0) < T_coh
            spk0(end+1:T_coh, 1) = 0;
        end
        [~, ~, ~, ~, ~, ~, ~, confC_spk, ~] = coherencycpb(breath_rs, spk0, params_coh);
        confC_spk_vals(kk) = confC_spk;

        fprintf('confC: dF/F=%.4f  spike=%.4f\n', confC_dff, confC_spk);

        %% ---- Per-ROI coherence loop ----
        validSess(kk) = true;

        % Preallocate session-level arrays
        sess_th_dff  = nan(nSpkROI, 1);
        sess_r_dff   = nan(nSpkROI, 1);
        sess_rlo_dff = nan(nSpkROI, 1);
        sess_rhi_dff = nan(nSpkROI, 1);
        sess_sig_dff = false(nSpkROI, 1);

        sess_th_spk  = nan(nSpkROI, 1);
        sess_r_spk   = nan(nSpkROI, 1);
        sess_rlo_spk = nan(nSpkROI, 1);
        sess_rhi_spk = nan(nSpkROI, 1);
        sess_sig_spk = false(nSpkROI, 1);
        sess_zerosp  = false(nSpkROI, 1);

        for jj = 1:nSpkROI
            roi_id = roiSpk_id(jj);

            %% -- dF/F coherence (coherencyc) --
            sig_c = dFF(:, roi_id) - mean(dFF(:, roi_id));
            [~, C_dff, phi_dff, ~, ~, ~, f_dff, ~, phistd_dff, Cerr_dff] = ...
                coherencyc(breath_rs, sig_c(1:T_coh), params_coh);
            f_dff = f_dff(:); C_dff = C_dff(:); phi_dff = phi_dff(:);

            mask_band = f_dff >= f_coh_band(1) & f_dff <= f_coh_band(2);
            if ~any(mask_band), mask_band = true(size(f_dff)); end

            sess_r_dff(jj)   = mean(C_dff(mask_band));
            sess_sig_dff(jj) = sess_r_dff(jj) > confC_dff;
            sess_th_dff(jj)  = angle(mean(exp(1i * (-phi_dff(mask_band)))));
            sess_rlo_dff(jj) = max(0,   mean(Cerr_dff(1, mask_band)));
            sess_rhi_dff(jj) = min(1,   mean(Cerr_dff(2, mask_band)));

            %% -- Spike train coherence (coherencycpb) --
            spk_train = caData(roi_id).spike_train;
            % Trim or pad spike_train to match T_coh
            if numel(spk_train) > T_coh
                spk_train = spk_train(1:T_coh);
            elseif numel(spk_train) < T_coh
                spk_train(end+1:T_coh, 1) = 0;
            end
            spk_train = spk_train(:);  % ensure column

            n_spk_in_window = sum(spk_train);
            if n_spk_in_window == 0
                % Zero spikes: coherence undefined
                sess_r_spk(jj)   = 0;
                sess_th_spk(jj)  = NaN;
                sess_sig_spk(jj) = false;
                sess_rlo_spk(jj) = 0;
                sess_rhi_spk(jj) = 0;
                sess_zerosp(jj)  = true;
                fprintf('  ROI %3d: zero spikes in coherence window, skipping spike coherence.\n', roi_id);
                continue;
            end

            [C_spk, phi_spk, ~, ~, ~, f_spk, zerosp_flag, ~, phistd_spk, Cerr_spk] = ...
                coherencycpb(breath_rs, spk_train, params_coh);
            f_spk = f_spk(:); C_spk = C_spk(:); phi_spk = phi_spk(:);

            if zerosp_flag
                % Chronux detected zero spikes (should not happen given check above, but guard)
                sess_r_spk(jj)   = 0;
                sess_th_spk(jj)  = NaN;
                sess_sig_spk(jj) = false;
                sess_rlo_spk(jj) = 0;
                sess_rhi_spk(jj) = 0;
                sess_zerosp(jj)  = true;
                fprintf('  ROI %3d: coherencycpb reported zerosp, skipping.\n', roi_id);
                continue;
            end

            mask_band_spk = f_spk >= f_coh_band(1) & f_spk <= f_coh_band(2);
            if ~any(mask_band_spk), mask_band_spk = true(size(f_spk)); end

            sess_r_spk(jj)   = mean(C_spk(mask_band_spk));
            sess_sig_spk(jj) = sess_r_spk(jj) > confC_spk;
            sess_th_spk(jj)  = angle(mean(exp(1i * (-phi_spk(mask_band_spk)))));
            sess_rlo_spk(jj) = max(0,   mean(Cerr_spk(1, mask_band_spk)));
            sess_rhi_spk(jj) = min(1,   mean(Cerr_spk(2, mask_band_spk)));
            sess_zerosp(jj)  = false;
        end

        %% ---- Append to master collectors ----
        all_th_dff     = [all_th_dff;     sess_th_dff];      %#ok<AGROW>
        all_r_dff      = [all_r_dff;      sess_r_dff];       %#ok<AGROW>
        all_rlo_dff    = [all_rlo_dff;    sess_rlo_dff];     %#ok<AGROW>
        all_rhi_dff    = [all_rhi_dff;    sess_rhi_dff];     %#ok<AGROW>
        all_sig_dff    = [all_sig_dff;    sess_sig_dff];     %#ok<AGROW>

        all_th_spk     = [all_th_spk;     sess_th_spk];     %#ok<AGROW>
        all_r_spk      = [all_r_spk;      sess_r_spk];      %#ok<AGROW>
        all_rlo_spk    = [all_rlo_spk;    sess_rlo_spk];    %#ok<AGROW>
        all_rhi_spk    = [all_rhi_spk;    sess_rhi_spk];    %#ok<AGROW>
        all_sig_spk    = [all_sig_spk;    sess_sig_spk];    %#ok<AGROW>
        all_zerosp     = [all_zerosp;     sess_zerosp];      %#ok<AGROW>

        all_roi_sess   = [all_roi_sess;   repmat(kk, nSpkROI, 1)]; %#ok<AGROW>
        all_roi_id     = [all_roi_id;     roiSpk_id(:)];    %#ok<AGROW>

        fprintf('Session %d: %d/%d dFF-sig, %d/%d spk-sig, %d zerosp\n', ...
            kk, sum(sess_sig_dff), nSpkROI, sum(sess_sig_spk), nSpkROI, sum(sess_zerosp));

    catch ME
        warning('Session %d (%s) failed: %s', kk, sessName, ME.message);
        fprintf('%s\n', getReport(ME, 'extended'));
        continue;
    end
end

%% 4 -- Validate we collected data
nTotalROIs = numel(all_r_dff);
assert(nTotalROIs > 0, 'No ROIs collected across any session.');
fprintf('\n===== TOTAL: %d ROIs from %d valid sessions =====\n', ...
    nTotalROIs, sum(validSess));

%% 5 -- Compute summary statistics
confC_dff_mean = mean(confC_dff_vals(validSess), 'omitnan');
confC_spk_mean = mean(confC_spk_vals(validSess), 'omitnan');

nSig_dff = sum(all_sig_dff);
nSig_spk = sum(all_sig_spk);
nZero    = sum(all_zerosp);

fprintf('dF/F coherence:  %d/%d significant (confC=%.4f)\n', nSig_dff, nTotalROIs, confC_dff_mean);
fprintf('Spike coherence: %d/%d significant (confC=%.4f), %d zerosp\n', nSig_spk, nTotalROIs, confC_spk_mean, nZero);

%% 6 -- Session colormap
validIdx   = find(validSess);
nValidSess = numel(validIdx);
if nValidSess <= 7
    sessColors = lines(max(nValidSess, 1));
else
    sessColors = turbo(nValidSess);
end

% Map each ROI's session index -> color row
sessColorMap = zeros(nSess, 3);
for ii = 1:nValidSess
    sessColorMap(validIdx(ii), :) = sessColors(ii, :);
end

%% 7 -- Figure: 2 polars + 1 scatter
fig_w = 18;  % cm
fig_h = 16;  % cm
fig = figure('Units','centimeters','Position',[2 2 fig_w fig_h], ...
    'Color','w', 'Name','Coherence: dFF vs Spike Train');

% ---- Polar 1: dF/F coherence ----
ax_pol_dff = polaraxes(fig, 'Units','normalized', ...
    'Position',[0.02 0.42 0.46 0.52]);
hold(ax_pol_dff, 'on');

th_circ = linspace(0, 2*pi, 360);
polarplot(ax_pol_dff, th_circ, repmat(confC_dff_mean, 1, 360), ...
    'k--', 'LineWidth', 1, 'DisplayName', sprintf('confC=%.3f', confC_dff_mean));

for rr = 1:nTotalROIs
    if all_sig_dff(rr)
        c_r = sessColorMap(all_roi_sess(rr), :);
    else
        c_r = [0 0 0];
    end
    polarplot(ax_pol_dff, all_th_dff(rr), all_r_dff(rr), 'o', ...
        'MarkerFaceColor', c_r, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 3, 'HandleVisibility','off');
end

ax_pol_dff.RLim = [0 1];
ax_pol_dff.ThetaZeroLocation = 'top';
ax_pol_dff.ThetaDir = 'clockwise';
ax_pol_dff.RAxisLocation = 180;
ax_pol_dff.FontSize = 7;
hold(ax_pol_dff, 'off');
title(ax_pol_dff, sprintf('dF/F Coherence (%d/%d sig)', nSig_dff, nTotalROIs), ...
    'Interpreter','none', 'FontSize', 8);

% ---- Polar 2: Spike train coherence ----
ax_pol_spk = polaraxes(fig, 'Units','normalized', ...
    'Position',[0.52 0.42 0.46 0.52]);
hold(ax_pol_spk, 'on');

polarplot(ax_pol_spk, th_circ, repmat(confC_spk_mean, 1, 360), ...
    'k--', 'LineWidth', 1, 'DisplayName', sprintf('confC=%.3f', confC_spk_mean));

for rr = 1:nTotalROIs
    if all_zerosp(rr)
        % Mark zero-spike ROIs with an X
        polarplot(ax_pol_spk, 0, 0, 'x', 'Color', [0.6 0.6 0.6], ...
            'MarkerSize', 4, 'HandleVisibility','off');
        continue;
    end
    if all_sig_spk(rr)
        c_r = sessColorMap(all_roi_sess(rr), :);
    else
        c_r = [0 0 0];
    end
    polarplot(ax_pol_spk, all_th_spk(rr), all_r_spk(rr), 'o', ...
        'MarkerFaceColor', c_r, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 3, 'HandleVisibility','off');
end

ax_pol_spk.RLim = [0 1];
ax_pol_spk.ThetaZeroLocation = 'top';
ax_pol_spk.ThetaDir = 'clockwise';
ax_pol_spk.RAxisLocation = 180;
ax_pol_spk.FontSize = 7;
hold(ax_pol_spk, 'off');
title(ax_pol_spk, sprintf('Spike Coherence (%d/%d sig, %d zerosp)', ...
    nSig_spk, nTotalROIs, nZero), 'Interpreter','none', 'FontSize', 8);

% ---- Scatter: dF/F coherence vs spike coherence ----
ax_scat = axes(fig, 'Units','normalized', 'Position',[0.12 0.06 0.35 0.30]);
hold(ax_scat, 'on');

% Only plot ROIs with valid spike coherence (exclude zerosp)
valid_mask = ~all_zerosp & ~isnan(all_th_spk);

for rr = 1:nTotalROIs
    if ~valid_mask(rr), continue; end
    is_either_sig = all_sig_dff(rr) || all_sig_spk(rr);
    if is_either_sig
        c_r = sessColorMap(all_roi_sess(rr), :);
        mfc = c_r;
    else
        c_r = [0 0 0];
        mfc = [0.7 0.7 0.7];
    end
    plot(ax_scat, all_r_dff(rr), all_r_spk(rr), 'o', ...
        'MarkerFaceColor', mfc, 'MarkerEdgeColor', c_r, 'MarkerSize', 4);
end

% Unity line
plot(ax_scat, [0 1], [0 1], 'k--', 'LineWidth', 0.5);

% Confidence thresholds
xline(ax_scat, confC_dff_mean, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(ax_scat, confC_spk_mean, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

xlim(ax_scat, [0 1]); ylim(ax_scat, [0 1]);
xlabel(ax_scat, 'dF/F Coherence (r)', 'FontSize', 8);
ylabel(ax_scat, 'Spike Train Coherence (r)', 'FontSize', 8);
title(ax_scat, 'dF/F vs Spike Coherence', 'Interpreter','none', 'FontSize', 8);
axis(ax_scat, 'square');
hold(ax_scat, 'off');

% Pearson correlation (valid ROIs only)
r_dff_valid = all_r_dff(valid_mask);
r_spk_valid = all_r_spk(valid_mask);
if numel(r_dff_valid) >= 3
    [rho, p_rho] = corr(r_dff_valid, r_spk_valid, 'Type', 'Pearson');
    text(ax_scat, 0.05, 0.95, sprintf('r=%.3f, p=%.1e\nn=%d', rho, p_rho, sum(valid_mask)), ...
        'Units','normalized', 'FontSize', 6, 'VerticalAlignment','top');
end

% ---- Session legend ----
ax_leg = axes(fig, 'Units','normalized', 'Position',[0.58 0.06 0.38 0.30], ...
    'Visible','off');
hold(ax_leg, 'on');
leg_entries = gobjects(nValidSess, 1);
for ii = 1:nValidSess
    sIdx = validIdx(ii);
    sName = strrep(sessNames{sIdx}, '_', ' ');
    nROI_sess = sum(all_roi_sess == sIdx);
    nSig_d = sum(all_sig_dff(all_roi_sess == sIdx));
    nSig_s = sum(all_sig_spk(all_roi_sess == sIdx));
    leg_entries(ii) = plot(ax_leg, NaN, NaN, 'o', ...
        'MarkerFaceColor', sessColors(ii,:), 'MarkerEdgeColor','k', ...
        'MarkerSize', 5, ...
        'DisplayName', sprintf('%s (n=%d, dFF:%d, spk:%d)', sName, nROI_sess, nSig_d, nSig_s));
end
hold(ax_leg, 'off');
legend(ax_leg, leg_entries, 'Location','northwest', 'FontSize', 5.5, ...
    'Box','off', 'Interpreter','none');

%% 8 -- Add super-title
[~, expName] = fileparts(char(inputPath));
sgtitle(fig, sprintf('Breathing Coherence: dF/F vs Spike Train  |  %s', ...
    strrep(expName, '_', ' ')), 'Interpreter','none', 'FontSize', 9);

%% 9 -- Save
if doSave
    outStem = fullfile(char(inputPath), 'breath_coherence_dFF_vs_spike');
    outPDF  = [outStem '.pdf'];
    outFig  = [outStem '.fig'];
    outMat  = [outStem '.mat'];

    exportgraphics(fig, outPDF, 'ContentType','vector', 'BackgroundColor','white');
    savefig(fig, outFig);

    % Save data for downstream use
    coh_results = struct();
    coh_results.all_th_dff    = all_th_dff;
    coh_results.all_r_dff     = all_r_dff;
    coh_results.all_rlo_dff   = all_rlo_dff;
    coh_results.all_rhi_dff   = all_rhi_dff;
    coh_results.all_sig_dff   = all_sig_dff;
    coh_results.all_th_spk    = all_th_spk;
    coh_results.all_r_spk     = all_r_spk;
    coh_results.all_rlo_spk   = all_rlo_spk;
    coh_results.all_rhi_spk   = all_rhi_spk;
    coh_results.all_sig_spk   = all_sig_spk;
    coh_results.all_zerosp    = all_zerosp;
    coh_results.all_roi_sess  = all_roi_sess;
    coh_results.all_roi_id    = all_roi_id;
    coh_results.confC_dff_vals = confC_dff_vals;
    coh_results.confC_spk_vals = confC_spk_vals;
    coh_results.confC_dff_mean = confC_dff_mean;
    coh_results.confC_spk_mean = confC_spk_mean;
    coh_results.sessNames      = sessNames;
    coh_results.validSess      = validSess;
    coh_results.f_breath_search = f_breath_search;
    coh_results.TW              = TW;
    coh_results.fwhm_factor     = fwhm_factor;
    coh_results.BaselineWinSec  = BaselineWinSec;
    coh_results.nDrop           = nDrop;
    save(outMat, 'coh_results', '-v7.3');

    fprintf('\nSaved:\n  %s\n  %s\n  %s\n', outPDF, outFig, outMat);
end

fprintf('\nDone.\n');
