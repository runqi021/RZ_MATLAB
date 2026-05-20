%% breath_time_summary_260329.m
%  Experiment-level TIME-domain summary: breath-triggered analysis
%  pooled across all FOVs.
%
%  Auto-detect mode:
%    - inputPath contains ca_spike_data.mat -> sub mode (single folder)
%    - Otherwise -> master mode (discover all subfolders with ca_spike_data.mat)
%
%  3 figures:
%    1) Peak-triggered breathing trace overlay (all sessions)
%    2) dF/F PETH heatmap — all spiking ROIs (sorted by peak latency)
%    3) dF/F PETH heatmap — sig ROIs only (sorted by peak latency)
%
%  Dependencies:
%    - helper.dFF_RZ
%    - Chronux (mtspectrumc, coherencyc)
%    - detect_session_fps.m

clear; close all; clc;

%% ========================= USER PARAMETERS =========================
inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx         = [];             % session indices to skip (master mode)
fps_img         = 30;             % imaging frame rate (Hz)
nDrop           = 30;             % frames to drop from start
TW              = 3;              % multitaper time-bandwidth (K = 2*TW-1 tapers)
f_breath_search = [1.5, 3];      % Hz (search band for breathing peak + coherence)
win_sec         = 0.5;            % half-window for PETH (seconds around breath peak)
prcLim          = [0.1 99.9];    % percentile clamp for dF/F heatmap
fmin            = 0.01;           % Hz (PSD lower bound)
fmax            = 14;             % Hz (PSD upper bound; [] -> fs/2)
% ====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

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

%% 2 -- Add Chronux to path
chronuxDir = fullfile(fileparts(mfilename('fullpath')), 'chronux_2_12');
if isfolder(chronuxDir) && ~contains(path, 'chronux')
    addpath(genpath(chronuxDir));
    fprintf('Added Chronux to path: %s\n', chronuxDir);
end

%% 3 -- Preallocate master collectors
all_breath_snip = [];       % [nSnippets x nTimeBins] — pooled breath snippets
all_breath_sess = [];       % session index per snippet

% Per-ROI PETH collector
all_roi_peth = struct('sessIdx',{}, 'roiId',{}, 'mean_peth',{}, ...
    'is_sig',{}, 'label',{}, 'dff_snips',{}, 'spike_times_per_peak',{}, ...
    'spike_t_dedup',{}, 'dt_nearest_per_peak',{});

% Session names and validity
sessNames = cell(nSess, 1);
validSess = false(nSess, 1);

% confC values per session
confC_vals = nan(nSess, 1);

% Track win_frames per session (for validation; should be constant if fs constant)
win_frames_per_sess = nan(nSess, 1);

%% 4 -- Main session loop (data collection)
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

        SAMload = load(SAM_path);
        F_raw = SAMload.F;
        F_raw(1:min(nDrop, size(F_raw,1)), :) = [];
        dFFout = helper.dFF_RZ(F_raw, 'FPS', fs);
        dFF = double(dFFout.dFF);     % [T x N]
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
        fs_breath = fs;  % breathing camera is frame-triggered 1:1 by calcium

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
        breath_dlc = breath_dlc - mean(breath_dlc);  % zero-mean; per-snippet [0,1] norm below

        nB = numel(breath_dlc);
        if nB >= T, breath_dlc = breath_dlc(1:T);
        else,       breath_dlc(end+1:T) = 0; end

        %% ---- Breathing PSD + FWHM (Chronux) ----
        params_breath.Fs     = fs_breath;
        params_breath.tapers = [TW, 2*TW - 1];
        params_breath.pad    = 0;
        params_breath.fpass  = [fmin, min(fmax_use, fs_breath/2)];
        params_breath.err    = [2, 0.05];
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
        fprintf('Breathing peak: %.3f Hz  FWHM: [%.3f, %.3f] Hz\n', ...
            f_breath_peak, f_fwhm(1), f_fwhm(2));

        %% ---- Compute PETH window ----
        b_frames     = sort(insp_onsets(:));
        win_frames   = round(win_sec * fs);

        % Common time grid (defined once, all sessions resampled to this)
        if ~exist('t_w', 'var') || isempty(t_w)
            nW  = 2 * round(win_sec * fps_img) + 1;  % reference grid at fallback fps
            t_w = linspace(-win_sec, win_sec, nW)';
        end

        win_frames_per_sess(kk) = win_frames;

        % Valid breath peaks (must have full window + edge margin within data)
        edge_margin = max(win_frames, nDrop);
        insp_valid = b_frames(b_frames > edge_margin & b_frames <= T - edge_margin);
        nInsp      = numel(insp_valid);
        fprintf('PETH window: +/-%.3f s (%d frames @ %dHz), %d valid onsets of %d total\n', ...
            win_sec, win_frames, fs, nInsp, numel(b_frames));

        if nInsp < 1
            warning('No valid breath peaks for PETH, skipping %s.', sessName);
            continue;
        end

        %% ---- Extract breath snippets (per-cycle [0,1] normalisation) ----
        t_local = (-win_frames:win_frames)' / fs;  % this session's time grid
        breath_snip_sess = zeros(nInsp, nW);
        for e = 1:nInsp
            c    = insp_valid(e);
            snip = breath_dlc(c - win_frames : c + win_frames);
            mn   = min(snip);  mx = max(snip);
            if mx > mn, snip = (snip - mn) / (mx - mn); end
            if numel(snip) == nW
                breath_snip_sess(e, :) = snip;
            else
                breath_snip_sess(e, :) = interp1(t_local, snip, t_w, 'linear', 'extrap');
            end
        end

        % Accumulate
        all_breath_snip = [all_breath_snip; breath_snip_sess];   %#ok<AGROW>
        all_breath_sess = [all_breath_sess; repmat(kk, nInsp, 1)]; %#ok<AGROW>

        %% ---- Chronux coherence: confC + per-ROI significance ----
        params_coh.Fs     = fs;
        params_coh.tapers = [TW, 2*TW - 1];
        params_coh.pad    = 0;
        params_coh.fpass  = f_breath_search;
        params_coh.err    = [2, 0.01];

        % Resample breath to calcium fs if needed
        if abs(fs_breath - fs) > 0.1
            breath_rs = resample(breath_dlc, round(fs), round(fs_breath));
        else
            breath_rs = breath_dlc;
        end
        T_coh     = min(T, numel(breath_rs));
        breath_rs = breath_rs(1:T_coh);

        % confC: compute once per session
        sig0 = dFF(1:T_coh, roiSpk_id(1)) - mean(dFF(1:T_coh, roiSpk_id(1)));
        [~, ~, ~, ~, ~, ~, ~, confC_val] = coherencyc(breath_rs, sig0, params_coh);
        confC_vals(kk) = confC_val;

        %% ---- Per-ROI: PETH + significance + per-cycle data ----
        for jj = 1:nSpkROI
            roi_id = roiSpk_id(jj);
            sig_c  = dFF(:, roi_id) - mean(dFF(:, roi_id));

            % Coherence for significance
            [~, C_coh, phi_coh, ~, ~, ~, f_coh, ~, phistd_coh, Cerr_coh] = ...
                coherencyc(breath_rs, sig_c(1:T_coh), params_coh);
            f_coh = f_coh(:); C_coh = C_coh(:);

            mask_fwhm = f_coh >= f_fwhm(1) & f_coh <= f_fwhm(2);
            if ~any(mask_fwhm), mask_fwhm = true(size(f_coh)); end
            r_k    = mean(C_coh(mask_fwhm));
            is_sig = r_k > confC_val;

            % Extract per-cycle dF/F snippets for this ROI
            dff_roi = dFF(:, roi_id);
            dff_snips = zeros(nInsp, nW);
            for e = 1:nInsp
                c = insp_valid(e);
                raw_snip = dff_roi(c - win_frames : c + win_frames);
                if numel(raw_snip) == nW
                    dff_snips(e, :) = raw_snip;
                else
                    dff_snips(e, :) = interp1(t_local, raw_snip, t_w, 'linear', 'extrap');
                end
            end
            mean_peth = mean(dff_snips, 1, 'omitnan');

            % Spike positions relative to each breath peak (in seconds)
            ca_ev = caData(roi_id).spike_idx(:);
            ca_ev = ca_ev(ca_ev >= 1 & ca_ev <= T);
            spike_times_per_peak = cell(nInsp, 1);
            for e = 1:nInsp
                c = insp_valid(e);
                spk_in_win = ca_ev(ca_ev >= c - win_frames & ca_ev <= c + win_frames);
                if ~isempty(spk_in_win)
                    spike_times_per_peak{e} = (spk_in_win - c) / fs;
                end
            end

            % Deduplicated spike times: each spike assigned to nearest peak only
            spike_t_dedup = zeros(numel(ca_ev), 1);  % preallocate max
            n_dedup = 0;
            for e = 1:numel(ca_ev)
                [min_dist, nearest_idx] = min(abs(double(insp_valid) - double(ca_ev(e))));
                if min_dist <= win_frames
                    n_dedup = n_dedup + 1;
                    spike_t_dedup(n_dedup) = (ca_ev(e) - insp_valid(nearest_idx)) / fs;
                end
            end
            spike_t_dedup = spike_t_dedup(1:n_dedup);

            % Signed dt to nearest calcium spike per breath peak (for sorted heatmap)
            dt_nearest_per_peak = nan(nInsp, 1);
            for e = 1:nInsp
                bp = insp_valid(e);
                before = ca_ev(ca_ev < bp);
                after  = ca_ev(ca_ev >= bp);
                dt_bef = Inf; dt_aft = Inf;
                if ~isempty(before), dt_bef = (bp - before(end)) / fs; end
                if ~isempty(after),  dt_aft = (after(1) - bp) / fs; end
                if dt_bef <= dt_aft
                    dt_nearest_per_peak(e) = -dt_bef;
                else
                    dt_nearest_per_peak(e) = dt_aft;
                end
            end

            % Store
            idx_r = numel(all_roi_peth) + 1;
            all_roi_peth(idx_r).sessIdx             = kk;
            all_roi_peth(idx_r).roiId               = roi_id;
            all_roi_peth(idx_r).mean_peth            = mean_peth;
            all_roi_peth(idx_r).is_sig               = is_sig;
            all_roi_peth(idx_r).label                = sprintf('S%dR%d', kk, roi_id);
            all_roi_peth(idx_r).dff_snips            = dff_snips;
            all_roi_peth(idx_r).spike_times_per_peak = spike_times_per_peak;
            all_roi_peth(idx_r).spike_t_dedup        = spike_t_dedup;
            all_roi_peth(idx_r).dt_nearest_per_peak  = dt_nearest_per_peak;

            sig_str = 'SIG'; if ~is_sig, sig_str = 'n.s.'; end
            fprintf('  ROI#%02d : r=%.3f  %s\n', roi_id, r_k, sig_str);
        end

        validSess(kk) = true;
        fprintf('Session %d complete: %d breath peaks, %d spiking ROIs\n', ...
            kk, nInsp, nSpkROI);

    catch ME
        fprintf('ERROR in %s:\n  %s\n', sessName, ME.message);
        for me_k = 1:numel(ME.stack)
            fprintf('  in %s (line %d)\n', ME.stack(me_k).name, ME.stack(me_k).line);
        end
    end
end  % end session loop

%% 5 -- Post-loop validation
nValidSess = sum(validSess);
fprintf('\n========== Collection complete ==========\n');
fprintf('Valid sessions: %d / %d\n', nValidSess, nSess);

assert(nValidSess > 0, 'No valid sessions to plot.');

nTotalSnip = size(all_breath_snip, 1);
nTotalROIs = numel(all_roi_peth);
nSigROIs   = sum([all_roi_peth.is_sig]);

fprintf('Total breath snippets: %d\n', nTotalSnip);
fprintf('Total spiking ROIs:    %d\n', nTotalROIs);
fprintf('Significant ROIs:      %d\n', nSigROIs);

% Mean confC across valid sessions
confC_mean_val = mean(confC_vals(validSess), 'omitnan');
fprintf('Mean confC:            %.4f\n', confC_mean_val);

% t_w and nW already defined in session loop

%% 6 -- Sort ROIs by peak latency of mean PETH
peth_mat_all = zeros(nTotalROIs, nW);
labels_all   = cell(nTotalROIs, 1);
for rr = 1:nTotalROIs
    peth_mat_all(rr, :) = all_roi_peth(rr).mean_peth;
    labels_all{rr}      = all_roi_peth(rr).label;
end

[~, peak_idx_all] = max(peth_mat_all, [], 2);
peak_lat_all      = t_w(peak_idx_all);
[~, sort_ord_all] = sort(peak_lat_all, 'ascend');

sig_mask = [all_roi_peth.is_sig];
sig_idx  = find(sig_mask);

% Sort order for sig ROIs only
if ~isempty(sig_idx)
    peth_mat_sig = peth_mat_all(sig_idx, :);
    [~, peak_idx_sig] = max(peth_mat_sig, [], 2);
    peak_lat_sig      = t_w(peak_idx_sig);
    [~, sort_ord_sig_local] = sort(peak_lat_sig, 'ascend');
    sort_ord_sig = sig_idx(sort_ord_sig_local);   % indices into all_roi_peth
end

%% 7 -- Build per-cycle heatmap matrices + spike overlay structs

% --- Helper: stack per-cycle dF/F + spike ticks for a given ROI order ---
%     Returns: big_mat, spike_overlay struct array, row_labels_pos, row_labels_txt
build_cycle_stack = @(roi_order) local_build_stack(all_roi_peth, roi_order, nW);

[big_mat_all, spk_ov_all, rlp_all, rlt_all] = build_cycle_stack(sort_ord_all);

if ~isempty(sig_idx)
    [big_mat_sig, spk_ov_sig, rlp_sig, rlt_sig] = build_cycle_stack(sort_ord_sig);
end

%% 7b -- Pool all spike times for PSTH
all_spike_t = [];
for rr = 1:nTotalROIs
    stp = all_roi_peth(rr).spike_times_per_peak;
    for ee = 1:numel(stp)
        if ~isempty(stp{ee})
            all_spike_t = [all_spike_t; stp{ee}(:)]; %#ok<AGROW>
        end
    end
end
fprintf('Pooled spike times: %d\n', numel(all_spike_t));

%% 8 -- Single figure: 2 columns x 3 rows
sq       = 3.17;      % square side (cm)  (3.8 / 1.2)
margin_l = 2.0;
margin_r = 1.6;       % room for colorbar
margin_b = 1.6;
margin_t = 1.2;
vgap     = 2.2;
hgap     = 2.4;

fig_w = margin_l + 2*sq + hgap + margin_r;
fig_h = margin_b + 4*sq + 3*vgap + margin_t;

row4_y = margin_b;
row3_y = row4_y + sq + vgap;
row2_y = row3_y + sq + vgap;
row1_y = row2_y + sq + vgap;

col1_x = margin_l;
col2_x = col1_x + sq + hgap;

figure('Color', 'w', 'Name', 'Breath Time Summary', ...
    'Units', 'centimeters', 'Position', [1 0.5 fig_w fig_h]);

%% ---- Panel 1: Peak-triggered breath trace overlay (normalised [0,1]) ----
ax1 = axes('Units', 'centimeters', 'Position', [col1_x row1_y sq sq]);
hold(ax1, 'on');
plot(ax1, t_w, all_breath_snip', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
grand_mean_breath = mean(all_breath_snip, 1, 'omitnan');
plot(ax1, t_w, grand_mean_breath, 'k-', 'LineWidth', 2);
xline(ax1, 0, 'r--', 'LineWidth', 1);
hold(ax1, 'off');
xlim(ax1, [t_w(1), t_w(end)]);
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Breathing (norm)');
title(ax1, sprintf('Breath overlay (%d sess, %d peaks)', nValidSess, nTotalSnip));
box(ax1, 'on');
pbaspect(ax1, [1 1 1]);

%% ---- Panel 2: Per-cycle dF/F heatmap + spike raster — ALL spiking ROIs ----
ax2 = axes('Units', 'centimeters', 'Position', [col1_x row2_y sq sq]);
hold(ax2, 'on');
imagesc(ax2, t_w, 1:size(big_mat_all,1), big_mat_all);
set(ax2, 'YDir', 'reverse');
colormap(ax2, flipud(gray(256)));
cHigh_all = prctile(big_mat_all(:), prcLim(2));
caxis(ax2, [0, cHigh_all]);
% Spike raster overlay
for sp = 1:numel(spk_ov_all)
    r = spk_ov_all(sp).row;
    t_spk = spk_ov_all(sp).tval;
    plot(ax2, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
end
xline(ax2, 0, 'r--', 'LineWidth', 1);
hold(ax2, 'off');
xlim(ax2, [t_w(1), t_w(end)]);
ylim(ax2, [0.5, size(big_mat_all,1)+0.5]);
set(ax2, 'YTick', rlp_all, 'YTickLabel', rlt_all, 'FontSize', 5);
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'ROI');
title(ax2, sprintf('Per-cycle dF/F + spikes — all ROIs (%d)', nTotalROIs));
cb2 = colorbar(ax2); ylabel(cb2, 'dF/F');
set(ax2, 'Units', 'centimeters', 'Position', [col1_x row2_y sq sq]);
cb2.Units = 'centimeters'; cb2.Position = [col1_x+sq+0.1, row2_y, 0.15, sq];

%% ---- Panel 3: Per-cycle dF/F heatmap + spike raster — SIG ROIs only ----
ax3 = axes('Units', 'centimeters', 'Position', [col1_x row3_y sq sq]);
if nSigROIs > 0
    hold(ax3, 'on');
    imagesc(ax3, t_w, 1:size(big_mat_sig,1), big_mat_sig);
    set(ax3, 'YDir', 'reverse');
    colormap(ax3, flipud(gray(256)));
    cHigh_sig = prctile(big_mat_sig(:), prcLim(2));
    caxis(ax3, [0, cHigh_sig]);
    % Spike raster overlay
    for sp = 1:numel(spk_ov_sig)
        r = spk_ov_sig(sp).row;
        t_spk = spk_ov_sig(sp).tval;
        plot(ax3, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
    end
    xline(ax3, 0, 'r--', 'LineWidth', 1);
    hold(ax3, 'off');
    xlim(ax3, [t_w(1), t_w(end)]);
    ylim(ax3, [0.5, size(big_mat_sig,1)+0.5]);
    set(ax3, 'YTick', rlp_sig, 'YTickLabel', rlt_sig, 'FontSize', 5);
    xlabel(ax3, 'Time (s)'); ylabel(ax3, 'ROI');
    title(ax3, sprintf('Per-cycle dF/F + spikes — sig ROIs (%d/%d)', nSigROIs, nTotalROIs));
    cb3 = colorbar(ax3); ylabel(cb3, 'dF/F');
    set(ax3, 'Units', 'centimeters', 'Position', [col1_x row3_y sq sq]);
    cb3.Units = 'centimeters'; cb3.Position = [col1_x+sq+0.1, row3_y, 0.15, sq];
else
    text(ax3, 0.5, 0.5, 'No sig ROIs', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title(ax3, 'PETH — sig ROIs (0)');
end

%% ---- Panel 4: Mean PETH heatmap — ALL spiking ROIs ----
peth_sorted_all_mat = peth_mat_all(sort_ord_all, :);
labels_sorted_all   = labels_all(sort_ord_all);

ax4 = axes('Units', 'centimeters', 'Position', [col2_x row1_y sq sq]);
imagesc(ax4, t_w, 1:nTotalROIs, peth_sorted_all_mat);
set(ax4, 'YDir', 'reverse');
colormap(ax4, flipud(gray(256)));
cHigh4 = prctile(peth_sorted_all_mat(:), prcLim(2));
caxis(ax4, [0, cHigh4]);
hold(ax4, 'on'); xline(ax4, 0, 'r--', 'LineWidth', 1); hold(ax4, 'off');
if nTotalROIs <= 40
    set(ax4, 'YTick', 1:nTotalROIs, 'YTickLabel', labels_sorted_all, 'FontSize', 5);
else
    ts4 = max(1, round(nTotalROIs/20));
    yt4 = 1:ts4:nTotalROIs;
    set(ax4, 'YTick', yt4, 'YTickLabel', labels_sorted_all(yt4), 'FontSize', 5);
end
xlabel(ax4, 'Time (s)'); ylabel(ax4, 'ROI');
title(ax4, sprintf('Mean PETH — all ROIs (%d)', nTotalROIs));
cb4 = colorbar(ax4); ylabel(cb4, 'dF/F');
set(ax4, 'Units', 'centimeters', 'Position', [col2_x row1_y sq sq]);
cb4.Units = 'centimeters'; cb4.Position = [col2_x+sq+0.1, row1_y, 0.15, sq];

%% ---- Panel 5: Mean PETH heatmap — SIG ROIs only ----
ax5 = axes('Units', 'centimeters', 'Position', [col2_x row2_y sq sq]);
if nSigROIs > 0
    peth_sorted_sig_mat = peth_mat_all(sort_ord_sig, :);
    labels_sorted_sig   = labels_all(sort_ord_sig);

    imagesc(ax5, t_w, 1:nSigROIs, peth_sorted_sig_mat);
    set(ax5, 'YDir', 'reverse');
    colormap(ax5, flipud(gray(256)));
    cHigh5 = prctile(peth_sorted_sig_mat(:), prcLim(2));
    caxis(ax5, [0, cHigh5]);
    hold(ax5, 'on'); xline(ax5, 0, 'r--', 'LineWidth', 1); hold(ax5, 'off');
    if nSigROIs <= 40
        set(ax5, 'YTick', 1:nSigROIs, 'YTickLabel', labels_sorted_sig, 'FontSize', 5);
    else
        ts5 = max(1, round(nSigROIs/20));
        yt5 = 1:ts5:nSigROIs;
        set(ax5, 'YTick', yt5, 'YTickLabel', labels_sorted_sig(yt5), 'FontSize', 5);
    end
    xlabel(ax5, 'Time (s)'); ylabel(ax5, 'ROI');
    title(ax5, sprintf('Mean PETH — sig ROIs (%d/%d)', nSigROIs, nTotalROIs));
    cb5 = colorbar(ax5); ylabel(cb5, 'dF/F');
    set(ax5, 'Units', 'centimeters', 'Position', [col2_x row2_y sq sq]);
    cb5.Units = 'centimeters'; cb5.Position = [col2_x+sq+0.1, row2_y, 0.15, sq];
else
    text(ax5, 0.5, 0.5, 'No sig ROIs', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title(ax5, 'Mean PETH — sig ROIs (0)');
end

%% ---- Panel 6: Spike rate PSTH (spikes/sec) ----
% Pool all spike times and compute firing rate histogram
dt_bin     = t_w(2) - t_w(1);
peth_edges = [t_w(:) - dt_bin/2; t_w(end) + dt_bin/2];
peth_counts = histcounts(all_spike_t, peth_edges);
peth_rate   = peth_counts / (nTotalSnip * dt_bin);  % spikes/sec

ax6 = axes('Units', 'centimeters', 'Position', [col2_x row3_y sq sq]);
hold(ax6, 'on');
bar(ax6, t_w, peth_rate, 1, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
xline(ax6, 0, 'r--', 'LineWidth', 1);
hold(ax6, 'off');
xlim(ax6, [t_w(1), t_w(end)]);
xlabel(ax6, 'Time (s)');
ylabel(ax6, 'Spike rate (Hz)');
title(ax6, sprintf('PSTH (%d ROIs, %d peaks)', nTotalROIs, nTotalSnip));
box(ax6, 'on');
pbaspect(ax6, [1 1 1]);

%% ---- Panel 7: Sorted per-cycle heatmap (by dt_nearest_ca) — ALL spiking ROIs ----
[sorted_all, spk_sort_all] = local_build_stack_sorted(all_roi_peth, sort_ord_all, nW);

ax7 = axes('Units', 'centimeters', 'Position', [col1_x row4_y sq sq]);
hold(ax7, 'on');
imagesc(ax7, t_w, 1:size(sorted_all,1), sorted_all);
set(ax7, 'YDir', 'reverse');
colormap(ax7, flipud(gray(256)));
cHigh7 = prctile(sorted_all(:), prcLim(2));
caxis(ax7, [0, cHigh7]);
% Spike raster overlay
for sp = 1:numel(spk_sort_all)
    r = spk_sort_all(sp).row;
    t_spk = spk_sort_all(sp).tval;
    plot(ax7, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
end
xline(ax7, 0, 'r--', 'LineWidth', 1);
hold(ax7, 'off');
xlim(ax7, [t_w(1), t_w(end)]);
ylim(ax7, [0.5, size(sorted_all,1)+0.5]);
xlabel(ax7, 'Time (s)'); ylabel(ax7, 'Cycle (sorted by dt)');
title(ax7, sprintf('Sorted by dt_{nearest} — all ROIs (%d)', nTotalROIs));
cb7 = colorbar(ax7); ylabel(cb7, 'dF/F');
set(ax7, 'Units', 'centimeters', 'Position', [col1_x row4_y sq sq]);
cb7.Units = 'centimeters'; cb7.Position = [col1_x+sq+0.1, row4_y, 0.15, sq];

%% ---- Panel 8: Sorted per-cycle heatmap (by dt_nearest_ca) — SIG ROIs only ----
ax8 = axes('Units', 'centimeters', 'Position', [col2_x row4_y sq sq]);
if nSigROIs > 0
    [sorted_sig, spk_sort_sig] = local_build_stack_sorted(all_roi_peth, sort_ord_sig, nW);

    hold(ax8, 'on');
    imagesc(ax8, t_w, 1:size(sorted_sig,1), sorted_sig);
    set(ax8, 'YDir', 'reverse');
    colormap(ax8, flipud(gray(256)));
    cHigh8 = prctile(sorted_sig(:), prcLim(2));
    caxis(ax8, [0, cHigh8]);
    % Spike raster overlay
    for sp = 1:numel(spk_sort_sig)
        r = spk_sort_sig(sp).row;
        t_spk = spk_sort_sig(sp).tval;
        plot(ax8, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
    end
    xline(ax8, 0, 'r--', 'LineWidth', 1);
    hold(ax8, 'off');
    xlim(ax8, [t_w(1), t_w(end)]);
    ylim(ax8, [0.5, size(sorted_sig,1)+0.5]);
    xlabel(ax8, 'Time (s)'); ylabel(ax8, 'Cycle (sorted by dt)');
    title(ax8, sprintf('Sorted by dt_{nearest} — sig ROIs (%d/%d)', nSigROIs, nTotalROIs));
    cb8 = colorbar(ax8); ylabel(cb8, 'dF/F');
    set(ax8, 'Units', 'centimeters', 'Position', [col2_x row4_y sq sq]);
    cb8.Units = 'centimeters'; cb8.Position = [col2_x+sq+0.1, row4_y, 0.15, sq];
else
    text(ax8, 0.5, 0.5, 'No sig ROIs', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title(ax8, 'Sorted by dt_{nearest} — sig ROIs (0)');
end

fprintf('\nDone. %d sessions, %d snippets, %d ROIs (%d sig).\n', ...
    nValidSess, nTotalSnip, nTotalROIs, nSigROIs);

%% ===== LOCAL FUNCTION: build per-cycle stacked matrix =====
function [big_mat, spike_overlay, row_labels_pos, row_labels_txt] = ...
        local_build_stack(all_roi_peth, roi_order, nW)
% Stack per-cycle dF/F snippets for ROIs in roi_order, with 2-row NaN gap
% between ROIs.  Collect spike overlay positions as struct array.

    big_mat        = [];
    spike_overlay  = struct('row',{}, 'tval',{});
    row_labels_pos = [];
    row_labels_txt = {};

    for ii = 1:numel(roi_order)
        rd       = all_roi_peth(roi_order(ii));
        dff_cyc  = rd.dff_snips;           % [nPeaks x nW]
        nPeaks   = size(dff_cyc, 1);
        row_start = size(big_mat, 1) + 1;
        big_mat   = [big_mat; dff_cyc; nan(2, nW)]; %#ok<AGROW>

        % Collect spike tick positions
        for e = 1:nPeaks
            spk_t = rd.spike_times_per_peak{e};
            if ~isempty(spk_t)
                row_in = row_start + e - 1;
                for s = 1:numel(spk_t)
                    spike_overlay(end+1).row  = row_in;     %#ok<AGROW>
                    spike_overlay(end).tval   = spk_t(s);
                end
            end
        end

        row_labels_pos(end+1) = row_start + nPeaks/2 - 0.5; %#ok<AGROW>
        row_labels_txt{end+1} = rd.label;                    %#ok<AGROW>
    end

    % Remove trailing 2-row NaN gap
    if size(big_mat, 1) >= 2
        big_mat = big_mat(1:end-2, :);
    end
end

%% ===== LOCAL FUNCTION: build sorted per-cycle stack (by dt_nearest_ca) =====
function [sorted_mat, spike_overlay] = local_build_stack_sorted(all_roi_peth, roi_order, nW)
% Pool all cycles from ROIs in roi_order, sort by dt_nearest_ca.
% No ROI grouping, no NaN gaps. Each row = one breath cycle.
% Returns sorted dF/F matrix and spike overlay struct.

    all_snips = [];
    all_spk_t = {};
    all_dt    = [];

    for ii = 1:numel(roi_order)
        rd = all_roi_peth(roi_order(ii));
        nP = size(rd.dff_snips, 1);
        all_snips = [all_snips; rd.dff_snips];          %#ok<AGROW>
        all_dt    = [all_dt; rd.dt_nearest_per_peak];   %#ok<AGROW>
        for e = 1:nP
            all_spk_t{end+1,1} = rd.spike_times_per_peak{e}; %#ok<AGROW>
        end
    end

    % Sort by dt_nearest ascending
    dt_sort = all_dt;
    dt_sort(isnan(dt_sort)) = Inf;
    [~, sort_idx] = sort(dt_sort, 'ascend');
    sorted_mat = all_snips(sort_idx, :);

    % Build spike overlay for sorted order
    spike_overlay = struct('row',{}, 'tval',{});
    for r = 1:numel(sort_idx)
        spk = all_spk_t{sort_idx(r)};
        if ~isempty(spk)
            for s = 1:numel(spk)
                spike_overlay(end+1).row = r;           %#ok<AGROW>
                spike_overlay(end).tval = spk(s);
            end
        end
    end
end
