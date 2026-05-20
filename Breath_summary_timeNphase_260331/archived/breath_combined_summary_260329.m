% breath_combined_summary_260329.m
%  Experiment-level combined TIME + PHASE summary figure.
%  Merges breath_time_summary_260329.m and breath_phase_summary_260328.m
%  into a single 5-column x 4-row figure (20 panels).
%
%  Layout:
%    Col1(time)         Col2(time)       Col3(phase)         Col4(phase)         Col5(phase)
%    R1 Breath overlay       Peak-trig avg    Breath overlay(phase)    Cycle-avg ALL       Polar coherence
%    R2 Per peak hmap ALL    Peak-trig sig    Per cycle hmap ALL       Cycle-avg SIG       Polar %active
%    R3 Per peak hmap SIG    Spike hist       Per cycle hmap SIG       Phase hist          Poincare
%    R4 Sorted peak ALL Sorted peak SIG  Sorted cycle ALL    Sorted cycle SIG    KDE
%
%  Auto-detect mode:
%    - inputPath contains ca_spike_data.mat -> sub mode (single folder)
%    - Otherwise -> master mode (discover all subfolders with ca_spike_data.mat)
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
TW              = 5;              % multitaper time-bandwidth (K = 2*TW-1 tapers)
f_breath_search = [1, 3.4];      % Hz (search band for breathing peak + coherence)
win_sec         = 0.5;            % half-window for peak-trig avg (seconds around breath peak)
prcLim          = [0.1 99.9];    % percentile clamp for dF/F heatmap
nPhaseBins      = 200;            % bins for phase-normalized breath cycles
fmin            = 0.01;           % Hz (PSD lower bound)
fmax            = 14;             % Hz (PSD upper bound; [] -> fs/2)
ibi_xlim        = [0 1];         % x/y limits for Poincare and KDE
insp_phase_lo   = 3*pi/2;        % inspiration start (rad)
insp_phase_hi   = 2*pi;          % inspiration end   (rad)
BaselineWinSec  = 20;            % dF/F sliding median baseline window (seconds)
doSave          = true;         % true = save PDF, false = just plot
% ====================================================================

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

%% 2 -- Add Chronux to path
chronuxDir = fullfile(fileparts(mfilename('fullpath')), 'chronux_2_12');
if isfolder(chronuxDir) && ~contains(path, 'chronux')
    addpath(genpath(chronuxDir));
    fprintf('Added Chronux to path: %s\n', chronuxDir);
end

%% 3 -- Preallocate master collectors

% --- TIME domain ---
all_breath_snip  = [];       % [nSnippets x nW]
all_breath_sess  = [];       % session index per snippet

% Combined per-ROI struct (time + phase fields)
all_roi_data = struct( ...
    'sessIdx',{}, 'roiId',{}, 'label',{}, 'is_sig',{}, ...
    'mean_peth',{}, 'dff_snips',{}, ...
    'spike_times_per_peak',{}, 'spike_t_dedup',{}, 'dt_nearest_per_peak',{}, ...
    'dff_phase',{}, 'nCycles',{}, 'spike_phase_per_cycle',{}, ...
    'th_k',{}, 'r_k',{}, 'pct_active',{});

% --- PHASE domain ---
all_cycle_mat   = [];       % [totalCycles x nPhaseBins]
all_cycle_sess  = [];
all_frame_phase = [];

% Polar / IBI collectors
all_th_k       = [];
all_r_k        = [];
all_rlo        = [];
all_rhi        = [];
all_dphi       = [];
all_roi_sess   = [];
all_pct_active = [];
all_is_sig_vec = [];     % logical vector aligned with all_th_k
all_spike_phase = {};
all_roi_id     = [];

ibi_all        = cell(nSess, 1);
confC_vals     = nan(nSess, 1);

% Spectra collectors
all_breath_psd   = cell(nSess, 1);
all_breath_freq  = cell(nSess, 1);
all_roi_psd      = {};
all_roi_psd_freq = {};

% Session metadata
sessNames      = cell(nSess, 1);
validSess      = false(nSess, 1);
win_frames_per_sess = nan(nSess, 1);
sess2fov       = zeros(nSess, 1);
nValidFov      = 0;

%% 4 -- Main session loop (merged data collection)
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
        F_raw = SAMload.F;
        %F_raw(1:min(nDrop, size(F_raw,1)), :) = [];
        dFFout = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
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

        % Adjust for dropped frames
        insp_onsets(insp_onsets <= nDrop) = [];
        insp_onsets = insp_onsets - nDrop;

        if numel(insp_onsets) < 2
            warning('Fewer than 2 breath peaks after drop, skipping %s.', sessName);
            continue;
        end

        % IBI
        peak_t_sec = insp_onsets(:) / fs_breath;
        ibi_all{kk} = diff(peak_t_sec);

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
        params_breath.err    = [2, 0.05];
        breath_z = breath_dlc / std(breath_dlc);
        [Sb, fbC, ~] = mtspectrumc(breath_z, params_breath);
        fbC = fbC(:); Sb = Sb(:);
        all_breath_psd{kk}  = Sb;  % linear power (convert to dB at plot time)
        all_breath_freq{kk} = fbC;

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

        %% ===== TIME COLLECTION =====

        %% ---- Peak-trig window ----
        b_frames     = sort(insp_onsets(:));
        win_frames   = round(win_sec * fs);

        % Common time grid (defined once, all sessions resampled to this)
        if ~exist('t_w', 'var') || isempty(t_w)
            nW  = 2 * round(win_sec * fps_img) + 1;
            t_w = linspace(-win_sec, win_sec, nW)';
        end

        win_frames_per_sess(kk) = win_frames;

        % Valid breath peaks (must have full window + edge margin)
        edge_margin = max(win_frames, nDrop);
        insp_valid = b_frames(b_frames > edge_margin & b_frames <= T - edge_margin);
        nInsp      = numel(insp_valid);
        fprintf('Peak-trig window: +/-%.3f s (%d frames @ %dHz), %d valid onsets of %d total\n', ...
            win_sec, win_frames, fs, nInsp, numel(b_frames));

        if nInsp < 1
            warning('No valid breath peaks for peak-trig avg, skipping %s.', sessName);
            continue;
        end

        %% ---- Extract breath snippets (per-cycle [0,1] normalisation) ----
        t_local = (-win_frames:win_frames)' / fs;
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

        all_breath_snip = [all_breath_snip; breath_snip_sess];   %#ok<AGROW>
        all_breath_sess = [all_breath_sess; repmat(kk, nInsp, 1)]; %#ok<AGROW>

        %% ===== PHASE COLLECTION =====

        %% ---- Phase-normalize each breathing cycle ----
        phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
        phase_1cyc = phase_1cyc(1:end-1);   % [0, 2pi)
        nPeaks_    = numel(insp_valid);
        nCyc_      = nPeaks_ - 1;
        cycle_mat  = nan(nCyc_, nPhaseBins);
        for c = 1:nCyc_
            i1 = insp_valid(c);
            i2 = insp_valid(c + 1);
            if i1 >= 1 && i2 <= T
                chunk = breath_dlc(i1:i2);
                orig_phase = linspace(0, 2*pi, numel(chunk));
                cycle_mat(c, :) = interp1(orig_phase, chunk, phase_1cyc, 'linear');
            end
        end
        valid_cyc = ~all(isnan(cycle_mat), 2);
        cycle_mat = cycle_mat(valid_cyc, :);
        nValidCyc = size(cycle_mat, 1);

        all_cycle_mat  = [all_cycle_mat;  cycle_mat];               %#ok<AGROW>
        all_cycle_sess = [all_cycle_sess; repmat(kk, nValidCyc, 1)]; %#ok<AGROW>

        fprintf('Breathing cycles: %d valid of %d total\n', nValidCyc, nCyc_);

        % Collect instantaneous phase of every frame in every cycle
        for c = 1:nCyc_
            bp1_f = insp_valid(c); bp2_f = insp_valid(c+1);
            L_f = bp2_f - bp1_f;
            if bp1_f >= 1 && bp2_f <= T && L_f > 0
                all_frame_phase = [all_frame_phase; ...
                    2*pi * ((0:L_f-1)' + 0.5) / L_f]; %#ok<AGROW>
            end
        end

        %% ===== COHERENCE (shared time+phase) =====
        params_coh.Fs     = fs;
        params_coh.tapers = [TW, 2*TW - 1];
        params_coh.pad    = 0;
        params_coh.fpass  = f_breath_search;
        params_coh.err    = [2, 0.01];

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

        %% ===== PER-ROI LOOP (time + phase combined) =====
        nValidFov = nValidFov + 1;
        sess2fov(kk) = nValidFov;

        % Trim raw F to match dFF length (dFF_RZ drops frames from front)
        F_trim = double(F_raw(end-T+1:end, :));

        for jj = 1:nSpkROI
            roi_id = roiSpk_id(jj);
            sig_c  = dFF(:, roi_id) - mean(dFF(:, roi_id));

            % --- Coherence for significance + polar data ---
            [~, C_coh, phi_coh, ~, ~, ~, f_coh, ~, phistd_coh, Cerr_coh] = ...
                coherencyc(breath_rs, sig_c(1:T_coh), params_coh);
            f_coh = f_coh(:); C_coh = C_coh(:); phi_coh = phi_coh(:);

            mask_fwhm = f_coh >= f_fwhm(1) & f_coh <= f_fwhm(2);
            if ~any(mask_fwhm), mask_fwhm = true(size(f_coh)); end
            r_k    = mean(C_coh(mask_fwhm));
            is_sig = r_k > confC_val;
            th_k   = angle(mean(exp(1i * (-phi_coh(mask_fwhm)))));
            dphi_k = mean(phistd_coh(mask_fwhm));
            rlo_k  = max(0,   mean(Cerr_coh(1, mask_fwhm)));
            rhi_k  = min(1,   mean(Cerr_coh(2, mask_fwhm)));

            % --- ROI spectrum (on derivative of raw F) ---
            params_roi_psd = params_breath;
            params_roi_psd.Fs = fs;
            params_roi_psd.fpass = [fmin, min(fmax_use, fs/2)];
            dsig = diff(F_trim(1:T_coh, roi_id)) * fs;
            [S_roi, f_roi] = mtspectrumc(dsig, params_roi_psd);
            idx_psd = numel(all_roi_psd) + 1;
            all_roi_psd{idx_psd}      = 10*log10(S_roi(:));
            all_roi_psd_freq{idx_psd} = f_roi(:);

            % --- TIME: per-cycle dF/F snippets ---
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

            % --- TIME: spike times relative to each breath peak ---
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

            % --- TIME: deduplicated spike times ---
            spike_t_dedup = zeros(numel(ca_ev), 1);
            n_dedup = 0;
            for e = 1:numel(ca_ev)
                [min_dist, nearest_idx] = min(abs(double(insp_valid) - double(ca_ev(e))));
                if min_dist <= win_frames
                    n_dedup = n_dedup + 1;
                    spike_t_dedup(n_dedup) = (ca_ev(e) - insp_valid(nearest_idx)) / fs;
                end
            end
            spike_t_dedup = spike_t_dedup(1:n_dedup);

            % --- TIME: signed dt to nearest calcium spike per breath peak ---
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

            % --- PHASE: %active cycles ---
            ca_ev_roi = ca_ev;
            nCyc_roi  = numel(insp_valid) - 1;
            active_count = 0;
            for cc = 1:nCyc_roi
                if any(ca_ev_roi >= insp_valid(cc) & ca_ev_roi < insp_valid(cc+1))
                    active_count = active_count + 1;
                end
            end
            pct_active = 100 * active_count / max(nCyc_roi, 1);

            % --- PHASE: spike phases for polar ---
            dt_l = nan(numel(ca_ev_roi), 1);
            dt_n = nan(numel(ca_ev_roi), 1);
            for ee = 1:numel(ca_ev_roi)
                c = ca_ev_roi(ee);
                bef = insp_valid(insp_valid <= c);
                aft = insp_valid(insp_valid >  c);
                if ~isempty(bef), dt_l(ee) = (c - bef(end)) / fs; end
                if ~isempty(aft), dt_n(ee) = (aft(1) - c)   / fs; end
            end
            kk_valid_sp = ~isnan(dt_l) & ~isnan(dt_n);
            spike_phase_roi = 2*pi * (dt_l(kk_valid_sp) + 0.5/fs) ./ ...
                (dt_l(kk_valid_sp) + dt_n(kk_valid_sp));

            % --- PHASE: dF/F per cycle ---
            nCyc_roi_ = numel(insp_valid) - 1;
            dff_phase_roi = nan(nCyc_roi_, nPhaseBins);
            for c_ = 1:nCyc_roi_
                i1_ = insp_valid(c_); i2_ = insp_valid(c_+1);
                if i1_ >= 1 && i2_ <= T
                    chunk_ = dff_roi(i1_:i2_);
                    orig_ph_ = linspace(0, 2*pi, numel(chunk_));
                    dff_phase_roi(c_,:) = interp1(orig_ph_, chunk_, phase_1cyc, 'linear');
                end
            end

            % --- PHASE: spike phases per cycle ---
            spk_per_cyc = cell(nCyc_roi_, 1);
            for c_ = 1:nCyc_roi_
                bp1_ = insp_valid(c_); bp2_ = insp_valid(c_+1);
                spk_ = ca_ev_roi(ca_ev_roi >= bp1_ & ca_ev_roi < bp2_);
                if ~isempty(spk_)
                    spk_per_cyc{c_} = 2*pi * (spk_ - bp1_ + 0.5) / (bp2_ - bp1_);
                end
            end

            % --- Store into combined struct ---
            idx_r = numel(all_roi_data) + 1;
            all_roi_data(idx_r).sessIdx              = kk;
            all_roi_data(idx_r).roiId                = roi_id;
            all_roi_data(idx_r).label                = sprintf('F#%02dR#%02d', sess2fov(kk), roi_id);
            all_roi_data(idx_r).is_sig               = is_sig;
            all_roi_data(idx_r).mean_peth             = mean_peth;
            all_roi_data(idx_r).dff_snips             = dff_snips;
            all_roi_data(idx_r).spike_times_per_peak  = spike_times_per_peak;
            all_roi_data(idx_r).spike_t_dedup         = spike_t_dedup;
            all_roi_data(idx_r).dt_nearest_per_peak   = dt_nearest_per_peak;
            all_roi_data(idx_r).dff_phase             = dff_phase_roi;
            all_roi_data(idx_r).nCycles               = nCyc_roi_;
            all_roi_data(idx_r).spike_phase_per_cycle = spk_per_cyc;
            all_roi_data(idx_r).th_k                  = th_k;
            all_roi_data(idx_r).r_k                   = r_k;
            all_roi_data(idx_r).pct_active            = pct_active;

            % --- Accumulate polar vectors ---
            all_th_k       = [all_th_k;       th_k];           %#ok<AGROW>
            all_r_k        = [all_r_k;        r_k];            %#ok<AGROW>
            all_rlo        = [all_rlo;        rlo_k];          %#ok<AGROW>
            all_rhi        = [all_rhi;        rhi_k];          %#ok<AGROW>
            all_dphi       = [all_dphi;       dphi_k];         %#ok<AGROW>
            all_roi_sess   = [all_roi_sess;   kk];             %#ok<AGROW>
            all_pct_active = [all_pct_active; pct_active];     %#ok<AGROW>
            all_is_sig_vec = [all_is_sig_vec; is_sig];         %#ok<AGROW>
            all_spike_phase = [all_spike_phase; {spike_phase_roi}]; %#ok<AGROW>
            all_roi_id     = [all_roi_id;     roi_id];         %#ok<AGROW>

            sig_str = 'SIG'; if ~is_sig, sig_str = 'n.s.'; end
            fprintf('  ROI#%02d : r=%.3f  th=%.2f rad  %s  %.1f%% active\n', ...
                roi_id, r_k, th_k, sig_str, pct_active);
        end

        validSess(kk) = true;
        fprintf('Session %d complete: %d breath peaks, %d cycles, %d spiking ROIs\n', ...
            kk, nInsp, nValidCyc, nSpkROI);

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

nTotalSnip   = size(all_breath_snip, 1);
nTotalROIs   = numel(all_roi_data);
nSigROIs     = sum([all_roi_data.is_sig]);
nTotalCycles = size(all_cycle_mat, 1);

fprintf('Total breath snippets: %d\n', nTotalSnip);
fprintf('Total breathing cycles:%d\n', nTotalCycles);
fprintf('Total spiking ROIs:    %d\n', nTotalROIs);
fprintf('Significant ROIs:      %d\n', nSigROIs);

confC_mean_val = mean(confC_vals(validSess), 'omitnan');
fprintf('Mean confC:            %.4f\n', confC_mean_val);

%% 6 -- Post-processing: build peak-trig avg matrix + unified alphabetical sort
peth_mat_all = zeros(nTotalROIs, nW);
labels_all   = cell(nTotalROIs, 1);
for rr = 1:nTotalROIs
    peth_mat_all(rr, :) = all_roi_data(rr).mean_peth;
    labels_all{rr}      = all_roi_data(rr).label;
end

% Unified alphabetical sort order for ALL unsorted plots
[~, sort_ord_all] = sort(labels_all);

sig_mask = [all_roi_data.is_sig];
sig_idx  = find(sig_mask);
nsig_idx = find(~sig_mask);
nNonSigROIs = numel(nsig_idx);

% Filter alphabetical order to sig/non-sig subsets
sort_ord_sig  = sort_ord_all(ismember(sort_ord_all, sig_idx));
sort_ord_nsig = sort_ord_all(ismember(sort_ord_all, nsig_idx));

fprintf('Non-sig ROIs:          %d\n', nNonSigROIs);

%% 7 -- TIME: build per-cycle stacks
[big_mat_all, spk_ov_all, rlp_all, rlt_all] = ...
    local_build_stack(all_roi_data, sort_ord_all, nW, 'dff_snips', 'spike_times_per_peak');

if nSigROIs > 0
    [big_mat_sig, spk_ov_sig, rlp_sig, rlt_sig] = ...
        local_build_stack(all_roi_data, sort_ord_sig, nW, 'dff_snips', 'spike_times_per_peak');
end

if nNonSigROIs > 0
    [big_mat_nsig, spk_ov_nsig, rlp_nsig, rlt_nsig] = ...
        local_build_stack(all_roi_data, sort_ord_nsig, nW, 'dff_snips', 'spike_times_per_peak');
end

%% 8 -- TIME: sorted stacks (by dt_nearest)
[sorted_all_t, spk_sort_all_t] = ...
    local_build_stack_sorted(all_roi_data, sort_ord_all, nW, ...
    'dff_snips', 'spike_times_per_peak', 'dt_nearest_per_peak');

if nSigROIs > 0
    [sorted_sig_t, spk_sort_sig_t] = ...
        local_build_stack_sorted(all_roi_data, sort_ord_sig, nW, ...
        'dff_snips', 'spike_times_per_peak', 'dt_nearest_per_peak');
end

if nNonSigROIs > 0
    [sorted_nsig_t, spk_sort_nsig_t] = ...
        local_build_stack_sorted(all_roi_data, sort_ord_nsig, nW, ...
        'dff_snips', 'spike_times_per_peak', 'dt_nearest_per_peak');
end

%% 9 -- TIME: pool spike times for spike hist (sig ROIs only)
all_spike_t = [];
for rr = 1:nTotalROIs
    if ~all_roi_data(rr).is_sig, continue; end
    stp = all_roi_data(rr).spike_times_per_peak;
    for ee = 1:numel(stp)
        if ~isempty(stp{ee})
            all_spike_t = [all_spike_t; stp{ee}(:)]; %#ok<AGROW>
        end
    end
end
fprintf('Pooled spike times (sig only): %d\n', numel(all_spike_t));

%% 10 -- PHASE post-processing
mean_cycle = mean(all_cycle_mat, 1, 'omitnan');
phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
phase_1cyc = phase_1cyc(1:end-1);
phase_full = [phase_1cyc, phase_1cyc + 2*pi];
cycle_dup  = [all_cycle_mat, all_cycle_mat];

% Normalize each cycle to [0,1] for phase overlay
for cc = 1:size(cycle_dup, 1)
    mn = min(cycle_dup(cc,:));
    mx = max(cycle_dup(cc,:));
    if mx > mn, cycle_dup(cc,:) = (cycle_dup(cc,:) - mn) / (mx - mn); end
end

%% 11 -- PHASE: build per-cycle phase stacks (grouped by ROI, alphabetical order)
mean_dff_phase_all = zeros(nTotalROIs, nPhaseBins);
for rr = 1:nTotalROIs
    valid_rows = ~all(isnan(all_roi_data(rr).dff_phase), 2);
    if any(valid_rows)
        mean_dff_phase_all(rr,:) = mean(all_roi_data(rr).dff_phase(valid_rows,:), 1, 'omitnan');
    end
end

% Per-cycle phase heatmap stacks (uses unified alphabetical sort_ord_*)
[ph_big_all, ph_spk_all, ph_rlp_all, ph_rlt_all] = ...
    local_build_stack_phase(all_roi_data, sort_ord_all, nPhaseBins);

if nSigROIs > 0
    [ph_big_sig, ph_spk_sig, ph_rlp_sig, ph_rlt_sig] = ...
        local_build_stack_phase(all_roi_data, sort_ord_sig, nPhaseBins);
end

%% 12 -- PHASE: cycle-avg heatmaps (one row per ROI, mean dff_phase, alphabetical order)
avg_phase_all = mean_dff_phase_all(sort_ord_all, :);
avg_labels_phase_all = labels_all(sort_ord_all);

if nSigROIs > 0
    avg_phase_sig = mean_dff_phase_all(sort_ord_sig, :);
    avg_labels_phase_sig = labels_all(sort_ord_sig);
end

%% 13 -- PHASE: sorted stacks (sort cycles by nearest spike phase)
[sorted_all_ph, spk_sort_all_ph] = ...
    local_build_stack_sorted_phase(all_roi_data, sort_ord_all, nPhaseBins);

if nSigROIs > 0
    [sorted_sig_ph, spk_sort_sig_ph] = ...
        local_build_stack_sorted_phase(all_roi_data, sort_ord_sig, nPhaseBins);
end

% Non-sig phase stacks
if nNonSigROIs > 0
    [ph_big_nsig, ph_spk_nsig, ph_rlp_nsig, ph_rlt_nsig] = ...
        local_build_stack_phase(all_roi_data, sort_ord_nsig, nPhaseBins);

    [sorted_nsig_ph, spk_sort_nsig_ph] = ...
        local_build_stack_sorted_phase(all_roi_data, sort_ord_nsig, nPhaseBins);

    avg_phase_nsig = mean_dff_phase_all(sort_ord_nsig, :);
    avg_labels_phase_nsig = labels_all(sort_ord_nsig);
end

%% 14 -- Phase spike data (sig ROIs only for histogram)
all_spike_ph_pool = {};   % sig only
for rr_ = 1:nTotalROIs
    if all_roi_data(rr_).is_sig
        spc_ = all_roi_data(rr_).spike_phase_per_cycle;
        for cc_ = 1:numel(spc_)
            if ~isempty(spc_{cc_})
                all_spike_ph_pool = [all_spike_ph_pool; {spc_{cc_}(:)}]; %#ok<AGROW>
            end
        end
    end
end
all_spike_ph_pool = vertcat(all_spike_ph_pool{:});
n_spk = numel(all_spike_ph_pool);

% Non-sig spike times (time domain) and phases
all_spike_t_nsig = [];
all_spike_ph_nsig = {};
for rr_ = 1:nTotalROIs
    if all_roi_data(rr_).is_sig, continue; end
    stp_ = all_roi_data(rr_).spike_times_per_peak;
    for ee_ = 1:numel(stp_)
        if ~isempty(stp_{ee_})
            all_spike_t_nsig = [all_spike_t_nsig; stp_{ee_}(:)]; %#ok<AGROW>
        end
    end
    spc_ = all_roi_data(rr_).spike_phase_per_cycle;
    for cc_ = 1:numel(spc_)
        if ~isempty(spc_{cc_})
            all_spike_ph_nsig = [all_spike_ph_nsig; {spc_{cc_}(:)}]; %#ok<AGROW>
        end
    end
end
if ~isempty(all_spike_ph_nsig)
    all_spike_ph_nsig = vertcat(all_spike_ph_nsig{:});
else
    all_spike_ph_nsig = [];
end

%% 15 -- Session colormap
validSessIdx = find(validSess);
cmap_sess = turbo(nValidSess);
sessColor = nan(nSess, 3);
for ii = 1:nValidSess
    sessColor(validSessIdx(ii), :) = cmap_sess(ii, :);
end

confC_mean = confC_mean_val;

% Breathing waveform for polar
breath_norm_polar = (mean_cycle - min(mean_cycle)) / (max(mean_cycle) - min(mean_cycle) + eps);

%% ===================================================================
%%  FIGURE: 5 columns x 4 rows = 20 panels
%% ===================================================================
sq_w     = 1.31 * 1.3 * 1.2 * 1.1;  % panel width  (cm)
sq_h     = 3.17 * 1.1;  % panel height (cm)
sq       = sq_w;
margin_l = 2.0;
margin_r = 1.6;
margin_b = 1.2;
margin_t = 1.0;
vgap     = 1.8;
hgap     = 2.0;

fig_w = margin_l + 8*sq_w + 7*hgap + margin_r;
fig_h = margin_b + 4*sq_h + 3*vgap + margin_t;

row4_y = margin_b;
row3_y = row4_y + sq_h + vgap;
row2_y = row3_y + sq_h + vgap;
row1_y = row2_y + sq_h + vgap;

col1_x  = margin_l;
col2_x  = col1_x  + sq_w + hgap;
col3_x  = col2_x  + sq_w + hgap;
col4_x  = col3_x  + sq_w + hgap;
col5_x  = col4_x  + sq_w + hgap;
col6_x  = col5_x  + sq_w + hgap;
col7_x  = col6_x  + sq_w + hgap;
col8_x  = col7_x  + sq_w + hgap;   % spectra

fig = figure('Color', 'w', 'Name', 'Breath Combined Summary', ...
    'Units', 'centimeters', 'Position', [0.5 0.5 fig_w fig_h], ...
    'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
    'PaperPosition', [0 0 fig_w fig_h]);
set(fig, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 7);

%% =====================================================================
%%  COLUMN 1 (TIME)
%% =====================================================================

%% ---- (R1,C1): Breath overlay ----
ax = axes(fig, 'Units','centimeters','Position',[col1_x row1_y sq_w sq_h]);
hold(ax,'on');
plot(ax, t_w, all_breath_snip', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
xline(ax, 0, 'k--', 'LineWidth', 1);
hold(ax,'off');
xlim(ax, [t_w(1), t_w(end)]);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Breathing (norm)');
th = title(ax, sprintf('Breath overlay\n%d cycles from %d FOVs', nTotalSnip, nValidSess));
th.Units = 'normalized'; th.Position(2) = 1.16;
box(ax,'on'); pbaspect(ax,[1 1 1]);

%% ---- (R2,C1): Per peak heatmap ALL ----
ax = axes(fig, 'Units','centimeters','Position',[col1_x row2_y sq_w sq_h]);
hold(ax,'on');
imagesc(ax, t_w, 1:size(big_mat_all,1), big_mat_all);
set(ax, 'YDir','reverse');
colormap(ax, flipud(gray(256)));
cH = prctile(big_mat_all(:), prcLim(2));
caxis(ax, [0, cH]);
% for sp = 1:numel(spk_ov_all)
%     r = spk_ov_all(sp).row; t_spk = spk_ov_all(sp).tval;
%     plot(ax, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
% end
xline(ax, 0, 'k--', 'LineWidth', 1);
hold(ax,'off');
xlim(ax, [t_w(1), t_w(end)]);
ylim(ax, [0.5, size(big_mat_all,1)+0.5]);
set(ax, 'YTick', rlp_all, 'YTickLabel', rlt_all, 'FontSize', 5);
xlabel(ax, 'Time (s)'); ylabel(ax, 'ROI');
title(ax, sprintf('Per peak heatmap ALL (%d)', nTotalROIs));
cb = colorbar(ax); ylabel(cb, 'dF/F');
set(ax, 'Units','centimeters','Position',[col1_x row2_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col1_x+sq_w+0.1, row2_y, 0.12, sq_h];

%% ---- (R3,C1): Per peak heatmap SIG ----
ax = axes(fig, 'Units','centimeters','Position',[col1_x row3_y sq_w sq_h]);
if nSigROIs > 0
    hold(ax,'on');
    imagesc(ax, t_w, 1:size(big_mat_sig,1), big_mat_sig);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    cH = prctile(big_mat_sig(:), prcLim(2));
    caxis(ax, [0, cH]);
    % for sp = 1:numel(spk_ov_sig)
    %     r = spk_ov_sig(sp).row; t_spk = spk_ov_sig(sp).tval;
    %     plot(ax, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
    % end
    xline(ax, 0, 'k--', 'LineWidth', 1);
    hold(ax,'off');
    xlim(ax, [t_w(1), t_w(end)]);
    ylim(ax, [0.5, size(big_mat_sig,1)+0.5]);
    set(ax, 'YTick', rlp_sig, 'YTickLabel', rlt_sig, 'FontSize', 5);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'ROI');
    title(ax, sprintf('Per peak heatmap SIG (%d/%d)', nSigROIs, nTotalROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col1_x row3_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col1_x+sq_w+0.1, row3_y, 0.12, sq_h];
else
    text(ax, 0.5,0.5, 'No sig ROIs', 'HorizontalAlignment','center', 'Units','normalized');
    title(ax, 'Per peak heatmap SIG (0)');
end

%% ---- (R4,C1): Sorted per peak hmap ALL (time) ----
ax = axes(fig, 'Units','centimeters','Position',[col1_x row4_y sq_w sq_h]);
hold(ax,'on');
imagesc(ax, t_w, 1:size(sorted_all_t,1), sorted_all_t);
set(ax, 'YDir','reverse');
colormap(ax, flipud(gray(256)));
cH = prctile(sorted_all_t(:), prcLim(2));
caxis(ax, [0, cH]);
% for sp = 1:numel(spk_sort_all_t)
%     r = spk_sort_all_t(sp).row; t_spk = spk_sort_all_t(sp).tval;
%     plot(ax, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
% end
xline(ax, 0, 'k--', 'LineWidth', 1);
hold(ax,'off');
xlim(ax, [t_w(1), t_w(end)]);
ylim(ax, [0.5, size(sorted_all_t,1)+0.5]);
xlabel(ax, 'Time (s)'); ylabel(ax, 'Cycle (sorted)');
title(ax, sprintf('Sorted per peak hmap ALL (%d)', nTotalROIs));
cb = colorbar(ax); ylabel(cb, 'dF/F');
set(ax, 'Units','centimeters','Position',[col1_x row4_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col1_x+sq_w+0.1, row4_y, 0.12, sq_h];

%% =====================================================================
%%  COLUMN 2 (TIME)
%% =====================================================================

%% ---- (R2,C2): Peak-trig avg hmap — ALL ----
peth_sorted_all_mat = peth_mat_all(sort_ord_all, :);
labels_sorted_all   = labels_all(sort_ord_all);

ax = axes(fig, 'Units','centimeters','Position',[col2_x row2_y sq_w sq_h]);
imagesc(ax, t_w, 1:nTotalROIs, peth_sorted_all_mat);
set(ax, 'YDir','reverse');
colormap(ax, flipud(gray(256)));
cH = prctile(peth_sorted_all_mat(:), prcLim(2));
caxis(ax, [0, cH]);
hold(ax,'on'); xline(ax, 0, 'k--', 'LineWidth', 1); hold(ax,'off');
if nTotalROIs <= 40
    set(ax, 'YTick', 1:nTotalROIs, 'YTickLabel', labels_sorted_all, 'FontSize', 5);
else
    ts = max(1, round(nTotalROIs/20));
    yt = 1:ts:nTotalROIs;
    set(ax, 'YTick', yt, 'YTickLabel', labels_sorted_all(yt), 'FontSize', 5);
end
xlabel(ax, 'Time (s)'); ylabel(ax, 'ROI');
title(ax, sprintf('Peak-trig avg hmap ALL (%d)', nTotalROIs));
cb = colorbar(ax); ylabel(cb, 'dF/F');
set(ax, 'Units','centimeters','Position',[col2_x row2_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col2_x+sq_w+0.1, row2_y, 0.12, sq_h];

%% ---- (R3,C2): Peak-trig avg hmap — SIG ----
ax = axes(fig, 'Units','centimeters','Position',[col2_x row3_y sq_w sq_h]);
if nSigROIs > 0
    peth_sorted_sig_mat = peth_mat_all(sort_ord_sig, :);
    labels_sorted_sig   = labels_all(sort_ord_sig);
    imagesc(ax, t_w, 1:nSigROIs, peth_sorted_sig_mat);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    cH = prctile(peth_sorted_sig_mat(:), prcLim(2));
    caxis(ax, [0, cH]);
    hold(ax,'on'); xline(ax, 0, 'k--', 'LineWidth', 1); hold(ax,'off');
    if nSigROIs <= 40
        set(ax, 'YTick', 1:nSigROIs, 'YTickLabel', labels_sorted_sig, 'FontSize', 5);
    else
        ts = max(1, round(nSigROIs/20));
        yt = 1:ts:nSigROIs;
        set(ax, 'YTick', yt, 'YTickLabel', labels_sorted_sig(yt), 'FontSize', 5);
    end
    xlabel(ax, 'Time (s)'); ylabel(ax, 'ROI');
    title(ax, sprintf('Peak-trig avg hmap SIG (%d/%d)', nSigROIs, nTotalROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col2_x row3_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col2_x+sq_w+0.1, row3_y, 0.12, sq_h];
else
    text(ax, 0.5,0.5, 'No sig ROIs', 'HorizontalAlignment','center', 'Units','normalized');
    title(ax, 'Peak-trig avg hmap SIG (0)');
end

%% ---- (R1,C2): Spike hist sig ----
dt_bin      = t_w(2) - t_w(1);
nBinsT      = numel(t_w);
peth_edges  = [t_w(:) - dt_bin/2; t_w(end) + dt_bin/2];
peth_counts = histcounts(all_spike_t, peth_edges);

ax = axes(fig, 'Units','centimeters','Position',[col2_x row1_y sq_w sq_h]);
hold(ax,'on');
bar(ax, t_w, peth_counts, 1, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
xline(ax, 0, 'k--', 'LineWidth', 1);
hold(ax,'off');
xlim(ax, [t_w(1), t_w(end)]);
xlabel(ax, 'Time (s)'); ylabel(ax, 'Spike count');
th = title(ax, sprintf('Spike hist sig (n=%d)', numel(all_spike_t)));
th.Units = 'normalized'; th.Position(2) = 1.08;
box(ax,'on'); pbaspect(ax,[1 1 1]);

%% ---- (R4,C2): Sorted per peak hmap SIG (time) ----
ax = axes(fig, 'Units','centimeters','Position',[col2_x row4_y sq_w sq_h]);
if nSigROIs > 0
    hold(ax,'on');
    imagesc(ax, t_w, 1:size(sorted_sig_t,1), sorted_sig_t);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    cH = prctile(sorted_sig_t(:), prcLim(2));
    caxis(ax, [0, cH]);
    % for sp = 1:numel(spk_sort_sig_t)
    %     r = spk_sort_sig_t(sp).row; t_spk = spk_sort_sig_t(sp).tval;
    %     plot(ax, [t_spk t_spk], [r-0.4 r+0.4], 'b', 'LineWidth', 1.5);
    % end
    xline(ax, 0, 'k--', 'LineWidth', 1);
    hold(ax,'off');
    xlim(ax, [t_w(1), t_w(end)]);
    ylim(ax, [0.5, size(sorted_sig_t,1)+0.5]);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Cycle (sorted)');
    title(ax, sprintf('Sorted hmap SIG (%d/%d)', nSigROIs, nTotalROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col2_x row4_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col2_x+sq_w+0.1, row4_y, 0.12, sq_h];
else
    text(ax, 0.5,0.5, 'No sig ROIs', 'HorizontalAlignment','center', 'Units','normalized');
    title(ax, 'Sorted hmap SIG (0)');
end

%% =====================================================================
%%  COLUMN 3 (PHASE)
%% =====================================================================

%% ---- (R1,C3): Breath overlay (phase) [0, 4pi] ----
ax = axes(fig, 'Units','centimeters','Position',[col3_x row1_y sq_w sq_h]);
hold(ax,'on');
patchColor = [0.85 0.85 0.85];
patchAlpha = 0.4;
patch(ax, [insp_phase_lo, insp_phase_hi, insp_phase_hi, insp_phase_lo], ...
    [-10 -10 10 10], patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
patch(ax, [insp_phase_lo+2*pi, insp_phase_hi+2*pi, insp_phase_hi+2*pi, insp_phase_lo+2*pi], ...
    [-10 -10 10 10], patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
plot(ax, phase_full, cycle_dup', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
xline(ax, 0, 'k--', 'LineWidth', 1);
xline(ax, 2*pi, 'k--', 'LineWidth', 1);
hold(ax,'off');
xlim(ax, [0, 4*pi]); ylim(ax, [0, 1]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Phase (rad)'); ylabel(ax, 'Breathing (norm)');
th = title(ax, sprintf('Breath overlay\n%d cycles from %d FOVs', nTotalCycles, nValidSess));
th.Units = 'normalized'; th.Position(2) = 1.16;
box(ax,'on'); pbaspect(ax,[1 1 1]);
yL = ax.YLim;
patches = findobj(ax, 'Type', 'patch');
for pp = 1:numel(patches), patches(pp).YData = [yL(1) yL(1) yL(2) yL(2)]; end

%% ---- (R3,C3): Per cycle heatmap SIG (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col3_x row3_y sq_w sq_h]);
if nSigROIs > 0
    ph_big_sig_dup = [ph_big_sig, ph_big_sig];
    cH = prctile(ph_big_sig_dup(:), prcLim(2));
    hold(ax,'on');
    imagesc(ax, phase_full, 1:size(ph_big_sig_dup,1), ph_big_sig_dup);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    caxis(ax, [0, cH]);
    % for sp_ = 1:numel(ph_spk_sig)
    %     r_ = ph_spk_sig(sp_).row; p_ = ph_spk_sig(sp_).phi;
    %     plot(ax, [p_ p_], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
    %     plot(ax, [p_+2*pi p_+2*pi], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
    % end
    xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
    hold(ax,'off');
    xlim(ax, [0, 4*pi]);
    ylim(ax, [0.5, size(ph_big_sig_dup,1)+0.5]);
    set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
        'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
    set(ax, 'YTick', ph_rlp_sig, 'YTickLabel', ph_rlt_sig, 'FontSize', 5);
    xlabel(ax, 'Phase'); ylabel(ax, 'ROI');
    title(ax, sprintf('Per cycle heatmap SIG (%d)', nSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col3_x row3_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col3_x+sq_w+0.1, row3_y, 0.12, sq_h];
else
    text(ax, 0.5,0.5, 'No sig ROIs', 'HorizontalAlignment','center', 'Units','normalized');
    title(ax, 'Per cycle heatmap SIG (0)');
end

%% ---- (R3,C4): Cycle-avg hmap SIG (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col4_x row3_y sq_w sq_h]);
if nSigROIs > 0
    avg_sig_dup = [avg_phase_sig, avg_phase_sig];
    imagesc(ax, phase_full, 1:nSigROIs, avg_sig_dup);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    cH = prctile(avg_sig_dup(:), prcLim(2));
    caxis(ax, [0, cH]);
    hold(ax,'on'); xline(ax, 2*pi, 'k--', 'LineWidth', 0.8); hold(ax,'off');
    set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
        'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
    if nSigROIs <= 40
        set(ax, 'YTick', 1:nSigROIs, 'YTickLabel', avg_labels_phase_sig, 'FontSize', 5);
    else
        ts = max(1, round(nSigROIs/20));
        yt = 1:ts:nSigROIs;
        set(ax, 'YTick', yt, 'YTickLabel', avg_labels_phase_sig(yt), 'FontSize', 5);
    end
    xlabel(ax, 'Phase'); ylabel(ax, 'ROI');
    title(ax, sprintf('Cycle-avg hmap SIG (%d)', nSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col4_x row3_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col4_x+sq_w+0.1, row3_y, 0.12, sq_h];
else
    text(ax, 0.5,0.5, 'No sig ROIs', 'HorizontalAlignment','center', 'Units','normalized');
    title(ax, 'Cycle-avg hmap SIG (0)');
end

%% ---- (R4,C3): Sorted per cycle hmap ALL (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col3_x row4_y sq_w sq_h]);
sorted_all_ph_dup = [sorted_all_ph, sorted_all_ph];
cH = prctile(sorted_all_ph_dup(:), prcLim(2));
hold(ax,'on');
imagesc(ax, phase_full, 1:size(sorted_all_ph_dup,1), sorted_all_ph_dup);
set(ax, 'YDir','reverse');
colormap(ax, flipud(gray(256)));
caxis(ax, [0, cH]);
% for sp_ = 1:numel(spk_sort_all_ph)
%     r_ = spk_sort_all_ph(sp_).row; p_ = spk_sort_all_ph(sp_).phi;
%     plot(ax, [p_ p_], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
%     plot(ax, [p_+2*pi p_+2*pi], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
% end
xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
hold(ax,'off');
xlim(ax, [0, 4*pi]);
ylim(ax, [0.5, size(sorted_all_ph_dup,1)+0.5]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Phase'); ylabel(ax, 'Cycle (sorted)');
title(ax, sprintf('Sorted hmap ALL (%d)', nTotalROIs));
cb = colorbar(ax); ylabel(cb, 'dF/F');
set(ax, 'Units','centimeters','Position',[col3_x row4_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col3_x+sq_w+0.1, row4_y, 0.12, sq_h];

%% =====================================================================
%%  COLUMN 4 (PHASE)
%% =====================================================================

%% ---- (R2,C3): Per cycle heatmap ALL (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col3_x row2_y sq_w sq_h]);
ph_big_all_dup = [ph_big_all, ph_big_all];
cH = prctile(ph_big_all_dup(:), prcLim(2));
hold(ax,'on');
imagesc(ax, phase_full, 1:size(ph_big_all_dup,1), ph_big_all_dup);
set(ax, 'YDir','reverse');
colormap(ax, flipud(gray(256)));
caxis(ax, [0, cH]);
% for sp_ = 1:numel(ph_spk_all)
%     r_ = ph_spk_all(sp_).row; p_ = ph_spk_all(sp_).phi;
%     plot(ax, [p_ p_], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
%     plot(ax, [p_+2*pi p_+2*pi], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
% end
xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
hold(ax,'off');
xlim(ax, [0, 4*pi]);
ylim(ax, [0.5, size(ph_big_all_dup,1)+0.5]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
set(ax, 'YTick', ph_rlp_all, 'YTickLabel', ph_rlt_all, 'FontSize', 5);
xlabel(ax, 'Phase'); ylabel(ax, 'ROI');
title(ax, sprintf('Per cycle heatmap ALL (%d)', nTotalROIs));
cb = colorbar(ax); ylabel(cb, 'dF/F');
set(ax, 'Units','centimeters','Position',[col3_x row2_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col3_x+sq_w+0.1, row2_y, 0.12, sq_h];

%% ---- (R2,C4): Cycle-avg hmap ALL (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col4_x row2_y sq_w sq_h]);
avg_all_dup = [avg_phase_all, avg_phase_all];
imagesc(ax, phase_full, 1:nTotalROIs, avg_all_dup);
set(ax, 'YDir','reverse');
colormap(ax, flipud(gray(256)));
cH = prctile(avg_all_dup(:), prcLim(2));
caxis(ax, [0, cH]);
hold(ax,'on'); xline(ax, 2*pi, 'k--', 'LineWidth', 0.8); hold(ax,'off');
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
if nTotalROIs <= 40
    set(ax, 'YTick', 1:nTotalROIs, 'YTickLabel', avg_labels_phase_all, 'FontSize', 5);
else
    ts = max(1, round(nTotalROIs/20));
    yt = 1:ts:nTotalROIs;
    set(ax, 'YTick', yt, 'YTickLabel', avg_labels_phase_all(yt), 'FontSize', 5);
end
xlabel(ax, 'Phase'); ylabel(ax, 'ROI');
title(ax, sprintf('Cycle-avg hmap ALL (%d)', nTotalROIs));
cb = colorbar(ax); ylabel(cb, 'dF/F');
set(ax, 'Units','centimeters','Position',[col4_x row2_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col4_x+sq_w+0.1, row2_y, 0.12, sq_h];

%% ---- (R1,C4): Spike probability sig + cosine fit ----
ax = axes(fig, 'Units','centimeters','Position',[col4_x row1_y sq_w sq_h]);
hold(ax,'on');
nBinsHist = 12;
edges_hist = linspace(0, 2*pi, nBinsHist + 1);
ctrs_hist  = (edges_hist(1:end-1) + edges_hist(2:end)) / 2;
cnt_hist   = histcounts(all_spike_ph_pool, edges_hist);
prob_hist  = cnt_hist / sum(cnt_hist);   % sums to 1 over one cycle

% Cosine fit on [0, 2pi]
mu_fit = angle(mean(exp(1i * all_spike_ph_pool)));
a_fit  = mean(prob_hist);
b_fit  = 2 * mean(prob_hist .* cos(ctrs_hist - mu_fit));

% Duplicate to [0, 4pi]
ctrs_dup  = [ctrs_hist, ctrs_hist + 2*pi];
prob_dup  = [prob_hist, prob_hist];

% Observed phase probability (gray line — uniform reference)
cnt_obs_ph  = histcounts(all_frame_phase, edges_hist);
prob_obs_ph = cnt_obs_ph / sum(cnt_obs_ph);
prob_obs_dup = [prob_obs_ph, prob_obs_ph];
plot(ax, ctrs_dup, prob_obs_dup, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

bar(ax, ctrs_dup, prob_dup, 1, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
x_fit = linspace(0, 4*pi, 500);
y_fit = a_fit + b_fit * cos(x_fit - mu_fit);
plot(ax, x_fit, y_fit, 'r-', 'LineWidth', 2);
xline(ax, 2*pi, 'k--', 'LineWidth', 1);
xlim(ax, [0, 4*pi]); ylim(ax, [0 0.5]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Phase (rad)'); ylabel(ax, 'P(spike)');
th = title(ax, sprintf('Spike probability sig (n=%d)', numel(all_spike_ph_pool)));
th.Units = 'normalized'; th.Position(2) = 1.08;
box(ax,'on'); pbaspect(ax,[1 1 1]);
hold(ax,'off');

%% ---- (R4,C4): Sorted per cycle hmap SIG (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col4_x row4_y sq_w sq_h]);
if nSigROIs > 0
    sorted_sig_ph_dup = [sorted_sig_ph, sorted_sig_ph];
    cH = prctile(sorted_sig_ph_dup(:), prcLim(2));
    hold(ax,'on');
    imagesc(ax, phase_full, 1:size(sorted_sig_ph_dup,1), sorted_sig_ph_dup);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    caxis(ax, [0, cH]);
    % for sp_ = 1:numel(spk_sort_sig_ph)
    %     r_ = spk_sort_sig_ph(sp_).row; p_ = spk_sort_sig_ph(sp_).phi;
    %     plot(ax, [p_ p_], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
    %     plot(ax, [p_+2*pi p_+2*pi], [r_-0.4 r_+0.4], 'b', 'LineWidth', 1.5);
    % end
    xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
    hold(ax,'off');
    xlim(ax, [0, 4*pi]);
    ylim(ax, [0.5, size(sorted_sig_ph_dup,1)+0.5]);
    set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
        'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
    xlabel(ax, 'Phase'); ylabel(ax, 'Cycle (sorted)');
    title(ax, sprintf('Sorted hmap SIG (%d)', nSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col4_x row4_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col4_x+sq_w+0.1, row4_y, 0.12, sq_h];
else
    text(ax, 0.5,0.5, 'No sig ROIs', 'HorizontalAlignment','center', 'Units','normalized');
    title(ax, 'Sorted hmap SIG (0)');
end

%% =====================================================================
%%  COLUMN 5 (PHASE: polars, Poincare, KDE)
%% =====================================================================

%% ---- (R1,C5): Polar coherence ----
polar_norm = [col5_x/fig_w, row1_y/fig_h, sq_w/fig_w, sq_h/fig_h];
ax_pol1 = polaraxes(fig, 'Position', polar_norm);
hold(ax_pol1,'on');
th_circ = linspace(0, 2*pi, 360);
polarplot(ax_pol1, th_circ, repmat(confC_mean, 1, 360), 'k--', 'LineWidth', 1);
polarplot(ax_pol1, [phase_1cyc, 2*pi], [breath_norm_polar, breath_norm_polar(1)], 'k-', 'LineWidth', 0.8);
for rr = 1:nTotalROIs
    if all_is_sig_vec(rr)
        c_r = sessColor(all_roi_sess(rr), :);
    else
        c_r = [0 0 0];
    end
    polarplot(ax_pol1, all_th_k(rr), all_r_k(rr), 'o', ...
        'MarkerFaceColor', c_r, 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5);
end
ax_pol1.RLim = [0, 1];
ax_pol1.ThetaZeroLocation = 'top';
ax_pol1.ThetaDir = 'clockwise';
ax_pol1.RAxisLocation = 180;
ax_pol1.FontSize = 6;
hold(ax_pol1,'off');
nSig = sum(all_is_sig_vec);
title(ax_pol1, sprintf('Coherence (%d/%d)', nSig, nTotalROIs), 'Interpreter','none');

%% ---- (R2,C5): Polar %active ----
polar_norm2 = [col5_x/fig_w, row2_y/fig_h, sq_w/fig_w, sq_h/fig_h];
ax_pol2 = polaraxes(fig, 'Position', polar_norm2);
hold(ax_pol2,'on');
polarplot(ax_pol2, [phase_1cyc, 2*pi], [breath_norm_polar, breath_norm_polar(1)], 'k-', 'LineWidth', 0.8);
for rr = 1:nTotalROIs
    if ~all_is_sig_vec(rr), continue; end
    c_r = sessColor(all_roi_sess(rr), :);
    r_pct = all_pct_active(rr) / 100;
    polarplot(ax_pol2, all_th_k(rr), r_pct, 'o', ...
        'MarkerFaceColor', c_r, 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5);
end
ax_pol2.RLim = [0, 1];
ax_pol2.ThetaZeroLocation = 'top';
ax_pol2.ThetaDir = 'clockwise';
ax_pol2.RAxisLocation = 180;
ax_pol2.FontSize = 6;
hold(ax_pol2,'off');
title(ax_pol2, sprintf('%%Active (%d sig)', nSig), 'Interpreter','none');

% Session colorbar for polars
ax_cb_p = axes(fig, 'Units','centimeters', ...
    'Position', [col5_x+sq_w+0.15, row2_y, 0.3, sq_h], 'Visible','off');
colormap(ax_cb_p, cmap_sess);
caxis(ax_cb_p, [1 nValidSess]);
cb_p = colorbar(ax_cb_p, 'Location','eastoutside');
cb_p.Label.String = 'Session';
cb_p.FontSize = 6;
if nValidSess > 1
    cb_p.Ticks = unique(round([ceil(nValidSess/2), nValidSess]));
else
    cb_p.Ticks = 1;
end

%% ---- (R3,C5): Poincare ----
ax = axes(fig, 'Units','centimeters','Position',[col5_x row3_y sq_w sq_h]);
hold(ax,'on');
for ks = 1:nSess
    if ~validSess(ks), continue; end
    ibi = ibi_all{ks};
    if isempty(ibi) || numel(ibi) < 2, continue; end
    c_s = sessColor(ks, :);
    ran1 = (rand(size(ibi)) - 0.5) .* 2 ./ fs ./ 2;
    ran2 = (rand(size(ibi)) - 0.5) .* 2 ./ fs ./ 2;
    scatter(ax, ibi(1:end-1)+ran1(1:end-1), ibi(2:end)+ran2(2:end), 6, c_s, ...
        'MarkerEdgeAlpha', 0.3);
end
plot(ax, ibi_xlim, ibi_xlim, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
xlim(ax, ibi_xlim); ylim(ax, ibi_xlim);
axis(ax, 'square');
xlabel(ax, 'IBI_n (s)'); ylabel(ax, 'IBI_{n+1} (s)');
th = title(ax, 'Poincare', 'Interpreter','none');
th.Units = 'normalized'; th.Position(2) = 1.04;
box(ax,'on');
colormap(ax, cmap_sess); caxis(ax, [1 nValidSess]);
cb = colorbar(ax); cb.Label.String = 'Session'; cb.FontSize = 6;
if nValidSess > 1
    cb.Ticks = unique(round([ceil(nValidSess/2), nValidSess]));
else
    cb.Ticks = 1;
end
set(ax, 'Units','centimeters','Position',[col5_x row3_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col5_x+sq_w+0.1, row3_y, 0.12, sq_h];

%% ---- (R4,C5): KDE ----
ax = axes(fig, 'Units','centimeters','Position',[col5_x row4_y sq_w sq_h]);
hold(ax,'on');
ibi_pooled = vertcat(ibi_all{validSess});
if ~isempty(ibi_pooled)
    xGrid = linspace(0, max(ibi_pooled)*1.5, 2000);
else
    xGrid = linspace(0, 2, 2000);
end
for ks = 1:nSess
    if ~validSess(ks), continue; end
    ibi = ibi_all{ks};
    ibi = ibi(isfinite(ibi) & ibi > 0);
    if numel(ibi) < 2, continue; end
    c_s = sessColor(ks, :);
    [fKDE, xi] = ksdensity(ibi, xGrid);
    plot(ax, xi, fKDE, '-', 'Color', c_s, 'LineWidth', 1.5);
end
xlim(ax, ibi_xlim);
xlabel(ax, 'IBI (s)'); ylabel(ax, 'Density');
th = title(ax, 'IBI KDE', 'Interpreter','none');
th.Units = 'normalized'; th.Position(2) = 1.04;
box(ax,'on'); axis(ax,'square');
colormap(ax, cmap_sess); caxis(ax, [1 nValidSess]);
cb = colorbar(ax); cb.Label.String = 'Session'; cb.FontSize = 6;
if nValidSess > 1
    cb.Ticks = unique(round([ceil(nValidSess/2), nValidSess]));
else
    cb.Ticks = 1;
end
set(ax, 'Units','centimeters','Position',[col5_x row4_y sq_w sq_h]);
cb.Units = 'centimeters'; cb.Position = [col5_x+sq_w+0.1, row4_y, 0.12, sq_h];

hold(ax,'off');

%% =====================================================================
%%  COLUMN 6: NON-SIG ROIs (time R1-R2, phase R3-R4)
%% =====================================================================

%% ---- (R1,C6): Per-cycle heatmap non-sig (time) ----
ax = axes(fig, 'Units','centimeters','Position',[col6_x row1_y sq_w sq_h]);
if nNonSigROIs > 0
    hold(ax,'on');
    imagesc(ax, t_w, 1:size(big_mat_nsig,1), big_mat_nsig);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    cH = prctile(big_mat_nsig(:), prcLim(2));
    caxis(ax, [0, cH]);
    xline(ax, 0, 'k--', 'LineWidth', 1);
    hold(ax,'off');
    xlim(ax, [t_w(1), t_w(end)]);
    ylim(ax, [0.5, size(big_mat_nsig,1)+0.5]);
    set(ax, 'YTick', rlp_nsig, 'YTickLabel', rlt_nsig, 'FontSize', 5);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'ROI');
    title(ax, sprintf('Per peak heatmap n.s. (%d)', nNonSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col6_x row1_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col6_x+sq_w+0.1, row1_y, 0.12, sq_h];
else
    title(ax, 'Per peak heatmap n.s. (0)');
end

%% ---- (R2,C6): Sorted per peak hmap non-sig (time) ----
ax = axes(fig, 'Units','centimeters','Position',[col6_x row2_y sq_w sq_h]);
if nNonSigROIs > 0
    hold(ax,'on');
    imagesc(ax, t_w, 1:size(sorted_nsig_t,1), sorted_nsig_t);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    cH = prctile(sorted_nsig_t(:), prcLim(2));
    caxis(ax, [0, cH]);
    xline(ax, 0, 'k--', 'LineWidth', 1);
    hold(ax,'off');
    xlim(ax, [t_w(1), t_w(end)]);
    ylim(ax, [0.5, size(sorted_nsig_t,1)+0.5]);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Cycle');
    title(ax, sprintf('Sorted hmap n.s. (%d)', nNonSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col6_x row2_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col6_x+sq_w+0.1, row2_y, 0.12, sq_h];
else
    title(ax, 'Sorted hmap n.s. (0)');
end

%% ---- (R3,C6): Per-cycle heatmap non-sig (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col6_x row3_y sq_w sq_h]);
if nNonSigROIs > 0
    ph_big_nsig_dup = [ph_big_nsig, ph_big_nsig];
    cH = prctile(ph_big_nsig_dup(:), prcLim(2));
    hold(ax,'on');
    imagesc(ax, phase_full, 1:size(ph_big_nsig_dup,1), ph_big_nsig_dup);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    caxis(ax, [0, cH]);
    xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
    hold(ax,'off');
    xlim(ax, [0, 4*pi]);
    ylim(ax, [0.5, size(ph_big_nsig_dup,1)+0.5]);
    set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], 'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
    set(ax, 'YTick', ph_rlp_nsig, 'YTickLabel', ph_rlt_nsig, 'FontSize', 5);
    xlabel(ax, 'Phase'); ylabel(ax, 'ROI');
    title(ax, sprintf('Per cycle heatmap n.s. (%d)', nNonSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col6_x row3_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col6_x+sq_w+0.1, row3_y, 0.12, sq_h];
else
    title(ax, 'Per cycle heatmap n.s. (0)');
end

%% ---- (R4,C6): Sorted per cycle hmap non-sig (phase) ----
ax = axes(fig, 'Units','centimeters','Position',[col6_x row4_y sq_w sq_h]);
if nNonSigROIs > 0
    sorted_nsig_ph_dup = [sorted_nsig_ph, sorted_nsig_ph];
    cH = prctile(sorted_nsig_ph_dup(:), prcLim(2));
    hold(ax,'on');
    imagesc(ax, phase_full, 1:size(sorted_nsig_ph_dup,1), sorted_nsig_ph_dup);
    set(ax, 'YDir','reverse');
    colormap(ax, flipud(gray(256)));
    caxis(ax, [0, cH]);
    xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
    hold(ax,'off');
    xlim(ax, [0, 4*pi]);
    ylim(ax, [0.5, size(sorted_nsig_ph_dup,1)+0.5]);
    set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], 'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
    xlabel(ax, 'Phase'); ylabel(ax, 'Cycle');
    title(ax, sprintf('Sorted hmap n.s. (%d)', nNonSigROIs));
    cb = colorbar(ax); ylabel(cb, 'dF/F');
    set(ax, 'Units','centimeters','Position',[col6_x row4_y sq_w sq_h]);
    cb.Units = 'centimeters'; cb.Position = [col6_x+sq_w+0.1, row4_y, 0.12, sq_h];
else
    title(ax, 'Sorted hmap n.s. (0)');
end

%% =====================================================================
%%  COLUMN 7: NON-SIG spike hist + phase hist
%% =====================================================================

%% ---- (R3,C7): Spike hist n.s. (time) ----
ax = axes(fig, 'Units','centimeters','Position',[col7_x row3_y sq_w sq_h]);
hold(ax,'on');
if ~isempty(all_spike_t_nsig)
    peth_counts_nsig = histcounts(all_spike_t_nsig, peth_edges);
    bar(ax, t_w, peth_counts_nsig, 1, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
end
xline(ax, 0, 'k--', 'LineWidth', 1);
hold(ax,'off');
xlim(ax, [t_w(1), t_w(end)]);
xlabel(ax, 'Time (s)'); ylabel(ax, 'Spike count');
th = title(ax, sprintf('Spike hist n.s. (n=%d)', numel(all_spike_t_nsig)));
th.Units = 'normalized'; th.Position(2) = 1.08;
box(ax,'on'); pbaspect(ax,[1 1 1]);

%% ---- (R4,C7): Phase hist n.s. ----
ax = axes(fig, 'Units','centimeters','Position',[col7_x row4_y sq_w sq_h]);
hold(ax,'on');
if ~isempty(all_spike_ph_nsig)
    cnt_nsig_ph = histcounts(all_spike_ph_nsig, edges_hist);
    prob_nsig_ph = cnt_nsig_ph / sum(cnt_nsig_ph);

    % Cosine fit
    mu_fit_ns = angle(mean(exp(1i * all_spike_ph_nsig)));
    a_fit_ns  = mean(prob_nsig_ph);
    b_fit_ns  = 2 * mean(prob_nsig_ph .* cos(ctrs_hist - mu_fit_ns));

    ctrs_dup_ns = [ctrs_hist, ctrs_hist + 2*pi];
    prob_dup_ns = [prob_nsig_ph, prob_nsig_ph];

    % Gray observed reference
    plot(ax, ctrs_dup, prob_obs_dup, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

    bar(ax, ctrs_dup_ns, prob_dup_ns, 1, 'FaceColor', [0.3 0.3 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    x_fit_ns = linspace(0, 4*pi, 500);
    y_fit_ns = a_fit_ns + b_fit_ns * cos(x_fit_ns - mu_fit_ns);
    plot(ax, x_fit_ns, y_fit_ns, 'r-', 'LineWidth', 2);
end
xline(ax, 2*pi, 'k--', 'LineWidth', 1);
xlim(ax, [0, 4*pi]); ylim(ax, [0 0.5]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Phase (rad)'); ylabel(ax, 'P(spike)');
th = title(ax, sprintf('Spike probability n.s. (n=%d)', numel(all_spike_ph_nsig)));
th.Units = 'normalized'; th.Position(2) = 1.08;
box(ax,'on'); pbaspect(ax,[1 1 1]);
hold(ax,'off');

%% ---- (R1,C7): Grand avg dF/F trace — SIG ROIs ----
ax = axes(fig, 'Units','centimeters','Position',[col7_x row1_y sq_w sq_h]);
hold(ax,'on');
if nSigROIs > 0
    grand_avg_sig = mean(avg_phase_sig, 1, 'omitnan');
    grand_avg_sig_dup = [grand_avg_sig, grand_avg_sig];
    plot(ax, phase_full, grand_avg_sig_dup, 'k-', 'LineWidth', 1.5);
end
xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
hold(ax,'off');
xlim(ax, [0, 4*pi]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], 'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Phase'); ylabel(ax, 'dF/F');
th = title(ax, sprintf('Avg dF/F SIG (%d)', nSigROIs));
th.Units = 'normalized'; th.Position(2) = 1.04;
box(ax,'on'); pbaspect(ax,[1 1 1]);

%% ---- (R2,C7): Grand avg dF/F trace — non-sig ROIs ----
ax = axes(fig, 'Units','centimeters','Position',[col7_x row2_y sq_w sq_h]);
hold(ax,'on');
if nNonSigROIs > 0
    grand_avg_nsig = mean(avg_phase_nsig, 1, 'omitnan');
    grand_avg_nsig_dup = [grand_avg_nsig, grand_avg_nsig];
    plot(ax, phase_full, grand_avg_nsig_dup, 'k-', 'LineWidth', 1.5);
end
xline(ax, 2*pi, 'k--', 'LineWidth', 0.8);
hold(ax,'off');
xlim(ax, [0, 4*pi]);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], 'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Phase'); ylabel(ax, 'dF/F');
th = title(ax, sprintf('Avg dF/F n.s. (%d)', nNonSigROIs));
th.Units = 'normalized'; th.Position(2) = 1.04;
box(ax,'on'); pbaspect(ax,[1 1 1]);

%% ---- (R1,C7 overlay): Per-ROI avg dF/F traces — SIG ----
% Add individual ROI traces (gray) behind the grand avg (black) on same panel
if nSigROIs > 0
    % Get the axes handle for R1,C7 (already created above)
    % Re-plot with per-ROI traces
    ax_r1c7 = findobj(fig, 'Type', 'axes', 'Position', [col7_x row1_y sq_w sq_h]);
    if ~isempty(ax_r1c7)
        hold(ax_r1c7, 'on');
        avg_sig_dup_lines = [avg_phase_sig, avg_phase_sig];
        for rr_ = 1:nSigROIs
            plot(ax_r1c7, phase_full, avg_sig_dup_lines(rr_,:), ...
                'Color', [0.7 0.7 0.7 0.4], 'LineWidth', 0.5);
        end
        % Re-plot grand avg on top
        grand_avg_sig_dup = [mean(avg_phase_sig, 1, 'omitnan'), mean(avg_phase_sig, 1, 'omitnan')];
        plot(ax_r1c7, phase_full, grand_avg_sig_dup, 'k-', 'LineWidth', 1.5);
        hold(ax_r1c7, 'off');
    end
end

%% ---- (R2,C7 overlay): Per-ROI avg dF/F traces — non-sig ----
if nNonSigROIs > 0
    ax_r2c7 = findobj(fig, 'Type', 'axes', 'Position', [col7_x row2_y sq_w sq_h]);
    if ~isempty(ax_r2c7)
        hold(ax_r2c7, 'on');
        avg_nsig_dup_lines = [avg_phase_nsig, avg_phase_nsig];
        for rr_ = 1:nNonSigROIs
            plot(ax_r2c7, phase_full, avg_nsig_dup_lines(rr_,:), ...
                'Color', [0.7 0.7 0.7 0.4], 'LineWidth', 0.5);
        end
        grand_avg_nsig_dup = [mean(avg_phase_nsig, 1, 'omitnan'), mean(avg_phase_nsig, 1, 'omitnan')];
        plot(ax_r2c7, phase_full, grand_avg_nsig_dup, 'k-', 'LineWidth', 1.5);
        hold(ax_r2c7, 'off');
    end
end

% (stackDFF columns removed)

%% =====================================================================
%%  COLUMN 8: Spectra (breathing + ROI)
%% =====================================================================

%% ---- (R1-R2,C10): Breathing PSD — mean ± std across sessions ----
tall_h10 = 2*sq_h + vgap;
ax = axes(fig, 'Units','centimeters','Position',[col8_x row2_y sq_w tall_h10]);
hold(ax,'on');

valid_psd_idx = find(~cellfun(@isempty, all_breath_psd));
if numel(valid_psd_idx) >= 1
    f_common = all_breath_freq{valid_psd_idx(1)};
    psd_mat = nan(numel(valid_psd_idx), numel(f_common));
    for ii = 1:numel(valid_psd_idx)
        fi = all_breath_freq{valid_psd_idx(ii)};
        si = all_breath_psd{valid_psd_idx(ii)};
        if numel(fi) == numel(f_common)
            psd_mat(ii,:) = si(:)';
        else
            psd_mat(ii,:) = interp1(fi, si, f_common, 'linear', NaN);
        end
    end
    m_lin = mean(psd_mat, 1, 'omitnan');
    s_lin = std(psd_mat, 0, 1, 'omitnan');
%     % Individual sessions in gray (commented out)
%     for ii = 1:size(psd_mat, 1)
%         plot(ax, f_common, 10*log10(psd_mat(ii,:)), '-', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5);
%     end
    fill(ax, [f_common(:); flipud(f_common(:))], ...
        10*log10([m_lin(:)+s_lin(:); flipud(m_lin(:)-s_lin(:))]), ...
        'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(ax, f_common, 10*log10(m_lin), 'k-', 'LineWidth', 1.5);
end
hold(ax,'off');
xlim(ax, [fmin, min(fmax, 15)]);
xlabel(ax, 'Freq (Hz)');
ylabel(ax, '$|\tilde{Z}(f)|^2 \; (\mathrm{Z}^2/\mathrm{Hz})$', 'Interpreter','latex');
th = title(ax, sprintf('Breath PSD (%d FOVs)', numel(valid_psd_idx)));
th.Units = 'normalized'; th.Position(2) = 1.02;
box(ax,'on');

%% ---- (R3-R4,C10): ROI PSD — all ROIs gray ----
ax = axes(fig, 'Units','centimeters','Position',[col8_x row4_y sq_w tall_h10]);
hold(ax,'on');
for ii = 1:numel(all_roi_psd)
    plot(ax, all_roi_psd_freq{ii}, all_roi_psd{ii}, '-', ...
        'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.5);
end
if numel(all_roi_psd) >= 1
    f_common_r = all_roi_psd_freq{1};
    roi_psd_mat = nan(numel(all_roi_psd), numel(f_common_r));
    for ii = 1:numel(all_roi_psd)
        fi = all_roi_psd_freq{ii};
        si = all_roi_psd{ii};
        if numel(fi) == numel(f_common_r)
            roi_psd_mat(ii,:) = si(:)';
        else
            roi_psd_mat(ii,:) = interp1(fi, si, f_common_r, 'linear', NaN);
        end
    end
    plot(ax, f_common_r, mean(roi_psd_mat, 1, 'omitnan'), 'k-', 'LineWidth', 1.5);
end
hold(ax,'off');
xlim(ax, [fmin, min(fmax, 15)]);
xlabel(ax, 'Freq (Hz)');
ylabel(ax, '$|\tilde{F}''(f)|^2 \; ([\mathrm{AU/s}]^2/\mathrm{Hz})$', 'Interpreter','latex');
th = title(ax, sprintf('F PSD (%d/%d ROIs)', nSigROIs, nTotalROIs));
th.Units = 'normalized'; th.Position(2) = 1.02;
box(ax,'on');

%% ---- Save PDF ----
if doSave
    outDir = fullfile(inputPath, 'breath_combined_summary');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    pdfPath = fullfile(outDir, 'breath_combined_summary.pdf');
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('Saved: %s\n', pdfPath);
end

%% ---- Save fov_map.mat ----
outDir_map = fullfile(inputPath, 'breath_combined_summary');
if ~exist(outDir_map, 'dir'), mkdir(outDir_map); end
validIdx_fm = find(validSess);
fov_map = cell(numel(validIdx_fm), 2);
for ff = 1:numel(validIdx_fm)
    kk_ = validIdx_fm(ff);
    fov_map{ff,1} = sessNames{kk_};
    sig_rois = [];
    for rr = 1:numel(all_roi_data)
        if all_roi_data(rr).sessIdx == kk_ && all_roi_data(rr).is_sig
            sig_rois = [sig_rois, all_roi_data(rr).roiId]; %#ok<AGROW>
        end
    end
    fov_map{ff,2} = sig_rois;
end
save(fullfile(outDir_map, 'fov_map.mat'), 'fov_map');
fprintf('Saved fov_map.mat (%d FOVs)\n', size(fov_map, 1));

%% =====================================================================
%%  FIGURE 2: Individual session breathing spectra (one panel per session)
%% =====================================================================
valid_psd_idx2 = find(~cellfun(@isempty, all_breath_psd));
nPsdSess = numel(valid_psd_idx2);
if nPsdSess > 0
    fig2_sq = 2.5;  % panel size (cm)
    fig2_sq_x = fig2_sq / 1.5 / 1.5;  % narrower panels
    fig2_gap = 0.8;
    fig2_vgap = 1.2;
    fig2_ml = 1.5; fig2_mr = 0.5; fig2_mb = 1.2; fig2_mt = 0.8;
    fig2_w = fig2_ml + nPsdSess * fig2_sq_x + (nPsdSess-1) * fig2_gap + fig2_mr;
    fig2_h = fig2_mb + 2*fig2_sq + fig2_vgap + fig2_mt;

    fig2 = figure('Color','w', 'Name','Breath PSD per session', ...
        'Units','centimeters', 'Position', [0.5 0.5 fig2_w fig2_h], ...
        'PaperUnits','centimeters', 'PaperSize', [fig2_w fig2_h], ...
        'PaperPosition', [0 0 fig2_w fig2_h]);
    set(fig2, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 7);

    row1_y2 = fig2_mb + fig2_sq + fig2_vgap;   % top row (dB)
    row2_y2 = fig2_mb;                          % bottom row (cumulative integral)

    % Compute shared y-limits across all sessions (dB row)
    all_dB_vals = [];
    for ii = 1:nPsdSess
        si = all_breath_psd{valid_psd_idx2(ii)};
        all_dB_vals = [all_dB_vals; 10*log10(si(:))]; %#ok<AGROW>
    end
    all_dB_vals = all_dB_vals(isfinite(all_dB_vals));
    yLo = floor(min(all_dB_vals)/5)*5;
    yHi = ceil(max(all_dB_vals)/5)*5;

    % Compute shared y-limits for cumulative integral (linear)
    all_cum_max = 0;
    for ii = 1:nPsdSess
        kk_ = valid_psd_idx2(ii);
        fi = all_breath_freq{kk_};
        si = all_breath_psd{kk_};
        cum_i = cumtrapz(fi(:), si(:));
        all_cum_max = max(all_cum_max, cum_i(end));
    end

    for ii = 1:nPsdSess
        kk_ = valid_psd_idx2(ii);
        fi = all_breath_freq{kk_};
        si = all_breath_psd{kk_};
        x0 = fig2_ml + (ii-1) * (fig2_sq_x + fig2_gap);

        % Row 1: PSD in dB
        ax2 = axes(fig2, 'Units','centimeters', 'Position', [x0 row1_y2 fig2_sq_x fig2_sq]);
        plot(ax2, fi, 10*log10(si), 'k-', 'LineWidth', 1);
        xlim(ax2, [fmin, min(fmax, 15)]);
        ylim(ax2, [yLo yHi]);
        box(ax2, 'on');
        title(ax2, sprintf('FOV#%d', ii), 'FontSize', 6);
        if ii == 1
            xlabel(ax2, 'Hz');
            ylabel(ax2, '$|\tilde{Z}(f)|^2 \; (\mathrm{Z}^2/\mathrm{Hz})$', 'Interpreter','latex');
        else
            set(ax2, 'YTickLabel', []);
        end

        % Row 2: cumulative integral (linear)
        cum_psd = cumtrapz(fi(:), si(:));
        ax3 = axes(fig2, 'Units','centimeters', 'Position', [x0 row2_y2 fig2_sq_x fig2_sq]);
        plot(ax3, fi, cum_psd, 'k-', 'LineWidth', 1);
        xlim(ax3, [fmin, min(fmax, 15)]);
        ylim(ax3, [0, all_cum_max * 1.05]);
        box(ax3, 'on');
        title(ax3, sprintf('\\int=%.2g', cum_psd(end)), 'FontSize', 6);
        if ii == 1
            xlabel(ax3, 'Hz'); ylabel(ax3, '\intPSD');
        else
            set(ax3, 'YTickLabel', []);
        end
    end

    if doSave
        outDir2 = fullfile(inputPath, 'breath_combined_summary');
        if ~exist(outDir2, 'dir'), mkdir(outDir2); end
        pdfPath2 = fullfile(outDir2, 'breath_psd_per_session.pdf');
        exportgraphics(fig2, pdfPath2, 'ContentType', 'vector', 'BackgroundColor', 'none');
        fprintf('Saved: %s\n', pdfPath2);
    end
end

fprintf('\nDone. %d sessions, %d snippets, %d cycles, %d ROIs (%d sig, %d n.s.).\n', ...
    nValidSess, nTotalSnip, nTotalCycles, nTotalROIs, nSigROIs, nNonSigROIs);

%% =====================================================================
%%  LOCAL FUNCTIONS
%% =====================================================================

%% ===== local_build_stack: per-cycle stacked matrix (time domain) =====
function [big_mat, spike_overlay, row_labels_pos, row_labels_txt] = ...
        local_build_stack(all_roi, roi_order, nW, field_snips, field_spk)
% Stack per-cycle dF/F snippets for ROIs in roi_order, with 2-row NaN gap.
% field_snips: name of the snippet field (e.g. 'dff_snips')
% field_spk:   name of the spike-times field (e.g. 'spike_times_per_peak')

    big_mat        = [];
    spike_overlay  = struct('row',{}, 'tval',{});
    row_labels_pos = [];
    row_labels_txt = {};

    for ii = 1:numel(roi_order)
        rd       = all_roi(roi_order(ii));
        dff_cyc  = rd.(field_snips);           % [nPeaks x nW]
        spk_data = rd.(field_spk);
        nPeaks   = size(dff_cyc, 1);

        % Filter to spiking cycles only
        spiking_mask = false(nPeaks, 1);
        for e = 1:nPeaks
            spiking_mask(e) = ~isempty(spk_data{e});
        end
        spiking_idx = find(spiking_mask);

        if isempty(spiking_idx)
            continue;
        end

        cyc_data  = dff_cyc(spiking_idx, :);
        row_start = size(big_mat, 1) + 1;
        big_mat   = [big_mat; cyc_data; nan(2, nW)]; %#ok<AGROW>

        for cc = 1:numel(spiking_idx)
            e = spiking_idx(cc);
            spk_t = spk_data{e};
            row_in = row_start + cc - 1;
            for s = 1:numel(spk_t)
                spike_overlay(end+1).row  = row_in;     %#ok<AGROW>
                spike_overlay(end).tval   = spk_t(s);
            end
        end

        row_labels_pos(end+1) = row_start + numel(spiking_idx)/2 - 0.5; %#ok<AGROW>
        row_labels_txt{end+1} = rd.label;                    %#ok<AGROW>
    end

    if size(big_mat, 1) >= 2
        big_mat = big_mat(1:end-2, :);
    end
end

%% ===== local_build_stack_sorted: sorted per-cycle stack (time domain) =====
function [sorted_mat, spike_overlay] = ...
        local_build_stack_sorted(all_roi, roi_order, ~, field_snips, field_spk, field_dt)
% Pool spiking cycles from ROIs in roi_order, sort by dt metric.

    all_snips = [];
    all_spk_t = {};
    all_dt    = [];

    for ii = 1:numel(roi_order)
        rd = all_roi(roi_order(ii));
        nP = size(rd.(field_snips), 1);
        spk_data = rd.(field_spk);
        for e = 1:nP
            if ~isempty(spk_data{e})  % only spiking cycles
                all_snips = [all_snips; rd.(field_snips)(e,:)]; %#ok<AGROW>
                all_dt    = [all_dt; rd.(field_dt)(e)];         %#ok<AGROW>
                all_spk_t{end+1,1} = spk_data{e};              %#ok<AGROW>
            end
        end
    end

    if isempty(all_snips)
        sorted_mat = zeros(0, size(all_roi(roi_order(1)).(field_snips), 2));
        spike_overlay = struct('row',{}, 'tval',{});
        return
    end

    dt_sort = all_dt;
    dt_sort(isnan(dt_sort)) = Inf;
    [~, sort_idx] = sort(dt_sort, 'ascend');
    sorted_mat = all_snips(sort_idx, :);

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

%% ===== local_build_stack_sorted_spiking: only cycles with spikes (time) =====
function [sorted_mat, spike_overlay] = ...
        local_build_stack_sorted_spiking(all_roi, roi_order, ~, field_snips, field_spk, field_dt)
% Same as local_build_stack_sorted but keeps only cycles with at least 1 spike.

    all_snips = [];
    all_spk_t = {};
    all_dt    = [];

    for ii = 1:numel(roi_order)
        rd = all_roi(roi_order(ii));
        nP = size(rd.(field_snips), 1);
        spk_data = rd.(field_spk);
        for e = 1:nP
            if ~isempty(spk_data{e})  % only keep cycles with spikes
                all_snips = [all_snips; rd.(field_snips)(e,:)]; %#ok<AGROW>
                all_dt    = [all_dt; rd.(field_dt)(e)];          %#ok<AGROW>
                all_spk_t{end+1,1} = spk_data{e};               %#ok<AGROW>
            end
        end
    end

    if isempty(all_snips)
        sorted_mat = [];
        spike_overlay = struct('row',{}, 'tval',{});
        return
    end

    [~, sort_idx] = sort(all_dt, 'ascend');
    sorted_mat = all_snips(sort_idx, :);

    spike_overlay = struct('row',{}, 'tval',{});
    for r = 1:numel(sort_idx)
        spk = all_spk_t{sort_idx(r)};
        for s = 1:numel(spk)
            spike_overlay(end+1).row = r;      %#ok<AGROW>
            spike_overlay(end).tval = spk(s);
        end
    end
end

%% ===== local_build_stack_phase: per-cycle stacked matrix (phase domain) =====
function [big_mat, spike_overlay, row_labels_pos, row_labels_txt] = ...
        local_build_stack_phase(all_roi, roi_order, nPhaseBins)
% Stack per-cycle phase-normalised dF/F for ROIs in roi_order (spiking cycles only).
% Uses dff_phase and spike_phase_per_cycle fields.

    big_mat        = [];
    spike_overlay  = struct('row',{}, 'phi',{});
    row_labels_pos = [];
    row_labels_txt = {};

    for ii = 1:numel(roi_order)
        rd      = all_roi(roi_order(ii));
        dff_ph  = rd.dff_phase;                  % [nCycles x nPhaseBins]
        spk_ph  = rd.spike_phase_per_cycle;
        valid_  = ~all(isnan(dff_ph), 2);
        valid_idx = find(valid_);

        % Filter to spiking cycles only (non-NaN rows with non-empty spike phase)
        spiking_idx = [];
        for cc_ = 1:numel(valid_idx)
            ci_ = valid_idx(cc_);
            if ci_ <= numel(spk_ph) && ~isempty(spk_ph{ci_})
                spiking_idx = [spiking_idx; cc_]; %#ok<AGROW>
            end
        end

        if isempty(spiking_idx)
            continue;
        end

        cyc_data  = dff_ph(valid_idx(spiking_idx), :);
        row_start = size(big_mat, 1) + 1;
        big_mat   = [big_mat; cyc_data; nan(2, nPhaseBins)]; %#ok<AGROW>

        for kk_ = 1:numel(spiking_idx)
            ci_ = valid_idx(spiking_idx(kk_));
            phi_ = spk_ph{ci_};
            row_in = row_start + kk_ - 1;
            for sp_ = 1:numel(phi_)
                spike_overlay(end+1).row = row_in;  %#ok<AGROW>
                spike_overlay(end).phi   = phi_(sp_);
            end
        end

        row_labels_pos(end+1) = row_start + numel(spiking_idx)/2 - 0.5; %#ok<AGROW>
        row_labels_txt{end+1} = rd.label;                               %#ok<AGROW>
    end

    if size(big_mat, 1) >= 2
        big_mat = big_mat(1:end-2, :);
    end
end

%% ===== local_build_stack_sorted_phase: sorted per-cycle stack (phase domain) =====
function [sorted_mat, spike_overlay] = ...
        local_build_stack_sorted_phase(all_roi, roi_order, ~)
% Pool spiking phase cycles from ROIs in roi_order, sort by nearest spike phase
% (minimum spike phase value in each cycle => closest to phase 0).

    all_snips  = [];
    all_spk_ph = {};
    all_dt_ph  = [];   % sorting metric: min spike phase per cycle (nearest to 0)

    for ii = 1:numel(roi_order)
        rd     = all_roi(roi_order(ii));
        dff_ph = rd.dff_phase;
        valid_ = ~all(isnan(dff_ph), 2);
        valid_idx = find(valid_);

        spk_ph = rd.spike_phase_per_cycle;

        for cc_ = 1:numel(valid_idx)
            ci_ = valid_idx(cc_);
            if ci_ <= numel(spk_ph) && ~isempty(spk_ph{ci_})  % only spiking cycles
                all_snips  = [all_snips;  dff_ph(ci_, :)];       %#ok<AGROW>
                all_dt_ph  = [all_dt_ph;  min(spk_ph{ci_})];    %#ok<AGROW>
                all_spk_ph = [all_spk_ph; {spk_ph{ci_}}];       %#ok<AGROW>
            end
        end
    end

    if isempty(all_snips)
        sorted_mat = [];
        spike_overlay = struct('row',{}, 'phi',{});
        return
    end

    [~, sort_idx] = sort(all_dt_ph, 'ascend');
    sorted_mat = all_snips(sort_idx, :);

    spike_overlay = struct('row',{}, 'phi',{});
    for r = 1:numel(sort_idx)
        spk = all_spk_ph{sort_idx(r)};
        if ~isempty(spk)
            for s = 1:numel(spk)
                spike_overlay(end+1).row = r;           %#ok<AGROW>
                spike_overlay(end).phi   = spk(s);
            end
        end
    end
end
