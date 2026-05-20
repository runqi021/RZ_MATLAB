% breath_polar_phase_diagnostic_260330.m
%  Breathing-only diagnostic: 6 panels comparing time, IBI, phase, and polar
%  representations of breathing waveforms.  No calcium/ROI data plotted.
%
%  Layout (2 rows x 4 columns):
%    Row 1: 1. Breath TIME overlay   2. Poincare (IBI)   3. IBI KDE    |
%    Row 2: 4. Phase overlay (Cart)  5. Polar waveform   6. Auto-coh   | 7. Breath PSD (tall)
%
%  Auto-detect mode:
%    - inputPath contains breath_peak_data.mat -> sub mode (single folder)
%    - Otherwise -> master mode (discover subfolders with breath_peak_data.mat)
%
%  Dependencies:
%    - Chronux (mtspectrumc, coherencyc)
%    - detect_session_fps.m
%  Note: confC is computed analytically from taper count, no calcium data needed.

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
% ======================================================================

%% ========================= USER PARAMETERS =========================
inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx         = [];             % session indices to skip (master mode)
fps_img         = 30;             % imaging frame rate (Hz)
nDrop           = 30;             % frames to drop from start
TW              = 5;              % multitaper time-bandwidth (K = 2*TW-1 tapers)
f_breath_search = [1.5, 3];      % Hz (search band for breathing peak; FWHM narrows for coherence)
nPhaseBins      = 200;            % bins for phase-normalized breath cycles
win_sec         = 0.5;            % half-window for peak-trig overlay (seconds)
ibi_xlim        = [0 2];          % x/y limits for Poincare and KDE
doSave          = true;           % true = save PDF, false = just plot
% ====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% 1 -- Auto-detect mode & discover sessions
if exist(fullfile(inputPath, '*breath_peak_data.mat'), 'file') == 2
    fprintf('Sub mode: processing single folder.\n');
    allMat = dir(fullfile(inputPath, '*breath_peak_data.mat'));
    isMasterMode = false; %#ok<NASGU>
else
    allMat = dir(fullfile(inputPath, '**', '*breath_peak_data.mat'));
    isMasterMode = true; %#ok<NASGU>

    if ~isempty(skipIdx)
        skipIdx = skipIdx(:)';
        keep = true(numel(allMat), 1);
        valid_skip = skipIdx(skipIdx >= 1 & skipIdx <= numel(allMat));
        keep(valid_skip) = false;
        allMat = allMat(keep);
    end
end

nSess = numel(allMat);
fprintf('Found %d session(s) with breath_peak_data.mat\n', nSess);
assert(nSess > 0, 'No breath_peak_data.mat found under %s', inputPath);

%% 2 -- Preallocate master collectors

% Time-domain breath snippets (plot 1)
all_breath_snip = [];

% Phase-domain breath cycles (plots 2, 3)
all_cycle_mat   = [];
all_cycle_sess  = [];

% IBI collectors (plots 2, 3)
ibi_all         = cell(nSess, 1);

% Breath PSD collectors (plot 7)
all_breath_psd  = cell(nSess, 1);
all_breath_freq = cell(nSess, 1);

% Auto-coherency (plot 6)
auto_r_vals     = nan(nSess, 1);
auto_th_vals    = nan(nSess, 1);
confC_vals      = nan(nSess, 1);

sessNames       = cell(nSess, 1);
validSess       = false(nSess, 1);
nValidFov       = 0;
t_w             = [];

%% 4 -- Main session loop
for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, sessName] = fileparts(folderPath);
    sessNames{kk} = sessName;

    fprintf('\n========== [%d/%d] %s ==========\n', kk, nSess, sessName);

    try
        %% ---- 4a. Auto-find required files ----
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

        %% ---- 4b. Auto-detect fps, get T from SAM ----
        fs = detect_session_fps(folderPath, fps_img);
        fs_breath = fs;

        SAMload = load(SAM_path, 'F');
        T = size(SAMload.F, 1);

        fprintf('T = %d frames\n', T);

        %% ---- 4c. Load breath peaks ----
        breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));
        insp_onsets = sort(breath_peak.insp_onset_idx(:));
        fp          = breath_peak.findpeak_params;

        insp_onsets(insp_onsets <= nDrop) = [];
        insp_onsets = insp_onsets - nDrop;

        if numel(insp_onsets) < 2
            warning('Fewer than 2 breath peaks after drop, skipping %s.', sessName);
            continue;
        end

        % IBI
        peak_t_sec = insp_onsets(:) / fs;
        ibi_all{kk} = diff(peak_t_sec);

        %% ---- 4d. Reconstruct breathing signal from DLC ----
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

        %% ---- 4e. Breathing PSD + FWHM (Chronux) ----
        params_breath.Fs     = fs_breath;
        params_breath.tapers = [TW, 2*TW - 1];
        params_breath.pad    = 0;
        params_breath.fpass  = [0.001, fs_breath/2];
        params_breath.err    = [2, 0.05];
        [Sb, fbC] = mtspectrumc(breath_dlc, params_breath);
        fbC = fbC(:); Sb = Sb(:);
        all_breath_psd{kk}  = Sb;
        all_breath_freq{kk} = fbC;

        mask_bp      = fbC >= f_breath_search(1) & fbC <= f_breath_search(2);
        [~, rel_idx] = max(Sb(mask_bp));
        idx_bp       = find(mask_bp, 1) + rel_idx - 1;
        half_max     = Sb(idx_bp) / 2;

        idx_lo = idx_bp;
        while idx_lo > 1          && Sb(idx_lo) > half_max, idx_lo = idx_lo - 1; end
        idx_hi = idx_bp;
        while idx_hi < numel(fbC) && Sb(idx_hi) > half_max, idx_hi = idx_hi + 1; end
        f_fwhm = [fbC(idx_lo), fbC(idx_hi)];
        f_fwhm = [max(f_fwhm(1), f_breath_search(1)), min(f_fwhm(2), f_breath_search(2))];
        freq_res = fbC(2) - fbC(1);
        min_bw   = 4 * freq_res;
        if diff(f_fwhm) < min_bw
            f_fwhm = [fbC(idx_bp) - min_bw/2, fbC(idx_bp) + min_bw/2];
            f_fwhm = [max(f_fwhm(1), f_breath_search(1)), min(f_fwhm(2), f_breath_search(2))];
        end
        fprintf('Breathing peak: %.3f Hz  FWHM: [%.3f, %.3f] Hz\n', ...
            fbC(idx_bp), f_fwhm(1), f_fwhm(2));

        %% ---- 4f. Peak-triggered breath snippets (plot 1) ----
        b_frames   = sort(insp_onsets(:));
        win_frames = round(win_sec * fs);

        % Common time grid (defined once)
        if isempty(t_w)
            nW  = 2 * round(win_sec * fps_img) + 1;
            t_w = linspace(-win_sec, win_sec, nW)';
        end
        nW = numel(t_w);

        edge_margin = max(win_frames, nDrop);
        insp_valid  = b_frames(b_frames > edge_margin & b_frames <= T - edge_margin);
        nInsp       = numel(insp_valid);

        if nInsp < 1
            warning('No valid breath peaks, skipping %s.', sessName);
            continue;
        end

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
        all_breath_snip = [all_breath_snip; breath_snip_sess]; %#ok<AGROW>

        %% ---- 4g. Phase-normalize each breathing cycle (plots 2, 3) ----
        phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
        phase_1cyc = phase_1cyc(1:end-1);
        nCyc_      = numel(insp_valid) - 1;
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

        %% ---- 4h. Auto-coherency of breathing (plot 4) ----
        params_coh.Fs     = fs;
        params_coh.tapers = [TW, 2*TW - 1];
        params_coh.pad    = 0;
        params_coh.fpass  = f_fwhm;
        params_coh.err    = [2, 0.01];

        if abs(fs_breath - fs) > 0.1
            breath_rs = resample(breath_dlc, round(fs), round(fs_breath));
        else
            breath_rs = breath_dlc;
        end
        T_coh     = min(T, numel(breath_rs));
        breath_rs = breath_rs(1:T_coh);

        % confC: analytical threshold from taper count and significance level
        K_tapers     = 2*TW - 1;
        p_conf       = params_coh.err(2);
        confC_val    = 1 - p_conf^(1/(K_tapers - 1));
        confC_vals(kk) = confC_val;

        % Auto-coherency
        [~, C_auto, phi_auto, ~, ~, ~, f_auto] = ...
            coherencyc(breath_rs, breath_rs, params_coh);
        f_auto = f_auto(:); C_auto = C_auto(:); phi_auto = phi_auto(:);

        mask_band = f_auto >= f_fwhm(1) & f_auto <= f_fwhm(2);
        if ~any(mask_band), mask_band = true(size(f_auto)); end
        auto_r_vals(kk)  = mean(C_auto(mask_band));
        auto_th_vals(kk) = angle(mean(exp(1i * (-phi_auto(mask_band)))));
        fprintf('  Auto-coh: r=%.4f  th=%.4f rad\n', auto_r_vals(kk), auto_th_vals(kk));

        %% ---- Mark session valid ----
        nValidFov = nValidFov + 1;
        validSess(kk) = true;
        fprintf('Session %d complete: %d peaks, %d cycles\n', kk, nInsp, nValidCyc);

    catch ME
        fprintf('ERROR in %s:\n  %s\n', sessName, ME.message);
        for me_k = 1:numel(ME.stack)
            fprintf('  in %s (line %d)\n', ME.stack(me_k).name, ME.stack(me_k).line);
        end
    end
end  % end session loop

%% 5 -- Post-processing
nValidSess = sum(validSess);
fprintf('\n========== Collection complete ==========\n');
fprintf('Valid sessions: %d / %d\n', nValidSess, nSess);
assert(nValidSess > 0, 'No valid sessions to plot.');

nTotalSnip   = size(all_breath_snip, 1);
nTotalCycles = size(all_cycle_mat, 1);
fprintf('Total breath snippets: %d\n', nTotalSnip);
fprintf('Total breathing cycles: %d\n', nTotalCycles);

confC_mean = mean(confC_vals(validSess), 'omitnan');

% Mean breathing cycle
mean_cycle = mean(all_cycle_mat, 1, 'omitnan');
phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
phase_1cyc = phase_1cyc(1:end-1);

% Normalize mean cycle to [0,1]
mean_cycle_01 = (mean_cycle - min(mean_cycle)) / (max(mean_cycle) - min(mean_cycle) + eps);

% Normalize each cycle to [0,1] for phase overlay (plot 2)
cycle_norm = all_cycle_mat;
for cc = 1:size(cycle_norm, 1)
    mn = min(cycle_norm(cc,:)); mx = max(cycle_norm(cc,:));
    if mx > mn, cycle_norm(cc,:) = (cycle_norm(cc,:) - mn) / (mx - mn); end
end

% Grand mean of [0,1]-normalized breath snippets (plot 1)
grand_mean_snip = mean(all_breath_snip, 1, 'omitnan');

% Mean auto-coherency point
auto_r_mean  = mean(auto_r_vals(validSess), 'omitnan');
auto_th_mean = angle(mean(exp(1i * auto_th_vals(validSess)), 'omitnan'));
fprintf('Mean auto-coherency: r=%.4f  th=%.4f rad\n', auto_r_mean, auto_th_mean);

%% 6 -- Session colors (for Poincare / KDE)
validSessIdx = find(validSess);
cmap_sess    = turbo(nValidSess);
sessColor    = nan(nSess, 3);
for ii = 1:nValidSess
    sessColor(validSessIdx(ii), :) = cmap_sess(ii, :);
end

%% 7 -- Figure: 6 panels (1 row x 6 columns)
sq       = 5;        % panel size (cm)
margin_l = 1.5;
margin_r = 1.0;
margin_b = 1.5;
margin_t = 2.0;
hgap     = 1.5;

nCols = 3;
nRows = 2;
vgap  = 1.5;

fig_w = margin_l + nCols*sq + (nCols-1)*hgap + margin_r;
fig_h = margin_b + nRows*sq + (nRows-1)*vgap + margin_t;

fig = figure('Units', 'centimeters', 'Position', [1 4 fig_w fig_h], ...
    'Color', 'w', 'Name', 'Breathing phase diagnostic');

col1_x = margin_l;
col2_x = col1_x + sq + hgap;
col3_x = col2_x + sq + hgap;
row1_y = margin_b + sq + vgap;   % top row
row2_y = margin_b;               % bottom row

%% ---- Plot 1: Breath overlay in TIME ----
ax1 = axes(fig, 'Units', 'centimeters', 'Position', [col1_x row1_y sq sq]);
hold(ax1, 'on');
plot(ax1, t_w, all_breath_snip', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
xline(ax1, 0, 'k--', 'LineWidth', 1);
hold(ax1, 'off');
xlim(ax1, [t_w(1), t_w(end)]);
ylim(ax1, [0, 1]);
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Breathing (norm)');
title(ax1, sprintf('Time overlay\n%d snippets, %d FOVs', nTotalSnip, nValidFov), 'FontSize', 7);
box(ax1, 'on'); pbaspect(ax1, [1 1 1]);

%% ---- Plot 2: Poincare (IBI_n vs IBI_{n+1}) ----
ax2 = axes(fig, 'Units', 'centimeters', 'Position', [col2_x row1_y sq sq]);
hold(ax2, 'on');
for ks = 1:nSess
    if ~validSess(ks), continue; end
    ibi = ibi_all{ks};
    if isempty(ibi) || numel(ibi) < 2, continue; end
    c_s = sessColor(ks, :);
    ran1 = (rand(size(ibi)) - 0.5) .* 2 ./ fs ./ 2;
    ran2 = (rand(size(ibi)) - 0.5) .* 2 ./ fs ./ 2;
    scatter(ax2, ibi(1:end-1)+ran1(1:end-1), ibi(2:end)+ran2(2:end), 6, c_s, ...
        'MarkerEdgeAlpha', 0.3);
end
plot(ax2, ibi_xlim, ibi_xlim, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
xlim(ax2, ibi_xlim); ylim(ax2, ibi_xlim);
axis(ax2, 'square');
xlabel(ax2, 'IBI_n (s)'); ylabel(ax2, 'IBI_{n+1} (s)');
title(ax2, 'Poincare', 'Interpreter','none', 'FontSize', 7);
box(ax2, 'on');
colormap(ax2, cmap_sess); caxis(ax2, [1 nValidSess]);

%% ---- Plot 3: IBI KDE ----
ax3 = axes(fig, 'Units', 'centimeters', 'Position', [col2_x row2_y sq sq]);
hold(ax3, 'on');
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
    plot(ax3, xi, fKDE, '-', 'Color', c_s, 'LineWidth', 1.5);
end
xlim(ax3, ibi_xlim);
xlabel(ax3, 'IBI (s)'); ylabel(ax3, 'Density');
title(ax3, 'IBI KDE', 'Interpreter','none', 'FontSize', 7);
box(ax3, 'on'); axis(ax3, 'square');

%% ---- Plot 4: Breath overlay in PHASE (Cartesian, [0, 2pi]) ----
ax4 = axes(fig, 'Units', 'centimeters', 'Position', [col1_x row2_y sq sq]);
hold(ax4, 'on');
phase_full = [phase_1cyc, phase_1cyc + 2*pi];
cycle_norm_dup = [cycle_norm, cycle_norm];
plot(ax4, phase_full, cycle_norm_dup', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.3);
%mean_01_dup = [mean_cycle_01, mean_cycle_01];
%plot(ax4, phase_full, mean_01_dup, 'k-', 'LineWidth', 1.5);
xline(ax4, 2*pi, 'k--', 'LineWidth', 0.8);
hold(ax4, 'off');
xlim(ax4, [0, 4*pi]);
ylim(ax4, [0, 1]);
set(ax4, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
    'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax4, 'Phase (rad)');
ylabel(ax4, 'Breathing (norm)');
title(ax4, sprintf('Phase overlay\n%d cycles, %d FOVs', nTotalCycles, nValidFov), 'FontSize', 7);
box(ax4, 'on'); pbaspect(ax4, [1 1 1]);

%% ---- Plot 5: Breath PSD (col 3, spans both rows) ----
tall_h = 2*sq + vgap;
ax7 = axes(fig, 'Units', 'centimeters', 'Position', [col3_x row2_y sq tall_h]);
hold(ax7, 'on');
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
    for ii = 1:size(psd_mat, 1)
        plot(ax7, f_common, 10*log10(psd_mat(ii,:)), '-', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.5);
    end
    m_lin = mean(psd_mat, 1, 'omitnan');
    % n_valid = sum(~isnan(psd_mat), 1);
    % sem_lin = std(psd_mat, 0, 1, 'omitnan') ./ sqrt(n_valid);
    plot(ax7, f_common, 10*log10(m_lin), 'k-', 'LineWidth', 1.5);
end
hold(ax7, 'off');
xlim(ax7, [0.01, 14]);
xlabel(ax7, 'Freq (Hz)');
ylabel(ax7, 'Power (dB)');
title(ax7, sprintf('Breath PSD (%d FOVs)', numel(valid_psd_idx)), 'FontSize', 7);
box(ax7, 'on');
set(ax7, 'XScale', 'log');

%% ---- Annotation ----
annotation(fig, 'textbox', [0.01, 0.94, 0.98, 0.05], ...
    'String', 'Row 1: time & Poincare    |    Row 2: phase & IBI KDE    |    Col 3: Chronux PSD', ...
    'EdgeColor', 'none', 'FontSize', 6, 'HorizontalAlignment', 'center', ...
    'FontAngle', 'italic', 'Interpreter', 'none');

%% 8 -- Save PDF
if doSave
    outDir  = inputPath;
    pdfPath = fullfile(outDir, 'breath_polar_phase_diagnostic.pdf');
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('Saved: %s\n', pdfPath);

    figPath = fullfile(outDir, 'breath_polar_phase_diagnostic.fig');
    savefig(fig, figPath);
    fprintf('Saved: %s\n', figPath);
end

fprintf('\nDone.\n');
