% whisker_phase_single.m
%  Single-file whisker DLC analysis: trace + spectrum + spectrogram.
%
%  Layout: (nWhiskers+1) columns x 3 rows
%    Each column = one whisker point, last column = average
%    Row 1: Time trace (detrended, z-scored)
%    Row 2: Multitaper PSD (Chronux)
%    Row 3: Multitaper spectrogram (Chronux mtspecgramc)
%
%  Dependencies:
%    - Chronux (mtspectrumc, mtspecgramc)
%    - Statistics Toolbox (pca)

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
% ======================================================================

%% ========================= USER PARAMETERS ===========================
% ---- Input file ----
dlc_csv = "C:\260407_KA_electro_dorsal_whisking\basler_20260407_164056_run003DLC_Resnet50_whisker_dorsalApr8shuffle1_snapshot_best-10.csv";

% ---- Acquisition ----
fps = 90;                    % video frame rate (Hz)

% ---- DLC bodyparts ----
point_ids       = [1 2 3 4]; % DLC bodypart IDs (1-indexed: 1=vL1, 2=vL2, ...)
pmin            = 0.5;       % likelihood threshold
fill_method     = 'linear';  % 'linear' or 'pchip'

% ---- Signal type ----
% 'pc1'  = first principal component of (x,y) per point (signed, max variance)
% 'y'    = Y coordinate (demeaned)
% 'x'    = X coordinate (demeaned)
signal_type     = 'pc1';

% ---- Multitaper PSD (Chronux) ----
TW              = 8;         % time-bandwidth product
fpass           = [0.1 44];  % frequency range (Hz)

% ---- Spectrogram (Chronux mtspecgramc) ----
spec_win_sec    = 0.5;       % window length (s)
spec_step_sec   = 0.05;      % step size (s)
spec_fpass      = [1 44];    % frequency range for spectrogram
spec_TW         = 3;

% ---- Plotting ----
psd_dB          = true;      % true = dB scale, false = linear
spec_clim_prc   = [1 99];    % percentile clamp for spectrogram color
trace_xlim      = [10 30];   % 20s representative window for row 1 (seconds)

% ---- Save ----
doSave          = true;
% ======================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% 1 -- Validate
fs  = fps;
[folderPath, csvStem] = fileparts(dlc_csv);
fprintf('CSV: %s\nfps: %d Hz\n', dlc_csv, fps);

%% 2 -- Load DLC CSV
data = readmatrix(dlc_csv, 'NumHeaderLines', 3);
nTriplets = floor((size(data, 2) - 1) / 3);
fprintf('%d frames, %d bodyparts\n', size(data, 1), nTriplets);

assert(all(point_ids >= 1 & point_ids <= nTriplets), ...
    'point_ids outside [1, %d]', nTriplets);

frame = data(:, 1);
t     = frame / fps;
T     = numel(t);

%% 3 -- Extract per-point x, y; likelihood filter + fill
nPts  = numel(point_ids);
X_all = nan(T, nPts);
Y_all = nan(T, nPts);

for k = 1:nPts
    pid = point_ids(k);
    cx  = 2 + 3 * (pid - 1);
    cy  = 3 + 3 * (pid - 1);
    cp  = 4 + 3 * (pid - 1);

    x = data(:, cx);
    y = data(:, cy);
    p = data(:, cp);

    bad = (p < pmin) | isnan(p);
    x(bad) = NaN;  y(bad) = NaN;
    x = fillmissing(x, fill_method);
    y = fillmissing(y, fill_method);

    X_all(:, k) = x;
    Y_all(:, k) = y;
end

%% 4 -- Build per-whisker signals [T x nPts], detrend + z-score
switch lower(signal_type)
    case 'y'
        RAW = Y_all;
    case 'x'
        RAW = X_all;
    case 'pc1'
        RAW = nan(T, nPts);
        for k = 1:nPts
            xy = [detrend(X_all(:,k)), detrend(Y_all(:,k))];
            [coeff, ~] = pca(xy);
            pc = xy * coeff(:,1);
            % fix sign ambiguity: align PC1 with Y direction
            if coeff(2,1) < 0, pc = -pc; end
            RAW(:,k) = pc;
        end
    otherwise
        error('signal_type must be ''pc1'', ''x'', or ''y''');
end

sig_each = nan(T, nPts);
for k = 1:nPts
    s = detrend(RAW(:, k) - mean(RAW(:, k)));
    sd = std(s);
    if sd > 0, s = s / sd; end
    sig_each(:, k) = s;
end

% Average across z-scored whiskers, re-normalize
sig_avg = mean(sig_each, 2);
sig_avg = sig_avg / std(sig_avg);

% All signals: columns 1..nPts = individual, column nPts+1 = average
sig_all = [sig_each, sig_avg];
nCols   = nPts + 1;

fprintf('%d whiskers + average, signal_type = %s, z-scored\n', nPts, signal_type);

%% 5 -- Multitaper PSD for each (on raw detrended signal, no filter)
params_mt.Fs     = fs;
params_mt.tapers = [TW, 2 * TW - 1];
params_mt.pad    = 0;
params_mt.fpass  = fpass;
params_mt.err    = [2, 0.05];

S_all     = [];
f_psd     = [];
S_err_all = cell(nCols, 1);

for k = 1:nCols
    [Sk, fk, Sek] = mtspectrumc(sig_all(:, k), params_mt);
    if isempty(S_all)
        S_all = nan(numel(Sk), nCols);
        f_psd = fk(:);
    end
    S_all(:, k)   = Sk(:);
    S_err_all{k}  = Sek;
end

%% 6 -- Spectrogram for each (Chronux mtspecgramc, on raw signal)
params_sg.Fs     = fs;
params_sg.tapers = [spec_TW, 2 * spec_TW - 1];
params_sg.pad    = 0;
params_sg.fpass  = spec_fpass;

movingwin = [spec_win_sec, spec_step_sec];

SG_all = cell(nCols, 1);
t_sg   = [];
f_sg   = [];

for k = 1:nCols
    [Sg, tg, fg] = mtspecgramc(sig_all(:, k), movingwin, params_sg);
    SG_all{k} = Sg;
    if isempty(t_sg)
        t_sg = tg(:);
        f_sg = fg(:);
    end
end

% Shared color limits
all_sg_dB = [];
for k = 1:nCols
    all_sg_dB = [all_sg_dB; 10*log10(SG_all{k}(:))]; %#ok<AGROW>
end
sg_clim = prctile(all_sg_dB, spec_clim_prc);

%% 7 -- Column labels from CSV header
col_labels = cell(1, nCols);
fid = fopen(dlc_csv, 'r');
fgetl(fid);
bp_line = fgetl(fid);
fclose(fid);
bp_parts = strsplit(bp_line, ',');

for k = 1:nPts
    pid = point_ids(k);
    col_idx = 2 + 3*(pid-1);
    if col_idx <= numel(bp_parts)
        col_labels{k} = bp_parts{col_idx};
    else
        col_labels{k} = sprintf('pt%d', pid);
    end
end
col_labels{nCols} = 'Average';

%% 8 -- Figure: 3 rows x nCols columns
sq_w     = 5.6;
sq_h     = 3.0;
margin_l = 1.8;
margin_r = 0.5;
margin_b = 1.2;
margin_t = 1.5;
hgap     = 1.2;
vgap     = 0.8;

fig_w = margin_l + nCols * sq_w + (nCols - 1) * hgap + margin_r;
fig_h = margin_b + 3 * sq_h + 2 * vgap + margin_t;

fig = figure('Units', 'centimeters', 'Position', [1 1 fig_w fig_h], ...
    'Color', 'w', 'Name', 'Whisker DLC: trace + PSD + spectrogram');

col_whisk = lines(nPts);
col_color = [col_whisk; 0 0 0];

for c = 1:nCols
    x0 = margin_l + (c - 1) * (sq_w + hgap);

    % ---- Row 1: time trace ----
    y1 = margin_b + 2 * (sq_h + vgap);
    ax_tr = axes(fig, 'Units', 'centimeters', 'Position', [x0 y1 sq_w sq_h]);
    plot(ax_tr, t, sig_all(:, c), 'Color', col_color(c,:), 'LineWidth', 0.3);
    xlim(ax_tr, trace_xlim);
    title(ax_tr, col_labels{c}, 'FontSize', 8, 'Color', col_color(c,:));
    box(ax_tr, 'on');
    if c == 1
        ylabel(ax_tr, sprintf('PC1 (z)'));
    else
        set(ax_tr, 'YTickLabel', []);
    end

    % ---- Row 2: PSD ----
    y2 = margin_b + 1 * (sq_h + vgap);
    ax_psd = axes(fig, 'Units', 'centimeters', 'Position', [x0 y2 sq_w sq_h]);
    hold(ax_psd, 'on');

    Se = S_err_all{c};
    Sv = S_all(:, c);
    if psd_dB
        plot(ax_psd, f_psd, 10*log10(Sv), 'Color', col_color(c,:), 'LineWidth', 1.0);
        fill(ax_psd, [f_psd; flipud(f_psd)], ...
            [10*log10(Se(1,:))'; flipud(10*log10(Se(2,:))')], ...
            col_color(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    else
        plot(ax_psd, f_psd, Sv, 'Color', col_color(c,:), 'LineWidth', 1.0);
        fill(ax_psd, [f_psd; flipud(f_psd)], ...
            [Se(1,:)'; flipud(Se(2,:)')], ...
            col_color(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end

    hold(ax_psd, 'off');
    set(ax_psd, 'XScale', 'log');
    xlim(ax_psd, fpass);
    box(ax_psd, 'on');
    if c == 1
        if psd_dB, ylabel(ax_psd, 'Power (dB)'); else, ylabel(ax_psd, 'Power'); end
    else
        set(ax_psd, 'YTickLabel', []);
    end
    if c == 1, xlabel(ax_psd, 'Hz'); end

    % ---- Row 3: Spectrogram ----
    y3 = margin_b;
    ax_sg = axes(fig, 'Units', 'centimeters', 'Position', [x0 y3 sq_w sq_h]);
    Sg_dB = 10 * log10(SG_all{c});
    imagesc(ax_sg, t_sg, f_sg, Sg_dB.');
    axis(ax_sg, 'xy');
    caxis(ax_sg, sg_clim);
    xlim(ax_sg, [t_sg(1) t_sg(end)]);
    ylim(ax_sg, spec_fpass);
    box(ax_sg, 'on');
    if c == 1
        ylabel(ax_sg, 'Freq (Hz)');
        xlabel(ax_sg, 'Time (s)');
    else
        set(ax_sg, 'YTickLabel', []);
    end
    if c == nCols
        cb = colorbar(ax_sg, 'Location', 'eastoutside');
        cb.Label.String = 'dB';
        cb.Label.FontSize = 6;
    end
end

% Supertitle
sessStr = regexprep(csvStem, 'DLC.*$', '');
sessStr = strrep(sessStr, '_', ' ');
annotation(fig, 'textbox', [0.01 0.95 0.98 0.04], ...
    'String', sprintf('%s  |  %s  |  %d whiskers + avg', sessStr, signal_type, nPts), ...
    'EdgeColor', 'none', 'FontSize', 8, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'none');

%% 9 -- Save
if doSave
    shortName = char(regexprep(csvStem, 'DLC.*$', ''));
    outStem = fullfile(folderPath, [shortName 'whiskerSpec']);

    pdfPath = char([outStem '.pdf']);
    print(fig, '-dpdf', '-bestfit', pdfPath);
    fprintf('Saved: %s\n', pdfPath);

    matPath = char([outStem '.mat']);
    results = struct( ...
        'dlc_csv',      dlc_csv, ...
        'fps',          fps, ...
        'signal_type',  signal_type, ...
        'point_ids',    point_ids, ...
        'col_labels',   {col_labels}, ...
        't',            t, ...
        'sig_each',     sig_each, ...
        'sig_avg',      sig_avg, ...
        'f_psd',        f_psd, ...
        'S_all',        S_all, ...
        'S_err_all',    {S_err_all}, ...
        't_sg',         t_sg, ...
        'f_sg',         f_sg, ...
        'SG_all',       {SG_all}, ...
        'params_mt',    params_mt, ...
        'spec_win',     movingwin);
    save(matPath, '-struct', 'results', '-v7.3');
    fprintf('Saved: %s\n', matPath);
end

fprintf('Done.\n');
