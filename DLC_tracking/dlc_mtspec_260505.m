% dlc_mtspec_260505.m
%   Multitaper PSD of one DLC bodypart axis (default: dot4 X) across all
%   runs found recursively under folderPath. Layout mimics row-1 of
%   breath_psd_per_session.m: one narrow panel per run, 10*log10(PSD) on
%   log-x, shared dB y-limits, only the first panel keeps axis labels.
%
%   Inputs:
%     <folderPath>/cam1_<ts>_runNNN/<...>DLC*shuffle*.csv  (one CSV per run)
%   Outputs (saved into folderPath):
%     mtspec_<bodypart>_<axis>_per_run.pdf
%     mtspec_<bodypart>_<axis>_per_run.mat
%
%   Dependencies: Chronux (mtspectrumc).

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

%% ========================= USER PARAMETERS ===========================
folderPath = "C:\Users\Admin\Desktop\260505_breathing_wt\cam1";
fps        = 50.7;            % Basler fps for this dataset
bodypart   = 'dot4';
axisName   = 'x';           % 'x' or 'y'
pmin       = 0.5;           % DLC likelihood threshold
fillMeth   = 'linear';      % fillmissing method

% Chronux multitaper
TW         = 6;
fpass      = [0.1 20];      % spectrum range to compute
xlim_plot  = [0.1 20];     % visible x-range (matches breath_psd_per_session)
peak_band  = [0.4 6];      % freq band for per-run peak detection

% Panel layout (cm) — wrap into a grid
panel_h     = 2.5;
panel_w     = 2.5 / 1.5 / 1.5;
panel_gap   = 0.4;       % horizontal gap
panel_vgap  = 1.0;       % vertical gap (a bit larger for run titles)
margin_l    = 1.5;
margin_r    = 0.5;
margin_b    = 1.2;
margin_t    = 0.8;
maxCols     = 20;        % panels per row
maxRows     = 15;        % cap on rows (warn + clip if exceeded)

doSave     = false;
% ======================================================================

%% FIND CSVs (recursive, exclude evaluation-results)
csv_pattern = '*DLC*shuffle*.csv';
Files = dir(fullfile(char(folderPath), '**', csv_pattern));
Files = Files(~contains({Files.folder}, 'evaluation-results'));
nRuns = numel(Files);
assert(nRuns > 0, 'No DLC CSVs found under %s', folderPath);

% Sort by run number
runNums = nan(nRuns, 1);
for k = 1:nRuns
    tok = regexp(Files(k).name, 'run(\d+)', 'tokens', 'once');
    if ~isempty(tok), runNums(k) = str2double(tok{1}); end
end
[~, ord] = sort(runNums);
Files   = Files(ord);
runNums = runNums(ord);
fprintf('Found %d CSVs (run %d ... %d)\n', nRuns, runNums(1), runNums(end));

%% LOAD + MULTITAPER PSD per run
params_mt.tapers = [TW, 2*TW - 1];
params_mt.pad    = 0;
params_mt.fpass  = fpass;
params_mt.err    = 0;

S_runs   = cell(nRuns, 1);
f_runs   = cell(nRuns, 1);
fps_runs = nan(nRuns, 1);

for k = 1:nRuns
    csv = fullfile(Files(k).folder, Files(k).name);

    % Per-run actual fps from timestamps.csv (Basler can drift far from configured fps)
    ts_csv = fullfile(Files(k).folder, 'timestamps.csv');
    if isfile(ts_csv)
        ts    = readmatrix(ts_csv);              % cols: frame_idx, camera_ts_ns, wall_s
        n_ts  = size(ts, 1);
        dur_s = (ts(end, 2) - ts(1, 2)) / 1e9;
        fps_run = (n_ts - 1) / dur_s;
    else
        fps_run = fps;
        warning('No timestamps.csv for run%03d, using fps=%g', runNums(k), fps);
    end
    fps_runs(k)  = fps_run;
    params_mt.Fs = fps_run;

    % Resolve column index for the requested bodypart from CSV header
    fid = fopen(csv, 'r');
    fgetl(fid);                     % scorer line
    bp_line = fgetl(fid);           % bodyparts line
    fclose(fid);
    bp_parts   = strsplit(strtrim(bp_line), ',');
    unique_bps = {};
    for c = 2:numel(bp_parts)
        nm = strtrim(bp_parts{c});
        if isempty(unique_bps) || ~strcmp(unique_bps{end}, nm)
            unique_bps{end+1} = nm; %#ok<AGROW>
        end
    end
    pid = find(strcmp(unique_bps, bodypart), 1);
    assert(~isempty(pid), 'bodypart %s missing in %s', bodypart, csv);

    data = readmatrix(csv, 'NumHeaderLines', 3);
    cx = 2 + 3*(pid - 1);
    cy = 3 + 3*(pid - 1);
    cp = 4 + 3*(pid - 1);
    if strcmpi(axisName, 'x')
        v = data(:, cx);
    else
        v = data(:, cy);
    end
    p = data(:, cp);

    % Likelihood mask -> fill -> demean
    bad = (p < pmin) | isnan(p);
    v(bad) = NaN;
    v = fillmissing(v, fillMeth, 'EndValues', 'nearest');
    v = v - mean(v, 'omitnan');

    [s, fk] = mtspectrumc(v, params_mt);
    S_runs{k} = s(:);
    f_runs{k} = fk(:);

    fprintf('  [%3d/%3d] run%03d  fps=%.2f  N=%d\n', k, nRuns, runNums(k), fps_run, numel(v));
end
fprintf('Per-run fps: mean=%.2f, min=%.2f, max=%.2f\n', ...
    mean(fps_runs, 'omitnan'), min(fps_runs), max(fps_runs));

%% Shared y-limits in dB across all runs
all_dB = [];
for k = 1:nRuns
    all_dB = [all_dB; 10*log10(S_runs{k})]; %#ok<AGROW>
end
all_dB = all_dB(isfinite(all_dB));
yLo = floor(min(all_dB) / 5) * 5;
yHi = ceil(max(all_dB) / 5) * 5;

%% Per-run peak frequency in peak_band, plus run start time from folder name
peak_freqs = nan(nRuns, 1);
run_times  = NaT(nRuns, 1);
for k = 1:nRuns
    fk     = f_runs{k};
    in_band = (fk >= peak_band(1)) & (fk <= peak_band(2));
    s_in   = S_runs{k}(in_band);
    f_in   = fk(in_band);
    [~, imax] = max(s_in);
    peak_freqs(k) = f_in(imax);

    [~, folderName] = fileparts(Files(k).folder);
    tok = regexp(folderName, '_(\d{8})_(\d{6})_run\d+', 'tokens', 'once');
    if ~isempty(tok)
        run_times(k) = datetime([tok{1} '_' tok{2}], 'InputFormat', 'yyyyMMdd_HHmmss');
    end
end

%% FIGURE: grid of narrow panels (maxCols per row), 10*log10 PSD, log-x
nRowsNeeded = ceil(nRuns / maxCols);
if nRowsNeeded > maxRows
    warning('nRuns=%d would need %d rows; clipping to first %d runs (%d rows × %d cols).', ...
        nRuns, nRowsNeeded, maxRows*maxCols, maxRows, maxCols);
    nRuns = maxRows * maxCols;
    nRowsNeeded = maxRows;
    Files   = Files(1:nRuns);
    runNums = runNums(1:nRuns);
    S_runs  = S_runs(1:nRuns);
end
nRowsUse = nRowsNeeded;

fig_w = margin_l + maxCols*panel_w + (maxCols - 1)*panel_gap + margin_r;
fig_h = margin_b + nRowsUse*panel_h + (nRowsUse - 1)*panel_vgap + margin_t;

fig1 = figure('Color', 'w', ...
    'Name', sprintf('PSD %s %s | %d runs', bodypart, axisName, nRuns), ...
    'Units', 'centimeters', 'Position', [0.5 0.5 fig_w fig_h], ...
    'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
    'PaperPosition', [0 0 fig_w fig_h]);
set(fig1, 'DefaultAxesFontSize', 6);

for k = 1:nRuns
    col = mod(k - 1, maxCols) + 1;
    row = ceil(k / maxCols);
    x0  = margin_l + (col - 1)*(panel_w + panel_gap);
    y0  = margin_b + (nRowsUse - row)*(panel_h + panel_vgap);

    ax = axes(fig1, 'Units', 'centimeters', 'Position', [x0 y0 panel_w panel_h]);
    plot(ax, f_runs{k}, 10*log10(S_runs{k}), 'k-', 'LineWidth', 0.8);
    xlim(ax, xlim_plot);  ylim(ax, [yLo yHi]);
    set(ax, 'XScale', 'log');  box(ax, 'on');
    title(ax, sprintf('run%03d', runNums(k)), 'FontSize', 5);

    is_first_col = (col == 1);
    is_last_row  = (row == nRowsUse) || (k > nRuns - maxCols);  % last row may be partial
    if is_first_col
        ylabel(ax, sprintf('%s_%s PSD (dB)', bodypart, axisName), 'Interpreter', 'none');
    else
        set(ax, 'YTickLabel', []);
    end
    if is_last_row
        xlabel(ax, 'Hz');
    else
        set(ax, 'XTickLabel', []);
    end
end

%% FIGURE 2: peak frequency over time (video timestamp from folder name)
fig2 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 22 9], 'PaperPositionMode', 'auto');
ax2 = axes(fig2);
plot(ax2, run_times, peak_freqs, 'ko-', 'LineWidth', 1, 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.30 0.50 0.80]);
xlabel(ax2, 'video start time');
ylabel(ax2, sprintf('peak frequency (Hz) in [%.1f, %.1f] Hz', peak_band(1), peak_band(2)));
title(ax2, sprintf('%s_%s peak freq vs time | n=%d runs', bodypart, axisName, nRuns), ...
    'Interpreter', 'none');
ylim(ax2, peak_band);
grid(ax2, 'on'); box(ax2, 'on');

%% SAVE
if doSave
    pdf1 = fullfile(char(folderPath), sprintf('mtspec_%s_%s_per_run.pdf',  bodypart, axisName));
    pdf2 = fullfile(char(folderPath), sprintf('mtspec_%s_%s_peakF_vs_t.pdf', bodypart, axisName));
    matp = fullfile(char(folderPath), sprintf('mtspec_%s_%s_per_run.mat',   bodypart, axisName));
    exportgraphics(fig1, pdf1, 'ContentType', 'vector', 'BackgroundColor', 'none');
    exportgraphics(fig2, pdf2, 'ContentType', 'vector', 'BackgroundColor', 'none');
    save(matp, 'S_runs', 'f_runs', 'fps_runs', 'runNums', 'run_times', 'peak_freqs', ...
        'peak_band', 'params_mt', 'folderPath', 'bodypart', 'axisName', '-v7.3');
    fprintf('\nSaved: %s\n', pdf1);
    fprintf('Saved: %s\n', pdf2);
    fprintf('Saved: %s\n', matp);
end

fprintf('Done — %d runs, %d freq bins (run 1).\n', nRuns, numel(f_runs{1}));
