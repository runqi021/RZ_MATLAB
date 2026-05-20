% whisker_LR_batch.m
%  Batch vL2 (L side) vs vR1 (R side) whisker comparison across all runs
%  in a folder. Signal = Y coordinate, demeaned. Two summary PDFs,
%  each a 5x6 grid of panels (one per run), vL2 in blue, vR1 in red.
%
%  Outputs (saved into folderPath):
%    whisker_<N>runs_traces.pdf   — Y demeaned, traceXlim window
%    whisker_<N>runs_spectra.pdf  — multitaper PSD, lin-lin
%    whisker_<N>runs_LR.mat       — cached signals + PSDs
%
%  Dependencies: Chronux (mtspectrumc).

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
% ======================================================================

%% ========================= USER PARAMETERS ===========================
folderPath  = "C:\260407_KA_electro_dorsal_whisking";

fps         = 90;

% Two whiskers to compare
L_pattern   = '*DLC_Resnet50_whisker_dorsalApr8*.csv';
L_point     = 'vL2';
R_pattern   = '*DLC_Resnet50_whisker_dorsal_RApr10*.csv';
R_point     = 'vR1';

pmin        = 0.5;
fillMeth    = 'linear';

% Multitaper PSD
TW          = 4;
fpass       = [5 45];

% Display
traceXlim   = [10 30];    % seconds, per-panel
nRowsGrid   = 5;
nColsGrid   = 6;

colL        = [0.00 0.45 0.85];   % blue
colR        = [0.85 0.10 0.10];   % red

doSave      = true;
% ======================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% 1 -- Pair L/R CSVs by run stem
L_files = dir(fullfile(folderPath, L_pattern));
R_files = dir(fullfile(folderPath, R_pattern));
assert(~isempty(L_files), 'No L CSVs in %s', folderPath);
assert(~isempty(R_files), 'No R CSVs in %s', folderPath);

getStem = @(n) regexprep(n, 'DLC_Resnet50.*$', '');
L_stems = arrayfun(@(f) string(getStem(f.name)), L_files);
R_stems = arrayfun(@(f) string(getStem(f.name)), R_files);

[run_stems, iL, iR] = intersect(L_stems, R_stems, 'stable');
nRuns = numel(run_stems);
fprintf('%d paired L/R runs\n', nRuns);
assert(nRuns > 0, 'No paired runs found.');

%% 2 -- Per run: load two CSVs, extract Y of the selected point, demean
L_y = cell(nRuns, 1);  R_y = cell(nRuns, 1);
t_all  = cell(nRuns, 1);
run_lbl = strings(nRuns, 1);
run_num = zeros(nRuns, 1);

for r = 1:nRuns
    fL = fullfile(L_files(iL(r)).folder, L_files(iL(r)).name);
    fR = fullfile(R_files(iR(r)).folder, R_files(iR(r)).name);

    [L_y{r}, t_all{r}] = local_point_y(fL, L_point, fps, pmin, fillMeth);
    [R_y{r}, ~       ] = local_point_y(fR, R_point, fps, pmin, fillMeth);

    run_lbl(r) = local_label(char(run_stems(r)));
    rn = regexp(char(run_stems(r)), 'run(\d+)', 'tokens', 'once');
    if ~isempty(rn), run_num(r) = str2double(rn{1}); end

    fprintf('[%2d/%2d] %s  (%d samples)\n', r, nRuns, run_lbl(r), numel(L_y{r}));
end

% Sort runs by run number so panels appear in order
[~, sortIdx] = sort(run_num);
L_y      = L_y(sortIdx);
R_y      = R_y(sortIdx);
t_all    = t_all(sortIdx);
run_lbl  = run_lbl(sortIdx);
run_stems = run_stems(sortIdx);

%% 3 -- Multitaper PSD per run per side (trim to common length)
params_mt.Fs     = fps;
params_mt.tapers = [TW, 2*TW - 1];
params_mt.pad    = 0;
params_mt.fpass  = fpass;
params_mt.err    = 0;

Lmin = min(min(cellfun(@numel, L_y)), min(cellfun(@numel, R_y)));

f_psd = [];
S_L = [];  S_R = [];
for r = 1:nRuns
    [sL, fk] = mtspectrumc(L_y{r}(1:Lmin), params_mt);
    [sR, ~ ] = mtspectrumc(R_y{r}(1:Lmin), params_mt);
    if isempty(f_psd)
        f_psd = fk(:);
        S_L = nan(numel(fk), nRuns);
        S_R = nan(numel(fk), nRuns);
    end
    S_L(:, r) = sL(:);
    S_R(:, r) = sR(:);
end

%% 4 -- Figure: traces grid
figT = figure('Units', 'centimeters', 'Position', [1 1 30 20], ...
    'Color', 'w', 'Name', sprintf('Traces %s vs %s (%d runs)', L_point, R_point, nRuns));

tl_T = tiledlayout(figT, nRowsGrid, nColsGrid, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_T, sprintf('Y demeaned — %s (blue) vs %s (red), %g-%gs', ...
    L_point, R_point, traceXlim(1), traceXlim(2)), 'FontSize', 10, 'FontWeight', 'bold');

for r = 1:nRuns
    ax = nexttile(tl_T);
    hold(ax, 'on');
    plot(ax, t_all{r}, L_y{r}, 'Color', colL, 'LineWidth', 0.3);
    plot(ax, t_all{r}, R_y{r}, 'Color', colR, 'LineWidth', 0.5);
    xlim(ax, traceXlim);
    title(ax, run_lbl(r), 'FontSize', 7);
    box(ax, 'on');
    if r > nRuns - nColsGrid
        xlabel(ax, 's');
    end
    if mod(r - 1, nColsGrid) == 0
        ylabel(ax, 'Y (px)');
    end
end

%% 5 -- Figure: PSD grid
figS = figure('Units', 'centimeters', 'Position', [1 1 30 20], ...
    'Color', 'w', 'Name', sprintf('PSDs %s vs %s (%d runs)', L_point, R_point, nRuns));

tl_S = tiledlayout(figS, nRowsGrid, nColsGrid, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_S, sprintf('Y PSD (lin-lin, TW=%g) — %s (blue) vs %s (red)', ...
    TW, L_point, R_point), 'FontSize', 10, 'FontWeight', 'bold');

for r = 1:nRuns
    ax = nexttile(tl_S);
    hold(ax, 'on');
    plot(ax, f_psd, S_L(:, r), 'Color', colL, 'LineWidth', 0.5);
    plot(ax, f_psd, S_R(:, r), 'Color', colR, 'LineWidth', 0.5);
    xlim(ax, fpass);
    title(ax, run_lbl(r), 'FontSize', 7);
    box(ax, 'on');
    if r > nRuns - nColsGrid
        xlabel(ax, 'Hz');
    end
    if mod(r - 1, nColsGrid) == 0
        ylabel(ax, 'px^2/Hz');
    end
end

%% 6 -- Save
if doSave
    tracesPdf  = fullfile(folderPath, sprintf('whisker_%druns_traces.pdf',  nRuns));
    spectraPdf = fullfile(folderPath, sprintf('whisker_%druns_spectra.pdf', nRuns));
    matPath    = fullfile(folderPath, sprintf('whisker_%druns_LR.mat',      nRuns));

    print(figT, '-dpdf', '-bestfit', char(tracesPdf));
    fprintf('Saved: %s\n', tracesPdf);
    print(figS, '-dpdf', '-bestfit', char(spectraPdf));
    fprintf('Saved: %s\n', spectraPdf);

    results = struct( ...
        'folderPath', char(folderPath), ...
        'run_stems',  {cellstr(run_stems)}, ...
        'run_lbl',    {cellstr(run_lbl)}, ...
        'fps',        fps, ...
        'L_point',    L_point, ...
        'R_point',    R_point, ...
        'L_y',        {L_y}, ...
        'R_y',        {R_y}, ...
        't_all',      {t_all}, ...
        'f_psd',      f_psd, ...
        'S_L',        S_L, ...
        'S_R',        S_R, ...
        'params_mt',  params_mt);
    save(matPath, '-struct', 'results', '-v7.3');
    fprintf('Saved: %s\n', matPath);
end

fprintf('Done (%d runs).\n', nRuns);

%% ==================== LOCAL FUNCTIONS ================================
function [y, tvec] = local_point_y(csv, point_name, fps, pmin, fillMeth)
% Return demeaned Y column of the requested DLC bodypart.
    bp_line = local_bodypart_line(csv);
    bp_parts = strsplit(strtrim(bp_line), ',');
    % Unique bodypart names in the order they appear (each name repeats 3x: x,y,p)
    unique_bps = {};
    for c = 2:numel(bp_parts)
        name = strtrim(bp_parts{c});
        if isempty(unique_bps) || ~strcmp(unique_bps{end}, name)
            unique_bps{end+1} = name; %#ok<AGROW>
        end
    end
    pid = find(strcmp(unique_bps, point_name), 1);
    assert(~isempty(pid), 'Bodypart %s not found in %s (found: %s)', ...
        point_name, csv, strjoin(unique_bps, ', '));

    data = readmatrix(csv, 'NumHeaderLines', 3);
    T    = size(data, 1);
    tvec = data(:, 1) / fps;

    cx = 2 + 3*(pid-1);  %#ok<NASGU>
    cy = 3 + 3*(pid-1);
    cp = 4 + 3*(pid-1);

    y = data(:, cy);
    p = data(:, cp);
    bad = (p < pmin) | isnan(p);
    y(bad) = NaN;
    y = fillmissing(y, fillMeth, 'EndValues', 'nearest');
    y = y - mean(y, 'omitnan');
end

function bp_line = local_bodypart_line(csv)
    fid = fopen(csv, 'r');
    fgetl(fid);               % scorer line
    bp_line = fgetl(fid);     % bodyparts line
    fclose(fid);
end

function lbl = local_label(stem)
% Build "runNNN [HH:MM]" label from a stem like basler_20260407_163056_run001
    rn = regexp(stem, 'run(\d+)', 'tokens', 'once');
    ts = regexp(stem, '_(\d{8})_(\d{6})_', 'tokens', 'once');
    if ~isempty(rn) && ~isempty(ts)
        hh = ts{2}(1:2);  mm = ts{2}(3:4);
        lbl = sprintf('run%s [%s:%s]', rn{1}, hh, mm);
    elseif ~isempty(rn)
        lbl = sprintf('run%s', rn{1});
    else
        lbl = stem;
    end
    lbl = string(lbl);
end
