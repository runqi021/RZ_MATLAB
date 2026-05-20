%% breath_cycle_overlay_260327.m
%  Overlay all breathing cycles normalized to phase [-2π, 2π].
%  Each peak-to-peak interval is one cycle interpolated to [0, 2π],
%  then duplicated to show [-2π, 2π].
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath  = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_2000f_00001";
FPS         = 30;
nPhaseBins  = 200;
saveFigs    = true;
% ──────────────────────────────────────────────────────────────────────────

%% ============================== LOAD ==============================

% --- *_dFF.mat (just for T and fs) ---
dffEntry = dir(fullfile(folderPath, '*_dFF.mat'));
assert(numel(dffEntry) == 1, "Expected 1 *_dFF.mat, found %d.", numel(dffEntry));
[~, stem] = fileparts(dffEntry(1).name);
D = load(fullfile(dffEntry.folder, dffEntry.name), 'dFF', 'params');
T = size(D.dFF, 1);
if isfield(D,'params') && isfield(D.params,'FPS')
    fs = double(D.params.FPS);
else
    fs = FPS;
    warning("params.FPS not found — using fallback FPS=%g", FPS);
end

% --- breath_peak_data.mat ---
bpFile = fullfile(folderPath, 'breath_peak_data.mat');
assert(isfile(bpFile), "breath_peak_data.mat not found.");
BP          = load(bpFile);
insp_onsets = sort(BP.insp_onset_idx(:));
fp          = BP.findpeak_params;

% --- DLC CSV → breath signal [T × 1] ---
csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end
assert(~isempty(csv_hits), "No DLC CSV found.");
dlc_raw  = readmatrix(fullfile(folderPath, csv_hits(1).name), 'NumHeaderLines', 3);
dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
dot_map  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
nDots    = numel(fp.dot_selection);
traces   = NaN(size(dlc_raw,1), nDots);
for d = 1:nDots
    di  = dot_map.(fp.dot_selection{d});
    xc  = dlc_raw(:, dot_cols(di,1));
    yc  = dlc_raw(:, dot_cols(di,2));
    pc  = dlc_raw(:, dot_cols(di,3));
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
        sig = fillmissing(sig, 'linear', 'EndValues','nearest');
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

%% ============================== PHASE NORMALIZE EACH CYCLE ==============================
b_frames = sort(insp_onsets(:));
nPeaks   = numel(b_frames);
nCycles  = nPeaks - 1;

phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
phase_1cyc = phase_1cyc(1:end-1);   % [0, 2π)

cycle_mat = nan(nCycles, nPhaseBins);
for c = 1:nCycles
    i1 = b_frames(c);
    i2 = b_frames(c + 1);
    if i1 >= 1 && i2 <= T
        chunk = breath_dlc(i1:i2);
        orig_phase = linspace(0, 2*pi, numel(chunk));
        cycle_mat(c, :) = interp1(orig_phase, chunk, phase_1cyc, 'linear');
    end
end

valid = ~all(isnan(cycle_mat), 2);
cycle_mat = cycle_mat(valid, :);
nValid = size(cycle_mat, 1);
mean_cycle = mean(cycle_mat, 1, 'omitnan');

% duplicate to [-2π, 2π]
phase_full = [phase_1cyc - 2*pi, phase_1cyc];
cycle_dup  = [cycle_mat, cycle_mat];
mean_dup   = [mean_cycle, mean_cycle];

Delta_k = diff(b_frames(logical([valid; false]))) / fs;
fprintf("%d peaks → %d valid cycles  mean period = %.3f s (%.1f Hz)\n", ...
    nPeaks, nValid, mean(Delta_k), 1/mean(Delta_k));

%% ============================== PLOT ==============================
fig = figure('Color', 'w', 'Name', 'breath_cycle_overlay');
set(fig, 'DefaultAxesFontSize', 10, 'DefaultTextFontSize', 11);

ax = axes(fig);
hold(ax, 'on');
plot(ax, phase_full, cycle_dup', 'Color', [0.75 0.75 0.75 0.3], 'LineWidth', 0.5);
plot(ax, phase_full, mean_dup, 'k', 'LineWidth', 2);
xline(ax, 0, 'r--', 'LineWidth', 1);
xline(ax, -2*pi, 'r--', 'LineWidth', 1);
xline(ax, 2*pi, 'r--', 'LineWidth', 1);
hold(ax, 'off');

xlim(ax, [-2*pi, 2*pi]);
set(ax, 'XTick', [-2*pi -pi 0 pi 2*pi], ...
    'XTickLabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
xlabel(ax, 'Phase (rad)');
ylabel(ax, 'Breathing (z-scored)');
title(ax, sprintf('%s  —  %d cycles overlaid', strrep(stem,'_','\_'), nValid), 'FontSize', 11);
box(ax, 'on');

drawnow;

%% ============================== SAVE ==============================
if saveFigs
    outDir = fullfile(folderPath, 'breath_trig_heatmap');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    outFig = fullfile(outDir, 'breath_cycle_overlay.fig');
    savefig(fig, outFig);

    outPdf = fullfile(outDir, 'breath_cycle_overlay.pdf');
    exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','none');

    fprintf("Saved:\n  %s\n  %s\n", outFig, outPdf);
end
