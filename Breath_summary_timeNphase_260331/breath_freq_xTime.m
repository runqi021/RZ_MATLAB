%% breath_freq_xTime.m
%  Per-session breathing peak frequency (Hz) vs raw-TIFF imaging start time.
%
%  For each session under inputPath that contains breath_peak_data.mat
%  (and a DLC CSV + raw .tif), this script:
%    1. Reconstructs the breathing trace from DLC + findpeak_params
%       (same recipe as Batch_breath_trig_combined_260328.m)
%    2. Computes a Chronux multitaper PSD
%    3. Picks the peak frequency inside f_breath_search
%    4. Reads SI.epoch = [YYYY MM DD HH MM SS] from the raw TIFF header
%
%  Output: one figure (peak Hz vs datetime, scatter), saved as .fig/.pdf/.mat
%  in inputPath/breath_freq_xTime/.
%
%  Auto-detect mode:
%    - inputPath itself contains breath_peak_data.mat       -> sub mode
%    - else, recurse and use every subfolder that has it    -> master mode
%
%  Dependencies: detect_session_fps, Chronux (mtspectrumc).
% -------------------------------------------------------------------------

clear; clc; close all;

%% PATH SETUP
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

%% USER PARAMETERS
inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
fps_img         = 30;            % fallback imaging fps
nDrop           = 30;            % frames to drop from start (matches batch)
TW              = 5;             % multitaper time-bandwidth (K = 2*TW-1 tapers)
fmin            = 0.01;          % Hz, PSD lower bound
fmax            = 14;            % Hz, PSD upper bound ([] -> fs/2)
f_breath_search = [1.5 3];     % Hz, peak search band
showLabels      = false;         % annotate each dot with session name
saveOutputs     = true;
% -------------------------------------------------------------------------

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% 1 -- Discover sessions
if exist(fullfile(inputPath, 'breath_peak_data.mat'), 'file') == 2
    fprintf('Sub mode: single folder.\n');
    allMat = dir(fullfile(inputPath, 'breath_peak_data.mat'));
else
    fprintf('Master mode: recursing under %s\n', inputPath);
    allMat = dir(fullfile(inputPath, '**', 'breath_peak_data.mat'));
end
nSess = numel(allMat);
assert(nSess > 0, 'No breath_peak_data.mat found under %s', inputPath);
fprintf('Found %d session(s).\n', nSess);

sessName  = strings(nSess, 1);
sessFolder= strings(nSess, 1);
peakFreq  = nan(nSess, 1);
sessTime  = NaT(nSess, 1);

%% 2 -- Per-session loop
for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, name] = fileparts(folderPath);
    sessName(kk)   = string(name);
    sessFolder(kk) = string(folderPath);
    fprintf('\n[%d/%d] %s\n', kk, nSess, name);

    % --- DLC CSV ---
    csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
    if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
    if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
    if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end
    if isempty(csv_hits)
        fprintf('  No DLC CSV, skipped.\n');
        continue
    end

    % --- breath_peak_data.mat (carries findpeak_params) ---
    BP = load(fullfile(allMat(kk).folder, allMat(kk).name));
    if ~isfield(BP, 'findpeak_params')
        fprintf('  No findpeak_params in breath_peak_data.mat, skipped.\n');
        continue
    end
    fp = BP.findpeak_params;

    % --- imaging fps (TIFF -> _meta.mat) ---
    fs = detect_session_fps(folderPath, fps_img);

    % --- Reconstruct breathing trace from DLC ---
    dlc_raw  = readmatrix(fullfile(folderPath, csv_hits(1).name), 'NumHeaderLines', 3);
    data_dlc = dlc_raw;
    data_dlc(1:min(nDrop, size(data_dlc,1)), :) = [];

    dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
    dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
    nActive  = numel(fp.dot_selection);
    traces   = NaN(size(data_dlc,1), nActive);
    for d = 1:nActive
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

    % --- Multitaper PSD ---
    if isempty(fmax), fmax_use = fs / 2; else, fmax_use = fmax; end
    params_breath        = struct();
    params_breath.Fs     = fs;
    params_breath.tapers = [TW, 2*TW - 1];
    params_breath.pad    = 0;
    params_breath.fpass  = [fmin, min(fmax_use, fs/2)];
    [Sb, fbC] = mtspectrumc(breath_dlc, params_breath);
    fbC = fbC(:); Sb = Sb(:);

    % --- Peak in search band ---
    mask_bp = fbC >= f_breath_search(1) & fbC <= f_breath_search(2);
    if ~any(mask_bp)
        fprintf('  No PSD bins in [%.2f %.2f] Hz, skipped.\n', f_breath_search);
        continue
    end
    [~, rel_idx] = max(Sb(mask_bp));
    idx_bp       = find(mask_bp, 1) + rel_idx - 1;
    peakFreq(kk) = fbC(idx_bp);

    % --- True imaging start time from raw TIFF epoch ---
    sessTime(kk) = read_tiff_epoch(folderPath);

    if isnat(sessTime(kk))
        fprintf('  peak=%.3f Hz   epoch=NaT (no TIFF metadata)\n', peakFreq(kk));
    else
        fprintf('  peak=%.3f Hz   epoch=%s\n', peakFreq(kk), datestr(sessTime(kk)));
    end
end

%% 3 -- Plot
valid = ~isnan(peakFreq) & ~isnat(sessTime);
nValid = sum(valid);
fprintf('\n%d/%d sessions usable for plot.\n', nValid, nSess);
assert(nValid > 0, 'No sessions had both a peak frequency and a TIFF epoch.');

[t_sorted, ord] = sort(sessTime(valid));
f_sorted        = peakFreq(valid);
f_sorted        = f_sorted(ord);
n_sorted        = sessName(valid);
n_sorted        = n_sorted(ord);

fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 18 9], ...
             'PaperUnits','centimeters', 'PaperSize',[18 9], ...
             'PaperPosition',[0 0 18 9]);
ax = axes(fig);
hold(ax, 'on');
plot(ax, t_sorted, f_sorted, 'o', ...
    'MarkerFaceColor', [0.20 0.45 0.80], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'LineWidth', 0.5);
if showLabels
    for ii = 1:nValid
        text(ax, t_sorted(ii), f_sorted(ii), "  " + n_sorted(ii), ...
            'FontSize', 7, 'Interpreter', 'none', ...
            'VerticalAlignment','middle', 'HorizontalAlignment','left');
    end
end
yline(ax, f_breath_search(1), 'k:', 'LineWidth', 0.5);
yline(ax, f_breath_search(2), 'k:', 'LineWidth', 0.5);
xlabel(ax, 'Imaging start time (raw TIFF epoch)');
ylabel(ax, 'Breathing peak frequency (Hz)');
title(ax, sprintf('Breathing peak frequency vs imaging time   (n=%d sessions)', nValid));
ax.FontName = 'Arial'; ax.FontSize = 9;
grid(ax, 'on'); box(ax, 'on');

%% 4 -- Save
if saveOutputs
    outDir = fullfile(inputPath, 'breath_freq_xTime');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    [~, parentName] = fileparts(inputPath);
    if isempty(parentName), parentName = 'session'; end
    outBase = char(fullfile(outDir, sprintf('breath_freq_xTime_%s', parentName)));

    savefig(fig, [outBase '.fig']);
    exportgraphics(fig, [outBase '.pdf'], ...
        'ContentType','vector', 'BackgroundColor','none');
    save([outBase '.mat'], 'sessName', 'sessFolder', 'peakFreq', 'sessTime', ...
        'inputPath', 'TW', 'f_breath_search', 'nDrop', 'fmin', 'fmax');
    fprintf('Saved:\n  %s.fig\n  %s.pdf\n  %s.mat\n', outBase, outBase, outBase);
end

fprintf('\nDone.\n');

%% ===================== LOCAL HELPERS =====================
function ep = read_tiff_epoch(folderPath)
% Read SI.epoch from the raw TIFF header. NaT if missing.
ep = NaT;
tif_hits = dir(fullfile(folderPath, '*.tif'));
if isempty(tif_hits), return; end
% Prefer raw TIFFs (no _preproc, _MC, _AVG suffix)
names  = {tif_hits.name};
is_raw = cellfun(@(n) ~contains(n, {'_preproc','_MC','_AVG'}, 'IgnoreCase', true), names);
order  = [find(is_raw), find(~is_raw)];
tif_hits = tif_hits(order);

for tt = 1:numel(tif_hits)
    tif_path = fullfile(tif_hits(tt).folder, tif_hits(tt).name);
    desc = '';
    try
        t = Tiff(tif_path, 'r');
        try
            desc = t.getTag('ImageDescription');
        catch
        end
        if isempty(desc)
            try
                desc = t.getTag('Software');
            catch
            end
        end
        t.close();
    catch
        try
            info1 = imfinfo(tif_path);
            info1 = info1(1);
            if isfield(info1,'ImageDescription'), desc = info1.ImageDescription; end
            if isempty(desc) && isfield(info1,'Software'), desc = info1.Software; end
        catch
        end
    end
    tok = regexp(desc, 'epoch\s*=\s*\[([^\]]+)\]', 'tokens', 'once');
    if ~isempty(tok)
        nums = str2double(strsplit(strtrim(tok{1})));
        if numel(nums) >= 6 && all(isfinite(nums(1:6)))
            ep = datetime(nums(1), nums(2), nums(3), ...
                          nums(4), nums(5), nums(6));
            return
        end
    end
end
end
