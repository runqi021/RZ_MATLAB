%% breath_trig_heatmap_260327.m
%  Per-ROI breath-triggered dF/F heatmap — two panels:
%    Left:  time-domain (sorted by nearest Ca spike)
%    Right: phase-normalized per cycle, duplicated to [-2π, 2π]
%  Saves .fig + .pdf per ROI into folderPath/breath_trig_heatmap/.
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath      = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_2000f_00001";
FPS             = 30;       % fallback imaging fps
prcLim          = [0.1 99.9]; % heatmap color percentile limits
nPhaseBins      = 100;      % phase resolution per cycle
signedSort      = true;     % false = sort by |dt|, true = signed dt (neg=Ca before breath)
saveFigs        = true;
% ──────────────────────────────────────────────────────────────────────────

%% ============================== LOAD ==============================

% --- *_dFF.mat ---
dffEntry = dir(fullfile(folderPath, '*_dFF.mat'));
assert(numel(dffEntry) == 1, "Expected 1 *_dFF.mat, found %d.", numel(dffEntry));
dffMatPath = fullfile(dffEntry.folder, dffEntry.name);
[~, stem]  = fileparts(dffMatPath);
D   = load(dffMatPath);
dFF = double(D.dFF);   % [T × N]
[T, N] = size(dFF);
if isfield(D,'params') && isfield(D.params,'FPS')
    fs = double(D.params.FPS);
else
    fs = FPS;
    warning("params.FPS not found — using fallback FPS=%g", FPS);
end

% --- ca_spike_data.mat ---
spikeFile = fullfile(folderPath, 'ca_spike_data.mat');
assert(isfile(spikeFile), "ca_spike_data.mat not found.");
Ks      = load(spikeFile);
roiList = find(logical(Ks.ifSpike(:)'));
caData  = Ks.roi_spikes;
assert(~isempty(roiList), "No spiking ROIs found.");

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

%% ============================== INSPIRATION TIMING ==============================
b_frames   = sort(insp_onsets(:));
nPeaks     = numel(b_frames);
Delta_k    = diff(b_frames) / fs;
win_sec    = mean(Delta_k) * 2.5;
win_frames = round(win_sec * fs);
t_snip     = (-win_frames:win_frames) / fs;

% phase grid for one cycle [0, 2π), duplicated to [-2π, 2π)
phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
phase_1cyc = phase_1cyc(1:end-1);   % drop duplicate endpoint
phase_full = [phase_1cyc - 2*pi, phase_1cyc];  % [-2π, 2π)

% valid onsets for phase: need prev and next peak within data
valid_ph_idx = [];
for b = 2:nPeaks-1
    if b_frames(b-1) >= 1 && b_frames(b+1) <= T
        valid_ph_idx(end+1) = b; %#ok<SAGROW>
    end
end
nValid_ph_all = numel(valid_ph_idx);

% phase-normalize spiking ROIs (for summary heatmap)
nSpk = numel(roiList);
avg_phase_all = nan(nSpk, nPhaseBins);
for ri = 1:nSpk
    roi = roiList(ri);
    dff_roi_all = dFF(:, roi);
    phase_snips_all = nan(nValid_ph_all, nPhaseBins);
    for b = 1:nValid_ph_all
        bi = valid_ph_idx(b);
        bp_curr = b_frames(bi);
        bp_next = b_frames(bi + 1);
        chunk = dff_roi_all(bp_curr:bp_next);
        orig_phase = linspace(0, 2*pi, numel(chunk));
        phase_snips_all(b, :) = interp1(orig_phase, chunk, phase_1cyc, 'linear');
    end
    avg_phase_all(ri, :) = mean(phase_snips_all, 1, 'omitnan');
end

fprintf("T=%d  N=%d  %d spiking ROIs  %d insp onsets\n", ...
    T, N, numel(roiList), nPeaks);

cmap_roi = lines(numel(roiList));

%% ============================== PER-ROI HEATMAP ==============================
for k = 1:numel(roiList)
    roi_id  = roiList(k);
    ca_ev   = caData(roi_id).spike_idx(:);
    dff_roi = dFF(:, roi_id);

    % =================================================================
    % LEFT: time-domain snippets
    % =================================================================
    nBreath  = nPeaks;
    snippets = nan(nBreath, 2*win_frames + 1);
    for b = 1:nBreath
        bp = b_frames(b);
        i1 = bp - win_frames;
        i2 = bp + win_frames;
        if i1 >= 1 && i2 <= T
            snippets(b, :) = dff_roi(i1:i2);
        end
    end
    valid_time     = ~all(isnan(snippets), 2);
    nValid_time    = sum(valid_time);
    valid_idx_time = find(valid_time);
    b_valid_time   = b_frames(valid_idx_time);

    % sort by dt_nearest_ca
    dt_near_time = nan(nValid_time, 1);
    for b = 1:nValid_time
        bp = b_valid_time(b);
        before = ca_ev(ca_ev <  bp);
        after  = ca_ev(ca_ev >= bp);
        dt_bef = Inf;  dt_aft = Inf;
        if ~isempty(before), dt_bef = (bp - before(end)) / fs; end
        if ~isempty(after),  dt_aft = (after(1) - bp)    / fs; end
        if signedSort
            if dt_bef <= dt_aft, dt_near_time(b) = -dt_bef;
            else,                dt_near_time(b) =  dt_aft; end
        else
            dt_near_time(b) = min(dt_bef, dt_aft);
        end
    end
    dt_sort = dt_near_time;
    dt_sort(isnan(dt_sort)) = Inf;
    [~, sort_t] = sort(dt_sort, 'ascend');

    snip_time_sorted = snippets(valid_idx_time(sort_t), :);
    cHigh_t = prctile(snip_time_sorted(:), prcLim(2));

    % =================================================================
    % RIGHT: phase-normalized per cycle, duplicated to [-2π, 2π]
    % =================================================================
    % valid onsets: need both previous and next peak within data
    valid_ph = false(nPeaks, 1);
    for b = 2:nPeaks-1
        if b_frames(b-1) >= 1 && b_frames(b+1) <= T
            valid_ph(b) = true;
        end
    end
    valid_idx_ph = find(valid_ph);
    nValid_ph    = numel(valid_idx_ph);
    b_valid_ph   = b_frames(valid_idx_ph);

    % interpolate each trial: prev→current = one cycle, current→next = one cycle
    phase_snips = nan(nValid_ph, nPhaseBins);
    for b = 1:nValid_ph
        bi   = valid_idx_ph(b);
        bp_prev = b_frames(bi - 1);
        bp_curr = b_frames(bi);
        bp_next = b_frames(bi + 1);
        % current→next cycle: frames bp_curr to bp_next, map to [0, 2π]
        chunk = dff_roi(bp_curr:bp_next);
        nPts  = numel(chunk);
        orig_phase = linspace(0, 2*pi, nPts);
        phase_snips(b, :) = interp1(orig_phase, chunk, phase_1cyc, 'linear');
    end

    % sort by dt_nearest_ca (same logic)
    dt_near_ph = nan(nValid_ph, 1);
    for b = 1:nValid_ph
        bp = b_valid_ph(b);
        before = ca_ev(ca_ev <  bp);
        after  = ca_ev(ca_ev >= bp);
        dt_bef = Inf;  dt_aft = Inf;
        if ~isempty(before), dt_bef = (bp - before(end)) / fs; end
        if ~isempty(after),  dt_aft = (after(1) - bp)    / fs; end
        if signedSort
            if dt_bef <= dt_aft, dt_near_ph(b) = -dt_bef;
            else,                dt_near_ph(b) =  dt_aft; end
        else
            dt_near_ph(b) = min(dt_bef, dt_aft);
        end
    end
    dt_sort_ph = dt_near_ph;
    dt_sort_ph(isnan(dt_sort_ph)) = Inf;
    [~, sort_p] = sort(dt_sort_ph, 'ascend');

    phase_sorted = phase_snips(sort_p, :);
    % duplicate: tile one cycle into [-2π, 2π]
    phase_sorted_dup = [phase_sorted, phase_sorted];
    cHigh_p = prctile(phase_sorted_dup(:), prcLim(2));

    % use same color scale for both panels
    cHigh = max(cHigh_t, cHigh_p);

    % =================================================================
    % FIGURE: two panels side by side
    % =================================================================
    fig = figure('Color', 'w', ...
        'Name', sprintf('breath_trig_heatmap ROI#%02d', roi_id));
    set(fig, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 10);
    sgtitle(fig, sprintf('%s  —  ROI #%d', strrep(stem,'_','\_'), roi_id), 'FontSize', 10);

    % --- Left: time-domain ---
    ax1 = subplot(1,2,1, 'Parent', fig);
    imagesc(ax1, t_snip, 1:nValid_time, snip_time_sorted);
    set(ax1, 'YDir', 'reverse');
    colormap(ax1, flipud(gray(256)));
    caxis(ax1, [0 cHigh]);
    hold(ax1, 'on');
    xline(ax1, 0, 'r--', 'LineWidth', 1.2);
    hold(ax1, 'off');
    xlim(ax1, [-win_sec, win_sec]);
    ylim(ax1, [0.5, nValid_time + 0.5]);
    xlabel(ax1, 'Time from breath peak (s)');
    ylabel(ax1, 'Breath # (sorted)');
    title(ax1, 'Time');
    axis(ax1, 'square');

    % --- Right: phase-normalized ---
    ax2 = subplot(1,2,2, 'Parent', fig);
    imagesc(ax2, phase_full, 1:nValid_ph, phase_sorted_dup);
    set(ax2, 'YDir', 'reverse');
    colormap(ax2, flipud(gray(256)));
    caxis(ax2, [0 cHigh]);
    hold(ax2, 'on');
    xline(ax2, 0, 'r--', 'LineWidth', 1.2);
    hold(ax2, 'off');
    xlim(ax2, [-2*pi, 2*pi]);
    ylim(ax2, [0.5, nValid_ph + 0.5]);
    set(ax2, 'XTick', [-2*pi -pi 0 pi 2*pi], ...
        'XTickLabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
    xlabel(ax2, 'Phase (rad)');
    set(ax2, 'YTickLabel', []);
    title(ax2, 'Phase-normalized');
    cb = colorbar(ax2);
    ylabel(cb, 'dF/F');
    axis(ax2, 'square');

    drawnow;

    % --- Save ---
    if saveFigs
        outDir = fullfile(folderPath, 'breath_trig_heatmap');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        outFig = fullfile(outDir, sprintf('breath_trig_heatmap_ROI%02d.fig', roi_id));
        savefig(fig, outFig);

        outPdf = fullfile(outDir, sprintf('breath_trig_heatmap_ROI%02d.pdf', roi_id));
        exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','none');

        fprintf("Saved ROI#%02d:\n  %s\n  %s\n", roi_id, outFig, outPdf);
    end
end

%% ============================== SUMMARY: ALL ROIs PHASE HEATMAP ==============================
avg_phase_dup = [avg_phase_all, avg_phase_all];  % duplicate to [-2π, 2π]

fig_sum = figure('Color', 'w', 'Name', 'breath_phase_all_ROIs');
set(fig_sum, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 10);

ax_sum = axes(fig_sum);
imagesc(ax_sum, phase_full, 1:nSpk, avg_phase_dup);
set(ax_sum, 'YDir', 'reverse');
colormap(ax_sum, flipud(gray(256)));
caxis(ax_sum, [0 prctile(avg_phase_dup(:), prcLim(2))]);
hold(ax_sum, 'on');
xline(ax_sum, 0, 'r--', 'LineWidth', 1.2);
hold(ax_sum, 'off');
xlim(ax_sum, [-2*pi, 2*pi]);
ylim(ax_sum, [0.5, nSpk + 0.5]);
set(ax_sum, 'XTick', [-2*pi -pi 0 pi 2*pi], ...
    'XTickLabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
set(ax_sum, 'YTick', 1:nSpk, 'YTickLabel', arrayfun(@(x) sprintf('#%d',x), roiList, 'Uni',false));
xlabel(ax_sum, 'Phase (rad)');
ylabel(ax_sum, 'ROI');
title(ax_sum, sprintf('%s  —  all ROIs mean phase dF/F', strrep(stem,'_','\_')), 'FontSize', 10);
cb_sum = colorbar(ax_sum);
ylabel(cb_sum, 'dF/F');

drawnow;

if saveFigs
    outDir = fullfile(folderPath, 'breath_trig_heatmap');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    outFig = fullfile(outDir, 'breath_phase_all_ROIs.fig');
    savefig(fig_sum, outFig);

    outPdf = fullfile(outDir, 'breath_phase_all_ROIs.pdf');
    exportgraphics(fig_sum, outPdf, 'ContentType','vector', 'BackgroundColor','none');

    fprintf("Saved summary:\n  %s\n  %s\n", outFig, outPdf);
end

fprintf("\nDone. Processed %d spiking ROIs + summary (%d total ROIs).\n", numel(roiList), N);
