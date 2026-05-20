%% breath_dt_scatter_260327.m
%  Per-ROI dt_last vs dt_next scatter with one-period duplication for
%  visual continuity across the cycle boundary.
%  Saves .fig + .pdf per ROI into folderPath/breath_dt_scatter/.
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_2000f_00001";
FPS        = 30;
saveFigs   = true;
% ──────────────────────────────────────────────────────────────────────────

%% ============================== LOAD ==============================
dffEntry = dir(fullfile(folderPath, '*_dFF.mat'));
assert(numel(dffEntry) == 1, "Expected 1 *_dFF.mat, found %d.", numel(dffEntry));
dffMatPath = fullfile(dffEntry.folder, dffEntry.name);
[~, stem]  = fileparts(dffMatPath);
D   = load(dffMatPath);
dFF = double(D.dFF);
[T, N] = size(dFF);
if isfield(D,'params') && isfield(D.params,'FPS')
    fs = double(D.params.FPS);
else
    fs = FPS;
    warning("params.FPS not found — using fallback FPS=%g", FPS);
end

spikeFile = fullfile(folderPath, 'ca_spike_data.mat');
assert(isfile(spikeFile), "ca_spike_data.mat not found.");
Ks      = load(spikeFile);
roiList = find(logical(Ks.ifSpike(:)'));
caData  = Ks.roi_spikes;
assert(~isempty(roiList), "No spiking ROIs found.");

bpFile = fullfile(folderPath, 'breath_peak_data.mat');
assert(isfile(bpFile), "breath_peak_data.mat not found.");
BP          = load(bpFile);
insp_onsets = sort(BP.insp_onset_idx(:));

b_frames = sort(insp_onsets(:));
Delta_k  = diff(b_frames) / fs;
T_mean   = mean(Delta_k);
T_std    = std(Delta_k);

cmap_roi = lines(numel(roiList));

fprintf("T=%d  N=%d  %d spiking ROIs  mean period=%.3f s\n", ...
    T, N, numel(roiList), T_mean);

%% ============================== PER-ROI SCATTER ==============================
for k = 1:numel(roiList)
    roi_id = roiList(k);
    col    = cmap_roi(k,:);
    ca_ev  = caData(roi_id).spike_idx(:);

    % dt_last / dt_next per spike
    dt_last = nan(numel(ca_ev), 1);
    dt_next = nan(numel(ca_ev), 1);
    for e = 1:numel(ca_ev)
        c   = ca_ev(e);
        bef = insp_onsets(insp_onsets <= c);
        aft = insp_onsets(insp_onsets >  c);
        if ~isempty(bef), dt_last(e) = (c - bef(end)) / fs; end
        if ~isempty(aft), dt_next(e) = (aft(1) - c)   / fs; end
    end
    keep = ~isnan(dt_last) & ~isnan(dt_next);
    xl = dt_last(keep);
    xn = dt_next(keep);
    nEv = numel(xl);

    % per-spike cycle length
    T_cyc = xl + xn;

    % --- FIGURE: two touching panels on the diagonal ---
    ax_max = T_mean;   % one mean breath cycle

    sq_cm  = 4.5;          % panel size (cm)
    ml = 1.8; mr = 0.4;    % margins
    mb = 1.0; mt = 1.2;
    fig_w = ml + 2*sq_cm + mr;
    fig_h = mb + 2*sq_cm + mt;

    fig = figure('Color', 'w', ...
        'Name', sprintf('breath_dt_scatter ROI#%02d', roi_id), ...
        'Units', 'centimeters', 'Position', [2 2 fig_w fig_h], ...
        'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
        'PaperPosition', [0 0 fig_w fig_h]);
    set(fig, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 10);
    sgtitle(fig, sprintf('%s  —  ROI #%d  (n=%d)', ...
        strrep(stem,'_','\_'), roi_id, nEv), 'FontSize', 10);

    % top-left panel:  [ml, mb+sq_cm, sq_cm, sq_cm]
    % bottom-right:    [ml+sq_cm, mb, sq_cm, sq_cm]
    positions = {[ml, mb + sq_cm, sq_cm, sq_cm], ...
                 [ml + sq_cm, mb, sq_cm, sq_cm]};

    for pp = 1:2
        ax = axes(fig, 'Units','centimeters', 'Position', positions{pp});
        hold(ax, 'on');

        % iso-period diagonals
        xv = linspace(0, ax_max, 100);
        for dv = T_mean + [-T_std, 0, T_std]
            plot(ax, xv, dv - xv, '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8);
        end

        % identity line
        plot(ax, [0 ax_max], [0 ax_max], '-', 'Color', [0.2 0.6 0.9], 'LineWidth', 1.2);

        % scatter
        scatter(ax, xl, xn, 12, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

        % diagonal crosshair
        mx = mean(xl);  my = mean(xn);
        s_diag = std((xl - xn) / sqrt(2));
        s_perp = std((xl + xn) / sqrt(2));
        d = [1 -1] / sqrt(2);
        p = [1  1] / sqrt(2);
        plot(ax, mx + s_diag*[-d(1) d(1)], my + s_diag*[-d(2) d(2)], ...
            '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 2);
        plot(ax, mx + s_perp*[-p(1) p(1)], my + s_perp*[-p(2) p(2)], ...
            '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 2);
        plot(ax, mx, my, 'o', 'Color', [0.85 0.33 0.1], ...
            'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerSize', 5);

        xlim(ax, [0, ax_max]);
        ylim(ax, [0, ax_max]);
        hold(ax, 'off');
        box(ax, 'on');

        if pp == 2
            xlabel(ax, 'dt_{last} (s)');
            set(ax, 'YTickLabel', []);
        else
            set(ax, 'XTickLabel', []);
            ylabel(ax, 'dt_{next} (s)');
        end
    end

    drawnow;

    % --- Save ---
    if saveFigs
        outDir = fullfile(folderPath, 'breath_dt_scatter');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        outFig = fullfile(outDir, sprintf('breath_dt_scatter_ROI%02d.fig', roi_id));
        savefig(fig, outFig);

        outPdf = fullfile(outDir, sprintf('breath_dt_scatter_ROI%02d.pdf', roi_id));
        exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','none');

        fprintf("Saved ROI#%02d:\n  %s\n  %s\n", roi_id, outFig, outPdf);
    end
end

fprintf("\nDone. Processed %d spiking ROIs.\n", numel(roiList));
