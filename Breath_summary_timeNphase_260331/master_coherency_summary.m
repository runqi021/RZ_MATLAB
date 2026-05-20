% master_coherency_summary.m
%  Collects coherence_data.mat from multiple experiment folders
%  (output of breath_combined_summary_260329.m) and plots a single
%  master polar plot — spike-train coherence, color-coded by experiment.
%
%  Usage: edit expFolders below, then run.

clear; close all; clc;

%% ========================= USER PARAMETERS =============================
expFolders = {
    "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing"
    "D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing"
    "C:\Users\Admin\Desktop\260330_sst_soma_g8s\phys\maybe_breathing"
    "D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed"
    "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing"
};

doSave = true;
outDir = "C:\Users\Admin\Desktop\RZ_MATLAB\Breath_summary_timeNphase_260331\master_summary";
% ========================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% 1 -- Collect coherence_data from each experiment
nExp = numel(expFolders);
all_th   = [];
all_r    = [];
all_rlo  = [];
all_rhi  = [];
all_dphi = [];
all_sig  = [];
all_exp  = [];   % experiment index per ROI
all_labels = {};
all_confC = [];
expNames = cell(nExp, 1);

for ee = 1:nExp
    matFile = fullfile(expFolders{ee}, 'breath_combined_summary', 'coherence_data.mat');
    if ~isfile(matFile)
        fprintf('WARNING: coherence_data.mat not found in %s, skipping.\n', expFolders{ee});
        continue;
    end

    cd = load(matFile);
    cd = cd.coherence_data;

    n = numel(cd.all_th_k_spk);
    all_th   = [all_th;   cd.all_th_k_spk];    %#ok<AGROW>
    all_r    = [all_r;    cd.all_r_k_spk];      %#ok<AGROW>
    all_rlo  = [all_rlo;  cd.all_rlo_spk];      %#ok<AGROW>
    all_rhi  = [all_rhi;  cd.all_rhi_spk];      %#ok<AGROW>
    all_dphi = [all_dphi; cd.all_dphi_spk];     %#ok<AGROW>
    all_sig  = [all_sig;  cd.all_is_sig_spk];   %#ok<AGROW>
    all_exp  = [all_exp;  repmat(ee, n, 1)];    %#ok<AGROW>
    all_labels = [all_labels; cd.labels_all];    %#ok<AGROW>
    all_confC = [all_confC; cd.confC_mean];      %#ok<AGROW>

    [~, folderName] = fileparts(expFolders{ee});
    if isempty(folderName)
        [~, folderName] = fileparts(fileparts(expFolders{ee}));
    end
    expNames{ee} = folderName;

    fprintf('Loaded %s: %d ROIs (%d sig)\n', folderName, n, sum(cd.all_is_sig_spk));
end

nTotal = numel(all_th);
nSig   = sum(all_sig);
confC_grand = mean(all_confC, 'omitnan');

fprintf('\n=== Master summary: %d ROIs across %d experiments (%d sig) ===\n', ...
    nTotal, nExp, nSig);

%% 2 -- Color map: one color per experiment
cmap_exp = lines(nExp);

%% 3 -- Master polar plot (spike-train coherence)
fig_w = 14; fig_h = 12;
fig = figure('Color', 'w', ...
    'Name', 'Master Coherency Summary (spike)', ...
    'Units', 'centimeters', 'Position', [2 2 fig_w fig_h], ...
    'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
    'PaperPosition', [0 0 fig_w fig_h]);

ax = polaraxes(fig, 'Position', [0.08 0.08 0.65 0.82]);
hold(ax, 'on');

% Confidence circle
th_circ = linspace(0, 2*pi, 360);
polarplot(ax, th_circ, repmat(confC_grand, 1, 360), 'k--', 'LineWidth', 1.2);

% Plot each ROI
for rr = 1:nTotal
    ee = all_exp(rr);
    if all_sig(rr)
        fc = cmap_exp(ee, :);
        ms = 4;
    else
        fc = [0.7 0.7 0.7];
        ms = 2.5;
    end
    polarplot(ax, all_th(rr), all_r(rr), 'o', ...
        'MarkerFaceColor', fc, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', ms, 'LineWidth', 0.5);
end

ax.RLim              = [0, 1];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir          = 'clockwise';
ax.RAxisLocation     = 180;
ax.FontSize          = 8;
hold(ax, 'off');
title(ax, sprintf('Spike coherence (%d/%d sig, confC=%.3f)', nSig, nTotal, confC_grand), ...
    'Interpreter', 'none', 'FontSize', 9);

% Legend (experiment colors)
ax_leg = axes(fig, 'Units', 'normalized', 'Position', [0.78 0.3 0.2 0.5], 'Visible', 'off');
hold(ax_leg, 'on');
for ee = 1:nExp
    plot(ax_leg, NaN, NaN, 'o', 'MarkerFaceColor', cmap_exp(ee,:), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6);
end
plot(ax_leg, NaN, NaN, 'o', 'MarkerFaceColor', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 4);
legend(ax_leg, [expNames; {'n.s.'}], ...
    'Interpreter', 'none', 'FontSize', 6, 'Location', 'west', 'Box', 'off');
hold(ax_leg, 'off');

drawnow;

%% 4 -- Save
if doSave
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    pdfPath = fullfile(outDir, 'master_coherency_summary.pdf');
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('Saved: %s\n', pdfPath);
end
