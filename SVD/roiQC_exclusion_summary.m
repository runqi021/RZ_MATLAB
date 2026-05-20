clear all; clc; close all;
%% roiQC_exclusion_summary.m
%  Aggregate exclusion data across all FOVs and generate CellPose tuning
%  suggestions + automatic QC threshold recommendations.
%
%  Reads roiQC_exclusion.mat files recursively under masterFolder.
%
%  Outputs (in masterFolder):
%    roiQC_master_summary.mat — aggregated data table
%    roiQC_master_summary.pdf — histograms + recommendations
%
%  RZ 2026-04

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
masterFolder = "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing";   % <-- set to top-level directory

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'mod')));
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% ======================================================================
%  LOAD ALL EXCLUSION FILES
%  ======================================================================
assert(masterFolder ~= "", 'Set masterFolder.');
exclHits = dir(fullfile(masterFolder, '**', 'roiQC_exclusion.mat'));
assert(~isempty(exclHits), 'No roiQC_exclusion.mat found under %s', masterFolder);

nFOV = numel(exclHits);
fprintf('Found %d FOV(s) with exclusion data.\n', nFOV);

% Aggregate all ROIs into flat arrays
all_fovID       = {};
all_roiLabel    = [];
all_p1_reject   = [];
all_p2_reject   = [];
all_final_keep  = [];
all_nPix        = [];
all_diam_px     = [];
all_diam_um     = [];
all_circ        = [];
all_peakSNR     = [];
all_skew        = [];
all_cp_diam     = [];

fov_summary = struct('fovID', {}, 'N_roi', {}, 'nKept', {}, 'nP1rej', {}, ...
    'nP2rej', {}, 'cp_diam', {});

for ff = 1:nFOV
    exclPath = fullfile(exclHits(ff).folder, exclHits(ff).name);
    E = load(exclPath);

    N = E.N_roi;
    fovID_str = string(E.fovID);

    all_fovID       = [all_fovID; repmat({char(fovID_str)}, N, 1)]; %#ok<AGROW>
    all_roiLabel    = [all_roiLabel;    E.roiLabels(:)];             %#ok<AGROW>
    all_p1_reject   = [all_p1_reject;   E.pass1_reject(:)];         %#ok<AGROW>
    all_final_keep  = [all_final_keep;  E.final_keep(:)];           %#ok<AGROW>
    all_nPix        = [all_nPix;        E.roi_nPix(:)];             %#ok<AGROW>
    all_diam_px     = [all_diam_px;     E.roi_diameter_px(:)];      %#ok<AGROW>
    all_diam_um     = [all_diam_um;     E.roi_diameter_um(:)];      %#ok<AGROW>
    all_circ        = [all_circ;        E.roi_circularity(:)];      %#ok<AGROW>
    all_peakSNR     = [all_peakSNR;     E.peakSNR(:)];              %#ok<AGROW>
    all_skew        = [all_skew;        E.skewness(:)];             %#ok<AGROW>

    if isfield(E, 'pass2_reject')
        all_p2_reject = [all_p2_reject; E.pass2_reject(:)]; %#ok<AGROW>
    else
        all_p2_reject = [all_p2_reject; false(N, 1)]; %#ok<AGROW>
    end

    cp_d = NaN;
    if isfield(E, 'cellpose_diameter') && isfinite(E.cellpose_diameter)
        cp_d = E.cellpose_diameter;
    end
    all_cp_diam = [all_cp_diam; repmat(cp_d, N, 1)]; %#ok<AGROW>

    nP1 = sum(E.pass1_reject);
    nP2 = 0;
    if isfield(E, 'pass2_reject'), nP2 = sum(E.pass2_reject); end
    nKept = sum(E.final_keep);

    fov_summary(ff).fovID   = char(fovID_str);
    fov_summary(ff).N_roi   = N;
    fov_summary(ff).nKept   = nKept;
    fov_summary(ff).nP1rej  = nP1;
    fov_summary(ff).nP2rej  = nP2;
    fov_summary(ff).cp_diam = cp_d;

    fprintf('  [%d] %-40s  ROIs=%d  kept=%d  P1rej=%d  P2rej=%d\n', ...
        ff, char(fovID_str), N, nKept, nP1, nP2);
end

N_total   = numel(all_final_keep);
N_kept    = sum(all_final_keep);
N_rejected = N_total - N_kept;

is_rejected = ~all_final_keep;

fprintf('\n========== AGGREGATE SUMMARY ==========\n');
fprintf('  Total ROIs:    %d across %d FOVs\n', N_total, nFOV);
fprintf('  Kept:          %d (%.1f%%)\n', N_kept, 100*N_kept/N_total);
fprintf('  Rejected:      %d (%.1f%%)\n', N_rejected, 100*N_rejected/N_total);

%% ======================================================================
%  CELLPOSE TUNING ANALYSIS
%  ======================================================================
fprintf('\n========== CELLPOSE TUNING SUGGESTIONS ==========\n');

% Size analysis
med_diam_kept = median(all_diam_um(~is_rejected), 'omitnan');
med_diam_rej  = median(all_diam_um(is_rejected),  'omitnan');
fprintf('  Diameter (um):  kept median=%.1f  |  rejected median=%.1f\n', ...
    med_diam_kept, med_diam_rej);

if isfinite(med_diam_rej) && med_diam_rej < med_diam_kept * 0.7
    fprintf('  >> Rejected ROIs are significantly SMALLER. Consider increasing CellPose diameter.\n');
elseif isfinite(med_diam_rej) && med_diam_rej > med_diam_kept * 1.3
    fprintf('  >> Rejected ROIs are significantly LARGER. Consider decreasing CellPose diameter.\n');
else
    fprintf('  >> Size distributions overlap. Diameter parameter appears reasonable.\n');
end

% Shape analysis
mean_circ_kept = mean(all_circ(~is_rejected), 'omitnan');
mean_circ_rej  = mean(all_circ(is_rejected),  'omitnan');
fprintf('  Circularity:    kept mean=%.3f  |  rejected mean=%.3f\n', ...
    mean_circ_kept, mean_circ_rej);

if isfinite(mean_circ_rej) && mean_circ_rej < mean_circ_kept * 0.7
    fprintf('  >> Rejected ROIs have LOW circularity (irregular shapes).\n');
    fprintf('     Consider lowering CellPose flow_threshold to reject irregular masks.\n');
end

% CellPose diameter used
cp_diams_used = unique(all_cp_diam(isfinite(all_cp_diam)));
if ~isempty(cp_diams_used)
    fprintf('  CellPose diameter(s) used: %s\n', mat2str(cp_diams_used'));
    med_kept_px = median(all_diam_px(~is_rejected), 'omitnan');
    fprintf('  Median kept ROI diameter:  %.1f px\n', med_kept_px);
    if abs(med_kept_px - cp_diams_used(1)) > 5
        fprintf('  >> Consider updating CellPose diameter to ~%.0f px to match actual cell sizes.\n', ...
            round(med_kept_px));
    end
end

%% ======================================================================
%  QC THRESHOLD DERIVATION
%  ======================================================================
fprintf('\n========== SUGGESTED QC THRESHOLDS ==========\n');

has_snr = any(isfinite(all_peakSNR));
has_skew = any(isfinite(all_skew));

if has_snr && N_rejected > 0
    % Find peakSNR threshold that best separates rejected from kept
    snr_kept = all_peakSNR(~is_rejected & isfinite(all_peakSNR));
    snr_rej  = all_peakSNR(is_rejected  & isfinite(all_peakSNR));

    if ~isempty(snr_rej) && ~isempty(snr_kept)
        thresholds = linspace(min(all_peakSNR(isfinite(all_peakSNR))), ...
                              max(all_peakSNR(isfinite(all_peakSNR))), 200);
        best_acc = 0;
        best_thr = thresholds(1);
        for ti = 1:numel(thresholds)
            thr = thresholds(ti);
            tp = sum(snr_rej < thr);     % correctly flagged rejected
            tn = sum(snr_kept >= thr);    % correctly passed kept
            acc = (tp + tn) / (numel(snr_rej) + numel(snr_kept));
            if acc > best_acc
                best_acc = acc;
                best_thr = thr;
            end
        end
        fprintf('  peakSNR threshold:  %.1f  (accuracy %.1f%%)\n', best_thr, 100*best_acc);
        fprintf('    Kept median=%.1f  |  Rejected median=%.1f\n', ...
            median(snr_kept), median(snr_rej));
    end
end

if has_skew && N_rejected > 0
    skew_kept = all_skew(~is_rejected & isfinite(all_skew));
    skew_rej  = all_skew(is_rejected  & isfinite(all_skew));

    if ~isempty(skew_rej) && ~isempty(skew_kept)
        thresholds = linspace(min(all_skew(isfinite(all_skew))), ...
                              max(all_skew(isfinite(all_skew))), 200);
        best_acc = 0;
        best_thr = thresholds(1);
        for ti = 1:numel(thresholds)
            thr = thresholds(ti);
            tp = sum(skew_rej < thr);
            tn = sum(skew_kept >= thr);
            acc = (tp + tn) / (numel(skew_rej) + numel(skew_kept));
            if acc > best_acc
                best_acc = acc;
                best_thr = thr;
            end
        end
        fprintf('  skewness threshold: %.2f  (accuracy %.1f%%)\n', best_thr, 100*best_acc);
        fprintf('    Kept median=%.2f  |  Rejected median=%.2f\n', ...
            median(skew_kept), median(skew_rej));
    end
end

% Diameter threshold
min_diam = NaN;
diam_kept = all_diam_um(~is_rejected);
diam_rej  = all_diam_um(is_rejected);
if ~isempty(diam_rej) && ~isempty(diam_kept)
    p5_kept  = prctile(diam_kept, 5);
    min_diam = max(p5_kept * 0.5, median(diam_rej));
    fprintf('  min diameter (um): %.1f  (5th pctile of kept = %.1f)\n', min_diam, p5_kept);
end

% Circularity threshold
p5_circ = NaN;
circ_kept = all_circ(~is_rejected);
circ_rej  = all_circ(is_rejected);
if ~isempty(circ_rej) && ~isempty(circ_kept)
    p5_circ = prctile(circ_kept, 5);
    fprintf('  min circularity:   %.2f  (5th pctile of kept = %.2f)\n', ...
        p5_circ, p5_circ);
end

%% ======================================================================
%  FIGURES
%  ======================================================================
tmpPDFs = {};

% ---- Figure 1: Morphometric distributions (rejected vs kept) ----
fig1 = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 14 8]);
t1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t1, sprintf('ROI QC Summary — %d FOVs, %d ROIs (%d rejected)', ...
    nFOV, N_total, N_rejected), 'FontSize', 14);

% Diameter histogram
ax1 = nexttile(t1);
hold(ax1, 'on');
edges = linspace(0, max(all_diam_um)*1.05, 30);
histogram(ax1, all_diam_um(~is_rejected), edges, 'FaceColor', [0.3 0.7 0.3], ...
    'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Kept');
histogram(ax1, all_diam_um(is_rejected), edges, 'FaceColor', [0.8 0.2 0.2], ...
    'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Rejected');
xlabel(ax1, 'Diameter (\mum)');
ylabel(ax1, 'Count');
title(ax1, 'ROI Diameter');
legend(ax1, 'Location', 'best');
box(ax1, 'off');

% Circularity histogram
ax2 = nexttile(t1);
hold(ax2, 'on');
edges_c = linspace(0, 1, 25);
histogram(ax2, all_circ(~is_rejected), edges_c, 'FaceColor', [0.3 0.7 0.3], ...
    'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Kept');
histogram(ax2, all_circ(is_rejected), edges_c, 'FaceColor', [0.8 0.2 0.2], ...
    'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Rejected');
xlabel(ax2, 'Circularity');
ylabel(ax2, 'Count');
title(ax2, 'ROI Circularity');
legend(ax2, 'Location', 'best');
box(ax2, 'off');

% Peak SNR histogram (if available)
ax3 = nexttile(t1);
if has_snr
    hold(ax3, 'on');
    snr_all = all_peakSNR(isfinite(all_peakSNR));
    edges_s = linspace(0, min(prctile(snr_all, 99), max(snr_all)), 30);
    histogram(ax3, all_peakSNR(~is_rejected & isfinite(all_peakSNR)), edges_s, ...
        'FaceColor', [0.3 0.7 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Kept');
    histogram(ax3, all_peakSNR(is_rejected & isfinite(all_peakSNR)), edges_s, ...
        'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Rejected');
    xlabel(ax3, 'Peak SNR');
    ylabel(ax3, 'Count');
    title(ax3, 'Peak SNR');
    legend(ax3, 'Location', 'best');
    box(ax3, 'off');
else
    text(ax3, 0.5, 0.5, 'No SNR data', 'HorizontalAlignment', 'center', 'FontSize', 12);
    axis(ax3, 'off');
end

% Skewness histogram (if available)
ax4 = nexttile(t1);
if has_skew
    hold(ax4, 'on');
    sk_all = all_skew(isfinite(all_skew));
    edges_k = linspace(min(sk_all), max(sk_all), 30);
    histogram(ax4, all_skew(~is_rejected & isfinite(all_skew)), edges_k, ...
        'FaceColor', [0.3 0.7 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Kept');
    histogram(ax4, all_skew(is_rejected & isfinite(all_skew)), edges_k, ...
        'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'DisplayName', 'Rejected');
    xlabel(ax4, 'Skewness');
    ylabel(ax4, 'Count');
    title(ax4, 'dFF Skewness');
    legend(ax4, 'Location', 'best');
    box(ax4, 'off');
else
    text(ax4, 0.5, 0.5, 'No skewness data', 'HorizontalAlignment', 'center', 'FontSize', 12);
    axis(ax4, 'off');
end

tmp1 = fullfile(tempdir, 'roiQC_summary_fig1.pdf');
exportgraphics(fig1, tmp1, 'ContentType', 'vector');
tmpPDFs{end+1} = tmp1;
close(fig1);

% ---- Figure 2: Per-FOV summary bar chart ----
fig2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 14 6]);
ax5 = axes(fig2);

fov_names = {fov_summary.fovID};
nROIs  = [fov_summary.N_roi];
nKept  = [fov_summary.nKept];
nP1rej = [fov_summary.nP1rej];
nP2rej = [fov_summary.nP2rej];

barData = [nKept; nP1rej; nP2rej]';
b = bar(ax5, barData, 'stacked');
b(1).FaceColor = [0.3 0.7 0.3];
b(2).FaceColor = [0.9 0.5 0.2];
b(3).FaceColor = [0.8 0.2 0.2];

set(ax5, 'XTick', 1:numel(fov_names), 'XTickLabel', fov_names, 'FontSize', 7);
xtickangle(ax5, 45);
ylabel(ax5, 'ROI count');
title(ax5, 'Per-FOV ROI Exclusion', 'FontSize', 12);
legend(ax5, {'Kept', 'Pass 1 reject', 'Pass 2 reject'}, 'Location', 'best');
box(ax5, 'off');

tmp2 = fullfile(tempdir, 'roiQC_summary_fig2.pdf');
exportgraphics(fig2, tmp2, 'ContentType', 'vector');
tmpPDFs{end+1} = tmp2;
close(fig2);

%% ======================================================================
%  SAVE OUTPUTS
%  ======================================================================
% Save aggregated data
summary = struct();
summary.nFOV          = nFOV;
summary.N_total       = N_total;
summary.N_kept        = N_kept;
summary.N_rejected    = N_rejected;
summary.fov_summary   = fov_summary;
summary.all_fovID     = all_fovID;
summary.all_roiLabel  = all_roiLabel;
summary.all_final_keep = all_final_keep;
summary.all_p1_reject  = all_p1_reject;
summary.all_p2_reject  = all_p2_reject;
summary.all_diam_um    = all_diam_um;
summary.all_diam_px    = all_diam_px;
summary.all_circ       = all_circ;
summary.all_peakSNR    = all_peakSNR;
summary.all_skew       = all_skew;

outMat = fullfile(masterFolder, 'roiQC_master_summary.mat');
save(outMat, '-struct', 'summary', '-v7.3');
fprintf('\nSaved: %s\n', outMat);

% Merge PDFs
if ~isempty(tmpPDFs)
    outPDF = fullfile(masterFolder, 'roiQC_master_summary.pdf');
    if isfile(outPDF), delete(outPDF); end
    append_pdfs_RZ(outPDF, tmpPDFs{:});
    fprintf('Saved: %s\n', outPDF);

    for ii = 1:numel(tmpPDFs)
        if isfile(tmpPDFs{ii}), delete(tmpPDFs{ii}); end
    end
end

% Save suggested thresholds
thresholds = struct();
if exist('best_thr', 'var')
    thresholds.peakSNR_threshold = best_thr;
end
thresholds.min_diameter_um   = min_diam;
thresholds.min_circularity   = p5_circ;
thresholds.notes = 'Derived from human QC decisions. Use as starting points, not absolute rules.';

thrMat = fullfile(masterFolder, 'roiQC_suggested_thresholds.mat');
save(thrMat, '-struct', 'thresholds');
fprintf('Saved: %s\n', thrMat);

fprintf('\n==================== SUMMARY COMPLETE ====================\n');
