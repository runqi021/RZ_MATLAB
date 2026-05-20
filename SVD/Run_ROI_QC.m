clear all; clc; close all;
%% Run_ROI_QC.m
%  Per-ROI quality control for 2-photon calcium imaging.
%
%  Tier 1: Trace-based SNR (from F_roi / dFF .mat files)
%  Tier 2: Motion dipole via per-ROI local SVD (from MC TIFF patches)
%
%  Outputs (in folderPath):
%    roiQC_output.mat  — all per-ROI scores, keep/reject flags, params
%    roiQC.pdf         — 3-page diagnostic report
%
%  RZ 2026-04-01

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
folderPath           = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_2000f_00001";      % <-- paste data subfolder
patchSize_um         = 40;      % local SVD patch size (microns), square
nLocalModes          = 30;      % truncated SVD modes per ROI
skewnessThreshold    = 0.3;     % dFF skewness below this -> suspect
peakSNR_Threshold    = 5;       % peak SNR below this -> suspect
motionScoreThreshold = 0.5;     % motion_score above this -> suspect
minROIPixels         = 10;      % minimum pixels for reliable plane fit
nTopROIs             = 5;       % worst ROIs to show in detail figure
savePDF              = true;

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'mod')));

%% ======================================================================
%  VALIDATE FOLDER
%  ======================================================================
assert(strlength(folderPath) > 0, 'Set folderPath at the top of the script.');
assert(isfolder(folderPath), 'Folder not found: %s', folderPath);

%% ======================================================================
%  DETECT FPS + PIXEL SIZE
%  ======================================================================
[fps, scan_meta] = detect_session_fps(folderPath);
fprintf('[ROI QC] Detected FPS = %d Hz\n', fps);

assert(isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um), ...
    'detect_session_fps did not return pixelSize_um. Need TIFF metadata.');
pixelSize_um = scan_meta.pixelSize_um;
fprintf('[ROI QC] Pixel size = %.4f um/px\n', pixelSize_um);

patchRadius_px = round((patchSize_um / 2) / pixelSize_um);
fprintf('[ROI QC] Patch radius = %d px  (patch = %d x %d px)\n', ...
    patchRadius_px, 2*patchRadius_px+1, 2*patchRadius_px+1);

%% ======================================================================
%  LOAD dFF DATA
%  ======================================================================
dff_hits = dir(fullfile(folderPath, '*_dFF.mat'));
assert(~isempty(dff_hits), 'No *_dFF.mat found in %s', folderPath);
dff_path = fullfile(dff_hits(1).folder, dff_hits(1).name);
fprintf('[ROI QC] Loading dFF:  %s\n', dff_path);

dffData = load(dff_path);
assert(isfield(dffData, 'F_roi'), '_dFF.mat missing F_roi.');
assert(isfield(dffData, 'dFF'),   '_dFF.mat missing dFF.');
assert(isfield(dffData, 'params'),'_dFF.mat missing params.');

F_roi = double(dffData.F_roi);   % [T_dff x N]
dFF   = double(dffData.dFF);     % [T_dff x N]
[T_dff, N_roi] = size(F_roi);

assert(isequal(size(dFF), [T_dff, N_roi]), ...
    'dFF size [%d x %d] does not match F_roi [%d x %d].', ...
    size(dFF,1), size(dFF,2), T_dff, N_roi);

tossFrames = 0;
if isfield(dffData.params, 'tossFrames')
    tossFrames = dffData.params.tossFrames;
end
fprintf('[ROI QC] F_roi: [%d x %d],  tossFrames = %d\n', T_dff, N_roi, tossFrames);

%% ======================================================================
%  LOAD cpSAM (maskL + avgPath)
%  ======================================================================
sam_hits = dir(fullfile(folderPath, '*_cpSAM_output.mat'));
assert(~isempty(sam_hits), 'No *_cpSAM_output.mat found in %s', folderPath);
sam_path = fullfile(sam_hits(1).folder, sam_hits(1).name);
fprintf('[ROI QC] Loading SAM:  %s\n', sam_path);

samData = load(sam_path);
assert(isfield(samData, 'maskL'), 'SAM file missing maskL.');
maskL = samData.maskL;   % [H x W] uint16 label image

% Average projection path for background display
avgPath = '';
if isfield(samData, 'avgPath') && (ischar(samData.avgPath) || isstring(samData.avgPath))
    avgPath = char(samData.avgPath);
end

%% ======================================================================
%  DISCOVER MC TIFF
%  ======================================================================
mc_hits = dir(fullfile(folderPath, '*_MC*.tif'));
if isempty(mc_hits)
    error('No *_MC*.tif found in %s', folderPath);
end

% Filter out AVG, QC, shifts TIFFs
names = {mc_hits.name};
keep_idx = true(size(names));
for ii = 1:numel(names)
    if contains(names{ii}, '_AVG', 'IgnoreCase', true) || ...
       contains(names{ii}, '_QC',  'IgnoreCase', true) || ...
       contains(names{ii}, '_shifts', 'IgnoreCase', true)
        keep_idx(ii) = false;
    end
end
mc_hits = mc_hits(keep_idx);
assert(~isempty(mc_hits), 'No valid MC TIFF after filtering AVG/QC/shifts in %s', folderPath);

% Prefer the TIFF with the most _MC occurrences (2-pass -> *_MC_MC.tif)
nMC = zeros(numel(mc_hits), 1);
for ii = 1:numel(mc_hits)
    nMC(ii) = numel(strfind(upper(mc_hits(ii).name), '_MC'));
end
[~, bestIdx] = max(nMC);
tifPath = fullfile(mc_hits(bestIdx).folder, mc_hits(bestIdx).name);
fprintf('[ROI QC] Using MC TIFF: %s\n', tifPath);

%% ======================================================================
%  LOAD MC TIFF
%  ======================================================================
fprintf('[ROI QC] Loading TIFF stack...\n');
tic;
tifInfo = imfinfo(tifPath);
H = tifInfo(1).Height;
W = tifInfo(1).Width;
T_tiff = numel(tifInfo);

Y_full = zeros(H, W, T_tiff, 'single');
for kk = 1:T_tiff
    Y_full(:,:,kk) = single(imread(tifPath, kk, 'Info', tifInfo));
end
fprintf('[ROI QC] Loaded %d frames [%d x %d] in %.1f s\n', T_tiff, H, W, toc);

%% ======================================================================
%  FRAME ALIGNMENT (TIFF vs dFF)
%  ======================================================================
if T_tiff > T_dff && tossFrames > 0
    fprintf('[ROI QC] TIFF has %d frames, dFF has %d. Skipping first %d toss frames from TIFF.\n', ...
        T_tiff, T_dff, tossFrames);
    Y_full = Y_full(:, :, tossFrames+1:end);
    T_tiff_aligned = size(Y_full, 3);
    assert(T_tiff_aligned == T_dff, ...
        'After toss, TIFF frames (%d) still != dFF frames (%d).', T_tiff_aligned, T_dff);
elseif T_tiff == T_dff
    fprintf('[ROI QC] TIFF and dFF frame counts match (%d).\n', T_dff);
else
    error('TIFF frames (%d) != dFF frames (%d) and tossFrames=%d does not resolve mismatch.', ...
        T_tiff, T_dff, tossFrames);
end
T = T_dff;

%% ======================================================================
%  TIER 1: TRACE-BASED SNR
%  ======================================================================
fprintf('[ROI QC] Computing Tier 1: Trace SNR...\n');

noise   = zeros(N_roi, 1);
skew    = zeros(N_roi, 1);
peakSNR = zeros(N_roi, 1);

% Detrend dFF before computing skewness and peakSNR
dFF_dt = detrend(dFF);

for i = 1:N_roi
    frame_diffs = diff(F_roi(:, i));
    noise(i) = 1.4826 * median(abs(frame_diffs)) / sqrt(2);

    skew(i) = skewness(dFF_dt(:, i));

    mean_F = mean(F_roi(:, i));
    if noise(i) > 0 && mean_F > 0
        peakSNR(i) = max(dFF_dt(:, i)) / (noise(i) / mean_F);
    else
        peakSNR(i) = 0;
    end
end

fprintf('[ROI QC] Tier 1 done. Mean skewness=%.3f, mean peakSNR=%.2f\n', ...
    mean(skew), mean(peakSNR));

%% ======================================================================
%  TIER 2: PER-ROI LOCAL SVD + DIPOLE METRIC
%  ======================================================================
fprintf('[ROI QC] Computing Tier 2: Per-ROI Local SVD...\n');

% Preallocate
motion_score  = NaN(N_roi, 1);
dipole_frac   = NaN(N_roi, nLocalModes);
dipole_angle  = NaN(N_roi, nLocalModes);
dom_dipole_angle = NaN(N_roi, 1);

% ROI centroids for patch extraction
props = regionprops(maskL, 'Centroid');

% Store local SVD results for Figure 3
Uloc_cell   = cell(N_roi, 1);
Sloc_cell   = cell(N_roi, 1);
patch_size_cell = cell(N_roi, 1);   % [nRows, nCols] of each patch
mask_local_cell = cell(N_roi, 1);

for i = 1:N_roi
    fprintf('[ROI QC] ROI %d/%d ...\n', i, N_roi);

    % Centroid (regionprops returns [col, row])
    cx = props(i).Centroid(1);
    cy = props(i).Centroid(2);

    % Patch bounds (clipped to image)
    r1 = max(1,  round(cy) - patchRadius_px);
    r2 = min(H,  round(cy) + patchRadius_px);
    c1 = max(1,  round(cx) - patchRadius_px);
    c2 = min(W,  round(cx) + patchRadius_px);

    nRows_patch = r2 - r1 + 1;
    nCols_patch = c2 - c1 + 1;
    nPix_patch  = nRows_patch * nCols_patch;

    % Local ROI mask within patch
    mask_local = (maskL(r1:r2, c1:c2) == i);
    nPix_roi = sum(mask_local(:));

    patch_size_cell{i} = [nRows_patch, nCols_patch];
    mask_local_cell{i} = mask_local;

    % Skip if too few ROI pixels
    if nPix_roi < minROIPixels
        fprintf('  -> Skipping ROI %d: only %d pixels (min=%d)\n', i, nPix_roi, minROIPixels);
        continue
    end

    % Extract patch and reshape to [nPix_patch x T]
    Y_patch = reshape(Y_full(r1:r2, c1:c2, :), [nPix_patch, T]);
    Y_patch = single(Y_patch);

    % Demean per pixel
    Y_patch = Y_patch - mean(Y_patch, 2);

    % Truncated SVD
    nModes_use = min(nLocalModes, min(nPix_patch, T) - 1);
    if nModes_use < 1
        continue
    end
    [U_loc, S_loc, ~] = svds(double(Y_patch), nModes_use);

    % Store for Figure 3
    Uloc_cell{i} = U_loc;
    Sloc_cell{i} = S_loc;

    % ROI pixel coordinates (centered)
    [rr, cc] = find(mask_local);
    x_c = cc - mean(cc);
    y_c = rr - mean(rr);

    % Design matrix for plane fit
    A = [x_c(:), y_c(:), ones(nPix_roi, 1)];

    % Per-mode dipole analysis
    s_loc = diag(S_loc);
    varExp_loc = s_loc.^2 / sum(s_loc.^2);

    for k = 1:nModes_use
        % Reshape eigenface to patch grid
        eigenface = reshape(U_loc(:, k), [nRows_patch, nCols_patch]);

        % Extract values at ROI pixels
        u_ik = eigenface(mask_local);

        % Plane fit
        coeffs = A \ u_ik(:);

        % Gradient component
        grad_vals = A(:, 1:2) * coeffs(1:2);

        dipole_frac(i, k)  = var(grad_vals) / max(var(u_ik(:)), eps);
        dipole_angle(i, k) = atan2(coeffs(2), coeffs(1));
    end

    % Weighted motion score
    motion_score(i) = sum(varExp_loc .* dipole_frac(i, 1:nModes_use)') / sum(varExp_loc);

    % Dominant dipole direction from highest-weighted mode
    weighted_dipole = varExp_loc .* dipole_frac(i, 1:nModes_use)';
    [~, dom_idx] = max(weighted_dipole);
    dom_dipole_angle(i) = dipole_angle(i, dom_idx);
end

fprintf('[ROI QC] Tier 2 done. Mean motion_score=%.4f\n', mean(motion_score, 'omitnan'));

%% ======================================================================
%  KEEP / REJECT DECISION
%  ======================================================================
reject_lowSNR   = peakSNR < peakSNR_Threshold;
reject_lowSkew  = skew < skewnessThreshold;
reject_motion   = motion_score > motionScoreThreshold;
reject_tooSmall = roi_nPix < minROIPixels;
reject_any      = reject_lowSNR | reject_lowSkew | reject_motion | reject_tooSmall;
keep            = ~reject_any;

%% ======================================================================
%  FIGURES
%  ======================================================================
fprintf('[ROI QC] Generating figures...\n');
tmpPDFs = {};

% ---- Figure 1: SNR Overview (2x2) ----
fig1 = figure('Visible','off', 'Color','w', 'Units','inches', 'Position',[1 1 12 9]);
t1 = tiledlayout(fig1, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(t1, 'ROI QC -- Tier 1: Trace SNR', 'Interpreter','none', 'FontSize',14);

% (1,1) Noise floor per ROI
ax1 = nexttile(t1);
[noise_sorted, noise_order] = sort(noise, 'descend');
bar(ax1, noise_sorted, 'FaceColor',[0.4 0.6 0.8], 'EdgeColor','none');
set(ax1, 'XTick', 1:N_roi, 'XTickLabel', noise_order, 'FontSize', 7);
xtickangle(ax1, 90);
xlabel(ax1, 'ROI #', 'Interpreter','none');
ylabel(ax1, 'Noise (a.u.)', 'Interpreter','none');
title(ax1, 'Noise floor per ROI (MAD estimator)', 'Interpreter','none');
xlim(ax1, [0.5, N_roi+0.5]);
set(ax1, 'Box','off');

% (1,2) Skewness per ROI
ax2 = nexttile(t1);
[skew_sorted, skew_order] = sort(skew, 'descend');
bar(ax2, skew_sorted, 'FaceColor',[0.5 0.7 0.4], 'EdgeColor','none');
hold(ax2, 'on');
yline(ax2, skewnessThreshold, '--r', 'LineWidth', 1.5);
set(ax2, 'XTick', 1:N_roi, 'XTickLabel', skew_order, 'FontSize', 7);
xtickangle(ax2, 90);
xlabel(ax2, 'ROI #', 'Interpreter','none');
ylabel(ax2, 'Skewness', 'Interpreter','none');
title(ax2, sprintf('dFF Skewness (threshold = %.2f)', skewnessThreshold), 'Interpreter','none');
xlim(ax2, [0.5, N_roi+0.5]);
set(ax2, 'Box','off');

% (2,1) Peak SNR per ROI
ax3 = nexttile(t1);
[peakSNR_sorted, snr_order] = sort(peakSNR, 'descend');
bar(ax3, peakSNR_sorted, 'FaceColor',[0.8 0.5 0.3], 'EdgeColor','none');
hold(ax3, 'on');
yline(ax3, peakSNR_Threshold, '--r', 'LineWidth', 1.5);
set(ax3, 'XTick', 1:N_roi, 'XTickLabel', snr_order, 'FontSize', 7);
xtickangle(ax3, 90);
xlabel(ax3, 'ROI #', 'Interpreter','none');
ylabel(ax3, 'Peak SNR', 'Interpreter','none');
title(ax3, sprintf('Peak SNR (threshold = %d)', peakSNR_Threshold), 'Interpreter','none');
xlim(ax3, [0.5, N_roi+0.5]);
set(ax3, 'Box','off');

% (2,2) Scatter: skewness vs peakSNR, colored by motion_score
ax4 = nexttile(t1);
ms_for_color = motion_score;
ms_for_color(isnan(ms_for_color)) = 0;
scatter(ax4, skew, peakSNR, 30, ms_for_color, 'filled', 'MarkerEdgeColor','none');
hold(ax4, 'on');
xline(ax4, skewnessThreshold, '--r', 'LineWidth', 1.2);
yline(ax4, peakSNR_Threshold, '--r', 'LineWidth', 1.2);
colormap(ax4, parula(256));
caxis(ax4, [0, max(max(ms_for_color), eps)]);
cb = colorbar(ax4);
cb.Label.String = 'motion_score';
cb.Label.Interpreter = 'none';
xlabel(ax4, 'Skewness', 'Interpreter','none');
ylabel(ax4, 'Peak SNR', 'Interpreter','none');
title(ax4, 'Skewness vs Peak SNR (color = motion score)', 'Interpreter','none');
set(ax4, 'Box','off');

if savePDF
    tmp1 = fullfile(tempdir, 'roiQC_fig1.pdf');
    exportgraphics(fig1, tmp1, 'ContentType','vector');
    tmpPDFs{end+1} = tmp1;
end
close(fig1);

% ---- Figure 2: Motion Dipole (1x2) ----
fig2 = figure('Visible','off', 'Color','w', 'Units','inches', 'Position',[1 1 12 5.5]);
t2 = tiledlayout(fig2, 1, 2, 'TileSpacing','compact', 'Padding','compact');
title(t2, 'ROI QC -- Tier 2: Motion Dipole', 'Interpreter','none', 'FontSize',14);

% Left: motion_score bar chart (sorted descending)
ax5 = nexttile(t2);
ms_plot = motion_score;
ms_plot(isnan(ms_plot)) = 0;
[ms_sorted, ms_order] = sort(ms_plot, 'descend');
bar(ax5, ms_sorted, 'FaceColor',[0.7 0.3 0.3], 'EdgeColor','none');
hold(ax5, 'on');
yline(ax5, motionScoreThreshold, '--r', 'LineWidth', 1.5);
set(ax5, 'XTick', 1:N_roi, 'XTickLabel', ms_order, 'FontSize', 7);
xtickangle(ax5, 90);
xlabel(ax5, 'ROI #', 'Interpreter','none');
ylabel(ax5, 'motion\_score', 'Interpreter','none');
title(ax5, sprintf('Motion Score (threshold = %.2f)', motionScoreThreshold), 'Interpreter','none');
xlim(ax5, [0.5, N_roi+0.5]);
set(ax5, 'Box','off');

% Right: Spatial map with avg image background
ax6 = nexttile(t2);
if ~isempty(avgPath) && isfile(avgPath)
    bgImg = mat2gray(double(imread(avgPath)));
else
    bgImg = zeros(size(maskL));
end
imagesc(ax6, bgImg); colormap(ax6, gray(256)); hold(ax6, 'on');
axis(ax6, 'image'); axis(ax6, 'off');
overlay_roi_values(ax6, maskL, ms_plot, 'hot', [0, max(max(ms_plot), eps)]);
cb2 = colorbar(ax6);
cb2.Label.String = 'motion\_score';
cb2.Label.Interpreter = 'none';
title(ax6, 'Motion Score Spatial Map', 'Interpreter','none');

if savePDF
    tmp2 = fullfile(tempdir, 'roiQC_fig2.pdf');
    exportgraphics(fig2, tmp2, 'ContentType','vector');
    tmpPDFs{end+1} = tmp2;
end
close(fig2);

% ---- Figure 3+: Eigenface Detail — 1 ROI per page, 30 modes (5x6) ----
[~, allSortedIdx] = sort(motion_score, 'descend', 'MissingPlacement','last');
nModesShow = min(nLocalModes, 30);
nRows_grid = 5;
nCols_grid = 6;

for pp = 1:N_roi
    roi_id = allSortedIdx(pp);
    U_loc  = Uloc_cell{roi_id};
    pSz    = patch_size_cell{roi_id};
    mLocal = mask_local_cell{roi_id};

    figR = figure('Visible','off', 'Color','w', 'Units','inches', ...
        'Position',[1 1 2.5*nCols_grid, 2.5*nRows_grid]);
    tR = tiledlayout(figR, nRows_grid, nCols_grid, 'TileSpacing','tight', 'Padding','compact');

    ms_val = motion_score(roi_id);
    if isnan(ms_val), ms_str = 'NaN'; else, ms_str = sprintf('%.3f', ms_val); end
    title(tR, sprintf('ROI #%d  |  motion\\_score = %s  |  skew = %.2f  |  peakSNR = %.1f', ...
        roi_id, ms_str, skew(roi_id), peakSNR(roi_id)), ...
        'Interpreter','none', 'FontSize',12);

    for k = 1:nRows_grid * nCols_grid
        axE = nexttile(tR);

        if isempty(U_loc) || k > size(U_loc, 2)
            axis(axE, 'off');
            continue
        end

        eigenface = reshape(U_loc(:, k), pSz);
        imagesc(axE, eigenface);
        hold(axE, 'on');

        % Overlay ROI boundary
        bnd = bwboundaries(mLocal, 4, 'noholes');
        for bb = 1:numel(bnd)
            plot(axE, bnd{bb}(:,2), bnd{bb}(:,1), '-y', 'LineWidth', 1.0);
        end

        % Symmetric clim centered at 0
        cLim = max(abs(eigenface(:)));
        if cLim == 0, cLim = 1; end
        caxis(axE, [-cLim, cLim]);
        colormap(axE, bluewhitered(256));
        axis(axE, 'image'); axis(axE, 'off');

        % Title with mode number and dipole fraction
        df_val = dipole_frac(roi_id, k);
        if isnan(df_val)
            title(axE, sprintf('mode %d', k), 'Interpreter','none', 'FontSize',7);
        else
            title(axE, sprintf('mode %d, dipole=%.2f', k, df_val), 'Interpreter','none', 'FontSize',7);
        end
    end

    if savePDF
        tmpR = fullfile(tempdir, sprintf('roiQC_roi%03d.pdf', roi_id));
        exportgraphics(figR, tmpR, 'ContentType','vector');
        tmpPDFs{end+1} = tmpR;
    end
    close(figR);
    fprintf('[ROI QC] Saved eigenface page for ROI #%d (%d/%d)\n', roi_id, pp, N_roi);
end

%% ======================================================================
%  SAVE .MAT
%  ======================================================================
results = struct();
results.noise        = noise;
results.skew         = skew;
results.peakSNR      = peakSNR;
results.motion_score = motion_score;
results.dipole_frac  = dipole_frac;
results.dipole_angle = dipole_angle;
results.dom_dipole_angle = dom_dipole_angle;
results.roi_nPix     = roi_nPix;

results.reject_lowSNR   = reject_lowSNR;
results.reject_lowSkew  = reject_lowSkew;
results.reject_motion   = reject_motion;
results.reject_tooSmall = reject_tooSmall;
results.reject_any      = reject_any;
results.keep            = keep;

results.params.folderPath           = folderPath;
results.params.fps                  = fps;
results.params.pixelSize_um         = pixelSize_um;
results.params.patchSize_um         = patchSize_um;
results.params.patchRadius_px       = patchRadius_px;
results.params.nLocalModes          = nLocalModes;
results.params.skewnessThreshold    = skewnessThreshold;
results.params.peakSNR_Threshold    = peakSNR_Threshold;
results.params.motionScoreThreshold = motionScoreThreshold;
results.params.minROIPixels         = minROIPixels;
results.params.T_dff                = T_dff;
results.params.N_roi                = N_roi;
results.params.tossFrames           = tossFrames;

outMat = fullfile(folderPath, 'roiQC_output.mat');
save(outMat, '-struct', 'results', '-v7.3');
fprintf('[ROI QC] Saved: %s\n', outMat);

%% ======================================================================
%  MERGE PDFs
%  ======================================================================
if savePDF && ~isempty(tmpPDFs)
    outPDF = fullfile(folderPath, 'roiQC.pdf');
    if isfile(outPDF), delete(outPDF); end
    append_pdfs_RZ(outPDF, tmpPDFs{:});
    fprintf('[ROI QC] Saved: %s\n', outPDF);

    % Clean up temp files
    for ii = 1:numel(tmpPDFs)
        if isfile(tmpPDFs{ii}), delete(tmpPDFs{ii}); end
    end
end

%% ======================================================================
%  CONSOLE SUMMARY
%  ======================================================================
[maxMS, maxMS_id] = max(motion_score);
nKeep = sum(keep);
nReject = sum(reject_any);

fprintf('\n==========================================================\n');
fprintf('  ROI QC SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Folder:          %s\n', folderPath);
fprintf('  FPS:             %d Hz\n', fps);
fprintf('  Frames:          %d\n', T_dff);
fprintf('  ROIs:            %d\n', N_roi);
fprintf('  Patch size:      %d um (%d px)\n', patchSize_um, 2*patchRadius_px+1);
fprintf('----------------------------------------------------------\n');
fprintf('  Tier 1 -- Trace SNR\n');
fprintf('    Mean skewness:     %.3f\n', mean(skew));
fprintf('    Mean peak SNR:     %.2f\n', mean(peakSNR));
fprintf('    Reject (low SNR):  %d / %d\n', sum(reject_lowSNR), N_roi);
fprintf('    Reject (low skew): %d / %d\n', sum(reject_lowSkew), N_roi);
fprintf('----------------------------------------------------------\n');
fprintf('  Tier 2 -- Motion Dipole\n');
fprintf('    Mean motion_score: %.4f\n', mean(motion_score, 'omitnan'));
fprintf('    Max motion_score:  %.4f  (ROI #%d)\n', maxMS, maxMS_id);
fprintf('    Reject (motion):   %d / %d\n', sum(reject_motion), N_roi);
fprintf('----------------------------------------------------------\n');
fprintf('  Summary\n');
fprintf('    Keep:    %d / %d  (%.1f%%)\n', nKeep, N_roi, 100*nKeep/N_roi);
fprintf('    Reject:  %d / %d\n', nReject, N_roi);
fprintf('==========================================================\n');
fprintf('[ROI QC] Done.\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function overlay_roi_values(ax, maskL, values, cmapType, clims)
%OVERLAY_ROI_VALUES  Draw filled ROI boundaries color-coded by a value vector.
%
%   cmapType: 'hot' for positive-only, 'bwr' for bluewhitered diverging

    N = numel(values);
    nColors = 256;

    if strcmp(cmapType, 'bwr')
        cmap = bluewhitered(nColors);
    else
        cmap = hot(nColors);
    end

    % Map values to colormap indices
    valNorm = (values - clims(1)) / (clims(2) - clims(1));
    valNorm = max(0, min(1, valNorm));
    cIdx    = round(valNorm * (nColors - 1)) + 1;

    % Get actual ROI labels from mask (may not be consecutive 1:N)
    labels = unique(maskL(:));
    labels(labels == 0) = [];

    % Build RGBA overlay
    [H, W] = size(maskL);
    ovR = zeros(H, W);
    ovG = zeros(H, W);
    ovB = zeros(H, W);
    ovA = zeros(H, W);

    for i = 1:min(N, numel(labels))
        roiMask = (maskL == labels(i));
        if ~any(roiMask, 'all'), continue; end

        c = cmap(cIdx(i), :);

        % Fill ROI interior with translucent color
        ovR(roiMask) = c(1);
        ovG(roiMask) = c(2);
        ovB(roiMask) = c(3);
        ovA(roiMask) = 0.55;

        % Draw boundary
        bnd = bwboundaries(roiMask, 4, 'noholes');
        for bb = 1:numel(bnd)
            plot(ax, bnd{bb}(:,2), bnd{bb}(:,1), '-', ...
                'Color', c, 'LineWidth', 1.2);
        end
    end

    ovRGB = cat(3, ovR, ovG, ovB);
    hImg = image(ax, ovRGB, 'AlphaData', ovA);
    uistack(hImg, 'top');

    % Force the axes colormap and clim for the colorbar
    colormap(ax, cmap);
    caxis(ax, clims);
end
