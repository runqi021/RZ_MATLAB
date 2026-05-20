clear all; clc; close all;
%% Run_ROI_motion_diagnostic.m
%  Quantify how much residual motion-correction displacement leaks into
%  each ROI's fluorescence trace.  Works entirely on pre-computed .mat
%  files — no TIFF loading.
%
%  Outputs (in folderPath):
%    motionDiag_output.mat  — all numeric results
%    motionDiag.pdf         — 3-page diagnostic report
%
%  RZ 2026-04-01

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
folderPath         = "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing\roi5_1400-1300-0_18lp_930_x4_512x512_3000f_00001";      % <-- paste data subfolder path here
corrThreshold      = 0.10;    % R^2 threshold to label an ROI "contaminated"
nTopROIs           = 5;       % number of worst ROIs in temporal figure
nSVDModes          = 10;      % SVD modes to classify as motion vs neural
modeCorr_threshold = 0.4;     % |corr| with shifts to call a mode "motion"
savePDF            = true;

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
%  DETECT FPS
%  ======================================================================
[fps, ~] = detect_session_fps(folderPath);
fprintf('Detected FPS = %d\n', fps);

%% ======================================================================
%  LOAD dFF DATA
%  ======================================================================
dff_hits = dir(fullfile(folderPath, '*_dFF.mat'));
assert(~isempty(dff_hits), 'No *_dFF.mat found in %s', folderPath);
dff_path = fullfile(dff_hits(1).folder, dff_hits(1).name);
fprintf('Loading dFF:  %s\n', dff_path);

dffData = load(dff_path);
assert(isfield(dffData, 'F_roi'), '_dFF.mat missing F_roi.');
assert(isfield(dffData, 'dFF'),   '_dFF.mat missing dFF.');
assert(isfield(dffData, 'params'),'_dFF.mat missing params.');

F_roi = double(dffData.F_roi);          % [T_dff x N]  post-toss raw F
[T_dff, N_roi] = size(F_roi);

tossFrames = 0;
if isfield(dffData.params, 'tossFrames')
    tossFrames = dffData.params.tossFrames;
end
fprintf('F_roi: [%d x %d],  tossFrames = %d\n', T_dff, N_roi, tossFrames);

%% ======================================================================
%  LOAD cpSAM (maskL + template fallback)
%  ======================================================================
sam_hits = dir(fullfile(folderPath, '*_cpSAM_output.mat'));
assert(~isempty(sam_hits), 'No *_cpSAM_output.mat found in %s', folderPath);
sam_path = fullfile(sam_hits(1).folder, sam_hits(1).name);
fprintf('Loading SAM:  %s\n', sam_path);

samData = load(sam_path);
assert(isfield(samData, 'maskL'), 'SAM file missing maskL.');
maskL = samData.maskL;   % [H x W] uint16 label image

%% ======================================================================
%  EXTRACT MC SHIFTS (sum across passes)
%  ======================================================================
[sx_total, sy_total, refImage] = extract_mc_shifts(folderPath);
T_shifts = numel(sx_total);
fprintf('MC shifts: T_shifts = %d  (across %s)\n', T_shifts, folderPath);

% --- Frame alignment: trim leading shift frames if MC ran before toss ---
if T_shifts > T_dff
    nTrim = T_shifts - T_dff;
    fprintf('Trimming first %d shift frames to align with F_roi (toss=%d).\n', nTrim, tossFrames);
    sx_total(1:nTrim) = [];
    sy_total(1:nTrim) = [];
    T_shifts = numel(sx_total);
end
assert(T_shifts == T_dff, ...
    'Frame mismatch after trim: T_shifts=%d, T_dff=%d', T_shifts, T_dff);

shift_mag = hypot(sx_total, sy_total);
t_sec = (0:T_dff-1)' / fps;

%% ======================================================================
%  SECTION A: PER-ROI MOTION CONTAMINATION (on raw F)
%  ======================================================================
fprintf('\n--- Section A: Per-ROI motion contamination ---\n');

rho_x   = zeros(N_roi, 1);
rho_y   = zeros(N_roi, 1);
rho_mag = zeros(N_roi, 1);
R2      = zeros(N_roi, 1);

X_reg = [ones(T_dff, 1), sx_total(:), sy_total(:)];

for i = 1:N_roi
    fi = F_roi(:, i);
    rho_x(i)   = corr(fi, sx_total(:));
    rho_y(i)   = corr(fi, sy_total(:));
    rho_mag(i) = corr(fi, shift_mag(:));

    % Multiple regression R^2
    b    = X_reg \ fi;
    yhat = X_reg * b;
    SS_res  = sum((fi - yhat).^2);
    SS_tot  = sum((fi - mean(fi)).^2);
    if SS_tot > 0
        R2(i) = 1 - SS_res / SS_tot;
    end
end

nContaminated = sum(R2 > corrThreshold);
fprintf('ROIs contaminated (R^2 > %.2f): %d / %d (%.1f%%)\n', ...
    corrThreshold, nContaminated, N_roi, 100*nContaminated/N_roi);
fprintf('Mean R^2 = %.4f,  Max R^2 = %.4f\n', mean(R2), max(R2));

%% ======================================================================
%  SECTION B: SVD ON ROI F MATRIX
%  ======================================================================
fprintf('\n--- Section B: SVD on F_roi matrix ---\n');

F_c = F_roi - mean(F_roi, 1);
nModes = min(nSVDModes, min(T_dff, N_roi));

[U, S_mat, ~] = svd(F_c, 'econ');

singVals = diag(S_mat);
varExp   = singVals.^2 / sum(singVals.^2);

svd_corr_x   = zeros(nModes, 1);
svd_corr_y   = zeros(nModes, 1);
svd_corr_mag = zeros(nModes, 1);
isMotion      = false(nModes, 1);

for k = 1:nModes
    temporal_k       = U(:, k);   % [T x 1] temporal mode
    svd_corr_x(k)   = corr(temporal_k, sx_total(:));
    svd_corr_y(k)   = corr(temporal_k, sy_total(:));
    svd_corr_mag(k) = corr(temporal_k, shift_mag(:));
    isMotion(k) = max(abs([svd_corr_x(k), svd_corr_y(k), svd_corr_mag(k)])) ...
                  > modeCorr_threshold;
end

motionVar_frac = sum(varExp(isMotion)) / sum(varExp(1:nModes));
fprintf('Motion modes: %d / %d  (%.1f%% of top-%d variance)\n', ...
    sum(isMotion), nModes, 100*motionVar_frac, nModes);

%% ======================================================================
%  SECTION C: FIGURES
%  ======================================================================
fprintf('\n--- Section C: Generating figures ---\n');
tmpPDFs = {};

% ---- Figure 1: Overview (2x2) ----
fig1 = figure('Visible','off', 'Color','w', 'Units','inches', 'Position',[1 1 11 8.5]);
t1 = tiledlayout(fig1, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(t1, 'ROI Motion Diagnostic — Overview', 'Interpreter','none', 'FontSize',14);

% (1,1) R^2 bar chart, sorted descending
ax1 = nexttile(t1);
[R2_sorted, sortIdx] = sort(R2, 'descend');
bar(ax1, R2_sorted, 'FaceColor',[0.3 0.3 0.8], 'EdgeColor','none');
hold(ax1, 'on');
yline(ax1, corrThreshold, '--r', 'LineWidth', 1.5);
xlabel(ax1, 'ROI (sorted)', 'Interpreter','none');
ylabel(ax1, 'R^2 (shift regression)', 'Interpreter','none');
title(ax1, 'R^2 per ROI (sorted)', 'Interpreter','none');
xlim(ax1, [0.5, N_roi+0.5]);
set(ax1, 'Box','off');

% (1,2) SVD singular value spectrum
ax2 = nexttile(t1);
hold(ax2, 'on');
for k = 1:nModes
    if isMotion(k)
        bar(ax2, k, 100*varExp(k), 'FaceColor','r', 'EdgeColor','none');
    else
        bar(ax2, k, 100*varExp(k), 'FaceColor',[0.2 0.4 0.8], 'EdgeColor','none');
    end
end
xlabel(ax2, 'SVD mode', 'Interpreter','none');
ylabel(ax2, 'Variance explained (%)', 'Interpreter','none');
title(ax2, 'SVD spectrum (red = motion)', 'Interpreter','none');
set(ax2, 'YScale','log', 'Box','off');
xlim(ax2, [0.5, nModes+0.5]);

% (2,1) Max |rho| per ROI bar chart
ax3 = nexttile(t1);
maxAbsRho = max(abs([rho_x, rho_y, rho_mag]), [], 2);
[maxRho_sorted, ~] = sort(maxAbsRho, 'descend');
bar(ax3, maxRho_sorted, 'FaceColor',[0.8 0.4 0.2], 'EdgeColor','none');
hold(ax3, 'on');
yline(ax3, sqrt(corrThreshold), '--r', 'LineWidth', 1.5);
xlabel(ax3, 'ROI (sorted)', 'Interpreter','none');
ylabel(ax3, 'max |rho|', 'Interpreter','none');
title(ax3, 'Max |Pearson r| per ROI (sorted)', 'Interpreter','none');
xlim(ax3, [0.5, N_roi+0.5]);
set(ax3, 'Box','off');

% (2,2) Text summary
ax4 = nexttile(t1);
axis(ax4, 'off');
summaryStr = { ...
    sprintf('N ROIs:            %d', N_roi), ...
    sprintf('FPS:               %d Hz', fps), ...
    sprintf('Frames (post-toss): %d', T_dff), ...
    '', ...
    sprintf('R^2 threshold:     %.2f', corrThreshold), ...
    sprintf('Contaminated ROIs: %d / %d (%.1f%%)', nContaminated, N_roi, 100*nContaminated/N_roi), ...
    sprintf('Mean R^2:          %.4f', mean(R2)), ...
    sprintf('Max  R^2:          %.4f', max(R2)), ...
    '', ...
    sprintf('SVD modes analyzed: %d', nModes), ...
    sprintf('Motion modes:       %d  (|r| > %.2f)', sum(isMotion), modeCorr_threshold), ...
    sprintf('Motion var frac:    %.1f%% of top-%d', 100*motionVar_frac, nModes), ...
    };
text(ax4, 0.05, 0.95, summaryStr, ...
    'VerticalAlignment','top', 'FontSize',11, 'FontName','FixedWidth', ...
    'Interpreter','none');

if savePDF
    tmp1 = fullfile(tempdir, 'motionDiag_fig1.pdf');
    exportgraphics(fig1, tmp1, 'ContentType','vector');
    tmpPDFs{end+1} = tmp1;
end
close(fig1);

% ---- Figure 2: Spatial Map (1x2) ----
fig2 = figure('Visible','off', 'Color','w', 'Units','inches', 'Position',[1 1 12 5.5]);
t2 = tiledlayout(fig2, 1, 2, 'TileSpacing','compact', 'Padding','compact');
title(t2, 'ROI Motion Diagnostic — Spatial Map', 'Interpreter','none', 'FontSize',14);

% Left: R^2 map (hot colormap)
ax5 = nexttile(t2);
if ~isempty(refImage)
    bgImg = mat2gray(double(refImage));
else
    bgImg = zeros(size(maskL));
end
imagesc(ax5, bgImg); colormap(ax5, gray(256)); hold(ax5, 'on');
axis(ax5, 'image'); axis(ax5, 'off');
overlay_roi_values(ax5, maskL, R2, 'hot', [0, max(max(R2), eps)]);
cb1 = colorbar(ax5);
cb1.Label.String = 'R^2';
cb1.Label.Interpreter = 'none';
title(ax5, 'R^2 (shift regression)', 'Interpreter','none');

% Right: signed rho_mag map (bluewhitered)
ax6 = nexttile(t2);
imagesc(ax6, bgImg); colormap(ax6, gray(256)); hold(ax6, 'on');
axis(ax6, 'image'); axis(ax6, 'off');
rhoLim = max(abs(rho_mag));
if rhoLim == 0, rhoLim = 1; end
overlay_roi_values(ax6, maskL, rho_mag, 'bwr', [-rhoLim, rhoLim]);
cb2 = colorbar(ax6);
cb2.Label.String = 'rho (shift mag)';
cb2.Label.Interpreter = 'none';
title(ax6, 'Pearson r (shift magnitude)', 'Interpreter','none');

if savePDF
    tmp2 = fullfile(tempdir, 'motionDiag_fig2.pdf');
    exportgraphics(fig2, tmp2, 'ContentType','vector');
    tmpPDFs{end+1} = tmp2;
end
close(fig2);

% ---- Figure 3: Temporal Traces (nTopROIs+1 rows) ----
nShow = min(nTopROIs, N_roi);
[~, topIdx] = sort(R2, 'descend');
topIdx = topIdx(1:nShow);

fig3 = figure('Visible','off', 'Color','w', 'Units','inches', ...
    'Position',[1 1 14 2.5*(nShow+1)]);
t3 = tiledlayout(fig3, nShow+1, 1, 'TileSpacing','tight', 'Padding','compact');
title(t3, 'ROI Motion Diagnostic — Temporal Traces', 'Interpreter','none', 'FontSize',14);

% Top panel: MC shifts
axShift = nexttile(t3);
hold(axShift, 'on');
plot(axShift, t_sec, sx_total, 'Color',[0.2 0.4 0.9], 'LineWidth',0.8);
plot(axShift, t_sec, sy_total, 'Color',[0.9 0.2 0.2], 'LineWidth',0.8);
plot(axShift, t_sec, shift_mag, 'Color','k', 'LineWidth',1.0);
ylabel(axShift, 'Shift (px)', 'Interpreter','none');
title(axShift, 'MC shifts: sx (blue), sy (red), |s| (black)', 'Interpreter','none');
legend(axShift, {'sx','sy','|s|'}, 'Location','northeast', 'Box','off');
set(axShift, 'Box','off');
axHandles = axShift;

% Per-ROI panels
for jj = 1:nShow
    roi_id = topIdx(jj);
    fi = F_roi(:, roi_id);
    b  = X_reg \ fi;
    yhat = X_reg * b;

    axR = nexttile(t3);
    hold(axR, 'on');
    plot(axR, t_sec, fi,   'Color','k',               'LineWidth',0.6);
    plot(axR, t_sec, yhat, 'Color',[0.9 0.15 0.15],   'LineWidth',1.0);
    ylabel(axR, 'F (a.u.)', 'Interpreter','none');
    title(axR, sprintf('ROI #%d,  R^2 = %.3f', roi_id, R2(roi_id)), ...
        'Interpreter','none');
    set(axR, 'Box','off');
    axHandles(end+1) = axR; %#ok<SAGROW>
end
xlabel(axHandles(end), 'Time (s)', 'Interpreter','none');
linkaxes(axHandles, 'x');
xlim(axHandles(1), [t_sec(1), t_sec(end)]);

if savePDF
    tmp3 = fullfile(tempdir, 'motionDiag_fig3.pdf');
    exportgraphics(fig3, tmp3, 'ContentType','vector');
    tmpPDFs{end+1} = tmp3;
end
close(fig3);

%% ======================================================================
%  SECTION D: SAVE
%  ======================================================================
fprintf('\n--- Section D: Saving results ---\n');

% .mat output
results = struct();
results.R2        = R2;
results.rho_x     = rho_x;
results.rho_y     = rho_y;
results.rho_mag   = rho_mag;
results.sortIdx_R2 = sortIdx;
results.nContaminated = nContaminated;

results.svd_varExp     = varExp(1:nModes);
results.svd_corr_x     = svd_corr_x;
results.svd_corr_y     = svd_corr_y;
results.svd_corr_mag   = svd_corr_mag;
results.svd_isMotion   = isMotion;
results.motionVar_frac = motionVar_frac;

results.params.folderPath         = folderPath;
results.params.corrThreshold      = corrThreshold;
results.params.nTopROIs           = nTopROIs;
results.params.nSVDModes          = nSVDModes;
results.params.modeCorr_threshold = modeCorr_threshold;
results.params.fps                = fps;
results.params.T_dff              = T_dff;
results.params.N_roi              = N_roi;
results.params.tossFrames         = tossFrames;

results.sx = sx_total;
results.sy = sy_total;

outMat = fullfile(folderPath, 'motionDiag_output.mat');
save(outMat, '-struct', 'results', '-v7.3');
fprintf('Saved: %s\n', outMat);

% Merge PDFs
if savePDF && ~isempty(tmpPDFs)
    outPDF = fullfile(folderPath, 'motionDiag.pdf');
    % Delete stale output if it exists so append_pdfs_RZ starts fresh
    if isfile(outPDF), delete(outPDF); end
    append_pdfs_RZ(outPDF, tmpPDFs{:});
    fprintf('Saved: %s\n', outPDF);

    % Clean up temp files
    for ii = 1:numel(tmpPDFs)
        if isfile(tmpPDFs{ii}), delete(tmpPDFs{ii}); end
    end
end

%% ======================================================================
%  CONSOLE SUMMARY
%  ======================================================================
fprintf('\n==========================================================\n');
fprintf('  ROI MOTION DIAGNOSTIC SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Folder:          %s\n', folderPath);
fprintf('  FPS:             %d Hz\n', fps);
fprintf('  Frames:          %d  (toss=%d)\n', T_dff, tossFrames);
fprintf('  ROIs:            %d\n', N_roi);
fprintf('----------------------------------------------------------\n');
fprintf('  R^2 threshold:   %.2f\n', corrThreshold);
fprintf('  Contaminated:    %d / %d  (%.1f%%)\n', ...
    nContaminated, N_roi, 100*nContaminated/N_roi);
fprintf('  Mean R^2:        %.4f\n', mean(R2));
fprintf('  Max  R^2:        %.4f  (ROI #%d)\n', max(R2), sortIdx(1));
fprintf('  Median |rho|:    %.4f\n', median(maxAbsRho));
fprintf('----------------------------------------------------------\n');
fprintf('  SVD modes:       %d analyzed,  %d classified as motion\n', nModes, sum(isMotion));
fprintf('  Motion var frac: %.1f%% of top-%d variance\n', 100*motionVar_frac, nModes);
fprintf('==========================================================\n');
fprintf('Done.\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function [sx, sy, refImg] = extract_mc_shifts(folder)
%EXTRACT_MC_SHIFTS  Load and sum MC shifts across all passes.
%
%   Strategy:
%     1. Find all *_MC_output.mat files (rigid MC passes). Each contains
%        mcOut.shifts. Sum shifts across passes for total displacement.
%     2. Fallback: if no _MC_output.mat, try *_MCinfo.mat (batch pipeline)
%        which has top-level sx, sy, r (last pass only).
%     3. Also check *_SPNR_output.mat (non-rigid pass). For SPNR, the
%        per-frame shift is the mean across patches.

    refImg = [];

    % --- Strategy 1: individual _MC_output.mat and _SPNR_output.mat ---
    mc_hits   = dir(fullfile(folder, '*_MC_output.mat'));
    spnr_hits = dir(fullfile(folder, '*_SPNR_output.mat'));
    all_hits  = [mc_hits; spnr_hits];

    if ~isempty(all_hits)
        sx = [];
        sy = [];
        for ff = 1:numel(all_hits)
            fpath = fullfile(all_hits(ff).folder, all_hits(ff).name);
            tmp = load(fpath, 'mcOut');
            assert(isfield(tmp, 'mcOut'), 'File %s missing mcOut.', fpath);
            mc = tmp.mcOut;

            [sx_pass, sy_pass] = shifts_from_mcOut(mc);
            fprintf('  Loaded shifts from: %s  (%d frames)\n', all_hits(ff).name, numel(sx_pass));

            if isempty(sx)
                sx = sx_pass;
                sy = sy_pass;
            else
                % Sum across passes — trim to common length
                nCommon = min(numel(sx), numel(sx_pass));
                sx = sx(1:nCommon) + sx_pass(1:nCommon);
                sy = sy(1:nCommon) + sy_pass(1:nCommon);
            end

            % Keep template from the last pass
            if isfield(mc, 'template') && ~isempty(mc.template)
                refImg = mc.template;
            end
        end
        return
    end

    % --- Strategy 2 (fallback): *_MCinfo.mat ---
    info_hits = dir(fullfile(folder, '*_MCinfo.mat'));
    if ~isempty(info_hits)
        fpath = fullfile(info_hits(1).folder, info_hits(1).name);
        tmp = load(fpath);
        fprintf('  [Fallback] Using MCinfo: %s\n', info_hits(1).name);

        if isfield(tmp, 'sx') && isfield(tmp, 'sy') && ~isempty(tmp.sx)
            sx = double(tmp.sx(:));
            sy = double(tmp.sy(:));
        elseif isfield(tmp, 'mcOut')
            [sx, sy] = shifts_from_mcOut(tmp.mcOut);
        else
            error('MCinfo file has neither sx/sy nor mcOut: %s', fpath);
        end

        if isfield(tmp, 'mcOut') && isfield(tmp.mcOut, 'template')
            refImg = tmp.mcOut.template;
        end
        return
    end

    error('No MC shift files (*_MC_output.mat, *_SPNR_output.mat, *_MCinfo.mat) found in:\n%s', folder);
end


function [sx, sy] = shifts_from_mcOut(mc)
%SHIFTS_FROM_MCOUT  Extract per-frame sx, sy from an mcOut struct.
%   For rigid MC, shifts(t).shifts is [1x2].
%   For SPNR (non-rigid), shifts(t).shifts is [nPatchY x nPatchX x 2]
%   (or similar multi-dim). We take the mean across patches.

    sh = mc.shifts;
    T  = numel(sh);
    sx = zeros(T, 1);
    sy = zeros(T, 1);

    for t = 1:T
        sft = sh(t).shifts;
        if isvector(sft) && numel(sft) == 2
            % Rigid: [1x2] or [2x1]
            sx(t) = sft(1);
            sy(t) = sft(2);
        else
            % Non-rigid: multi-patch, take FOV-mean shift
            % NoRMCorre convention: shifts are stored with dim order
            % matching the image dims, plus one trailing dim for x/y.
            % For 2D: [nPatchY x nPatchX x 2]
            sx(t) = mean(sft(:,:,1), 'all');
            if size(sft, 3) >= 2
                sy(t) = mean(sft(:,:,2), 'all');
            end
        end
    end
end


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
