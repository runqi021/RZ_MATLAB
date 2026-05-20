clear all; clc; close all;
%% Run_ROI_QC_Pass1.m
%  Two-pass ROI QC -- Pass 1: Static avg-projection grid.
%
%  For each FOV: generates a per-ROI montage grid (gamma-corrected avg
%  projection crop with yellow ROI boundary), a Page 0 overview, and a
%  per-ROI sequential review GUI for keep/reject decisions.
%
%  Two input modes:
%    masterFolder: scan all subfolders for *_cpSAM_output.mat
%    folderPath:   single FOV
%
%  Outputs (per FOV folder):
%    roiQC_pass1_montage.png  -- per-ROI grid image (with rejection marks)
%    roiQC_exclusion.mat      -- exclusion data + morphometrics
%
%  RZ 2026-04

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
masterFolder    = "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing";      % multi-experiment: scans all subfolders
folderPath      = "";      % single FOV mode (ignored if masterFolder set)

patchSize_um    = 30;      % crop half-width per ROI (microns)
gamma_val       = 0.6;     % gamma correction for display
clip_lo         = 0.4;     % brightness clip low percentile
clip_hi         = 99.5;    % brightness clip high percentile
interactiveMode = true;    % true = sequential per-ROI GUI; false = text input

skipExisting    = false;    % skip FOVs with existing roiQC_exclusion.mat

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'mod')));
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% ======================================================================
%  DISCOVER FOVs
%  ======================================================================
samHits = discover_cpSAM(masterFolder, folderPath);
nFOV = numel(samHits);
fprintf('\n=== Found %d FOV(s) with cpSAM output ===\n', nFOV);
for i = 1:nFOV
    fprintf('  [%d] %s\n', i, fullfile(samHits(i).folder, samHits(i).name));
end
fprintf('\n');

%% ======================================================================
%  BATCH LOOP
%  ======================================================================
nProcessed = 0;
nSkipped   = 0;
nFailed    = 0;
failedList = {};

for ff = 1:nFOV
    samPath   = fullfile(samHits(ff).folder, samHits(ff).name);
    fovFolder = samHits(ff).folder;
    [~, fovID] = fileparts(fovFolder);

    fprintf('\n========== [%d/%d] %s ==========\n', ff, nFOV, fovID);

    % Skip if exclusion already exists
    exclPath = fullfile(fovFolder, 'roiQC_exclusion.mat');
    if skipExisting && isfile(exclPath)
        fprintf('  [SKIP] roiQC_exclusion.mat already exists.\n');
        nSkipped = nSkipped + 1;
        continue
    end

    try
        %% ---- Load cpSAM ----
        samData = load(samPath);
        assert(isfield(samData, 'maskL'), 'cpSAM file missing maskL.');
        maskL = samData.maskL;

        props = regionprops(maskL, 'Centroid', 'Area', 'Perimeter');
        validROI = find(arrayfun(@(p) p.Area > 0, props));
        N_roi = numel(validROI);
        fprintf('  ROIs: %d\n', N_roi);

        if N_roi == 0
            fprintf('  No ROIs found. Skipping.\n');
            continue
        end

        %% ---- Detect FPS + pixel size ----
        [fps, scan_meta] = detect_session_fps(fovFolder);
        if isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um)
            pixelSize_um = scan_meta.pixelSize_um;
        else
            pixelSize_um = 1.7778;
            fprintf(2, '  WARNING: pixelSize_um not found. Using fallback %.4f um/px.\n', pixelSize_um);
        end
        crop_px = round(patchSize_um / pixelSize_um / 2);
        patchSz = 2 * crop_px + 1;
        fprintf('  Patch: %d um -> %d px half-width -> %d x %d tile\n', ...
            patchSize_um, crop_px, patchSz, patchSz);

        %% ---- Load avg projection ----
        avgImg = load_avg_projection(samData, fovFolder);
        [H, W] = size(avgImg);
        fprintf('  Avg projection: [%d x %d]\n', H, W);

        %% ---- Morphometrics ----
        roi_nPix        = zeros(N_roi, 1);
        roi_diameter_px = zeros(N_roi, 1);
        roi_diameter_um = zeros(N_roi, 1);
        roi_circularity = zeros(N_roi, 1);
        roi_centroid    = zeros(N_roi, 2);

        for i = 1:N_roi
            rid = validROI(i);
            roi_nPix(i)        = props(rid).Area;
            roi_diameter_px(i) = 2 * sqrt(props(rid).Area / pi);
            roi_diameter_um(i) = roi_diameter_px(i) * pixelSize_um;
            perim = props(rid).Perimeter;
            if perim > 0
                roi_circularity(i) = 4 * pi * props(rid).Area / perim^2;
            end
            roi_centroid(i,:) = [props(rid).Centroid(2), props(rid).Centroid(1)]; % [row, col]
        end

        %% ---- Optional: load dFF for SNR metrics ----
        peakSNR_vec = NaN(N_roi, 1);
        skew_vec    = NaN(N_roi, 1);

        dff_hits = dir(fullfile(fovFolder, '*_dFF.mat'));
        if ~isempty(dff_hits)
            dffData = load(fullfile(dff_hits(1).folder, dff_hits(1).name), 'F_roi', 'dFF');
            if isfield(dffData, 'F_roi') && isfield(dffData, 'dFF')
                F_roi  = double(dffData.F_roi);
                dFF_dt = detrend(double(dffData.dFF));
                nROI_dff = min(N_roi, size(F_roi, 2));
                for i = 1:nROI_dff
                    noise_i = 1.4826 * median(abs(diff(F_roi(:, i)))) / sqrt(2);
                    skew_vec(i) = skewness(dFF_dt(:, i));
                    mean_F = mean(F_roi(:, i));
                    if noise_i > 0 && mean_F > 0
                        peakSNR_vec(i) = max(dFF_dt(:, i)) / (noise_i / mean_F);
                    end
                end
                fprintf('  SNR metrics computed from dFF.\n');
            end
        end

        %% ---- Build montage grid ----
        [montageRGB, gridInfo] = build_roi_grid(avgImg, maskL, validROI, props, ...
            crop_px, patchSz, gamma_val, clip_lo, clip_hi);

        %% ---- Build individual ROI tiles ----
        roiTiles = build_roi_tiles(avgImg, maskL, validROI, props, ...
            crop_px, patchSz, gamma_val, clip_lo, clip_hi);

        %% ---- Load previous exclusion if exists ----
        prev_reject = [];
        prev_toss   = false;
        if isfile(exclPath)
            prev_excl = load(exclPath, 'pass1_reject', 'toss_fov');
            if isfield(prev_excl, 'pass1_reject') && numel(prev_excl.pass1_reject) == N_roi
                prev_reject = prev_excl.pass1_reject;
                if isfield(prev_excl, 'toss_fov')
                    prev_toss = prev_excl.toss_fov;
                end
                fprintf('  Loaded previous Pass 1 decisions (%d rejected).\n', sum(prev_reject));
            end
        end

        %% ---- Display + interactive exclusion ----
        rejected      = false(N_roi, 1);
        reject_reason = repmat("", N_roi, 1);
        toss_fov      = false;

        if interactiveMode
            % Metric vectors for display (NaN-safe)
            metrics = struct('diameter_um', roi_diameter_um, ...
                             'circularity', roi_circularity, ...
                             'peakSNR',     peakSNR_vec, ...
                             'skewness',    skew_vec);

            [rejected, reject_reason, toss_fov] = sequential_static_review( ...
                roiTiles, montageRGB, gridInfo, validROI, N_roi, ...
                patchSz, fovID, metrics, prev_reject, prev_toss);
        else
            % --- Text-input fallback ---
            baseFigSz = 900;
            if N_roi > 80, figScale = 1.5; else, figScale = 1.0; end
            figSz = round(baseFigSz * figScale);

            fig_m = figure('Color', 'k', 'Visible', 'on', ...
                'Name', sprintf('Pass 1 QC - %s (%d ROIs)', fovID, N_roi), ...
                'Position', [50 50 figSz figSz]);
            ax_m = axes(fig_m, 'Units', 'normalized', 'Position', [0.02 0.04 0.96 0.92]);
            image(ax_m, montageRGB);
            axis(ax_m, 'image'); axis(ax_m, 'off');
            title(ax_m, sprintf('%s - %d ROIs', strrep(fovID,'_',' '), N_roi), ...
                'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');

            montagePng = fullfile(fovFolder, 'roiQC_pass1_montage.png');
            exportgraphics(fig_m, montagePng, 'Resolution', 200);
            fprintf('  Montage saved: %s\n', montagePng);
            close(fig_m);

            fprintf('  Review the montage and enter ROI numbers to EXCLUDE.\n');
            excludeStr = input('  Exclude ROIs (comma-separated, Enter to keep all): ', 's');
            if ~isempty(strtrim(excludeStr))
                excludeNums = str2num(excludeStr); %#ok<ST2NM>
                for ei = 1:numel(excludeNums)
                    idx = find(validROI == excludeNums(ei), 1);
                    if ~isempty(idx)
                        rejected(idx) = true;
                        reject_reason(idx) = "manual";
                    end
                end
            end
        end

        % Save montage PNG with rejection marks
        montagePng = fullfile(fovFolder, 'roiQC_pass1_montage.png');
        save_marked_montage(montagePng, montageRGB, gridInfo, validROI, ...
            N_roi, patchSz, rejected, toss_fov);

        %% ---- Save exclusion data ----
        excl = struct();
        excl.fovID           = string(fovID);
        excl.folderPath      = string(fovFolder);
        excl.N_roi           = N_roi;
        excl.roiLabels       = validROI(:);
        excl.toss_fov        = toss_fov;
        excl.pass1_reject    = rejected;
        excl.pass1_reason    = reject_reason;
        excl.pass2_reject    = false(N_roi, 1);
        excl.pass2_reason    = repmat("", N_roi, 1);
        excl.final_keep      = ~rejected;

        excl.roi_nPix        = roi_nPix;
        excl.roi_diameter_px = roi_diameter_px;
        excl.roi_diameter_um = roi_diameter_um;
        excl.roi_circularity = roi_circularity;
        excl.roi_centroid    = roi_centroid;

        excl.peakSNR         = peakSNR_vec;
        excl.skewness        = skew_vec;

        excl.pixelSize_um       = pixelSize_um;
        excl.cellpose_diameter  = NaN;
        if isfield(samData, 'opt') && isstruct(samData.opt) && isfield(samData.opt, 'Diameter')
            excl.cellpose_diameter = samData.opt.Diameter;
        end
        excl.pass1_datetime  = datetime('now');
        excl.pass2_datetime  = NaT;

        save(exclPath, '-struct', 'excl', '-v7.3');

        nKeep = sum(~rejected);
        fprintf('  Exclusion saved: %s\n', exclPath);
        fprintf('  Kept: %d/%d (%.0f%%)  |  Rejected: %d\n', ...
            nKeep, N_roi, 100*nKeep/N_roi, sum(rejected));

        nProcessed = nProcessed + 1;

    catch ME
        fprintf(2, '\n*** FAILED [%d/%d]: %s\n%s\n', ff, nFOV, fovID, ME.message);
        for si = 1:numel(ME.stack)
            fprintf(2, '  %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
        end
        nFailed = nFailed + 1;
        failedList{end+1} = fovID; %#ok<SAGROW>
    end
end

%% ======================================================================
%  CONSOLE SUMMARY
%  ======================================================================
fprintf('\n==================== PASS 1 COMPLETE ====================\n');
fprintf('  Processed: %d\n', nProcessed);
fprintf('  Skipped:   %d\n', nSkipped);
fprintf('  Failed:    %d\n', nFailed);
fprintf('  Total:     %d\n', nFOV);
if nFailed > 0
    fprintf('  Failed FOVs:\n');
    for i = 1:numel(failedList)
        fprintf('    %s\n', failedList{i});
    end
end
fprintf('=========================================================\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function samHits = discover_cpSAM(masterFolder, folderPath)
%DISCOVER_CPSAM  Find all *_cpSAM_output.mat files.
    if masterFolder ~= ""
        samHits = dir(fullfile(masterFolder, '**', '*_cpSAM_output.mat'));
    elseif folderPath ~= ""
        samHits = dir(fullfile(folderPath, '*_cpSAM_output.mat'));
        if isempty(samHits)
            % Try one level of recursion
            samHits = dir(fullfile(folderPath, '**', '*_cpSAM_output.mat'));
        end
    else
        error('Set either masterFolder or folderPath.');
    end
    assert(~isempty(samHits), 'No *_cpSAM_output.mat found.');
    fprintf('[discover] Found %d cpSAM file(s)\n', numel(samHits));
end


function avgImg = load_avg_projection(samData, fovFolder)
%LOAD_AVG_PROJECTION  Load avg projection from cpSAM avgPath or MC TIFF.
    avgImg = [];

    % Try cpSAM avgPath first
    if isfield(samData, 'avgPath')
        ap = char(samData.avgPath);
        if isfile(ap)
            avgImg = double(imread(ap));
            fprintf('  Avg loaded from cpSAM avgPath.\n');
            return
        end
    end

    % Try *_AVG_for_CP.tif in folder
    avg_hits = dir(fullfile(fovFolder, '*_AVG_for_CP.tif'));
    if ~isempty(avg_hits)
        avgImg = double(imread(fullfile(avg_hits(1).folder, avg_hits(1).name)));
        fprintf('  Avg loaded from AVG_for_CP.tif.\n');
        return
    end

    % Fall back: stream from MC TIFF
    tifPath = discover_mc_tiff(fovFolder);
    fprintf('  Streaming avg from MC TIFF...\n');
    tic;
    info = imfinfo(tifPath);
    T = numel(info);
    H = info(1).Height; W = info(1).Width;
    acc = zeros(H, W, 'double');
    for k = 1:T
        acc = acc + double(imread(tifPath, k, 'Info', info));
    end
    avgImg = acc / T;
    fprintf('  Streamed %d frames in %.1f s.\n', T, toc);
end


function tifPath = discover_mc_tiff(fovFolder)
%DISCOVER_MC_TIFF  Find the best MC TIFF in a folder.
    mc_hits = dir(fullfile(fovFolder, '*_MC*.tif'));
    assert(~isempty(mc_hits), 'No *_MC*.tif found in %s', fovFolder);

    % Filter out AVG, QC, shifts, SELROI variants
    names = {mc_hits.name};
    keep = true(size(names));
    for ii = 1:numel(names)
        if contains(names{ii}, '_AVG', 'IgnoreCase', true) || ...
           contains(names{ii}, '_QC',  'IgnoreCase', true) || ...
           contains(names{ii}, '_shifts', 'IgnoreCase', true) || ...
           contains(names{ii}, '_SELROI', 'IgnoreCase', true)
            keep(ii) = false;
        end
    end
    mc_hits = mc_hits(keep);
    assert(~isempty(mc_hits), 'No valid MC TIFF after filtering in %s', fovFolder);

    % Prefer TIFF with most _MC occurrences (2-pass -> *_MC_MC.tif)
    nMC = zeros(numel(mc_hits), 1);
    for ii = 1:numel(mc_hits)
        nMC(ii) = numel(strfind(upper(mc_hits(ii).name), '_MC'));
    end
    [~, bestIdx] = max(nMC);
    tifPath = fullfile(mc_hits(bestIdx).folder, mc_hits(bestIdx).name);
end


function [montageRGB, gridInfo] = build_roi_grid(avgImg, maskL, validROI, props, ...
        crop_px, patchSz, gamma_val, clip_lo, clip_hi)
%BUILD_ROI_GRID  Assemble per-ROI crop tiles into a montage RGB image.
    N_roi = numel(validROI);
    [H, W] = size(avgImg);

    nCols  = ceil(sqrt(N_roi));
    nRows  = ceil(N_roi / nCols);
    border = 2;
    tileW  = patchSz + border;
    tileH  = patchSz + border;
    montageW = nCols * tileW + border;
    montageH = nRows * tileH + border;

    montageRGB = zeros(montageH, montageW, 3, 'uint8');

    tileRow = zeros(N_roi, 1);
    tileCol = zeros(N_roi, 1);

    for i = 1:N_roi
        rid = validROI(i);
        cx = props(rid).Centroid(1);  % col
        cy = props(rid).Centroid(2);  % row

        % Crop bounds (clipped to image edges)
        r1 = max(1,  round(cy) - crop_px);
        r2 = min(H,  round(cy) + crop_px);
        c1 = max(1,  round(cx) - crop_px);
        c2 = min(W,  round(cx) + crop_px);

        % Brightness normalization + gamma (matching batch_breath_trig_combined)
        rp = avgImg(r1:r2, c1:c2);
        lo = prctile(rp(:), clip_lo);
        hi = prctile(rp(:), clip_hi);
        if hi <= lo, hi = lo + eps; end
        rp = max(0, min(1, (rp - lo) / (hi - lo)));
        rp = rp .^ gamma_val;

        % Yellow ROI boundary overlay
        mc    = maskL(r1:r2, c1:c2);
        perim = bwperim(mc == rid);
        rpRGB = repmat(rp, [1 1 3]);
        for ch = 1:3
            pl = rpRGB(:,:,ch);
            pl(perim) = [1 1 0] * (ch == [1;2;3]);
            rpRGB(:,:,ch) = pl;
        end

        % Pad to uniform tile size (edge ROIs may be clipped)
        tile = zeros(patchSz, patchSz, 3);
        pH = r2 - r1 + 1;
        pW = c2 - c1 + 1;
        tile(1:pH, 1:pW, :) = rpRGB;

        % Grid position
        gi = ceil(i / nCols);
        gj = mod(i - 1, nCols) + 1;
        tileRow(i) = gi;
        tileCol(i) = gj;

        y0 = (gi - 1) * tileH + border + 1;
        x0 = (gj - 1) * tileW + border + 1;

        montageRGB(y0:y0+patchSz-1, x0:x0+patchSz-1, :) = uint8(255 * tile);
    end

    gridInfo.nCols   = nCols;
    gridInfo.nRows   = nRows;
    gridInfo.border  = border;
    gridInfo.tileW   = tileW;
    gridInfo.tileH   = tileH;
    gridInfo.tileRow = tileRow;
    gridInfo.tileCol = tileCol;
end


function roiTiles = build_roi_tiles(avgImg, maskL, validROI, props, ...
        crop_px, patchSz, gamma_val, clip_lo, clip_hi)
%BUILD_ROI_TILES  Pre-build individual ROI crop tiles as a cell array.
%  Same crop/gamma/boundary logic as build_roi_grid but stored per-ROI
%  for large single-ROI display in sequential review.
%
%  Returns roiTiles{i} = [patchSz x patchSz x 3] uint8 RGB.
    N_roi = numel(validROI);
    [H, W] = size(avgImg);
    roiTiles = cell(N_roi, 1);

    for i = 1:N_roi
        rid = validROI(i);
        cx = props(rid).Centroid(1);  % col
        cy = props(rid).Centroid(2);  % row

        % Crop bounds (clipped to image edges)
        r1 = max(1,  round(cy) - crop_px);
        r2 = min(H,  round(cy) + crop_px);
        c1 = max(1,  round(cx) - crop_px);
        c2 = min(W,  round(cx) + crop_px);

        % Brightness normalization + gamma
        rp = avgImg(r1:r2, c1:c2);
        lo = prctile(rp(:), clip_lo);
        hi = prctile(rp(:), clip_hi);
        if hi <= lo, hi = lo + eps; end
        rp = max(0, min(1, (rp - lo) / (hi - lo)));
        rp = rp .^ gamma_val;

        % Yellow ROI boundary overlay
        mc    = maskL(r1:r2, c1:c2);
        perim = bwperim(mc == rid);
        rpRGB = repmat(rp, [1 1 3]);
        for ch = 1:3
            pl = rpRGB(:,:,ch);
            pl(perim) = [1 1 0] * (ch == [1;2;3]);
            rpRGB(:,:,ch) = pl;
        end

        % Pad to uniform tile size (edge ROIs may be clipped)
        tile = zeros(patchSz, patchSz, 3);
        pH = r2 - r1 + 1;
        pW = c2 - c1 + 1;
        tile(1:pH, 1:pW, :) = rpRGB;

        roiTiles{i} = uint8(255 * tile);
    end
end


function [rejected, reject_reason, toss_fov] = sequential_static_review( ...
        roiTiles, montageRGB, gridInfo, validROI, N_roi, patchSz, fovID, ...
        metrics, prev_reject, prev_toss)
%SEQUENTIAL_STATIC_REVIEW  Per-ROI sequential review + grid confirmation.
%
%  Phase 1 (sequential): shows one ROI tile at a time, large view.
%    K = keep, R = reject, U = undo, Left/Right = navigate,
%    T = toss FOV (toggle), Enter = go to grid confirmation, Esc = cancel.
%
%  Phase 2 (grid confirmation): shows the full montage with white X marks
%    on rejected ROIs. L-click = toggle. Enter = confirm, Esc = back.
%
%  Decision array: 1=keep (default), -1=reject.

    rejected      = false(N_roi, 1);
    reject_reason = repmat("", N_roi, 1);
    toss_fov      = false;

    % Decision: 1=keep (default), -1=reject
    decision = ones(N_roi, 1);

    % Load previous decisions if available
    if ~isempty(prev_reject) && numel(prev_reject) == N_roi
        decision(prev_reject) = -1;
        toss_fov = prev_toss;
        fprintf('  Resuming with previous decisions.\n');
    end

    % Figure sizing (same as old grid figure)
    baseFigSz = 900;
    if N_roi > 80, figScale = 1.5; else, figScale = 1.0; end
    figSz = round(baseFigSz * figScale);

    % White cross pointer for grid phase
    crossCData = NaN(16, 16);
    crossCData(8, :) = 2; crossCData(:, 8) = 2;  % 2 = white in pointer colormap
    crossCData(8, 8) = 2;

    fig = figure('Color', 'k', 'Visible', 'on', ...
        'Name', sprintf('Pass 1 QC - %s (%d ROIs)', fovID, N_roi), ...
        'Position', [50 50 figSz figSz], ...
        'MenuBar', 'none', 'ToolBar', 'none');

    fig.UserData = '';
    fig.KeyPressFcn = @(~, evt) set(gcbf, 'UserData', evt.Key);

    cur = 1;          % current ROI index (1-based into N_roi)
    phase = 1;        % 1 = sequential, 2 = grid confirmation
    needRedraw = true;
    done = false;
    cancelled = false;

    while ~done && isvalid(fig)
        if needRedraw
            clf(fig);
            fig.Color = 'k';
            fig.UserData = '';

            if phase == 1
                fig.Pointer = 'arrow';
            else
                set(fig, 'Pointer', 'custom', ...
                    'PointerShapeCData', crossCData, ...
                    'PointerShapeHotSpot', [8 8]);
            end

            if phase == 1
                % ========== PHASE 1: Sequential per-ROI view ==========
                rid = validROI(cur);

                % Main axes: show ROI tile large
                ax_img = axes(fig, 'Units', 'normalized', ...
                    'Position', [0.05 0.12 0.90 0.78]);
                imshow(roiTiles{cur}, 'Parent', ax_img);
                axis(ax_img, 'image');

                % Decision color for title
                if decision(cur) == 1
                    decStr = 'KEEP'; decColor = [0.2 0.9 0.2];
                elseif decision(cur) == -1
                    decStr = 'REJECT'; decColor = [1.0 0.3 0.3];
                else
                    decStr = 'KEEP'; decColor = [0.2 0.9 0.2];
                end

                % Title with ROI# and decision
                if toss_fov
                    titleStr = sprintf('TOSS FOV | ROI #%d (%d/%d) [%s]', ...
                        rid, cur, N_roi, decStr);
                    titleColor = [1.0 0.3 0.3];
                else
                    titleStr = sprintf('ROI #%d (%d/%d) [%s]', ...
                        rid, cur, N_roi, decStr);
                    titleColor = decColor;
                end
                title(ax_img, titleStr, 'Color', titleColor, ...
                    'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

                % Metrics text at bottom of image (NaN-safe)
                metricParts = {};
                if isfinite(metrics.diameter_um(cur))
                    metricParts{end+1} = sprintf('Diam: %.1f um', metrics.diameter_um(cur)); %#ok<AGROW>
                end
                if isfinite(metrics.circularity(cur))
                    metricParts{end+1} = sprintf('Circ: %.2f', metrics.circularity(cur)); %#ok<AGROW>
                end
                if isfinite(metrics.peakSNR(cur))
                    metricParts{end+1} = sprintf('SNR: %.1f', metrics.peakSNR(cur)); %#ok<AGROW>
                end
                if isfinite(metrics.skewness(cur))
                    metricParts{end+1} = sprintf('Skew: %.2f', metrics.skewness(cur)); %#ok<AGROW>
                end
                if ~isempty(metricParts)
                    metricStr = strjoin(metricParts, '  |  ');
                    annotation(fig, 'textbox', [0.05 0.04 0.90 0.06], ...
                        'String', metricStr, 'Color', [0.7 0.7 0.7], ...
                        'EdgeColor', 'none', 'FontSize', 9, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FitBoxToText', 'off');
                end

                % Status bar at top
                nDecided = sum(decision ~= 0);
                nReject  = sum(decision == -1);
                nKeep    = N_roi - nReject;
                statusStr = sprintf('%d/%d decided (%d keep, %d reject)  |  [K]eep [R]eject [U]ndo  |  Left/Right=nav  |  [T]oss FOV  |  Enter=confirm  |  Esc=cancel', ...
                    nDecided, N_roi, nKeep, nReject);
                annotation(fig, 'textbox', [0.0 0.93 1.0 0.07], ...
                    'String', statusStr, 'Color', [0.8 0.8 0.8], ...
                    'EdgeColor', 'none', 'FontSize', 8, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FitBoxToText', 'off');

            else
                % ========== PHASE 2: Grid confirmation ==========
                draw_grid_view(fig, montageRGB, gridInfo, validROI, ...
                    N_roi, patchSz, decision, toss_fov);
            end

            needRedraw = false;
            drawnow;
        end

        if phase == 1
            % ---- Phase 1: poll keypress ----
            pause(0.05);
            key = fig.UserData;
            if isempty(key), continue; end
            fig.UserData = '';

            switch key
                case 'k'  % Keep
                    decision(cur) = 1;
                    if cur < N_roi, cur = cur + 1; end
                    needRedraw = true;

                case 'r'  % Reject
                    decision(cur) = -1;
                    if cur < N_roi, cur = cur + 1; end
                    needRedraw = true;

                case 'u'  % Undo -> keep
                    decision(cur) = 1;
                    needRedraw = true;

                case 'rightarrow'
                    if cur < N_roi, cur = cur + 1; needRedraw = true; end

                case 'leftarrow'
                    if cur > 1, cur = cur - 1; needRedraw = true; end

                case 't'  % Toss FOV toggle
                    toss_fov = ~toss_fov;
                    if toss_fov
                        decision(:) = -1;
                    else
                        decision(:) = 1;
                    end
                    needRedraw = true;

                case 'return'  % Enter -> go to grid confirmation
                    phase = 2;
                    needRedraw = true;

                case 'escape'  % Cancel all
                    cancelled = true;
                    done = true;
            end

        else
            % ---- Phase 2: click-to-toggle + keys (white pointer) ----
            try
                w = waitforbuttonpress;
            catch
                done = true; continue  % figure closed
            end
            if w == 1  % key press
                ch = double(get(fig, 'CurrentCharacter'));
                if isempty(ch) || ch == 13  % Enter = confirm
                    done = true;
                elseif ch == 27  % Esc = back to sequential
                    phase = 1;
                    needRedraw = true;
                end
            else  % mouse click
                cp = get(gca, 'CurrentPoint');
                px = round(cp(1,1)); py = round(cp(1,2));
                gj_click = floor((px - gridInfo.border - 1) / gridInfo.tileW) + 1;
                gi_click = floor((py - gridInfo.border - 1) / gridInfo.tileH) + 1;
                if gj_click >= 1 && gj_click <= gridInfo.nCols && ...
                   gi_click >= 1 && gi_click <= gridInfo.nRows
                    roi_idx = (gi_click - 1) * gridInfo.nCols + gj_click;
                    if roi_idx >= 1 && roi_idx <= N_roi
                        if decision(roi_idx) == -1
                            decision(roi_idx) = 1;
                        else
                            decision(roi_idx) = -1;
                        end
                        toss_fov = false;
                        needRedraw = true;
                    end
                end
            end
        end
    end

    if isvalid(fig), close(fig); end

    % Apply decisions
    if cancelled
        % Escape = cancel: treat everything as keep
        fprintf('  Review cancelled. All ROIs kept.\n');
        toss_fov = false;
    else
        for i = 1:N_roi
            if decision(i) == -1
                rejected(i) = true;
                if toss_fov
                    reject_reason(i) = "toss_fov";
                else
                    reject_reason(i) = "manual";
                end
            end
        end
        nRej = sum(decision == -1);
        fprintf('  Review: %d rejected, %d kept\n', nRej, N_roi - nRej);
    end
end


function draw_grid_view(fig, montageRGB, gridInfo, validROI, ...
        N_roi, patchSz, decision, toss_fov)
%DRAW_GRID_VIEW  Draw the montage grid with white X on rejected, click-to-toggle.
    fig.Color = [0.15 0.15 0.15];

    ax_grid = axes(fig, 'Units', 'normalized', ...
        'Position', [0.02 0.06 0.96 0.86]);
    image(ax_grid, montageRGB);
    axis(ax_grid, 'image'); axis(ax_grid, 'off');
    hold(ax_grid, 'on');

    for i = 1:N_roi
        rid_i = validROI(i);
        gi = gridInfo.tileRow(i);
        gj = gridInfo.tileCol(i);
        x_center = (gj - 1) * gridInfo.tileW + gridInfo.border + 1 + patchSz / 2;
        y_top    = (gi - 1) * gridInfo.tileH + gridInfo.border + 1;

        text(ax_grid, x_center, y_top + 2, sprintf('%d', rid_i), ...
            'Color', 'y', 'FontSize', 7, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

        if decision(i) == -1
            x0 = (gj - 1) * gridInfo.tileW + gridInfo.border + 0.5;
            y0 = (gi - 1) * gridInfo.tileH + gridInfo.border + 0.5;
            xw = patchSz;
            yh = patchSz;
            patch(ax_grid, [x0 x0+xw x0+xw x0], [y0 y0 y0+yh y0+yh], 'w', ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(ax_grid, [x0 x0+xw], [y0 y0+yh], 'w-', 'LineWidth', 2);
            plot(ax_grid, [x0 x0+xw], [y0+yh y0], 'w-', 'LineWidth', 2);
        end
    end

    if toss_fov
        [mH, mW, ~] = size(montageRGB);
        patch(ax_grid, [0.5 mW+0.5 mW+0.5 0.5], [0.5 0.5 mH+0.5 mH+0.5], 'r', ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end
    hold(ax_grid, 'off');

    nReject = sum(decision == -1);
    if toss_fov
        titleStr = sprintf('TOSS FOV - all %d excluded | L-click=toggle | Enter=save | Esc=back', N_roi);
        titleColor = [1.0 0.3 0.3];
    else
        titleStr = sprintf('Reject: %d  Keep: %d  |  L-click=toggle  |  Enter=save  |  Esc=back', ...
            nReject, N_roi - nReject);
        titleColor = [1 1 1];
    end
    title(ax_grid, titleStr, 'Color', titleColor, ...
        'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
end


function save_marked_montage(outPath, montageRGB, gridInfo, validROI, ...
        N_roi, patchSz, rejected, toss_fov)
%SAVE_MARKED_MONTAGE  Save montage PNG with rejection marks and ROI labels.
    fig_s = figure('Color', 'k', 'Visible', 'off', ...
        'Position', [50 50 900 900]);
    ax_s = axes(fig_s, 'Units', 'normalized', 'Position', [0.02 0.04 0.96 0.92]);
    image(ax_s, montageRGB);
    axis(ax_s, 'image'); axis(ax_s, 'off');
    hold(ax_s, 'on');

    % Yellow ROI labels + white X on rejected
    for i = 1:N_roi
        rid = validROI(i);
        gi = gridInfo.tileRow(i);
        gj = gridInfo.tileCol(i);
        x_center = (gj - 1) * gridInfo.tileW + gridInfo.border + 1 + patchSz / 2;
        y_top    = (gi - 1) * gridInfo.tileH + gridInfo.border + 1;

        text(ax_s, x_center, y_top + 2, sprintf('%d', rid), ...
            'Color', 'y', 'FontSize', 7, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

        if rejected(i)
            x0 = (gj - 1) * gridInfo.tileW + gridInfo.border + 0.5;
            y0 = (gi - 1) * gridInfo.tileH + gridInfo.border + 0.5;
            xw = patchSz;
            yh = patchSz;
            patch(ax_s, [x0 x0+xw x0+xw x0], [y0 y0 y0+yh y0+yh], 'w', ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(ax_s, [x0 x0+xw], [y0 y0+yh], 'w-', 'LineWidth', 2);
            plot(ax_s, [x0 x0+xw], [y0+yh y0], 'w-', 'LineWidth', 2);
        end
    end

    % Red wash if toss_fov
    if toss_fov
        [mH, mW, ~] = size(montageRGB);
        patch(ax_s, [0.5 mW+0.5 mW+0.5 0.5], [0.5 0.5 mH+0.5 mH+0.5], 'r', ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end
    hold(ax_s, 'off');

    exportgraphics(fig_s, outPath, 'Resolution', 200);
    fprintf('  Montage saved: %s\n', outPath);
    close(fig_s);
end


