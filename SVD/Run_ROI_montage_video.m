clear all; clc; close all;
%% Run_ROI_montage_video.m
%  Generate a montage AVI of all ROI patches from the MC TIFF.
%  Each ROI's local patch is tiled in a grid, playing at imaging FPS.
%
%  Two input modes:
%    masterFolder: batch — scan all subfolders for *_cpSAM_output.mat
%    folderPath:   single FOV
%
%  Outputs (per FOV, in Fmovie_perROI/):
%    all_ROIs_montage.mp4       — montage video
%    montage_ROI_labels.png     — reference image with ROI numbers
%
%  Skips FOVs where all_ROIs_montage.mp4 already exists.

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
masterFolder = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";      % multi-experiment: scans all subfolders
folderPath   = "";      % single FOV mode (ignored if masterFolder set)

patchSize_um = 40;      % patch size per ROI (microns), square

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'mod')));

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

    % Skip if montage already exists
    outAVI = fullfile(fovFolder, 'Fmovie_perROI', 'all_ROIs_montage.mp4');
    if isfile(outAVI)
        fprintf('  [SKIP] all_ROIs_montage.mp4 already exists.\n');
        nSkipped = nSkipped + 1;
        continue
    end

    try
        generate_montage(fovFolder, samPath, patchSize_um);
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
fprintf('\n==================== MONTAGE VIDEO COMPLETE ====================\n');
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
fprintf('================================================================\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function samHits = discover_cpSAM(masterFolder, folderPath)
    if masterFolder ~= ""
        samHits = dir(fullfile(masterFolder, '**', '*_cpSAM_output.mat'));
    elseif folderPath ~= ""
        samHits = dir(fullfile(folderPath, '*_cpSAM_output.mat'));
        if isempty(samHits)
            samHits = dir(fullfile(folderPath, '**', '*_cpSAM_output.mat'));
        end
    else
        error('Set either masterFolder or folderPath.');
    end
    assert(~isempty(samHits), 'No *_cpSAM_output.mat found.');
    fprintf('[discover] Found %d cpSAM file(s)\n', numel(samHits));
end


function generate_montage(folderPath, samPath, patchSize_um)
%GENERATE_MONTAGE  Generate montage AVI for a single FOV.
    %% Detect FPS + pixel size
    [fps, scan_meta] = detect_session_fps(folderPath);
    fprintf('  FPS = %d Hz\n', fps);

    assert(isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um), ...
        'detect_session_fps did not return pixelSize_um.');
    pixelSize_um = scan_meta.pixelSize_um;
    patchRadius_px = round((patchSize_um / 2) / pixelSize_um);
    patchSz = 2 * patchRadius_px + 1;

    %% Load cpSAM
    samData = load(samPath, 'maskL');
    assert(isfield(samData, 'maskL'), 'SAM file missing maskL.');
    maskL = samData.maskL;

    props = regionprops(maskL, 'Centroid');
    N_roi = numel(props);
    fprintf('  ROIs: %d\n', N_roi);
    if N_roi == 0
        fprintf('  No ROIs. Skipping.\n');
        return
    end

    %% Discover + load MC TIFF
    tifPath = discover_mc_tiff(folderPath);
    fprintf('  MC TIFF: %s\n', tifPath);
    fprintf('  Loading TIFF stack...\n');
    tic;
    tifInfo = imfinfo(tifPath);
    H = tifInfo(1).Height;
    W = tifInfo(1).Width;
    T = numel(tifInfo);

    Y_full = zeros(H, W, T, 'single');
    for kk = 1:T
        Y_full(:,:,kk) = single(imread(tifPath, kk, 'Info', tifInfo));
    end
    fprintf('  Loaded %d frames [%d x %d] in %.1f s\n', T, H, W, toc);

    %% Patch bounds
    roi_patches = struct('r1',{}, 'r2',{}, 'c1',{}, 'c2',{});
    for i = 1:N_roi
        cx = props(i).Centroid(1);
        cy = props(i).Centroid(2);
        roi_patches(i).r1 = max(1,  round(cy) - patchRadius_px);
        roi_patches(i).r2 = min(H,  round(cy) + patchRadius_px);
        roi_patches(i).c1 = max(1,  round(cx) - patchRadius_px);
        roi_patches(i).c2 = min(W,  round(cx) + patchRadius_px);
    end

    %% Grid layout
    nCols_grid = ceil(sqrt(N_roi));
    nRows_grid = ceil(N_roi / nCols_grid);
    border = 2;
    tileW = patchSz + border;
    tileH = patchSz + border;
    montageW = nCols_grid * tileW + border;
    montageH = nRows_grid * tileH + border;

    %% Per-ROI brightness limits (all frames)
    roi_lo = zeros(N_roi, 1);
    roi_hi = zeros(N_roi, 1);
    for i = 1:N_roi
        p = roi_patches(i);
        patch_all = Y_full(p.r1:p.r2, p.c1:p.c2, :);
        roi_lo(i) = double(prctile(patch_all(:), 0.5));
        roi_hi(i) = double(prctile(patch_all(:), 99.5));
        if roi_hi(i) <= roi_lo(i), roi_hi(i) = roi_lo(i) + 1; end
    end

    %% Bitmap font (5x3 per digit)
    digit_font = zeros(5, 3, 10);
    digit_font(:,:,1)  = [1 1 1; 1 0 1; 1 0 1; 1 0 1; 1 1 1]; % 0
    digit_font(:,:,2)  = [0 1 0; 1 1 0; 0 1 0; 0 1 0; 1 1 1]; % 1
    digit_font(:,:,3)  = [1 1 1; 0 0 1; 1 1 1; 1 0 0; 1 1 1]; % 2
    digit_font(:,:,4)  = [1 1 1; 0 0 1; 1 1 1; 0 0 1; 1 1 1]; % 3
    digit_font(:,:,5)  = [1 0 1; 1 0 1; 1 1 1; 0 0 1; 0 0 1]; % 4
    digit_font(:,:,6)  = [1 1 1; 1 0 0; 1 1 1; 0 0 1; 1 1 1]; % 5
    digit_font(:,:,7)  = [1 1 1; 1 0 0; 1 1 1; 1 0 1; 1 1 1]; % 6
    digit_font(:,:,8)  = [1 1 1; 0 0 1; 0 1 0; 0 1 0; 0 1 0]; % 7
    digit_font(:,:,9)  = [1 1 1; 1 0 1; 1 1 1; 1 0 1; 1 1 1]; % 8
    digit_font(:,:,10) = [1 1 1; 1 0 1; 1 1 1; 0 0 1; 1 1 1]; % 9

    label_masks = cell(N_roi, 1);
    for i = 1:N_roi
        digits = num2str(i) - '0';
        glyph = [];
        for d = 1:numel(digits)
            if d > 1, glyph = [glyph, zeros(5, 1)]; end %#ok<AGROW>
            glyph = [glyph, digit_font(:,:, digits(d) + 1)]; %#ok<AGROW>
        end
        label_masks{i} = logical(glyph);
    end

    %% ROI boundary masks
    roi_perim_cell = cell(N_roi, 1);
    for i = 1:N_roi
        p = roi_patches(i);
        mask_local = (maskL(p.r1:p.r2, p.c1:p.c2) == i);
        perim_local = bwperim(mask_local, 4);
        perim_tile = false(patchSz, patchSz);
        pH = p.r2 - p.r1 + 1;
        pW = p.c2 - p.c1 + 1;
        perim_tile(1:pH, 1:pW) = perim_local;
        roi_perim_cell{i} = perim_tile;
    end

    %% Write montage AVI
    movieDir = fullfile(folderPath, 'Fmovie_perROI');
    if ~exist(movieDir, 'dir'), mkdir(movieDir); end

    outAVI = fullfile(movieDir, 'all_ROIs_montage.mp4');
    vw = VideoWriter(outAVI, 'MPEG-4');
    vw.FrameRate = fps;
    vw.Quality = 95;
    open(vw);

    fprintf('  Writing montage (%d x %d, %d frames @ %d fps)...\n', ...
        montageW, montageH, T, fps);
    tic;
    for f = 1:T
        montage_frame = zeros(montageH, montageW, 3, 'uint8');

        for i = 1:N_roi
            p = roi_patches(i);
            ptch = double(Y_full(p.r1:p.r2, p.c1:p.c2, f));

            ptch = (ptch - roi_lo(i)) / (roi_hi(i) - roi_lo(i));
            ptch = uint8(255 * max(0, min(1, ptch)));

            tile = zeros(patchSz, patchSz, 'uint8');
            pH = p.r2 - p.r1 + 1;
            pW = p.c2 - p.c1 + 1;
            tile(1:pH, 1:pW) = ptch;

            gi = ceil(i / nCols_grid);
            gj = mod(i - 1, nCols_grid) + 1;
            y0 = (gi - 1) * tileH + border + 1;
            x0 = (gj - 1) * tileW + border + 1;

            for ch = 1:3
                montage_frame(y0:y0+patchSz-1, x0:x0+patchSz-1, ch) = tile;
            end

            % Yellow boundary
            perim = roi_perim_cell{i};
            yr = y0:y0+patchSz-1;
            xr = x0:x0+patchSz-1;
            montage_frame(yr, xr, 1) = montage_frame(yr, xr, 1) .* uint8(~perim) + uint8(perim)*255;
            montage_frame(yr, xr, 2) = montage_frame(yr, xr, 2) .* uint8(~perim) + uint8(perim)*255;
            montage_frame(yr, xr, 3) = montage_frame(yr, xr, 3) .* uint8(~perim);

            % Bitmap label (yellow)
            lbl = label_masks{i};
            lh = size(lbl, 1);
            lw = size(lbl, 2);
            if lh <= patchSz && lw <= patchSz
                for lr = 1:lh
                    for lc = 1:lw
                        if lbl(lr, lc)
                            montage_frame(y0+lr-1, x0+lc-1, 1) = 255;
                            montage_frame(y0+lr-1, x0+lc-1, 2) = 255;
                            montage_frame(y0+lr-1, x0+lc-1, 3) = 0;
                        end
                    end
                end
            end
        end

        writeVideo(vw, montage_frame);
    end
    close(vw);
    fprintf('  Saved: %s (%.1f s)\n', outAVI, toc);

    % Reference label image — average projection across all frames
    avg_montage = zeros(montageH, montageW, 3, 'uint8');
    for i = 1:N_roi
        p = roi_patches(i);
        ptch = mean(double(Y_full(p.r1:p.r2, p.c1:p.c2, :)), 3);
        ptch = (ptch - roi_lo(i)) / (roi_hi(i) - roi_lo(i));
        ptch = uint8(255 * max(0, min(1, ptch)));

        tile = zeros(patchSz, patchSz, 'uint8');
        pH = p.r2 - p.r1 + 1;
        pW = p.c2 - p.c1 + 1;
        tile(1:pH, 1:pW) = ptch;

        gi = ceil(i / nCols_grid);
        gj = mod(i - 1, nCols_grid) + 1;
        y0 = (gi - 1) * tileH + border + 1;
        x0 = (gj - 1) * tileW + border + 1;

        for ch = 1:3
            avg_montage(y0:y0+patchSz-1, x0:x0+patchSz-1, ch) = tile;
        end

        % Yellow boundary
        perim = roi_perim_cell{i};
        yr = y0:y0+patchSz-1;
        xr = x0:x0+patchSz-1;
        avg_montage(yr, xr, 1) = avg_montage(yr, xr, 1) .* uint8(~perim) + uint8(perim)*255;
        avg_montage(yr, xr, 2) = avg_montage(yr, xr, 2) .* uint8(~perim) + uint8(perim)*255;
        avg_montage(yr, xr, 3) = avg_montage(yr, xr, 3) .* uint8(~perim);

        % Bitmap label (yellow)
        lbl = label_masks{i};
        lh = size(lbl, 1);
        lw = size(lbl, 2);
        if lh <= patchSz && lw <= patchSz
            for lr = 1:lh
                for lc = 1:lw
                    if lbl(lr, lc)
                        avg_montage(y0+lr-1, x0+lc-1, 1) = 255;
                        avg_montage(y0+lr-1, x0+lc-1, 2) = 255;
                        avg_montage(y0+lr-1, x0+lc-1, 3) = 0;
                    end
                end
            end
        end
    end
    labelPng = fullfile(movieDir, 'montage_ROI_labels.png');
    imwrite(avg_montage, labelPng);
    fprintf('  Labels (avg proj): %s\n', labelPng);
end


function tifPath = discover_mc_tiff(fovFolder)
    mc_hits = dir(fullfile(fovFolder, '*_MC*.tif'));
    assert(~isempty(mc_hits), 'No *_MC*.tif found in %s', fovFolder);
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
    nMC = zeros(numel(mc_hits), 1);
    for ii = 1:numel(mc_hits)
        nMC(ii) = numel(strfind(upper(mc_hits(ii).name), '_MC'));
    end
    [~, bestIdx] = max(nMC);
    tifPath = fullfile(mc_hits(bestIdx).folder, mc_hits(bestIdx).name);
end
