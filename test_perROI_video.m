%% test_perROI_video.m
%  Generate a per-ROI MP4 combining:
%    Top:    cropped calcium movie patch (grayscale+yellow ROI perimeter)
%    Bottom: F trace (white on black) with a moving yellow time-indicator line
%
%  Test script — hardcoded folder and ROI index.

clear; clc; close all;

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
folderPath = 'D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing\roi5_1400-1200-20_20lp_930_x4.4_512x512_6000f_00001';
roiList    = [3 7];       % empty = all ROIs; or e.g. [1 3 5] for specific ROIs
patchHalf_um = 20;        % half-width of crop window (microns)
gamma      = 0.6;         % display gamma for movie patch
videoQuality = 95;        % MPEG-4 quality (0-100)

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

%% ======================================================================
%  1) DETECT FPS & PIXEL SIZE
%  ======================================================================
[fps, scan_meta] = detect_session_fps(folderPath);
fprintf('FPS = %d Hz\n', fps);

assert(isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um), ...
    'detect_session_fps did not return a valid pixelSize_um.');
pixelSize_um = scan_meta.pixelSize_um;
patchHalf_px = round(patchHalf_um / pixelSize_um);
fprintf('Pixel size = %.4f um  ->  patch half-width = %d px\n', pixelSize_um, patchHalf_px);

%% ======================================================================
%  2) LOAD cpSAM OUTPUT (maskL + F)
%  ======================================================================
samHits = dir(fullfile(folderPath, '*_cpSAM_output.mat'));
assert(~isempty(samHits), 'No *_cpSAM_output.mat found in:\n  %s', folderPath);
samPath = fullfile(samHits(1).folder, samHits(1).name);
fprintf('Loading cpSAM: %s\n', samHits(1).name);

samData = load(samPath, 'maskL', 'F');
assert(isfield(samData, 'maskL'), 'cpSAM file missing maskL.');
assert(isfield(samData, 'F'),     'cpSAM file missing F.');
maskL = samData.maskL;
F_all = samData.F;   % [T x N]

[T_F, N_roi] = size(F_all);
fprintf('F: [%d x %d]  (T x N)\n', T_F, N_roi);

if isempty(roiList), roiList = 1:N_roi; end

%% ======================================================================
%  3) COMPUTE dF/F
%  ======================================================================
fprintf('Computing dF/F with helper.dFF_RZ ...\n');
outDFF = helper.dFF_RZ(F_all, 'FPS', fps);
dFF = outDFF.dFF;          % [T2 x N]
t_dff = outDFF.t_dff;      % [T2 x 1]
T_dff = size(dFF, 1);

%% ======================================================================
%  4) DISCOVER & LOAD MC TIFF
%  ======================================================================
tifPath = discover_mc_tiff(folderPath);
fprintf('MC TIFF: %s\n', tifPath);

fprintf('Loading TIFF stack...\n');
tic;
tifInfo = imfinfo(tifPath);
H = tifInfo(1).Height;
W = tifInfo(1).Width;
T_tif = numel(tifInfo);
fprintf('  TIFF: [%d x %d x %d]\n', H, W, T_tif);

% Sanity: TIFF frames must match F rows
assert(T_tif == T_F, ...
    'Frame count mismatch: TIFF has %d frames but cpSAM F has %d rows.', T_tif, T_F);

% Load full stack (single precision)
Y = zeros(H, W, T_tif, 'single');
for kk = 1:T_tif
    Y(:,:,kk) = single(imread(tifPath, kk, 'Info', tifInfo));
end
fprintf('  Loaded in %.1f s\n', toc);

%% ======================================================================
%  5) LOOP OVER ROIs
%  ======================================================================
props = regionprops(maskL, 'Centroid');
movieDir = fullfile(folderPath, 'Fmovie_perROI');
if ~exist(movieDir, 'dir'), mkdir(movieDir); end

fprintf('\n=== Rendering %d ROI videos ===\n', numel(roiList));
for ri = 1:numel(roiList)
roiIdx = roiList(ri);
fprintf('\n--- ROI %d (%d/%d) ---\n', roiIdx, ri, numel(roiList));

cx = round(props(roiIdx).Centroid(1));
cy = round(props(roiIdx).Centroid(2));

r1 = max(1,  cy - patchHalf_px);
r2 = min(H,  cy + patchHalf_px);
c1 = max(1,  cx - patchHalf_px);
c2 = min(W,  cx + patchHalf_px);
pH = r2 - r1 + 1;
pW = c2 - c1 + 1;

fprintf('ROI %d centroid: (%d, %d)  patch: rows [%d:%d] cols [%d:%d] = %dx%d px\n', ...
    roiIdx, cx, cy, r1, r2, c1, c2, pH, pW);

% Extract all patch frames
patch_all = Y(r1:r2, c1:c2, :);

% Percentile clipping limits (across all frames of this ROI patch)
lo = double(prctile(patch_all(:), 0.5));
hi = double(prctile(patch_all(:), 99.5));
if hi <= lo, hi = lo + 1; end
fprintf('  Clip limits: [%.1f, %.1f]\n', lo, hi);

% ROI perimeter mask (local patch coordinates)
mask_local = (maskL(r1:r2, c1:c2) == roiIdx);
perim_local = bwperim(mask_local, 4);

% F trace for this ROI
F_roi = F_all(:, roiIdx);
% Upscale movie patch to ~200px for visibility
movieSz = 200;  % target size for the movie panel (square)
scaleFactor = movieSz / max(pH, pW);

% Upscale perimeter mask too
perim_up = imresize(double(perim_local), scaleFactor, 'nearest') > 0.5;
upH = size(perim_up, 1);
upW = size(perim_up, 2);

% Trace panel: same height as movie, 4x wider
traceH = upH;
traceW = upW * 4;

% Render F trace as a white line on black background image (thick line)
traceBG = zeros(traceH, traceW, 3, 'uint8');

F_min = min(F_roi);
F_max = max(F_roi);
if F_max == F_min, F_max = F_min + 1; end

yMargin = round(traceH * 0.08);
yLo = traceH - yMargin;
yHi = yMargin + 1;

x_px = linspace(1, traceW, T_tif);
y_px = yLo - (F_roi - F_min) / (F_max - F_min) * (yLo - yHi);
y_px = round(y_px);
y_px = max(1, min(traceH, y_px));

% Draw thick white trace (±1 pixel width)
lineThick = 1;  % half-thickness
for tt = 1:T_tif-1
    x0 = round(x_px(tt));
    x1 = round(x_px(tt+1));
    y0 = y_px(tt);
    y1 = y_px(tt+1);
    nSteps = max(abs(x1-x0), abs(y1-y0));
    if nSteps == 0, nSteps = 1; end
    for s = 0:nSteps
        xi = round(x0 + (x1-x0) * s / nSteps);
        yi = round(y0 + (y1-y0) * s / nSteps);
        for dy = -lineThick:lineThick
            yy = max(1, min(traceH, yi + dy));
            xx = max(1, min(traceW, xi));
            traceBG(yy, xx, :) = 255;
        end
    end
end

fprintf('Movie panel: [%d x %d]  Trace panel: [%d x %d]\n', upH, upW, traceH, traceW);

%% ======================================================================
%  7) WRITE VIDEO
%  ======================================================================
outMP4 = fullfile(movieDir, sprintf('ROI_%02d_movie.mp4', roiIdx));
vw = VideoWriter(outMP4, 'MPEG-4');
vw.FrameRate = fps;
vw.Quality   = videoQuality;
open(vw);

% Layout: [movie | gap | trace]  side by side
gapW = 4;
totalW = upW + gapW + traceW;
totalH = upH;

% H.264 requires even dimensions
if mod(totalW, 2) ~= 0, totalW = totalW + 1; end
if mod(totalH, 2) ~= 0, totalH = totalH + 1; end

fprintf('Writing video: [%d x %d], %d frames @ %d fps\n', totalW, totalH, T_tif, fps);
fprintf('Output: %s\n', outMP4);
tic;

for ff = 1:T_tif
    % --- Left panel: upscaled movie patch ---
    ptch = double(patch_all(:,:,ff));
    ptch = (ptch - lo) / (hi - lo);
    ptch = max(0, min(1, ptch));
    ptch = ptch .^ gamma;

    % Upscale
    ptch_up = imresize(ptch, [upH upW], 'bilinear');
    ptch_u8 = uint8(255 * ptch_up);
    movieRGB = repmat(ptch_u8, [1 1 3]);

    % Burn in yellow ROI perimeter
    movieRGB(:,:,1) = movieRGB(:,:,1) .* uint8(~perim_up) + uint8(perim_up) * 255;
    movieRGB(:,:,2) = movieRGB(:,:,2) .* uint8(~perim_up) + uint8(perim_up) * 255;
    movieRGB(:,:,3) = movieRGB(:,:,3) .* uint8(~perim_up);

    % --- Right panel: F trace + time indicator ---
    traceFrame = traceBG;

    % Yellow vertical line (2px wide)
    xLine = round(x_px(ff));
    xLine = max(2, min(traceW-1, xLine));
    traceFrame(:, xLine-1:xLine+1, 1) = 255;
    traceFrame(:, xLine-1:xLine+1, 2) = 255;
    traceFrame(:, xLine-1:xLine+1, 3) = 0;

    % --- Assemble full frame ---
    fullFrame = zeros(totalH, totalW, 3, 'uint8');
    fullFrame(1:upH, 1:upW, :) = movieRGB;
    fullFrame(1:traceH, upW+gapW+1:upW+gapW+traceW, :) = traceFrame;

    writeVideo(vw, fullFrame);

    if mod(ff, 500) == 0
        fprintf('  frame %d / %d\n', ff, T_tif);
    end
end

close(vw);
elapsed = toc;
fprintf('Done! %.1f s (%.1f fps encoding rate)\n', elapsed, T_tif/elapsed);
fprintf('Saved: %s\n', outMP4);

end  % for roiIdx
fprintf('\n=== All %d ROI videos complete ===\n', numel(roiList));

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function tifPath = discover_mc_tiff(fovFolder)
%DISCOVER_MC_TIFF  Find the MC TIFF with the most "_MC" occurrences.
%  Excludes _AVG, _QC, _shifts, _SELROI files and non-.tif files.
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
