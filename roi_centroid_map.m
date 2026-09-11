%% roi_centroid_map.m
% Compute every ROI centroid across all FOVs in a dataset folder, transform
% to shared stage coordinates (um), and scatter them to eyeball clustering.
%
% Quick exploratory step toward "cell registration": a single physical cell
% can appear as multiple ROIs across repeated/overlapping FOV recordings.
% Plotting all centroids in one stage-coordinate frame shows whether ROIs
% from different FOVs land at the same physical location.
%
% Standalone, plot-only. Reads:
%   <fov>/*_cpSAM_output.mat  -> maskL (HxW labeled mask)
%   detect_session_fps()      -> motorPosition [x y z] um, pixelSize_um
%
% Runqi Zhang / 2026

clear; close all;

%% ---------------- USER PARAMS ----------------
datasetPath  = "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing";
xSign        = +1;     % flip if image-X is mirrored vs stage-X
ySign        = +1;     % flip if image-Y is mirrored vs stage-Y (verify w/ known nudge)
fallback_fps = 30;
PixelSizeBase = 1.7778; % um/px at zoom=1 (fallback if pixelSize_um missing)

% --- overlay options ---
showImages = true;     % overlay each FOV's AVG projection (tinted per FOV) in stage coords
showDots   = true;     % also draw centroid dots on top of the images
imgClip    = [2 99.5]; % percentile clip for per-image contrast normalization
imgGamma   = 0.8;      % <1 lifts dim pixels in the transparency mask

%% ---------------- PATH SETUP ----------------
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);   % so detect_session_fps resolves

%% ---------------- DISCOVER FOVs ----------------
sam_hits = dir(fullfile(datasetPath, '**', '*_cpSAM_output.mat'));
assert(~isempty(sam_hits), 'No *_cpSAM_output.mat found under %s', datasetPath);
nFov = numel(sam_hits);
fprintf('Found %d FOV(s) under %s\n\n', nFov, datasetPath);

% accumulators (one row per ROI)
all_x = []; all_y = []; all_z = []; all_area = []; all_fov = [];
fov_names = strings(nFov,1);
fov_z     = nan(nFov,1);
fov_keep  = false(nFov,1);
fov_extent = nan(nFov,4);   % [xc yc halfW halfH] in um for extent rectangles
fov_img    = cell(nFov,1);  % normalized AVG image [0..1] per FOV (for overlay)
fov_imgXY  = nan(nFov,4);   % [xL xR yTopRow yBotRow] image-edge coords in um
summary = strings(nFov,1);

%% ---------------- PER-FOV LOOP ----------------
for ff = 1:nFov
    folderPath = sam_hits(ff).folder;
    samPath    = fullfile(folderPath, sam_hits(ff).name);
    [~, fovName] = fileparts(folderPath);
    fov_names(ff) = string(fovName);

    % --- labeled mask -> centroids ---
    S = load(samPath, 'maskL');
    if ~isfield(S, 'maskL') || isempty(S.maskL)
        fprintf(2, '[skip] %s: no maskL\n', fovName); continue;
    end
    maskL = S.maskL;
    [H, W] = size(maskL);
    labels = setdiff(unique(maskL(:)), 0);
    if isempty(labels)
        fprintf(2, '[skip] %s: no ROIs in maskL\n', fovName); continue;
    end
    props = regionprops(maskL, 'Centroid', 'Area');
    props = props(labels);                       % align rows with F columns
    cents = reshape([props.Centroid], 2, []).';  % [N x 2] = [x_col, y_row]
    areas_px = [props.Area].';                    % [N x 1]

    % --- metadata (shared stage frame + per-FOV pixel size) ---
    [~, sm] = detect_session_fps(folderPath, fallback_fps);
    if ~isfield(sm,'motorPosition') || numel(sm.motorPosition) < 3
        fprintf(2, '[skip] %s: no motorPosition\n', fovName); continue;
    end
    motor = sm.motorPosition(:).';
    px_um = NaN;
    if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um > 0
        px_um = sm.pixelSize_um;
    elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor > 0
        px_um = PixelSizeBase / sm.zoomFactor;
    end
    if ~isfinite(px_um)
        fprintf(2, '[skip] %s: no pixelSize\n', fovName); continue;
    end

    % --- transform to stage coords (um) ---
    cx_um = motor(1) + xSign * (cents(:,1) - W/2) * px_um;
    cy_um = motor(2) + ySign * (cents(:,2) - H/2) * px_um;
    cz_um = motor(3) * ones(size(cx_um));
    area_um2 = areas_px * px_um^2;

    % --- accumulate ---
    all_x   = [all_x;   cx_um];      %#ok<AGROW>
    all_y   = [all_y;   cy_um];      %#ok<AGROW>
    all_z   = [all_z;   cz_um];      %#ok<AGROW>
    all_area= [all_area; area_um2];  %#ok<AGROW>
    all_fov = [all_fov; ff*ones(size(cx_um))]; %#ok<AGROW>

    fov_z(ff)     = motor(3);
    fov_keep(ff)  = true;
    fov_extent(ff,:) = [motor(1), motor(2), (W/2)*px_um, (H/2)*px_um];

    % --- AVG projection image, normalized, for stage-coord overlay ---
    if showImages
        avgIm = load_avg_image(samPath, folderPath);
        if ~isempty(avgIm)
            avgIm = double(avgIm);
            lo = prctile(avgIm(:), imgClip(1));
            hi = prctile(avgIm(:), imgClip(2));
            if hi <= lo, hi = lo + 1; end
            fov_img{ff} = min(max((avgIm - lo) / (hi - lo), 0), 1);
            % edge coords: col/row centers 1 and W/H mapped via the same transform
            fov_imgXY(ff,:) = [ motor(1) + xSign*(1 - W/2)*px_um, ...
                                motor(1) + xSign*(W - W/2)*px_um, ...
                                motor(2) + ySign*(1 - H/2)*px_um, ...
                                motor(2) + ySign*(H - H/2)*px_um ];
        end
    end

    zoomStr = '';
    if isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor), zoomStr = sprintf('%gx', sm.zoomFactor); end
    summary(ff) = sprintf('  %-48s N=%-4d  XYZ=[%.0f %.0f %.0f]  zoom=%-5s  px=%.4f um', ...
        fovName, numel(labels), motor(1), motor(2), motor(3), zoomStr, px_um);
end

%% ---------------- SUMMARY ----------------
fprintf('FOV summary:\n');
for ff = 1:nFov
    if strlength(summary(ff)) > 0, fprintf('%s\n', summary(ff)); end
end
fprintf('\nTotal ROIs plotted: %d across %d FOV(s)\n', numel(all_x), nnz(fov_keep));
assert(~isempty(all_x), 'No ROIs accumulated — nothing to plot.');

%% ---------------- VISUALIZE (2D XY, colored by FOV, Z annotated) ----------------
keepIdx = find(fov_keep);
cmap = lines(numel(keepIdx));

figure('Color','w','Position',[100 100 900 800]); hold on;
legH = gobjects(numel(keepIdx),1);
legTxt = strings(numel(keepIdx),1);

% Pass 1: tinted AVG images (drawn first so dots/rectangles sit on top).
% Overlap regions blend colors -> candidate same-cell across FOVs.
if showImages
    for ii = 1:numel(keepIdx)
        ff = keepIdx(ii);
        if isempty(fov_img{ff}), continue; end
        I = fov_img{ff};
        tint = reshape(cmap(ii,:), 1, 1, 3);
        rgb  = I .* tint;                 % grayscale modulates the FOV hue
        alpha = I .^ imgGamma;            % dark background -> transparent
        e = fov_imgXY(ff,:);
        image('XData', e(1:2), 'YData', e(3:4), 'CData', rgb, ...
              'AlphaData', alpha, 'AlphaDataMapping', 'none');
    end
end

% Pass 2: extent rectangles + centroid dots
for ii = 1:numel(keepIdx)
    ff = keepIdx(ii);
    sel = (all_fov == ff);
    % faint FOV extent rectangle (shows spatial overlap between FOVs)
    e = fov_extent(ff,:);
    rectangle('Position', [e(1)-e(3), e(2)-e(4), 2*e(3), 2*e(4)], ...
        'EdgeColor', [cmap(ii,:) 0.5], 'LineStyle','-', 'LineWidth', 0.75);
    % ROI centroids (legend handle taken from here either way)
    if showDots
        legH(ii) = scatter(all_x(sel), all_y(sel), 36, cmap(ii,:), 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor','k', 'LineWidth',0.25);
    else
        legH(ii) = plot(NaN, NaN, 's', 'MarkerFaceColor', cmap(ii,:), ...
            'MarkerEdgeColor','none', 'MarkerSize', 8);
    end
    legTxt(ii) = sprintf('%s  (z=%g, n=%d)', shortName(fov_names(ff)), fov_z(ff), nnz(sel));
end

axis equal; grid on; box on;
xlabel('stage X (\mum)'); ylabel('stage Y (\mum)');
[~, dsName] = fileparts(datasetPath);
title(sprintf('ROI centroids in stage coords — %s  (%d ROIs, %d FOVs)', ...
    dsName, numel(all_x), numel(keepIdx)), 'Interpreter','none');
legend(legH, legTxt, 'Location','eastoutside', 'Interpreter','none', 'FontSize',8);
set(gca,'YDir','normal');

%% ---------------- local helpers ----------------
function s = shortName(name)
% trim long FOV folder names for the legend
name = char(name);
if numel(name) > 34, s = [name(1:31) '...']; else, s = name; end
end

function im = load_avg_image(samPath, folderPath)
% Return the AVG projection for this FOV as a 2D array (or [] if none found).
% Priority: avgPath stored in cpSAM_output.mat -> *_AVG_for_CP.tif -> *_AVG*.tif
im = [];
try
    S = load(samPath, 'avgPath');
    if isfield(S,'avgPath') && ~isempty(S.avgPath) && isfile(S.avgPath)
        im = imread(char(S.avgPath));
    end
catch
end
if isempty(im)
    cand = dir(fullfile(folderPath, '*_AVG_for_CP.tif'));
    if isempty(cand), cand = dir(fullfile(folderPath, '*_AVG*.tif')); end
    if ~isempty(cand)
        im = imread(fullfile(cand(1).folder, cand(1).name));
    end
end
if ~isempty(im) && size(im,3) > 1, im = im(:,:,1); end  % first channel only
end
