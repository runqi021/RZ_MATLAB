%% roi_pair_morph_quant.m
% Quantify cross-FOV ROI pairs for cell registration: for every pair of ROIs
% from DIFFERENT FOVs that are spatially near each other, compute
%   (1) spatial distance in shared stage coords (um), and
%   (2) morphology / appearance similarity,
% then scatter them so same-cell clusters separate from chance neighbors.
%
% Within-FOV pairs are NOT considered (Cellpose already gives distinct cells).
%
% Appearance similarity is OFFSET-TOLERANT: peak normalized cross-correlation
% (normxcorr2) of the two AVG patches over a search window, on a common um
% grid (FOVs have different zoom). This stays high for a true same-cell pair
% even when motor-position error leaves the centroids ~20 um apart.
%
% Ground-truth same-cell pairs (read off the map) are overlaid to validate.
%
% Standalone, plot-only. Reuses the loading logic of roi_centroid_map.m.
% Runqi Zhang / 2026

clear; close all;

%% ---------------- USER PARAMS ----------------
datasetPath  = "D:\260721_Sert_soma_G8s\phys";
xSign        = +1;
ySign        = +1;
fallback_fps = 30;
PixelSizeBase = 1.7778;

maxPairDist_um = 150;     % only score cross-FOV pairs within this XY distance
maxPairDz_um   = 40;     % ...and within this stage-Z difference
patchHalf_um   = 11;      % template half-size (um) around centroid
searchHalf_um  = 12;     % search half-size (um) on the other ROI (>= patchHalf + max offset)
targetUmPerPx  = 0.5;    % common resampling grid for patches
imgClip        = [5 99.5];

%% ---------------- PATH SETUP ----------------
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

%% ---------------- LOAD ROIs + per-FOV AVG images ----------------
sam_hits = dir(fullfile(datasetPath, '**', '*_cpSAM_output.mat'));
assert(~isempty(sam_hits), 'No *_cpSAM_output.mat found under %s', datasetPath);
nFov = numel(sam_hits);

% per-ROI table
R = struct('x',{},'y',{},'z',{},'fov',{},'roi',{},'lab',{},'area',{},'ecc',{},'sol',{}, ...
           'cx_px',{},'cy_px',{});
fov_img   = cell(nFov,1);   % normalized AVG image (double)
fov_maskL = cell(nFov,1);   % labeled mask (uint) for ROI-outline overlay
fov_pxum  = nan(nFov,1);
fov_name  = strings(nFov,1);

for ff = 1:nFov
    folderPath = sam_hits(ff).folder;
    samPath    = fullfile(folderPath, sam_hits(ff).name);
    [~, fov_name(ff)] = fileparts(folderPath);

    S = load(samPath, 'maskL');
    if ~isfield(S,'maskL') || isempty(S.maskL), continue; end
    maskL = S.maskL; [H, W] = size(maskL);
    labels = setdiff(unique(maskL(:)), 0);
    if isempty(labels), continue; end
    fov_maskL{ff} = maskL;

    props = regionprops(maskL, 'Centroid','Area','Eccentricity','Solidity');
    props = props(labels);

    [~, sm] = detect_session_fps(folderPath, fallback_fps);
    if ~isfield(sm,'motorPosition') || numel(sm.motorPosition) < 3, continue; end
    motor = sm.motorPosition(:).';
    px_um = NaN;
    if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um>0
        px_um = sm.pixelSize_um;
    elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor>0
        px_um = PixelSizeBase / sm.zoomFactor;
    end
    if ~isfinite(px_um), continue; end
    fov_pxum(ff) = px_um;

    avgIm = load_avg_image(samPath, folderPath);
    if ~isempty(avgIm)
        avgIm = double(avgIm);
        lo = prctile(avgIm(:), imgClip(1)); hi = prctile(avgIm(:), imgClip(2));
        if hi <= lo, hi = lo + 1; end
        fov_img{ff} = min(max((avgIm - lo)/(hi - lo), 0), 1);  % normalized [0,1]
    end

    for k = 1:numel(props)
        c = props(k).Centroid;     % [x_col, y_row]
        R(end+1) = struct( ...
            'x',   motor(1) + xSign*(c(1)-W/2)*px_um, ...
            'y',   motor(2) + ySign*(c(2)-H/2)*px_um, ...
            'z',   motor(3), ...
            'fov', ff, ...
            'roi', k, ...                       % within-FOV ROI index (F column)
            'lab', labels(k), ...               % label value in maskL
            'area',props(k).Area * px_um^2, ...
            'ecc', props(k).Eccentricity, ...
            'sol', props(k).Solidity, ...
            'cx_px', c(1), 'cy_px', c(2)); %#ok<SAGROW>
    end
end
assert(~isempty(R), 'No ROIs loaded.');
xs=[R.x].'; ys=[R.y].'; zs=[R.z].'; fov=[R.fov].'; area=[R.area].';
roiIdx=[R.roi].'; labArr=[R.lab].'; ecc=[R.ecc].'; sol=[R.sol].'; cxpx=[R.cx_px].'; cypx=[R.cy_px].';
nROI = numel(R);
fprintf('Loaded %d ROIs across %d FOVs.\n', nROI, nnz(~isnan(fov_pxum)));

%% ---------------- ENUMERATE CROSS-FOV CANDIDATE PAIRS ----------------
% NOTE: distance is used ONLY as a cheap prefilter so we don't run normxcorr2
% on all ~N^2 pairs. It is NOT a scoring feature (per request). To ignore
% distance entirely, set maxPairDist_um = Inf (slower: scores every pair).
% empty [] means "no gate" (search ALL cross-FOV pairs)
if isempty(maxPairDist_um), maxPairDist_um = Inf; end
if isempty(maxPairDz_um),   maxPairDz_um   = Inf; end

D = squareform(pdist([xs ys]));               % XY distance matrix (um)
[ia, ib] = find(triu(true(nROI), 1));         % all i<j
keep = (D(sub2ind([nROI nROI], ia, ib)) <= maxPairDist_um) ...
     & (fov(ia) ~= fov(ib)) ...
     & (abs(zs(ia) - zs(ib)) <= maxPairDz_um);
ia = ia(keep); ib = ib(keep);
nPair = numel(ia);
fprintf('Cross-FOV candidate pairs within %.0f um (Dz<=%.0f): %d\n', ...
    maxPairDist_um, maxPairDz_um, nPair);

pair_dist = D(sub2ind([nROI nROI], ia, ib));
pair_dz   = abs(zs(ia) - zs(ib));
pair_ncc  = nan(nPair,1);                     % offset-tolerant appearance sim
pair_iou  = nan(nPair,1);                     % mask IoU after NCC alignment

for p = 1:nPair
    i = ia(p); j = ib(p);
    % appearance NCC + the two ROI masks in the aligned frame -> IoU
    [pair_ncc(p), ~, ~, mA, mB] = pair_appearance( ...
        fov_img{fov(i)}, cxpx(i), cypx(i), fov_pxum(fov(i)), ...
        fov_img{fov(j)}, cxpx(j), cypx(j), fov_pxum(fov(j)), ...
        patchHalf_um, searchHalf_um, targetUmPerPx, ...
        fov_maskL{fov(i)}, labArr(i), fov_maskL{fov(j)}, labArr(j));
    if ~isempty(mA) && ~isempty(mB)
        uni = nnz(mA | mB);
        if uni > 0, pair_iou(p) = nnz(mA & mB) / uni; end
    end
end

%% ---------------- NCC vs IoU scatter (top 400 IoU black) ----------------
% Each point = one cross-FOV pair. x = appearance NCC, y = mask IoU.
% Top 400 by IoU shown black (the candidate same-cell matches).
nBlack = min(400, nPair);
[~, iouOrder] = sort(pair_iou, 'descend', 'MissingPlacement','last');
isTop = false(nPair,1); isTop(iouOrder(1:nBlack)) = true;

figure('Color','w','Position',[120 120 760 680]);
scatter(pair_ncc(~isTop), pair_iou(~isTop), 18, [.6 .6 .6], 'filled', 'MarkerFaceAlpha',0.4); hold on;
scatter(pair_ncc(isTop),  pair_iou(isTop),  30, 'k', 'filled');
xlabel('appearance: peak NCC'); ylabel('mask IoU (aligned)');
title(sprintf('IoU vs NCC — %d pairs (top %d IoU black)', nPair, nBlack));
legend({'rest', sprintf('top %d IoU', nBlack)}, 'Location','northwest');
grid on; box on;

%% ---------------- local helpers ----------------
function s = shortName(name)
% trim long FOV folder names for tables/titles
name = char(name);
if numel(name) > 28, s = [name(1:25) '...']; else, s = name; end
end

function im = load_avg_image(samPath, folderPath)
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
    if ~isempty(cand), im = imread(fullfile(cand(1).folder, cand(1).name)); end
end
if ~isempty(im) && size(im,3) > 1, im = im(:,:,1); end
end

function [ncc, tmpl, matched, tmplMask, matchedMask] = pair_appearance( ...
        imgA,cxA,cyA,pxA, imgB,cxB,cyB,pxB, halfT, halfS, target, maskA,labA, maskB,labB)
% Offset-tolerant appearance similarity between two ROIs.
%   tmpl    = small template patch around A (common um grid)
%   matched = the sub-crop of B's larger search patch at the NCC peak shift
%             -> tmpl and matched are EXACTLY what produced the score, same size.
%   tmplMask/matchedMask = each ROI's own binary mask on the SAME crop/shift
%             (only built when requested and masks are supplied).
ncc = NaN; tmpl = []; matched = []; tmplMask = []; matchedMask = [];
tmpl = extract_patch(imgA, cxA, cyA, halfT, pxA, target);
srch = extract_patch(imgB, cxB, cyB, halfS, pxB, target);
if isempty(tmpl) || isempty(srch) || any(size(srch) <= size(tmpl)), tmpl=[]; return; end
C = normxcorr2(tmpl - mean(tmpl(:)), srch - mean(srch(:)));
[ncc, idx] = max(C(:));
[yp, xp] = ind2sub(size(C), idx);
% normxcorr2 peak (yp,xp) = where template's lower-right corner aligns in 'valid' frame
r0 = yp - size(tmpl,1) + 1;  c0 = xp - size(tmpl,2) + 1;
r0 = min(max(r0,1), size(srch,1)-size(tmpl,1)+1);
c0 = min(max(c0,1), size(srch,2)-size(tmpl,2)+1);
matched = srch(r0:r0+size(tmpl,1)-1, c0:c0+size(tmpl,2)-1);

if nargout > 3 && nargin >= 15 && ~isempty(maskA) && ~isempty(maskB)
    tmplMask    = extract_patch(double(maskA==labA), cxA, cyA, halfT, pxA, target) > 0.5;
    srchMask    = extract_patch(double(maskB==labB), cxB, cyB, halfS, pxB, target) > 0.5;
    matchedMask = srchMask(r0:r0+size(tmpl,1)-1, c0:c0+size(tmpl,2)-1);
end
end

function draw_outline(mask, col)
% overlay binary-mask boundary on the current axes (image already drawn)
if isempty(mask) || ~any(mask(:)), return; end
hold on;
B = bwboundaries(mask);
for b = 1:numel(B)
    plot(B{b}(:,2), B{b}(:,1), '-', 'Color', col, 'LineWidth', 1);
end
hold off;
end

function patch = extract_patch(img, cx_px, cy_px, half_um, px_um, target_um_px)
% Crop a +/- half_um window around (cx,cy) and resample to target_um_px grid.
patch = [];
if isempty(img) || ~isfinite(px_um), return; end
rpx = ceil(half_um / px_um);
imgp = padarray(img, [rpx rpx], 'replicate');     % so crop is always valid
cc = round(cx_px) + rpx; rr = round(cy_px) + rpx;  % shifted center in padded img
crop = imgp(rr-rpx:rr+rpx, cc-rpx:cc+rpx);
n = round(2*half_um / target_um_px) + 1;           % output side in common grid
patch = imresize(crop, [n n]);
end
