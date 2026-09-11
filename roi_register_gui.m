%% roi_register_gui.m
% Cell registration across FOVs by mask-IoU, with clustering + manual proofing GUI.
%
% Pipeline:
%   1. Load every ROI (centroid in stage um, mask, AVG image, metadata).
%   2. For cross-FOV candidate pairs, align by appearance NCC and measure mask IoU.
%   3. Keep edges with IoU >= iouThresh (0.68) -> graph -> connected components.
%      Each component = one candidate cell (it may have been captured 2,3,4...x).
%   4. Flag components that contain >1 ROI from the SAME FOV (impossible for one
%      cell -> a chaining/over-merge error to fix by hand).
%   5. GUI: step through each multi-ROI cluster, click ROIs to include/exclude,
%      reject, or accept. Save -> roi_cell_registration.mat (per-ROI cell_id).
%
% Default = treat each auto-cluster as one cell; you only EDIT the exceptions.
%
% Standalone. Runqi Zhang / 2026

clear; close all;

%% ---------------- USER PARAMS ----------------
datasetPath   = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
iouThresh     = 0.68;    % IoU >= this -> candidate same-cell edge
xSign         = +1;
ySign         = +1;
fallback_fps  = 30;
PixelSizeBase = 1.7778;

maxPairDist_um = 200;    % spatial prefilter (tractability only, not a score)
maxPairDz_um   = 50;
patchHalf_um   = 9;      % template half-size (um) for NCC alignment
searchHalf_um  = 11;     % search half-size (um)
targetUmPerPx  = 0.5;    % common resampling grid
imgClip        = [5 99.5];
ctxHalf_um     = 13;     % half-size (um) of the patch shown per ROI in the GUI

outPath = fullfile(datasetPath, 'roi_cell_registration.mat');

%% ---------------- PATH SETUP ----------------
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);

%% ---------------- LOAD ROIs + per-FOV AVG images ----------------
sam_hits = dir(fullfile(datasetPath, '**', '*_cpSAM_output.mat'));
assert(~isempty(sam_hits), 'No *_cpSAM_output.mat found under %s', datasetPath);
nFov = numel(sam_hits);

R = struct('x',{},'y',{},'z',{},'fov',{},'roi',{},'lab',{},'cx_px',{},'cy_px',{});
fov_img   = cell(nFov,1);
fov_maskL = cell(nFov,1);
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

    props = regionprops(maskL, 'Centroid');
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
        fov_img{ff} = min(max((avgIm - lo)/(hi - lo), 0), 1);
    end

    for k = 1:numel(props)
        c = props(k).Centroid;
        R(end+1) = struct( ...
            'x', motor(1) + xSign*(c(1)-W/2)*px_um, ...
            'y', motor(2) + ySign*(c(2)-H/2)*px_um, ...
            'z', motor(3), 'fov', ff, 'roi', k, 'lab', labels(k), ...
            'cx_px', c(1), 'cy_px', c(2)); %#ok<SAGROW>
    end
end
assert(~isempty(R), 'No ROIs loaded.');
xs=[R.x].'; ys=[R.y].'; zs=[R.z].'; fov=[R.fov].';
roiIdx=[R.roi].'; labArr=[R.lab].'; cxpx=[R.cx_px].'; cypx=[R.cy_px].';
nROI = numel(R);
fprintf('Loaded %d ROIs across %d FOVs.\n', nROI, nnz(~isnan(fov_pxum)));

%% ---------------- IoU EDGES (cross-FOV, aligned) ----------------
D = squareform(pdist([xs ys]));
[ia, ib] = find(triu(true(nROI), 1));
keep = (D(sub2ind([nROI nROI], ia, ib)) <= maxPairDist_um) ...
     & (fov(ia) ~= fov(ib)) & (abs(zs(ia)-zs(ib)) <= maxPairDz_um);
ia = ia(keep); ib = ib(keep);
nPair = numel(ia);

eU = []; eV = []; eW = [];
for p = 1:nPair
    i = ia(p); j = ib(p);
    [~, ~, ~, mA, mB] = pair_appearance( ...
        fov_img{fov(i)}, cxpx(i), cypx(i), fov_pxum(fov(i)), ...
        fov_img{fov(j)}, cxpx(j), cypx(j), fov_pxum(fov(j)), ...
        patchHalf_um, searchHalf_um, targetUmPerPx, ...
        fov_maskL{fov(i)}, labArr(i), fov_maskL{fov(j)}, labArr(j));
    if isempty(mA) || isempty(mB), continue; end
    uni = nnz(mA | mB);
    if uni == 0, continue; end
    iou = nnz(mA & mB) / uni;
    if iou >= iouThresh
        eU(end+1) = i; eV(end+1) = j; eW(end+1) = iou; %#ok<SAGROW>
    end
end
fprintf('%d cross-FOV pairs scored; %d edges at IoU >= %.2f\n', nPair, numel(eW), iouThresh);

%% ---------------- CLUSTER (connected components) ----------------
if isempty(eU)
    comp = (1:nROI)';
else
    G = graph(eU, eV, eW, nROI);
    comp = conncomp(G)';
end
% group ROIs by component; keep only multi-ROI clusters for review
[uc, ~, ic] = unique(comp);
clusters = {};
for u = 1:numel(uc)
    mem = find(ic == u);
    if numel(mem) >= 2, clusters{end+1} = mem(:)'; end %#ok<SAGROW>
end
nClust = numel(clusters);
% conflict flag: >1 ROI from same FOV within a cluster
conflict = false(1, nClust);
for c = 1:nClust
    fv = fov(clusters{c});
    conflict(c) = numel(fv) ~= numel(unique(fv));
end
fprintf('%d multi-ROI clusters (%d with same-FOV conflicts to check).\n', ...
    nClust, nnz(conflict));

%% ---------------- LAUNCH GUI ----------------
S = struct();
S.clusters = clusters; S.conflict = conflict; S.cur = 1; S.nClust = nClust;
S.group = cellfun(@(m) ones(size(m)), clusters, 'uni', 0);  % default: all -> group 1 (one cell)
S.fov=fov; S.roiIdx=roiIdx; S.fov_name=fov_name; S.labArr=labArr;
S.cxpx=cxpx; S.cypx=cypx; S.fov_img=fov_img; S.fov_maskL=fov_maskL; S.fov_pxum=fov_pxum;
S.ctxHalf_um=ctxHalf_um; S.targetUmPerPx=targetUmPerPx;
S.nROI=nROI; S.xs=xs; S.ys=ys; S.zs=zs; S.comp=comp; S.outPath=outPath;
S.iouThresh=iouThresh; S.R=R;

if nClust == 0
    fprintf('No multi-ROI clusters to review. Saving singletons only.\n');
    save_registration(S);
    return
end

fig = figure('Color','w','Name','ROI registration proofing','NumberTitle','off', ...
    'Position',[80 80 1300 800], 'KeyPressFcn',@(~,e) on_key(e));
guidata(fig, S);
redraw(fig);


%% ============================ GUI CALLBACKS ============================
function redraw(fig)
S = guidata(fig); clf(fig);
c   = S.cur; mem = S.clusters{c}; grp = S.group{c};
n   = numel(mem);

% --- status bar: conflict = a group (>0) with a repeated FOV ---
warnStr = '';
for v = unique(grp(grp>0))'
    fv = S.fov(mem(grp==v));
    if numel(fv) ~= numel(unique(fv))
        warnStr = '   *** SAME-FOV CONFLICT within a group ***'; break;
    end
end
nCells = numel(unique(grp(grp>0)));
uicontrol(fig,'Style','text','Units','normalized','Position',[0.02 0.93 0.96 0.05], ...
    'BackgroundColor','w','FontSize',11,'HorizontalAlignment','left', ...
    'String',sprintf('Cluster %d / %d   |   %d ROIs   ->   %d cell(s), %d out%s', ...
        c, S.nClust, n, nCells, nnz(grp==0), warnStr), ...
    'ForegroundColor', tern(~isempty(warnStr),[0.8 0 0],[0 0 0]));

% --- member patches, outline colored by group ---
nCols = min(6, n); nRows = ceil(n/nCols);
x0=0.03; y0=0.14; ww=0.94; hh=0.76;
for k = 1:n
    rr = floor((k-1)/nCols); cc = mod(k-1,nCols);
    ax = axes(fig,'Position',[x0+cc*ww/nCols, y0+(nRows-1-rr)*hh/nRows, ...
        0.92*ww/nCols, 0.86*hh/nRows]); %#ok<LAXES>
    [patch, mask] = roi_patch(S, mem(k));
    if ~isempty(patch), imagesc(ax, patch); caxis(ax,[0 1]); end
    axis(ax,'image','off'); colormap(ax, flipud(gray));
    hold(ax,'on');
    oc = group_color(grp(k));
    if ~isempty(mask)
        B = bwboundaries(mask);
        for b=1:numel(B), plot(ax,B{b}(:,2),B{b}(:,1),'-','Color',oc,'LineWidth',1.8); end
    end
    if grp(k)==0, gtag='out'; else, gtag=sprintf('g%d',grp(k)); end
    title(ax, sprintf('%s | F%d r%d', gtag, S.fov(mem(k)), S.roiIdx(mem(k))), ...
        'FontSize',8, 'Color', oc, 'Interpreter','none');
    set(ax,'ButtonDownFcn',@(~,~) cycle_group(fig,k));
    set(allchild(ax),'HitTest','off');
end

% --- buttons ---
btn = @(x,str,cb) uicontrol(fig,'Style','pushbutton','Units','normalized', ...
    'Position',[x 0.03 0.13 0.07],'String',str,'FontSize',10,'Callback',cb);
btn(0.03,'<< Prev',       @(~,~) nav(fig,-1));
btn(0.17,'All OUT',       @(~,~) set_all(fig,0));
btn(0.31,'Reset (1 cell)',@(~,~) set_all(fig,1));
btn(0.45,'Next >>',       @(~,~) nav(fig,+1));
btn(0.73,'SAVE & FINISH', @(~,~) finish(fig));
uicontrol(fig,'Style','text','Units','normalized','Position',[0.59 0.02 0.13 0.08], ...
    'BackgroundColor','w','FontSize',8, ...
    'String','click a cell: cycle group  g1->g2->...->out. same group = same cell.');
end

function cycle_group(fig,k)
S = guidata(fig); c = S.cur; maxG = numel(S.clusters{c});
g = S.group{c}(k) + 1; if g > maxG, g = 0; end   % wrap past last group -> out
S.group{c}(k) = g; guidata(fig,S); redraw(fig);
end

function set_all(fig,v)
S = guidata(fig); S.group{S.cur}(:) = v; guidata(fig,S); redraw(fig);
end

function nav(fig,d)
S = guidata(fig); S.cur = min(max(S.cur+d,1),S.nClust); guidata(fig,S); redraw(fig);
end

function on_key(e)
fig = gcf;
switch e.Key
    case 'rightarrow', nav(fig,+1);
    case 'leftarrow',  nav(fig,-1);
end
end

function finish(fig)
S = guidata(fig);
save_registration(S);
fprintf('Saved registration to %s\n', S.outPath);
close(fig);
end

function save_registration(S)
% Build per-ROI cell_id: start each ROI as its own cell, then merge the
% INCLUDED members of each reviewed cluster into one shared id.
cellId = (1:S.nROI)';
if isfield(S,'clusters')
    for c = 1:numel(S.clusters)
        mem = S.clusters{c}; grp = S.group{c};
        for v = unique(grp(grp>0))'          % each positive group -> its own cell
            gm = mem(grp==v);
            if numel(gm) >= 2, cellId(gm) = gm(1); end
        end
        % grp==0 ("out") members keep their own singleton id
    end
end
[~,~,cellId] = unique(cellId);   % relabel 1..K contiguous

reg = struct('fov',num2cell(S.fov), 'fov_name',num2cell(S.fov_name(S.fov)), ...
    'roi',num2cell(S.roiIdx), 'x',num2cell(S.xs), 'y',num2cell(S.ys), ...
    'z',num2cell(S.zs), 'cell_id',num2cell(cellId));
nCells = numel(unique(cellId));
counts = accumarray(cellId,1);
fprintf('Registration: %d ROIs -> %d cells (max captures = %d).\n', ...
    S.nROI, nCells, max(counts));
save(S.outPath, 'reg', 'cellId', 'iouThresh_used', '-v7'); %#ok<NASGU>
iouThresh_used = S.iouThresh; %#ok<NASGU>
save(S.outPath, 'reg', 'cellId', 'iouThresh_used');
end

%% ============================ helpers ============================
function [patch, mask] = roi_patch(S, idx)
patch = extract_patch(S.fov_img{S.fov(idx)}, S.cxpx(idx), S.cypx(idx), ...
    S.ctxHalf_um, S.fov_pxum(S.fov(idx)), S.targetUmPerPx);
mask  = extract_patch(double(S.fov_maskL{S.fov(idx)}==S.labArr(idx)), ...
    S.cxpx(idx), S.cypx(idx), S.ctxHalf_um, S.fov_pxum(S.fov(idx)), S.targetUmPerPx) > 0.5;
end

function out = tern(cond,a,b), if cond, out=a; else, out=b; end, end

function col = group_color(g)
% group 0 ("out") = gray; groups 1.. = distinct colors
if g == 0, col = [0.6 0.6 0.6]; return; end
cmap = lines(8); col = cmap(mod(g-1,8)+1, :);
end

function [ncc, tmpl, matched, tmplMask, matchedMask] = pair_appearance( ...
        imgA,cxA,cyA,pxA, imgB,cxB,cyB,pxB, halfT, halfS, target, maskA,labA, maskB,labB)
ncc = NaN; tmpl = []; matched = []; tmplMask = []; matchedMask = [];
tmpl = extract_patch(imgA, cxA, cyA, halfT, pxA, target);
srch = extract_patch(imgB, cxB, cyB, halfS, pxB, target);
if isempty(tmpl) || isempty(srch) || any(size(srch) <= size(tmpl)), tmpl=[]; return; end
C = normxcorr2(tmpl - mean(tmpl(:)), srch - mean(srch(:)));
[ncc, idx] = max(C(:));
[yp, xp] = ind2sub(size(C), idx);
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

function patch = extract_patch(img, cx_px, cy_px, half_um, px_um, target_um_px)
patch = [];
if isempty(img) || ~isfinite(px_um), return; end
rpx = ceil(half_um / px_um);
imgp = padarray(img, [rpx rpx], 'replicate');
cc = round(cx_px) + rpx; rr = round(cy_px) + rpx;
crop = imgp(rr-rpx:rr+rpx, cc-rpx:cc+rpx);
n = round(2*half_um / target_um_px) + 1;
patch = imresize(crop, [n n]);
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
