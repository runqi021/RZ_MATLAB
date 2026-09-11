function roi_ncc_round1_toss_260801()
%% ROUND 1 of 3 - TOSS bad observations from the auto groups.
% Click a patch to toss / un-toss it. Tossing writes NA (-1) into the shared
% assignment vector grpOf in roi_match_curated.mat -- the data is NOT deleted,
% just excluded from downstream analysis. Patches are cached once, so it's fast.
%
% CONTROLS: left-click tile = toss/un-toss (red X = tossed)
%           n / right = next cell    p / left = prev cell    s = save    q/Enter = save+quit
%
% Shared assignment grpOf(r):  >0 cell id | 0 ungrouped | -1 tossed
% Runqi Zhang / 2026-07-24
close all;

%% ---- params ----
addpath(fileparts(mfilename('fullpath')));      % so ncc_round_cfg_260801 resolves
cellcfg     = ncc_round_cfg_260801();                 % <-- EDIT PATHS THERE, not here
datasetPath = cellcfg.datasetPath;
outDirName  = cellcfg.outDirName;
tileSize = 130;  nCols = 6;

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));  addpath(repoRoot);
outDir   = fullfile(datasetPath, outDirName);
[roi, fov_name, fov_folder, fov_pxum, params, cellMembers, memberNCC] = load_results(outDir);
nccThr   = params.nccMergeThr;
nRoiTot  = numel(roi.fov);  nCells = numel(cellMembers);
[grpOf, homeOf] = load_or_init_state(outDir, cellMembers, nRoiTot);
[pC, mC, ~] = get_patches(outDir, roi, fov_folder, fov_pxum, params, nRoiTot);   % cache-first (fast)

hFig = figure('Color','w','Name','Round 1: toss','Position',[100 70 nCols*tileSize+140 740]);
gi = 1;
while true
    mem = cellMembers{gi};
    tossed = grpOf(mem) < 0;
    mapper = draw_toss(hFig, gi, nCells, mem, tossed, memberNCC{gi}, pC, mC, roi, fov_name, nccThr, tileSize, nCols);
    [xx, yy, btn] = ginput(1);
    if isempty(btn), break; end
    if btn == 1
        idx = mapper(xx, yy);
        if idx >= 1 && idx <= numel(mem)
            r = mem(idx);
            if grpOf(r) < 0, grpOf(r) = homeOf(r);          % un-toss -> back to its cell
            else, homeOf(r) = grpOf(r); grpOf(r) = -1; end  % toss
        end
    else
        switch lower(char(btn))
            case {'n', char(29)}, gi = min(gi+1, nCells);
            case {'p', char(28)}, gi = max(gi-1, 1);
            case 's', save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers); fprintf('saved\n');
            case {'q', char(13)}, break;
        end
    end
end
save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers);
fprintf('Round 1 done: %d observations tossed. -> %s\n', nnz(grpOf<0), fullfile(outDir,'roi_match_curated.mat'));
if ishandle(hFig), close(hFig); end
end

%% ================= shared helpers (copied into each round) =================
function [roi, fov_name, fov_folder, fov_pxum, params, cellMembers, memberNCC] = load_results(outDir)
resPath = fullfile(outDir, 'roi_match_results.mat');
assert(isfile(resPath), 'Not found: %s (run roi_pair_morph_match_260727.m first)', resPath);
S = load(resPath, 'match');  M = S.match;
roi=M.roi; fov_name=M.fov_name; fov_folder=M.fov_folder; fov_pxum=M.fov_pxum;
params=M.params; cellMembers=M.cellMembers; memberNCC=M.memberNCC;
end

function [fov_raw, fov_maskL] = reload_fovs(roi, fov_folder)
nFov = numel(fov_folder);  fov_raw = cell(nFov,1);  fov_maskL = cell(nFov,1);
fprintf('Reloading FOV images/masks...\n');
for f = unique(roi.fov(:)).'
    folderPath = char(fov_folder(f));
    sh = dir(fullfile(folderPath, '*_cpSAM_output.mat'));
    if isempty(sh), continue; end
    Q = load(fullfile(sh(1).folder, sh(1).name), 'maskL');
    if isfield(Q,'maskL'), fov_maskL{f} = Q.maskL; end
    avgIm = load_avg_image(fullfile(sh(1).folder, sh(1).name), folderPath);
    if ~isempty(avgIm), fov_raw{f} = double(avgIm); end
end
end
function I = norm_img(avgIm, clip, gamma)
lo = prctile(avgIm(:), clip(1)); hi = prctile(avgIm(:), clip(2)); if hi <= lo, hi = lo + 1; end
I = min(max((avgIm - lo)/(hi - lo), 0), 1); if gamma ~= 1, I = I .^ gamma; end
end

function [pC, mC, Umat] = get_patches(outDir, roi, fov_folder, fov_pxum, params, nRoiTot)
% Cache-first (display look already baked). Fallback: DISPLAY patches use
% dispClip+gamma; SCORING Umat uses linear imgClip.
cp = fullfile(outDir,'roi_match_patches.mat');
if isfile(cp)
    C = load(cp,'patchCache');
    if isfield(C,'patchCache') && isfield(C.patchCache,'pAll') && numel(C.patchCache.pAll)==nRoiTot
        pc = C.patchCache; pC = pc.pAll; mC = pc.mAll; Umat = double(pc.Umat);
        for r = 1:numel(pC), if ~isempty(pC{r}), pC{r} = double(pC{r}); end, end
        fprintf('Loaded patch cache (%d ROIs).\n', nRoiTot); return;
    end
end
fprintf('No patch cache; reloading FOVs + extracting (slower)...\n');
viewHalf = getfielddef(params,'viewHalf_um',16); target = params.targetUmPerPx; patchHalf = params.patchHalf_um;
imgClip  = getfielddef(params,'imgClip',[5 99.5]);
dispClip = getfielddef(params,'dispClip',[1 99.9]); dispGamma = getfielddef(params,'dispGamma',0.6);
[fov_raw, fov_maskL] = reload_fovs(roi, fov_folder);
fov_img = cell(numel(fov_raw),1);
for f = 1:numel(fov_raw)
    if ~isempty(fov_raw{f}), fov_img{f} = norm_img(fov_raw{f}, imgClip, 1); end
end
nPix = (round(2*patchHalf/target)+1)^2;
pC = cell(nRoiTot,1); mC = cell(nRoiTot,1); Umat = nan(nRoiTot,nPix);
for r = 1:nRoiTot
    f = roi.fov(r); if f<1||f>numel(fov_raw)||isempty(fov_raw{f}), continue; end
    pC{r} = disp_norm_patch(extract_patch(fov_raw{f}, roi.cx_px(r), roi.cy_px(r), viewHalf, fov_pxum(f), target), dispClip, dispGamma);
    mC{r} = extract_patch(double(fov_maskL{f}==roi.lab(r)), roi.cx_px(r), roi.cy_px(r), viewHalf, fov_pxum(f), target) > 0.5;
    t = extract_patch(fov_img{f}, roi.cx_px(r), roi.cy_px(r), patchHalf, fov_pxum(f), target);
    if ~isempty(t), v = t(:)-mean(t(:)); nv = norm(v); if nv>0, Umat(r,:) = (v/nv).'; end, end
end
end

function [grpOf, homeOf] = load_or_init_state(outDir, cellMembers, nRoiTot)
cp = fullfile(outDir,'roi_match_curated.mat');
if isfile(cp)
    C = load(cp,'curated');
    if isfield(C,'curated') && isfield(C.curated,'grpOf') && numel(C.curated.grpOf)==nRoiTot
        grpOf = C.curated.grpOf;
        if isfield(C.curated,'homeOf') && numel(C.curated.homeOf)==nRoiTot, homeOf = C.curated.homeOf; else, homeOf = grpOf; end
        fprintf('Resuming from existing curation.\n');
        return;
    elseif isfield(C,'curated')
        warning('roi_match_curated.mat ROI count mismatch (ROI set changed) - starting curation fresh.');
    end
end
grpOf = zeros(nRoiTot,1);
for gi = 1:numel(cellMembers), grpOf(cellMembers{gi}) = gi; end
homeOf = grpOf;
end

function save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers)
% Curated partition = every NON-tossed ROI. grpOf>0 form their cells; each
% remaining ungrouped ROI (grpOf==0) is emitted as its OWN 1x cell so ALL
% detected ROIs are kept (tossed grpOf<0 excluded). grpOf saved unchanged
% (0 stays 'pending' so the rounds can be re-run).
posIds = sort(unique(grpOf(grpOf>0)));
groups = cell(0,1);
for g = 1:numel(posIds)
    m = find(grpOf==posIds(g)); [~,o] = sort(roi.fov(m)); groups{end+1,1} = m(o); %#ok<AGROW>
end
ungr = find(grpOf==0);
for k = 1:numel(ungr), groups{end+1,1} = ungr(k); end %#ok<AGROW>
curated = struct('grpOf',grpOf,'homeOf',homeOf, ...
    'curatedCellMembers',{groups},'original',{cellMembers},'roi',roi,'fov_name',fov_name);
save(fullfile(outDir,'roi_match_curated.mat'),'curated','-v7.3');
rows = {};
for g = 1:numel(groups)
    m = groups{g};
    for c = 1:numel(m)
        rr = m(c);
        rows(end+1,:) = {g, numel(m), char(fov_name(roi.fov(rr))), roi.roi(rr), ...
                         roi.x(rr), roi.y(rr), roi.z(rr)}; %#ok<AGROW>
    end
end
if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', ...
        {'cell_id','cell_size','fov_name','roi_index','x_um','y_um','z_um'});
    writetable(T, fullfile(outDir,'roi_match_curated.csv'));
end
end

function mapper = draw_toss(hFig, gi, nCells, mem, tossed, mn, pC, mC, roi, fov_name, nccThr, TS, nCols)
figure(hFig); clf(hFig);
N = numel(mem); nc = max(min(N,nCols),1); nr = max(ceil(N/nc),1); GAP = 8;
Wc = nc*(TS+GAP)+GAP; Hc = nr*(TS+GAP)+GAP;
img = ones(Hc,Wc)*0.15; mimg = false(Hc,Wc); box = zeros(N,4);
for c = 1:N
    r = mem(c); pv = pC{r}; mv = mC{r};
    row = ceil(c/nc); col = c-(row-1)*nc;
    x0 = GAP+(col-1)*(TS+GAP); y0 = GAP+(row-1)*(TS+GAP); box(c,:) = [x0 y0 x0+TS-1 y0+TS-1];
    if ~isempty(pv)
        img(y0:y0+TS-1,x0:x0+TS-1) = imresize(pv,[TS TS]);
        mimg(y0:y0+TS-1,x0:x0+TS-1) = imresize(mv,[TS TS]) > 0.5;
    end
end
ax = axes('Parent',hFig); imagesc(ax,img); colormap(ax,gray); axis(ax,'image'); axis(ax,'off'); hold(ax,'on');
B = bwboundaries(mimg); for b = 1:numel(B), plot(ax,B{b}(:,2),B{b}(:,1),'y-','LineWidth',1); end
for c = 1:N
    r = mem(c); x0 = box(c,1); y0 = box(c,2); x1 = box(c,3); y1 = box(c,4);
    tcol = 'w'; if mn(c) < nccThr, tcol = [1 .4 .4]; end
    text(ax, x0+3, y0+11, sprintf('%s #%d', fovLabel(fov_name(roi.fov(r))), roi.roi(r)), ...
        'Color',tcol,'FontSize',7,'Interpreter','none','FontWeight','bold');
    text(ax, x0+3, y1-6, sprintf('r=%.2f', mn(c)), 'Color',tcol,'FontSize',7);
    if tossed(c)
        plot(ax,[x0 x1],[y0 y1],'r-','LineWidth',2); plot(ax,[x0 x1],[y1 y0],'r-','LineWidth',2);
        plot(ax,[x0 x0 x1 x1 x0],[y0 y1 y1 y0 y0],'r-','LineWidth',2);
    end
end
title(ax, sprintf('ROUND 1 TOSS   cell %d/%d   %dx (keep %d)   [click=toss | n/p | s | q]', ...
    gi, nCells, N, nnz(~tossed)), 'Interpreter','none','FontSize',9);
drawnow;
    function idx = m_(x,y)
        col_ = floor((x-GAP)/(TS+GAP))+1; row_ = floor((y-GAP)/(TS+GAP))+1;
        idx = (row_-1)*nc + col_;
        if col_<1 || col_>nc || idx<1 || idx>N, idx = 0; end
    end
mapper = @m_;
end

function v = getfielddef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function s = fovLabel(name)
name = char(name);
if numel(name) > 34, s = [name(1:20) '..' name(end-11:end)]; else, s = name; end
end
function im = load_avg_image(samPath, folderPath)
im = [];
try
    S = load(samPath, 'avgPath');
    if isfield(S,'avgPath') && ~isempty(S.avgPath) && isfile(S.avgPath), im = imread(char(S.avgPath)); end
catch
end
if isempty(im)
    cand = dir(fullfile(folderPath, '*_AVG_for_CP.tif'));
    if isempty(cand), cand = dir(fullfile(folderPath, '*_AVG*.tif')); end
    if ~isempty(cand), im = imread(fullfile(cand(1).folder, cand(1).name)); end
end
if ~isempty(im) && size(im,3) > 1, im = im(:,:,1); end
end
function p = disp_norm_patch(p, clip, gamma)
% PER-ROI-PATCH display normalization (display-only; never touches matching).
if isempty(p), return; end
lo = prctile(p(:), clip(1)); hi = prctile(p(:), clip(2)); if hi <= lo, hi = lo + 1; end
p = min(max((p - lo)/(hi - lo), 0), 1); if gamma ~= 1, p = p .^ gamma; end
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
