function roi_review4_ungrouped_260724()
%% ROUND 4 of 4 - place NEVER-GROUPED observations, in two phases.
% PHASE 1 (triage): grid of all ungrouped ROIs; click the junk to TOSS.
% PHASE 2 (AUTO-SUGGEST): for each surviving ungrouped ROI its BEST-matching
%   existing cell (medoid NCC) is proposed when NCC > assignThr; you just
%   confirm. Shows the NCC-vs-IoU scatter + histogram first, then steps
%   high->low. Whatever you SKIP stays ungrouped and is kept as its own 1x
%   cell on save (nothing is dropped).
%
% PHASE 1: left-click = toss/un-toss | n/p = page | d = done->phase 2 | q = save+quit
% PHASE 2: y = assign | n/space = skip | x = toss | p = back | s = save | q/Enter = save+quit
%
% Shared assignment grpOf(r):  >0 cell id | 0 ungrouped | -1 tossed
% Runqi Zhang / 2026-07-24
close all;

%% ---- params ----
datasetPath = "D:\260721_Sert_soma_G8s\phys";
outDirName  = 'roi_match_out_260724';
maxAssignDist_um = 200;  % phase 2: only match an ungrouped ROI to cells within this
assignThr        = 0.3;  % phase 2: only suggest an assignment with medoid NCC > this
assignTile = 200;                                  % phase 2: tile size for the ROI|cell view
gridTile   = 95;   gridCols   = 8;  gridRows = 5;   % phase 1: triage grid

repoRoot = fileparts(fileparts(mfilename('fullpath')));  addpath(repoRoot);
outDir   = fullfile(datasetPath, outDirName);
[roi, fov_name, fov_folder, fov_pxum, params, cellMembers, ~] = load_results(outDir);
nRoiTot  = numel(roi.fov);
[grpOf, homeOf] = load_or_init_state(outDir, cellMembers, nRoiTot);
[pC, mC, Umat] = get_patches(outDir, roi, fov_folder, fov_pxum, params, nRoiTot);   % cache-first (fast)

hFig = figure('Color','w','Name','Round 4: ungrouped','Position',[80 60 gridCols*gridTile+40 gridRows*gridTile+120]);

%% ================= PHASE 1: triage-toss ungrouped =================
ung0 = find(grpOf == 0);            % snapshot of never-grouped ROIs at entry
perPage = gridCols*gridRows;  page = 1;  goPhase2 = false;
while ~isempty(ung0)
    nPages = max(ceil(numel(ung0)/perPage),1);  page = min(max(page,1),nPages);
    pageIdx = ung0((page-1)*perPage+1 : min(page*perPage, numel(ung0)));
    tossed  = grpOf(pageIdx) < 0;
    mapper = draw_grid_toss(hFig, page, nPages, pageIdx, tossed, pC, mC, roi, fov_name, gridTile, gridCols, gridRows);
    [xx, yy, btn] = ginput(1);
    if isempty(btn), break; end
    if btn == 1
        ti = mapper(xx, yy);
        if ti>=1 && ti<=numel(pageIdx)
            r = pageIdx(ti);
            if grpOf(r)==0, grpOf(r) = -1; else, grpOf(r) = 0; end   % toss / un-toss
        end
    else
        switch lower(char(btn))
            case {'n', char(29)}, page = page + 1;
            case {'p', char(28)}, page = page - 1;
            case 'd', goPhase2 = true; break;
            case 's', save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers); fprintf('saved\n');
            case 'q', break;
        end
    end
end
if isempty(ung0), goPhase2 = true; end

%% ================= PHASE 2: AUTO-SUGGEST assignments (NCC > assignThr) =================
if goPhase2
    [rlist, clist, nccV, iouV, distV] = suggest_assign(grpOf, roi, Umat, mC, maxAssignDist_um);
    fprintf('Ungrouped ROIs with a candidate cell: %d  |  NCC>%.2f: %d\n', ...
        numel(rlist), assignThr, nnz(nccV>assignThr));
    draw_merge_scatter(nccV, iouV, assignThr);        % best-match NCC vs IoU + histogram
    keep = find(nccV > assignThr);
    [~, ord] = sort(nccV(keep), 'descend');  order = keep(ord);
    qi = 1;
    while qi <= numel(order)
        s = order(qi);  r = rlist(s);  cid = clist(s);
        mem = find(grpOf==cid);
        if grpOf(r)~=0 || isempty(mem) || any(roi.fov(mem)==roi.fov(r))
            qi = qi + 1; continue;      % r already handled / cell gone / FOV now conflicts
        end
        medC = medoid_of(cid, grpOf, Umat);
        draw_suggest(hFig, qi, numel(order), r, cid, numel(mem), nccV(s), iouV(s), distV(s), ...
            medC, pC, mC, roi, fov_name, assignTile);
        [~, ~, btn] = ginput(1);
        if isempty(btn), break; end
        switch lower(char(btn))
            case 'y', grpOf(r) = cid; qi = qi + 1;              % assign
            case {'n', ' ', char(29)}, qi = qi + 1;             % skip -> stays 1x cell on save
            case 'x', grpOf(r) = -1; qi = qi + 1;               % toss
            case {'p', char(28)}, qi = max(qi-1, 1);
            case 's', save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers); fprintf('saved\n');
            case {'q', char(13)}, break;
        end
    end
end
save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers);
fprintf('Round 4 done. -> %s\n', fullfile(outDir,'roi_match_curated.mat'));
if ishandle(hFig), close(hFig); end
end

%% ---- assign / grid helpers ----
function [rlist, clist, nccV, iouV, distV] = suggest_assign(grpOf, roi, Umat, mC, maxDist)
% For each ungrouped ROI, its BEST-matching existing cell (medoid NCC), among
% cells within maxDist that DON'T already contain the ROI's FOV (one-per-FOV).
ids = sort(unique(grpOf(grpOf>0)));  n = numel(ids);
med = zeros(n,1); cx = zeros(n,1); cy = zeros(n,1);
nFov = max([1; roi.fov(:)]);  cellFov = false(n, nFov);
for k = 1:n
    m = find(grpOf==ids(k));
    med(k) = medoid_of(ids(k), grpOf, Umat);
    cx(k) = mean(roi.x(m));  cy(k) = mean(roi.y(m));
    cellFov(k, unique(roi.fov(m))) = true;
end
ung = find(grpOf==0);
rlist = zeros(0,1); clist = zeros(0,1); nccV = zeros(0,1); iouV = zeros(0,1); distV = zeros(0,1);
for ri = 1:numel(ung)
    r = ung(ri);  ur = Umat(r,:);  if any(isnan(ur)), continue; end
    fr = roi.fov(r);  rx = roi.x(r);  ry = roi.y(r);
    best = -Inf; bk = 0; bd = NaN; biou = NaN;
    for k = 1:n
        if cellFov(k, fr), continue; end                       % cell already has this FOV
        d = hypot(rx-cx(k), ry-cy(k));  if d > maxDist, continue; end
        uk = Umat(med(k),:);  if any(isnan(uk)), continue; end
        ncc = ur * uk.';
        if ncc > best
            best = ncc; bk = k; bd = d;  biou = NaN;
            A = mC{r}; Bm = mC{med(k)};
            if ~isempty(A) && ~isempty(Bm) && isequal(size(A),size(Bm))
                u = nnz(A|Bm); if u>0, biou = nnz(A&Bm)/u; end
            end
        end
    end
    if bk > 0
        rlist(end+1,1)=r; clist(end+1,1)=ids(bk); nccV(end+1,1)=best; iouV(end+1,1)=biou; distV(end+1,1)=bd; %#ok<AGROW>
    end
end
end

function r = medoid_of(id, grpOf, Umat)
mem = find(grpOf==id);
if numel(mem) == 1, r = mem; return; end
G = Umat(mem,:) * Umat(mem,:)';
sc = nan(numel(mem),1);
for c = 1:numel(mem), o = true(numel(mem),1); o(c) = false; sc(c) = mean(G(c,o),'omitnan'); end
[~, mi] = max(sc); r = mem(mi);
end

function mapper = draw_grid_toss(hFig, page, nPages, pageIdx, tossed, pC, mC, roi, fov_name, TS, nCols, nRows)
figure(hFig); clf(hFig);
K = numel(pageIdx); GAP = 6;
Wc = nCols*(TS+GAP)+GAP; Hc = nRows*(TS+GAP)+GAP;
img = ones(Hc,Wc)*0.15; mimg = false(Hc,Wc); box = zeros(K,4);
for c = 1:K
    r = pageIdx(c); row = ceil(c/nCols); col = c-(row-1)*nCols;
    x0 = GAP+(col-1)*(TS+GAP); y0 = GAP+(row-1)*(TS+GAP); box(c,:) = [x0 y0 x0+TS-1 y0+TS-1];
    pv = pC{r}; mv = mC{r};
    if ~isempty(pv)
        img(y0:y0+TS-1,x0:x0+TS-1) = imresize(pv,[TS TS]);
        mimg(y0:y0+TS-1,x0:x0+TS-1) = imresize(mv,[TS TS]) > 0.5;
    end
end
ax = axes('Parent',hFig); imagesc(ax,img); colormap(ax,gray); axis(ax,'image'); axis(ax,'off'); hold(ax,'on');
B = bwboundaries(mimg); for b = 1:numel(B), plot(ax,B{b}(:,2),B{b}(:,1),'y-','LineWidth',0.5); end
for c = 1:K
    r = pageIdx(c); x0 = box(c,1); y0 = box(c,2); x1 = box(c,3); y1 = box(c,4);
    text(ax, x0+2, y0+9, sprintf('%s #%d', fovLabel(fov_name(roi.fov(r))), roi.roi(r)), ...
        'Color','w','FontSize',6,'Interpreter','none');
    if tossed(c)
        plot(ax,[x0 x1],[y0 y1],'r-','LineWidth',2); plot(ax,[x0 x1],[y1 y0],'r-','LineWidth',2);
        plot(ax,[x0 x0 x1 x1 x0],[y0 y1 y1 y0 y0],'r-','LineWidth',2);
    end
end
title(ax, sprintf('ROUND 4 PHASE 1 (toss ungrouped)   page %d/%d   %d shown, %d tossed   [click=toss | n/p | d=done | q]', ...
    page, nPages, K, nnz(tossed)), 'Interpreter','none','FontSize',9);
drawnow;
    function idx = m_(x,y)
        col_ = floor((x-GAP)/(TS+GAP))+1; row_ = floor((y-GAP)/(TS+GAP))+1;
        idx = (row_-1)*nCols + col_;
        if col_<1 || col_>nCols || idx<1 || idx>K, idx = 0; end
    end
mapper = @m_;
end

function draw_merge_scatter(nccV, iouV, thr)
figure('Color','w','Name','Round 4: ROI-cell NCC','Position',[80 300 940 380]);
hi = nccV > thr;
subplot(1,2,1);
scatter(nccV(~hi), iouV(~hi), 18, [.6 .6 .6], 'filled', 'MarkerFaceAlpha',0.4); hold on;
scatter(nccV(hi),  iouV(hi),  30, 'r', 'filled');
xline(thr,'r--'); xlabel('best ROI-cell NCC (medoid)'); ylabel('mask IoU'); grid on; box on;
title(sprintf('%d ungrouped with a match  (%d > %.2f)', numel(nccV), nnz(hi), thr));
subplot(1,2,2);
if ~isempty(nccV)
    histogram(nccV, linspace(min([0;nccV]),1,40), 'FaceColor',[.3 .3 .8]); hold on;
    xline(thr,'r--','LineWidth',1);
end
xlabel('best ROI-cell NCC'); ylabel('# ungrouped ROIs'); title('NCC distribution'); grid on; box on;
drawnow;
end

function draw_suggest(hFig, qi, nQ, r, cid, nc, ncc, iou, dist, medC, pC, mC, roi, fov_name, TS)
figure(hFig); clf(hFig);
GAP = 10;  Wc = 2*TS + 3*GAP;  Hc = TS + 2*GAP;
img = ones(Hc,Wc)*0.15; mimg = false(Hc,Wc);
placeone(r, GAP);  placeone(medC, 2*GAP+TS);
ax = axes('Parent',hFig); imagesc(ax,img); colormap(ax,gray); axis(ax,'image'); axis(ax,'off'); hold(ax,'on');
B = bwboundaries(mimg); for k = 1:numel(B), plot(ax,B{k}(:,2),B{k}(:,1),'y-','LineWidth',1); end
text(ax, GAP+3, GAP+13, sprintf('ungrouped ROI  %s #%d', fovLabel(fov_name(roi.fov(r))), roi.roi(r)), ...
    'Color','c','FontSize',9,'FontWeight','bold','Interpreter','none');
text(ax, 2*GAP+TS+3, GAP+13, sprintf('-> cell %d (%dx)', cid, nc), ...
    'Color','c','FontSize',9,'FontWeight','bold','Interpreter','none');
title(ax, sprintf('ASSIGN? %d/%d   NCC=%.3f  IoU=%.2f  d=%.0fum    [y=assign  n=skip  x=toss  p=back  s=save  q=quit]', ...
    qi, nQ, ncc, iou, dist), 'Interpreter','none','FontSize',10);
drawnow;
    function placeone(rr, x0)
        pv = pC{rr};  y0 = GAP;
        if isempty(pv), return; end
        img(y0:y0+TS-1, x0:x0+TS-1)  = imresize(pv,[TS TS]);
        if ~isempty(mC{rr}), mimg(y0:y0+TS-1, x0:x0+TS-1) = imresize(mC{rr},[TS TS]) > 0.5; end
    end
end

%% ================= shared helpers =================
function [roi, fov_name, fov_folder, fov_pxum, params, cellMembers, memberNCC] = load_results(outDir)
resPath = fullfile(outDir, 'roi_match_results.mat');
assert(isfile(resPath), 'Not found: %s (run roi_pair_morph_match_260724.m first)', resPath);
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
ids = sort(unique(grpOf(grpOf>0)));
curatedCellMembers = cell(numel(ids),1);
for g = 1:numel(ids)
    m = find(grpOf==ids(g)); [~,o] = sort(roi.fov(m)); curatedCellMembers{g} = m(o);
end
curated = struct('grpOf',grpOf,'homeOf',homeOf,'cellIds',ids, ...
    'curatedCellMembers',{curatedCellMembers},'original',{cellMembers},'roi',roi,'fov_name',fov_name);
save(fullfile(outDir,'roi_match_curated.mat'),'curated','-v7.3');
rows = {};
for g = 1:numel(ids)
    m = curatedCellMembers{g};
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
