function roi_review3_merge_260727()
%% ROUND 3 of 4 - MERGE same-cell splits, NCC-suggested.
% Runs a SECOND round of NCC, this time BETWEEN cells (medoid vs medoid, only
% for cells that don't already share a FOV), plots the NCC-vs-IoU scatter + NCC
% distribution (like the first pass), then steps you through the suggested
% merges (NCC > mergeSuggestThr) ranked high-to-low.
%
% CONTROLS (review window):
%   y ... MERGE this pair (2nd cell -> 1st)   n / space / right ... skip
%   p / left ... back      s ... save      q / Enter ... save+quit
%
% Shared assignment grpOf(r):  >0 cell id | 0 ungrouped | -1 tossed
% Runqi Zhang / 2026-07-24
close all;

%% ---- params ----
addpath(fileparts(mfilename('fullpath')));      % so cell_cfg_260727 resolves
cellcfg     = cell_cfg_260727();                 % <-- EDIT PATHS THERE, not here
datasetPath = cellcfg.datasetPath;
outDirName  = cellcfg.outDirName;
maxMergeDist_um = 300;   % only consider merging cells whose centroids are within this
mergeSuggestThr = 0.4;   % only review candidate merges with medoid NCC > this
tileSize = 200;

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));  addpath(repoRoot);
outDir   = fullfile(datasetPath, outDirName);
[roi, fov_name, fov_folder, fov_pxum, params, cellMembers, ~] = load_results(outDir);
nRoiTot  = numel(roi.fov);
[grpOf, homeOf] = load_or_init_state(outDir, cellMembers, nRoiTot);
[pC, mC, Umat] = get_patches(outDir, roi, fov_folder, fov_pxum, params, nRoiTot);   % cache-first (fast)

%% ---- 2nd round of NCC: cell-vs-cell (medoid), disjoint-FOV only ----
[pairs, nccV, iouV, distV] = suggest_merges(grpOf, roi, Umat, mC, maxMergeDist_um);
fprintf('Cell-vs-cell candidates (disjoint FOV, <=%.0fum): %d  |  NCC>%.2f: %d\n', ...
    maxMergeDist_um, size(pairs,1), mergeSuggestThr, nnz(nccV>mergeSuggestThr));
draw_merge_scatter(nccV, iouV, mergeSuggestThr);        % scatter + NCC histogram (non-blocking)

%% ---- review suggestions, high NCC first ----
keep = find(nccV > mergeSuggestThr);
[~, ord] = sort(nccV(keep), 'descend');  order = keep(ord);
hFig = figure('Color','w','Name','Round 3: suggested merges','Position',[120 80 2*tileSize+80 tileSize+150]);
qi = 1;
while qi <= numel(order)
    p = order(qi);  a = pairs(p,1);  b = pairs(p,2);
    ma = find(grpOf==a);  mb = find(grpOf==b);
    if isempty(ma) || isempty(mb) || any(ismember(roi.fov(ma), roi.fov(mb)))
        qi = qi + 1; continue;      % a cell was merged away / now shares a FOV
    end
    medA = medoid_of(a, grpOf, Umat);  medB = medoid_of(b, grpOf, Umat);
    draw_pair(hFig, qi, numel(order), a, b, numel(ma), numel(mb), nccV(p), iouV(p), distV(p), ...
        medA, medB, pC, mC, roi, fov_name, tileSize);
    [~, ~, btn] = ginput(1);
    if isempty(btn), break; end
    switch lower(char(btn))
        case 'y', grpOf(grpOf==b) = a; qi = qi + 1;             % merge B into A
        case {'n', ' ', char(29)}, qi = qi + 1;
        case {'p', char(28)}, qi = max(qi-1, 1);
        case 's', save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers); fprintf('saved\n');
        case {'q', char(13)}, break;
    end
end
save_state(outDir, grpOf, homeOf, roi, fov_name, cellMembers);
fprintf('Round 3 done. %d cells remain. -> %s\n', numel(unique(grpOf(grpOf>0))), fullfile(outDir,'roi_match_curated.mat'));
if ishandle(hFig), close(hFig); end
end

%% ---- round-3 helpers ----
function r = medoid_of(id, grpOf, Umat)
mem = find(grpOf==id);
if numel(mem) == 1, r = mem; return; end
G = Umat(mem,:) * Umat(mem,:)';
sc = nan(numel(mem),1);
for c = 1:numel(mem), o = true(numel(mem),1); o(c) = false; sc(c) = mean(G(c,o),'omitnan'); end
[~, mi] = max(sc); r = mem(mi);
end

function [pairs, nccV, iouV, distV] = suggest_merges(grpOf, roi, Umat, mC, maxDist)
% cell-vs-cell NCC/IoU on medoid representatives, only for cells with DISJOINT
% FOV sets (merging same-FOV cells would break one-ROI-per-FOV) and within maxDist.
ids = sort(unique(grpOf(grpOf>0)));  n = numel(ids);
med = zeros(n,1); cx = zeros(n,1); cy = zeros(n,1);
nFov = max([1; roi.fov(:)]);  cellFov = false(n, nFov);
for k = 1:n
    m = find(grpOf==ids(k));
    med(k) = medoid_of(ids(k), grpOf, Umat);
    cx(k) = mean(roi.x(m));  cy(k) = mean(roi.y(m));
    cellFov(k, unique(roi.fov(m))) = true;
end
pairs = zeros(0,2); nccV = zeros(0,1); iouV = zeros(0,1); distV = zeros(0,1);
for a = 1:n-1
    ua = Umat(med(a),:);  if any(isnan(ua)), continue; end
    for b = a+1:n
        if any(cellFov(a,:) & cellFov(b,:)), continue; end      % share a FOV -> not mergeable
        d = hypot(cx(a)-cx(b), cy(a)-cy(b));
        if d > maxDist, continue; end
        ub = Umat(med(b),:);  if any(isnan(ub)), continue; end
        ncc = ua * ub.';
        iou = NaN;  A = mC{med(a)};  Bm = mC{med(b)};
        if ~isempty(A) && ~isempty(Bm) && isequal(size(A),size(Bm))
            u = nnz(A|Bm);  if u>0, iou = nnz(A&Bm)/u; end
        end
        pairs(end+1,:) = [ids(a) ids(b)]; nccV(end+1,1) = ncc; iouV(end+1,1) = iou; distV(end+1,1) = d; %#ok<AGROW>
    end
end
end

function draw_merge_scatter(nccV, iouV, thr)
figure('Color','w','Name','Round 3: cell-cell NCC','Position',[80 300 940 380]);
hi = nccV > thr;
subplot(1,2,1);
scatter(nccV(~hi), iouV(~hi), 18, [.6 .6 .6], 'filled', 'MarkerFaceAlpha',0.4); hold on;
scatter(nccV(hi),  iouV(hi),  30, 'r', 'filled');
xline(thr,'r--'); xlabel('cell-cell NCC (medoid)'); ylabel('mask IoU'); grid on; box on;
title(sprintf('%d candidate cell-pairs  (%d > %.2f)', numel(nccV), nnz(hi), thr));
subplot(1,2,2);
if ~isempty(nccV)
    histogram(nccV, linspace(min([0;nccV]),1,40), 'FaceColor',[.3 .3 .8]); hold on;
    xline(thr,'r--','LineWidth',1);
end
xlabel('cell-cell NCC'); ylabel('# pairs'); title('NCC distribution'); grid on; box on;
drawnow;
end

function draw_pair(hFig, qi, nQ, a, b, na, nb, ncc, iou, dist, medA, medB, pC, mC, roi, fov_name, TS)
figure(hFig); clf(hFig);
GAP = 10;  Wc = 2*TS + 3*GAP;  Hc = TS + 2*GAP;
img = ones(Hc,Wc)*0.15; mimg = false(Hc,Wc);
placeone(medA, GAP);  placeone(medB, 2*GAP+TS);
ax = axes('Parent',hFig); imagesc(ax,img); colormap(ax,gray); axis(ax,'image'); axis(ax,'off'); hold(ax,'on');
B = bwboundaries(mimg); for k = 1:numel(B), plot(ax,B{k}(:,2),B{k}(:,1),'y-','LineWidth',1); end
text(ax, GAP+3, GAP+13, sprintf('cell %d (%dx)  %s', a, na, fovLabel(fov_name(roi.fov(medA)))), ...
    'Color','c','FontSize',9,'FontWeight','bold','Interpreter','none');
text(ax, 2*GAP+TS+3, GAP+13, sprintf('cell %d (%dx)  %s', b, nb, fovLabel(fov_name(roi.fov(medB)))), ...
    'Color','c','FontSize',9,'FontWeight','bold','Interpreter','none');
title(ax, sprintf('MERGE? %d/%d   NCC=%.3f  IoU=%.2f  d=%.0fum    [y=merge  n=skip  p=back  s=save  q=quit]', ...
    qi, nQ, ncc, iou, dist), 'Interpreter','none','FontSize',10);
drawnow;
    function placeone(r, x0)
        pv = pC{r};  y0 = GAP;
        if isempty(pv), return; end
        img(y0:y0+TS-1, x0:x0+TS-1)  = imresize(pv,[TS TS]);
        if ~isempty(mC{r}), mimg(y0:y0+TS-1, x0:x0+TS-1) = imresize(mC{r},[TS TS]) > 0.5; end
    end
end

%% ================= shared helpers =================
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
