function cell_activity_summary_260726()
%% cell_activity_summary_260726
% For each matched CELL (a group of same-cell ROIs across recordings), pool the
% activity into ONE summary figure:
%   LEFT  = the cell's average-projection patch (medoid observation) + outline.
%   RIGHT = one stacked row per recording: that recording's dF/F for this cell
%           (black) with its breathing trace overlaid (blue), on a seconds axis.
% One figure per cell, saved to <datasetPath>\cell_id_test\.
%
% Grouping comes from roi_match_curated.mat if present (your curation), else
% from roi_match_results.mat (auto grouping). dF/F from *_ch1_dFF.mat, breathing
% from breath_peak_pc1.mat (fallback breath_pc1.mat); breath is length-aligned
% to dF/F by dropping its front toss frames. fps via detect_session_fps.
%
% Runqi Zhang / 2026-07-26
close all;

%% ---- params ----
datasetPath = "D:\260721_Sert_soma_G8s\phys";
outDirName  = 'roi_match_out_260724';
outSub      = 'cell_id_test';
minCellSize = 1;       % 1 = INCLUDE ungrouped-untossed ROIs as their own 1x cell (full
                       %   inventory - every non-tossed ROI gets a figure). Set 2 for
                       %   pooled multi-recordings only. (Tossed grpOf<0 are always excluded.)
maxCells    = inf;      % cap for a test run (set Inf for all cells)
leftHalf_um = 35;      % half-size (um) of the LEFT avg-projection patch
showFirst   = 6;       % open this many figures on screen; the rest saved silently
fallback_fps = 30;

repoRoot = fileparts(fileparts(mfilename('fullpath')));  addpath(repoRoot);
outDir = fullfile(datasetPath, outDirName);
assert(isfile(fullfile(outDir,'roi_match_results.mat')), ...
    'roi_match_results.mat not found in %s (run the matcher first)', outDir);
M = load(fullfile(outDir,'roi_match_results.mat'),'match');  match = M.match;
roi = match.roi;  fov_name = match.fov_name;  fov_folder = match.fov_folder;  fov_pxum = match.fov_pxum;
params = match.params;
target    = params.targetUmPerPx;
patchHalf = params.patchHalf_um;
dispClip  = getfielddef(params,'dispClip',[1 99.9]);
dispGamma = getfielddef(params,'dispGamma',0.6);

% grouping: prefer curated. Rebuild cells from the per-ROI group vector so
% completeness is guaranteed regardless of how the file was saved:
%   curated grpOf:  >0 = a cell, 0 = ungrouped-untossed (-> own 1x cell), <0 = tossed (excluded)
%   results  grp :  per-ROI group id (singletons already their own group)
cp = fullfile(outDir,'roi_match_curated.mat');
cells = cell(0,1);
if isfile(cp)
    C = load(cp,'curated');  g = C.curated.grpOf;  src = 'curated';
    posIds = sort(unique(g(g>0)));
    for i = 1:numel(posIds), cells{end+1,1} = find(g==posIds(i)); end %#ok<AGROW>
    ung = find(g==0);                                   % ungrouped-untossed -> 1x cells
    for i = 1:numel(ung), cells{end+1,1} = ung(i); end  %#ok<AGROW>
else
    g = match.grp;  src = 'auto (results)';
    ids = sort(unique(g));
    for i = 1:numel(ids), cells{end+1,1} = find(g==ids(i)); end %#ok<AGROW>
end
sizes = cellfun(@numel, cells);
sel = find(sizes >= minCellSize);
[~, o] = sort(sizes(sel), 'descend');  sel = sel(o);
if numel(sel) > maxCells
    fprintf('Capping to %d of %d cells (set maxCells=Inf for all).\n', maxCells, numel(sel));
    sel = sel(1:maxCells);
end
fprintf('Grouping source: %s. Plotting %d cells (size>=%d).\n', src, numel(sel), minCellSize);
assert(~isempty(sel), 'No cells with size >= %d.', minCellSize);

outSubDir = fullfile(datasetPath, outSub);
if ~exist(outSubDir,'dir'), mkdir(outSubDir);
else, old = dir(fullfile(outSubDir,'cell_*.png')); for q=1:numel(old), delete(fullfile(old(q).folder,old(q).name)); end
end

%% ---- preload per-FOV data for the FOVs we need ----
neededFov = unique(roi.fov(vertcat(cells{sel})));
nF = numel(fov_folder);
DFF = cell(nF,1);  BR = cell(nF,1);  FPS = nan(nF,1);  AVG = cell(nF,1);  MASK = cell(nF,1);
fprintf('Loading dFF + breathing for %d FOVs...\n', numel(neededFov));
for fi = neededFov(:).'
    folder = char(fov_folder(fi));
    dd = dir(fullfile(folder,'*_ch1_dFF.mat'));
    if ~isempty(dd), S = load(fullfile(dd(1).folder,dd(1).name),'dFF'); if isfield(S,'dFF'), DFF{fi} = S.dFF; end, end
    BR{fi} = load_breath(folder);
    [FPS(fi), ~] = detect_session_fps(folder, fallback_fps);
    sh = dir(fullfile(folder,'*_cpSAM_output.mat'));
    if ~isempty(sh)
        Q = load(fullfile(sh(1).folder,sh(1).name),'maskL');
        if isfield(Q,'maskL'), MASK{fi} = Q.maskL; end
        AVG{fi} = load_avg_image(fullfile(sh(1).folder,sh(1).name), folder);
    end
end

%% ---- one figure per cell ----
hOff = figure('Color','w','Visible','off');
for si = 1:numel(sel)
    ci  = sel(si);
    mem = cells{ci};
    [~, mo] = sort(roi.fov(mem));  mem = mem(mo);      % order rows by FOV
    N   = numel(mem);
    md  = cell_medoid(mem, roi, AVG, fov_pxum, patchHalf, target);

    onScreen = si <= showFirst;
    if onScreen, hf = figure('Color','w','Visible','on'); else, hf = hOff; clf(hf); end
    set(hf, 'Position', [80 60 1150 max(360, 90*N+90)]);

    % --- LEFT: avg-projection patch of the medoid ---
    f = roi.fov(md);  px = fov_pxum(f);
    pv = disp_norm_patch(extract_patch(AVG{f}, roi.cx_px(md), roi.cy_px(md), leftHalf_um, px, target), dispClip, dispGamma);
    mv = extract_patch(double(MASK{f}==roi.lab(md)), roi.cx_px(md), roi.cy_px(md), leftHalf_um, px, target) > 0.5;
    axL = axes('Parent',hf,'Position',[0.03 0.30 0.28 0.44]);
    if ~isempty(pv), imagesc(axL, pv); end
    colormap(axL,gray); axis(axL,'image'); axis(axL,'off'); hold(axL,'on');
    draw_outline(axL, mv, 'y');
    title(axL, sprintf('cell %d  (%dx)\nmedoid: %s #%d', ci, N, fovLabel(fov_name(f)), roi.roi(md)), ...
        'Interpreter','none','FontSize',9);

    % --- RIGHT: stacked dF/F + breathing, one row per recording ---
    x0 = 0.37;  w = 0.60;  yTop = 0.95;  yBot = 0.07;  h = (yTop-yBot)/N;
    for i = 1:N
        r = mem(i);  f = roi.fov(r);  k = roi.roi(r);  fps = FPS(f);
        ax = axes('Parent',hf,'Position',[x0, yTop-i*h, w, h*0.86]);  hold(ax,'on');
        dff = [];
        if ~isempty(DFF{f}) && k>=1 && k<=size(DFF{f},2), dff = DFF{f}(:,k); end
        if isempty(dff)
            text(ax,0.5,0.5,'(no dF/F)','Units','normalized','HorizontalAlignment','center','Color',[.6 .6 .6]);
            axis(ax,'off');
        else
            br = BR{f};
            nToss = numel(br) - numel(dff);
            if ~isempty(br) && nToss > 0, br = br(nToss+1:end); end   % drop front toss frames
            n = numel(dff); if ~isempty(br), n = min(n, numel(br)); end
            dff = dff(1:n);  tt = (0:n-1)/fps;
            if ~isempty(br)
                brs = rescale(br(1:n), min(dff), max(dff));           % scale breath into dF/F range
                plot(ax, tt, brs, '-', 'Color', [0.2 0.45 1 0.55], 'LineWidth', 0.6);
            end
            plot(ax, tt, dff, 'k', 'LineWidth', 0.6);
            xlim(ax, [0 max(tt)]);  box(ax,'off');
            set(ax,'YTick',[],'FontSize',7);
            ylabel(ax, sprintf('%s\n#%d', fovLabel(fov_name(f)), k), 'Interpreter','none','FontSize',6,'Rotation',0, ...
                'HorizontalAlignment','right','VerticalAlignment','middle');
            if i < N, set(ax,'XTick',[]); else, xlabel(ax,'time (s)','FontSize',8); end
        end
    end

    exportgraphics(hf, fullfile(outSubDir, sprintf('cell_%03d_%02dx.png', ci, N)), 'Resolution', 130);
end
if ishandle(hOff), close(hOff); end
fprintf('Saved %d per-cell figures to %s\n', numel(sel), outSubDir);
end

%% ================= helpers =================
function br = load_breath(folder)
br = [];
bp = fullfile(folder,'breath_peak_pc1.mat');
if isfile(bp), S = load(bp); if isfield(S,'breath'), br = S.breath(:); return; end, end
bp = fullfile(folder,'breath_pc1.mat');
if isfile(bp), S = load(bp,'breathTrace'); if isfield(S,'breathTrace'), br = S.breathTrace(:); end, end
end

function md = cell_medoid(mem, roi, AVG, fov_pxum, patchHalf, target)
% medoid = member whose avg-patch has the highest mean NCC to the others
n = numel(mem);
if n == 1, md = mem; return; end
P = cell(n,1);
for i = 1:n
    r = mem(i); f = roi.fov(r);
    p = extract_patch(AVG{f}, roi.cx_px(r), roi.cy_px(r), patchHalf, fov_pxum(f), target);
    if ~isempty(p), v = p(:)-mean(p(:)); nv = norm(v); if nv>0, P{i} = v/nv; end, end
end
sc = nan(n,1);
for i = 1:n
    if isempty(P{i}), continue; end
    s = 0; c = 0;
    for j = 1:n, if j~=i && ~isempty(P{j}), s = s + P{i}.'*P{j}; c = c + 1; end, end
    if c > 0, sc(i) = s/c; end
end
[~, mi] = max(sc);  if isempty(mi) || isnan(sc(mi)), mi = 1; end
md = mem(mi);
end

function p = disp_norm_patch(p, clip, gamma)
if isempty(p), return; end
lo = prctile(p(:), clip(1)); hi = prctile(p(:), clip(2)); if hi <= lo, hi = lo + 1; end
p = min(max((p - lo)/(hi - lo), 0), 1); if gamma ~= 1, p = p .^ gamma; end
end

function draw_outline(ax, mask, col)
if isempty(mask) || ~any(mask(:)), return; end
B = bwboundaries(mask);
for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', col, 'LineWidth', 1); end
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
im = double(im);
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
