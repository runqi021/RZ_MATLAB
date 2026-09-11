%% roi_pair_morph_match_260727.m
% Cross-FOV ROI MATCHING for cell registration.
%
% Motivation: the same physical patch is sometimes imaged REPEATEDLY as
% several FOVs. Those repeats should be treated as "one ROI seen N times",
% not N different ROIs. This script finds, for every pair of ROIs from
% DIFFERENT FOVs, an offset-tolerant appearance match and a mask overlap,
% then renders a GALLERY of the top matches side by side so you can eyeball
% which pairs are the same cell.
%
% Within-FOV pairs are NOT considered (Cellpose already gives distinct cells).
%
% Appearance match = zero-shift normalized cross-correlation of the two AVG
% patches (each cut CENTERED on its own ROI centroid, resampled to a common um
% grid so different-zoom FOVs are comparable). Because both patches are already
% centroid-centered, corr2(patch_i,patch_j) == the NCC peak; no search window
% needed. patchHalf_um sets how much surrounding context the match sees.
%
% RE-ZERO / TOSS: if the stage was re-zeroed partway through a session, the
% motor coordinates of FOVs before vs after no longer share an origin, so
% their stage distance is meaningless. For now we TOSS such FOVs manually:
% put their folder path (or a distinctive substring of it) in tossFolders
% and they are dropped BEFORE loading. Each gallery tile also prints the
% pair distance in um, so an obvious visual match sitting at an implausibly
% large distance flags a re-zeroed FOV you should add to the toss list.
% (Future experiments: don't re-zero -> absolute coords stay comparable.)
%
% Standalone, plot-only. Reuses the loading logic of roi_pair_morph_quant.m.
% Runqi Zhang / 2026-07-24

clear; close all;

%% ---------------- USER PARAMS ----------------
addpath(fileparts(mfilename('fullpath')));      % so cell_cfg_260727 resolves
cellcfg     = cell_cfg_260727();                 % <-- EDIT PATHS THERE, not here
datasetPath = cellcfg.datasetPath;
xSign        = +1;
ySign        = -1;   % MEASURED 2026-08-30 on 260829: image row runs AGAINST stage y.
                     %   +1 gave 15.0% of ROIs a cross-FOV partner <10um (median NN 27.1um);
                     %   -1 gives 43.7% (median 11.7um), 2.92x, over 26 overlapping FOV pairs.
                     %   Same defect that halved Vgat/0730 (see reference_stage_axis_convention).
fallback_fps = 30;
PixelSizeBase = 1.7778;

% --- manual FOV toss (re-zeroed / bad FOVs) ---
% Each entry is a folder path OR a distinctive substring of one. Any FOV
% whose *_cpSAM_output.mat folder contains an entry (case-insensitive) is
% dropped before loading. Leave empty to keep all FOVs.
tossFolders = [ ...
    "ro1_2x_x1018y1130_18lp_3000f_00001"   % basename, NOT a full path: the data was
    ...  % "D:\260721_Sert_soma_G8s\phys\...\some_fov_folder"
    ];

% --- pairing / matching ---
maxPairDist_um = 300;    % XY prefilter: only score cross-FOV pairs within this (Inf/[] = all)
maxPairDz_um   = 30;     % ...and within this stage-Z difference (Inf/[] = all)
patchHalf_um   = 35;     % half-size (um) of the patch cut around each centroid = how much
                         %   context the NCC sees. ~11 = soma only; raise (~18) to pull in
                         %   neighbouring cells so adjacent similar somata stop confusing.
targetUmPerPx  = 0.5;    % common resampling grid (um/px). ~1.5 = optical lower limit and
                         %   coarser than every FOV's native pixel, so ALL patches are
                         %   DOWN-sampled to one baseline -> no invented detail, less noise,
                         %   different zooms made comparable. (0.5 upsampled -> noise emphasis.)
imgClip        = [1 99.5];
dispClip       = [1 99.9];   % DISPLAY-only clip, computed PER-ROI-PATCH (each patch stretched to
                             %   its own 1-99.9 pct so every cell pops regardless of FOV brightness)
dispGamma      = 0.6;        % DISPLAY-only gamma to lift dim structure. Matching stays LINEAR
                             %   (gamma would break NCC's intensity-invariance) -> affects only
                             %   what you SEE, never the scores.

% --- gallery ---
nTop           = 12;     % top pairs shown per gallery

% --- gradual-round grouping (same-cell tracks) ---
nccMergeThr    = 0.60;   % merge two ROIs into one cell if peak NCC >= this
cellSortBy     = 'confidence';  % order cells/figures by: 'confidence' (worst intra-group
                                %   NCC first = review these) or 'size' (largest group first)

% --- output / saving ---
saveOutputs     = true;              % save results .mat + membership CSV + per-cell PNGs
outDirName  = cellcfg.outDirName;
viewHalf_um     = patchHalf_um;      % DISPLAY patch half-size = COMPARE size, so what you see in
                                     %   figures/review GUIs is exactly what NCC scored. Set to a
                                     %   larger number if you want figures to show extra context.
showFirstNcells = 6;                 % open this many per-cell figures on screen; rest saved silently
savePerCellPng  = false;             % write one PNG per cell. SLOW when there are many cells
                                     %   (exportgraphics per cell). The 3 review GUIs render cells
                                     %   interactively, so leave OFF unless you want static previews.

%% ---------------- PATH SETUP ----------------
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);

outDir = '';
if saveOutputs
    outDir = fullfile(datasetPath, outDirName);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    else
        % rerun = regenerate THIS script's own products. Explicitly does NOT touch
        % roi_match_curated.* (your hand-curation) — delete that yourself if you want
        % to restart curation from a fresh auto-grouping.
        old = [dir(fullfile(outDir,'cell_*.png')); dir(fullfile(outDir,'pairs_top_*.png')); ...
               dir(fullfile(outDir,'ncc_vs_iou_scatter.png')); ...
               dir(fullfile(outDir,'roi_match_results.mat')); ...
               dir(fullfile(outDir,'roi_match_cells.csv')); ...
               dir(fullfile(outDir,'roi_match_patches.mat'))];
        for q = 1:numel(old), delete(fullfile(old(q).folder, old(q).name)); end
        if isfile(fullfile(outDir,'roi_match_curated.mat'))
            fprintf('NOTE: kept existing roi_match_curated.mat (review GUIs resume from it). Delete it to restart curation.\n');
        end
    end
end

%% ---------------- LOAD ROIs + per-FOV AVG images ----------------
sam_hits = dir(fullfile(datasetPath, '**', '*_cpSAM_output.mat'));
assert(~isempty(sam_hits), 'No *_cpSAM_output.mat found under %s', datasetPath);

% --- apply manual toss BEFORE loading ---
tossMask = false(numel(sam_hits),1);
for ff = 1:numel(sam_hits)
    tossMask(ff) = is_tossed(sam_hits(ff).folder, tossFolders);
end
if any(tossMask)
    fprintf('Tossed %d FOV(s) by manual toss list:\n', nnz(tossMask));
    tossedFolders = unique(string({sam_hits(tossMask).folder}));
    for t = 1:numel(tossedFolders), fprintf('   - %s\n', tossedFolders(t)); end
end
sam_hits = sam_hits(~tossMask);
assert(~isempty(sam_hits), 'All FOVs were tossed — nothing to load.');
nFov = numel(sam_hits);

% per-ROI table
R = struct('x',{},'y',{},'z',{},'fov',{},'roi',{},'lab',{},'area',{},'ecc',{},'sol',{}, ...
           'cx_px',{},'cy_px',{});
fov_img   = cell(nFov,1);   % MATCHING-normalized AVG image (linear, for NCC)
fov_raw   = cell(nFov,1);   % RAW AVG image (source for PER-ROI-PATCH display normalization)
fov_maskL = cell(nFov,1);   % labeled mask (uint) for ROI-outline overlay
fov_pxum  = nan(nFov,1);
fov_name  = strings(nFov,1);
fov_folder = strings(nFov,1);   % full folder path (so the reviewer can reload avg/mask)

for ff = 1:nFov
    folderPath = sam_hits(ff).folder;
    samPath    = fullfile(folderPath, sam_hits(ff).name);
    [~, fov_name(ff)] = fileparts(folderPath);
    fov_folder(ff) = string(folderPath);

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
        fov_raw{ff} = avgIm;                                        % raw (for per-patch display norm)
        lo = prctile(avgIm(:), imgClip(1)); hi = prctile(avgIm(:), imgClip(2));
        if hi <= lo, hi = lo + 1; end
        fov_img{ff} = min(max((avgIm - lo)/(hi - lo), 0), 1);       % matching-normalized (linear)
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
fprintf('Loaded %d ROIs across %d FOVs (after toss).\n', nROI, nnz(~isnan(fov_pxum)));

%% ---------------- ENUMERATE CROSS-FOV CANDIDATE PAIRS ----------------
% Distance is a cheap prefilter so we don't run normxcorr2 on all ~N^2 pairs.
% It is NOT a scoring feature. Set maxPairDist_um = Inf/[] to score every pair.
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
assert(nPair > 0, 'No cross-FOV candidate pairs — loosen maxPairDist_um / maxPairDz_um.');

pair_dist = D(sub2ind([nROI nROI], ia, ib));
pair_dz   = abs(zs(ia) - zs(ib));
pair_ncc  = nan(nPair,1);                     % offset-tolerant appearance sim
pair_iou  = nan(nPair,1);                     % mask IoU after NCC alignment

tScore = tic;
% --- VECTORIZED zero-shift NCC (replaces ~83k FFT-based normxcorr2 calls) ---
% Each patch is extracted CENTERED on its own ROI centroid, so a true same
% cell is already centered in BOTH patches; the offset search normxcorr2 did
% was only +/-1 um with these params. So precompute one unit-norm, mean-
% subtracted template vector + one mask vector per ROI, then score every pair
% as a batched dot product / mask overlap.
%   sum(unit_i .* unit_j)  ==  corr2(patch_i, patch_j)  ==  peak NCC at 0 shift.
% This drops the parallel pool entirely (no parpool startup) and runs in secs.
nPix = (round(2*patchHalf_um/targetUmPerPx) + 1)^2;
Umat = nan(nROI, nPix);          % unit-norm mean-subtracted template vectors
Mmat = false(nROI, nPix);        % ROI mask vectors (same um grid)
for r = 1:nROI
    f = fov(r); img = fov_img{f}; px = fov_pxum(f); mk = fov_maskL{f};
    t = extract_patch(img, cxpx(r), cypx(r), patchHalf_um, px, targetUmPerPx);
    if ~isempty(t)
        v = t(:) - mean(t(:)); nv = norm(v);
        if nv > 0, Umat(r,:) = (v / nv).'; end
    end
    if ~isempty(mk)
        mm = extract_patch(double(mk==labArr(r)), cxpx(r), cypx(r), patchHalf_um, px, targetUmPerPx) > 0.5;
        Mmat(r,:) = mm(:).';
    end
end

CH = 10000;                      % pair-chunk size (caps transient memory)
for s0 = 1:CH:nPair
    e0 = min(s0+CH-1, nPair);
    A = ia(s0:e0); B = ib(s0:e0);
    pair_ncc(s0:e0) = sum(Umat(A,:) .* Umat(B,:), 2);        % NaN if either patch invalid
    inter = sum(Mmat(A,:) & Mmat(B,:), 2);
    uni   = sum(Mmat(A,:) | Mmat(B,:), 2);
    iou = inter ./ uni; iou(uni == 0) = NaN;
    pair_iou(s0:e0) = iou;
end
fprintf('Pair scoring: %.1f s (%d pairs).\n', toc(tScore), nPair);

%% ---------------- MATCHED-ROI GALLERIES (the thing to examine) ----------------
% For each ranking measure, show the top nTop cross-FOV pairs as EXACTLY what
% the score compared: template A | best-aligned crop of B, each with its own
% ROI outline overlaid. Title carries NCC, IoU and the stage distance (um).
% Rank by NCC and by IoU separately so you can see what each selects.
measures = struct( ...
    'name', {'appearance NCC', 'mask IoU'}, ...
    'vals', {pair_ncc,         pair_iou}, ...
    'dir',  {'descend',        'descend'});

nShow = min(nTop, nPair);
for m = 1:numel(measures)
    [~, order] = sort(measures(m).vals, measures(m).dir, 'MissingPlacement','last');
    top = order(1:nShow);

    % --- console table ---
    fprintf('\nTop %d cross-FOV pairs by %s:\n', nShow, measures(m).name);
    fprintf('  rank   NCC    IoU   d(um) | FOV_a #roi  <->  FOV_b #roi\n');
    for r = 1:nShow
        p = top(r); i = ia(p); j = ib(p);
        fprintf('  %3d  %.3f  %.3f  %5.0f | %s #%d  <->  %s #%d\n', ...
            r, pair_ncc(p), pair_iou(p), pair_dist(p), ...
            shortName(fov_name(fov(i))), roiIdx(i), shortName(fov_name(fov(j))), roiIdx(j));
    end

    % --- patch gallery ---
    hfPair = figure('Color','w','Position',[40+40*m 40 560 90*nShow]);
    tl = tiledlayout(nShow, 2, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('Top %d cross-FOV matches by %s', nShow, measures(m).name), ...
        'Interpreter','none');
    for r = 1:nShow
        p = top(r); i = ia(p); j = ib(p);
        % show each ROI's own centroid-centered patch (exactly what was scored)
        [pa, ma] = view_patch(fov_raw{fov(i)}, fov_maskL{fov(i)}, cxpx(i), cypx(i), labArr(i), fov_pxum(fov(i)), patchHalf_um, targetUmPerPx);
        [pb, mb] = view_patch(fov_raw{fov(j)}, fov_maskL{fov(j)}, cxpx(j), cypx(j), labArr(j), fov_pxum(fov(j)), patchHalf_um, targetUmPerPx);
        pa = disp_norm_patch(pa, dispClip, dispGamma);  pb = disp_norm_patch(pb, dispClip, dispGamma);
        nexttile; if ~isempty(pa), imagesc(pa); end; axis image off; colormap(gca,gray);
        draw_outline(ma, 'y');
        title(sprintf('#%d  %s #%d', r, shortName(fov_name(fov(i))), roiIdx(i)), ...
            'Interpreter','none','FontSize',7);
        nexttile; if ~isempty(pb), imagesc(pb); end; axis image off; colormap(gca,gray);
        draw_outline(mb, 'y');
        title(sprintf('%s #%d | NCC=%.2f IoU=%.2f d=%.0fum', ...
            shortName(fov_name(fov(j))), roiIdx(j), pair_ncc(p), pair_iou(p), pair_dist(p)), ...
            'Interpreter','none','FontSize',7);
    end
    if saveOutputs
        exportgraphics(hfPair, fullfile(outDir, sprintf('pairs_top_by_%s.png', ...
            matlab.lang.makeValidName(measures(m).name))), 'Resolution',150);
    end
end

%% ---------------- NCC vs IoU scatter (top nTop IoU black) ----------------
% Each point = one cross-FOV pair. Same-cell pairs cluster at high NCC + high IoU.
nBlack = min(nTop, nPair);
[~, iouOrder] = sort(pair_iou, 'descend', 'MissingPlacement','last');
isTop = false(nPair,1); isTop(iouOrder(1:nBlack)) = true;

hfScatter = figure('Color','w','Position',[120 120 760 680]);
scatter(pair_ncc(~isTop), pair_iou(~isTop), 18, [.6 .6 .6], 'filled', 'MarkerFaceAlpha',0.4); hold on;
scatter(pair_ncc(isTop),  pair_iou(isTop),  30, 'k', 'filled');
xlabel('appearance: peak NCC'); ylabel('mask IoU (aligned)');
title(sprintf('IoU vs NCC — %d pairs (top %d IoU black)', nPair, nBlack));
legend({'rest', sprintf('top %d IoU', nBlack)}, 'Location','northwest');
grid on; box on;
if saveOutputs, exportgraphics(hfScatter, fullfile(outDir, 'ncc_vs_iou_scatter.png'), 'Resolution',150); end

%% ---------------- GRADUAL-ROUND GROUPING: same-cell tracks ----------------
% A patch can be recorded 3-6x across FOVs; pairwise scores alone won't say
% "these 5 ROIs are ONE cell." Collapse pairs into groups by greedy single-
% link agglomeration on NCC (morphology), highest score first ("gradual
% rounds"), with a hard ONE-ROI-PER-FOV constraint (a cell can't be two ROIs
% in the same FOV). Each resulting group = one physical cell + its repeats.
parent  = (1:nROI)';
uf_rank = zeros(nROI,1);
fovMask = false(nROI, nFov);                 % FOVs occupied by each group (tracked at root)
for r = 1:nROI, fovMask(r, fov(r)) = true; end

[~, ord] = sort(pair_ncc, 'descend', 'MissingPlacement','last');
nMerged = 0;
for k = 1:nPair
    p = ord(k);
    if isnan(pair_ncc(p)) || pair_ncc(p) < nccMergeThr, break; end   % sorted: rest are below thr
    ri = uf_find(parent, ia(p));
    rj = uf_find(parent, ib(p));
    if ri == rj, continue; end
    if any(fovMask(ri,:) & fovMask(rj,:)), continue; end             % share a FOV -> reject merge
    if uf_rank(ri) < uf_rank(rj), t = ri; ri = rj; rj = t; end       % union by rank
    parent(rj)    = ri;
    fovMask(ri,:) = fovMask(ri,:) | fovMask(rj,:);
    if uf_rank(ri) == uf_rank(rj), uf_rank(ri) = uf_rank(ri) + 1; end
    nMerged = nMerged + 1;
end

root = zeros(nROI,1);
for r = 1:nROI, root(r) = uf_find(parent, r); end
[~, ~, grp] = unique(root);
groupSize = accumarray(grp, 1);
multi = find(groupSize >= 2);
[~, so] = sort(groupSize(multi), 'descend');
multi = multi(so);
fprintf('\nGrouping @ NCC>=%.2f: %d merges -> %d cells, %d seen >1x (max %dx).\n', ...
    nccMergeThr, nMerged, numel(groupSize), numel(multi), max([1; groupSize(multi)]));

% --- ordered member list per cell (sorted by FOV), reused below ---
cellMembers = cell(numel(multi),1);
for gi = 1:numel(multi)
    mem = find(grp == multi(gi));
    [~, mo] = sort(fov(mem));
    cellMembers{gi} = mem(mo);
end

% --- intra-group consistency (grouping is single-link on PAIRS; here we score
%     every member against every other member of its cell). memberNCC{gi}(c) =
%     that observation's mean NCC to the rest; cellMinNCC = the weakest internal
%     pair = the honest "is this really one cell?" score; cellMedNCC = typical. ---
memberNCC  = cell(numel(multi),1);
cellMinNCC = nan(numel(multi),1);
cellMedNCC = nan(numel(multi),1);
for gi = 1:numel(multi)
    mem = cellMembers{gi};  N = numel(mem);
    G = Umat(mem,:) * Umat(mem,:)';           % N×N all-pairs NCC (diag ~1)
    offMask = ~eye(N,'logical');
    cellMinNCC(gi) = min(G(offMask), [], 'omitnan');
    cellMedNCC(gi) = median(G(offMask), 'omitnan');
    mn = nan(N,1);
    for c = 1:N
        o = true(N,1); o(c) = false;
        mn(c) = mean(G(c,o), 'omitnan');       % this member's mean NCC to the rest
    end
    memberNCC{gi} = mn;
end

% --- order cells for display/review ---
if strcmpi(cellSortBy, 'confidence')
    [~, so2] = sort(cellMinNCC, 'ascend', 'MissingPlacement','last');   % shakiest first
else
    [~, so2] = sort(cellfun(@numel, cellMembers), 'descend');           % largest first
end
multi = multi(so2);  cellMembers = cellMembers(so2);
memberNCC = memberNCC(so2);  cellMinNCC = cellMinNCC(so2);  cellMedNCC = cellMedNCC(so2);

% --- console listing (with intra-group NCC) ---
fprintf('Cells ordered by %s. minNCC = weakest internal pair (review low ones).\n', cellSortBy);
for gi = 1:numel(multi)
    mem = cellMembers{gi};
    parts = arrayfun(@(r) sprintf('%s#%d', fovLabel(fov_name(fov(r))), roiIdx(r)), ...
        mem, 'UniformOutput', false);
    fprintf('  cell %2d (%dx | min=%.2f med=%.2f): %s\n', ...
        gi, numel(mem), cellMinNCC(gi), cellMedNCC(gi), strjoin(parts, '  '));
end

%% ---------------- SAVE RESULTS (.mat + membership CSV) ----------------
if saveOutputs
    params = struct('datasetPath',char(datasetPath), ...
        'maxPairDist_um',maxPairDist_um,'maxPairDz_um',maxPairDz_um, ...
        'patchHalf_um',patchHalf_um, ...
        'targetUmPerPx',targetUmPerPx,'nccMergeThr',nccMergeThr);
    params.tossFolders = cellstr(tossFolders);

    % assign vector fields one-by-one (struct('x',vec,...) would make a struct ARRAY)
    roiT = struct(); roiT.x=xs; roiT.y=ys; roiT.z=zs; roiT.fov=fov; roiT.roi=roiIdx;
    roiT.lab=labArr; roiT.area=area; roiT.ecc=ecc; roiT.sol=sol; roiT.cx_px=cxpx; roiT.cy_px=cypx;
    pairsT = struct(); pairsT.i=ia; pairsT.j=ib; pairsT.ncc=pair_ncc; pairsT.iou=pair_iou; pairsT.dist=pair_dist;

    params.viewHalf_um = viewHalf_um;  params.imgClip = imgClip;
    params.dispClip = dispClip;  params.dispGamma = dispGamma;
    match = struct();
    match.params = params;  match.fov_name = fov_name;  match.fov_folder = fov_folder;  match.fov_pxum = fov_pxum;
    match.roi = roiT;       match.pairs = pairsT;
    match.grp = grp;        match.groupSize = groupSize;  match.multi = multi;
    match.cellMembers = cellMembers;
    match.memberNCC = memberNCC;  match.cellMinNCC = cellMinNCC;  match.cellMedNCC = cellMedNCC;
    save(fullfile(outDir,'roi_match_results.mat'), 'match', '-v7.3');

    % membership CSV: one row per (cell, observation), incl. consistency scores
    rows = {};
    for gi = 1:numel(multi)
        mem = cellMembers{gi};  mn = memberNCC{gi};
        for c = 1:numel(mem)
            r = mem(c);
            rows(end+1,:) = {gi, numel(mem), cellMinNCC(gi), cellMedNCC(gi), mn(c), ...
                             char(fov_name(fov(r))), roiIdx(r), xs(r), ys(r), zs(r)}; %#ok<SAGROW>
        end
    end
    if ~isempty(rows)
        T = cell2table(rows, 'VariableNames', ...
            {'cell_id','cell_size','cell_min_ncc','cell_med_ncc','member_ncc', ...
             'fov_name','roi_index','x_um','y_um','z_um'});
        writetable(T, fullfile(outDir,'roi_match_cells.csv'));
    end
    fprintf('Saved results (.mat) + membership CSV to %s\n', outDir);

    % ---- patch cache: lets the 3 review GUIs load instantly (no FOV reload / re-extract) ----
    % Store each ROI's DISPLAY patch + mask (viewHalf grid) and the SCORING unit-vectors
    % Umat (patchHalf grid). single precision keeps the file small.
    pAll = cell(nROI,1);  mAll = cell(nROI,1);
    for r = 1:nROI
        f = fov(r);
        pAll{r} = single(disp_norm_patch(extract_patch(fov_raw{f}, cxpx(r), cypx(r), viewHalf_um, fov_pxum(f), targetUmPerPx), dispClip, dispGamma));
        mAll{r} = extract_patch(double(fov_maskL{f}==labArr(r)), cxpx(r), cypx(r), viewHalf_um, fov_pxum(f), targetUmPerPx) > 0.5;
    end
    patchCache = struct('pAll',{pAll}, 'mAll',{mAll}, 'Umat',single(Umat), 'nROI',nROI, ...
        'viewHalf_um',viewHalf_um, 'targetUmPerPx',targetUmPerPx, 'patchHalf_um',patchHalf_um);
    save(fullfile(outDir,'roi_match_patches.mat'), 'patchCache', '-v7.3');
    fprintf('Saved patch cache (%d ROIs) to %s\n', nROI, outDir);
end

%% ---------------- ONE FIGURE PER CELL (all saved; first few shown) ----------------
% Each figure = one physical cell, its N observations tiled (<=6 per row),
% on a slightly larger view window so surrounding context is visible.
% SPEEDUP: reuse ONE hidden figure for the saved-only cells (creating 170+
% fresh figures is the slow part); only the first showFirstNcells get their
% own on-screen window.
if savePerCellPng
    tFig = tic;
    hfSave = [];
    if saveOutputs && numel(multi) > showFirstNcells
        hfSave = figure('Color','w','Visible','off');
    end
    for gi = 1:numel(multi)
        mem = cellMembers{gi};
        N   = numel(mem);
        onScreen = gi <= showFirstNcells;
        if onScreen
            hf = figure('Color','w','Visible','on');
        else
            hf = hfSave; clf(hf);
        end
        nc = min(N, 6); nr = ceil(N/nc);
        set(hf, 'Position', [100 100 210*nc 230*nr]);
        mn = memberNCC{gi};
        tl = tiledlayout(hf, nr, nc, 'TileSpacing','compact','Padding','compact');
        title(tl, sprintf('cell %d  —  %dx  —  minNCC %.2f   medNCC %.2f', ...
            gi, N, cellMinNCC(gi), cellMedNCC(gi)), 'Interpreter','none');
        for c = 1:N
            r = mem(c); f = fov(r); px = fov_pxum(f);
            pv = disp_norm_patch(extract_patch(fov_raw{f}, cxpx(r), cypx(r), viewHalf_um, px, targetUmPerPx), dispClip, dispGamma);
            mv = extract_patch(double(fov_maskL{f}==labArr(r)), cxpx(r), cypx(r), viewHalf_um, px, targetUmPerPx) > 0.5;
            ax = nexttile(tl);
            if ~isempty(pv), imagesc(ax, pv); end
            axis(ax,'image'); axis(ax,'off'); colormap(ax, gray);
            % weak members (low mean NCC to the rest) get a red outline + red title
            if mn(c) < nccMergeThr, outCol = 'r'; txtCol = 'r'; else, outCol = 'y'; txtCol = 'k'; end
            draw_outline(mv, outCol);
            title(ax, sprintf('%s\n#%d  r=%.2f', fovLabel(fov_name(f)), roiIdx(r), mn(c)), ...
                'Interpreter','none','FontSize',7, 'Color', txtCol);
        end
        minPct = cellMinNCC(gi); if ~isfinite(minPct), minPct = 0; end
        if saveOutputs
            exportgraphics(hf, fullfile(outDir, ...
                sprintf('cell_%03d_%02dx_min%02d.png', gi, N, round(100*minPct))), 'Resolution',110);
        end
    end
    if ~isempty(hfSave), close(hfSave); end
    if saveOutputs && ~isempty(multi)
        fprintf('Saved %d per-cell figures in %.1f s to %s\n', numel(multi), toc(tFig), outDir);
    end
end

%% ---------------- local helpers ----------------
function r = uf_find(parent, x)
% union-find root (no path compression; nROI is small)
r = x;
while parent(r) ~= r, r = parent(r); end
end

function tf = is_tossed(folderPath, tossList)
% true if folderPath contains any entry of tossList (case-insensitive, / == \)
tf = false;
if isempty(tossList), return; end
fp = lower(strrep(char(folderPath), '/', '\'));
tossList = string(tossList);
for t = 1:numel(tossList)
    tt = lower(strrep(char(tossList(t)), '/', '\'));
    if ~isempty(tt) && contains(fp, tt), tf = true; return; end
end
end

function s = shortName(name)
% trim long FOV folder names for tables/titles
name = char(name);
if numel(name) > 28, s = [name(1:25) '...']; else, s = name; end
end

function s = fovLabel(name)
% keep HEAD (zoom) + TAIL (condition/run) so near-identical FOV names stay
% distinguishable, e.g. roi7_1.7x_..O2_00001 vs roi7_1.8x_..roomAir_00002
name = char(name);
if numel(name) > 34, s = [name(1:20) '..' name(end-11:end)]; else, s = name; end
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

function [patch, mask] = view_patch(img, maskL, cx, cy, lab, px, half_um, target)
% Centroid-centered image patch + its ROI mask on the common um grid.
% This is EXACTLY what the zero-shift NCC compares for this ROI.
patch = extract_patch(img, cx, cy, half_um, px, target);
mask  = [];
if ~isempty(maskL)
    mask = extract_patch(double(maskL==lab), cx, cy, half_um, px, target) > 0.5;
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

function p = disp_norm_patch(p, clip, gamma)
% PER-ROI-PATCH display normalization: stretch this patch to its own clip
% percentiles then apply gamma. Display-only (never touches matching).
if isempty(p), return; end
lo = prctile(p(:), clip(1)); hi = prctile(p(:), clip(2)); if hi <= lo, hi = lo + 1; end
p = min(max((p - lo)/(hi - lo), 0), 1); if gamma ~= 1, p = p .^ gamma; end
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
