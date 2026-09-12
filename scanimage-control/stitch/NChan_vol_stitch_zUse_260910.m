function NChan_vol_stitch_zUse_260910()
clc;

rawDir = "D:\Data\RZ\test260911";
refChID = 3;        % registration reference. USE THE VESSEL CHANNEL -- on 260804
                    % ch3 gives NCC 0.94 while ch1 gives 0.27-0.49.
chUse   = [];       % channels to WRITE OUT. [] = every ch<N>\avgz folder found,
                    % same convention as tile_deinterleave_avgz_260910. Or name
                    % them explicitly, e.g. [1 3]. refChID is always included.
                    % Registration is ALWAYS done on refChID alone; this only
                    % controls which channels get stitched with those coords.
zUse_um = []; zUseProj = [];
%zUse_um = [-300 20]; zUseProj = [];
rescaled_um = 1;
stageSource = "auto";

%%
step_um         = 400; um_per_px_guess = 1.4815;

placeInit    = "stage"; row_increases_up = true;

flipX = false; flipY = false;

SearchRadiusPx   = 120;
SearchPerpPx     = 40;
HighPassSigmaPx  = 0;
TemplateFrac     = 0.60;
MinTemplatePx    = 24;
NccMin           = 0.45;
ExclRadiusPx     = 15;
MarginMin        = 0.02;
subpixel         = true;
ScaleWarnFrac    = 0.02;

PriorWeight      = 0.05;
DarkContrast     = 0.025;
DarkPriorBoost   = 200;
HuberK           = 2.5;
IrlsIters        = 6;

useLinearBlending = true;
rowsPerStrip      = 64;
nRows_tile = [];  nCols_tile = [];

writeStitched = true;
compareTo     = "";

refAvgzDir   = fullfile(rawDir, sprintf("ch%d", refChID),   "avgz");

% Discover channels from the folders that actually exist, rather than assuming
% two. The deinterleaver writes one ch<N>\avgz per saved channel and the count
% is not known until a header is read, so scanning is the only honest way.
chFound = discover_channels(rawDir);
if isempty(chFound)
    error("No ch<N>%savgz folders under %s -- run the deinterleaver first.", filesep, rawDir);
end
if isempty(chUse)
    chIDs = chFound;
else
    chIDs = intersect(chUse(:).', chFound, 'stable');
    gone  = setdiff(chUse(:).', chFound);
    if ~isempty(gone)
        warning("requested channel(s) [%s] have no avgz folder -- skipped.", num2str(gone));
    end
end
if ~ismember(refChID, chFound)
    error("refChID = %d has no avgz folder (found [%s]).", refChID, num2str(chFound));
end
chIDs = unique([refChID, chIDs], 'stable');   % ref always written
otherChIDs   = setdiff(chIDs, refChID, 'stable');
otherAvgzDirs = arrayfun(@(c) fullfile(rawDir, sprintf("ch%d", c), "avgz"), ...
                         otherChIDs, 'uni', 0);
fprintf("[chan] found [%s], registering on ch%d, writing [%s]\n", ...
        num2str(chFound), refChID, num2str(chIDs));
outDir = fullfile(rawDir, "matlab_stitch");
if ~isfolder(outDir), mkdir(outDir); end

coordMat = fullfile(outDir, sprintf("stitch_coords_from_ch%d.mat", refChID));

if ~isfolder(refAvgzDir), error("Missing folder: %s", refAvgzDir); end
for m = 1:numel(otherAvgzDirs)
    if ~isfolder(otherAvgzDirs{m})
        error("Missing folder: %s", otherAvgzDirs{m});
    end
end

refFilesAll0 = list_tifs(refAvgzDir);
if isempty(refFilesAll0), error("No TIFFs in %s", refAvgzDir); end

N0 = numel(refFilesAll0);
okCR = false(N0,1); okTI = false(N0,1);
col0 = nan(N0,1); row0 = nan(N0,1); tileIdx0 = nan(N0,1);
for i = 1:N0
    [c,r,tf] = parse_colrow_from_name(refFilesAll0(i));
    if tf, okCR(i) = true; col0(i) = c; row0(i) = r; end
    [k,tf2] = parse_tileidx_from_name(refFilesAll0(i));
    if tf2, okTI(i) = true; tileIdx0(i) = k; end
end
nCR = nnz(okCR); nTI = nnz(okTI);
if     nCR>0 && nTI==0, mode = "colrow";
elseif nTI>0 && nCR==0, mode = "tileidx";
elseif nCR==0 && nTI==0
    error("No files match either naming: need colXX/rowYY OR tile_### in filenames.");
else
    error("Mixed naming in %s -- remove stray files / use one convention.", refAvgzDir);
end
fprintf("[mode] %s\n", mode);

[tileFileRef, tileFilesOther, tileAtRC, nRows, nCols, c_of_t, rB_of_t] = ...
    build_order(mode, refFilesAll0, okCR, okTI, col0, row0, tileIdx0, ...
                otherAvgzDirs, otherChIDs, row_increases_up, nRows_tile, nCols_tile);

% tileFiles{k} is the tile list for channel chIDs(k); the ref sits at index 1.
tileFiles = [{tileFileRef}, tileFilesOther];
nTiles = numel(tileFileRef);
fprintf("[grid] nRows=%d nCols=%d nTiles=%d\n", nRows, nCols, nTiles);

[H0, W0, inClass] = read_first_page_size_class(tileFileRef(1));
nZ0 = count_tiff_pages_file(tileFileRef(1));
fprintf("[tile] size=%dx%d, nZ=%d, class=%s\n", H0, W0, nZ0, inClass);

ACQ = read_acq_params(rawDir, tileFileRef(1));
if ACQ.ok
    fprintf("[acq]  %s\n", ACQ.src);
    fprintf("[acq]  zoom %g, FOV %.2f um, %.4f um/px, channels %s, %g slices x %g frames\n", ...
            ACQ.zoom, ACQ.fovUm, ACQ.umPerPx, mat2str(ACQ.channelSave(:).'), ...
            ACQ.nSlices, ACQ.framesPerSlice);
    if ~isempty(ACQ.channelSave) && ~ismember(refChID, ACQ.channelSave)
        warning("refChID = %d is not in the acquired channelSave %s.", ...
                refChID, mat2str(ACQ.channelSave(:).'));
    end
    if isfinite(ACQ.umPerPx) && abs(ACQ.umPerPx - um_per_px_guess)/ACQ.umPerPx > 1e-3
        fprintf("[acq]  seeding um_per_px from the header: %.4f (typed value was %.4f)\n", ...
                ACQ.umPerPx, um_per_px_guess);
        um_per_px_guess = ACQ.umPerPx;
    end
else
    fprintf("[acq]  no raw tile in %s -- zoom/z/channels unknown, using the typed settings\n", rawDir);
end

if isempty(zUse_um)
    if isempty(zUseProj)
        zUse = 1:nZ0;
        fprintf("[zUse] all %d pages (no range set)\n", nZ0);
    else
        zUse = zUseProj(:).';
        zUse = zUse(zUse >= 1 & zUse <= nZ0);
        assert(~isempty(zUse), "zUseProj is empty after clipping to 1..%d.", nZ0);
        fprintf("[zUse] pages %d..%d (%d of %d) from zUseProj -- PAGE INDICES, " + ...
                "depth in um unknown\n", min(zUse), max(zUse), numel(zUse), nZ0);
    end
else
    if ACQ.ok && ~isempty(ACQ.zs)
        zs = ACQ.zs;  zsrc = ACQ.src;
    else
        [zs, zsrc] = find_stack_zs(rawDir, tileFileRef(1));
    end
    assert(~isempty(zs), ['zUse_um is in MICRONS and needs SI.hStackManager.zs, which ' ...
        'lives only in a raw tile -- none found in %s. Put one raw tile there, or ' ...
        'use zUseProj (page indices) instead.'], rawDir);
    assert(numel(zs) == nZ0, ['the raw header lists %d slices but the avgz tile has ' ...
        '%d pages -- they are not one page per slice, so a micron range cannot be ' ...
        'mapped onto pages.'], numel(zs), nZ0);
    lo = min(zUse_um); hi = max(zUse_um);
    zUse = find(zs >= lo & zs <= hi).';
    assert(~isempty(zUse), ['no slice falls in [%g %g] um; the stack spans %g to %g ' ...
        'um in steps of %g.'], lo, hi, zs(1), zs(end), median(diff(zs)));
    fprintf("[zUse] %g to %g um -> pages %d..%d (%d of %d), z %g to %g um  [zs from %s]\n", ...
            lo, hi, min(zUse), max(zUse), numel(zUse), nZ0, ...
            zs(min(zUse)), zs(max(zUse)), zsrc);
end

for t = 1:nTiles
    [Hk,Wk,clsK] = read_first_page_size_class(tileFileRef(t));
    if Hk~=H0 || Wk~=W0 || clsK~=string(inClass)
        error("Tile size/class mismatch: %s", tileFileRef(t));
    end
end

sumProj = cell(nTiles,1);
for t = 1:nTiles
    tt = Tiff(tileFileRef(t),'r');
    S = zeros(H0,W0,'double');
    for iz = zUse, tt.setDirectory(iz); S = S + double(tt.read()); end
    tt.close();
    if HighPassSigmaPx > 0
        S = S - imgaussfilt(S, HighPassSigmaPx);
    end
    sumProj{t} = S;
end
if HighPassSigmaPx > 0
    fprintf("[proj] high-pass sigma %g px before registration\n", HighPassSigmaPx);
end

[stageX, stageY, haveStage, stageSrc] = parse_stage_tokens(tileFileRef, rawDir, stageSource);
step_px = step_um / um_per_px_guess;

if placeInit == "stage" && haveStage
    x0 = (stageX - min(stageX)) / um_per_px_guess;
    y0 = (max(stageY) - stageY)  / um_per_px_guess;
    fprintf("[init] stage from %s, %.4f um/px -> step ~%.1f px\n", ...
            stageSrc, um_per_px_guess, step_px);
else
    if placeInit == "stage"
        warning("no _x###_y### tokens in the filenames -- falling back to the grid lattice");
    end
    x0 = zeros(nTiles,1); y0 = zeros(nTiles,1);
    for t = 1:nTiles
        rTop = nRows - rB_of_t(t) + 1;
        x0(t) = (c_of_t(t)-1) * step_px;
        y0(t) = (rTop-1)      * step_px;
    end
    fprintf("[init] uniform grid, step %.1f px\n", step_px);
end
if flipX, x0 = max(x0) - x0; end
if flipY, y0 = max(y0) - y0; end

idxGrid = zeros(nRows,nCols);
for t = 1:nTiles, idxGrid(rB_of_t(t), c_of_t(t)) = t; end

E = struct('i',{},'j',{},'dxG',{},'dyG',{},'axis',{});
for rB = 1:nRows
    for c = 1:nCols
        i = idxGrid(rB,c);
        if c < nCols
            j = idxGrid(rB,c+1);
            E(end+1) = struct('i',i,'j',j,'dxG',x0(j)-x0(i),'dyG',y0(j)-y0(i),'axis','x');
        end
        if rB < nRows
            j = idxGrid(rB+1,c);
            E(end+1) = struct('i',i,'j',j,'dxG',x0(j)-x0(i),'dyG',y0(j)-y0(i),'axis','y');
        end
    end
end
nE = numel(E);
fprintf("[reg] %d edges, search +/-%d px along / +/-%d across\n", nE, SearchRadiusPx, SearchPerpPx);

dxE = zeros(nE,1); dyE = zeros(nE,1); okE = false(nE,1);
ncc = nan(nE,1); margin = nan(nE,1); why = strings(nE,1);
t0 = tic;
for e = 1:nE
    [dxE(e), dyE(e), okE(e), ncc(e), margin(e), why(e)] = refine_edge_ncc( ...
        sumProj{E(e).i}, sumProj{E(e).j}, E(e).dxG, E(e).dyG, ...
        SearchRadiusPx, SearchPerpPx, TemplateFrac, MinTemplatePx, ...
        NccMin, ExclRadiusPx, MarginMin, subpixel);
    if ~okE(e), dxE(e) = E(e).dxG; dyE(e) = E(e).dyG; end
end
fprintf("[reg] %d/%d edges accepted (%.0f s)\n", nnz(okE), nE, toc(t0));
qn = prctile(ncc, [5 25 50 75 95]);
fprintf("[reg] NCC over all %d edges: p5 %.2f  p25 %.2f  med %.2f  p75 %.2f  p95 %.2f\n", ...
        nE, qn(1), qn(2), qn(3), qn(4), qn(5));
if nnz(~okE) > 0
    [u,~,ix] = unique(why(~okE));
    cnt = accumarray(ix,1);
    for k = 1:numel(u), fprintf("      rejected %3d x  %s\n", cnt(k), u(k)); end
end

isx = strcmp({E.axis},'x').';
resid = hypot(dxE - [E.dxG].', dyE - [E.dyG].');
fprintf("[reg] residual vs guess: median %.1f px, max %.1f px (accepted edges)\n", ...
        median(resid(okE)), max(resid(okE)));

if haveStage
    si = [E.i].'; sj = [E.j].';
    dsx = abs(stageX(sj) - stageX(si));  dsy = abs(stageY(sj) - stageY(si));
    uX = dsx(okE &  isx) ./ abs(dxE(okE &  isx));
    uY = dsy(okE & ~isx) ./ abs(dyE(okE & ~isx));
    uAll = [uX; uY];
    if ~isempty(uAll)
        fprintf("[scale] measured %.4f um/px (x %.4f, y %.4f) vs guess %.4f\n", ...
                median(uAll), median(uX), median(uY), um_per_px_guess);
        if abs(median(uAll)/um_per_px_guess - 1) > ScaleWarnFrac
            warning(['um_per_px_guess is %.1f%% off. Set um_per_px_guess = %.4f ' ...
                     'for a tighter search next run.'], ...
                     100*(median(uAll)/um_per_px_guess - 1), median(uAll));
        end
    end
end

if haveStage && exist('uX','var') && ~isempty(uX) && ~isempty(uY)
    sxPx = median(uX);  syPx = median(uY);
    x0 = (stageX - min(stageX)) / sxPx;
    y0 = (max(stageY) - stageY) / syPx;
    if flipX, x0 = max(x0) - x0; end
    if flipY, y0 = max(y0) - y0; end
    fprintf('[prior] rebuilt at the measured scale: x %.4f, y %.4f um/px\n', sxPx, syPx);
end

contentIdx = zeros(nTiles,1);
for t = 1:nTiles
    q = prctile(sumProj{t}(:), [50 99.5]);
    contentIdx(t) = (q(2) - q(1)) / max(1, abs(q(1)));
end
isDark = contentIdx < DarkContrast;
priorW = repmat(PriorWeight, nTiles, 1);
priorW(isDark) = PriorWeight * DarkPriorBoost;
fprintf('[dark] contrast index: min %.3f  med %.3f  max %.3f   -> %d tile(s) below %.3f pinned\n', ...
        min(contentIdx), median(contentIdx), max(contentIdx), nnz(isDark), DarkContrast);
for t = find(isDark).'
    fprintf('         %s  (%.3f)\n', tileFileRef(t), contentIdx(t));
end

loop_closure_report(E, dxE, dyE, okE, idxGrid, nRows, nCols);

ii = [E(okE).i].'; jj = [E(okE).j].';
assert(~isempty(ii), 'every edge was rejected -- nothing to solve');

report_connectivity(nTiles, ii, jj, tileFileRef);

w = 1 + 19 * max(0, (ncc(okE) - NccMin) / max(eps, 1 - NccMin));

[x, y, solveInfo] = solve_robust(nTiles, ii, jj, dxE(okE), dyE(okE), w, ...
                                 x0, y0, priorW, HuberK, IrlsIters);
x = x - min(x);  y = y - min(y);
fprintf("[solve] IRLS %d iters, %d/%d edges downweighted below half\n", ...
        IrlsIters, solveInfo.nDown, numel(ii));

pr = hypot((x(jj)-x(ii)) - dxE(okE), (y(jj)-y(ii)) - dyE(okE));
fprintf("[solve] post-solve edge residual: median %.2f px, 95th %.2f, max %.2f\n", ...
        median(pr), prctile(pr,95), max(pr));

Hn = H0; Wn = W0;
if ~exist('sxPx','var'), sxPx = um_per_px_guess; syPx = um_per_px_guess; end
x_raw = x;  y_raw = y;
if ~isempty(rescaled_um) && rescaled_um > 0 && haveStage && exist('sxPx','var')
    fx_s = sxPx / rescaled_um;
    fy_s = syPx / rescaled_um;
    Wn = max(1, round(W0 * fx_s));
    Hn = max(1, round(H0 * fy_s));
    x  = x * fx_s;   y = y * fy_s;
    x_um = x;  y_um = y;
    fprintf(['[scale] mosaic resampled to %g um/px: tile %dx%d -> %dx%d ' ...
             '(x %.4f, y %.4f um/px in)\n'], rescaled_um, H0, W0, Hn, Wn, sxPx, syPx);
elseif ~isempty(rescaled_um) && rescaled_um > 0
    warning(['rescaled_um is set but no measured scale is available (no stage ' ...
             'positions, or no accepted edges) -- writing at native pixel size.']);
end

if ~exist('x_um','var')
    x_um = x_raw; y_um = y_raw;
end
outW = ceil(max(x) + Wn);  outH = ceil(max(y) + Hn);
writtenSpace = "raw";
if ~isempty(rescaled_um) && (Hn ~= H0 || Wn ~= W0), writtenSpace = "rescaled"; end

fprintf("[stitch] mosaic %dx%d\n", outH, outW);

% tileFileRef is kept because stitch_manual_gui_260910.m:198 reads it.
% tileFiles/chIDs are the N-channel form; tileFileOther is a legacy alias for
% the first non-ref channel so any older consumer still loads.
if numel(tileFiles) > 1, tileFileOther = tileFiles{2}; else, tileFileOther = strings(0,1); end
save(coordMat, 'x','y','tileFileRef','tileFileOther','tileFiles','chIDs', ...
     'refChID','c_of_t','rB_of_t', ...
     'nRows','nCols','zUse','zUse_um','zUseProj','stageSrc', ...
     'x_raw','y_raw','H0','W0', ...
     'x_um','y_um','Hn','Wn','rescaled_um','writtenSpace', ...
     'outH','outW','step_um','um_per_px_guess','sxPx','syPx', ...
     'row_increases_up','flipX','flipY','mode','E','dxE','dyE','okE','ncc','margin','pr');
fprintf("[save] %s\n", coordMat);

cur = fullfile(outDir, compareTo);
if compareTo ~= "" && isfile(cur)
    compare_to_curated(cur, x, y, tileFileRef, W0, H0);
end

qc_figure(outDir, E, dxE, dyE, okE, ncc, pr, isx, x, y, x0, y0, W0, H0, refChID);

if ~writeStitched
    fprintf("[skip] writeStitched=false -- coords and QC only\n"); return
end
rawTag = sprintf("raw_%gum", round(mean([sxPx syPx]), 1));
SPACES = struct('tag', {rawTag}, 'x', {x_raw}, 'y', {y_raw}, 'H', {H0}, 'W', {W0});
if ~isempty(rescaled_um) && (Hn ~= H0 || Wn ~= W0)
    SPACES(2) = struct('tag', sprintf("rescaled_%gum", rescaled_um), ...
                       'x', x_um, 'y', y_um, 'H', Hn, 'W', Wn);
end

for sp = 1:numel(SPACES)
    SP = SPACES(sp);
    xs = SP.x - min(SP.x);  ys = SP.y - min(SP.y);
    oW = ceil(max(xs) + SP.W);  oH = ceil(max(ys) + SP.H);
    fprintf("[write] %s space: tile %dx%d -> mosaic %dx%d\n", SP.tag, SP.H, SP.W, oH, oW);

    if useLinearBlending, wTile = linear_blend_weights(SP.H, SP.W);
    else,                 wTile = ones(SP.H, SP.W, 'double'); end

    place = struct([]);
    for t = 1:nTiles
        xInt = floor(xs(t)); yInt = floor(ys(t));
        fx = xs(t) - xInt;   fy = ys(t) - yInt;
        place(t).fx = fx;  place(t).fy = fy;
        place(t).rr = (yInt+1):(yInt+SP.H);
        place(t).cc = (xInt+1):(xInt+SP.W);
        if subpixel && (abs(fx)>1e-6 || abs(fy)>1e-6)
            place(t).w2 = shift_bilinear(wTile, fx, fy);
        else
            place(t).w2 = wTile;
        end
    end

    for ci = 1:numel(chIDs)
        cc     = chIDs(ci);
        tfiles = tileFiles{ci};
        fn = fullfile(outDir, sprintf("stitched_ch%d_avgz_%s.tif", cc, SP.tag));
        if isfile(fn), delete(fn); end
        write_stitched_from_tilefiles(fn, tfiles, place, oH, oW, SP.H, SP.W, zUse, ...
            inClass, rowsPerStrip, ...
            sprintf("ch%d, coords from ch%d, %s space; zUse=[%s] um = pages %d..%d; mode=%s", ...
                    cc, refChID, SP.tag, num2str(zUse_um), min(zUse), max(zUse), mode));
        fprintf("[DONE] %s\n", fn);
    end
end
end

function [dx, dy, ok, peak, margin, why] = refine_edge_ncc(A, B, dxG, dyG, ...
        Ralong, Rperp, tfrac, tmin, nccMin, exclR, marginMin, subpix)

dx = dxG; dy = dyG; ok = false; peak = NaN; margin = NaN; why = "";

swapped = false;
horiz = abs(dxG) >= abs(dyG);
if (horiz && dxG < 0) || (~horiz && dyG < 0)
    [A, B] = deal(B, A);  dxG = -dxG;  dyG = -dyG;  swapped = true;
end

H = size(A,1); W = size(A,2);
if horiz
    ovGuess = W - round(dxG);
    tw = max(tmin, round(tfrac * ovGuess));
    if tw >= W || ovGuess <= tmin, why = "overlap too small"; return; end
    cs0 = max(1, W - tw - Ralong + 1);
    As  = A(:, cs0:W);
    vM  = min(Rperp, floor((H - tmin)/2));
    rT  = (1+vM):(H-vM);
    T   = B(rT, 1:tw);
else
    ovGuess = H - round(dyG);
    th = max(tmin, round(tfrac * ovGuess));
    if th >= H || ovGuess <= tmin, why = "overlap too small"; return; end
    rs0 = max(1, H - th - Ralong + 1);
    As  = A(rs0:H, :);
    hM  = min(Rperp, floor((W - tmin)/2));
    cT  = (1+hM):(W-hM);
    T   = B(1:th, cT);
end

if any(size(T) < tmin) || any(size(As) < size(T)), why = "patch too small"; return; end
if std(T(:)) < eps || std(As(:)) < eps, why = "flat patch"; return; end

C  = normxcorr2(T, As);
[th_, tw_] = size(T);
[ah_, aw_] = size(As);
Cv = C(th_:ah_, tw_:aw_);
if isempty(Cv), why = "no valid overlap"; return; end

[peak, k] = max(Cv(:));
[fr, fc]  = ind2sub(size(Cv), k);

msk = true(size(Cv));
msk(max(1,fr-exclR):min(end,fr+exclR), max(1,fc-exclR):min(end,fc+exclR)) = false;
if nnz(msk) > 10, margin = peak - max(Cv(msk)); else, margin = Inf; end

if subpix
    [sr, sc] = subpix2(Cv, fr, fc);
else
    sr = 0; sc = 0;
end

if horiz
    dyN = (fr + sr) - rT(1);
    dxN = cs0 + (fc + sc) - 2;
else
    dxN = (fc + sc) - cT(1);
    dyN = rs0 + (fr + sr) - 2;
end

if abs(dxN - dxG) > Ralong + 1 || abs(dyN - dyG) > max(Ralong,Rperp) + 1
    why = "peak outside the search window"; return
end
if ~isfinite(peak) || peak < nccMin, why = sprintf("ncc < %.2f", nccMin); return; end
if ~isfinite(margin) || margin < marginMin
    why = sprintf("not unique (margin < %.3f)", marginMin); return
end

dx = dxN; dy = dyN; ok = true;
if swapped, dx = -dx; dy = -dy; end
end

function [sr, sc] = subpix2(C, r, c)
sr = 0; sc = 0;
if r>1 && r<size(C,1)
    d = C(r-1,c) - 2*C(r,c) + C(r+1,c);
    if abs(d) > 1e-12, sr = 0.5*(C(r-1,c) - C(r+1,c))/d; end
end
if c>1 && c<size(C,2)
    d = C(r,c-1) - 2*C(r,c) + C(r,c+1);
    if abs(d) > 1e-12, sc = 0.5*(C(r,c-1) - C(r,c+1))/d; end
end
sr = max(-1,min(1,sr));  sc = max(-1,min(1,sc));
end

function comp = report_connectivity(n, ii, jj, names)
lab = zeros(n,1); comp = 0;
adj = cell(n,1);
for e = 1:numel(ii)
    adj{ii(e)}(end+1) = jj(e);  adj{jj(e)}(end+1) = ii(e);
end
for s = 1:n
    if lab(s), continue; end
    comp = comp + 1;  lab(s) = comp;  stack = s;
    while ~isempty(stack)
        v = stack(end); stack(end) = [];
        for u = adj{v}
            if ~lab(u), lab(u) = comp; stack(end+1) = u; end
        end
    end
end
sz = accumarray(lab, 1);
[big, ~] = max(sz);
iso = find(cellfun(@isempty, adj));
fprintf('[graph] %d component(s); largest holds %d of %d tiles\n', comp, big, n);
if ~isempty(iso)
    fprintf('[graph] %d tile(s) with NO accepted edge -- pinned to the stage:\n', numel(iso));
    for k = 1:numel(iso)
        fprintf('          %s\n', names(iso(k)));
    end
end
end

function report_connectivity_old(n, ii, jj)
seen = false(n,1); seen(1) = true;
adj = cell(n,1);
for e = 1:numel(ii)
    adj{ii(e)}(end+1) = jj(e);  adj{jj(e)}(end+1) = ii(e);
end
stack = 1;
while ~isempty(stack)
    v = stack(end); stack(end) = [];
    for u = adj{v}
        if ~seen(u), seen(u) = true; stack(end+1) = u; end
    end
end
if ~all(seen)
    fprintf(['[graph] %d of %d tiles are NOT edge-connected to tile 1 -- those ' ...
             'rest on the stage prior alone\n'], nnz(~seen), n);
else
    fprintf('[graph] all %d tiles edge-connected\n', n);
end
end

function [x, y, info] = solve_robust(n, ii, jj, dx, dy, w0, px, py, priorW, huberK, nIter)
px = px(:) - min(px);  py = py(:) - min(py);
w  = w0(:);
info = struct('nDown', 0, 'sigma', NaN);
for it = 1:nIter
    [x, y] = solve_once(n, ii, jj, dx, dy, w, px, py, priorW);
    rx = (x(jj)-x(ii)) - dx;   ry = (y(jj)-y(ii)) - dy;
    r  = hypot(rx, ry);
    sig = 1.4826 * median(abs(r - median(r)));
    if sig < 1e-6, sig = max(1e-6, median(r)); end
    f = min(1, huberK * sig ./ max(r, eps));
    w = w0(:) .* f;
    info.sigma = sig;
end
info.nDown = nnz(f < 0.5);
end

function [x, y] = solve_once(n, ii, jj, dx, dy, w, px, py, priorW)
m = numel(dx);
rows = 2*m + 2*n;
A = spalloc(rows, 2*n, 4*m + 2*n);
b = zeros(rows,1);  ww = zeros(rows,1);
r = 0;
for e = 1:m
    r = r + 1;
    A(r, 2*jj(e)-1) = 1;  A(r, 2*ii(e)-1) = -1;  b(r) = dx(e);  ww(r) = w(e);
    r = r + 1;
    A(r, 2*jj(e))   = 1;  A(r, 2*ii(e))   = -1;  b(r) = dy(e);  ww(r) = w(e);
end
for t = 1:n
    r = r + 1;  A(r, 2*t-1) = 1;  b(r) = px(t);  ww(r) = priorW(t);
    r = r + 1;  A(r, 2*t)   = 1;  b(r) = py(t);  ww(r) = priorW(t);
end
W = spdiags(sqrt(ww), 0, rows, rows);
p = lsqminnorm(W*A, W*b);
x = p(1:2:end);  y = p(2:2:end);
end

function loop_closure_report(E, dxE, dyE, okE, idxGrid, nRows, nCols)
key = containers.Map('KeyType','char','ValueType','double');
for e = 1:numel(E)
    key(sprintf('%d_%d', E(e).i, E(e).j)) = e;
end
err = [];
for rB = 1:nRows-1
    for c = 1:nCols-1
        A = idxGrid(rB,c);  B = idxGrid(rB,c+1);
        C = idxGrid(rB+1,c); D = idxGrid(rB+1,c+1);
        kAB = sprintf('%d_%d',A,B); kCD = sprintf('%d_%d',C,D);
        kAC = sprintf('%d_%d',A,C); kBD = sprintf('%d_%d',B,D);
        if ~(isKey(key,kAB) && isKey(key,kCD) && isKey(key,kAC) && isKey(key,kBD)), continue; end
        e1 = key(kAB); e2 = key(kBD); e3 = key(kAC); e4 = key(kCD);
        if ~all(okE([e1 e2 e3 e4])), continue; end
        ex = (dxE(e1) + dxE(e2)) - (dxE(e3) + dxE(e4));
        ey = (dyE(e1) + dyE(e2)) - (dyE(e3) + dyE(e4));
        err(end+1) = hypot(ex, ey);
    end
end
if isempty(err)
    fprintf("[loop] no closed 2x2 cycle has all four edges accepted\n");
else
    fprintf("[loop] %d closed cycles: median %.2f px, 95th %.2f, max %.2f\n", ...
            numel(err), median(err), prctile(err,95), max(err));
end
end

function [x, y] = solve_global_positions_weighted(n, ii, jj, dx, dy, w)
m = numel(dx);
A = spalloc(2*m + 2, 2*n, 4*m + 2);
b = zeros(2*m + 2, 1);  ww = zeros(2*m + 2, 1);
row = 0;
for e = 1:m
    row = row + 1;
    A(row, 2*jj(e)-1) = 1;  A(row, 2*ii(e)-1) = -1;  b(row) = dx(e);  ww(row) = w(e);
    row = row + 1;
    A(row, 2*jj(e))   = 1;  A(row, 2*ii(e))   = -1;  b(row) = dy(e);  ww(row) = w(e);
end
A(end-1, 1) = 1; b(end-1) = 0; ww(end-1) = 50;
A(end,   2) = 1; b(end)   = 0; ww(end)   = 50;
Wmat = spdiags(sqrt(ww), 0, numel(ww), numel(ww));
p = lsqminnorm(Wmat*A, Wmat*b);
x = p(1:2:end);  y = p(2:2:end);
end

function compare_to_curated(curFile, x, y, tileFileRef, W0, H0)
C = load(curFile);
if ~isfield(C,'x') || ~isfield(C,'files'), fprintf("[check] %s has no x/files\n", curFile); return; end
[~, bnA] = cellfun(@(s) fileparts(char(s)), num2cell(tileFileRef), 'uni', 0);
[~, bnC] = cellfun(@(s) fileparts(char(s)), num2cell(C.files),     'uni', 0);
[tf, loc] = ismember(bnA, bnC);
if nnz(tf) < 3, fprintf("[check] could not match tiles to the curated set\n"); return; end
dx = (x(tf) - C.x(loc(tf)));  dy = (y(tf) - C.y(loc(tf)));
dx = dx - median(dx);  dy = dy - median(dy);
d  = hypot(dx, dy);
fprintf(['[check] vs curated (%d tiles, common translation removed):\n' ...
         '        median %.1f px, 95th %.1f px, max %.1f px  (tile is %dx%d)\n'], ...
         nnz(tf), median(d), prctile(d,95), max(d), H0, W0);
end

function qc_figure(outDir, E, dxE, dyE, okE, ncc, pr, isx, x, y, x0, y0, W0, H0, refCh)
COL = [42 120 214; 235 104 52]/255;
INK2 = [0.32 0.32 0.31];
try, opengl('software'); catch, end
f = figure('Color','w','Units','pixels','Position',[40 40 1400 760], 'Visible','off');
tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl,1);
histogram(ax, ncc(okE), 24, 'FaceColor',COL(1,:),'EdgeColor','none'); grid(ax,'on');
xlabel(ax,'peak NCC'); ylabel(ax,'edges'); title(ax,'edge quality','FontWeight','normal');

ax = nexttile(tl,2); hold(ax,'on');
r = hypot(dxE-[E.dxG].', dyE-[E.dyG].');
histogram(ax, r(okE &  isx), 24, 'FaceColor',COL(1,:),'EdgeColor','none','DisplayName','x');
histogram(ax, r(okE & ~isx), 24, 'FaceColor',COL(2,:),'EdgeColor','none','DisplayName','y');
grid(ax,'on'); xlabel(ax,'|refined - guess| (px)'); ylabel(ax,'edges');
title(ax,'how wrong the guess was','FontWeight','normal');
lg=legend(ax); lg.Box='off'; hold(ax,'off');

ax = nexttile(tl,3);
histogram(ax, pr, 24, 'FaceColor',COL(1,:),'EdgeColor','none'); grid(ax,'on');
xlabel(ax,'post-solve edge residual (px)'); ylabel(ax,'edges');
title(ax,'solve consistency','FontWeight','normal');

ax = nexttile(tl,4); hold(ax,'on');
plot(ax, x0-min(x0), y0-min(y0), 'o', 'Color',[.7 .7 .7], 'MarkerSize',4, 'DisplayName','initial');
plot(ax, x, y, '.', 'Color',COL(1,:), 'MarkerSize',12, 'DisplayName','solved');
set(ax,'YDir','reverse'); axis(ax,'equal'); grid(ax,'on');
xlabel(ax,'x (px)'); ylabel(ax,'y (px)'); title(ax,'tile centres','FontWeight','normal');
lg=legend(ax,'Location','best'); lg.Box='off'; hold(ax,'off');

ax = nexttile(tl,5); hold(ax,'on');
for e = 1:numel(E)
    c = COL(2,:); if okE(e), c = COL(1,:); end
    plot(ax, [x(E(e).i) x(E(e).j)] + W0/2, [y(E(e).i) y(E(e).j)] + H0/2, '-', 'Color', c, 'LineWidth',1);
end
set(ax,'YDir','reverse'); axis(ax,'equal'); grid(ax,'on');
xlabel(ax,'x (px)'); ylabel(ax,'y (px)');
title(ax, sprintf('edges: %d accepted (blue), %d rejected (orange)', nnz(okE), nnz(~okE)), ...
      'FontWeight','normal'); hold(ax,'off');

ax = nexttile(tl,6); hold(ax,'on');
sx = arrayfun(@(e) e.dxG, E).';  sy = arrayfun(@(e) e.dyG, E).';
scatter(ax, hypot(sx,sy), hypot(dxE,dyE), 14, COL(1,:), 'filled', 'MarkerFaceAlpha',0.5);
xl = xlim(ax); plot(ax, xl, xl, '-', 'Color', INK2); grid(ax,'on');
xlabel(ax,'guessed step (px)'); ylabel(ax,'refined step (px)');
title(ax,'guess vs measured','FontWeight','normal'); hold(ax,'off');

title(tl, sprintf('auto-stitch QC  --  ch%d  --  %d/%d edges accepted', refCh, nnz(okE), numel(E)), ...
      'FontWeight','bold','Interpreter','none');
base = fullfile(outDir, sprintf('autostitch_qc_ch%d', refCh));
exportgraphics(f, base + ".png", 'Resolution', 170, 'BackgroundColor','white');
close(f);
fprintf("[qc]   %s.png\n", base);
end

function [tileFileRef, tileFilesOther, tileAtRC, nRows, nCols, c_of_t, rB_of_t] = ...
    build_order(mode, refFilesAll0, okCR, okTI, col0, row0, tileIdx0, ...
                otherAvgzDirs, otherChIDs, row_increases_up, nRows_tile, nCols_tile)
% otherAvgzDirs is a CELL of avgz folders, one per non-reference channel, and
% tileFilesOther comes back as a matching cell of tile lists. The pairing rule is
% identical for every channel, so it is applied in a loop rather than written
% once for a hard-coded "other".
if mode == "colrow"
    refFilesAll = refFilesAll0(okCR);
    col = col0(okCR); row = row0(okCR);
    uCol = sort(unique(col(:))); uRow = sort(unique(row(:)));
    nCols = numel(uCol); nRows = numel(uRow);
    cBin = zeros(numel(refFilesAll),1); rBin = zeros(numel(refFilesAll),1);
    for i=1:numel(refFilesAll)
        cBin(i) = find(uCol==col(i),1);  rBin(i) = find(uRow==row(i),1);
    end
    if row_increases_up, rB = rBin; else, rB = nRows - rBin + 1; end
    tileAtRC = nan(nRows, nCols);
    for i=1:numel(refFilesAll)
        if ~isnan(tileAtRC(rB(i), cBin(i)))
            error("Duplicate mapping at (rB,c)=(%d,%d).", rB(i), cBin(i));
        end
        tileAtRC(rB(i), cBin(i)) = i;
    end
    if any(isnan(tileAtRC(:))), error("Missing tiles: some (row,col) absent."); end

    nTiles = nRows*nCols;
    tileOrder = zeros(nTiles,1); k = 0;
    for rBin_here = 1:nRows
        if row_increases_up, rB_here = rBin_here; else, rB_here = nRows - rBin_here + 1; end
        if mod(rBin_here,2)==1, colsSnake = 1:nCols; else, colsSnake = nCols:-1:1; end
        for c_here = colsSnake
            k = k + 1; tileOrder(k) = tileAtRC(rB_here, c_here);
        end
    end
    tileFileRef = refFilesAll(tileOrder);

    tileFilesOther = cell(1, numel(otherAvgzDirs));
    for m = 1:numel(otherAvgzDirs)
        tf_m = strings(nTiles,1);
        for t=1:nTiles
            [~, bn, ext] = fileparts(tileFileRef(t));
            cand = fullfile(otherAvgzDirs{m}, bn + ext);
            if isfile(cand), tf_m(t) = cand;
            else
                [c0,r0,tf] = parse_colrow_from_name(tileFileRef(t));
                if ~tf, error("Cannot parse col/row from %s", tileFileRef(t)); end
                tf_m(t) = find_by_colrow(otherAvgzDirs{m}, c0, r0);
            end
        end
        if any(tf_m==""), error("Some ch%d tiles missing.", otherChIDs(m)); end
        tileFilesOther{m} = tf_m;
    end
    inv = zeros(numel(refFilesAll),1); inv(tileOrder) = 1:nTiles;
    tileAtRC = arrayfun(@(v) inv(v), tileAtRC);
else
    refFilesAll = refFilesAll0(okTI);
    tileIdx = tileIdx0(okTI);
    [tileIdxSorted, ord] = sort(tileIdx(:), 'ascend');
    refFilesAll = refFilesAll(ord);
    if numel(unique(tileIdxSorted)) ~= numel(tileIdxSorted)
        error("Duplicate tile indices in REF channel.");
    end
    nTiles = numel(refFilesAll);
    [nRows, nCols] = infer_grid_from_ntiles(nTiles, nRows_tile, nCols_tile);
    if nRows*nCols ~= nTiles
        error("nRows*nCols != nTiles (%d*%d != %d).", nRows, nCols, nTiles);
    end
    tileAtRC = nan(nRows,nCols);
    for t = 1:nTiles
        rB_here  = ceil(t / nCols);
        posInRow = mod(t-1, nCols) + 1;
        if mod(rB_here,2)==1, c_here = posInRow; else, c_here = nCols - posInRow + 1; end
        tileAtRC(rB_here, c_here) = t;
    end
    tileFileRef = refFilesAll;
    tileFilesOther = cell(1, numel(otherAvgzDirs));
    for m = 1:numel(otherAvgzDirs)
        otherFilesAll0 = list_tifs(otherAvgzDirs{m});
        K0 = numel(otherFilesAll0); okO = false(K0,1); tileIdxO0 = nan(K0,1);
        for i=1:K0
            [k,tf] = parse_tileidx_from_name(otherFilesAll0(i));
            if tf, okO(i) = true; tileIdxO0(i) = k; end
        end
        otherFilesAll = otherFilesAll0(okO); tileIdxO = tileIdxO0(okO);
        [tileIdxOSorted, ordO] = sort(tileIdxO(:), 'ascend');
        otherFilesAll = otherFilesAll(ordO);
        if ~isequal(tileIdxOSorted(:), tileIdxSorted(:))
            error("Tile index mismatch between ch%d and ch(ref).", otherChIDs(m));
        end
        tileFilesOther{m} = otherFilesAll;
    end
end

c_of_t = zeros(numel(tileFileRef),1); rB_of_t = zeros(numel(tileFileRef),1);
for rr = 1:nRows
    for cc = 1:nCols
        t = tileAtRC(rr,cc);  c_of_t(t) = cc;  rB_of_t(t) = rr;
    end
end
end

function ch = discover_channels(rawDir)
% Channel IDs that actually have a ch<N>\avgz folder under rawDir, ascending.
% Matches the deinterleaver's output layout; the number of channels is not
% knowable in advance, so it is read off the filesystem rather than assumed.
ch = [];
d = dir(char(rawDir));
d = d([d.isdir]);
for i = 1:numel(d)
    tok = regexp(d(i).name, '^ch(\d+)$', 'tokens', 'once');
    if isempty(tok), continue; end
    if isfolder(fullfile(rawDir, d(i).name, "avgz"))
        ch(end+1) = str2double(tok{1}); %#ok<AGROW>
    end
end
ch = sort(unique(ch));
end

function A = read_acq_params(rawDir, refTile)
A = struct('ok', false, 'src', "", 'zs', [], 'zoom', NaN, 'fovUm', NaN, ...
           'umPerPx', NaN, 'channelSave', [], 'nSlices', NaN, 'framesPerSlice', NaN);
c = dir(fullfile(rawDir, "*.tif"));
if isempty(c), return; end
files = string(fullfile({c.folder}, {c.name}));
names = string({c.name});

pick = 1;
[cr, rr, tf] = parse_colrow_from_name(refTile);
if tf
    hit = find(contains(names, sprintf("col%02d_row%02d", cr, rr)), 1);
    if ~isempty(hit), pick = hit; end
end
try
    t = Tiff(files(pick), 'r'); s = string(t.getTag('Software')); t.close();
catch
    return
end
A.src            = names(pick);
A.zs             = si_vec(s, "SI.hStackManager.zs");
A.zoom           = si_num(s, "SI.hRoiManager.scanZoomFactor");
A.channelSave    = si_vec(s, "SI.hChannels.channelSave");
A.nSlices        = si_num(s, "SI.hStackManager.numSlices");
A.framesPerSlice = si_num(s, "SI.hStackManager.framesPerSlice");
px  = si_num(s, "SI.hRoiManager.pixelsPerLine");
fov = si_mat(s, "SI.hRoiManager.imagingFovUm");
if ~isempty(fov) && size(fov,2) >= 2 && isfinite(px) && px > 0
    A.fovUm   = max(fov(:,1)) - min(fov(:,1));
    A.umPerPx = A.fovUm / px;
end
A.ok = true;
end

function v = si_num(s, key)
v = si_vec(s, key);
if isempty(v), v = NaN; else, v = v(1); end
end

function v = si_vec(s, key)
v = [];
tok = regexp(s, key + "\s*=\s*([^\n\r]*)", 'tokens', 'once');
if isempty(tok), return; end
raw = erase(string(tok{1}), ["[", "]", "{", "}", "'"]);
raw = replace(replace(raw, ";", " "), ",", " ");
d = str2double(split(strtrim(raw)));
v = d(~isnan(d));
end

function M = si_mat(s, key)
M = [];
tok = regexp(s, key + "\s*=\s*([^\n\r]*)", 'tokens', 'once');
if isempty(tok), return; end
rows = split(erase(string(tok{1}), ["[", "]"]), ";");
for i = 1:numel(rows)
    d = str2double(split(strtrim(replace(rows(i), ",", " ")))); d = d(~isnan(d)).';
    if i == 1, M = zeros(numel(rows), numel(d)); end
    M(i, 1:numel(d)) = d;
end
end

function [zs, src] = find_stack_zs(rawDir, refTile)
zs = []; src = "";
c = dir(fullfile(rawDir, "*.tif"));
if isempty(c), return; end
files = string(fullfile({c.folder}, {c.name}));

pick = 1;
[cr, rr, tf] = parse_colrow_from_name(refTile);
if tf
    hit = find(contains(files, sprintf("col%02d_row%02d", cr, rr)), 1);
    if ~isempty(hit), pick = hit; end
end
zs = read_zs(files(pick));
[~, nm, ex] = fileparts(files(pick));  src = nm + ex;

if numel(files) > 1
    other = read_zs(files(mod(pick, numel(files)) + 1));
    if ~isempty(other) && ~isequal(numel(other), numel(zs))
        warning(['two raw tiles report different slice counts (%d vs %d) -- the ' ...
                 'stack was not uniform across tiles, so a micron range is only ' ...
                 'right for the tile it was read from.'], numel(zs), numel(other));
    end
end
end

function zs = read_zs(f)
zs = [];
try
    t = Tiff(f, 'r'); s = string(t.getTag('Software')); t.close();
catch
    return
end
tok = regexp(s, "SI.hStackManager.zs\s*=\s*([^\n\r]*)", 'tokens', 'once');
if isempty(tok), return; end
raw = erase(string(tok{1}), ["[", "]"]);
raw = replace(replace(raw, ";", " "), ",", " ");
v = str2double(split(strtrim(raw)));
zs = v(~isnan(v));
end

function [sx, sy, have, src] = parse_stage_tokens(files, rawDir, want)
n = numel(files); sx = nan(n,1); sy = nan(n,1); src = "filename";

if want ~= "filename"
    hx = nan(n,1); hy = nan(n,1);
    raws = dir(fullfile(rawDir, "*.tif"));
    if ~isempty(raws)
        rawNames = string({raws.name});
        rawPaths = string(fullfile({raws.folder}, {raws.name}));
        for k = 1:n
            [~, bn] = fileparts(files(k));
            stem = regexprep(bn, '_ch\d+_avgz$', '');
            hit  = find(startsWith(rawNames, stem + "."), 1);
            if isempty(hit)
                [c, r, tf] = parse_colrow_from_name(files(k));
                if tf
                    hit = find(contains(rawNames, sprintf("col%02d_row%02d", c, r)), 1);
                end
            end
            if isempty(hit), continue; end
            mp = read_motor_position(rawPaths(hit));
            if numel(mp) >= 2, hx(k) = mp(1); hy(k) = mp(2); end
        end
    end
    if all(isfinite(hx)) && all(isfinite(hy))
        sx = hx; sy = hy; src = "header (SI.hMotors.motorPosition)";
        have = true; return
    end
    if want == "header"
        error(['stageSource = "header" but motorPosition could not be read for ' ...
               'every tile from raw files in %s (%d of %d found).'], ...
               rawDir, nnz(isfinite(hx)), n);
    end
    if ~isempty(raws)
        fprintf("[stage] header incomplete (%d/%d) -- falling back to filenames\n", ...
                nnz(isfinite(hx)), n);
    end
end

for k = 1:n
    [~, bn] = fileparts(files(k));
    a = regexp(bn, '_x(-?\d+(?:\.\d+)?)', 'tokens', 'once');
    b = regexp(bn, '_y(-?\d+(?:\.\d+)?)', 'tokens', 'once');
    if ~isempty(a), sx(k) = str2double(a{1}); end
    if ~isempty(b), sy(k) = str2double(b{1}); end
end
have = all(isfinite(sx)) && all(isfinite(sy));
end

function mp = read_motor_position(f)
mp = [];
try
    t = Tiff(f, 'r'); s = string(t.getTag('Software')); t.close();
catch
    return
end
tok = regexp(s, "SI.hMotors.motorPosition\s*=\s*([^\n\r]*)", 'tokens', 'once');
if isempty(tok), return; end
raw = erase(string(tok{1}), ["[", "]"]);
raw = replace(replace(raw, ";", " "), ",", " ");
v = str2double(split(strtrim(raw)));
mp = v(~isnan(v));
end

function files = list_tifs(folder)
c = [dir(fullfile(folder,"*.tif")); dir(fullfile(folder,"*.tiff"))];
files = strings(numel(c),1);
for k=1:numel(c), files(k) = string(fullfile(c(k).folder, c(k).name)); end
end

function [col,row,tf] = parse_colrow_from_name(fullpath)
[~, bn, ~] = fileparts(fullpath);
tokC = regexp(bn, 'col(\d+)', 'tokens', 'once');
tokR = regexp(bn, 'row(\d+)', 'tokens', 'once');
tf = ~isempty(tokC) && ~isempty(tokR);
if ~tf, col = NaN; row = NaN;
else,   col = str2double(tokC{1}); row = str2double(tokR{1}); end
end

function [idx, tf] = parse_tileidx_from_name(fullpath)
[~, bn, ~] = fileparts(fullpath);
tok = regexp(bn, '(?i)tile[_-]?(\d+)', 'tokens', 'once');
tf = ~isempty(tok);
if ~tf, idx = NaN; else, idx = str2double(tok{1}); end
end

function [nRows, nCols] = infer_grid_from_ntiles(nTiles, nRows_user, nCols_user)
if ~isempty(nRows_user) && ~isempty(nCols_user)
    nRows = nRows_user; nCols = nCols_user; return
elseif ~isempty(nRows_user)
    nRows = nRows_user; nCols = round(nTiles/nRows); return
elseif ~isempty(nCols_user)
    nCols = nCols_user; nRows = round(nTiles/nCols); return
end
best = [1 nTiles]; bestScore = inf;
for r = 1:floor(sqrt(nTiles))
    if mod(nTiles, r)==0
        c = nTiles/r; r1 = min(r,c); c1 = max(r,c);
        if abs(c1-r1) < bestScore, bestScore = abs(c1-r1); best = [r1 c1]; end
    end
end
nRows = best(1); nCols = best(2);
end

function fn = find_by_colrow(folder, col, row)
c = dir(fullfile(folder, sprintf("*col%02d*row%02d*.tif*", col, row)));
if isempty(c), c = dir(fullfile(folder, sprintf("*col%d*row%d*.tif*", col, row))); end
if isempty(c), error("No file for col%d row%d in %s", col, row, folder); end
fn = string(fullfile(c(1).folder, c(1).name));
end

function [H, W, cls] = read_first_page_size_class(fn)
t = Tiff(fn,'r'); A = t.read(); t.close();
H = size(A,1); W = size(A,2); cls = string(class(A));
end

function n = count_tiff_pages_file(fn)
t = Tiff(fn,'r'); n = 1;
while ~t.lastDirectory(), t.nextDirectory(); n = n + 1; end
t.close();
end

function w = linear_blend_weights(H, W)
[xg, yg] = meshgrid(1:W, 1:H);
w = double(min(min(xg-1, W-xg), min(yg-1, H-yg)) + 1);
end

function B = shift_bilinear(A, fx, fy)
[H, W] = size(A);
[X, Y] = meshgrid(1:W, 1:H);
B = interp2(X, Y, A, X-fx, Y-fy, 'linear', 0);
end

function write_stitched_from_tilefiles(outFile, tileFiles, place, outH, outW, H0, W0, zUse, inClass, rowsPerStrip, desc)
nT = numel(tileFiles);
tw = Tiff(outFile, 'w8');
cleanup = onCleanup(@() tw.close());
tag.ImageLength = outH;  tag.ImageWidth = outW;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.BitsPerSample = 16;  tag.SamplesPerPixel = 1;
tag.RowsPerStrip = min(outH, rowsPerStrip);
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.Software = 'NChan_vol_stitch_zUse_260910';
tag.ImageDescription = char(desc);
if strcmpi(inClass,'int16'), tag.SampleFormat = Tiff.SampleFormat.Int;
else,                        tag.SampleFormat = Tiff.SampleFormat.UInt; end

rd = cell(nT,1);
for t = 1:nT, rd{t} = Tiff(tileFiles(t), 'r'); end
cl2 = onCleanup(@() cellfun(@(o) o.close(), rd));

for zi = 1:numel(zUse)
    num = zeros(outH, outW, 'double');
    den = zeros(outH, outW, 'double');
    for t = 1:nT
        rd{t}.setDirectory(zUse(zi));
        P = double(rd{t}.read());
        if size(P,1) ~= H0 || size(P,2) ~= W0
            P = imresize(P, [H0 W0], 'bicubic');
        end
        if abs(place(t).fx)>1e-6 || abs(place(t).fy)>1e-6
            P = shift_bilinear(P, place(t).fx, place(t).fy);
        end
        rr = place(t).rr;  cc = place(t).cc;
        num(rr,cc) = num(rr,cc) + P .* place(t).w2;
        den(rr,cc) = den(rr,cc) + place(t).w2;
    end
    den(den==0) = 1;
    M = num ./ den;
    tw.setTag(tag);
    if strcmpi(inClass,'int16'), tw.write(int16(round(max(min(M,32767),-32768))));
    else,                        tw.write(uint16(round(max(min(M,65535),0)))); end
    if zi < numel(zUse), tw.writeDirectory(); end
end
end
