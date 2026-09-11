function NChan_vol_stitch_zUse_260910()
%NCHAN_VOL_STITCH_ZUSE_260910  Auto tile stitcher -- bounded NCC, self-calibrating.
%
% Drop-in replacement for NChan_vol_stitch_zUse.m.  The file discovery, the two
% naming conventions, the channel pairing, zUse, the blending and the TIFF writer
% are all as they were.  What changed is the registration, which is why the old
% one failed.
%
% ============================ WHY THE OLD ONE FAILED ==========================
% It cut the correlation patch AT THE INITIAL GUESS and phase-correlated that
% patch.  A patch of width P can only represent shifts of +/- P/2 before the FFT
% wraps, so the guess had to be right to within half the overlap it predicted --
% the one thing it cannot be, since the guess is what is being corrected.
%
% Measured on 260909_ChAT_g8m_Shiverer\map:
%   um_per_px_guess 0.8889195  ->  step 450.0 px,  patch 512-450 = 62 px wide
%   unambiguous range           ->  +/- 31 px
%   true step (144 curated pairs) 418.75 px  ->  correction needed -31.2 px
%   i.e. 101% of the range.  Every horizontal edge aliases or fails psrMin, then
%   falls back to the wrong guess at weight 0.2, and the global solve is dragged
%   onto a 450 px lattice that is 31 px per step too coarse -- 250 px of error
%   accumulated across 9 columns.
%
% The scale guess was simply wrong: 0.8889 is the ScanImage header value, and on
% this rig the stage and the header disagree by ~6% (0.9462 measured).  But
% hard-coding 0.9462 would only move the cliff, so:
%
% ============================== WHAT IS DIFFERENT =============================
% 1. THE TEMPLATE IS CUT INDEPENDENTLY OF THE GUESS.  A strip is taken from the
%    LEADING EDGE of the second tile, sized as a fraction of the EXPECTED overlap,
%    and searched for inside a generous strip of the first tile.  The guess then
%    only has to land within SearchRadiusPx, not within half a patch.
% 2. normxcorr2, NOT raw phase correlation.  Linear (zero-padded) rather than
%    circular, so there is no wraparound to alias into; normalised, so a
%    brightness difference between tiles cannot masquerade as a shift; and the
%    peak value is a correlation coefficient, which is an interpretable accept
%    threshold rather than a unitless PSR.
% 3. THE SEARCH IS EXPLICITLY BOUNDED and the bound is reported.  A residual that
%    wants to leave the window is a FAILURE, not a silently accepted wrap.
% 4. INITIAL PLACEMENT FROM THE STAGE when the filenames carry _x###_y### tokens,
%    instead of a uniform step*index lattice.  Handles a non-uniform or skipped
%    step, and is what the manual GUI already does.
% 5. SELF-CALIBRATION.  After the pass, um/px is re-measured from the successful
%    edges against the stage travel, and reported.  If it disagrees with the
%    guess by more than ScaleWarnFrac you are told, with the number to paste back.
% 6. FAILED EDGES ARE DROPPED, not kept at low weight.  A wrong edge held at
%    weight 0.2 still pulls; 40 of them pull hard and in the same direction.
%    Connectivity is checked instead, and a disconnected grid is an error.
% 7. IT REPORTS. Per-edge outcome, residual distribution, post-solve edge
%    residuals, and a QC figure -- so the next failure is diagnosable.
%
% Runqi Zhang / 2026-09-10

clc;

%% ========================= USER SETTINGS =========================
rawDir = "E:\260804_shiver_dbh\map";   % contains ch*/avgz

refChID   = 3; otherChID = 1;
zUse_um     = [-300 0];   % in micron, e.g. [-300 0] -> use 0-300 micron under surface. [] = whole stack. Takes precedence over zUseProj.
zUseProj    = [];         % raw PAGE INDICES, only used when zUse_um is [] -- the fallback for a folder with no raw tile.
rescaled_um = 1;          % output pixel size in um. [] = keep the native pixels.
                          % Resampled with the MEASURED per-axis scale, so one
                          % output pixel is one micron and the mosaic is isotropic
                          % -- which also removes the ~6% x-vs-y difference the
                          % raw pixels carry (1.618 vs 1.524 um/px on 260804).
stageSource = "auto";     % "auto" | "header" | "filename"
                          % "header" = SI.hMotors.motorPosition from the raw tile
                          % (authoritative).  "auto" uses it when every tile has
                          % a raw file and falls back to the filename tokens.

% 260804: zoom 1.2, FOV 758.55 um, 1.4815 um/px, 600 um step -> 107 px overlap.
% 260909: zoom 2,   FOV 455.13 um, 0.8889 um/px, 400 um step ->  93 px overlap.
step_um         = 400;        % nominal stage step; only a fallback now
um_per_px_guess = 1.4815;     % the ScanImage HEADER value for zoom 1.2.
% This is a SEED, not a calibration, and it does not have to be right -- the
% bounded search below tolerates any error up to SearchRadiusPx, and the script
% re-measures the true value from the accepted edges and prints it.
%
% Seeding it with a hand-measured number would make the auto-stitcher depend on
% a manual stitch having already been done, which defeats the point. For the
% record the curated 260909 map gives ~0.9462 um/px (400 um of stage travel over
% 418.75 px, 144 pairs), i.e. the header understates the field by ~6% on this
% rig -- but that number came FROM a curation and is only used here to check the
% self-calibration lands in the right place, never to seed it.

placeInit    = "stage";       % "stage" (filename _x###_y### tokens) | "grid"
% WHAT THE FILENAME TOKENS ARE.  `_x-2700_y-1350` is the ACTUAL motor position,
% not the commanded one: auto_acq_n_shutterOFF.m moves, then reads
% hMotors.motorPosition and prints it with %.0f.  So the only difference from
% SI.hMotors.motorPosition in the raw header is rounding to whole microns --
% measured at 0.49 um worst case over 36 tiles of 260804 (sd 0.30).
%
% That makes the filenames a sound source, which matters because an avgz tile is
% a processed TIFF and carries no SI metadata at all; the header is only
% available where a raw tile was kept.  Either way it is just the SEED for a
% search that then measures the step from the images.
row_increases_up = true;      % grid mode only: filename row index goes UP

flipX = false;
flipY = false;

% ---- registration ----
SearchRadiusPx   = 120;       % bound on the residual ALONG the step direction
SearchPerpPx     = 40;        % bound ACROSS it
HighPassSigmaPx  = 0;        % 0 = off.  See the note at the projection step.
TemplateFrac     = 0.60;      % template width as a fraction of expected overlap
MinTemplatePx    = 24;
NccMin           = 0.45;      % accept: peak correlation coefficient
ExclRadiusPx     = 15;        % uniqueness test: ignore this much around the peak
MarginMin        = 0.02;      % accept: peak minus best rival outside that radius
% NO peak-to-sidelobe gate.  PSR is a statistic for PHASE correlation, whose
% surface is a sharp spike on flat noise.  A normxcorr2 surface is smooth and
% broad, so even a perfect match sits only a few SD above its own neighbourhood:
% carrying the old psrMin = 4 across rejected 142 of 144 edges that had ALREADY
% passed NCC >= 0.3.  The peak VALUE is the quality measure for NCC; how far it
% beats the best rival peak is the uniqueness measure.
subpixel         = true;
ScaleWarnFrac    = 0.02;      % warn if the measured um/px is this far off

% ---- global solve ----
PriorWeight      = 0.05;      % weight of the stage prior per tile, vs edge ~1-20
DarkContrast     = 0.025;     % contrast index below which a tile is "nothing"
DarkPriorBoost   = 200;       % prior weight multiplier for those -- effectively pinned
HuberK           = 2.5;       % IRLS: residuals past this many robust SDs are cut
IrlsIters        = 6;

% ---- blending / output ----
useLinearBlending = true;
rowsPerStrip      = 64;
nRows_tile = [];  nCols_tile = [];    % tile_### mode grid override

writeStitched = true;         % false = coords + QC only (fast, for tuning)
compareTo     = "";           % "stitch_manual_coords.mat", or "" to skip
%% =================================================================

refAvgzDir   = fullfile(rawDir, sprintf("ch%d", refChID),   "avgz");
otherAvgzDir = fullfile(rawDir, sprintf("ch%d", otherChID), "avgz");
outDir = fullfile(rawDir, "matlab_stitch");
if ~isfolder(outDir), mkdir(outDir); end



coordMat = fullfile(outDir, sprintf("stitch_coords_from_ch%d.mat", refChID));

if ~isfolder(refAvgzDir),   error("Missing folder: %s", refAvgzDir); end
if ~isfolder(otherAvgzDir), error("Missing folder: %s", otherAvgzDir); end

%% ========================= DISCOVER FILES (unchanged) =========================
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

[tileFileRef, tileFileOther, tileAtRC, nRows, nCols, c_of_t, rB_of_t] = ...
    build_order(mode, refFilesAll0, okCR, okTI, col0, row0, tileIdx0, ...
                otherAvgzDir, otherChID, row_increases_up, nRows_tile, nCols_tile);
nTiles = numel(tileFileRef);
fprintf("[grid] nRows=%d nCols=%d nTiles=%d\n", nRows, nCols, nTiles);

%% ========================= READ SIZE / PAGES =========================
[H0, W0, inClass] = read_first_page_size_class(tileFileRef(1));
nZ0 = count_tiff_pages_file(tileFileRef(1));
fprintf("[tile] size=%dx%d, nZ=%d, class=%s\n", H0, W0, nZ0, inClass);

% ---- acquisition parameters, from ONE raw tile ------------------------------
% Zoom, field of view, pixel size, saved channels and the z axis are properties
% of the ACQUISITION, not of a tile, so a single raw file is enough.  That is the
% normal case when a map is copied off the rig for stitching: all the avgz tiles
% plus one raw tile kept for its metadata.  Per-tile STAGE POSITION is separate
% and needs every raw file -- keeping the two apart is what stops a folder with
% one raw tile from falling back to a typed-in um/px it does not need to.
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

% ---- zUse: MICRONS of motor z, resolved to page indices ---------------------
% zUse is a depth RANGE in microns, not page numbers, so the same setting means
% the same tissue on a stack with a different step size or start.  It is matched
% against SI.hStackManager.zs, the per-slice motor z -- the same convention
% stitch_manual_gui_fast_260818 uses.  Motor z runs NEGATIVE into the brain and
% the ventral surface is where z was zeroed, so [-300 0] is the top 300 um.
% Order does not matter: [0 -300] is the same slab, and both ends are inclusive.
%
% THE AVGZ TILES CARRY NO METADATA AT ALL -- no Software tag, no ImageDescription
% (checked).  So zs has to come from a RAW tile beside ch*/, and if there is none
% a micron range cannot be honoured.  That is an error rather than a silent
% fallback to all pages: quietly projecting 540 um when 300 was asked for changes
% the picture completely and nothing in the output would say so.
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
        zs = ACQ.zs;  zsrc = ACQ.src;          % already read, do not re-open
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

%% ========================= SUM PROJECTIONS =========================
sumProj = cell(nTiles,1);
for t = 1:nTiles
    tt = Tiff(tileFileRef(t),'r');
    S = zeros(H0,W0,'double');
    for iz = zUse, tt.setDirectory(iz); S = S + double(tt.read()); end
    tt.close();
    % HIGH-PASS BEFORE REGISTRATION.  Each avgz tile here is itself 56 z planes,
    % so summing all of them gives a thick projection dominated by SLOW shading --
    % vignetting, depth falloff, a bright surface vessel smeared through the
    % stack.  NCC is computed over the whole template, so that low-frequency
    % component supplies most of the correlation and swamps the cell-sized detail
    % that actually localises the match: the peak goes broad and wanders.
    % Subtracting a Gaussian blur leaves only structure at the scale of somata and
    % vessels, which is what carries the registration information.
    if HighPassSigmaPx > 0
        S = S - imgaussfilt(S, HighPassSigmaPx);
    end
    sumProj{t} = S;
end
if HighPassSigmaPx > 0
    fprintf("[proj] high-pass sigma %g px before registration\n", HighPassSigmaPx);
end

%% ========================= INITIAL PLACEMENT =========================
[stageX, stageY, haveStage, stageSrc] = parse_stage_tokens(tileFileRef, rawDir, stageSource);
step_px = step_um / um_per_px_guess;

if placeInit == "stage" && haveStage
    % Measured rig convention: image column runs along +stage x, image row along
    % -stage y.  Using the per-tile stage position rather than index*step means a
    % non-uniform or skipped step is carried exactly instead of averaged away.
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

%% ========================= EDGES =========================
idxGrid = zeros(nRows,nCols);
for t = 1:nTiles, idxGrid(rB_of_t(t), c_of_t(t)) = t; end

E = struct('i',{},'j',{},'dxG',{},'dyG',{},'axis',{});
for rB = 1:nRows
    for c = 1:nCols
        i = idxGrid(rB,c);
        if c < nCols
            j = idxGrid(rB,c+1);
            E(end+1) = struct('i',i,'j',j,'dxG',x0(j)-x0(i),'dyG',y0(j)-y0(i),'axis','x'); %#ok<AGROW>
        end
        if rB < nRows
            j = idxGrid(rB+1,c);
            E(end+1) = struct('i',i,'j',j,'dxG',x0(j)-x0(i),'dyG',y0(j)-y0(i),'axis','y'); %#ok<AGROW>
        end
    end
end
nE = numel(E);
fprintf("[reg] %d edges, search +/-%d px along / +/-%d across\n", nE, SearchRadiusPx, SearchPerpPx);

%% ========================= PAIRWISE NCC =========================
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
% The NCC distribution over ALL edges, not only the accepted ones.  Without it a
% threshold is picked blind, and "the gate is too strict" and "the data is bad"
% look identical from the accept count alone.
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

%% ========================= SELF-CALIBRATION =========================
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

%% ========================= GLOBAL SOLVE =========================
% ---- REBUILD THE PRIOR AT THE MEASURED SCALE -------------------------------
% The prior started life at um_per_px_guess, which is the header value and is
% wrong here by 8%.  A prior that is stretched 8% is not a mild bias: it pulls
% against every correct edge, and pinning a dark tile to it puts that tile
% confidently in the wrong place.  Now that the edges have measured the scale,
% the prior is rebuilt on it -- PER AXIS, because x and y disagree by ~6% (the
% resonant and galvo axes are not calibrated together).
%
% This does not make the prior circular.  It is still the STAGE that says where
% each tile is; the edges only supply the microns-to-pixels conversion, which is
% one number per axis over the whole mosaic, not a per-tile adjustment.
if haveStage && exist('uX','var') && ~isempty(uX) && ~isempty(uY)
    sxPx = median(uX);  syPx = median(uY);
    x0 = (stageX - min(stageX)) / sxPx;
    y0 = (max(stageY) - stageY) / syPx;
    if flipX, x0 = max(x0) - x0; end
    if flipY, y0 = max(y0) - y0; end
    fprintf('[prior] rebuilt at the measured scale: x %.4f, y %.4f um/px\n', sxPx, syPx);
end

% ---- DARK TILES ------------------------------------------------------------
% A tile with no structure has no opinion about where it belongs, and it should
% not be given one.  Contrast index = (p99.5 - p50)/p50 on the reference
% projection: how far the bright tail rises above the tile's own background,
% which is scale-free, so it compares tiles of different exposure.  Anything
% below DarkContrast is pinned hard to its stage position -- typically the map's
% corners, where the field ran off the tissue.  Pinning is better than leaving
% them free: a dark tile usually still produces ONE spurious edge somewhere, and
% a single unopposed bad edge will drag an otherwise unconstrained tile far.
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

% LOOP CLOSURE, before solving anything.  Every 2x2 block of tiles is a cycle:
% going A->B->D must land where A->C->D lands.  The deviation needs no ground
% truth and no solve, and it is the cleanest evidence that a set of edges is
% mutually inconsistent -- which is the failure mode a per-edge quality score
% cannot see, because each edge can look confident on its own and still be wrong.
loop_closure_report(E, dxE, dyE, okE, idxGrid, nRows, nCols);

ii = [E(okE).i].'; jj = [E(okE).j].';
assert(~isempty(ii), 'every edge was rejected -- nothing to solve');

% CONNECTIVITY IS REPORTED, NOT REQUIRED.  A map this size covers millimetres and
% much of it is background: on the 260909 9x9 the median edge NCC was 0.09, not
% because the gate was strict but because a pair of blank strips carries nothing
% to register.  Demanding one edge-connected graph makes an empty corner fatal.
% The stage prior is what covers those tiles -- image evidence where it exists,
% the recorded stage position where it does not -- so an isolated tile is a note,
% not a failure.
report_connectivity(nTiles, ii, jj, tileFileRef);

% Weight by the correlation coefficient: 1 at the accept threshold, 20 at a
% perfect match.  PSR is gone, so it cannot set the weights either.
w = 1 + 19 * max(0, (ncc(okE) - NccMin) / max(eps, 1 - NccMin));

% ROBUST SOLVE.  81 tiles are 162 unknowns and 144 edges are 288 equations, so
% the system is overdetermined 1.8x -- and that redundancy is the whole reason a
% bad edge is survivable.  Plain least squares spreads one bad edge over the
% whole mosaic; IRLS lets the solution itself decide which edges are lying and
% demote them, re-solving until the weights settle.
%
% The gauge comes from a WEAK PRIOR pulling every tile toward its stage position
% rather than from anchoring tile 1.  Anchoring one tile fixes the origin but
% leaves the mosaic free to drift over long distances -- errors accumulate along
% whatever path the edges happen to form.  A weak prior on all 81 costs nothing
% where the edges agree and quietly holds the far corner where they do not.
[x, y, solveInfo] = solve_robust(nTiles, ii, jj, dxE(okE), dyE(okE), w, ...
                                 x0, y0, priorW, HuberK, IrlsIters);
x = x - min(x);  y = y - min(y);
fprintf("[solve] IRLS %d iters, %d/%d edges downweighted below half\n", ...
        IrlsIters, solveInfo.nDown, numel(ii));

% Post-solve residual: how well the solved positions honour the edges they were
% built from.  This is the number that says whether the stitch is trustworthy --
% the pre-solve residual only says how wrong the guess was.
pr = hypot((x(jj)-x(ii)) - dxE(okE), (y(jj)-y(ii)) - dyE(okE));
fprintf("[solve] post-solve edge residual: median %.2f px, 95th %.2f, max %.2f\n", ...
        median(pr), prctile(pr,95), max(pr));

%% ================= RESAMPLE THE MOSAIC TO rescaled_um ==========================
% The tiles are acquired at ~1.6 um/px and the two axes differ by ~6% (x 1.618,
% y 1.524 here), so a native-pixel mosaic is neither round-numbered nor square.
% Rescaling to rescaled_um with the MEASURED per-axis scale fixes both at once: the
% output is isotropic, and one pixel is one micron, so any distance read off it
% is already in microns.
%
% The measured scale is used, not the header's, because the header is 8% off --
% rescaling by it would stamp that error permanently into the pixels.
% BOTH COORDINATE SPACES ARE KEPT.  The solve happens in RAW pixels -- NCC runs
% on the original tiles and nothing is resampled before the measurement -- so the
% raw-pixel solution is the primary result.  The rescaled one is derived from it.
% Saving only the space that happened to be written left a consumer unable to
% tell which one it had.
%   x_raw / y_raw  with H0, W0      -> mosaic at the true acquired pixel size,
%                                      original pixels, no interpolation
%   x_um  / y_um   with Hn, Wn      -> mosaic at rescaled_um, isotropic, one
%                                      pixel = one micron
Hn = H0; Wn = W0;
if ~exist('sxPx','var'), sxPx = um_per_px_guess; syPx = um_per_px_guess; end
x_raw = x;  y_raw = y;
if ~isempty(rescaled_um) && rescaled_um > 0 && haveStage && exist('sxPx','var')
    fx_s = sxPx / rescaled_um;          % output px per input px, per axis
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

if ~exist('x_um','var')      % no rescale requested: the two spaces coincide
    x_um = x_raw; y_um = y_raw;
end
outW = ceil(max(x) + Wn);  outH = ceil(max(y) + Hn);
writtenSpace = "raw";
if ~isempty(rescaled_um) && (Hn ~= H0 || Wn ~= W0), writtenSpace = "rescaled"; end

% The filename says which pixel space the file is in.  A mosaic at the acquired
% pixel size and one resampled to 1 um/px look identical in a viewer but differ
% in every distance read off them, so they must not be able to share a name --
% and with distinct names both can sit in the folder at once.
fprintf("[stitch] mosaic %dx%d\n", outH, outW);

save(coordMat, 'x','y','tileFileRef','tileFileOther','c_of_t','rB_of_t', ...
     'nRows','nCols','zUse','zUse_um','zUseProj','stageSrc', ...
     'x_raw','y_raw','H0','W0', ...          % raw-pixel space (the solve)
     'x_um','y_um','Hn','Wn','rescaled_um','writtenSpace', ...  % rescaled space
     'outH','outW','step_um','um_per_px_guess','sxPx','syPx', ...
     'row_increases_up','flipX','flipY','mode','E','dxE','dyE','okE','ncc','margin','pr');
fprintf("[save] %s\n", coordMat);

%% ========================= CHECK AGAINST CURATION =========================
cur = fullfile(outDir, compareTo);
if compareTo ~= "" && isfile(cur)
    compare_to_curated(cur, x, y, tileFileRef, W0, H0);
end

qc_figure(outDir, E, dxE, dyE, okE, ncc, pr, isx, x, y, x0, y0, W0, H0, refChID);

%% ========================= WRITE =========================
if ~writeStitched
    fprintf("[skip] writeStitched=false -- coords and QC only\n"); return
end
% BOTH SPACES ARE WRITTEN, for every channel.  They answer different questions:
% "_raw" keeps the acquired pixels untouched -- nothing interpolated, which is
% what anything quantitative needs -- while "_rescale" is isotropic at
% rescaled_um, so a distance can be read straight off the image.  One solve
% produces both; only the rendering differs, so there is no reason to force a
% choice up front and a re-run to change it.
% The tag carries the PIXEL SIZE, not just the word: "raw_1.6um" and
% "rescaled_1um" say what a pixel of that file actually is, which is the one
% thing you need when setting the scale in ImageJ and the one thing a bare
% "raw" / "rescale" does not tell you.  The raw figure is the MEASURED pixel,
% rounded to 0.1 um -- the axes differ by ~6%, so it is a label, not a spec;
% sxPx / syPx in the coords file carry the exact per-axis values.
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

    for cc = [refChID otherChID]
        if cc == refChID, tfiles = tileFileRef; else, tfiles = tileFileOther; end
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

%% ============================ REGISTRATION ==================================
function [dx, dy, ok, peak, margin, why] = refine_edge_ncc(A, B, dxG, dyG, ...
        Ralong, Rperp, tfrac, tmin, nccMin, exclR, marginMin, subpix)
%REFINE_EDGE_NCC  Bounded normalised cross-correlation of one tile pair.
%
% THE TEMPLATE IS CUT FROM THE LEADING EDGE OF B AND ITS SIZE COMES FROM THE
% EXPECTED OVERLAP, NOT FROM THE GUESSED POSITION.  That is the whole fix: the
% old code cut its patch at the guess, so the representable shift range shrank
% exactly as the guess got worse, and a guess wrong by more than half the
% predicted overlap could never be corrected.  Here the guess only selects how
% big a strip to take and where to centre the search window.

dx = dxG; dy = dyG; ok = false; peak = NaN; margin = NaN; why = "";

% Normalise the edge so the step is non-negative along its dominant axis; the
% answer is un-swapped at the end.  Avoids four near-identical index blocks.
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
Cv = C(th_:ah_, tw_:aw_);          % fully-overlapping positions only
if isempty(Cv), why = "no valid overlap"; return; end

[peak, k] = max(Cv(:));
[fr, fc]  = ind2sub(size(Cv), k);

% Uniqueness: how far the peak beats the best RIVAL peak outside an exclusion
% radius.  On a smooth NCC surface the immediate neighbourhood is nearly as high
% as the peak by construction, so the rival has to be sought beyond it.
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

% The bound is enforced, not assumed: a peak that wants to sit outside the
% window is a rejection.  Silently accepting it is how the old code turned a
% wrapped alias into a confident wrong answer.
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
%REPORT_CONNECTIVITY  Component structure of the accepted-edge graph.
%
% Reports the LARGEST component, not connectivity to tile 1.  Tile 1 is the
% acquisition-order start, i.e. a CORNER, and a corner is the likeliest tile in
% the whole map to be empty -- on 260804 it is the dimmest tile in the reference
% channel.  Measuring the graph from it said "35 of 36 disconnected" about a
% mosaic whose other 35 tiles were perfectly connected to each other.
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
            if ~lab(u), lab(u) = comp; stack(end+1) = u; end %#ok<AGROW>
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
% A grid held together by too few accepted edges solves to something, and that
% something is a superposition of rigid pieces floating relative to each other.
% lsqminnorm will not complain; this does.
seen = false(n,1); seen(1) = true;
adj = cell(n,1);
for e = 1:numel(ii)
    adj{ii(e)}(end+1) = jj(e);  adj{jj(e)}(end+1) = ii(e);
end
stack = 1;
while ~isempty(stack)
    v = stack(end); stack(end) = [];
    for u = adj{v}
        if ~seen(u), seen(u) = true; stack(end+1) = u; end %#ok<AGROW>
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
%SOLVE_ROBUST  Weighted least squares over the tile graph, with IRLS and a prior.
%
% Unknowns are the 2n tile positions; each edge contributes two equations
% x_j - x_i = dx and y_j - y_i = dy; each tile contributes two weak equations
% x_t = priorX_t, y_t = priorY_t.  The prior both fixes the gauge and stops
% long-range drift, so no tile has to be anchored.
%
% IRLS: solve, measure each edge's residual against the solution, rescale the
% weights by a Huber factor, repeat.  An edge that disagrees with the consensus
% of its neighbours loses its vote instead of bending the mosaic around itself.
px = px(:) - min(px);  py = py(:) - min(py);
w  = w0(:);
info = struct('nDown', 0, 'sigma', NaN);
for it = 1:nIter
    [x, y] = solve_once(n, ii, jj, dx, dy, w, px, py, priorW);
    rx = (x(jj)-x(ii)) - dx;   ry = (y(jj)-y(ii)) - dy;
    r  = hypot(rx, ry);
    % Robust scale from the residuals themselves: the MAD, not the SD, so a few
    % badly wrong edges cannot inflate the very threshold meant to catch them.
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
%LOOP_CLOSURE_REPORT  Consistency of every 2x2 cycle of tiles.
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
        err(end+1) = hypot(ex, ey); %#ok<AGROW>
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

%% ============================== REPORTING ===================================
function compare_to_curated(curFile, x, y, tileFileRef, W0, H0)
C = load(curFile);
if ~isfield(C,'x') || ~isfield(C,'files'), fprintf("[check] %s has no x/files\n", curFile); return; end
% stitch_manual_gui_fast_260818 folds the nudge into x already.
[~, bnA] = cellfun(@(s) fileparts(char(s)), num2cell(tileFileRef), 'uni', 0);
[~, bnC] = cellfun(@(s) fileparts(char(s)), num2cell(C.files),     'uni', 0);
[tf, loc] = ismember(bnA, bnC);
if nnz(tf) < 3, fprintf("[check] could not match tiles to the curated set\n"); return; end
dx = (x(tf) - C.x(loc(tf)));  dy = (y(tf) - C.y(loc(tf)));
dx = dx - median(dx);  dy = dy - median(dy);      % placement is gauge-free
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
% outDir is a STRING, so [base '.png'] builds a 1x2 STRING ARRAY instead of
% concatenating, and exportgraphics reports "unrecognized inputs". Use +.
base = fullfile(outDir, sprintf('autostitch_qc_ch%d', refCh));
exportgraphics(f, base + ".png", 'Resolution', 170, 'BackgroundColor','white');
close(f);
fprintf("[qc]   %s.png\n", base);
end

%% ========================= ORDER / NAMING (unchanged logic) =================
function [tileFileRef, tileFileOther, tileAtRC, nRows, nCols, c_of_t, rB_of_t] = ...
    build_order(mode, refFilesAll0, okCR, okTI, col0, row0, tileIdx0, ...
                otherAvgzDir, otherChID, row_increases_up, nRows_tile, nCols_tile)
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

    tileFileOther = strings(nTiles,1);
    for t=1:nTiles
        [~, bn, ext] = fileparts(tileFileRef(t));
        cand = fullfile(otherAvgzDir, bn + ext);
        if isfile(cand), tileFileOther(t) = cand;
        else
            [c0,r0,tf] = parse_colrow_from_name(tileFileRef(t));
            if ~tf, error("Cannot parse col/row from %s", tileFileRef(t)); end
            tileFileOther(t) = find_by_colrow(otherAvgzDir, c0, r0);
        end
    end
    if any(tileFileOther==""), error("Some ch%d tiles missing.", otherChID); end
    % remap tileAtRC to acquisition-order indices
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
    otherFilesAll0 = list_tifs(otherAvgzDir);
    K0 = numel(otherFilesAll0); okO = false(K0,1); tileIdxO0 = nan(K0,1);
    for i=1:K0
        [k,tf] = parse_tileidx_from_name(otherFilesAll0(i));
        if tf, okO(i) = true; tileIdxO0(i) = k; end
    end
    otherFilesAll = otherFilesAll0(okO); tileIdxO = tileIdxO0(okO);
    [tileIdxOSorted, ordO] = sort(tileIdxO(:), 'ascend');
    otherFilesAll = otherFilesAll(ordO);
    if ~isequal(tileIdxOSorted(:), tileIdxSorted(:))
        error("Tile index mismatch between channels.");
    end
    tileFileOther = otherFilesAll;
end

c_of_t = zeros(numel(tileFileRef),1); rB_of_t = zeros(numel(tileFileRef),1);
for rr = 1:nRows
    for cc = 1:nCols
        t = tileAtRC(rr,cc);  c_of_t(t) = cc;  rB_of_t(t) = rr;
    end
end
end

function A = read_acq_params(rawDir, refTile)
%READ_ACQ_PARAMS  Acquisition parameters from ONE raw tile.
%
% z positions, zoom, field of view, pixel size, saved channels, slice and frame
% counts.  All are properties of the acquisition, so one raw file carries them
% for the whole map -- which is what makes a folder of avgz tiles plus a single
% raw tile fully self-describing.  Deliberately separate from the per-tile stage
% position, which needs every raw file.
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
    % FOV form, not PIX_BASE/zoom: it already carries scanAngleMultiplier and
    % any scan rotation, so it stays right on a non-square or shifted field.
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
    M(i, 1:numel(d)) = d; %#ok<AGROW>
end
end

function [zs, src] = find_stack_zs(rawDir, refTile)
%FIND_STACK_ZS  Per-slice motor z, from a raw tile beside the ch*/ folders.
%
% Prefers the raw tile matching the reference tile's own col/row, so the z axis
% belongs to a tile that is actually in the mosaic.  Falls back to any raw tile
% in the folder: the stack parameters come from one acquisition script and are
% identical across tiles, which is checked by comparing against a second file
% when one exists.
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
%PARSE_STAGE_TOKENS  Per-tile stage position, from the header if it exists.
%
% SI.hMotors.motorPosition in the RAW tile is the authoritative record of where
% the stage went.  The _x###_y### in a filename is the SAME quantity, just
% rounded: auto_acq_n_shutterOFF.m moves, then reads hMotors.motorPosition and
% prints %.0f.  Measured on 260804 the two agree to within 0.49 um over 36 tiles
% (sd 0.30) -- rounding alone.  The header is preferred only for that last
% fraction of a micron; a stalled or nudged stage IS recorded either way.
n = numel(files); sx = nan(n,1); sy = nan(n,1); src = "filename";

if want ~= "filename"
    hx = nan(n,1); hy = nan(n,1);
    raws = dir(fullfile(rawDir, "*.tif"));
    if ~isempty(raws)
        rawNames = string({raws.name});
        rawPaths = string(fullfile({raws.folder}, {raws.name}));
        for k = 1:n
            [~, bn] = fileparts(files(k));
            % avgz basename is the raw basename plus _ch<N>_avgz
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

%% ============================== IO / BLEND ==================================
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
        % Resample to the output pixel size BEFORE the sub-pixel shift, because
        % the shift is expressed in OUTPUT pixels -- doing it the other way round
        % would apply an output-scale offset to input-scale data.
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
