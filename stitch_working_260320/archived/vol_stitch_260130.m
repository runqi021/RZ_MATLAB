stitch_snake_7x7_imagejstyle_scanimage_meta

%%
function stitch_snake_7x7_imagejstyle_scanimage_meta()
% stitch_snake_7x7_imagejstyle_scanimage_meta
% ------------------------------------------------------------
% ScanImage-aware snake tiling stitcher (with metadata):
% - Finds tile_01 ... tile_NN under rawDir (recursive)
% - Reads ScanImage metadata from ImageDescription/Software tag:
%     * saved channels (e.g., [1 3])
%     * framesPerSlice (T), numSlices (Z) and z-step / zs
%     * imagingFovUm + pixelsPerLine/linesPerFrame -> um_per_px
%     * scanFrameRate / scanFramePeriod -> finterval
% - Uses one reference channel to compute a registration projection per tile
% - Pairwise phase correlation on neighbor overlaps, then global LSQ
% - Stitches and writes an ImageJ hyperstack TIFF (C x Z x T)
%
% Output: rawDir\matlab_stitch\stitched_hyperstack_scanimage.tif

%% ---------------- USER ----------------
rawDir = "E:\live_251104_definitely_not_sst_tdtomato_fitc\raw\groupZ";
outDir = fullfile(rawDir, "matlab_stitch");
outTif = fullfile(outDir, "stitched_hyperstack_scanimage.tif");

% If you want to force grid dims, set these; otherwise inferred from #tiles
nRows = [];     % e.g. 9
nCols = [];     % e.g. 9

% If you want to force nominal step size (um) for initial guesses (optional)
step_um_force = 400;   % set [] to not force; used only if stage positions not parsed

% optional flips (if mosaic mirrored)
flipX = false;
flipY = false;

% Which saved channel to use for registration projection:
%  - 'first' uses the first saved channel in metadata (usually best)
%  - or set to an actual channel ID (e.g. 1 or 3)
regChannel = "first";

% registration controls
minOverlapFrac = 0.05;
psrMin = 5;
maxResidualPx = 150;
subpixel = true;

% speed: how many time frames to sample per z when building reg projections
% (use Inf to use all)
regSampleT = 50;

% blending
useLinearBlending = true;

% TIFF writing
rowsPerStrip = 64;

%% ---------------- FIND TILE FILES ----------------
if ~isfolder(rawDir), error("Not a folder: %s", rawDir); end
if ~isfolder(outDir), mkdir(outDir); end
if isfile(outTif), delete(outTif); end

tileFile = find_all_tile_files(rawDir);
nTiles = numel(tileFile);
if nTiles < 2, error("Found <2 tiles under %s", rawDir); end
fprintf("[stitch] Found %d tiles.\n", nTiles);

% infer grid if not provided
if isempty(nRows) || isempty(nCols)
    [nRows, nCols] = infer_grid_dims(nTiles);
    fprintf("[grid] inferred nRows=%d, nCols=%d\n", nRows, nCols);
else
    if nRows*nCols ~= nTiles
        error("nRows*nCols (%d) != nTiles (%d)", nRows*nCols, nTiles);
    end
end

%% ---------------- READ META FROM FIRST TILE ----------------
meta0 = read_scanimage_metadata(tileFile(1));

% image size / class from first page
tt0 = Tiff(tileFile(1), 'r');
tt0.setDirectory(1);
I0 = tt0.read();
[H0,W0] = size(I0);
inClass = class(I0);
tt0.close();

% um_per_px from ScanImage if possible
um_per_px = NaN;
[um_per_px, um_per_px_x, um_per_px_y] = infer_um_per_px(meta0, H0, W0);
if isfinite(um_per_px)
    fprintf("[meta] um/px ~ %.6f (x=%.6f, y=%.6f)\n", um_per_px, um_per_px_x, um_per_px_y);
else
    warning("[meta] could not infer um/px from metadata; you MUST set it manually in infer_um_per_px().");
end

% saved channels
chSave = infer_saved_channels(meta0);
nCh = numel(chSave);
fprintf("[meta] saved channels: [%s] (nCh=%d)\n", num2str(chSave), nCh);

% Z/T dimensions
[nZ, nT, zStep_um, finterval_s, orderMode] = infer_zt(meta0, tileFile(1), nCh);
fprintf("[meta] Z=%d, T=%d, zStep=%.4f um, finterval=%.6f s, order=%s\n", nZ, nT, zStep_um, finterval_s, orderMode);

% choose reg channel index in saved channel list
if isstring(regChannel) && regChannel=="first"
    regChID = chSave(1);
else
    regChID = double(regChannel);
end
regChIdx = find(chSave==regChID, 1);
if isempty(regChIdx)
    error("Requested regChannel=%g not in saved channels [%s].", regChID, num2str(chSave));
end
fprintf("[reg] using channel ID %d (index %d in saved list)\n", regChID, regChIdx);

%% ---------------- OPEN TILES + BUILD REG PROJECTIONS ----------------
tile = struct([]);
stageXY = nan(nTiles,2);

for t = 1:nTiles
    fn = tileFile(t);
    tt = Tiff(fn, 'r');

    tile(t).t = tt;
    tile(t).fn = fn;

    % verify size and nPages
    tt.setDirectory(1);
    I1 = tt.read();
    if ~isequal(size(I1), [H0 W0])
        error("Tile %d size mismatch.", t);
    end

    nPages = count_tiff_pages(tt);
    tile(t).nPages = nPages;

    % parse per-tile metadata (for stage positions etc.)
    meta = read_scanimage_metadata(fn);
    tile(t).meta = meta;

    stageXY(t,:) = infer_stage_xy(meta); % [x y] in motor units (usually um)

    % build registration projection from (regChIdx), across Z and sampled T
    tile(t).sumProj = build_reg_projection(tt, H0, W0, nCh, nZ, nT, regChIdx, orderMode, regSampleT);

    if t==1
        fprintf("[tile] size=%dx%d, class=%s, pages=%d\n", H0,W0,inClass,nPages);
    end
end
cleanupTiles = onCleanup(@() close_all_tiles(tile)); %#ok<NASGU>

%% ---------------- GRID: snake order, start left-bottom ----------------
idxRC = zeros(nRows, nCols);
idx = 1;
for rB = 1:nRows
    if mod(rB,2)==1
        cols = 1:nCols;        % left->right
    else
        cols = nCols:-1:1;     % right->left
    end
    for c = cols
        idxRC(rB,c) = idx;
        idx = idx + 1;
    end
end

%% ---------------- INITIAL POSITIONS (prefer stageXY if valid) ----------------
useStage = all(isfinite(stageXY(:))) && (range(stageXY(:,1))>0) && (range(stageXY(:,2))>0) && isfinite(um_per_px);

x0 = zeros(nTiles,1);
y0 = zeros(nTiles,1);

if useStage
    % stage to pixels (relative)
    x_um = stageXY(:,1) - min(stageXY(:,1));
    y_um = stageXY(:,2) - min(stageXY(:,2));

    x0 = x_um / um_per_px;
    y0 = y_um / um_per_px;

    fprintf("[init] using stage positions for initial guesses.\n");
else
    % fallback to nominal step (snake grid)
    if isempty(step_um_force)
        error("No valid stage XY from metadata and step_um_force is empty. Provide step_um_force.");
    end
    if ~isfinite(um_per_px)
        error("No valid um_per_px. Need ScanImage FOV meta or set um_per_px manually.");
    end

    step_px = step_um_force / um_per_px;
    for rB = 1:nRows
        rTop = nRows - rB + 1;
        for c = 1:nCols
            tt = idxRC(rB,c);
            x0(tt) = (c-1) * step_px;
            y0(tt) = (rTop-1) * step_px;
        end
    end
    fprintf("[init] using snake grid with step_um_force=%.3f um.\n", step_um_force);
end

if flipX, x0 = max(x0) - x0; end
if flipY, y0 = max(y0) - y0; end

%% ---------------- BUILD NEIGHBOR EDGES (4-neighborhood in grid) ----------------
edges = []; % [i j dxGuess dyGuess]
for rB = 1:nRows
    for c = 1:nCols
        i = idxRC(rB,c);
        if c < nCols
            j = idxRC(rB,c+1);
            edges(end+1,:) = [i j (x0(j)-x0(i)) (y0(j)-y0(i))]; %#ok<AGROW>
        end
        if rB < nRows
            j = idxRC(rB+1,c); % neighbor above
            edges(end+1,:) = [i j (x0(j)-x0(i)) (y0(j)-y0(i))]; %#ok<AGROW>
        end
    end
end
fprintf("[reg] edges=%d\n", size(edges,1));

%% ---------------- PAIRWISE REG ON REG PROJ ----------------
ii = edges(:,1); jj = edges(:,2);
dxG = edges(:,3); dyG = edges(:,4);

dxE = dxG; dyE = dyG;
wE  = 0.2 * ones(size(dxG));    % fallback if registration fails

for e = 1:numel(dxG)
    i = ii(e); j = jj(e);

    [dx, dy, ok, psr] = refine_shift_phasecorr(tile(i).sumProj, tile(j).sumProj, dxG(e), dyG(e), ...
        minOverlapFrac, psrMin, maxResidualPx, subpixel);

    if ok
        dxE(e) = dx;
        dyE(e) = dy;
        wE(e)  = min(20, max(1, psr));
    end
end
fprintf("[reg] ok edges=%d/%d\n", nnz(wE>1), numel(wE));

%% ---------------- GLOBAL OPT (weighted LSQ) ----------------
[x, y] = solve_global_positions_weighted(nTiles, ii, jj, dxE, dyE, wE);
x = x - min(x);
y = y - min(y);

outW = ceil(max(x) + W0);
outH = ceil(max(y) + H0);
fprintf("[stitch] mosaic=%dx%d, C=%d, Z=%d, T=%d\n", outH, outW, nCh, nZ, nT);

%% ---------------- PRECOMPUTE BLENDING & TILE PLACEMENT ----------------
if useLinearBlending
    wTile = linear_blend_weights(H0, W0);
else
    wTile = ones(H0,W0,'double');
end

place = struct([]);
for t = 1:nTiles
    xi = x(t); yi = y(t);
    xInt = floor(xi); yInt = floor(yi);
    fx = xi - xInt; fy = yi - yInt;

    r0 = yInt + 1;
    c0 = xInt + 1;
    rr = r0:(r0+H0-1);
    cc = c0:(c0+W0-1);

    if subpixel && (abs(fx) > 1e-6 || abs(fy) > 1e-6)
        w2 = shift_bilinear_cached(wTile, fx, fy);
    else
        w2 = wTile;
    end

    place(t).xInt = xInt; place(t).yInt = yInt;
    place(t).fx = fx; place(t).fy = fy;
    place(t).rr = rr; place(t).cc = cc;
    place(t).w2 = w2;
end

%% ---------------- OUTPUT TIFF SETUP (ImageJ Hyperstack) ----------------
[bitsPerSample, sampleFormat] = class_to_tiff_format(inClass);
bytesPerPixel = bitsPerSample/8;

nImages = nCh * nZ * nT;
estBytes = double(outH)*double(outW)*double(nImages)*double(bytesPerPixel);
useBigTiff = estBytes > 3.5e9;

if useBigTiff
    tout = Tiff(outTif, 'w8');
else
    tout = Tiff(outTif, 'w');
end
cleanupOut = onCleanup(@() tout.close()); %#ok<NASGU>

tag.ImageLength = outH;
tag.ImageWidth  = outW;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.BitsPerSample = bitsPerSample;
tag.SampleFormat = sampleFormat;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = rowsPerStrip;
tag.Software = 'MATLAB snake stitch (ScanImage meta; phasecorr; linear blend)';

% ImageJ hyperstack description
desc = sprintf([ ...
    'ImageJ=1.54\n' ...
    'images=%d\n' ...
    'channels=%d\n' ...
    'slices=%d\n' ...
    'frames=%d\n' ...
    'hyperstack=true\n' ...
    'mode=grayscale\n' ...
    'loop=false\n'], nImages, nCh, nZ, nT);

if isfinite(zStep_um) && zStep_um~=0
    desc = desc + sprintf("spacing=%.6f\nunit=um\n", zStep_um);
end
if isfinite(finterval_s) && finterval_s>0
    desc = desc + sprintf("finterval=%.9f\n", finterval_s);
end
tag.ImageDescription = char(desc);

%% ---------------- STITCH AND WRITE (ImageJ order: C fastest, then Z, then T) ----------------
imgCounter = 0;

for tIdx = 1:nT
    fprintf("[write] T %d/%d\n", tIdx, nT);
    for zIdx = 1:nZ
        for cIdx = 1:nCh
            imgCounter = imgCounter + 1;

            acc  = zeros(outH, outW, 'double');
            wacc = zeros(outH, outW, 'double');

            for tileIdx = 1:nTiles
                page = scanimage_page_index(cIdx, zIdx, tIdx, nCh, nZ, nT, orderMode);

                tile(tileIdx).t.setDirectory(page);
                I = double(tile(tileIdx).t.read());

                fx = place(tileIdx).fx; fy = place(tileIdx).fy;

                if subpixel && (abs(fx) > 1e-6 || abs(fy) > 1e-6)
                    I2 = shift_bilinear_cached(I, fx, fy);
                    w2 = place(tileIdx).w2;
                else
                    I2 = I;
                    w2 = wTile;
                end

                rr = place(tileIdx).rr;
                cc = place(tileIdx).cc;

                acc(rr,cc)  = acc(rr,cc)  + I2 .* w2;
                wacc(rr,cc) = wacc(rr,cc) + w2;
            end

            out = acc ./ max(wacc, eps);
            out = cast(round(out), inClass);

            tout.setTag(tag);
            tout.write(out);
            if imgCounter < nImages
                tout.writeDirectory();
            end
        end
    end
end

fprintf("\n[DONE] %s\n", outTif);

end

%% ===================== helpers =====================

function files = find_all_tile_files(root)
% find all tile_XX tiffs; sorts by XX
cands = [dir(fullfile(root,"**","tile_*.tif")); dir(fullfile(root,"**","tile_*.tiff"))];
if isempty(cands), error("No tile_*.tif(f) under %s", root); end
names = string({cands.name});
nums = nan(numel(cands),1);
for k=1:numel(cands)
    tok = regexp(names(k), "tile_(\d+)", "tokens", "once");
    if ~isempty(tok), nums(k)=str2double(tok{1}); end
end
[~,ord] = sort(nums);
cands = cands(ord);
files = strings(numel(cands),1);
for k=1:numel(cands)
    files(k) = string(fullfile(cands(k).folder, cands(k).name));
end
end

function [r,c] = infer_grid_dims(nTiles)
% pick factor pair closest to square
f = factor_pairs(nTiles);
d = abs(f(:,1)-f(:,2));
[~,i] = min(d);
r = f(i,1); c = f(i,2);
end

function f = factor_pairs(n)
f = [];
for a=1:floor(sqrt(n))
    if mod(n,a)==0
        b = n/a;
        f(end+1,:) = [a b]; %#ok<AGROW>
    end
end
% prefer r>=c? keep as found; user can swap by setting nRows/nCols
end

function meta = read_scanimage_metadata(tiffFile)
% tries ImageDescription and Software tags; returns struct with map meta.kv(key)=valueStr
meta = struct();
meta.kv = containers.Map('KeyType','char','ValueType','char');

tt = Tiff(tiffFile,'r');
tt.setDirectory(1);

txt = "";
try
    d = tt.getTag('ImageDescription');
    if ~isempty(d), txt = txt + string(d) + newline; end
catch
end
try
    s = tt.getTag('Software');
    if ~isempty(s), txt = txt + string(s) + newline; end
catch
end
tt.close();

% normalize separators: ScanImage sometimes uses "↵"
txt = replace(txt, char(8629), newline); % ↵
txt = replace(txt, "\r\n", newline);
txt = replace(txt, "\r", newline);

lines = splitlines(txt);
for i=1:numel(lines)
    L = strtrim(lines(i));
    if startsWith(L,"SI.")
        % expect "key = value"
        parts = split(L, "=");
        if numel(parts)>=2
            key = strtrim(parts(1));
            val = strtrim(join(parts(2:end), "="));
            meta.kv(char(key)) = char(val);
        end
    end
end
end

function chSave = infer_saved_channels(meta)
chSave = [];
if isKey(meta.kv,'SI.hChannels.channelSave')
    chSave = parse_numeric_vector(meta.kv('SI.hChannels.channelSave'));
elseif isKey(meta.kv,'SI.hChannels.channelsActive')
    chSave = parse_numeric_vector(meta.kv('SI.hChannels.channelsActive'));
end
chSave = chSave(:).';
if isempty(chSave)
    % last resort: assume 1 channel
    chSave = 1;
end
end

function [nZ, nT, zStep_um, finterval_s, orderMode] = infer_zt(meta, tiffFile, nCh)
% infer nZ, nT, z step, and ordering assumption
zStep_um = NaN;
finterval_s = NaN;

% frame interval
if isKey(meta.kv,'SI.hRoiManager.scanFramePeriod')
    finterval_s = str2double(meta.kv('SI.hRoiManager.scanFramePeriod'));
elseif isKey(meta.kv,'SI.hRoiManager.scanFrameRate')
    fr = str2double(meta.kv('SI.hRoiManager.scanFrameRate'));
    if isfinite(fr) && fr>0, finterval_s = 1/fr; end
end

% z info
nZ = NaN;
if isKey(meta.kv,'SI.hStackManager.numSlices')
    nZ = str2double(meta.kv('SI.hStackManager.numSlices'));
end

zs = [];
if isKey(meta.kv,'SI.hStackManager.zs')
    zs = parse_numeric_vector(meta.kv('SI.hStackManager.zs'));
    if ~isempty(zs)
        nZ = numel(zs);
        if numel(zs) >= 2
            zStep_um = median(diff(zs));
        end
    end
end

if isKey(meta.kv,'SI.hStackManager.stackZStepSize') && ~isfinite(zStep_um)
    zStep_um = str2double(meta.kv('SI.hStackManager.stackZStepSize'));
end
if isfinite(zStep_um), zStep_um = abs(zStep_um); end

% frames per slice
nT = NaN;
if isKey(meta.kv,'SI.hStackManager.framesPerSlice')
    nT = str2double(meta.kv('SI.hStackManager.framesPerSlice'));
end

% ordering mode
% If fastZ is enabled, often order is volume-interleaved (T major, Z minor).
orderMode = "sliceMajor"; % default: z blocks
if isKey(meta.kv,'SI.hFastZ.enable')
    en = lower(string(meta.kv('SI.hFastZ.enable')));
    if en=="true" || en=="1"
        orderMode = "volumeInterleaved";
    end
end

% final fallback from page count
tt = Tiff(tiffFile,'r');
nPages = count_tiff_pages(tt);
tt.close();

if ~isfinite(nZ) || nZ<1
    nZ = 1;
end
if ~isfinite(nT) || nT<1
    nFramesTotal = floor(nPages / nCh);
    if orderMode=="sliceMajor"
        nT = floor(nFramesTotal / nZ);
    else
        % volumeInterleaved: frames are volumes; each volume has nZ slices
        nT = floor(nFramesTotal / nZ);
    end
end

% guard
nFramesTotal = floor(nPages / nCh);
if nZ*nT > nFramesTotal
    % if inconsistent, collapse z to 1 and put everything in T
    warning("[meta] Z*T > nFramesTotal from pages. Falling back to Z=1.");
    nZ = 1;
    nT = nFramesTotal;
end
end

function [um_per_px, um_per_px_x, um_per_px_y] = infer_um_per_px(meta, H, W)
um_per_px = NaN; um_per_px_x = NaN; um_per_px_y = NaN;

% pixels
ppl = NaN; lpf = NaN;
if isKey(meta.kv,'SI.hRoiManager.pixelsPerLine')
    ppl = str2double(meta.kv('SI.hRoiManager.pixelsPerLine'));
else
    ppl = W;
end
if isKey(meta.kv,'SI.hRoiManager.linesPerFrame')
    lpf = str2double(meta.kv('SI.hRoiManager.linesPerFrame'));
else
    lpf = H;
end

% FOV corners in um
if isKey(meta.kv,'SI.hRoiManager.imagingFovUm')
    f = parse_numeric_matrix_4x2(meta.kv('SI.hRoiManager.imagingFovUm'));
    if ~isempty(f)
        width_um  = max(f(:,1)) - min(f(:,1));
        height_um = max(f(:,2)) - min(f(:,2));
        um_per_px_x = width_um / ppl;
        um_per_px_y = height_um / lpf;
        um_per_px = mean([um_per_px_x um_per_px_y]);
    end
end
end

function xy = infer_stage_xy(meta)
% returns [x y] from SI.hMotors.motorPosition if present
xy = [NaN NaN];
if isKey(meta.kv,'SI.hMotors.motorPosition')
    v = parse_numeric_vector(meta.kv('SI.hMotors.motorPosition'));
    if numel(v) >= 2
        xy = v(1:2);
    end
end
end

function page = scanimage_page_index(cIdx, zIdx, tIdx, nCh, nZ, nT, orderMode)
% returns 1-based TIFF directory index
% Assumes channels are saved sequentially per "frame".
% Two common patterns:
%  1) sliceMajor: all frames of z1, then all frames of z2, ...
%        frameLinear = (z-1)*nT + t
%  2) volumeInterleaved: for each timepoint(volume) t, slices z=1..nZ
%        frameLinear = (t-1)*nZ + z
if orderMode=="volumeInterleaved"
    frameLinear = (tIdx-1)*nZ + zIdx;
else
    frameLinear = (zIdx-1)*nT + tIdx;
end
page = (frameLinear-1)*nCh + cIdx;
end

function sumProj = build_reg_projection(tt, H, W, nCh, nZ, nT, regChIdx, orderMode, regSampleT)
% Average over time (sampled) within each z, then sum across z.
if isinf(regSampleT) || regSampleT>=nT
    tList = 1:nT;
else
    tList = unique(round(linspace(1,nT,regSampleT)));
end

sumProj = zeros(H,W,'double');
for z=1:nZ
    acc = zeros(H,W,'double');
    for tIdx = tList
        page = scanimage_page_index(regChIdx, z, tIdx, nCh, nZ, nT, orderMode);
        tt.setDirectory(page);
        acc = acc + double(tt.read());
    end
    acc = acc / numel(tList);
    sumProj = sumProj + acc;
end
end

function v = parse_numeric_vector(s)
% extracts all numbers (supports brackets, semicolons)
tok = regexp(char(s), '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
v = str2double(tok);
end

function M = parse_numeric_matrix_4x2(s)
% for strings like: [-56.8 -56.8; 56.8 -56.8; 56.8 56.8; -56.8 56.8]
v = parse_numeric_vector(s);
if numel(v) >= 8
    M = reshape(v(1:8), [2 4]).'; % 4x2
else
    M = [];
end
end

function n = count_tiff_pages(t)
t.setDirectory(1);
n = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    n = n + 1;
end
t.setDirectory(1);
end

function close_all_tiles(tile)
for k = 1:numel(tile)
    try, tile(k).t.close(); catch, end %#ok<CTCH>
end
end

function [dx, dy, ok, psr] = refine_shift_phasecorr(A, B, dxGuess, dyGuess, minOverlapFrac, psrMin, maxResidualPx, subpixel)
H = size(A,1); W = size(A,2);
ovW = W - abs(dxGuess);
ovH = H - abs(dyGuess);
if ovW < minOverlapFrac*W || ovH < minOverlapFrac*H
    dx = dxGuess; dy = dyGuess; ok = false; psr = NaN; return;
end
[Af, Bf] = overlap_patches(A, B, dxGuess, dyGuess);

Af = Af - mean(Af(:));
Bf = Bf - mean(Bf(:));
w = hann2(size(Af,1), size(Af,2));
Af = Af .* w;
Bf = Bf .* w;

[resDx, resDy, psr] = phasecorr_shift(Af, Bf, subpixel);

if ~isfinite(psr) || psr < psrMin || abs(resDx) > maxResidualPx || abs(resDy) > maxResidualPx
    dx = dxGuess; dy = dyGuess; ok = false; return;
end
dx = dxGuess + resDx;
dy = dyGuess + resDy;
ok = true;
end

function [fixedPatch, movingPatch] = overlap_patches(A, B, dx, dy)
H = size(A,1); W = size(A,2);
dx = round(dx); dy = round(dy);

if dx >= 0
    a_c = (dx+1):W;
    b_c = 1:(W-dx);
else
    a_c = 1:(W+dx);
    b_c = (-dx+1):W;
end

if dy >= 0
    a_r = (dy+1):H;
    b_r = 1:(H-dy);
else
    a_r = 1:(H+dy);
    b_r = (-dy+1):H;
end

fixedPatch  = A(a_r, a_c);
movingPatch = B(b_r, b_c);
end

function [dx, dy, psr] = phasecorr_shift(F, M, subpixel)
H = size(F,1); W = size(F,2);
FF = fft2(F);
FM = fft2(M);
R = FF .* conj(FM);
R = R ./ max(abs(R), eps);
c = real(ifft2(R));

[peakVal, idx] = max(c(:));
[py, px] = ind2sub(size(c), idx);

dy0 = py - 1;
dx0 = px - 1;
if dy0 > H/2, dy0 = dy0 - H; end
if dx0 > W/2, dx0 = dx0 - W; end

dx = dx0;
dy = dy0;

mask = true(size(c));
win = 5;
r1 = max(1,py-win); r2 = min(H,py+win);
c1 = max(1,px-win); c2 = min(W,px+win);
mask(r1:r2, c1:c2) = false;
mu = mean(c(mask));
sd = std(c(mask));
psr = (peakVal - mu) / max(sd, eps);

if subpixel
    dx = dx0 + subpix_parabola(c(py, wrap(px-1,W)), c(py,px), c(py, wrap(px+1,W)));
    dy = dy0 + subpix_parabola(c(wrap(py-1,H),px), c(py,px), c(wrap(py+1,H),px));
end
end

function o = subpix_parabola(vL, v0, vR)
den = (vL - 2*v0 + vR);
if abs(den) < 1e-12
    o = 0;
else
    o = 0.5 * (vL - vR) / den;
end
end

function k = wrap(k, n)
k = mod(k-1, n) + 1;
end

function w = hann2(H, W)
wy = hann(H);
wx = hann(W);
w = wy * wx.';
end

function [x, y] = solve_global_positions_weighted(n, ii, jj, dx, dy, w)
m = numel(dx);
A = spalloc(2*m + 2, 2*n, 4*m + 2);
b = zeros(2*m + 2, 1);
ww = zeros(2*m + 2, 1);

row = 0;
for e = 1:m
    i = ii(e); j = jj(e);

    row = row + 1;
    A(row, 2*j-1) = 1;
    A(row, 2*i-1) = -1;
    b(row) = dx(e);
    ww(row) = w(e);

    row = row + 1;
    A(row, 2*j) = 1;
    A(row, 2*i) = -1;
    b(row) = dy(e);
    ww(row) = w(e);
end

% anchor tile 1 at (0,0)
A(end-1, 1) = 1; b(end-1) = 0; ww(end-1) = 50;
A(end,   2) = 1; b(end)   = 0; ww(end)   = 50;

Wmat = spdiags(sqrt(ww), 0, numel(ww), numel(ww));
p = lsqminnorm(Wmat*A, Wmat*b);

x = p(1:2:end);
y = p(2:2:end);
end

function w = linear_blend_weights(H, W)
[x, y] = meshgrid(1:W, 1:H);
dx = min(x-1, W-x);
dy = min(y-1, H-y);
w = min(dx, dy) + 1;
w = double(w);
end

function I2 = shift_bilinear_cached(I, fx, fy)
% cached grids for speed
persistent X Y lastH lastW
[H,W] = size(I);
if isempty(X) || isempty(Y) || lastH~=H || lastW~=W
    [X,Y] = meshgrid(1:W, 1:H);
    lastH = H; lastW = W;
end
I2 = interp2(X, Y, double(I), X - fx, Y - fy, 'linear', 0);
end

function [bitsPerSample, sampleFormat] = class_to_tiff_format(cls)
switch cls
    case 'uint8'
        bitsPerSample = 8;  sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint16'
        bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint32'
        bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.UInt;
    case 'int16'
        bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.Int;
    case 'single'
        bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'double'
        bitsPerSample = 64; sampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        error("Unsupported input class for TIFF writing: %s", cls);
end
end
