clear; clc; close all;
stitch_avgz_by_colrow_twochan_applycoords_SNAKEROW

function stitch_avgz_by_colrow_twochan_applycoords_SNAKEROW()
% stitch_avgz_by_colrow_twochan_applycoords_SNAKEROW
% ------------------------------------------------------------
% WORKING "nometa" stitcher adapted for your ScanImage naming + acquisition:
%
% 1) Tiles are AVGZ stacks (Z pages) already saved per tile in:
%      rawDir\ch3\avgz\*.tif
%      rawDir\ch1\avgz\*.tif
%
% 2) Tile positions are inferred from filename "colXX" and "rowYY".
%
% 3) Acquisition order is SNAKE ACROSS ROWS (row-major snake):
%      row00: col00 -> colNN
%      row01: colNN -> col00
%      row02: col00 -> colNN ...
%
% 4) Stitches refCh (ch3) FIRST to compute coordinates (x,y),
%    saves coords, then applies SAME coords to otherCh (ch1).
%
% 5) Lets you choose zUse to drop dirty end z (applies to both channels).
%
% Outputs:
%   - stitched_ch3_avgz.tif
%   - stitched_ch1_avgz.tif
%   - stitch_coords_from_ch3.mat
% ------------------------------------------------------------

%% ========================= USER SETTINGS =========================
rawDir = "Z:\p1yao\RZ\in vivo\sert_soma_g8s_retake";  % contains ch*/avgz

refChID   = 3;
otherChID = 1;

refAvgzDir   = fullfile(rawDir, sprintf("ch%d", refChID),   "avgz");
otherAvgzDir = fullfile(rawDir, sprintf("ch%d", otherChID), "avgz");

outDir = fullfile(rawDir, "matlab_stitch");
if ~isfolder(outDir), mkdir(outDir); end

outRef   = fullfile(outDir, sprintf("stitched_ch%d_avgz.tif", refChID));
outOther = fullfile(outDir, sprintf("stitched_ch%d_avgz.tif", otherChID));
coordMat = fullfile(outDir, sprintf("stitch_coords_from_ch%d.mat", refChID));

% Initial guess (same role as your working nometa script)
step_um = 600;
um_per_px_guess = 1.4815;

% IMPORTANT:
% If filename row index increases as you move "up" physically, set true.
row_increases_up = true;

% optional flips if mosaic mirrored
flipX = false;
flipY = false;

% registration
minOverlapFrac = 0.05;
psrMin = 5;
maxResidualPx = 150;
subpixel = true;

% blending
useLinearBlending = true;

% TIFF writing
rowsPerStrip = 64;

% ---------------- Z slice selection ----------------
% [] keeps all. Example: zUse = 2:(nZ-2)
zUse = 1:30;

%% ========================= DISCOVER FILES =========================
if ~isfolder(refAvgzDir),   error("Missing folder: %s", refAvgzDir); end
if ~isfolder(otherAvgzDir), error("Missing folder: %s", otherAvgzDir); end

refFilesAll = list_tifs(refAvgzDir);
if isempty(refFilesAll), error("No TIFFs in %s", refAvgzDir); end

% Parse col,row for each ref file
col = nan(numel(refFilesAll),1);
row = nan(numel(refFilesAll),1);
ok  = false(numel(refFilesAll),1);

for i = 1:numel(refFilesAll)
    [c,r,tf] = parse_colrow_from_name(refFilesAll(i));
    if tf
        col(i) = c; row(i) = r; ok(i) = true;
    end
end

refFilesAll = refFilesAll(ok);
col = col(ok);
row = row(ok);

if isempty(refFilesAll)
    error("No files matched colXX/rowYY in %s", refAvgzDir);
end

% Unique sorted bins -> infer grid size
uCol = sort(unique(col(:)));
uRow = sort(unique(row(:)));
nCols = numel(uCol);
nRows = numel(uRow);
nTiles = nRows*nCols;

fprintf("[grid] inferred nRows=%d nCols=%d nTiles=%d\n", nRows, nCols, nTiles);

% Map each file into (cBin, rBin) in [1..nCols],[1..nRows]
cBin = zeros(numel(refFilesAll),1);
rBin = zeros(numel(refFilesAll),1);
for i=1:numel(refFilesAll)
    cBin(i) = find(uCol==col(i),1);
    rBin(i) = find(uRow==row(i),1);
end

% Convert rBin to rB (bottom-based row index for geometry)
% rB=1 bottom, rB=nRows top
if row_increases_up
    rB = rBin;                 % smallest row is bottom
else
    rB = nRows - rBin + 1;     % smallest row is top
end

% tileAtRC(rB,c) holds index into refFilesAll (NOT tile index)
tileAtRC = nan(nRows, nCols);
for i=1:numel(refFilesAll)
    if ~isnan(tileAtRC(rB(i), cBin(i)))
        error("Duplicate mapping at (rB,c)=(%d,%d). Filenames not unique?", rB(i), cBin(i));
    end
    tileAtRC(rB(i), cBin(i)) = i;
end
if any(isnan(tileAtRC(:)))
    error("Missing tiles: some (row,col) missing in filenames.");
end

%% ========================= ACQUISITION ORDER: SNAKE ACROSS ROWS =========================
% Acquisition is row-major snake over *filename row order*:
% row00: col00->colNN, row01: colNN->col00, ...
% That means we loop rBin = 1..nRows, and alternate column direction.
tileOrder = zeros(nTiles,1);
k = 0;

for rBin_here = 1:nRows
    % map to bottom-based row index used in tileAtRC
    if row_increases_up
        rB_here = rBin_here;
    else
        rB_here = nRows - rBin_here + 1;
    end

    if mod(rBin_here,2)==1
        colsSnake = 1:nCols;        % row00 left->right
    else
        colsSnake = nCols:-1:1;     % row01 right->left
    end

    for c_here = colsSnake
        k = k + 1;
        tileOrder(k) = tileAtRC(rB_here, c_here); % index into refFilesAll
    end
end

tileFileRef = refFilesAll(tileOrder);

% Match other channel files by identical basename if possible; else fallback by col/row
tileFileOther = strings(nTiles,1);
for t=1:nTiles
    [~, bn, ext] = fileparts(tileFileRef(t));
    cand = fullfile(otherAvgzDir, bn + ext);
    if isfile(cand)
        tileFileOther(t) = cand;
    else
        % fallback by parsed col/row from ref
        [c0,r0,tf] = parse_colrow_from_name(tileFileRef(t));
        if ~tf, error("Cannot parse col/row from %s", tileFileRef(t)); end
        tileFileOther(t) = find_by_colrow(otherAvgzDir, c0, r0);
    end
end
if any(tileFileOther=="")
    error("Some ch%d tiles missing in %s", otherChID, otherAvgzDir);
end

fprintf("[order] using row-major snake acquisition order\n");
disp("First 5 tiles:");
disp(tileFileRef(1:min(5,end)));

%% ========================= READ SIZE/CLASS/NZ =========================
[H0, W0, inClass] = read_first_page_size_class(tileFileRef(1));
nZ0 = count_tiff_pages_file(tileFileRef(1));
fprintf("[tile] size=%dx%d, nZ=%d, class=%s\n", H0, W0, nZ0, inClass);

% Validate all tiles consistent + finalize zUse
if isempty(zUse)
    zUse = 1:nZ0;
else
    zUse = zUse(:).';
    zUse = zUse(zUse>=1 & zUse<=nZ0);
    if isempty(zUse), error("zUse empty after clipping."); end
end
fprintf("[zUse] %d slices: [%d .. %d] (nZ=%d)\n", numel(zUse), min(zUse), max(zUse), nZ0);

for t=1:nTiles
    [Hk,Wk,clsK] = read_first_page_size_class(tileFileRef(t));
    if Hk~=H0 || Wk~=W0 || clsK~=string(inClass)
        error("Tile size/class mismatch: %s", tileFileRef(t));
    end
    nZk = count_tiff_pages_file(tileFileRef(t));
    if nZk~=nZ0
        error("Tile nZ mismatch: %s (nZ=%d vs %d)", tileFileRef(t), nZk, nZ0);
    end
end

%% ========================= BUILD SUM PROJ (REF) =========================
% No need to keep file handles open; just compute sum projections and close.
sumProj = cell(nTiles,1);
for t=1:nTiles
    tt = Tiff(tileFileRef(t),'r');
    c = onCleanup(@() tt.close()); %#ok<NASGU>

    S = zeros(H0,W0,'double');
    for iz = zUse
        tt.setDirectory(iz);
        S = S + double(tt.read());
    end
    sumProj{t} = S;
end

%% ========================= INITIAL GUESSES x0,y0 =========================
step_px = step_um / um_per_px_guess;

% Create inverse maps: for each refFilesAll index -> (rB,c)
inv_rB = nan(numel(refFilesAll),1);
inv_c  = nan(numel(refFilesAll),1);
for rr=1:nRows
    for cc=1:nCols
        ii = tileAtRC(rr,cc);
        inv_rB(ii) = rr;
        inv_c(ii)  = cc;
    end
end

% For each tile index t (acquisition order), get its grid coords (rB,c)
c_of_t  = zeros(nTiles,1);
rB_of_t = zeros(nTiles,1);
for t=1:nTiles
    ii = tileOrder(t);       % index into refFilesAll
    c_of_t(t)  = inv_c(ii);
    rB_of_t(t) = inv_rB(ii);
end

x0 = zeros(nTiles,1);
y0 = zeros(nTiles,1);
for t=1:nTiles
    c = c_of_t(t);
    rB_here = rB_of_t(t);              % bottom-based
    rTop = nRows - rB_here + 1;        % convert to top-down image coords
    x0(t) = (c-1)    * step_px;
    y0(t) = (rTop-1) * step_px;
end
if flipX, x0 = max(x0) - x0; end
if flipY, y0 = max(y0) - y0; end

%% ========================= BUILD NEIGHBOR EDGES (4-neighborhood) =========================
% Build idxGrid(rB,c) = tile index t
idxGrid = zeros(nRows,nCols);
for t=1:nTiles
    idxGrid(rB_of_t(t), c_of_t(t)) = t;
end

edges = []; % [i j dxGuess dyGuess]
for rB_here = 1:nRows
    for c_here = 1:nCols
        i = idxGrid(rB_here, c_here);

        if c_here < nCols
            j = idxGrid(rB_here, c_here+1);
            edges(end+1,:) = [i j (x0(j)-x0(i)) (y0(j)-y0(i))]; %#ok<AGROW>
        end
        if rB_here < nRows
            j = idxGrid(rB_here+1, c_here); % neighbor above
            edges(end+1,:) = [i j (x0(j)-x0(i)) (y0(j)-y0(i))]; %#ok<AGROW>
        end
    end
end
fprintf("[reg] edges=%d\n", size(edges,1));

%% ========================= PAIRWISE REG ON REF SUM PROJ =========================
ii = edges(:,1); jj = edges(:,2);
dxG = edges(:,3); dyG = edges(:,4);

dxE = dxG; dyE = dyG;
wE  = 0.2 * ones(size(dxG));

for e = 1:numel(dxG)
    i = ii(e); j = jj(e);
    [dx, dy, ok, psr] = refine_shift_phasecorr(sumProj{i}, sumProj{j}, dxG(e), dyG(e), ...
        minOverlapFrac, psrMin, maxResidualPx, subpixel);

    if ok
        dxE(e) = dx;
        dyE(e) = dy;
        wE(e)  = min(20, max(1, psr));
    end
end
fprintf("[reg] ok edges=%d/%d\n", nnz(wE>1), numel(wE));

%% ========================= GLOBAL OPT (coords) =========================
[x, y] = solve_global_positions_weighted(nTiles, ii, jj, dxE, dyE, wE);
x = x - min(x);
y = y - min(y);

outW = ceil(max(x) + W0);
outH = ceil(max(y) + H0);
fprintf("[stitch] mosaic=%dx%d, zUsed=%d/%d\n", outH, outW, numel(zUse), nZ0);

save(coordMat, 'x','y','tileFileRef','tileFileOther','c_of_t','rB_of_t','nRows','nCols','zUse','H0','W0','outH','outW', ...
    'step_um','um_per_px_guess','row_increases_up','flipX','flipY');
fprintf("[save] %s\n", coordMat);

%% ========================= BLENDING WEIGHTS =========================
if useLinearBlending
    wTile = linear_blend_weights(H0, W0);
else
    wTile = ones(H0,W0,'double');
end

%% ========================= PRECOMPUTE TILE PLACEMENTS =========================
place = struct([]);
for t = 1:nTiles
    xi = x(t); yi = y(t);
    xInt = floor(xi); yInt = floor(yi);
    fx = xi - xInt; fy = yi - yInt;

    r0 = yInt + 1; c0 = xInt + 1;
    rr = r0:(r0+H0-1);
    cc = c0:(c0+W0-1);

    if subpixel && (abs(fx)>1e-6 || abs(fy)>1e-6)
        w2 = shift_bilinear(wTile, fx, fy);
    else
        w2 = wTile;
    end

    place(t).fx = fx;
    place(t).fy = fy;
    place(t).rr = rr;
    place(t).cc = cc;
    place(t).w2 = w2;
end

%% ========================= STITCH REF (ch3) =========================
if isfile(outRef), delete(outRef); end
write_stitched_from_tilefiles(outRef, tileFileRef, place, outH, outW, H0, W0, zUse, inClass, rowsPerStrip, ...
    sprintf("coords from ch%d; row-snake; zUse=[%d..%d]", refChID, min(zUse), max(zUse)));
fprintf("[DONE] %s\n", outRef);

%% ========================= STITCH OTHER (ch1) SAME coords =========================
if isfile(outOther), delete(outOther); end
write_stitched_from_tilefiles(outOther, tileFileOther, place, outH, outW, H0, W0, zUse, inClass, rowsPerStrip, ...
    sprintf("coords from ch%d applied to ch%d; row-snake; zUse=[%d..%d]", refChID, otherChID, min(zUse), max(zUse)));
fprintf("[DONE] %s\n", outOther);

end

%% ===================== helpers =====================

function files = list_tifs(folder)
c = [dir(fullfile(folder,"*.tif")); dir(fullfile(folder,"*.tiff"))];
files = strings(numel(c),1);
for k=1:numel(c)
    files(k) = string(fullfile(c(k).folder, c(k).name));
end
end

function [col,row,tf] = parse_colrow_from_name(fullpath)
[~, bn, ~] = fileparts(fullpath);
tokC = regexp(bn, 'col(\d+)', 'tokens', 'once');
tokR = regexp(bn, 'row(\d+)', 'tokens', 'once');
tf = ~isempty(tokC) && ~isempty(tokR);
if ~tf
    col = NaN; row = NaN;
else
    col = str2double(tokC{1});
    row = str2double(tokR{1});
end
end

function fn = find_by_colrow(folder, col, row)
pat1 = sprintf("*col%02d*row%02d*.tif*", col, row);
c = dir(fullfile(folder, pat1));
if isempty(c)
    % if col/row are not zero-padded in names
    pat2 = sprintf("*col%d*row%d*.tif*", col, row);
    c = dir(fullfile(folder, pat2));
end
if isempty(c)
    fn = "";
else
    fn = string(fullfile(c(1).folder, c(1).name));
end
end

function [H, W, cls] = read_first_page_size_class(fn)
tt = Tiff(fn,'r');
tt.setDirectory(1);
I = tt.read();
tt.close();
[H,W] = size(I);
cls = class(I);
end

function n = count_tiff_pages_file(fn)
tt = Tiff(fn,'r');
tt.setDirectory(1);
n = 1;
while ~tt.lastDirectory()
    tt.nextDirectory();
    n = n + 1;
end
tt.close();
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
    a_c = (dx+1):W; b_c = 1:(W-dx);
else
    a_c = 1:(W+dx); b_c = (-dx+1):W;
end

if dy >= 0
    a_r = (dy+1):H; b_r = 1:(H-dy);
else
    a_r = 1:(H+dy); b_r = (-dy+1):H;
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

dx = dx0; dy = dy0;

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
if abs(den) < 1e-12, o = 0; else, o = 0.5 * (vL - vR) / den; end
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
[xg, yg] = meshgrid(1:W, 1:H);
dx = min(xg-1, W-xg);
dy = min(yg-1, H-yg);
w = min(dx, dy) + 1;
w = double(w);
end

function I2 = shift_bilinear(I, fx, fy)
[H,W] = size(I);
[xg, yg] = meshgrid(1:W, 1:H);
I2 = interp2(xg, yg, double(I), xg - fx, yg - fy, 'linear', 0);
end

function [bitsPerSample, sampleFormat] = class_to_tiff_format(cls)
switch char(cls)
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
        error("Unsupported TIFF class: %s", cls);
end
end

function write_stitched_from_tilefiles(outTif, tileFiles, place, outH, outW, H0, W0, zUse, inClass, rowsPerStrip, softwareStr)
nTiles = numel(tileFiles);

% open all tiles (avgz stacks are small; this is fine)
tileT = cell(nTiles,1);
for t=1:nTiles
    tileT{t} = Tiff(tileFiles(t),'r');
end
cleanup = onCleanup(@() cellfun(@(tt) safe_close(tt), tileT)); %#ok<NASGU>

[bitsPerSample, sampleFormat] = class_to_tiff_format(inClass);

nZout = numel(zUse);
bytesPerPixel = bitsPerSample/8;
estBytes = double(outH) * double(outW) * double(nZout) * double(bytesPerPixel);
useBigTiff = estBytes > 3.5e9;

if isfile(outTif), delete(outTif); end
if useBigTiff
    tout = Tiff(outTif, 'w8');
else
    tout = Tiff(outTif, 'w');
end
cleanup2 = onCleanup(@() tout.close()); %#ok<NASGU>

tag.ImageLength = outH;
tag.ImageWidth  = outW;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.BitsPerSample = bitsPerSample;
tag.SampleFormat = sampleFormat;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = rowsPerStrip;
tag.Software = softwareStr;

for kz = 1:nZout
    iz = zUse(kz);
    fprintf("[write] z %d/%d (src %d)\n", kz, nZout, iz);

    acc  = zeros(outH, outW, 'double');
    wacc = zeros(outH, outW, 'double');

    for t=1:nTiles
        tileT{t}.setDirectory(iz);
        I = double(tileT{t}.read());

        fx = place(t).fx; fy = place(t).fy;
        if (abs(fx)>1e-6 || abs(fy)>1e-6)
            I2 = shift_bilinear(I, fx, fy);
        else
            I2 = I;
        end

        rr = place(t).rr;
        cc = place(t).cc;
        w2 = place(t).w2;

        acc(rr,cc)  = acc(rr,cc)  + I2 .* w2;
        wacc(rr,cc) = wacc(rr,cc) + w2;
    end

    out = acc ./ max(wacc, eps);
    out = cast(round(out), inClass);

    tout.setTag(tag);
    tout.write(out);
    if kz < nZout
        tout.writeDirectory();
    end
end
end

function safe_close(tt)
try, tt.close(); catch, end %#ok<CTCH>
end
