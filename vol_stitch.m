stitch_snake_7x7_imagejstyle_nometa

%%
function stitch_snake_7x7_imagejstyle_nometa()
% Grid/Collection Stitching style (no metadata):
% - tiles: tile_01 ... tile_49 (7x7)
% - snake order, starting from left-bottom
% - initial spacing: 500 um (converted to pixels using um_per_px_guess)
% - registration on SUM projection (sum over Z pages)
% - apply solved XY to entire tile stack
% - fusion: linear blending weights (distance-to-border)
%
% Output: rawDir\matlab_stitch\stitched_stack_imagejstyle_nometa.tif

%% ---------------- USER ----------------
rawDir = "E:\live_251104_definitely_not_sst_tdtomato_fitc\raw\groupZ";
outDir = fullfile(rawDir, "matlab_stitch");
outTif = fullfile(outDir, "stitched_stack_imagejstyle_nometa.tif");

nRows = 6; nCols = 6;
nTiles = nRows * nCols;

% initial condition
step_um = 600;

% IMPORTANT (since no metadata): microns per pixel
um_per_px_guess = 1.4815;   % <-- SET THIS (e.g., 1.48 for ~758um FOV @ 512px)

% optional flips if mosaic mirrored
flipX = false;
flipY = false;

% registration
minOverlapFrac = 0.05;    % overlap needed (fraction of width/height) to attempt registration
psrMin = 5;               % peak-to-sidelobe ratio threshold (lower if noisy)
maxResidualPx = 150;       % reject crazy residuals (px)
subpixel = true;          % subpixel peak refine + apply subpixel shifts (interp2)

% blending
useLinearBlending = true; % "linear blending" like ImageJ option :contentReference[oaicite:1]{index=1}

% TIFF writing
rowsPerStrip = 64;

%% ---------------- FIND TILE FILES ----------------
if ~isfolder(rawDir), error("Not a folder: %s", rawDir); end
if ~isfolder(outDir), mkdir(outDir); end
if isfile(outTif), delete(outTif); end

tileFile = strings(nTiles,1);
for i = 1:nTiles
    tileFile(i) = find_tile_file(rawDir, i);
end
fprintf("[stitch] Found all %d tiles.\n", nTiles);

%% ---------------- OPEN TILES, GET SIZE + nZ, BUILD SUM PROJ ----------------
tile = struct([]);
for t = 1:nTiles
    fn = tileFile(t);
    tt = Tiff(fn, 'r');
    % keep handle open (we'll reuse later)
    tile(t).t = tt;
    tile(t).fn = fn;

    % read first slice to get size/class
    tt.setDirectory(1);
    I0 = tt.read();
    [H,W] = size(I0);

    % count pages (nZ)
    nZ = count_tiff_pages(tt);

    tile(t).H = H;
    tile(t).W = W;
    tile(t).nZ = nZ;

    % sum projection over Z (in double)
    sumProj = zeros(H,W,'double');
    for iz = 1:nZ
        tt.setDirectory(iz);
        sumProj = sumProj + double(tt.read());
    end
    tile(t).sumProj = sumProj;

    if t==1
        inClass = class(I0);
        fprintf("[tile] size=%dx%d, nZ=%d, class=%s\n", H,W,nZ,inClass);
    else
        if tile(t).H~=tile(1).H || tile(t).W~=tile(1).W || tile(t).nZ~=tile(1).nZ
            error("Tiles differ in size or nZ. Split into separate folders.");
        end
    end
end
cleanupTiles = onCleanup(@() close_all_tiles(tile)); %#ok<NASGU>

H0 = tile(1).H; W0 = tile(1).W; nZ0 = tile(1).nZ;

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

% initial positions (pixels)
step_px = step_um / um_per_px_guess;
x0 = zeros(nTiles,1);
y0 = zeros(nTiles,1);

for rB = 1:nRows
    rTop = nRows - rB + 1;  % convert bottom-based to top-based image rows
    for c = 1:nCols
        t = idxRC(rB,c);
        x0(t) = (c-1) * step_px;
        y0(t) = (rTop-1) * step_px;
    end
end

if flipX, x0 = max(x0) - x0; end
if flipY, y0 = max(y0) - y0; end

%% ---------------- BUILD NEIGHBOR EDGES (4-neighborhood) ----------------
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

%% ---------------- PAIRWISE REG ON SUM PROJ ----------------
ii = edges(:,1); jj = edges(:,2);
dxG = edges(:,3); dyG = edges(:,4);

dxE = dxG; dyE = dyG;
wE  = 0.2 * ones(size(dxG));    % low weight fallback if registration fails

for e = 1:numel(dxG)
    i = ii(e); j = jj(e);

    [dx, dy, ok, psr] = refine_shift_phasecorr(tile(i).sumProj, tile(j).sumProj, dxG(e), dyG(e), ...
        minOverlapFrac, psrMin, maxResidualPx, subpixel);

    if ok
        dxE(e) = dx;
        dyE(e) = dy;
        wE(e)  = min(20, max(1, psr));  % correlation confidence weight
    end
end
fprintf("[reg] ok edges=%d/%d\n", nnz(wE>1), numel(wE));

%% ---------------- GLOBAL OPT (weighted LSQ) ----------------
[x, y] = solve_global_positions_weighted(nTiles, ii, jj, dxE, dyE, wE);

x = x - min(x);
y = y - min(y);

outW = ceil(max(x) + W0);
outH = ceil(max(y) + H0);
fprintf("[stitch] mosaic=%dx%d, nZ=%d\n", outH, outW, nZ0);

%% ---------------- BLENDING WEIGHTS ----------------
if useLinearBlending
    wTile = linear_blend_weights(H0, W0);
else
    wTile = ones(H0,W0,'double');
end

%% ---------------- OUTPUT TIFF SETUP ----------------
% infer class from first tile, first page
tile(1).t.setDirectory(1);
img0 = tile(1).t.read();
inClass = class(img0);

[bitsPerSample, sampleFormat] = class_to_tiff_format(inClass);

bytesPerPixel = bitsPerSample/8;
estBytes = double(outH) * double(outW) * double(nZ0) * double(bytesPerPixel);
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
tag.Software = 'MATLAB snake stitch (no meta; sum-proj reg; linear blend)';

%% ---------------- STITCH SLICE-BY-SLICE ----------------
for iz = 1:nZ0
    fprintf("[write] z %d/%d\n", iz, nZ0);

    acc  = zeros(outH, outW, 'double');
    wacc = zeros(outH, outW, 'double');

    for t = 1:nTiles
        tile(t).t.setDirectory(iz);
        I = double(tile(t).t.read());

        xi = x(t); yi = y(t);
        xInt = floor(xi); yInt = floor(yi);
        fx = xi - xInt; fy = yi - yInt;

        if subpixel && (abs(fx) > 1e-6 || abs(fy) > 1e-6)
            I2 = shift_bilinear(I, fx, fy);
            w2 = shift_bilinear(wTile, fx, fy);
        else
            I2 = I;
            w2 = wTile;
        end

        r0 = yInt + 1;
        c0 = xInt + 1;
        rr = r0:(r0+H0-1);
        cc = c0:(c0+W0-1);

        acc(rr,cc)  = acc(rr,cc)  + I2 .* w2;
        wacc(rr,cc) = wacc(rr,cc) + w2;
    end

    out = acc ./ max(wacc, eps);
    out = cast(round(out), inClass);

    tout.setTag(tag);
    tout.write(out);
    if iz < nZ0
        tout.writeDirectory();
    end
end

fprintf("\n[DONE] %s\n", outTif);

end

%% ===================== helpers =====================

function fn = find_tile_file(root, idx)
% find tile_XX.tif/f in root (recursive)
pat = sprintf("tile_%02d", idx);
cands = [dir(fullfile(root,"**",pat+"*.tif")); dir(fullfile(root,"**",pat+"*.tiff"))];
if isempty(cands)
    error("Missing file for %s under %s", pat, root);
end
% pick first (if duplicates exist, you should clean the folder)
fn = string(fullfile(cands(1).folder, cands(1).name));
end

function n = count_tiff_pages(t)
% count directories without imfinfo
t.setDirectory(1);
n = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    n = n + 1;
end
% return to first
t.setDirectory(1);
end

function close_all_tiles(tile)
for k = 1:numel(tile)
    try, tile(k).t.close(); catch, end %#ok<CTCH>
end
end

function [dx, dy, ok, psr] = refine_shift_phasecorr(A, B, dxGuess, dyGuess, minOverlapFrac, psrMin, maxResidualPx, subpixel)
% Align B to A. Returns refined shift (dx,dy) for B relative to A.

H = size(A,1); W = size(A,2);

% overlap size from guess
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
% distance-to-border weights (linear ramp)
[x, y] = meshgrid(1:W, 1:H);
dx = min(x-1, W-x);
dy = min(y-1, H-y);
w = min(dx, dy) + 1;
w = double(w);
end

function I2 = shift_bilinear(I, fx, fy)
[H,W] = size(I);
[x, y] = meshgrid(1:W, 1:H);
I2 = interp2(x, y, double(I), x - fx, y - fy, 'linear', 0);
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
