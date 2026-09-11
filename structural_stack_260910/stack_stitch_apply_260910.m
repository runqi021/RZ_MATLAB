function out = stack_stitch_apply_260910(coordFile, varargin)
%STACK_STITCH_APPLY_260910  Render a stitched VOLUME from curated tile coords.
%
%   out = stack_stitch_apply_260910(coordFile)
%   out = stack_stitch_apply_260910(coordFile, 'Channels', [1 3], ...)
%
% Takes the `stack_stitch_coords.mat` written by stack_stitch_gui_260910 and
% renders the full-resolution mosaic for every requested channel.  Separate from
% the GUI on purpose: the placement is the only thing that needs a human, so
% once it is curated this runs headless and can be re-run for another channel
% without touching the coords.
%
% BLENDING IS A LINEAR FEATHER, NOT A LAST-WRITER-WINS PASTE.  Each tile carries
% a weight that ramps from ~0 at its own edge to 1 `Feather` microns inside, and
% the output is the weighted mean.  A hard paste puts a visible step everywhere
% the tiles disagree -- and they do disagree, because the two were acquired
% minutes apart at different stage positions with the same power ramp, so the
% overlap strip sits at a different distance from each tile's optical axis.
%
% THE WEIGHTS ARE 2-D.  They depend on x and y only, never on z, so the
% denominator is one [H x W] plane rather than a whole volume -- which is what
% keeps this inside RAM.  The numerator is single [H x W x Z]: about 1.2 GB for
% a 3-tile strip, so it is accumulated one tile at a time and divided plane by
% plane on the way out.
%
% Z OFFSETS ARE IN WHOLE PLANES.  The volumes are 1 um isotropic, so a plane IS
% a micron and no interpolation is needed or wanted: resampling in z to chase a
% sub-micron offset would blur every tile to fix a shift smaller than the axial
% PSF.
%
% OPTIONS
%   'Channels'  []      SI channel numbers ([] = every channel in the tiles)
%   'Feather'   30      blend ramp width, um
%   'OutDir'    ''      default <session root>\matlab_stitch
%   'SaveProj'  true    also write max and avg XY projections of the mosaic
%
% Runqi Zhang / 2026-09-10

p = inputParser;
p.addParameter('Channels', []);
p.addParameter('Feather',  30, @(x)isscalar(x) && x >= 0);
p.addParameter('OutDir',   '');
p.addParameter('SaveProj', true);
p.parse(varargin{:});
opt = p.Results;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);

C = load(coordFile);
assert(isfield(C, 'stitch'), '%s is not a stack_stitch coords file', coordFile);
St = C.stitch;

outDir = opt.OutDir;
if isempty(outDir), outDir = fullfile(St.root, 'matlab_stitch'); end
if ~isfolder(outDir), mkdir(outDir); end

chList = opt.Channels;
if isempty(chList), chList = St.channels; end

n  = numel(St.metaFiles);
pt = St.pxOut;

fprintf('\n================ stack_stitch_apply_260910 ================\n');
fprintf('coords : %s\n', coordFile);
fprintf('tiles  : %d   voxel %g um   feather %g um\n', n, pt, opt.Feather);

%% ---------------------------- canvas geometry ------------------------------
% x/y are the top-left corner of each tile in output pixels; z is a whole-plane
% offset.  Shift everything so the mosaic starts at (1,1,1).
x = St.x(:) + St.manX(:);
y = St.y(:) + St.manY(:);
z = St.manZ(:);
h = St.tileH(:); w = St.tileW(:); d = St.tileD(:);

x = round(x - min(x));  y = round(y - min(y));  z = round(z - min(z));
W = max(x + w);  H = max(y + h);  D = max(z + d);
fprintf('canvas : %d x %d x %d px  (%g x %g x %g um)\n', H, W, D, H*pt, W*pt, D*pt);
for i = 1:n
    fprintf('  tile %d  %-22s  x %5d  y %5d  z %4d\n', i, St.names{i}, x(i), y(i), z(i));
end

%% -------------------------------- weights ----------------------------------
% One ramp per tile SIZE (they are usually identical, so this is computed once
% and reused).  max(...,1e-3) keeps a pixel that only one tile covers from
% getting a zero denominator at the very edge of the mosaic.
den = zeros(H, W, 'single');
Wt  = cell(n,1);
featherPx = max(1, opt.Feather / pt);
for i = 1:n
    Wt{i} = feather_ramp(h(i), w(i), featherPx);
    den(y(i)+(1:h(i)), x(i)+(1:w(i))) = den(y(i)+(1:h(i)), x(i)+(1:w(i))) + Wt{i};
end
den(den == 0) = 1;                      % untouched canvas stays 0/1 = 0

%% ------------------------------ per channel --------------------------------
out = struct('ch', {}, 'file', {});
for ch = chList(:).'
    fprintf('\n--- channel %d ---\n', ch);
    num = zeros(H, W, D, 'single');
    t0 = tic;
    for i = 1:n
        tif = tile_stack_path(St.metaFiles{i}, ch);
        V   = single(read_tiff_stack(tif));
        assert(size(V,1) == h(i) && size(V,2) == w(i) && size(V,3) == d(i), ...
            'tile %d is [%s] but the coords say [%d %d %d]', i, ...
            num2str(size(V)), h(i), w(i), d(i));
        rr = y(i)+(1:h(i));  cc = x(i)+(1:w(i));  zz = z(i)+(1:d(i));
        num(rr,cc,zz) = num(rr,cc,zz) + V .* Wt{i};
        fprintf('  tile %d placed (%.0f s)\n', i, toc(t0));
        clear V
    end

    fn = fullfile(outDir, sprintf('stitched_ch%d_volume.tif', ch));
    write_tiff_int16_planes(fn, num, den);
    dd = dir(fn);
    fprintf('  %s  (%.2f GB)\n', fn, dd.bytes/2^30);
    out(end+1) = struct('ch', ch, 'file', fn); %#ok<AGROW>

    if opt.SaveProj
        num = num ./ den;
        mx = max(num, [], 3);  av = mean(num, 3);
        imwrite(uint16(max(0, mx)), fullfile(outDir, sprintf('stitched_ch%d_maxproj.tif', ch)));
        imwrite(uint16(max(0, av)), fullfile(outDir, sprintf('stitched_ch%d_avgproj.tif', ch)));
        fprintf('  + max / avg XY projections\n');
    end
    clear num
end

save(fullfile(outDir, 'stack_stitch_applied.mat'), 'St', 'x', 'y', 'z', 'H', 'W', 'D', 'opt');
fprintf('\ndone -> %s\n', outDir);
end

% =============================================================================
function R = feather_ramp(h, w, fpx)
%FEATHER_RAMP  1 in the middle, ramping to ~0 at the tile border.
% Separable (the min of an x ramp and a y ramp) rather than a true distance
% transform: on a rectangle they agree everywhere except within fpx of a corner,
% and it needs no toolbox.
cx = min(1:w, w:-1:1);  cy = min(1:h, h:-1:1);
rx = min(1, single(cx) / fpx);
ry = min(1, single(cy) / fpx);
R  = max(min(ry(:), rx(:).'), 1e-3);
end

function tif = tile_stack_path(metaFile, ch)
S = load(metaFile, 'meta');
ci = find(S.meta.channels == ch, 1);
assert(~isempty(ci), 'channel %d not in %s', ch, metaFile);
tif = S.meta.outTif{ci};
if ~isfile(tif)
    [~, nm, ex] = fileparts(tif);
    tif = fullfile(fileparts(metaFile), [nm ex]);
end
assert(isfile(tif), 'missing tile stack %s', tif);
end

function V = read_tiff_stack(fn)
t = Tiff(fn, 'r'); c = onCleanup(@() close(t));
A = t.read();  n = 1;
while ~t.lastDirectory(), t.nextDirectory(); n = n + 1; end
V = zeros(size(A,1), size(A,2), n, 'like', A);
t.setDirectory(1);
for k = 1:n
    V(:,:,k) = t.read();
    if k < n, t.nextDirectory(); end
end
end

function write_tiff_int16_planes(fn, num, den)
%WRITE_TIFF_INT16_PLANES  num ./ den, one plane at a time, as a BigTIFF.
% Dividing on the way out rather than in place avoids a second full-size copy of
% a volume that is already over a gigabyte.
[ny, nx, nz] = size(num);
if isfile(fn), delete(fn); end
t = Tiff(fn, 'w8');  c = onCleanup(@() close(t));
tag.ImageLength = ny;  tag.ImageWidth = nx;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.BitsPerSample = 16;  tag.SampleFormat = Tiff.SampleFormat.Int;
tag.SamplesPerPixel = 1;
tag.RowsPerStrip        = ny;   % ONE strip per page.  MATLAB's Tiff writer pads
                                % the final strip out to a full RowsPerStrip
                                % instead of writing a short one, so a 455-row
                                % page at 64 rows/strip carries 57 rows of zeros
                                % -- measured at 52 237 B/page, 11%% of the file.
                                % A single strip has no partial strip to pad.

tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.Software = 'stack_stitch_apply_260910';
for k = 1:nz
    P = num(:,:,k) ./ den;
    t.setTag(tag);
    t.write(int16(round(max(min(P, 32767), -32768))));
    if k < nz, t.writeDirectory(); end
end
end
