function stitch_apply_coords_volume(datasetRoot, chan, zUse)
%STITCH_APPLY_COORDS_VOLUME  Apply manual stitch coords to a z-stack -> VOLUME.
%
%   stitch_manual_gui solves tile positions and writes matlab_stitch/
%   stitch_manual_coords.mat, but only saves a 2D AVG-PROJECTION. This applies
%   those exact pixel positions to the multi-page tiles and writes a stitched
%   MULTI-PAGE VOLUME (Z preserved) -- the input you need for an X-Z / depth
%   projection (see phys_zStack/vessel_XZ_projection_260529.m).
%
%   Blends each Z slice independently with the GUI's linear-blend weights.
%
%   stitch_apply_coords_volume                      % 260522 root, same chan as coords, all Z
%   stitch_apply_coords_volume(root)                % explicit root
%   stitch_apply_coords_volume(root, 3)             % render ch3 (swap token if coords were ch1)
%   stitch_apply_coords_volume(root, 3, 2:28)       % subset of Z slices
%
%   Output: matlab_stitch/stitched_ch<chan>_volume.tif  (multi-page, BigTIFF if needed)
%
%   Standalone. MATLAB R2021a+. Runqi Zhang / 2026.

if nargin<1 || isempty(datasetRoot)
    datasetRoot = "C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\map";
end
if nargin<3, zUse = []; end
datasetRoot = char(datasetRoot);

coordDir = fullfile(datasetRoot,'matlab_stitch');
coordMat = fullfile(coordDir,'stitch_manual_coords.mat');
assert(isfile(coordMat), 'Coords not found: %s\nRun stitch_manual_gui first.', coordMat);

p = load(coordMat);
need = {'x','y','outH','outW','H0','W0','files'};
for k=1:numel(need)
    assert(isfield(p,need{k}), 'Coords missing field "%s" -- re-save from the GUI.', need{k});
end
x=p.x(:); y=p.y(:); outH=p.outH; outW=p.outW; H0=p.H0; W0=p.W0;
srcFiles = string(p.files(:));
nTiles = numel(srcFiles);

% --- channel the coords were solved on (from filename token) ---
tok = regexp(char(srcFiles(1)), 'ch(\d+)', 'tokens', 'once');
srcCh = []; if ~isempty(tok), srcCh = str2double(tok{1}); end
if nargin<2 || isempty(chan)
    chan = srcCh;
    assert(~isempty(chan), 'Could not parse channel from coords; pass chan explicitly.');
end

% --- map tiles to requested channel (swap token if different) ---
if isempty(srcCh) || chan==srcCh
    tileFiles = srcFiles;
else
    srcTok = sprintf('ch%d', srcCh);  dstTok = sprintf('ch%d', chan);
    tileFiles = strings(nTiles,1);
    for t=1:nTiles
        f = char(srcFiles(t));
        g = strrep(f, [filesep srcTok filesep], [filesep dstTok filesep]); % ch1\ -> ch3\
        g = regexprep(g, ['_' srcTok '_'], ['_' dstTok '_']);              % _ch1_ -> _ch3_
        assert(isfile(g), 'Destination tile not found:\n  %s\n(from %s)', g, f);
        tileFiles(t) = string(g);
    end
end

% --- probe class + nZ, finalize zUse ---
I0 = imread(char(tileFiles(1)),1); inClass = class(I0);
nZ0 = count_tiff_pages(tileFiles(1));
for t=1:nTiles
    nZk = count_tiff_pages(tileFiles(t));
    assert(nZk==nZ0, 'Tile nZ mismatch: %s (nZ=%d vs %d)', tileFiles(t), nZk, nZ0);
end
if isempty(zUse), zUse = 1:nZ0; else
    zUse = zUse(:).'; zUse = zUse(zUse>=1 & zUse<=nZ0);
    assert(~isempty(zUse), 'zUse empty after clipping to 1..%d.', nZ0);
end
nZout = numel(zUse);
fprintf('[vol] %d tiles, mosaic %dx%d, Z=%d/%d, class=%s\n', ...
        nTiles, outH, outW, nZout, nZ0, inClass);

% --- open all tiles ---
T = cell(nTiles,1);
for t=1:nTiles, T{t} = Tiff(char(tileFiles(t)),'r'); end
cleanT = onCleanup(@() cellfun(@safeClose, T)); %#ok<NASGU>

% --- precompute placement windows + blend weights ---
wT  = linear_blend_weights(H0,W0);
rr = cell(nTiles,1); cc = cell(nTiles,1);
for t=1:nTiles
    r0=round(y(t)); c0=round(x(t));
    r0=max(0,min(r0,outH-H0)); c0=max(0,min(c0,outW-W0));
    rr{t}=r0+1:r0+H0; cc{t}=c0+1:c0+W0;
end

% --- output TIFF (BigTIFF if large) ---
out = fullfile(coordDir, sprintf('stitched_ch%d_volume.tif', chan));
if isfile(out), delete(out); end
[bps,sf] = class_to_tiff_format(inClass);
estBytes = double(outH)*double(outW)*double(nZout)*double(bps/8);
tout = Tiff(out, ternary(estBytes>3.5e9,'w8','w'));
cleanO = onCleanup(@() tout.close()); %#ok<NASGU>
tag.ImageLength=outH; tag.ImageWidth=outW;
tag.Photometric=Tiff.Photometric.MinIsBlack; tag.SamplesPerPixel=1;
tag.BitsPerSample=bps; tag.SampleFormat=sf;
tag.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
tag.Compression=Tiff.Compression.None; tag.RowsPerStrip=64;
tag.Software='stitch_apply_coords_volume';

for kz=1:nZout
    iz = zUse(kz);
    acc = zeros(outH,outW,'double'); wac = zeros(outH,outW,'double');
    for t=1:nTiles
        T{t}.setDirectory(iz);
        img = double(T{t}.read());
        acc(rr{t},cc{t}) = acc(rr{t},cc{t}) + img.*wT;
        wac(rr{t},cc{t}) = wac(rr{t},cc{t}) + wT;
    end
    M = cast(round(acc./max(wac,eps)), inClass);
    tout.setTag(tag); tout.write(M);
    if kz<nZout, tout.writeDirectory(); end
    fprintf('  z %d/%d (src %d)\n', kz, nZout, iz);
end
fprintf('[DONE] %s\n', out);
end

%% ===================== helpers =====================
function n = count_tiff_pages(file)
tf=Tiff(char(file),'r'); c=onCleanup(@() safeClose(tf)); %#ok<NASGU>
n=1; while ~tf.lastDirectory(), tf.nextDirectory(); n=n+1; end
end

function w=linear_blend_weights(H,W)
[xg,yg]=meshgrid(1:W,1:H);
w=double(min(min(xg-1,W-xg),min(yg-1,H-yg))+1);
end

function [bps,sf]=class_to_tiff_format(cls)
switch char(cls)
    case 'uint8',  bps=8;  sf=Tiff.SampleFormat.UInt;
    case 'uint16', bps=16; sf=Tiff.SampleFormat.UInt;
    case 'int16',  bps=16; sf=Tiff.SampleFormat.Int;
    case 'single', bps=32; sf=Tiff.SampleFormat.IEEEFP;
    case 'double', bps=64; sf=Tiff.SampleFormat.IEEEFP;
    otherwise, error('Unsupported class: %s',cls);
end
end

function v=ternary(cond,a,b), if cond, v=a; else, v=b; end, end
function safeClose(t), try t.close(); catch, end, end
