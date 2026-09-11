function stitch_apply_coords_chan(datasetRoot, srcCh, dstCh)
%STITCH_APPLY_COORDS_CHAN  Re-render a stitched mosaic for another channel
%   using positions already solved by stitch_manual_gui (single-channel, ch1).
%
%   stitch_manual_gui aligns ONE channel and writes matlab_stitch/
%   stitch_manual_coords.mat (x,y per tile, outH/outW, H0/W0, ch1 file list).
%   This applies those exact pixel positions to a different channel's tiles,
%   matched 1:1 by swapping the channel token in each filename, and writes
%   stitched_ch<dstCh>_avgproj.tif. Blend logic mirrors the GUI's blendFull().
%
%   stitch_apply_coords_chan                       % defaults: 260522 root, ch1 -> ch3
%   stitch_apply_coords_chan(root, 1, 3)           % explicit
%
%   Standalone. MATLAB R2021a+. Runqi Zhang / 2026.

if nargin<1 || isempty(datasetRoot)
    datasetRoot = "C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\map";
end
if nargin<2 || isempty(srcCh), srcCh = 1; end
if nargin<3 || isempty(dstCh), dstCh = 3; end
datasetRoot = char(datasetRoot);

coordDir = fullfile(datasetRoot,'matlab_stitch');
coordMat = fullfile(coordDir,'stitch_manual_coords.mat');
assert(isfile(coordMat), 'Coords not found: %s', coordMat);

p = load(coordMat);
need = {'x','y','outH','outW','H0','W0','files'};
for k=1:numel(need)
    assert(isfield(p,need{k}), 'Coords missing field "%s" — re-save from the GUI.', need{k});
end
x=p.x(:); y=p.y(:); outH=p.outH; outW=p.outW; H0=p.H0; W0=p.W0;
srcFiles = string(p.files(:));
nTiles = numel(srcFiles);

% --- Map source-channel tiles to destination channel by filename token ---
srcTok = sprintf('ch%d', srcCh);
dstTok = sprintf('ch%d', dstCh);
dstFiles = strings(nTiles,1);
for t=1:nTiles
    f = char(srcFiles(t));
    g = strrep(f, [filesep srcTok filesep], [filesep dstTok filesep]); % ch1\ -> ch3\ folder
    g = regexprep(g, ['_' srcTok '_'], ['_' dstTok '_']);              % _ch1_ -> _ch3_ name
    assert(isfile(g), 'Destination tile not found:\n  %s\n(from %s)', g, f);
    dstFiles(t) = string(g);
end

% --- Probe class from first dest tile ---
I0 = imread(dstFiles(1),1);
inClass = class(I0);

% --- Build avg projections (mean over all Z pages) and blend ---
fprintf('Rendering ch%d mosaic (%d tiles, %dx%d)...\n', dstCh, nTiles, outH, outW);
wT  = linear_blend_weights(H0,W0);
acc = zeros(outH,outW,'double');
wac = zeros(outH,outW,'double');
for t=1:nTiles
    img = avg_proj(dstFiles(t), H0, W0);
    r0=round(y(t)); c0=round(x(t));
    r0=max(0,min(r0,outH-H0)); c0=max(0,min(c0,outW-W0));
    rr=r0+1:r0+H0; cc=c0+1:c0+W0;
    acc(rr,cc)=acc(rr,cc)+img.*wT;
    wac(rr,cc)=wac(rr,cc)+wT;
    if mod(t,10)==0, fprintf('  %d/%d\n', t, nTiles); end
end
M = acc./max(wac,eps);

out = fullfile(coordDir, sprintf('stitched_ch%d_avgproj.tif', dstCh));
write_single_tiff(out, cast(round(M),inClass), inClass);
fprintf('Saved: %s\n', out);
end

% ===================== helpers (copied from stitch_manual_gui) =====================
function img = avg_proj(file, H0, W0)
tf=Tiff(char(file),'r'); cln=onCleanup(@() safeClose(tf)); %#ok<NASGU>
nZ=1; while ~tf.lastDirectory(), tf.nextDirectory(); nZ=nZ+1; end
acc=zeros(H0,W0,'double');
for k=1:nZ, tf.setDirectory(k); acc=acc+double(tf.read()); end
img=acc/nZ;
end

function w=linear_blend_weights(H,W)
[xg,yg]=meshgrid(1:W,1:H);
w=double(min(min(xg-1,W-xg),min(yg-1,H-yg))+1);
end

function write_single_tiff(outFile,img,cls)
[bps,sf]=class_to_tiff_format(cls);
if isfile(outFile), delete(outFile); end
t=Tiff(outFile,'w'); cln=onCleanup(@() safeClose(t)); %#ok<NASGU>
tag.ImageLength=size(img,1); tag.ImageWidth=size(img,2);
tag.Photometric=Tiff.Photometric.MinIsBlack; tag.SamplesPerPixel=1;
tag.BitsPerSample=bps; tag.SampleFormat=sf;
tag.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
tag.Compression=Tiff.Compression.None; tag.RowsPerStrip=64;
tag.Software='stitch_apply_coords_chan avgproj';
t.setTag(tag); t.write(img);
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

function safeClose(t), try t.close(); catch, end, end
