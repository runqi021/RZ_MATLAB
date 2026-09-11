function stitch_apply_coords_gamma(datasetRoot, ch, gammaVal)
%STITCH_APPLY_COORDS_GAMMA  Gamma-corrected avg-proj mosaic, reusing saved coords.
%   Same placement as stitch_manual_gui / stitch_apply_coords_chan, but each
%   Z page is GAMMA-corrected BEFORE the average projection:
%       norm = page / globalMax            (shared scale across ALL tiles)
%       page_g = norm.^gammaVal * globalMax
%   globalMax is the max over every page of every tile (one pre-pass), so the
%   gamma mapping is identical for all tiles and seams stay matched.
%   gammaVal<1 brightens dim somata. Output keeps the input class & intensity
%   scale. Writes stitched_ch<ch>_avgproj_gamma<g>.tif (linear mosaic untouched).
%
%   stitch_apply_coords_gamma                       % defaults: 260522 root, ch1, g=0.6
%   stitch_apply_coords_gamma(root, 1, 0.6)
%
%   Standalone. MATLAB R2021a+. Runqi Zhang / 2026.

if nargin<1 || isempty(datasetRoot)
    datasetRoot = "C:\Users\Admin\Desktop\Vglut2";
end
if nargin<2 || isempty(ch),       ch = 1;        end
if nargin<3 || isempty(gammaVal), gammaVal = 0.6; end
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

% --- Map saved file list (whatever channel was aligned) to requested channel ---
srcFiles = string(p.files(:));
nTiles   = numel(srcFiles);
dstTok   = sprintf('ch%d', ch);
files    = strings(nTiles,1);
for t=1:nTiles
    f = char(srcFiles(t));
    g = regexprep(f, ['([\\/])ch\d+([\\/])'], ['$1' dstTok '$2']); % folder token
    g = regexprep(g, '_ch\d+_', ['_' dstTok '_']);                 % filename token
    assert(isfile(g), 'Tile not found:\n  %s', g);
    files(t) = string(g);
end

inClass = class(imread(files(1),1));

% --- Pass 1: global max over all pages of all tiles ---
fprintf('Pass 1/2: scanning global max (%d tiles)...\n', nTiles);
globalMax = 0;
for t=1:nTiles
    globalMax = max(globalMax, tile_stack_max(files(t)));
end
assert(globalMax>0, 'All-zero data — cannot normalize.');
fprintf('  globalMax = %g\n', globalMax);

% --- Pass 2: gamma per page -> avg proj -> blend at saved positions ---
fprintf('Pass 2/2: gamma %.3g avg-proj + blend (%dx%d)...\n', gammaVal, outH, outW);
wT  = linear_blend_weights(H0,W0);
acc = zeros(outH,outW,'double');
wac = zeros(outH,outW,'double');
for t=1:nTiles
    img = avg_proj_gamma(files(t), H0, W0, globalMax, gammaVal); % gamma-then-mean
    r0=round(y(t)); c0=round(x(t));
    r0=max(0,min(r0,outH-H0)); c0=max(0,min(c0,outW-W0));
    rr=r0+1:r0+H0; cc=c0+1:c0+W0;
    acc(rr,cc)=acc(rr,cc)+img.*wT;
    wac(rr,cc)=wac(rr,cc)+wT;
    if mod(t,10)==0, fprintf('  %d/%d\n', t, nTiles); end
end
M = acc./max(wac,eps);

gstr = strrep(sprintf('%.2g',gammaVal),'.','p');
out  = fullfile(coordDir, sprintf('stitched_ch%d_avgproj_gamma%s.tif', ch, gstr));
write_single_tiff(out, cast(round(M),inClass), inClass);
fprintf('Saved: %s\n', out);
end

% ===================== helpers =====================
function m = tile_stack_max(file)
tf=Tiff(char(file),'r'); cln=onCleanup(@() safeClose(tf)); %#ok<NASGU>
m=0; m=max(m,double(max(tf.read(),[],'all')));
while ~tf.lastDirectory(), tf.nextDirectory(); m=max(m,double(max(tf.read(),[],'all'))); end
end

function img = avg_proj_gamma(file, H0, W0, gmax, gammaVal)
tf=Tiff(char(file),'r'); cln=onCleanup(@() safeClose(tf)); %#ok<NASGU>
acc=zeros(H0,W0,'double'); nZ=0;
done=false;
while ~done
    pg = double(tf.read())/gmax;          % normalize to [0,1] on shared scale
    acc = acc + (pg.^gammaVal)*gmax;      % gamma, then back to intensity units
    nZ = nZ+1;
    if tf.lastDirectory(), done=true; else, tf.nextDirectory(); end
end
img = acc/nZ;
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
tag.Software='stitch_apply_coords_gamma avgproj';
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
