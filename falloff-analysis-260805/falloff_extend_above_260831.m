function falloff_extend_above_260831(dataDir, roiTag, depthMin, channel)
%FALLOFF_EXTEND_ABOVE_260831  Add the above-pia planes to a cached falloff volume.
%
% falloff_projections_260805.m builds its grid with
%     depth = (0 : dzNat : max(...))
% so every plane ABOVE the surface is discarded and the cached volume starts at
% the pia. This prepends those planes back, writing a SEPARATE cache so the
% original is untouched:
%     falloff_proj_<tag>_vol.mat  ->  falloff_proj_<tag>_vol_from<D>.mat
%
% IT DOES NOT RE-READ THE WHOLE ROI. Above-surface depths exist only in the
% SHALLOWEST stack, so only that file is opened, and only the slices covering
% [depthMin, -1]. That is ~10 slices instead of ~13 GB across four stacks.
%
% HOW FAR UP EACH ROI CAN GO is set by acquisition, not by this script:
%     260728_vglut2 roi1_z-10-100   ->  -10 um
%     260804_shiver roi2_2x_z-50-150 -> -50 um
% asking for more than the stack holds is an error rather than a silent clamp.
%
% THE PROCESSING MATCHES build_volume EXACTLY, or the prepended planes would sit
% on a different scale from the rest of the volume:
%   * average framesPerSlice frames, taking the saved channel's interleaved page
%   * bicubic resample to PIX_TARGET (1 um), the same imresize call
%   * pctPlane per plane from the SI power ramp (constant here: pzAdjust false)
%   * pages walked with nextDirectory, never setDirectory(t,k) -- the latter is
%     quadratic over these 23k-IFD files
%
% Runqi Zhang / 2026-08-31

%% ===================== USER-EDITABLE =====================
DATA_DIR  = 'C:\fall-off\260728_vglut2_vessel';
ROI_TAG   = 'roi1';
DEPTH_MIN = -10;           % um; negative = above the pia
CHANNEL   = [];
PIX_BASE  = 1.7778;
PIX_TARGET= 1.0;
%% =========================================================
if nargin >= 1 && ~isempty(dataDir), DATA_DIR  = dataDir; end
if nargin >= 2 && ~isempty(roiTag),  ROI_TAG   = roiTag;  end
if nargin >= 3 && ~isempty(depthMin),DEPTH_MIN = depthMin;end
if nargin >= 4,                      CHANNEL   = channel; end
assert(DEPTH_MIN < 0, 'DEPTH_MIN must be negative (above the surface)');

volFile = fullfile(DATA_DIR, sprintf('falloff_proj_%s_vol.mat', ROI_TAG));
assert(isfile(volFile), 'no base cache: %s', volFile);
outFile = fullfile(DATA_DIR, sprintf('falloff_proj_%s_vol_from%d.mat', ROI_TAG, DEPTH_MIN));

%% --- find the shallowest stack ------------------------------------------------
f = dir(fullfile(DATA_DIR, sprintf('%s*_z*-*_*.tif', ROI_TAG)));
assert(~isempty(f), 'no stacks for %s in %s', ROI_TAG, DATA_DIR);
dSt = nan(numel(f),1);
for i = 1:numel(f)
    t = regexp(f(i).name, '_z(-?\d+)-(-?\d+)_', 'tokens', 'once');
    if ~isempty(t), dSt(i) = str2double(t{1}); end
end
[~, iSh] = min(dSt);
tp = fullfile(f(iSh).folder, f(iSh).name);
fprintf('shallowest stack: %s  (starts at %g um)\n', f(iSh).name, dSt(iSh));
assert(dSt(iSh) <= DEPTH_MIN, ...
    ['%s only reaches %g um above the pia; %g um was asked for. ' ...
     'That depth was never acquired.'], f(iSh).name, dSt(iSh), DEPTH_MIN);

%% --- metadata from the header (fast; no imfinfo over 23k IFDs) ---------------
t = Tiff(tp,'r');
sw=''; try, sw=t.getTag('Software'); catch, end
im=''; try, im=t.getTag('ImageDescription'); catch, end
H = t.getTag('ImageLength'); W = t.getTag('ImageWidth');
txt = [sw newline im];
gv = @(k) str2double(regexp(txt,[k '\s*=\s*([-\d\.eE]+)'],'tokens','once'));
n    = gv('SI\.hStackManager\.numSlices');
fpsl = gv('SI\.hStackManager\.framesPerSlice');
dz   = abs(gv('SI\.hStackManager\.stackZStepSize'));
zoom = gv('SI\.hRoiManager\.scanZoomFactor');
ppl  = gv('SI\.hRoiManager\.pixelsPerLine');
P0   = gv('SI\.hBeams\.stackStartPower');
P1   = gv('SI\.hBeams\.stackEndPower');
pzOn = ~isempty(regexpi(txt,'SI\.hBeams\.pzAdjust\s*=\s*true','once'));
chS  = str2num(char(regexp(txt,'SI\.hChannels\.channelSave\s*=\s*([^\n]*)','tokens','once'))); %#ok<ST2NM>
nCh  = max(numel(chS),1);
chIdx = 1;
if nCh > 1
    assert(~isempty(CHANNEL), ['%s saves %d channels -- pass CHANNEL, or the ' ...
        'wrong label is read without an error'], f(iSh).name, nCh);
    chIdx = find(chS == CHANNEL, 1);
    assert(~isempty(chIdx), 'channel %g not in channelSave %s', CHANNEL, mat2str(chS));
end
assert(~pzOn || abs(P1-P0) < 1e-9, ...
    'ramped power on the shallow stack -- pctZ would need the exp ramp, not implemented here');
if ~isfinite(ppl) || ppl <= 0, ppl = W; end
pxUm = PIX_BASE * (512/ppl) / zoom;
d    = dSt(iSh) + (0:n-1)*dz;               % depth of each slice, 1 um steps
fprintf('  %d slices, %d frames/slice, %d ch, %g um steps, %.4f um/px, power %g%%\n', ...
        n, fpsl, nCh, dz, pxUm, P0);

want = find(d >= DEPTH_MIN & d < 0);
assert(~isempty(want), 'no slices between %g and 0', DEPTH_MIN);
fprintf('  slices %d-%d cover depth %g..%g um\n', want(1), want(end), d(want(1)), d(want(end)));

%% --- read those slices --------------------------------------------------------
sc    = pxUm / PIX_TARGET;
probe = imresize(zeros(H,W,'single'), sc, 'bicubic');
Vnew  = zeros(size(probe,1), size(probe,2), numel(want), 'single');

% Sequential walk: setDirectory(t,k) is quadratic over these files.
pageWanted = zeros(numel(want)*fpsl, 1);
c = 0;
for q = 1:numel(want)
    iz = want(q);
    for fr = 1:fpsl
        c = c + 1;
        pageWanted(c) = ((iz-1)*fpsl + fr-1)*nCh + chIdx;   % 1-based page
    end
end
maxPage = max(pageWanted);
acc = zeros(H, W, numel(want), 'single');
setDirectory(t, 1);
pg = 1;
tic;
while true
    hit = find(pageWanted == pg, 1);
    if ~isempty(hit)
        q = ceil(hit / fpsl);
        acc(:,:,q) = acc(:,:,q) + single(t.read());
    end
    if pg >= maxPage, break; end
    nextDirectory(t); pg = pg + 1;
end
t.close();
for q = 1:numel(want)
    Vnew(:,:,q) = imresize(acc(:,:,q)/fpsl, sc, 'bicubic');
end
fprintf('  read %d pages in %.1f s\n', maxPage, toc);

%% --- prepend to the cached volume ---------------------------------------------
S = load(volFile);
assert(abs(S.pxUm - pxUm) < 1e-6, 'pxUm mismatch: cache %g, stack %g', S.pxUm, pxUm);
assert(isequal(size(S.V,1), size(Vnew,1)) && isequal(size(S.V,2), size(Vnew,2)), ...
    'in-plane size mismatch: cache %s, new %s', mat2str(size(S.V,[1 2])), mat2str(size(Vnew,[1 2])));

V        = cat(3, Vnew, S.V);
depth    = [d(want).'; S.depth];
pctPlane = [repmat(P0, numel(want), 1); S.pctPlane];
DATA_DIR_ = S.DATA_DIR; ROI_TAG_ = S.ROI_TAG; %#ok<NASGU>
assert(all(diff(depth) > 0), 'depth is not strictly increasing after the prepend');

save(outFile, 'V','depth','pctPlane','pxUm','DATA_DIR_','ROI_TAG_','-v7.3');
fprintf('\nwrote %s\n  %d planes, %g..%g um (was %g..%g)\n', outFile, ...
        numel(depth), depth(1), depth(end), S.depth(1), S.depth(end));
end
