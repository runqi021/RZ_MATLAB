function meta = struct_stack_prep_260910(tifPath, varargin)
%STRUCT_STACK_PREP_260910  Structural 3-D stack -> group average -> 1 um isotropic.
%
%   meta = struct_stack_prep_260910(tifPath)
%   meta = struct_stack_prep_260910(tifPath, 'Name', value, ...)
%
% Stage 1 of the structural-map pipeline.  Reads a raw ScanImage slow z-stack
% ONCE, averages the frames acquired at each z, resamples the volume so one
% voxel is one cubic micron, and writes one TIFF per saved channel.
%
% EVERYTHING COMES FROM THE HEADER.  Nothing about geometry is assumed:
% pixel size, slice count, frames per slice, z step and the channel list are all
% read from the SI metadata, so a single-channel stack and a two-channel one go
% down the same path and neither needs an edit.  The three numbers that matter:
%
%   pixel size  = objectiveResolution * 2*|imagingFovDeg| / pixelsPerLine
%                 (equivalently PIX_BASE/zoom; the FOV form is used because it
%                 survives scanAngleMultiplier != 1, which the zoom form does not)
%   dz          = median(diff(SI.hStackManager.zs))
%   page order  = ((slice-1)*framesPerSlice + frame-1)*nCh + chIdx
%
% THE PAGE ORDER IS THE TRAP.  Channels are interleaved PER FRAME, not per
% slice, so a stride that assumes per-slice interleaving still returns an image
% of the right size and never errors -- it just silently averages the other
% label.  chIdx is the POSITION in channelSave, while the name of the channel is
% its SI channel NUMBER (here [1 3], so channel 3 is the second plane).  Outputs
% are named by the SI number.
%
% WHY 1 um ISOTROPIC.  The raw voxel is 0.889 x 0.889 x 1 um, so any 3-D view
% of the raw array is stretched ~12% in z and every distance read off a side
% view is wrong by that much.  Lateral resampling is bicubic; the axial axis is
% only touched (pchip) when the step is not already 1 um.
%
% Averaging is done in SINGLE.  30 frames of int16 sum to at most ~1e6, well
% inside single's 2^24 exactly-representable integer range, so the mean is exact
% and the accumulator is half the size of a double one.
%
% OUTPUT  (beside the input TIFF)
%   <stem>_ch<N>_avg1um.tif    int16, [ny nx nz], 1 um isotropic, one per channel
%   <stem>_stackmeta.mat       meta struct below + per-plane percentiles
%
% int16 rather than uint16 because ScanImage digitises to int16 and the dark
% level sits near zero: clamping negatives away would put a floor on the noise
% and bias every background estimate made downstream off this file.
%
% OPTIONS
%   'PixTarget'   1.0    output voxel, um
%   'ZSurface'    0      motor z of the pia, um.  depth = ZSurface - z.
%                        Stacks are routinely started above the pia, so the
%                        shallowest planes get NEGATIVE depth and that is
%                        correct, not a bug.  NOTHING is shifted by this beyond
%                        the depth axis written into the meta file.
%   'Channels'    []     SI channel numbers to keep ([] = all saved channels)
%   'Overwrite'   false  re-read even if the outputs already exist
%
% Runqi Zhang / 2026-09-10

%% ------------------------------- options ----------------------------------
p = inputParser;
p.addParameter('PixTarget', 1.0,   @(x)isscalar(x) && x > 0);
p.addParameter('ZSurface',  0,     @isscalar);
p.addParameter('Channels',  [],    @isnumeric);
p.addParameter('Overwrite', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);

assert(isfile(tifPath), 'Not a file: %s', tifPath);
[dataDir, stem] = fileparts(tifPath);

fprintf('\n=================== struct_stack_prep_260910 ===================\n');
fprintf('file : %s\n', tifPath);

%% ---------------------------- read the header ------------------------------
meta = read_si_header(tifPath);
meta.pixTarget = opt.PixTarget;
meta.zSurface  = opt.ZSurface;

fprintf('\n--- ScanImage header ---\n');
fprintf('  SI version        %s.%s\n', meta.siVersionMajor, meta.siVersionMinor);
fprintf('  frame             %d x %d px, zoom %g\n', meta.linesPerFrame, meta.pixelsPerLine, meta.zoom);
fprintf('  pixel size        %.4f um/px  (objectiveResolution %g um/deg)\n', meta.pxUm, meta.objectiveResolution);
fprintf('  slices            %d,  framesPerSlice %d\n', meta.nSlices, meta.framesPerSlice);
fprintf('  z                 %g -> %g um motor, step %g um\n', meta.zs(1), meta.zs(end), meta.dz);
fprintf('  channels saved    %s   (%d channel(s))\n', mat2str(meta.channelSave(:).'), meta.nCh);
fprintf('  power             pzAdjust %d, %g -> %g %%, Lz %g, zPowerReference %g\n', ...
        meta.pzAdjust, meta.stackStartPower, meta.stackEndPower, meta.lengthConstants, meta.zPowerReference);
fprintf('  motorPosition     %s\n', mat2str(meta.motorPosition));
fprintf('  epoch             %s\n', meta.epochStr);
fprintf('  expected pages    %d  (%d slices x %d frames x %d ch)\n', ...
        meta.nPages, meta.nSlices, meta.framesPerSlice, meta.nCh);
fprintf('  file size         %.2f GiB  (pixel data %.2f GiB + %.0f MB of IFDs, %.0f B/page)\n', ...
        meta.fileBytes/2^30, meta.pixelBytes/2^30, meta.overheadBytes/2^20, meta.overheadPerPage);
% Compare against the EXPECTED byte count, not bytes/frameBytes: every page
% carries its own IFD and tag block (~2 kB here), so the naive division invents
% pages that are not there and cries truncation on a perfectly complete file.
if meta.fileBytes < meta.pixelBytes
    warning(['file is %.2f GiB but %d pages need %.2f GiB of pixel data -- ' ...
             'the stack is TRUNCATED; only what is present will be read.'], ...
             meta.fileBytes/2^30, meta.nPages, meta.pixelBytes/2^30);
elseif meta.overheadPerPage > 8192
    warning(['%.0f bytes/page of non-pixel data is more than an IFD should need -- ' ...
             'check the header page count (%d) against the file.'], ...
             meta.overheadPerPage, meta.nPages);
end

sib = dir(fullfile(dataDir, [regexprep(stem,'_\d+$','') '_*.tif']));
if numel(sib) > 1
    warning('%d files share this prefix; only %s is read.', numel(sib), [stem '.tif']);
end

%% -------------------------- which channels to do ---------------------------
chWanted = opt.Channels;
if isempty(chWanted), chWanted = meta.channelSave(:).'; end
chIdx = arrayfun(@(c) find(meta.channelSave(:).' == c, 1), chWanted, 'uni', 0);
assert(all(~cellfun(@isempty, chIdx)), 'Requested channel not in channelSave %s', mat2str(meta.channelSave(:).'));
chIdx = cell2mat(chIdx);

outTif = arrayfun(@(c) fullfile(dataDir, sprintf('%s_ch%d_avg1um.tif', stem, c)), ...
                  chWanted, 'uni', 0);
outMat = fullfile(dataDir, sprintf('%s_stackmeta.mat', stem));

if ~opt.Overwrite && all(cellfun(@isfile, outTif)) && isfile(outMat)
    fprintf('\noutputs already exist -- nothing to do (pass ''Overwrite'',true to force)\n');
    S = load(outMat, 'meta'); meta = S.meta; return
end

%% ------------------------------ read + average -----------------------------
% One pass over the file.  Pages are walked with nextDirectory: setDirectory(t,k)
% re-seeks from the start of the IFD chain every call and is quadratic over
% tens of thousands of pages (minutes -> hours on a stack this size).
fprintf('\n--- reading %d pages ---\n', meta.nPages);
t0 = tic;
ny = meta.linesPerFrame; nx = meta.pixelsPerLine;
nz = meta.nSlices; nf = meta.framesPerSlice; nc = meta.nCh;

ACC = zeros(ny, nx, nz, numel(chIdx), 'single');
keep = false(1, nc); keep(chIdx) = true;
slot = zeros(1, nc); slot(chIdx) = 1:numel(chIdx);

t = Tiff(tifPath, 'r');
cleanupT = onCleanup(@() close(t));
nRead = 0; done = false;
for k = 1:nz
    for fme = 1:nf
        for c = 1:nc
            if keep(c)
                ACC(:,:,k,slot(c)) = ACC(:,:,k,slot(c)) + single(t.read());
            end
            nRead = nRead + 1;
            if nRead == meta.nPages, done = true; break; end
            try
                t.nextDirectory();
            catch
                done = true;
                warning('file ended after %d of %d pages (slice %d of %d)', ...
                        nRead, meta.nPages, k, nz);
                break
            end
        end
        if done, break; end
    end
    if done, break; end
    if mod(k, 50) == 0
        fprintf('  slice %4d / %d   (%.0f s)\n', k, nz, toc(t0));
    end
end
clear cleanupT
ACC = ACC / nf;
nzGot = min(nz, ceil(nRead / (nf*nc)));
if nzGot < nz
    ACC = ACC(:,:,1:nzGot,:);
end
fprintf('  read %d pages in %.0f s -> %d slices averaged (%d frames each)\n', ...
        nRead, toc(t0), nzGot, nf);

zs    = meta.zs(1:nzGot);
depth = meta.zSurface - zs;             % um below pia; negative above it

%% ------------------------- resample to 1 um isotropic ----------------------
pt = opt.PixTarget;
nxNew = max(1, round(nx * meta.pxUm / pt));
nyNew = max(1, round(ny * meta.pxUm / pt));
fprintf('\n--- resampling to %g um isotropic ---\n', pt);
fprintf('  lateral  %.4f -> %g um/px : %dx%d -> %dx%d (bicubic)\n', meta.pxUm, pt, ny, nx, nyNew, nxNew);

dNew = depth;
needZ = abs(abs(meta.dz) - pt) > 1e-6;
if needZ
    dNew = (depth(1) : pt : depth(end)).';
    fprintf('  axial    %.4f -> %g um    : %d -> %d planes (pchip)\n', abs(meta.dz), pt, nzGot, numel(dNew));
else
    fprintf('  axial    already %g um    : %d planes unchanged\n', pt, nzGot);
end

V = cell(numel(chIdx), 1);
for ci = 1:numel(chIdx)
    A = imresize(ACC(:,:,:,ci), [nyNew nxNew], 'bicubic');   % acts on the 1st two dims
    if needZ
        sz = size(A);
        A  = reshape(A, [], sz(3)).';
        A  = interp1(depth, A, dNew, 'pchip');               % keep single: double ~4x RAM
        A  = reshape(A.', [sz(1) sz(2) numel(dNew)]);
    end
    V{ci} = A;
    fprintf('  channel %d -> [%d %d %d]\n', chWanted(ci), size(A,1), size(A,2), size(A,3));
end
clear ACC A

%% --------------------------- per-plane statistics --------------------------
% Cheap, and it is what locates the tissue surface later without another read.
% top/bot 5% of each plane, exactly as in the falloff pipeline, so a number here
% means the same thing as a number there.
fprintf('\n--- per-plane top/bot 5%% ---\n');
stats = struct('ch', {}, 'top5', {}, 'bot5', {}, 'med', {});
for ci = 1:numel(chIdx)
    A = V{ci}; nzz = size(A,3);
    top5 = zeros(nzz,1); bot5 = zeros(nzz,1); med = zeros(nzz,1);
    for k = 1:nzz
        v = sort(reshape(A(:,:,k), [], 1));
        m = max(1, round(0.05*numel(v)));
        bot5(k) = mean(v(1:m));  top5(k) = mean(v(end-m+1:end));  med(k) = v(round(end/2));
    end
    stats(ci).ch = chWanted(ci); stats(ci).top5 = top5; stats(ci).bot5 = bot5; stats(ci).med = med;

    % Surface marker, REPORTED ONLY -- nothing is shifted by it.  On a
    % power-ramped stack the raw (top-bot) peak is meaningless (power climbs
    % with depth by construction), so the half-rise of the normalised curve is
    % used: the blurred water/tissue edge where the beam enters.
    d  = top5 - bot5;
    dS = smoothdata(d, 'movmean', max(3, round(7/pt)));
    dN = (dS - min(dS)) / max(eps, max(dS) - min(dS));
    iH = find(dN >= 0.5, 1, 'first');
    stats(ci).surfHalfRise = dNew(max(1, iH));
    [~, iP] = max(dS);
    stats(ci).surfPeak = dNew(iP);
    fprintf('  ch%d  (top-bot) half-rise at depth %+.1f um, peak at %+.1f um  [reported, NOT applied]\n', ...
            chWanted(ci), stats(ci).surfHalfRise, stats(ci).surfPeak);
end

%% -------------------------------- save -------------------------------------
fprintf('\n--- saving ---\n');
for ci = 1:numel(chIdx)
    write_tiff_int16(outTif{ci}, V{ci});
    d = dir(outTif{ci});
    fprintf('  %s  (%.2f GB)\n', outTif{ci}, d.bytes/2^30);
end

meta.depth     = dNew;
meta.zsUsed    = zs;
meta.nSlicesRead = nzGot;
meta.pxOut     = pt;
meta.size      = [nyNew nxNew numel(dNew)];
meta.channels  = chWanted;
meta.outTif    = outTif;
meta.planeStats = stats;
meta.builtOn   = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
save(outMat, 'meta', '-v7.3');
fprintf('  %s\n', outMat);
fprintf('\ndone in %.0f s\n', toc(t0));
end

% =============================================================================
function meta = read_si_header(tifPath)
%READ_SI_HEADER  Everything the pipeline needs, straight out of the TIFF.

t = Tiff(tifPath, 'r');
soft = t.getTag('Software');
try, desc = t.getTag('ImageDescription'); catch, desc = ''; end
t.close();

s = string(soft);
s = replace(s, char(8593), newline);
s = replace(s, "↵", newline);

gn = @(key) local_num(s, key);
gs = @(key) local_str(s, key);

meta.tifPath          = tifPath;
meta.siVersionMajor   = char(strip(gs("SI.VERSION_MAJOR"), "'"));
meta.siVersionMinor   = char(strip(gs("SI.VERSION_MINOR"), "'"));
meta.linesPerFrame    = gn("SI.hRoiManager.linesPerFrame");
meta.pixelsPerLine    = gn("SI.hRoiManager.pixelsPerLine");
meta.zoom             = gn("SI.hRoiManager.scanZoomFactor");
meta.objectiveResolution = gn("SI.objectiveResolution");
meta.framesPerSlice   = gn("SI.hStackManager.framesPerSlice");
meta.nSlices          = gn("SI.hStackManager.numSlices");
meta.stackZStepSize   = gn("SI.hStackManager.stackZStepSize");
meta.stackZStartPos   = gn("SI.hStackManager.stackZStartPos");
meta.stackZEndPos     = gn("SI.hStackManager.stackZEndPos");
meta.zPowerReference  = gn("SI.hStackManager.zPowerReference");
meta.pzAdjust         = gn("SI.hBeams.pzAdjust");
meta.stackStartPower  = gn("SI.hBeams.stackStartPower");
meta.stackEndPower    = gn("SI.hBeams.stackEndPower");
meta.lengthConstants  = gn("SI.hBeams.lengthConstants");
meta.scanFrameRate    = gn("SI.hRoiManager.scanFrameRate");
meta.channelSave      = local_vec(s, "SI.hChannels.channelSave");
meta.motorPosition    = local_vec(s, "SI.hMotors.motorPosition");
meta.zs               = local_vec(s, "SI.hStackManager.zs");
meta.softwareRaw      = char(s);

if isempty(meta.channelSave), meta.channelSave = 1; end
meta.nCh = numel(meta.channelSave);

% Pixel size from the FOV, not from zoom.  imagingFovUm already carries
% scanAngleMultiplier and any scan rotation, so it stays right when a stack is
% acquired with a non-square or shifted field, where PIX_BASE/zoom quietly is not.
fovUm = local_mat(s, "SI.hRoiManager.imagingFovUm");
if ~isempty(fovUm) && size(fovUm,2) >= 2
    widthUm = max(fovUm(:,1)) - min(fovUm(:,1));
else
    fovDeg  = local_mat(s, "SI.hRoiManager.imagingFovDeg");
    widthUm = (max(fovDeg(:,1)) - min(fovDeg(:,1))) * meta.objectiveResolution;
end
meta.fovUm = widthUm;
meta.pxUm  = widthUm / meta.pixelsPerLine;

if isempty(meta.zs)
    meta.zs = (meta.stackZStartPos : meta.stackZStepSize : meta.stackZEndPos).';
end
meta.zs = meta.zs(:);
if numel(meta.zs) > meta.nSlices, meta.zs = meta.zs(1:meta.nSlices); end
meta.dz = median(diff(meta.zs));

meta.nPages    = meta.nSlices * meta.framesPerSlice * meta.nCh;
d = dir(tifPath); meta.fileBytes = d.bytes;
meta.pixelBytes     = meta.nPages * meta.linesPerFrame * meta.pixelsPerLine * 2;
meta.overheadBytes  = meta.fileBytes - meta.pixelBytes;
meta.overheadPerPage = meta.overheadBytes / max(1, meta.nPages);

ep = local_vec(string(desc), "epoch");
if numel(ep) >= 6
    meta.epochStr = sprintf('%04d-%02d-%02d %02d:%02d:%05.2f', ep(1),ep(2),ep(3),ep(4),ep(5),ep(6));
else
    meta.epochStr = '(none)';
end
end

function v = local_num(s, key)
raw = local_str(s, key);
if strcmpi(raw, "true"),  v = 1; return; end     % SI writes booleans as words,
if strcmpi(raw, "false"), v = 0; return; end     % which str2double turns to NaN
v = local_vec(s, key);
if isempty(v), v = NaN; else, v = v(1); end
end

function out = local_str(s, key)
tok = regexp(s, key + "\s*=\s*([^\n\r]*)", 'tokens', 'once');
if isempty(tok), out = ""; else, out = strip(string(tok{1})); end
end

function v = local_vec(s, key)
raw = local_str(s, key);
v = [];
if raw == "", return; end
raw = erase(raw, ["[", "]", "{", "}", "'"]);
raw = replace(raw, ";", " ");
raw = replace(raw, ",", " ");
v = str2double(split(strtrim(raw)));
v = v(~isnan(v));
end

function M = local_mat(s, key)
raw = local_str(s, key);
M = [];
if raw == "", return; end
raw = erase(raw, ["[", "]"]);
rows = split(raw, ";");
M = zeros(numel(rows), 0);
for i = 1:numel(rows)
    v = str2double(split(strtrim(replace(rows(i), ",", " "))));
    v = v(~isnan(v)).';
    if i == 1, M = zeros(numel(rows), numel(v)); end
    M(i, 1:numel(v)) = v; %#ok<AGROW>
end
end

% =============================================================================
function write_tiff_int16(fn, V)
%WRITE_TIFF_INT16  Multi-page int16 TIFF, one page per z plane.
%
% int16 preserves ScanImage's signed digitiser scale, so the dark level stays
% where it was and a background estimate made off this file is not floored at 0.

V = int16(round(max(min(V, 32767), -32768)));
[ny, nx, nz] = size(V);

if isfile(fn), delete(fn); end
t = Tiff(fn, 'w8');            % BigTIFF: these volumes go past the 4 GB mark
c = onCleanup(@() close(t));
tag.ImageLength         = ny;
tag.ImageWidth          = nx;
tag.Photometric         = Tiff.Photometric.MinIsBlack;
tag.BitsPerSample       = 16;
tag.SampleFormat        = Tiff.SampleFormat.Int;
tag.SamplesPerPixel     = 1;
tag.RowsPerStrip        = ny;   % ONE strip per page.  MATLAB's Tiff writer pads
                                % the final strip out to a full RowsPerStrip
                                % instead of writing a short one, so a 455-row
                                % page at 64 rows/strip carries 57 rows of zeros
                                % -- measured at 52 237 B/page, 11%% of the file.
                                % A single strip has no partial strip to pad.
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression         = Tiff.Compression.None;
tag.Software            = 'struct_stack_prep_260910';

for k = 1:nz
    t.setTag(tag);
    t.write(V(:,:,k));
    if k < nz, t.writeDirectory(); end
end
end
