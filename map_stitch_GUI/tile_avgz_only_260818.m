clear; clc; close all;
avgz_only_all_files

function avgz_only_all_files()
% avgz_only_all_files
% ------------------------------------------------------------
% AVGZ-ONLY variant of tile_deinterleave_N_groupZavg_260131.m.
%
% For EVERY ScanImage TIFF under rawDir, writes ONLY the avgZ (mean over T)
% per channel -- it does NOT write the deinterleaved Z*T stacks. Those run
% ~1.4 GB per tile per channel and are not used by stitch_manual_gui or
% NChan_vol_stitch_zUse, which read avgz only.
%
% Output folder structure:
%   rawDir\ch<id>\avgz\<base>_ch<id>_avgz.tif      (nZ pages)
%
% Pages are streamed IN FILE ORDER (nextDirectory), not random-accessed with
% setDirectory -- for a 5400-page 2.8 GB tile that is the difference between
% a linear scan and re-walking the IFD chain 5400 times.
%
% Standalone. MATLAB R2021a+. Runqi Zhang / 2026.
% ------------------------------------------------------------

%% ========================= USER SETTINGS =========================
rawDir = "C:\260824_Vglut2-soma-g8s_vagotomized\map";

% Channel IDs from SI.hChannels.channelSave
channels = [1];

recursive    = false;   % true if tiles live in subfolders
strictMeta   = true;    % require pages == nCh*nZ*nT
rowsPerStrip = 64;
maxFiles     = Inf;     % process at most N tiles this run (Inf = all).
                        % Re-running skips finished avgz, so this is resumable.

% If raw is int16, convert to uint16 safely (no wrapping)
forceUint16 = true;
int16Offset = 32768;

rawDir   = string(rawDir);
channels = channels(:).';
assert(~isempty(channels), "channels must be non-empty");

%% ========================= DISCOVER FILES =========================
if ~isfolder(rawDir), error("Not a folder: %s", rawDir); end

if recursive
    files = list_all_tiffs_recursive(rawDir);
else
    files = list_all_tiffs_nonrecursive(rawDir);
end
files = unique(files);

% Do NOT process our own outputs
for cc = channels
    files = files(~contains(lower(files), lower(filesep + "ch" + string(cc) + filesep)));
end
files = files(~contains(lower(files), lower(filesep + "avgz" + filesep)));
files = files(~contains(lower(files), "stitched"));
files = files(~contains(lower(files), "cache"));

if isempty(files), error("No TIFFs found in %s", rawDir); end
fprintf("[scan] %d files\n", numel(files));

%% ========================= OUTPUT FOLDERS =========================
nReqCh = numel(channels);
outZ   = strings(nReqCh,1);
for k = 1:nReqCh
    outZ(k) = fullfile(rawDir, sprintf("ch%d", channels(k)), "avgz");
    if ~isfolder(outZ(k)), mkdir(outZ(k)); end
end

%% ========================= MAIN LOOP =========================
nOK = 0; nSkip = 0; tAll = tic;

for k = 1:numel(files)
    if nOK >= maxFiles
        fprintf("\n[stop] maxFiles=%d reached; re-run to continue.\n", maxFiles); break;
    end
    fn = files(k);
    [~, baseName, ext] = fileparts(fn);
    fprintf("\n[%4d/%4d] %s%s\n", k, numel(files), baseName, ext);

    % --- fast ScanImage meta from FIRST IFD only ---
    try
        meta = read_scanimage_meta_fast_firstIFD(fn);
    catch ME
        nSkip = nSkip + 1; fprintf("  [skip] meta read failed: %s\n", ME.message); continue;
    end
    if ~isKey(meta.kv,'SI.hChannels.channelSave')
        nSkip = nSkip + 1; fprintf("  [skip] no SI.hChannels.channelSave\n"); continue;
    end

    chSave = parse_numeric_vector(meta.kv('SI.hChannels.channelSave'));
    chSave = chSave(:).';
    if isempty(chSave), nSkip = nSkip + 1; fprintf("  [skip] empty channelSave\n"); continue; end
    nCh = numel(chSave);

    chIdx = nan(nReqCh,1);
    bad = false;
    for c = 1:nReqCh
        j = find(chSave==channels(c), 1);
        if isempty(j)
            fprintf("  [skip] ch%d not in channelSave=[%s]\n", channels(c), num2str(chSave));
            bad = true; break;
        end
        chIdx(c) = j;
    end
    if bad, nSkip = nSkip + 1; continue; end

    [nZ, nT, orderMode] = infer_zt_fast(meta);
    if ~isfinite(nZ) || ~isfinite(nT) || nZ<1 || nT<1
        nSkip = nSkip + 1; fprintf("  [skip] bad nZ/nT\n"); continue;
    end

    if strictMeta
        try
            nPages = count_tiff_pages_file(fn);
        catch ME
            nSkip = nSkip + 1; fprintf("  [skip] count pages failed: %s\n", ME.message); continue;
        end
        expected = nCh*nZ*nT;
        if nPages ~= expected
            error("STRICT META mismatch: pages=%d expected=%d (=nCh*nZ*nT) for %s", nPages, expected, baseName);
        end
    end

    outAz = strings(nReqCh,1);
    for c = 1:nReqCh
        outAz(c) = fullfile(outZ(c), baseName + sprintf("_ch%d_avgz.tif", channels(c)));
    end

    fprintf("  meta: channelSave=[%s]  nZ=%d  nT=%d  order=%s\n", num2str(chSave), nZ, nT, orderMode);

    % If all avgz already exist with the right page count, skip
    if all(isfile(outAz))
        try
            if all(arrayfun(@(f) count_tiff_pages_file(f)==nZ, outAz))
                fprintf("  [skip] avgz already exist\n"); continue;
            end
        catch
            % fallthrough -> rewrite
        end
    end

    try
        tT = tic;
        stream_avgz_multichan(fn, outAz, chIdx, nCh, nZ, nT, orderMode, ...
                              rowsPerStrip, forceUint16, int16Offset);
        nOK = nOK + 1;
        fprintf("  [ok] avgz written (%.1f s)\n", toc(tT));
    catch ME
        nSkip = nSkip + 1; fprintf("  [FAIL] %s\n", ME.message);
    end
end

fprintf("\n[DONE] ok=%d  skipped=%d  (%.1f min)\n", nOK, nSkip, toc(tAll)/60);
for c = 1:nReqCh
    fprintf("  ch%d avgz : %s\n", channels(c), outZ(c));
end
end

%% ========================= AVGZ WORKER (N channels, sequential pages) =========================
function stream_avgz_multichan(rawFn, outAz, chIdx, nCh, nZ, nT, orderMode, ...
                               rowsPerStrip, forceUint16, int16Offset)
% Streams every page once in file order and accumulates the mean over T for
% each requested channel, writing one avgz page per Z.

tr = Tiff(rawFn,'r'); c0 = onCleanup(@() tr.close()); %#ok<NASGU>
tr.setDirectory(1);
I0 = tr.read();
[H,W] = size(I0);
inClassRaw = class(I0);
needShift  = forceUint16 && strcmpi(inClassRaw,'int16');
if needShift, outClass = 'uint16'; else, outClass = inClassRaw; end
[bitsPerSample, sampleFormat] = class_to_tiff_format(outClass);

nReq = numel(outAz);
for c = 1:nReq
    if isfile(outAz(c)), delete(outAz(c)); end
end

tw = cell(nReq,1);
for c = 1:nReq
    tw{c} = Tiff(outAz(c),'w');
end
cl = onCleanup(@() close_all_tiffs(tw)); %#ok<NASGU>

tag = struct();
tag.ImageLength         = H;
tag.ImageWidth          = W;
tag.Photometric         = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel     = 1;
tag.BitsPerSample       = bitsPerSample;
tag.SampleFormat        = sampleFormat;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression         = Tiff.Compression.None;
tag.RowsPerStrip        = rowsPerStrip;

% Carry the ScanImage header through to the avgz copy, so downstream tools
% (stitch_manual_gui_fast_260818) can read motorPosition / imagingFovUm from
% the avgz file itself instead of having to hunt down the raw tile.
try
    siTxt = tr.getTag('Software');
    if ~isempty(siTxt), tag.Software = siTxt; end
catch
end
try
    imTxt = tr.getTag('ImageDescription');
    if ~isempty(imTxt), tag.ImageDescription = imTxt; end
catch
end

% Precompute, for every page, which requested channel it belongs to (0 = skip)
% and which Z it lands on. Cheap (nCh*nZ*nT integers) and keeps the inner
% loop free of per-page index arithmetic.
nPages  = nCh*nZ*nT;
pageReq = zeros(nPages,1);   % index into chIdx, or 0
pageZ   = zeros(nPages,1);
for c = 1:nReq
    for z = 1:nZ
        for t = 1:nT
            p = scanimage_page_index(chIdx(c), z, t, nCh, nZ, nT, orderMode);
            pageReq(p) = c;
            pageZ(p)   = z;
        end
    end
end

acc     = zeros(H,W,nReq,'double');
written = zeros(nReq,1);
curZ    = 0;

for p = 1:nPages
    if p > 1, tr.nextDirectory(); end          % sequential scan, no IFD re-walk
    c = pageReq(p);
    if c == 0, continue; end

    z = pageZ(p);
    if z ~= curZ
        if curZ > 0
            [acc, written] = flush_z(tw, tag, acc, written, nT, nZ, outClass);
        end
        curZ = z;
    end

    I = double(tr.read());
    if needShift, I = min(max(I + double(int16Offset), 0), 65535); end
    acc(:,:,c) = acc(:,:,c) + I;

    if mod(p, nCh*nT) == 0
        fprintf("  [read] z %d/%d\n", z, nZ);
    end
end
flush_z(tw, tag, acc, written, nT, nZ, outClass);   % last Z
end

function [acc, written] = flush_z(tw, tag, acc, written, nT, nZ, outClass)
for c = 1:numel(tw)
    if written(c) > 0, tw{c}.writeDirectory(); end
    tw{c}.setTag(tag);
    tw{c}.write(cast(round(acc(:,:,c) ./ max(nT,1)), outClass));
    written(c) = written(c) + 1;
    if written(c) > nZ, error("wrote more than nZ=%d avgz pages", nZ); end
end
acc(:) = 0;
end

function close_all_tiffs(tw)
for c = 1:numel(tw)
    try, tw{c}.close(); catch, end
end
end

%% ========================= UTILITIES =========================
function files = list_all_tiffs_nonrecursive(root)
cands = [dir(fullfile(root,"*.tif")); dir(fullfile(root,"*.tiff"))];
files = strings(numel(cands),1);
for k=1:numel(cands)
    files(k) = string(fullfile(cands(k).folder, cands(k).name));
end
end

function files = list_all_tiffs_recursive(root)
cands = [dir(fullfile(root,"**","*.tif")); dir(fullfile(root,"**","*.tiff"))];
files = strings(numel(cands),1);
for k=1:numel(cands)
    files(k) = string(fullfile(cands(k).folder, cands(k).name));
end
end

function meta = read_scanimage_meta_fast_firstIFD(tiffFile)
t = Tiff(tiffFile,'r');
c = onCleanup(@() t.close()); %#ok<NASGU>

txt = "";
try, txt = txt + string(t.getTag('Software')) + newline; catch, end
try, txt = txt + string(t.getTag('ImageDescription')) + newline; catch, end

txt = replace(txt, char(8629), newline);
txt = replace(txt, "\r\n", newline);
txt = replace(txt, "\r", newline);

meta = struct();
meta.kv = containers.Map('KeyType','char','ValueType','char');

lines = splitlines(txt);
for i = 1:numel(lines)
    L = strtrim(lines(i));
    if startsWith(L,"SI.")
        eq = strfind(L,'=');
        if ~isempty(eq)
            k = strtrim(extractBefore(L, eq(1)));
            v = strtrim(extractAfter(L,  eq(1)));
            meta.kv(char(k)) = char(v);
        end
    end
end
end

function v = parse_numeric_vector(s)
tok = regexp(char(s), '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
v = str2double(tok);
end

function [nZ, nT, orderMode] = infer_zt_fast(meta)
orderMode = "sliceMajor";
if isKey(meta.kv,'SI.hFastZ.enable')
    en = lower(string(meta.kv('SI.hFastZ.enable')));
    if en=="true" || en=="1"
        orderMode = "volumeInterleaved";
    end
end

nT = NaN;
if isKey(meta.kv,'SI.hStackManager.framesPerSlice')
    nT = str2double(meta.kv('SI.hStackManager.framesPerSlice'));
end

nZ = NaN;
if isKey(meta.kv,'SI.hStackManager.zs')
    zs = parse_numeric_vector(meta.kv('SI.hStackManager.zs'));
    if ~isempty(zs), nZ = numel(zs); end
elseif isKey(meta.kv,'SI.hStackManager.numSlices')
    nZ = str2double(meta.kv('SI.hStackManager.numSlices'));
end
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

function page = scanimage_page_index(cIdx, zIdx, tIdx, nCh, nZ, nT, orderMode)
if orderMode=="volumeInterleaved"
    frameLinear = (tIdx-1)*nZ + zIdx;
else
    frameLinear = (zIdx-1)*nT + tIdx;
end
page = (frameLinear-1)*nCh + cIdx;
end

function [bitsPerSample, sampleFormat] = class_to_tiff_format(cls)
switch cls
    case 'uint8',  bitsPerSample = 8;  sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint16', bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint32', bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.UInt;
    case 'int16',  bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.Int;
    case 'single', bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'double', bitsPerSample = 64; sampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise, error("Unsupported TIFF class: %s", cls);
end
end
