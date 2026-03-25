clear; clc; close all;
deinterleave_and_avgz_all_files

function deinterleave_and_avgz_all_files()
% deinterleave_and_avgz_all_files
% ------------------------------------------------------------
% For EVERY ScanImage TIFF under rawDir:
%   - Save DEINTERLEAVED per channel: Z*T pages (z-major, then t)
%   - Save AVGZ (mean over T) per channel: Z pages
%
% Output folder structure:
%   rawDir\ch<ch1_id>\                (deinterleave Z*T)
%   rawDir\ch<ch1_id>\avgz\           (avgZ Z)
%   rawDir\ch<ch2_id>\                (deinterleave Z*T)
%   rawDir\ch<ch2_id>\avgz\
%
% Reads pixels with Tiff.setDirectory() (streaming, low RAM).
% ------------------------------------------------------------

%% ========================= USER SETTINGS =========================
rawDir = "D:\RUNQI\260319_shiver_cal-590_enblock";

% Channel IDs from SI.hChannels.channelSave
%   Two-channel:  channels = [1 3]
%   Single-channel: channels = 1   (or 3, etc.)
channels = 2; %[2];

recursive   = false;   % true if your tiles are in subfolders
strictMeta  = true;    % require pages == nCh*nZ*nT
rowsPerStrip = 64;

% If raw is int16, convert to uint16 safely (no wrapping)
forceUint16 = true;
int16Offset = 32768;

nReqCh = numel(channels);
if nReqCh == 2
    ch1_id = channels(1);
    ch2_id = channels(2);
elseif nReqCh == 1
    ch1_id = channels(1);
else
    error("channels must have 1 or 2 elements");
end

%% ========================= DISCOVER FILES =========================
if ~isfolder(rawDir), error("Not a folder: %s", rawDir); end

if recursive
    files = list_all_tiffs_recursive(rawDir);
else
    files = list_all_tiffs_nonrecursive(rawDir);
end
files = unique(files);

% Do NOT process your own outputs
for cc = channels(:).'
    files = files(~contains(lower(files), lower(filesep + "ch" + string(cc) + filesep)));
end
files = files(~contains(lower(files), lower(filesep + "avgz" + filesep)));
files = files(~contains(lower(files), "stitched"));
files = files(~contains(lower(files), "cache"));

if isempty(files)
    error("No TIFFs found in %s", rawDir);
end
fprintf("[scan] %d files\n", numel(files));

%% ========================= OUTPUT FOLDERS =========================
outCh1  = fullfile(rawDir, sprintf("ch%d", ch1_id));
outCh1Z = fullfile(outCh1, "avgz");
ensure_dir(outCh1); ensure_dir(outCh1Z);
if nReqCh == 2
    outCh2  = fullfile(rawDir, sprintf("ch%d", ch2_id));
    outCh2Z = fullfile(outCh2, "avgz");
    ensure_dir(outCh2); ensure_dir(outCh2Z);
end

%% ========================= MAIN LOOP =========================
nOK = 0; nSkip = 0;

for k = 1:numel(files)
    fn = files(k);
    [~, baseName, ext] = fileparts(fn);

    fprintf("\n[%4d/%4d] %s%s\n", k, numel(files), baseName, ext);

    % --- fast ScanImage meta from FIRST IFD only ---
    try
        meta = read_scanimage_meta_fast_firstIFD(fn);
    catch ME
        nSkip = nSkip + 1;
        fprintf("  [skip] meta read failed: %s\n", ME.message);
        continue;
    end

    if ~isKey(meta.kv,'SI.hChannels.channelSave')
        nSkip = nSkip + 1;
        fprintf("  [skip] no SI.hChannels.channelSave (not ScanImage?)\n");
        continue;
    end

    chSave = parse_numeric_vector(meta.kv('SI.hChannels.channelSave'));
    chSave = chSave(:).';
    if isempty(chSave)
        nSkip = nSkip + 1;
        fprintf("  [skip] empty channelSave\n");
        continue;
    end

    nCh = numel(chSave);
    ch1_idx = find(chSave==ch1_id, 1);
    if isempty(ch1_idx)
        nSkip = nSkip + 1;
        fprintf("  [skip] ch%d not in channelSave=[%s]\n", ch1_id, num2str(chSave));
        continue;
    end
    if nReqCh == 2
        ch2_idx = find(chSave==ch2_id, 1);
        if isempty(ch2_idx)
            nSkip = nSkip + 1;
            fprintf("  [skip] ch%d not in channelSave=[%s]\n", ch2_id, num2str(chSave));
            continue;
        end
    end

    [nZ, nT, orderMode] = infer_zt_fast(meta);
    if ~isfinite(nZ) || ~isfinite(nT) || nZ<1 || nT<1
        nSkip = nSkip + 1;
        fprintf("  [skip] bad nZ/nT\n");
        continue;
    end

    % --- strict page count check ---
    if strictMeta
        try
            nPages = count_tiff_pages_file(fn);
        catch ME
            nSkip = nSkip + 1;
            fprintf("  [skip] count pages failed: %s\n", ME.message);
            continue;
        end
        expected = nCh*nZ*nT;
        if nPages ~= expected
            error("STRICT META mismatch: pages=%d expected=%d (=nCh*nZ*nT)", nPages, expected);
        end
    end

    % --- output paths ---
    needDeint = (nCh > 1);  % skip deinterleave when raw file is already single-channel
    if needDeint
        outDe1 = fullfile(outCh1,  baseName + sprintf("_ch%d_deint.tif", ch1_id));
    end
    outAz1 = fullfile(outCh1Z, baseName + sprintf("_ch%d_avgz.tif",  ch1_id));
    if nReqCh == 2
        outDe2 = fullfile(outCh2,  baseName + sprintf("_ch%d_deint.tif", ch2_id));
        outAz2 = fullfile(outCh2Z, baseName + sprintf("_ch%d_avgz.tif",  ch2_id));
    end

    fprintf("  meta: channelSave=[%s]  nZ=%d  nT=%d  order=%s\n", num2str(chSave), nZ, nT, orderMode);

    % If all outputs already exist with correct page counts, skip
    if nReqCh == 2
        allExist = isfile(outDe1) && isfile(outDe2) && isfile(outAz1) && isfile(outAz2);
    elseif needDeint
        allExist = isfile(outDe1) && isfile(outAz1);
    else
        allExist = isfile(outAz1);
    end
    if allExist
        try
            ok1 = count_tiff_pages_file(outAz1)==nZ;
            if needDeint
                ok1 = ok1 && count_tiff_pages_file(outDe1)==nZ*nT;
            end
            if nReqCh == 2
                ok1 = ok1 && count_tiff_pages_file(outDe2)==nZ*nT && count_tiff_pages_file(outAz2)==nZ;
            end
            if ok1
                fprintf("  [skip] outputs already exist\n");
                continue;
            end
        catch
            % fallthrough -> rewrite
        end
    end

    % --- do the work (streaming, starts saving now) ---
    try
        if nReqCh == 2
            stream_deint_and_avgz_twochan(fn, outDe1, outAz1, ch1_idx, ...
                                              outDe2, outAz2, ch2_idx, ...
                                              nCh, nZ, nT, orderMode, rowsPerStrip, forceUint16, int16Offset);
        elseif needDeint
            stream_deint_and_avgz_onechan(fn, outDe1, outAz1, ch1_idx, ...
                                              nCh, nZ, nT, orderMode, rowsPerStrip, forceUint16, int16Offset);
        else
            stream_avgz_only(fn, outAz1, ...
                             nZ, nT, orderMode, rowsPerStrip, forceUint16, int16Offset);
        end
        nOK = nOK + 1;
        fprintf("  [ok] wrote %s\n", ternary(needDeint, "deint+avgz", "avgz"));
    catch ME
        nSkip = nSkip + 1;
        fprintf("  [FAIL] %s\n", ME.message);
    end
end

fprintf("\n[DONE] ok=%d  skipped=%d\n", nOK, nSkip);
if nReqCh == 2 || nCh > 1
    fprintf("  ch%d deint: %s\n", ch1_id, outCh1);
end
fprintf("  ch%d avgz : %s\n", ch1_id, outCh1Z);
if nReqCh == 2
    fprintf("  ch%d deint: %s\n", ch2_id, outCh2);
    fprintf("  ch%d avgz : %s\n", ch2_id, outCh2Z);
end

end

%% ========================= CORE WORKER =========================
function stream_deint_and_avgz_twochan(rawFn, outDe1, outAz1, c1_idx, outDe2, outAz2, c2_idx, ...
                                       nCh, nZ, nT, orderMode, rowsPerStrip, forceUint16, int16Offset)
% Reads each (z,t) once and writes:
%  - deinterleave (z-major, then t) for ch1 and ch2
%  - avgz page per z for ch1 and ch2

% open input
tr = Tiff(rawFn,'r');  c0 = onCleanup(@() tr.close());
tr.setDirectory(1);
I0 = tr.read();
[H,W] = size(I0);
inClassRaw = class(I0);

% decide output class
if forceUint16 && strcmpi(inClassRaw,'int16')
    outClass = 'uint16';
else
    outClass = inClassRaw;
end
[bitsPerSample, sampleFormat] = class_to_tiff_format(outClass);

% recreate outputs cleanly
if isfile(outDe1), delete(outDe1); end
if isfile(outDe2), delete(outDe2); end
if isfile(outAz1), delete(outAz1); end
if isfile(outAz2), delete(outAz2); end

td1 = Tiff(outDe1,'w'); c1 = onCleanup(@() td1.close());
td2 = Tiff(outDe2,'w'); c2 = onCleanup(@() td2.close());
ta1 = Tiff(outAz1,'w'); c3 = onCleanup(@() ta1.close());
ta2 = Tiff(outAz2,'w'); c4 = onCleanup(@() ta2.close());

tag = struct();
tag.ImageLength = H;
tag.ImageWidth  = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.BitsPerSample = bitsPerSample;
tag.SampleFormat = sampleFormat;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = rowsPerStrip;

% counters
pDe = 0;

for zIdx = 1:nZ
    acc1 = zeros(H,W,'double');
    acc2 = zeros(H,W,'double');

    for tIdx = 1:nT
        % ---- channel 1 ----
        p1 = scanimage_page_index(c1_idx, zIdx, tIdx, nCh, nZ, nT, orderMode);
        tr.setDirectory(p1);
        I1 = double(tr.read());
        if forceUint16 && strcmpi(inClassRaw,'int16')
            I1 = I1 + double(int16Offset);
            I1 = min(max(I1,0),65535);
        end
        acc1 = acc1 + I1;

        % ---- channel 2 ----
        p2 = scanimage_page_index(c2_idx, zIdx, tIdx, nCh, nZ, nT, orderMode);
        tr.setDirectory(p2);
        I2 = double(tr.read());
        if forceUint16 && strcmpi(inClassRaw,'int16')
            I2 = I2 + double(int16Offset);
            I2 = min(max(I2,0),65535);
        end
        acc2 = acc2 + I2;

        % ---- write deinterleave pages (z-major, then t) ----
        pDe = pDe + 1;

        td1.setTag(tag); td1.write(cast(round(I1), outClass));
        if pDe < nZ*nT, td1.writeDirectory(); end

        td2.setTag(tag); td2.write(cast(round(I2), outClass));
        if pDe < nZ*nT, td2.writeDirectory(); end
    end

    % ---- write avgz page for this z ----
    A1 = acc1 ./ max(nT,1);
    A2 = acc2 ./ max(nT,1);

    ta1.setTag(tag); ta1.write(cast(round(A1), outClass));
    if zIdx < nZ, ta1.writeDirectory(); end

    ta2.setTag(tag); ta2.write(cast(round(A2), outClass));
    if zIdx < nZ, ta2.writeDirectory(); end

    fprintf("  [write] z %d/%d\n", zIdx, nZ);
end
end

%% ========================= SINGLE-CHANNEL WORKER =========================
function stream_deint_and_avgz_onechan(rawFn, outDe1, outAz1, c1_idx, ...
                                       nCh, nZ, nT, orderMode, rowsPerStrip, forceUint16, int16Offset)
% Single-channel version: extracts one channel (deinterleave if nCh>1),
% writes deinterleaved TIFF + avgZ TIFF.

% open input
tr = Tiff(rawFn,'r');  c0 = onCleanup(@() tr.close());
tr.setDirectory(1);
I0 = tr.read();
[H,W] = size(I0);
inClassRaw = class(I0);

% decide output class
if forceUint16 && strcmpi(inClassRaw,'int16')
    outClass = 'uint16';
else
    outClass = inClassRaw;
end
[bitsPerSample, sampleFormat] = class_to_tiff_format(outClass);

% recreate outputs cleanly
if isfile(outDe1), delete(outDe1); end
if isfile(outAz1), delete(outAz1); end

td1 = Tiff(outDe1,'w'); c1 = onCleanup(@() td1.close());
ta1 = Tiff(outAz1,'w'); c2 = onCleanup(@() ta1.close());

tag = struct();
tag.ImageLength = H;
tag.ImageWidth  = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.BitsPerSample = bitsPerSample;
tag.SampleFormat = sampleFormat;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = rowsPerStrip;

pDe = 0;

for zIdx = 1:nZ
    acc1 = zeros(H,W,'double');

    for tIdx = 1:nT
        p1 = scanimage_page_index(c1_idx, zIdx, tIdx, nCh, nZ, nT, orderMode);
        tr.setDirectory(p1);
        I1 = double(tr.read());
        if forceUint16 && strcmpi(inClassRaw,'int16')
            I1 = I1 + double(int16Offset);
            I1 = min(max(I1,0),65535);
        end
        acc1 = acc1 + I1;

        pDe = pDe + 1;
        td1.setTag(tag); td1.write(cast(round(I1), outClass));
        if pDe < nZ*nT, td1.writeDirectory(); end
    end

    % avgz page for this z
    A1 = acc1 ./ max(nT,1);
    ta1.setTag(tag); ta1.write(cast(round(A1), outClass));
    if zIdx < nZ, ta1.writeDirectory(); end

    fprintf("  [write] z %d/%d\n", zIdx, nZ);
end
end

%% ========================= AVGZ-ONLY WORKER (single channel, nCh==1) =========================
function stream_avgz_only(rawFn, outAz1, ...
                          nZ, nT, orderMode, rowsPerStrip, forceUint16, int16Offset)
% No deinterleaving needed (nCh==1). Just compute avgZ per Z slice.

tr = Tiff(rawFn,'r');  c0 = onCleanup(@() tr.close());
tr.setDirectory(1);
I0 = tr.read();
[H,W] = size(I0);
inClassRaw = class(I0);

if forceUint16 && strcmpi(inClassRaw,'int16')
    outClass = 'uint16';
else
    outClass = inClassRaw;
end
[bitsPerSample, sampleFormat] = class_to_tiff_format(outClass);

if isfile(outAz1), delete(outAz1); end
ta1 = Tiff(outAz1,'w'); c1 = onCleanup(@() ta1.close());

tag = struct();
tag.ImageLength = H;
tag.ImageWidth  = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.BitsPerSample = bitsPerSample;
tag.SampleFormat = sampleFormat;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = rowsPerStrip;

for zIdx = 1:nZ
    acc = zeros(H,W,'double');

    for tIdx = 1:nT
        p = scanimage_page_index(1, zIdx, tIdx, 1, nZ, nT, orderMode);
        tr.setDirectory(p);
        I = double(tr.read());
        if forceUint16 && strcmpi(inClassRaw,'int16')
            I = I + double(int16Offset);
            I = min(max(I,0),65535);
        end
        acc = acc + I;
    end

    A = acc ./ max(nT,1);
    ta1.setTag(tag); ta1.write(cast(round(A), outClass));
    if zIdx < nZ, ta1.writeDirectory(); end

    fprintf("  [write] z %d/%d\n", zIdx, nZ);
end
end

%% ========================= UTILITIES =========================
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function ensure_dir(p)
if ~isfolder(p), mkdir(p); end
end

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
c = onCleanup(@() t.close());

txt = "";
try, txt = txt + string(t.getTag('Software')) + newline; catch, end
try, txt = txt + string(t.getTag('ImageDescription')) + newline; catch, end

txt = replace(txt, char(8629), newline); % ↵
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
