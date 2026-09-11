function tile_deinterleave_avgz_260910()
%TILE_DEINTERLEAVE_AVGZ_260910  Deinterleave + avgz, sequential-read version.
%
% Same job and same outputs as tile_deinterleave_N_groupZavg_260131.m:
%   rawDir\ch<N>\        <base>_ch<N>_deint.tif   Z*T pages, z-major then t
%   rawDir\ch<N>\avgz\   <base>_ch<N>_avgz.tif    Z pages, mean over T
%
% ===================== WHERE THE TIME ACTUALLY GOES =====================
% MEASURED, after two wrong guesses.  On one 260804 tile (2 ch x 54 z x 50 t =
% 5400 pages of 512x512) written to E:, which does 129 MB/s sequential:
%
%   READ is NOT the problem.  setDirectory(p) per page vs sequential
%   nextDirectory(), alternated three times so neither gets a warmer cache:
%       setDirectory  29.8 / 29.5 / 29.8 s
%       nextDirectory 27.8 / 27.5 / 27.7 s      -> 1.08x, i.e. nothing.
%   The textbook objection to setDirectory -- that it rewalks the IFD chain and
%   is therefore quadratic -- does not bite here, because libtiff caches the
%   directory offsets it has already visited and the pages are read in ASCENDING
%   order.  It would bite on random or backward access.
%
%   WRITING is the problem.  One 2700-page 1.32 GB channel file:
%       raw fwrite, no TIFF      10.5 s   129 MB/s
%       Tiff, one strip/page    103.8 s    13 MB/s
%       Tiff, 64-row strips     169.2 s     8 MB/s   <- the ORIGINAL setting
%       Tiff as BigTIFF         134.7 s    10 MB/s
%       hand-laid TIFF           22.6 s    60 MB/s
%   MATLAB's Tiff class costs 10-16x the raw byte rate.
%
% SO THE SPEEDUP COMES FROM TWO CHANGES, both on the write path:
%   1. RowsPerStrip = ImageLength instead of 64            1.63x
%   2. the deint file written directly, IFDs appended last 4.6x
% Sequential reading is kept because it is free and slightly faster, not because
% it was ever the bottleneck.
%
% Also removed: count_tiff_pages_file walks every IFD (~2-4 s here) and the old
% script called it up to five times per tile, once for the strict check and once
% per existing output, just to decide whether to skip.  Now at most once, and the
% skip check reads the IFD table via imfinfo rather than the pixels.
%
% Measured and NOT a factor, so left alone: double vs single accumulation is
% 1.09 vs 1.12 s per 200 pages.
%
% THE DEINT FILE IS THE EXPENSIVE OUTPUT: 2 x 1.32 GB per tile against 2 x 28 MB
% for the avgz.  writeDeint = false is ~4x faster overall and is all the
% stitcher consumes.
%
% Runqi Zhang / 2026-09-10

%% ========================= USER SETTINGS =========================
rawDir      = "E:\260804_shiver_dbh\map";
chUse       = [];      % [] = EVERY channel the file saved, read per tile from
                       % SI.hChannels.channelSave.  Give a list ([1 3]) only to
                       % force a subset; an entry that was not actually saved is
                       % dropped with a note instead of failing the tile.

writeDeint  = false;   % avgz only.  The deint file is a regenerable cache --
                       % the raw TIFF is still on disk -- and nothing downstream
                       % reads it: the stitcher and both GUIs open ch<N>/avgz only.
                       % 35 s and 57 MB per tile, against 53 s and 2.5 GB.
maxFiles    = [];      % [] = every tile; a number caps it, for a test run
outRoot     = "";      % "" = write into rawDir; set a path to redirect (testing)

recursive   = false;
strictMeta  = true;    % verify pages == nCh*nZ*nT before doing any work
forceUint16 = true;    % int16 -> uint16 with an offset, no wrapping
int16Offset = 32768;
%% =================================================================

if ~isfolder(rawDir), error("Not a folder: %s", rawDir); end
if outRoot == "", outRoot = rawDir; end

files = list_tiffs(rawDir, recursive);
% Skip this script's own outputs.  Matching ANY ch<N> folder rather than two
% named ones keeps the filter independent of how many channels there turn out
% to be -- which is not known until each tile's header is read.
files = files(cellfun(@isempty, regexpi(cellstr(files), ['\' filesep 'ch\d+\' filesep], 'once')));
files = files(~contains(lower(files), lower(filesep + "avgz" + filesep)));
files = files(~contains(lower(files), "stitched"));
files = files(~contains(lower(files), "cache"));
files = unique(files);
if isempty(files), error("No TIFFs found in %s", rawDir); end
if ~isempty(maxFiles), files = files(1:min(maxFiles, numel(files))); end
fprintf("[scan] %d file(s)\n", numel(files));

nOK = 0; nSkip = 0; tAll = tic;

for k = 1:numel(files)
    fn = files(k);
    [~, baseName, ext] = fileparts(fn);
    fprintf("\n[%3d/%3d] %s%s\n", k, numel(files), baseName, ext);

    try
        meta = read_si_meta(fn);
    catch ME
        nSkip = nSkip + 1; fprintf("  [skip] meta read failed: %s\n", ME.message); continue
    end
    if ~isKey(meta,'SI.hChannels.channelSave')
        nSkip = nSkip + 1; fprintf("  [skip] no channelSave (not ScanImage?)\n"); continue
    end

    % channelSave IS the channel count: its length sets the page cycle, and its
    % order sets each channel's position within that cycle.  Nothing here is
    % hard-coded to two.
    chSave = num_vec(meta('SI.hChannels.channelSave')); chSave = chSave(:).';
    nCh = numel(chSave);
    if isempty(chUse)
        chIDs = chSave;
    else
        chIDs = intersect(chUse(:).', chSave, 'stable');
        gone  = setdiff(chUse(:).', chSave);
        if ~isempty(gone)
            fprintf("  [note] requested ch [%s] not saved in this tile -- ignored\n", num2str(gone));
        end
    end
    if isempty(chIDs)
        nSkip = nSkip + 1;
        fprintf("  [skip] none of the requested channels are in channelSave [%s]\n", num2str(chSave));
        continue
    end
    cIdx = arrayfun(@(id) find(chSave == id, 1), chIDs);  % position in the page cycle
    nSel = numel(chIDs);

    [nZ, nT, orderMode] = infer_zt(meta);
    if ~isfinite(nZ) || ~isfinite(nT) || nZ < 1 || nT < 1
        nSkip = nSkip + 1; fprintf("  [skip] bad nZ/nT\n"); continue
    end
    fprintf("  channelSave [%s]  writing ch [%s]  nZ %d  nT %d  order %s\n", ...
            num2str(chSave), num2str(chIDs), nZ, nT, orderMode);

    outDe = strings(nSel,1);  outAz = strings(nSel,1);
    for m = 1:nSel
        dCh = fullfile(outRoot, sprintf("ch%d", chIDs(m)));
        dZ  = fullfile(dCh, "avgz");
        ensure_dir(dCh); ensure_dir(dZ);
        outDe(m) = fullfile(dCh, baseName + sprintf("_ch%d_deint.tif", chIDs(m)));
        outAz(m) = fullfile(dZ,  baseName + sprintf("_ch%d_avgz.tif",  chIDs(m)));
    end

    % Existence check by SIZE, not by walking every IFD of every output.  The
    % old version called count_tiff_pages_file four times here, each a full pass
    % over the file -- i.e. it could spend longer deciding to skip than the
    % actual work would have taken on the avgz.
    want = cellstr(outAz(:));  wantPages = repmat(nZ, 1, nSel);
    if writeDeint
        want = [want; cellstr(outDe(:))];
        wantPages = [wantPages repmat(nZ*nT, 1, nSel)];
    end
    if all(cellfun(@isfile, want)) && outputs_look_complete(want, wantPages)
        fprintf("  [skip] outputs already present and the right size\n"); continue
    end

    try
        t0 = tic;
        one_tile(fn, nCh, nZ, nT, orderMode, cIdx, strictMeta, ...
                 outDe, outAz, writeDeint, forceUint16, int16Offset);
        nOK = nOK + 1;
        fprintf("  [ok] %.1f s\n", toc(t0));
    catch ME
        nSkip = nSkip + 1; fprintf("  [FAIL] %s\n", ME.message);
    end
end

fprintf("\n[DONE] ok %d, skipped %d, total %.1f s\n", nOK, nSkip, toc(tAll));
end

%% ============================ CORE ============================
function one_tile(rawFn, nCh, nZ, nT, orderMode, cIdx, strictMeta, ...
                  outDe, outAz, writeDeint, forceUint16, int16Offset)
% cIdx holds each wanted channel's POSITION IN THE PAGE CYCLE (1..nCh), and
% outDe/outAz hold one filename per entry.  Everything below is written over
% numel(cIdx), so one channel or four cost the same code.

nSel = numel(cIdx);

tr = Tiff(rawFn, 'r');  cR = onCleanup(@() tr.close());
I0 = tr.read();  [H, W] = size(I0);  inClass = class(I0);

nPages = nCh*nZ*nT;
if strictMeta
    n = count_pages(rawFn);
    if n ~= nPages
        error("page count %d != nCh*nZ*nT = %d", n, nPages);
    end
end

doOffset = forceUint16 && strcmpi(inClass,'int16');
if doOffset, outClass = 'uint16'; else, outClass = inClass; end

% SEQUENTIAL ORDER ONLY.  The reader is advanced with nextDirectory, so pages
% must be consumed in ascending order -- which is exactly the order the loop
% wants them in for a slice-major stack.  A volume-interleaved stack needs its
% deinterleaved pages in an order the file does not supply, so that case goes
% to the indexed path rather than being silently written out of order.
seqOK = (orderMode == "sliceMajor");
if writeDeint && ~seqOK
    warning(['volume-interleaved stack: the deinterleaved page order cannot be ' ...
             'produced by one sequential pass, so this tile uses the slow indexed ' ...
             'path.  avgz alone (writeDeint = false) stays fast.']);
    indexed_path(rawFn, outDe, outAz, nCh, nZ, nT, orderMode, cIdx, ...
                 H, W, outClass, doOffset, int16Offset);
    return
end

% WRITING.  MATLAB's Tiff class costs 10-16x the raw byte rate on this data,
% measured on one 2700-page 1.32 GB channel file written to E: (drive does
% 129 MB/s):
%       raw fwrite              10.5 s   129 MB/s
%       Tiff, one strip/page   103.8 s    13 MB/s
%       Tiff, 64-row strips    169.2 s     8 MB/s   <- the original setting
%       Tiff BigTIFF           134.7 s    10 MB/s
%       hand-laid TIFF          22.6 s    60 MB/s
% So the deint file is written directly: pixel pages streamed back to back, all
% IFDs appended at the end.  4.6x on the file that dominates the run.
classicMax = 4e9;                       % classic TIFF offsets are uint32
deintBytes = 8 + nZ*nT*H*W*2 + nZ*nT*130;
useFast    = deintBytes < classicMax;
if writeDeint && ~useFast
    fprintf("  [note] deint would exceed 4 GB -- using the slower BigTIFF writer\n");
end

az = cell(nSel,1);  de = cell(nSel,1);  td = cell(nSel,1);  tag = [];
for m = 1:nSel, az{m} = ftw_begin(outAz(m), H, W, outClass); end
if writeDeint
    if useFast
        for m = 1:nSel, de{m} = ftw_begin(outDe(m), H, W, outClass); end
    else
        tag = base_tag(H, W, outClass);
        for m = 1:nSel, td{m} = Tiff(outDe(m), 'w8'); end
    end
end

% page-cycle position -> output slot; 0 means this channel is not wanted and its
% page is stepped over without being read.
slotOf = zeros(1, nCh);  slotOf(cIdx) = 1:nSel;

acc  = zeros(H, W, nSel, 'single');
page = 0;  nDe = 0;
tr.setDirectory(1);

for z = 1:nZ
    acc(:) = 0;
    for t = 1:nT
        for c = 1:nCh
            page = page + 1;
            sl = slotOf(c);
            if sl > 0
                I = tr.read();
                if doOffset, I = uint16(int32(I) + int16Offset); end
                acc(:,:,sl) = acc(:,:,sl) + single(I);
                if writeDeint
                    if useFast
                        de{sl} = ftw_write(de{sl}, I);
                    else
                        td{sl}.setTag(tag);  td{sl}.write(I);
                        if nDe + 1 < nZ*nT, td{sl}.writeDirectory(); end
                    end
                end
            end
            if page < nPages, tr.nextDirectory(); end
        end
        if writeDeint, nDe = nDe + 1; end
    end
    for m = 1:nSel
        az{m} = ftw_write(az{m}, cast(round(acc(:,:,m) / nT), outClass));
    end
end

for m = 1:nSel, ftw_end(az{m}); end
if writeDeint
    if useFast
        for m = 1:nSel, ftw_end(de{m}); end
    else
        for m = 1:nSel, td{m}.close(); end
    end
end
end

function indexed_path(rawFn, outDe, outAz, nCh, nZ, nT, orderMode, cIdx, ...
                      H, W, outClass, doOffset, int16Offset)
% setDirectory per page: the only way to get deinterleaved output in a usable
% order out of a volume-interleaved stack.  Slow, and only reached when
% writeDeint is on -- avgz alone never needs it.
%
% This used to be called with filenames in parameters that expected already-open
% Tiff objects, and one argument short of its own signature, so the first
% volume-interleaved tile with writeDeint on would have errored.  Dormant,
% because writeDeint defaults to false.  It opens its own handles now.
nSel = numel(cIdx);
tr  = Tiff(rawFn, 'r');  cR = onCleanup(@() tr.close());
tag = base_tag(H, W, outClass);

td = cell(nSel,1);  ta = cell(nSel,1);
for m = 1:nSel
    td{m} = Tiff(outDe(m), deint_mode(H, W, nZ*nT));
    ta{m} = Tiff(outAz(m), 'w');
end

acc = zeros(H, W, nSel, 'single');  pDe = 0;
for z = 1:nZ
    acc(:) = 0;
    for t = 1:nT
        for m = 1:nSel
            tr.setDirectory(page_index(cIdx(m), z, t, nCh, nZ, nT, orderMode));
            I = tr.read();
            if doOffset, I = uint16(int32(I) + int16Offset); end
            acc(:,:,m) = acc(:,:,m) + single(I);
            td{m}.setTag(tag);  td{m}.write(I);
            if pDe + 1 < nZ*nT, td{m}.writeDirectory(); end
        end
        pDe = pDe + 1;
    end
    for m = 1:nSel
        ta{m}.setTag(tag);  ta{m}.write(cast(round(acc(:,:,m) / nT), outClass));
        if z < nZ, ta{m}.writeDirectory(); end
    end
end
for m = 1:nSel, td{m}.close();  ta{m}.close(); end
end

%% ========================== HELPERS ==========================
function tag = base_tag(H, W, cls)
tag = struct();
tag.ImageLength = H;  tag.ImageWidth = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
[tag.BitsPerSample, tag.SampleFormat] = class_fmt(cls);
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
% ONE strip per page.  MATLAB's writer pads a partial final strip out to a full
% RowsPerStrip, so a height that is not a multiple of it silently inflates every
% page; a single strip has no partial strip to pad and needs fewer offsets.
tag.RowsPerStrip = H;
end

function m = deint_mode(H, W, nPages)
% Classic TIFF tops out at 4 GB of offsets; switch to BigTIFF before that.
if H*W*2*nPages > 3.5e9, m = 'w8'; else, m = 'w'; end
end

function tf = outputs_look_complete(files, wantPages)
tf = true;
for k = 1:numel(files)
    d = dir(files{k});
    if isempty(d), tf = false; return; end
    info = imfinfo(files{k});          % reads the IFD table, not the pixels
    if numel(info) ~= wantPages(k), tf = false; return; end
end
end

function ensure_dir(p), if ~isfolder(p), mkdir(p); end, end

function files = list_tiffs(root, recursive)
if recursive
    c = [dir(fullfile(root,"**","*.tif")); dir(fullfile(root,"**","*.tiff"))];
else
    c = [dir(fullfile(root,"*.tif")); dir(fullfile(root,"*.tiff"))];
end
files = strings(numel(c),1);
for k = 1:numel(c), files(k) = string(fullfile(c(k).folder, c(k).name)); end
end

function kv = read_si_meta(f)
t = Tiff(f,'r'); c = onCleanup(@() t.close());
txt = "";
try, txt = txt + string(t.getTag('Software')) + newline; catch, end
try, txt = txt + string(t.getTag('ImageDescription')) + newline; catch, end
txt = replace(txt, char(8629), newline);
kv = containers.Map('KeyType','char','ValueType','char');
L = splitlines(txt);
for i = 1:numel(L)
    s = strtrim(L(i));
    if startsWith(s, "SI.")
        e = strfind(s, '=');
        if ~isempty(e)
            kv(char(strtrim(extractBefore(s, e(1))))) = char(strtrim(extractAfter(s, e(1))));
        end
    end
end
end

function v = num_vec(s)
tok = regexp(char(s), '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
v = str2double(tok);
end

function [nZ, nT, orderMode] = infer_zt(kv)
orderMode = "sliceMajor";
if isKey(kv,'SI.hFastZ.enable')
    e = lower(string(kv('SI.hFastZ.enable')));
    if e == "true" || e == "1", orderMode = "volumeInterleaved"; end
end
nT = NaN;
if isKey(kv,'SI.hStackManager.framesPerSlice')
    nT = str2double(kv('SI.hStackManager.framesPerSlice'));
end
nZ = NaN;
if isKey(kv,'SI.hStackManager.zs')
    z = num_vec(kv('SI.hStackManager.zs'));
    if ~isempty(z), nZ = numel(z); end
end
if ~isfinite(nZ) && isKey(kv,'SI.hStackManager.numSlices')
    nZ = str2double(kv('SI.hStackManager.numSlices'));
end
end

function n = count_pages(f)
t = Tiff(f,'r'); c = onCleanup(@() t.close());
n = 1;
while ~t.lastDirectory(), t.nextDirectory(); n = n + 1; end
end

function p = page_index(cIdx, z, t, nCh, nZ, nT, orderMode)
if orderMode == "volumeInterleaved"
    frame = (t-1)*nZ + z;
else
    frame = (z-1)*nT + t;
end
p = (frame-1)*nCh + cIdx;
end

function [bits, fmt] = class_fmt(cls)
switch cls
    case 'uint8',  bits = 8;  fmt = Tiff.SampleFormat.UInt;
    case 'uint16', bits = 16; fmt = Tiff.SampleFormat.UInt;
    case 'uint32', bits = 32; fmt = Tiff.SampleFormat.UInt;
    case 'int16',  bits = 16; fmt = Tiff.SampleFormat.Int;
    case 'single', bits = 32; fmt = Tiff.SampleFormat.IEEEFP;
    case 'double', bits = 64; fmt = Tiff.SampleFormat.IEEEFP;
    otherwise, error("Unsupported TIFF class: %s", cls);
end
end


%% ==================== DIRECT TIFF WRITER ====================
% Classic little-endian TIFF, uncompressed, one strip per page, every page the
% same size.  Pixels are streamed back to back from byte 8 and all IFDs are
% appended at the end, so nothing ever seeks backwards except the 4-byte patch
% of the first-IFD pointer at close.
%
% Only valid while the file stays under 4 GB -- classic TIFF offsets are uint32.
% The caller checks that and falls back to Tiff('w8') when it would not.

function w = ftw_begin(fn, H, W, cls)
if isfile(fn), delete(fn); end
w.fid = fopen(fn, 'w');
if w.fid < 0, error("cannot open %s for writing", fn); end
w.H = H; w.W = W; w.cls = cls; w.n = 0; w.fn = fn;
fwrite(w.fid, [73 73], 'uint8');   % 'II'
fwrite(w.fid, 42, 'uint16');       % classic TIFF
fwrite(w.fid, 0, 'uint32');        % first-IFD offset, patched in ftw_end
end

function w = ftw_write(w, P)
% TRANSPOSED.  MATLAB is column-major, TIFF stores raster rows.  Writing P as-is
% yields a file whose metadata is entirely valid and whose every page is silently
% transposed -- invisible on a square tile until you compare pixels, which is
% exactly how it got past the first test here.
w.n = w.n + 1;
fwrite(w.fid, P.', w.cls);
end

function ftw_end(w)
if isempty(w) || w.fid < 0, return; end
H = w.H; W = w.W; n = w.n; nBytes = H*W*2;
switch w.cls
    case 'uint16', sfmt = 1;
    case 'int16',  sfmt = 2;
    otherwise, error("ftw: unsupported class %s", w.cls);
end
ifd0 = ftell(w.fid);
ifdSize = 2 + 10*12 + 4;
for k = 1:n
    fwrite(w.fid, 10, 'uint16');
    ftw_tag(w.fid, 256, 4, W);                    % ImageWidth
    ftw_tag(w.fid, 257, 4, H);                    % ImageLength
    ftw_tag(w.fid, 258, 3, 16);                   % BitsPerSample
    ftw_tag(w.fid, 259, 3, 1);                    % Compression = none
    ftw_tag(w.fid, 262, 3, 1);                    % BlackIsZero
    ftw_tag(w.fid, 273, 4, 8 + (k-1)*nBytes);     % StripOffsets
    ftw_tag(w.fid, 277, 3, 1);                    % SamplesPerPixel
    ftw_tag(w.fid, 278, 4, H);                    % RowsPerStrip = whole page
    ftw_tag(w.fid, 279, 4, nBytes);               % StripByteCounts
    ftw_tag(w.fid, 339, 3, sfmt);                 % SampleFormat
    if k < n, nxt = ifd0 + k*ifdSize; else, nxt = 0; end
    fwrite(w.fid, nxt, 'uint32');
end
fseek(w.fid, 4, 'bof');
fwrite(w.fid, ifd0, 'uint32');
fclose(w.fid);
end

function ftw_abort(w)
if ~isempty(w) && isstruct(w) && w.fid >= 0
    fclose(w.fid);
    if isfile(w.fn), delete(w.fn); end
end
end

function ftw_tag(fid, tag, type, val)
fwrite(fid, tag, 'uint16'); fwrite(fid, type, 'uint16'); fwrite(fid, 1, 'uint32');
if type == 3      % a SHORT lives in the first 2 bytes of the 4-byte value field
    fwrite(fid, val, 'uint16'); fwrite(fid, 0, 'uint16');
else
    fwrite(fid, val, 'uint32');
end
end
