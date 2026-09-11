function tile_deinterleave_avgz_260910()

rawDir      = "E:\260804_shiver_dbh\map";
chUse       = [];
writeDeint  = false; 

%%
maxFiles    = [];    outRoot     = "";
strictMeta  = true;  recursive   = false;
forceUint16 = true;  int16Offset = 32768;

if ~isfolder(rawDir), error("Not a folder: %s", rawDir); end
if outRoot == "", outRoot = rawDir; end

files = list_tiffs(rawDir, recursive);
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
    cIdx = arrayfun(@(id) find(chSave == id, 1), chIDs);
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

function one_tile(rawFn, nCh, nZ, nT, orderMode, cIdx, strictMeta, ...
                  outDe, outAz, writeDeint, forceUint16, int16Offset)

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

seqOK = (orderMode == "sliceMajor");
if writeDeint && ~seqOK
    warning(['volume-interleaved stack: the deinterleaved page order cannot be ' ...
             'produced by one sequential pass, so this tile uses the slow indexed ' ...
             'path.  avgz alone (writeDeint = false) stays fast.']);
    indexed_path(rawFn, outDe, outAz, nCh, nZ, nT, orderMode, cIdx, ...
                 H, W, outClass, doOffset, int16Offset);
    return
end

classicMax = 4e9;
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

function tag = base_tag(H, W, cls)
tag = struct();
tag.ImageLength = H;  tag.ImageWidth = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
[tag.BitsPerSample, tag.SampleFormat] = class_fmt(cls);
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = H;
end

function m = deint_mode(H, W, nPages)
if H*W*2*nPages > 3.5e9, m = 'w8'; else, m = 'w'; end
end

function tf = outputs_look_complete(files, wantPages)
tf = true;
for k = 1:numel(files)
    d = dir(files{k});
    if isempty(d), tf = false; return; end
    info = imfinfo(files{k});
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

function w = ftw_begin(fn, H, W, cls)
if isfile(fn), delete(fn); end
w.fid = fopen(fn, 'w');
if w.fid < 0, error("cannot open %s for writing", fn); end
w.H = H; w.W = W; w.cls = cls; w.n = 0; w.fn = fn;
fwrite(w.fid, [73 73], 'uint8');
fwrite(w.fid, 42, 'uint16');
fwrite(w.fid, 0, 'uint32');
end

function w = ftw_write(w, P)
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
    ftw_tag(w.fid, 256, 4, W);
    ftw_tag(w.fid, 257, 4, H);
    ftw_tag(w.fid, 258, 3, 16);
    ftw_tag(w.fid, 259, 3, 1);
    ftw_tag(w.fid, 262, 3, 1);
    ftw_tag(w.fid, 273, 4, 8 + (k-1)*nBytes);
    ftw_tag(w.fid, 277, 3, 1);
    ftw_tag(w.fid, 278, 4, H);
    ftw_tag(w.fid, 279, 4, nBytes);
    ftw_tag(w.fid, 339, 3, sfmt);
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
if type == 3
    fwrite(fid, val, 'uint16'); fwrite(fid, 0, 'uint16');
else
    fwrite(fid, val, 'uint32');
end
end
