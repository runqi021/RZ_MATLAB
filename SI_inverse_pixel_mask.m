clear; close all; clc;
%%
fn="C:\Users\Admin\Downloads\shi_250115\test\preBotZ_1200_400_z50-700_15-55lp_00001.tif";

% si_dwell_correct_tiff
% - Reads ScanImage header from info(1).Software
% - Determines:
%     (1) interleave (page-interleaved channels) from header + page count
%     (2) raw type (int16 vs uint16) from TIFF SampleFormat/BitsPerSample
% - Applies pixel dwell correction using SI.hScan2D.mask / SI.hScan2D.sampleRate
% - Writes corrected multipage TIFF with SAME sample format as raw (int16 or uint16)

info = imfinfo(fn);
assert(~isempty(info), "imfinfo failed: %s", fn);

% ---- header from Software (your case) ----
hdr = string(info(1).Software);
hdr = replace(hdr, "↵", newline);

% ---- dwell fields ----
sampleRate = getSI_scalar(hdr, "SI.hScan2D.sampleRate");   % Hz
mask       = getSI_array (hdr, "SI.hScan2D.mask");         % supports [..;..;..]
assert(numel(mask) >= 16, "Parsed mask length=%d (should be ~512).", numel(mask));

% ---- channel saved / interleave inference ----
channelsActive = getSI_numlist_maybe(hdr, "SI.hChannels.channelsActive"); % e.g. 2
channelSave    = getSI_numlist_maybe(hdr, "SI.hChannels.channelSave");    % e.g. 2
savedCh        = inferSavedChannels(channelsActive, channelSave);

nPages = numel(info);

% ---- TIFF storage type ----
[sf, bits, spp] = readTiffNumericTags(fn);
rawClass = inferRawClassFromTags(sf, bits);
if rawClass == ""
    % fallback
    rawClass = string(class(imread(fn, 1, "Info", info)));
end

isInterleavedPages = (numel(savedCh) > 1) && (spp == 1) && (mod(nPages, numel(savedCh)) == 0);

fprintf("File: %s\n", fn);
fprintf("Raw class inferred: %s (SampleFormat=%s, Bits=%d)\n", rawClass, mat2str(sf), bits);
fprintf("Saved channels inferred: [%s]\n", strjoin(string(savedCh), " "));
fprintf("SamplesPerPixel=%d, pages=%d, InterleavedPages=%d\n", spp, nPages, isInterleavedPages);

% ---- dwell gain along X ----
H = info(1).Height;
W = info(1).Width;

mask_mapped = mapMaskToWidth(mask, W);                 % W x 1 or 1 x W
mask_row    = mask_mapped(:)';                         % 1 x W

gain_inv = mask_row ./ mean(mask_row);                 % 1 x W (mean-normalized inverse)
gx = reshape(single(gain_inv), 1, W);                  % 1 x W for broadcasting

fprintf("inverse gain (mask/mean): min=%.3f max=%.3f\n", min(gain_inv), max(gain_inv));

% ---- output ----
[inDir, base, ext] = fileparts(fn);
outPath = fullfile(inDir, base + "_maskInv" + ext);

% If file exists, try delete; if locked, write temp then move at end
tmpPath = outPath;
useTemp = false;
if exist(outPath, "file")
    try
        delete(outPath);
    catch
        tmpPath = fullfile(inDir, base + "_dwellFlat_TMP" + ext);
        if exist(tmpPath, "file"), delete(tmpPath); end
        useTemp = true;
    end
end

%% ---- write with Tiff (supports int16) ----
t = openTiffWithRetry(tmpPath);

tag.ImageLength = H;
tag.ImageWidth  = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = spp; % usually 1
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.BitsPerSample = 16;

if rawClass == "int16"
    tag.SampleFormat = Tiff.SampleFormat.Int;
elseif rawClass == "uint16"
    tag.SampleFormat = Tiff.SampleFormat.UInt;
else
    % If you ever hit this, you can extend, but your data are 16-bit
    t.close();
    error("Unsupported rawClass for this writer: %s", rawClass);
end

% Keep Software metadata (optional, but useful)
tag.Software = char(string(info(1).Software));

% RowsPerStrip
try
    tag.RowsPerStrip = info(1).RowsPerStrip;
catch
    tag.RowsPerStrip = min(H, 128);
end

for k = 1:nPages
    I = imread(fn, k, "Info", info);

    % dwell correction along X
    Y = single(I) .* gx;  % HxW (or HxWxC if spp>1; most SI uses spp=1)

    % cast back to SAME type
    if rawClass == "int16"
        Y = int16(max(-32768, min(32767, round(Y))));
    else
        Y = uint16(max(0, min(65535, round(Y))));
    end

    t.setTag(tag);
    t.write(Y);

    if k < nPages
        t.writeDirectory();
    end

    if mod(k,200) == 0 || k == nPages
        fprintf("Wrote %d / %d pages\n", k, nPages);
    end
end
t.close();

% Move temp into final if needed
if useTemp
    try
        movefile(tmpPath, outPath, "f");
    catch
        warning("Wrote temp file but could not move it. Temp is at: %s", tmpPath);
        outPath = tmpPath;
    end
end

fprintf("DONE. Output: %s\n", outPath);

%%
figure; histogram(Y);
figure; histogram(I);
%%

%% ===================== helpers =====================

function t = openTiffWithRetry(path)
% Always BigTIFF for safety.
nTry = 10;
for i = 1:nTry
    try
        t = Tiff(path, 'w8');
        return;
    catch
        pause(0.1);
    end
end
error("Failed to open output TIFF for writing: %s (is it locked?)", path);
end

function rhs = getSI_rhs(h, key)
pat = key + " = ";
idx = strfind(h, pat);
assert(~isempty(idx), "Key not found: %s", key);
s = extractAfter(h, idx(1) + strlength(pat) - 1);
nl = strfind(s, newline);
if isempty(nl), rhs = s; else, rhs = extractBefore(s, nl(1)); end
rhs = string(rhs);
end

function v = getSI_scalar(h, key)
rhs = getSI_rhs(h, key);
v = str2double(strtrim(rhs));
assert(~isnan(v), "Failed to parse scalar for %s", key);
end

function a = getSI_array(h, key)
rhs = strtrim(getSI_rhs(h, key));
assert(startsWith(rhs, "[") && contains(rhs, "]"), "Expected bracket array for %s", key);
body = extractBetween(rhs, "[", "]");
if iscell(body), body = body{1}; end
body = string(body);
% critical: convert ';' and ',' into spaces so sscanf reads all values
body = replace(body, ";", " ");
body = replace(body, ",", " ");
body = replace(body, newline, " ");
body = strip(body);
a = sscanf(char(body), "%f").';
end

function nums = getSI_numlist_maybe(h, key)
rhs = strtrim(getSI_rhs(h, key));
if startsWith(rhs, "[") && contains(rhs, "]")
    body = extractBetween(rhs, "[", "]");
    if iscell(body), body = body{1}; end
    body = replace(string(body), ";", " ");
    nums = sscanf(char(body), "%f").';
else
    v = str2double(rhs);
    if isnan(v), nums = []; else, nums = v; end
end
end

function saved = inferSavedChannels(chAct, chSave)
% In your case channelsActive=2 and channelSave=2 -> saved channel is [2]
if isempty(chSave)
    saved = chAct;
else
    saved = chSave;
end
saved = unique(saved(saved >= 1 & isfinite(saved)));
if isempty(saved), saved = 1; end
end

function maskOut = mapMaskToWidth(maskIn, Wout)
maskIn = double(maskIn(:));
Win = numel(maskIn);
assert(Win >= 2, "mask length %d too small.", Win);
if Win == Wout
    maskOut = maskIn;
else
    xi = linspace(1, Win, Wout);
    maskOut = interp1(1:Win, maskIn, xi, "linear", "extrap");
    maskOut = max(maskOut, eps);
end
end

function [sf, bits, spp] = readTiffNumericTags(f)
sf = []; bits = NaN; spp = 1;
try
    tt = Tiff(f,'r');
    try sf   = tt.getTag('SampleFormat');    catch, sf = []; end
    try bits = tt.getTag('BitsPerSample');   catch, bits = NaN; end
    try spp  = tt.getTag('SamplesPerPixel'); catch, spp = 1; end
    tt.close();
catch
end
end

function cls = inferRawClassFromTags(sf, bits)
cls = "";
if ~isempty(sf) && ~isnan(bits)
    if sf == 2 && bits == 16, cls = "int16";  return; end
    if sf == 1 && bits == 16, cls = "uint16"; return; end
end
end
