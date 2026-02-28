bg_subtract_firstpeak_perZ_drop1_cropLR20

%%
function bg_subtract_firstpeak_perZ_drop1_cropLR20()
% bg_subtract_firstpeak_perZ_drop1_cropLR20
% ------------------------------------------------------------
% For a Z-stack TIFF with 50 frames per Z:
%   - DROP the first frame of EACH Z group  -> write 49 frames/Z
%   - ALSO crop 20 pixels from LEFT and RIGHT of every frame (W -> W-40)
%   - For each Z, estimate bg(z) as LEFT-MOST peak of smoothed histogram
%     (same findpeaks logic as your code), using ONLY kept frames and ONLY
%     cropped pixels.
%   - Subtract bg(z) from all kept frames in that Z group.
%   - Write output TIFF preserving sample format (int16 stays int16, etc).
%
% Output pages: nZ * (framesPerZ-1)
% Output width: W0 - 40
% ------------------------------------------------------------

%% ---------------- USER SETTINGS ----------------
tiffPath = "C:\Users\Admin\Downloads\shi_250115\wt_fitc\live_not_sst_fitc_1.2x_930nm_7-26lp_z220_-10_100f_col01_row01_x-900_y101_00001.tif" ;

framesPerZ    = 100;
dropFirstPerZ = 1;          % drop 1st frame in each Z group
cropLR        = 20;         % <-- crop 20 px on left and 20 px on right

subsample   = 4;
nBins       = 512;
smoothW     = 9;
minPromFrac = 0.01;

clipNegative = false;

zStart_um = 0;
zEnd_um   = 220;
zStep_um  = 10;

outSuffix = "_minusBgPerZ_drop1_cropLR20";

%% ---------------- CHECK INPUT + METADATA ----------------
assert(isfile(tiffPath), "File not found: %s", tiffPath);
info   = imfinfo(tiffPath);
nPages = numel(info);
H0     = info(1).Height;
W0     = info(1).Width;

fprintf("[bgPerZ] %s\n", tiffPath);
fprintf("[bgPerZ] Pages: %d\n", nPages);

if mod(nPages, framesPerZ) ~= 0
    error("nPages (%d) not divisible by framesPerZ (%d).", nPages, framesPerZ);
end
nZ = nPages / framesPerZ;

keepPerZ = framesPerZ - dropFirstPerZ;
if keepPerZ <= 0
    error("dropFirstPerZ=%d too large (framesPerZ=%d).", dropFirstPerZ, framesPerZ);
end

% crop columns
cols = (1+cropLR):(W0-cropLR);
if isempty(cols) || numel(cols) < 10
    error("cropLR=%d too large for width=%d.", cropLR, W0);
end
W = numel(cols);
H = H0;

[inClass, bits, sf] = detect_tiff_class_from_tags(tiffPath);
fprintf("[bgPerZ] TIFF tags: class=%s, Bits=%d, SampleFormat=%s\n", inClass, bits, sf);

% z vector (QC only)
z_expected = (zStart_um:zStep_um:zEnd_um);
if numel(z_expected) == nZ
    z_um = z_expected(:);
else
    z_um = linspace(zStart_um, zEnd_um, nZ).';
end

%% ---------------- PASS 1: ESTIMATE bg(z) USING KEPT FRAMES + CROPPED IMG ----------------
bg = zeros(nZ,1);
bg_nSamples = zeros(nZ,1);

tIn = Tiff(tiffPath, 'r');
cIn = onCleanup(@() tIn.close());

for iz = 1:nZ
    vals = [];

    firstKeep_f = dropFirstPerZ + 1;
    page0 = (iz-1)*framesPerZ + firstKeep_f;

    % bounds from first kept frame (cropped)
    setDirectory(tIn, page0);
    I0  = double(tIn.read());
    I0  = I0(:, cols);
    I0s = I0(1:subsample:end, 1:subsample:end);
    lo0 = prctile(I0s(:), 0.5);
    hi0 = prctile(I0s(:), 99.5);
    if hi0 <= lo0, hi0 = lo0 + 1; end

    % gather samples from kept frames (cropped)
    for f = firstKeep_f:framesPerZ
        pageIdx = (iz-1)*framesPerZ + f;
        setDirectory(tIn, pageIdx);
        I  = double(tIn.read());
        I  = I(:, cols);
        Is = I(1:subsample:end, 1:subsample:end);
        vals = [vals; Is(:)]; %#ok<AGROW>
    end

    vals = vals(isfinite(vals));
    bg_nSamples(iz) = numel(vals);

    lo = min(lo0, prctile(vals, 0.5));
    hi = max(hi0, prctile(vals, 99.5));
    if hi <= lo, hi = lo + 1; end

    edges = linspace(lo, hi, nBins+1);
    counts = histcounts(vals, edges);
    bin_centers = (edges(1:end-1) + edges(2:end))/2;

    counts_smooth = smoothdata(counts, 'gaussian', smoothW);

    [pks, locs] = findpeaks(counts_smooth, bin_centers, ...
        'MinPeakProminence', max(counts_smooth)*minPromFrac);

    if isempty(pks)
        warning("[bgPerZ] iz=%d: no peaks; fallback to lo.", iz);
        bg(iz) = lo;
    else
        [~, idx_min] = min(locs);
        bg(iz) = locs(idx_min);
    end

    if mod(iz, max(1,floor(nZ/10)))==0 || iz==nZ
        fprintf("[bgPerZ] iz=%d/%d bg=%.2f (n=%d)\n", iz, nZ, bg(iz), bg_nSamples(iz));
    end
end

%% ---------------- PASS 2: APPLY bg(z) TO KEPT FRAMES, CROP LR, WRITE ----------------
[folder,name,ext] = fileparts(tiffPath);
outPath = fullfile(folder, name + outSuffix + ext);
csvPath = fullfile(folder, name + "_bgPerZ_drop1_cropLR20.csv");

if isfile(outPath), delete(outPath); end

tOut = Tiff(outPath, 'w');
cOut = onCleanup(@() tOut.close());

tag = make_tiff_tag_custom(H, W, inClass, bits, sf);

pageOut = 0;
nPagesOut = nZ * keepPerZ;

for iz = 1:nZ
    b = bg(iz);
    firstKeep_f = dropFirstPerZ + 1;

    for f = firstKeep_f:framesPerZ
        pageIdx = (iz-1)*framesPerZ + f;
        setDirectory(tIn, pageIdx);
        I  = tIn.read();
        I  = I(:, cols);              % crop LR

        Id = double(I) - b;

        if clipNegative
            Id(Id < 0) = 0;
        end

        Iout = cast_like(Id, inClass);

        pageOut = pageOut + 1;
        tOut.setTag(tag);
        tOut.write(Iout);

        if pageOut < nPagesOut
            tOut.writeDirectory();
        end
    end

    if mod(iz, max(1,floor(nZ/10)))==0 || iz==nZ
        fprintf("[write] z %d/%d\n", iz, nZ);
    end
end

fprintf("Wrote: %s\n", outPath);
fprintf("Output pages: %d (= %d Z * %d kept)\n", nPagesOut, nZ, keepPerZ);
fprintf("Output size: %d x %d (cropped LR %d px each side)\n", H, W, cropLR);

%% ---------------- SAVE CSV + QC PLOT ----------------
z_expected = (zStart_um:zStep_um:zEnd_um);
if numel(z_expected) == nZ, z_um = z_expected(:); end

Ttbl = table((1:nZ).', z_um, bg, bg_nSamples, ...
    'VariableNames', {'z_index','z_um','bg_firstPeak','n_samples'});
writetable(Ttbl, csvPath);
fprintf("Saved: %s\n", csvPath);

figure('Color','w');
plot(z_um, bg, '-o', 'LineWidth', 1);
xlabel('Z (\mum)'); ylabel('Background (first peak, ADU)');
title(sprintf('Per-Z bg (first peak), drop1 frame, cropLR=%d', cropLR));
grid on;

end

%% ================= HELPERS =================
function [cls,bits,fmtStr] = detect_tiff_class_from_tags(fp)
t = Tiff(fp,'r');
bits = t.getTag('BitsPerSample');
sf   = t.getTag('SampleFormat');
t.close();

switch sf
    case Tiff.SampleFormat.UInt
        fmtStr='UInt';
        if bits==8, cls='uint8';
        elseif bits==16, cls='uint16';
        elseif bits==32, cls='uint32';
        else, cls='uint16';
        end
    case Tiff.SampleFormat.Int
        fmtStr='Int';
        if bits==16, cls='int16';
        elseif bits==32, cls='int32';
        else, cls='int16';
        end
    case Tiff.SampleFormat.IEEEFP
        fmtStr='IEEEFP';
        cls='single';
    otherwise
        fmtStr='Unknown';
        cls='int16';
end
end

function tag = make_tiff_tag_custom(H, W, cls, bits, fmtStr)
tag.ImageLength = H;
tag.ImageWidth  = W;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = min(H, 256);

tag.BitsPerSample = bits;
switch fmtStr
    case 'UInt'
        tag.SampleFormat = Tiff.SampleFormat.UInt;
    case 'Int'
        tag.SampleFormat = Tiff.SampleFormat.Int;
    case 'IEEEFP'
        tag.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        if startsWith(cls,'u'), tag.SampleFormat = Tiff.SampleFormat.UInt;
        elseif startsWith(cls,'int'), tag.SampleFormat = Tiff.SampleFormat.Int;
        else, tag.SampleFormat = Tiff.SampleFormat.IEEEFP;
        end
end
end

function out = cast_like(x, cls)
switch cls
    case 'int16'
        x = min(max(x, -32768), 32767);
        out = int16(round(x));
    case 'int32'
        x = min(max(x, double(intmin('int32'))), double(intmax('int32')));
        out = int32(round(x));
    case 'uint16'
        x = min(max(x, 0), 65535);
        out = uint16(round(x));
    case 'uint8'
        x = min(max(x, 0), 255);
        out = uint8(round(x));
    case 'single'
        out = single(x);
    otherwise
        x = min(max(x, -32768), 32767);
        out = int16(round(x));
end
end
