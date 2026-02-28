groupZ_avg_every_100frames_tiffconvention

%%
function groupZ_avg_every_100frames_tiffconvention()
% Average every 100 frames -> 1 frame (groupZ), and write using Tiff tags
% (same convention as your stitcher). Supports int16.
%
% Output: <name>_groupZ.tif in the same folder.

inDir = "E:\live_251104_definitely_not_sst_tdtomato_fitc\raw";
framesPerGroup = 100;

% TIFF writing
rowsPerStrip = 64;

tifs = [dir(fullfile(inDir, "*.tif")); dir(fullfile(inDir, "*.tiff"))];
assert(~isempty(tifs), "No TIFF files found in: %s", inDir);

for k = 1:numel(tifs)
    inFile = fullfile(tifs(k).folder, tifs(k).name);
    [~, base, ext] = fileparts(inFile);
    outTif = fullfile(tifs(k).folder, base + "_groupZ" + ext);

    % ---- open input with Tiff (like your stitcher) ----
    tin = Tiff(inFile, 'r');
    cleanupIn = onCleanup(@() tin.close()); %#ok<NASGU>

    % read first page for size/class
    tin.setDirectory(1);
    I0 = tin.read();
    [H,W] = size(I0);
    inClass = class(I0);

    % count pages
    nF = count_tiff_pages(tin);
    assert(mod(nF, framesPerGroup) == 0, ...
        "File %s has %d frames; not divisible by %d.", tifs(k).name, nF, framesPerGroup);

    nZ = nF / framesPerGroup;

    fprintf("\n[%d/%d] %s | %d frames -> %d Z slices (avg %d) | %dx%d | %s\n", ...
        k, numel(tifs), tifs(k).name, nF, nZ, framesPerGroup, H, W, inClass);

    % ---- output TIFF setup (same tag convention) ----
    [bitsPerSample, sampleFormat] = class_to_tiff_format(inClass);
    bytesPerPixel = bitsPerSample/8;

    estBytes = double(H) * double(W) * double(nZ) * double(bytesPerPixel);
    useBigTiff = estBytes > 3.5e9;

    if isfile(outTif), delete(outTif); end
    if useBigTiff
        tout = Tiff(outTif, 'w8');
    else
        tout = Tiff(outTif, 'w');
    end
    cleanupOut = onCleanup(@() tout.close()); %#ok<NASGU>

    tag.ImageLength = H;
    tag.ImageWidth  = W;
    tag.Photometric = Tiff.Photometric.MinIsBlack;
    tag.SamplesPerPixel = 1;
    tag.BitsPerSample = bitsPerSample;
    tag.SampleFormat = sampleFormat;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tag.Compression = Tiff.Compression.None;
    tag.RowsPerStrip = rowsPerStrip;
    tag.Software = 'MATLAB groupZ avg (Tiff tags; like stitcher)';

    % integer clipping bounds (supports int16 negatives)
    if isinteger(I0)
        vmin = double(intmin(inClass));
        vmax = double(intmax(inClass));
    end

    % ---- process group-by-group and write page-by-page ----
    for iz = 1:nZ
        fprintf("[write] z %d/%d\n", iz, nZ);

        idx0 = (iz-1)*framesPerGroup + 1;
        idx1 = iz*framesPerGroup;

        acc = zeros(H, W, 'double');
        for t = idx0:idx1
            tin.setDirectory(t);
            acc = acc + double(tin.read());
        end

        out = acc / framesPerGroup;

        % cast back to original unit/class
        if isinteger(I0)
            out = min(max(out, vmin), vmax);
            out = cast(round(out), inClass);
        else
            out = cast(out, inClass);
        end

        tout.setTag(tag);
        tout.write(out);

        if iz < nZ
            tout.writeDirectory();
        end
    end

    fprintf("[DONE] %s\n", outTif);
end

fprintf("\nAll done.\n");
end

%% ===================== helpers (same style as yours) =====================

function n = count_tiff_pages(t)
t.setDirectory(1);
n = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    n = n + 1;
end
t.setDirectory(1);
end

function [bitsPerSample, sampleFormat] = class_to_tiff_format(cls)
switch cls
    case 'uint8'
        bitsPerSample = 8;  sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint16'
        bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint32'
        bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.UInt;
    case 'int8'
        bitsPerSample = 8;  sampleFormat = Tiff.SampleFormat.Int;
    case 'int16'
        bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.Int;
    case 'int32'
        bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.Int;
    case 'single'
        bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'double'
        bitsPerSample = 64; sampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        error("Unsupported TIFF class: %s", cls);
end
end
