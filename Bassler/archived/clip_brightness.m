clip40_rescale_overwrite_all_avi_fast

%%
function clip40_rescale_overwrite_all_avi_fast()
% Faster: chunked read + LUT mapping, overwrite AVI

inDir   = "D:\RZ\260112_whisking_KA_dorsal_20nL\250112_whisking_KA_dorsal_20nL-RZ-2026-01-12\videos";
clipMax = 40;
scale   = 255 / clipMax;

files = dir(fullfile(inDir, "*.avi"));
if isempty(files), error("No AVI files found in: %s", inDir); end
fprintf("Found %d AVI files.\n", numel(files));

% ---- LUT for uint8 input: y = min(x,clipMax)*scale ----
% if your frames are uint8 grayscale already, this makes mapping O(1)
lut = uint8( min((0:255), clipMax) * scale );

chunkN = 300; % frames per chunk; tune (200-800 usually good)

for k = 1:numel(files)
    inPath = fullfile(files(k).folder, files(k).name);
    fprintf("[%d/%d] %s\n", k, numel(files), files(k).name);

    vr = VideoReader(inPath);
    tmpPath = fullfile(files(k).folder, ".__tmp__" + string(files(k).name));

    vw = VideoWriter(tmpPath, "Motion JPEG AVI");
    vw.FrameRate = vr.FrameRate;
    open(vw);

    % Best-effort frame count (may be empty/NaN depending on codec)
    try
        nFrames = vr.NumFrames; %#ok<VIDREAD>
    catch
        nFrames = floor(vr.Duration * vr.FrameRate);
    end
    if ~isfinite(nFrames) || nFrames < 1
        nFrames = floor(vr.Duration * vr.FrameRate);
    end

    i = 1;
    while i <= nFrames
        j = min(i + chunkN - 1, nFrames);

        % read chunk: returns HxWxCxN
        F = read(vr, [i j]);

        % ---- ensure grayscale uint8 ----
        if ndims(F) == 4 && size(F,3) == 3
            % fast rgb2gray over 4D: do weighted sum in single then uint8
            R = single(F(:,:,1,:)); G = single(F(:,:,2,:)); B = single(F(:,:,3,:));
            Fg = uint8(0.2989*R + 0.5870*G + 0.1140*B);
        else
            Fg = F;
            if ~isa(Fg, "uint8")
                Fg = im2uint8(Fg);
            end
        end

        % ---- LUT mapping (vectorized) ----
        % lut is 1x256, MATLAB uses 1-based indexing
        out = reshape(lut(double(Fg)+1), size(Fg));

        % ---- write chunk ----
        % writeVideo expects HxWx1xN or HxWxN? It accepts HxWxN for grayscale in many versions.
        % To be safe: write frame-by-frame within the chunk (still much faster than readFrame loop)
        nC = size(out, 4);
        for t = 1:nC
            writeVideo(vw, out(:,:,:,t));
        end

        i = j + 1;
    end

    close(vw);
    vr = []; clear vr;
    pause(0.05);

    % Windows-safe overwrite
    ok = false;
    for attempt = 1:25
        try
            if isfile(inPath), delete(inPath); end
            movefile(tmpPath, inPath, "f");
            ok = true;
            break;
        catch ME
            pause(0.1);
            if attempt == 25
                fprintf(2, "FAILED to overwrite: %s\nTemp kept at: %s\n", inPath, tmpPath);
                rethrow(ME);
            end
        end
    end

    if ok, fprintf("  overwritten OK\n"); end
end

fprintf("Done. All AVIs clipped at %d and rescaled in-place.\n", clipMax);
end
