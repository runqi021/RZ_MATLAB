function meta = tiffinfo(tiffPath)
%TIFFINFO  Print key ScanImage acquisition parameters from a (raw) TIFF.
%
%   meta = tiffinfo(tiffPath)
%
%   Reads the ScanImage header (info(1).Software / ImageDescription) and
%   prints + returns the most-used acquisition parameters:
%       fps                    frame rate (Hz)
%       motorPosition          [x y z] stage coordinate (um)
%       zoom                   scanZoomFactor
%       pixelSize_um           micron / pixel (FOV/pixels, fallback base/zoom)
%       laserPower_pct         Pockels power(s) (%)
%       channelSave            channels written to disk
%       channelInputRange      per-channel ADC input range (V), the "threshold"
%
%   NOTE: only RAW ScanImage TIFFs carry this metadata. Processed TIFFs
%   (_minusDark / _MC / _AVG) have none -- use the raw file or _meta.mat.
%
%   See also: detect_session_fps (the pipeline's cached metadata reader).

    info = imfinfo(tiffPath);

    meta_str = "";
    if isfield(info(1),'Software') && ~isempty(info(1).Software)
        meta_str = meta_str + string(info(1).Software) + newline;
    end
    if isfield(info(1),'ImageDescription') && ~isempty(info(1).ImageDescription)
        meta_str = meta_str + string(info(1).ImageDescription) + newline;
    end
    assert(strlength(strtrim(meta_str)) > 0, ...
        'No ScanImage metadata (Software/ImageDescription) in:\n%s', tiffPath);

    g  = @(k) get_scalar(meta_str, k);
    gv = @(k) get_vector(meta_str, k);

    meta = struct();
    meta.source_tif = string(tiffPath);

    % --- fps ---
    meta.fps = g("SI.hRoiManager.scanFrameRate");
    if isnan(meta.fps)
        p = g("SI.hRoiManager.scanFramePeriod");
        if ~isnan(p) && p > 0, meta.fps = 1/p; end
    end

    % --- XYZ stage coordinate (um) ---
    meta.motorPosition = gv("SI.hMotors.motorPosition");

    % --- zoom ---
    meta.zoom = g("SI.hRoiManager.scanZoomFactor");

    % --- pixel size (um): prefer FOV(um)/pixels, fall back to base/zoom ---
    fovUm     = gv("SI.hRoiManager.imagingFovUm");      % 4 corners, [x y] flattened
    pxPerLine = g("SI.hRoiManager.pixelsPerLine");
    meta.pixelSize_um = NaN;
    if numel(fovUm) >= 8 && ~isnan(pxPerLine) && pxPerLine > 0
        xs = fovUm(1:2:end);
        meta.pixelSize_um = (max(xs) - min(xs)) / pxPerLine;
    elseif ~isnan(meta.zoom) && meta.zoom > 0
        meta.pixelSize_um = 1.7778 / meta.zoom;         % PixelSizeBase fallback
    end

    % --- laser power (%) ---
    meta.laserPower_pct = gv("SI.hBeams.powers");

    % --- channels + input ranges ("threshold", V) ---
    meta.channelSave       = gv("SI.hChannels.channelSave");
    meta.channelInputRange = get_cell_ranges(meta_str, "SI.hChannels.channelInputRange");

    % ---------------- pretty print ----------------
    [~, name] = fileparts(tiffPath);
    fprintf('\n===== %s =====\n', name);
    fprintf('  fps           : %.4f Hz\n', meta.fps);
    fprintf('  stage XYZ     : [%.2f, %.2f, %.2f] um\n', meta.motorPosition);
    fprintf('  zoom          : %gx\n', meta.zoom);
    fprintf('  pixel size    : %.4f um/px\n', meta.pixelSize_um);
    fprintf('  laser power   : %s %%\n', num2str(meta.laserPower_pct(:)'));
    if isempty(meta.channelSave)
        fprintf('  channels      : (none parsed)\n');
    else
        for c = meta.channelSave(:)'
            if c <= numel(meta.channelInputRange) && ~isempty(meta.channelInputRange{c})
                r = meta.channelInputRange{c};
                fprintf('  ch%-2d range    : [%g, %g] V\n', c, r(1), r(end));
            else
                fprintf('  ch%-2d range    : (not found)\n', c);
            end
        end
    end
    fprintf('\n');
end

% ======================= local helpers =======================
function v = get_scalar(meta, key)
    v = NaN;
    pat = key + "\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)";
    tok = regexp(meta, pat, "tokens", "once");
    if isempty(tok), return; end
    vv = str2double(string(tok{1}));
    if isfinite(vv), v = vv; end
end

function vec = get_vector(meta, key)
    % Parses "key = [ ... ]" (also tolerates a bare scalar) into a row vector.
    vec = [];
    pat = key + "\s*=\s*(\[[^\]]*\]|[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)";
    tok = regexp(meta, pat, "tokens", "once");
    if isempty(tok), return; end
    nums = regexp(string(tok{1}), "[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", "match");
    if isempty(nums), return; end
    vec = str2double(nums);
    vec = vec(isfinite(vec));
end

function ranges = get_cell_ranges(meta, key)
    % Parses "key = {[-0.5 0.5] [-0.5 0.5] ...}" into a cell of numeric vectors.
    ranges = {};
    pat = key + "\s*=\s*\{([^}]*)\}";
    tok = regexp(meta, pat, "tokens", "once");
    if isempty(tok), return; end
    blocks = regexp(string(tok{1}), "\[([^\]]*)\]", "tokens");
    for i = 1:numel(blocks)
        nums = str2double(regexp(blocks{i}{1}, ...
            "[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", "match"));
        ranges{end+1} = nums(isfinite(nums)); %#ok<AGROW>
    end
end
