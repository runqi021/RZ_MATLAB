function [fps, scan_meta] = detect_session_fps(folderPath, fallback_fps)
% DETECT_SESSION_FPS  Get imaging frame rate for a session folder.
%
%   [fps, scan_meta] = detect_session_fps(folderPath)
%   [fps, scan_meta] = detect_session_fps(folderPath, fallback_fps)
%
%   Priority chain (TIFF is ground truth):
%     1. *_meta.mat with source='tif_imfinfo' or 'qc' -> trust it
%     2. *.tif with ScanImage header -> parse, save/overwrite _meta.mat
%     3. *_dFF.mat with params.fps   -> last resort before fallback
%     4. Nothing found               -> return fallback_fps (default 30)
%
%   Any existing _meta.mat from non-TIFF sources (e.g. 'dFF.mat') is
%   overwritten when a TIFF with valid metadata is found.

if nargin < 2 || isempty(fallback_fps)
    fallback_fps = 30;
end
assert(isfolder(folderPath), 'Folder does not exist: %s', folderPath);

% =========================================================================
%  Step 1 — check for existing *_meta.mat from TIFF/QC (trusted sources)
% =========================================================================
meta_hits = dir(fullfile(folderPath, '*_meta.mat'));
if ~isempty(meta_hits)
    meta_path = fullfile(meta_hits(1).folder, meta_hits(1).name);
    scan_meta = load(meta_path);
    if isfield(scan_meta, 'fps') && isnumeric(scan_meta.fps) && scan_meta.fps > 0
        % Trust if source is TIFF/QC, or if saved by dffQC pipeline (has source_tif)
        if (isfield(scan_meta, 'source') && any(strcmp(scan_meta.source, {'tif_imfinfo', 'qc'}))) ...
                || isfield(scan_meta, 'source_tif')
            fps = scan_meta.fps;
            return
        end
        % Otherwise: untrusted source — fall through
        fprintf('[detect_fps] Found _meta.mat but untrusted — verifying against TIFF.\n');
    end
end

% =========================================================================
%  Step 2 — read ScanImage metadata from TIFFs (ground truth)
% =========================================================================
tif_hits = dir(fullfile(folderPath, '*.tif'));
if ~isempty(tif_hits)
    % Prefer raw TIFFs (no _preproc, _MC, _AVG suffix)
    names = {tif_hits.name};
    is_raw = cellfun(@(n) ~contains(n, {'_preproc','_MC','_AVG'},'IgnoreCase',true), names);
    order = [find(is_raw), find(~is_raw)];
    tif_hits = tif_hits(order);

    for tt = 1:numel(tif_hits)
        tif_path = fullfile(tif_hits(tt).folder, tif_hits(tt).name);
        try
            info = imfinfo(tif_path);
            meta_str = build_meta_string(info(1));
            if strlength(strtrim(meta_str)) == 0, continue; end

            raw_fps = parse_scanimage_scalar(meta_str, "SI.hRoiManager.scanFrameRate");
            if isempty(raw_fps) || ~isfinite(raw_fps) || raw_fps <= 0
                period = parse_scanimage_scalar(meta_str, "SI.hRoiManager.scanFramePeriod");
                if ~isempty(period) && isfinite(period) && period > 0
                    raw_fps = 1 / period;
                end
            end

            if ~isempty(raw_fps) && isfinite(raw_fps) && raw_fps > 0
                fps = round(raw_fps);
                scan_meta = build_meta_struct(meta_str, fps, raw_fps, tif_path);
                scan_meta.source = 'tif_imfinfo';

                % Save (or overwrite untrusted) _meta.mat
                [~, tif_stem] = fileparts(tif_hits(tt).name);
                cache_path = fullfile(folderPath, [tif_stem, '_meta.mat']);
                save(cache_path, '-struct', 'scan_meta');
                fprintf('[detect_fps] Saved %s from TIFF (fps=%d)\n', cache_path, fps);

                % Also delete any other _meta.mat that had wrong source
                for mm = 1:numel(meta_hits)
                    old_path = fullfile(meta_hits(mm).folder, meta_hits(mm).name);
                    if ~strcmp(old_path, cache_path) && isfile(old_path)
                        delete(old_path);
                        fprintf('[detect_fps] Removed stale %s\n', old_path);
                    end
                end
                return
            end
        catch ME
            fprintf(2, '[detect_fps] WARNING: failed to read TIFF %s: %s\n', ...
                tif_path, ME.message);
        end
    end
end

% =========================================================================
%  Step 3 — fall back to *_dFF.mat params.fps (not fully trusted)
% =========================================================================
dff_hits = dir(fullfile(folderPath, '*_dFF.mat'));
if ~isempty(dff_hits)
    dff_path = fullfile(dff_hits(1).folder, dff_hits(1).name);
    tmp = load(dff_path, 'params');
    if isfield(tmp, 'params') && isfield(tmp.params, 'fps') ...
            && isnumeric(tmp.params.fps) && tmp.params.fps > 0
        fps = tmp.params.fps;
        scan_meta = struct('fps', fps, 'source', 'dFF.mat');
        fprintf(2, '[detect_fps] WARNING: using _dFF.mat fps=%d (no TIFF metadata found). May be wrong for old QC output.\n', fps);
        return
    end
end

% =========================================================================
%  Step 4 — fallback
% =========================================================================
fps = fallback_fps;
scan_meta = struct('fps', fps, 'source', 'fallback');
fprintf(2, '[detect_fps] WARNING: no metadata found in %s, using fallback fps=%d\n', ...
    folderPath, fps);

end  % detect_session_fps


% =========================================================================
%  Local helpers
% =========================================================================

function scan_meta = build_meta_struct(meta_str, fps, raw_fps, tif_path)
% Build full metadata struct matching Batch_dffQC_260325.m output.
PixelSizeBase = 1.7778;

zoomFactor = parse_scanimage_scalar(meta_str, "SI.hRoiManager.scanZoomFactor");
if ~isempty(zoomFactor) && isfinite(zoomFactor) && zoomFactor > 0
    pixelSize_um = PixelSizeBase / zoomFactor;
else
    pixelSize_um = NaN;
end

scan_meta = struct();
scan_meta.fps                = fps;
scan_meta.scanFrameRate_raw  = raw_fps;
scan_meta.zoomFactor         = zoomFactor;
scan_meta.pixelSize_um       = pixelSize_um;
scan_meta.laserPower_pct     = parse_scanimage_vector(meta_str, "SI.hBeams.powers");
scan_meta.channelSave        = parse_scanimage_vector(meta_str, "SI.hChannels.channelSave");
scan_meta.channelInputRanges = parse_scanimage_vector(meta_str, "SI.hChannels.channelInputRanges");
scan_meta.motorPosition      = parse_scanimage_vector(meta_str, "SI.hMotors.motorPosition");
scan_meta.framesPerSlice     = parse_scanimage_scalar(meta_str, "SI.hStackManager.framesPerSlice");
scan_meta.numSlices          = parse_scanimage_scalar(meta_str, "SI.hStackManager.numSlices");
scan_meta.source_tif         = tif_path;
end

function meta_str = build_meta_string(info1)
meta_str = "";
if isfield(info1, 'Software') && ~isempty(info1.Software)
    meta_str = meta_str + string(info1.Software) + newline;
end
if isfield(info1, 'ImageDescription') && ~isempty(info1.ImageDescription)
    meta_str = meta_str + string(info1.ImageDescription) + newline;
end
end

function v = parse_scanimage_scalar(metaStr, key)
v = [];
pat = key + "\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
vv = str2double(string(tok{1}));
if isfinite(vv), v = vv; end
end

function vec = parse_scanimage_vector(metaStr, key)
vec = [];
pat = key + "\s*=\s*([^\r\n]+)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
rhs = strtrim(string(tok{1}));
nums = regexp(rhs, "[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", "match");
if ~isempty(nums)
    vec = str2double(string(nums));
    vec = vec(isfinite(vec));
end
end
