function generate_perROI_videos(folderPath, samPath, patchHalf_um)
%GENERATE_PERROI_VIDEOS  Per-ROI MP4: movie patch + F trace with time indicator.
%
%   generate_perROI_videos(folderPath, samPath, patchHalf_um)
%
%   Outputs (in folderPath/Fmovie_perROI/):
%     ROI_01_movie.mp4, ROI_02_movie.mp4, ...
%
%   Each video: [movie patch (left) | F trace with yellow time line (right)]
%   Black background, white trace, yellow ROI perimeter + time indicator.

if nargin < 3, patchHalf_um = 15; end
videoQuality = 98;

%% Detect FPS + pixel size
[fps, scan_meta] = detect_session_fps(folderPath);
fprintf('  FPS = %d Hz\n', fps);

assert(isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um), ...
    'detect_session_fps did not return a valid pixelSize_um.');
pixelSize_um = scan_meta.pixelSize_um;
patchHalf_px = round(patchHalf_um / pixelSize_um);

%% Load cpSAM + compute dFF
samData = load(samPath, 'maskL', 'F');
assert(isfield(samData, 'maskL'), 'cpSAM file missing maskL.');
assert(isfield(samData, 'F'), 'cpSAM file missing F.');
maskL = samData.maskL;
F_all = double(samData.F);
[T_F, N_roi] = size(F_all);
fprintf('  F: [%d x %d], %d ROIs\n', T_F, N_roi, N_roi);

if N_roi == 0, fprintf('  No ROIs. Skipping.\n'); return; end

dFFout = helper.dFF_RZ(F_all, 'FPS', fps);
dFF_all = dFFout.dFF;  % [T x N]

%% Load MC TIFF
tifPath = discover_mc_tiff(folderPath);
fprintf('  MC TIFF: %s\n', tifPath);
fprintf('  Loading TIFF stack...\n');
tic;
tifInfo = imfinfo(tifPath);
H = tifInfo(1).Height;
W = tifInfo(1).Width;
T_tif = numel(tifInfo);

Y = zeros(H, W, T_tif, 'single');
for kk = 1:T_tif
    Y(:,:,kk) = single(imread(tifPath, kk, 'Info', tifInfo));
end
fprintf('  Loaded %d frames [%d x %d] in %.1f s\n', T_tif, H, W, toc);

assert(T_tif == T_F, 'Frame mismatch: TIFF=%d, F=%d.', T_tif, T_F);

%% Setup
props = regionprops(maskL, 'Centroid');
movieDir = fullfile(folderPath, 'Fmovie_perROI');
if ~exist(movieDir, 'dir'), mkdir(movieDir); end

fprintf('  Rendering %d ROI videos...\n', N_roi);

for roiIdx = 1:N_roi
    outMP4 = fullfile(movieDir, sprintf('ROI_%02d_movie.mp4', roiIdx));
    if isfile(outMP4)
        fprintf('  [SKIP] ROI %d (exists)\n', roiIdx);
        continue;
    end

    cx = round(props(roiIdx).Centroid(1));
    cy = round(props(roiIdx).Centroid(2));
    r1 = max(1, cy - patchHalf_px);  r2 = min(H, cy + patchHalf_px);
    c1 = max(1, cx - patchHalf_px);  c2 = min(W, cx + patchHalf_px);
    pH = r2 - r1 + 1;  pW = c2 - c1 + 1;

    patch_all = Y(r1:r2, c1:c2, :);

    mask_local = (maskL(r1:r2, c1:c2) == roiIdx);

    % Contrast from ROI-masked pixels only — captures the cell's own dynamic range
    roi_vals = patch_all(repmat(mask_local, [1 1 T_tif]));
    lo = double(prctile(roi_vals, 30));
    hi = double(prctile(roi_vals, 99.9));
    if hi <= lo, hi = lo + 1; end
    perim_local = bwperim(mask_local, 4);
    dFF_roi = dFF_all(:, roiIdx);
    t_sec = (0:T_tif-1)' / fps;

    % Upscale movie patch to ~300px
    movieSz = 300;
    scaleFactor = movieSz / max(pH, pW);
    perim_up = imresize(double(perim_local), scaleFactor, 'nearest') > 0.5;
    upH = size(perim_up, 1);  upW = size(perim_up, 2);

    % Trace panel: same height, 4x wider + right margin for scale bar
    scaleMargin = 70;  % px reserved on right for scale bar
    traceH = upH;
    traceW = upW * 4 + scaleMargin;
    traceDrawW = traceW - scaleMargin;  % trace only drawn in this region

    % Pre-render dFF trace (white on black, thick line)
    traceBG = zeros(traceH, traceW, 3, 'uint8');
    dFF_min = min(dFF_roi);  dFF_max = max(dFF_roi);
    if dFF_max == dFF_min, dFF_max = dFF_min + 1; end
    dFF_range = dFF_max - dFF_min;
    yMargin = round(traceH * 0.08);
    yLo = traceH - yMargin;  yHi = yMargin + 1;
    plotH = yLo - yHi;
    x_px = linspace(1, traceDrawW, T_tif);
    y_px = round(yLo - (dFF_roi - dFF_min) / dFF_range * plotH);
    y_px = max(1, min(traceH, y_px));

    lineThick = 2;  % half-thickness ±2px
    for tt = 1:T_tif-1
        x0 = round(x_px(tt));  x1 = round(x_px(tt+1));
        y0 = y_px(tt);  y1 = y_px(tt+1);
        nSteps = max(abs(x1-x0), abs(y1-y0));
        if nSteps == 0, nSteps = 1; end
        for s = 0:nSteps
            xi = max(1, min(traceDrawW, round(x0 + (x1-x0)*s/nSteps)));
            yi = round(y0 + (y1-y0)*s/nSteps);
            for dy = -lineThick:lineThick
                yy = max(1, min(traceH, yi + dy));
                traceBG(yy, xi, :) = 255;
            end
        end
    end

    % 20% dFF scale bar (white, in right margin area, vertically centered)
    scaleBar_dFF = 0.20;
    barH_px = round(scaleBar_dFF / dFF_range * plotH);
    barH_px = max(barH_px, 5);
    barX = traceDrawW + round(scaleMargin / 2);
    barYcenter = round((yHi + yLo) / 2);
    barY1 = barYcenter - round(barH_px / 2);
    barY2 = barY1 + barH_px;
    barY1 = max(1, barY1);  barY2 = min(traceH, barY2);
    % Vertical bar (3px wide)
    for bx = max(1, barX-1):min(traceW, barX+1)
        traceBG(barY1:barY2, bx, :) = 255;
    end
    % Label "20%" above bar
    scaleFont = make_digit_font();
    scaleGlyph = render_text('0.2', scaleFont);
    scaleGlyph = imresize(double(scaleGlyph), 3, 'nearest') > 0.5;
    sgh = size(scaleGlyph, 1);  sgw = size(scaleGlyph, 2);
    sly = max(1, barY1 - sgh - 3);
    slx = max(1, barX - round(sgw / 2));
    if sly + sgh - 1 <= traceH && slx + sgw - 1 <= traceW
        for ch = 1:3
            plane = traceBG(sly:sly+sgh-1, slx:slx+sgw-1, ch);
            plane(scaleGlyph) = 255;
            traceBG(sly:sly+sgh-1, slx:slx+sgw-1, ch) = plane;
        end
    end

    % Pre-render timestamp glyphs (bitmap font)
    digitFont = make_digit_font();

    % Video dimensions
    gapW = 4;
    totalW = upW + gapW + traceW;  totalH = upH;
    if mod(totalW, 2) ~= 0, totalW = totalW + 1; end
    if mod(totalH, 2) ~= 0, totalH = totalH + 1; end

    vw = VideoWriter(outMP4, 'MPEG-4');
    vw.FrameRate = fps;  vw.Quality = videoQuality;
    open(vw);

    for ff = 1:T_tif
        ptch = double(patch_all(:,:,ff));
        ptch = max(0, min(1, (ptch - lo) / (hi - lo)));
        ptch_up = imresize(ptch, [upH upW], 'bilinear');
        ptch_u8 = uint8(255 * ptch_up);
        movieRGB = repmat(ptch_u8, [1 1 3]);

        movieRGB(:,:,1) = movieRGB(:,:,1) .* uint8(~perim_up) + uint8(perim_up)*255;
        movieRGB(:,:,2) = movieRGB(:,:,2) .* uint8(~perim_up) + uint8(perim_up)*255;
        movieRGB(:,:,3) = movieRGB(:,:,3) .* uint8(~perim_up);

        traceFrame = traceBG;

        % Burn timestamp into trace panel (top-left, white, 2x size)
        tStr = sprintf('%.1fs', t_sec(ff));
        tGlyph = render_text(tStr, digitFont);
        tGlyph = imresize(double(tGlyph), 3, 'nearest') > 0.5;
        gh = size(tGlyph, 1);  gw = size(tGlyph, 2);
        ty = 3;  tx = 5;
        if tx + gw - 1 <= traceW && ty + gh - 1 <= traceH
            for ch = 1:3
                plane = traceFrame(ty:ty+gh-1, tx:tx+gw-1, ch);
                plane(tGlyph) = 255;
                traceFrame(ty:ty+gh-1, tx:tx+gw-1, ch) = plane;
            end
        end
        xLine = max(2, min(traceW-1, round(x_px(ff))));
        traceFrame(:, xLine-1:xLine+1, 1) = 255;
        traceFrame(:, xLine-1:xLine+1, 2) = 255;
        traceFrame(:, xLine-1:xLine+1, 3) = 0;

        fullFrame = zeros(totalH, totalW, 3, 'uint8');
        fullFrame(1:upH, 1:upW, :) = movieRGB;
        fullFrame(1:traceH, upW+gapW+1:upW+gapW+traceW, :) = traceFrame;
        writeVideo(vw, fullFrame);
    end
    close(vw);
    fprintf('  ROI %d/%d saved\n', roiIdx, N_roi);
end
fprintf('  All %d ROI videos done.\n', N_roi);
end


function tifPath = discover_mc_tiff(fovFolder)
    mc_hits = dir(fullfile(fovFolder, '*_MC*.tif'));
    assert(~isempty(mc_hits), 'No *_MC*.tif found in %s', fovFolder);
    names = {mc_hits.name};
    keep = true(size(names));
    for ii = 1:numel(names)
        if contains(names{ii}, '_AVG', 'IgnoreCase', true) || ...
           contains(names{ii}, '_QC',  'IgnoreCase', true) || ...
           contains(names{ii}, '_shifts', 'IgnoreCase', true) || ...
           contains(names{ii}, '_SELROI', 'IgnoreCase', true)
            keep(ii) = false;
        end
    end
    mc_hits = mc_hits(keep);
    assert(~isempty(mc_hits), 'No valid MC TIFF after filtering in %s', fovFolder);
    nMC = zeros(numel(mc_hits), 1);
    for ii = 1:numel(mc_hits)
        nMC(ii) = numel(strfind(upper(mc_hits(ii).name), '_MC'));
    end
    [~, bestIdx] = max(nMC);
    tifPath = fullfile(mc_hits(bestIdx).folder, mc_hits(bestIdx).name);
end


function digitFont = make_digit_font()
%MAKE_DIGIT_FONT  7x5 bitmap font for digits 0-9, '.', 's'
    digitFont = struct();
    digitFont.chars = '0123456789.s%';
    g = zeros(7, 5, 13);
    g(:,:,1)  = [0 1 1 1 0; 1 0 0 0 1; 1 0 0 1 1; 1 0 1 0 1; 1 1 0 0 1; 1 0 0 0 1; 0 1 1 1 0]; % 0
    g(:,:,2)  = [0 0 1 0 0; 0 1 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 1 1 1 0]; % 1
    g(:,:,3)  = [0 1 1 1 0; 1 0 0 0 1; 0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0; 1 1 1 1 1]; % 2
    g(:,:,4)  = [0 1 1 1 0; 1 0 0 0 1; 0 0 0 0 1; 0 0 1 1 0; 0 0 0 0 1; 1 0 0 0 1; 0 1 1 1 0]; % 3
    g(:,:,5)  = [0 0 0 1 0; 0 0 1 1 0; 0 1 0 1 0; 1 0 0 1 0; 1 1 1 1 1; 0 0 0 1 0; 0 0 0 1 0]; % 4
    g(:,:,6)  = [1 1 1 1 1; 1 0 0 0 0; 1 1 1 1 0; 0 0 0 0 1; 0 0 0 0 1; 1 0 0 0 1; 0 1 1 1 0]; % 5
    g(:,:,7)  = [0 1 1 1 0; 1 0 0 0 0; 1 0 0 0 0; 1 1 1 1 0; 1 0 0 0 1; 1 0 0 0 1; 0 1 1 1 0]; % 6
    g(:,:,8)  = [1 1 1 1 1; 0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0]; % 7
    g(:,:,9)  = [0 1 1 1 0; 1 0 0 0 1; 1 0 0 0 1; 0 1 1 1 0; 1 0 0 0 1; 1 0 0 0 1; 0 1 1 1 0]; % 8
    g(:,:,10) = [0 1 1 1 0; 1 0 0 0 1; 1 0 0 0 1; 0 1 1 1 1; 0 0 0 0 1; 0 0 0 0 1; 0 1 1 1 0]; % 9
    g(:,:,11) = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0]; % .
    g(:,:,12) = [0 1 1 1 0; 1 0 0 0 0; 0 1 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 0 0 0 1; 0 1 1 1 0]; % s
    g(:,:,13) = [1 0 0 1 0; 0 0 1 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 1 0 0 0; 0 1 0 0 1; 0 0 0 1 0]; % %  (slash with dots)
    digitFont.glyphs = logical(g);
end


function img = render_text(str, digitFont)
%RENDER_TEXT  Render string into a logical bitmap using digitFont.
    img = false(7, 0);
    for k = 1:numel(str)
        idx = find(digitFont.chars == str(k), 1);
        if isempty(idx)
            img = [img, false(7, 3)]; %#ok<AGROW> % space for unknown
        else
            if k > 1, img = [img, false(7, 1)]; end %#ok<AGROW> % 1px gap
            img = [img, digitFont.glyphs(:,:,idx)]; %#ok<AGROW>
        end
    end
end
