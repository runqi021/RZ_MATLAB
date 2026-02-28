%%
% make_selROI_overlay_and_stackDFF.m
% ------------------------------------------------------------
% Outputs (ALL keep original movie intensity units; NO 0-255 scaling):
%   1) <movie>_SELROI_overlay_uint16.tif   : uint16 grayscale multipage TIFF with WHITE ROI outlines + ROI numbers
%   2) <movie>_SELROI_proj_uint16.tif      : uint16 grayscale projection (mean) with WHITE outlines + ROI numbers
%   3) <movie>_SELROI_mask_relabel.tif     : uint16 label mask where selected ROIs are relabeled 1..Nsel
%   4) <movie>_SELROI_stackDFF.png         : stacked dFF plot for selected ROIs, labeled 1..Nsel,
%                                           traces scaled by gain, wider spacing, 20% dFF scale bar
%
% Relabel rule:
%   original CP label == n_sel(i)  --->  new label == i
% So mask label i matches dFF_select(:, i)

clear; clc;

%% ---------------- USER PATHS ----------------
SAM_mat  = "E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\fov2_-200_L_30lp_1100nm_5x_00003_MC_cpSAM_output.mat";
movieTif = "E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\fov2_-200_L_30lp_1100nm_5x_00003_MC.tif";

maskTif  = "E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\fov2_-200_L_30lp_1100nm_5x_00003_MC_AVG_for_CP_cp_masks.tif";
maskNpy  = "E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\fov2_-200_L_30lp_1100nm_5x_00003_MC_AVG_for_CP_seg.npy"; % optional fallback

% selected ROI IDs (ORIGINAL labels in CP mask + columns in SAM F)
n_sel = [2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25, 26, 27];

% acquisition
fps = 30;
tossFrames = 30;     % toss first second (30 frames @ 30 Hz)

% outline thickness (pixels)
edgeThickness = 1;   % 1 or 2

% ---------------- ROI LABEL SETTINGS (ADJUST HERE) ----------------
labelFontName  = "Arial";
labelFontSize  = 16;         % <--- change this for label size
labelPadPx     = 8;          % <--- distance OUTSIDE ROI bbox
labelSide      = "right";    % <--- "right" | "left" | "top" | "bottom"
labelThickenPx = 0;          % 0 = no thickening. (keep 0 if you hate bold)

% ---------------- dFF PLOT SETTINGS (ADJUST HERE) ----------------
traceGain    = 2.5;  % <--- scale each trace amplitude by this factor (visual only)
dffScaleBar  = 0.20; % 20% dFF (scale bar text shows 20%; bar height follows traceGain)
minOffset    = 1.00; % minimum spacing in plotted units
offsetMult   = 12;   % spacing = max(minOffset, offsetMult * robustScale)

%% ---------------- OUTPUTS ----------------
[outDir, baseName] = fileparts(movieTif);
outOverlayTif   = fullfile(outDir, baseName + "_SELROI_overlay_uint16.tif");
outProjTif      = fullfile(outDir, baseName + "_SELROI_proj_uint16.tif");
outRelabelMask  = fullfile(outDir, baseName + "_SELROI_mask_relabel.tif");
outStackPng     = fullfile(outDir, baseName + "_SELROI_stackDFF.png");

%% ---------------- LOAD SAM F and compute dFF ----------------
S = load(SAM_mat);
assert(isfield(S, "F"), "SAM mat does not contain field F.");
F_raw = S.F;

F = F_raw;
if tossFrames > 0
    assert(size(F,1) > tossFrames, "Not enough frames to toss.");
    F(1:tossFrames, :) = [];
end

dFFout = helper.dFF_RZ(F);
dFF = dFFout.dFF;

% select columns by original ROI ids -> columns become 1..Nsel
dFF_select = dFF(:, n_sel);
nSel = numel(n_sel);

%% ---------------- LOAD CP LABEL MASK ----------------
L = load_label_mask(maskTif, maskNpy);  % uint16 2D label image

%% ---------------- RELABEL MASK: keep only n_sel and map to 1..Nsel ----------------
Lsel = zeros(size(L), "uint16");
for i = 1:nSel
    Lsel(L == uint16(n_sel(i))) = uint16(i);
end

% warn if some selected labels are absent in mask
counts = zeros(nSel,1);
for i=1:nSel, counts(i) = nnz(Lsel==i); end
missing = find(counts==0);
if ~isempty(missing)
    warning("Some selected ROIs had 0 pixels in the mask after relabeling (new IDs): %s", mat2str(missing'));
end

imwrite(Lsel, outRelabelMask, "tif", "Compression","lzw");

%% ---------------- CHECK MOVIE SIZE + FRAME RANGE ----------------
infoV = imfinfo(movieTif);
nFramesMovie = numel(infoV);
H = infoV(1).Height; W = infoV(1).Width;

assert(all(size(Lsel) == [H W]), ...
    "Mask size (%dx%d) does not match movie (%dx%d).", size(Lsel,1), size(Lsel,2), H, W);

% match video length to dFF length (after toss)
nT = size(dFF_select, 1);
startFrame = tossFrames + 1;
endFrame = min(nFramesMovie, startFrame + nT - 1);
if endFrame < (startFrame + nT - 1)
    warning("Movie shorter than dFF length. Using frames %d..%d (movie end).", startFrame, endFrame);
end
nOutFrames = endFrame - startFrame + 1;

%% ---------------- PRECOMPUTE OUTLINES + ROI NUMBER TEXT MASK ----------------
edgeAll = false(H,W);

for i = 1:nSel
    bw = (Lsel == i);
    if any(bw(:))
        e  = bwperim(bw, 8);
        if edgeThickness > 1
            e = imdilate(e, strel("disk", edgeThickness-1, 0));
        end
        edgeAll = edgeAll | e;
    end
end

% --- label positions OUTSIDE each ROI (bbox-based, top-left-ish anchor) ---
labelXY = nan(nSel,2); % [x,y] in 1-based image coords
for i = 1:nSel
    bw = (Lsel == i);
    if ~any(bw(:)), continue; end

    rp = regionprops(bw, "BoundingBox");
    bb = rp(1).BoundingBox; % [x y w h]

    switch lower(string(labelSide))
        case "right"
            x = bb(1) + bb(3) + labelPadPx;
            y = bb(2);
        case "left"
            x = bb(1) - labelPadPx;   % (we'll shift by string width inside renderer)
            y = bb(2);
        case "top"
            x = bb(1);
            y = bb(2) - labelPadPx;
        case "bottom"
            x = bb(1);
            y = bb(2) + bb(4) + labelPadPx;
        otherwise
            error("labelSide must be right/left/top/bottom");
    end

    labelXY(i,:) = [x y];
end

labelStr = arrayfun(@(k) sprintf("%d", k), 1:nSel, "UniformOutput", false);

% render Arial, non-bold labels into a logical mask
textMask = make_text_mask_java(H, W, labelXY, labelStr, labelFontName, labelFontSize, labelSide);

if labelThickenPx > 0
    textMask = imdilate(textMask, strel("disk", labelThickenPx, 0));
end

markAll = edgeAll | textMask;

%% ---------------- SAVE PROJECTION (uint16, original units) ----------------
projN = min(200, nOutFrames);
proj = compute_mean_frame_tiff_uint16(movieTif, startFrame, startFrame + projN - 1);

projOverlay = proj;
projOverlay(markAll) = max(proj(:));   % WHITE within original unit range (no scaling)
imwrite(projOverlay, outProjTif, "tif", "Compression","lzw");

%% ---------------- WRITE OVERLAY MOVIE (uint16 grayscale, original units) ----------------
if exist(outOverlayTif, "file"), delete(outOverlayTif); end

tIn = Tiff(movieTif, "r");
fprintf("Writing overlay TIFF (uint16 grayscale): %s\n", outOverlayTif);

for k = 1:nOutFrames
    dirIdx = startFrame + (k-1);
    tIn.setDirectory(dirIdx);
    frame = tIn.read();
    if ~isa(frame, "uint16"), frame = uint16(frame); end

    out = frame;
    out(markAll) = max(frame(:));      % WHITE outline + labels in this frame's own unit range

    if k == 1
        imwrite(out, outOverlayTif, "tif", "Compression","lzw");
    else
        imwrite(out, outOverlayTif, "tif", "WriteMode","append", "Compression","lzw");
    end

    if mod(k, 200) == 0
        fprintf("  %d / %d frames\n", k, nOutFrames);
    end
end
tIn.close();

%% ---------------- STACKED dFF PLOT (scaled traces + wider spacing + 20% bar) ----------------
figS = figure("Color","w");
stackDFF_relabel_with_scalebar(dFF_select, fps, dffScaleBar, traceGain, minOffset, offsetMult);
title(sprintf("stackDFF (selected ROIs relabeled 1..%d) | gain=%.2g | orig labels: %s", ...
    nSel, traceGain, mat2str(n_sel)));
exportgraphics(figS, outStackPng, "Resolution", 200);
close(figS);

fprintf("\nDONE.\n- Overlay movie:   %s\n- Projection:      %s\n- Relabeled mask:  %s\n- stackDFF plot:   %s\n", ...
    outOverlayTif, outProjTif, outRelabelMask, outStackPng);

%% ======================== LOCAL FUNCTIONS ================================

function L = load_label_mask(maskTif, maskNpy)
    if isfile(maskTif)
        infoM = imfinfo(maskTif);
        if numel(infoM) > 1
            warning("Mask TIF has %d pages; using FIRST page.", numel(infoM));
        end
        L = imread(maskTif, 1);
    elseif isfile(maskNpy)
        assert(exist("readNPY","file")==2, ...
            "maskNpy provided but readNPY not found. Use maskTif or install npy-matlab (readNPY).");
        L = readNPY(maskNpy);
    else
        error("No mask found. Checked:\n%s\n%s", maskTif, maskNpy);
    end

    if ndims(L) > 2, L = L(:,:,1); end
    if ~isa(L, "uint16"), L = uint16(L); end
end

function meanU16 = compute_mean_frame_tiff_uint16(tifPath, firstIdx, lastIdx)
    info = imfinfo(tifPath);
    lastIdx = min(lastIdx, numel(info));
    H = info(1).Height; W = info(1).Width;

    t = Tiff(tifPath, "r");
    acc = zeros(H, W, "double");
    n = 0;

    for ii = firstIdx:lastIdx
        t.setDirectory(ii);
        fr = t.read();
        acc = acc + double(fr);
        n = n + 1;
    end
    t.close();

    m = round(acc / max(n,1));  % original units, rounded
    m(m < 0) = 0;
    m(m > double(intmax("uint16"))) = double(intmax("uint16"));
    meanU16 = uint16(m);
end

function stackDFF_relabel_with_scalebar(dFF, fps, barAmp, gain, minOffset, offsetMult)
    [T, N] = size(dFF);
    t = (0:T-1) / fps;

    % baseline removal per trace
    d = dFF;
    for i=1:N
        d(:,i) = d(:,i) - median(d(:,i), "omitnan");
    end

    % apply gain (visual)
    d = gain * d;

    % spacing
    s = mad(d(:), 1);
    if ~isfinite(s) || s==0
        s = std(d(:), "omitnan");
    end
    if ~isfinite(s) || s==0
        s = 0.05 * gain;
    end
    offset = max(minOffset, offsetMult * s);

    hold on;
    for i=1:N
        y0 = (i-1)*offset;
        plot(t, d(:,i) + y0, "k", "LineWidth", 0.8);
        text(t(1), y0, sprintf("%d", i), ...
            "Color","k","FontSize",9, "FontName","Arial","FontWeight","normal", ...
            "HorizontalAlignment","left", "VerticalAlignment","bottom");
    end

    yticks((0:N-1)*offset);
    yticklabels(string(1:N));
    xlabel("Time (s)");
    ylabel("ROI (relabelled)");
    box off;

    % scale bar (20% dFF, but plotted height scales with gain)
    barH = barAmp * gain;
    x0 = t(1) + 0.92*(t(end)-t(1));
    y0 = 0.10*offset;
    plot([x0 x0], [y0 y0+barH], "k", "LineWidth", 2);
    text(x0, y0+barH/2, sprintf("  %d%% dFF", round(barAmp*100)), ...
        "Color","k","FontSize",10, "FontName","Arial","FontWeight","normal", ...
        "HorizontalAlignment","left", "VerticalAlignment","middle");
end

function textMask = make_text_mask_java(H, W, xyTopLeft, labelStr, fontName, fontSize, labelSide)
% make_text_mask_java
% Renders text into a BYTE grayscale Java BufferedImage, returns logical mask.
% Fixes the Color(...) range error by using Color.white / Color.black.

    import java.awt.image.BufferedImage
    import java.awt.Color
    import java.awt.Font
    import java.awt.RenderingHints

    fontSize = int32(round(fontSize));

    bi = BufferedImage(W, H, BufferedImage.TYPE_BYTE_GRAY);
    g  = bi.createGraphics();

    % background black
    g.setColor(Color.black);
    g.fillRect(0,0,W,H);

    % font: Arial, PLAIN (non-bold)
    g.setFont(Font(char(fontName), Font.PLAIN, fontSize));
    g.setColor(Color.white);

    % no anti-alias (thin but clean)
    g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,      RenderingHints.VALUE_ANTIALIAS_OFF);

    fm = g.getFontMetrics();
    ascent = int32(fm.getAscent());
    height = int32(fm.getHeight());

    labelSide = lower(string(labelSide));

    for i = 1:numel(labelStr)
        xy = xyTopLeft(i,:);
        if any(~isfinite(xy)), continue; end
        s = char(labelStr{i});
        sw = int32(fm.stringWidth(s));

        % MATLAB 1-based -> Java 0-based
        x = int32(round(xy(1) - 1));
        yTop = int32(round(xy(2) - 1));

        % if left-side labeling, shift by string width so it sits left of bbox
        if labelSide == "left"
            x = x - sw;
        end

        % clamp to image bounds
        if x < 0, x = 0; end
        if yTop < 0, yTop = 0; end
        if x + sw >= W, x = max(int32(0), int32(W) - sw - 1); end
        if yTop + height >= H, yTop = max(int32(0), int32(H) - height - 1); end

        % drawString uses BASELINE y
        y = yTop + ascent;
        g.drawString(s, x, y);
    end

    g.dispose();

    % Extract bytes
    r  = bi.getRaster();
    db = r.getDataBuffer();
    data = int8(db.getData());         % int8[]
    u8 = typecast(data(:), 'uint8');   % 0..255
    img = reshape(u8, [W, H])';        % -> [H W]

    textMask = img > 0;
end
