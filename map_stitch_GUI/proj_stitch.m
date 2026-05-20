close all; clear; clc;
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
% proj_stitch.m
% Ch1 max projection + Ch3 avg projection from stitched Z-stacks.
% Ch1 uses proj gamma, Ch3 uses per-Z gamma.
% Single figure: ch1 left, ch3 right.

%% ======== 1) User parameters — always set ========
stitchFolder = 'C:\Users\Admin\Desktop\chat_260209\matlab_stitch';
repeatFromMeta = false;    % true = load all params from proj_meta.mat, ignore below

z_use_ch1 = 4:21;   % which Z slices to include ([] = all)
z_use_ch3 = [];      % which Z slices to include ([] = all)

pLo_ch1 = 95;      pHi_ch1 = 100;   % ch1 percentile clip
pLo_ch3 = 75;       pHi_ch3 = 96;   % ch3 percentile clip

smooth_ch1 = 0.6;      % ch1 gaussian sigma in pixels ([] = no smoothing)
gamma_ch1 = 0.3;       % ch1 proj gamma ([] = no gamma)
gamma_ch3_perz = 0.3;  % ch3 per-Z gamma ([] = no gamma)
gamma_ch3_proj = 0.5;   % ch3 proj gamma, applied after per-Z ([] = no gamma)

doSave = true;      % true = save TIFFs; false = just plot

scaleBar_um = 500;         % scale bar length (um)
scaleBar_side = "left";   % "right" or "left"
rotateCW = 5;              % clockwise rotation: 0, 90, 180, or 270 degrees
flipLR = false;            % flip horizontally (left-right mirror)
crop_x_um = [0 3600];   % [x_start, x_size] in um from left; [] = no x crop
crop_y_um = [0 3600];    % [y_start, y_size] in um from bottom; [] = no y crop

%% ======== 2) Acquisition parameters — auto-detected if raw TIFF in parent ========
um_per_px_base = 1.7778;   % um/px at zoom=1 (microscope constant, always needed)
zoom = 1.2;                % fallback if no raw TIFF found
range_ch1 = [-1 1];   % fallback channel 1 input range [V]
range_ch3 = [-0.25 0.25]; % fallback channel 3 input range [V]

%% ======== Repeat from saved meta ========
if repeatFromMeta
    metaFile = fullfile(stitchFolder, 'proj_meta.mat');
    assert(isfile(metaFile), 'No proj_meta.mat in %s', stitchFolder);
    m = load(metaFile);
    z_use_ch1 = m.z_use_ch1;  z_use_ch3 = m.z_use_ch3;
    pLo_ch1 = m.pLo_ch1;  pHi_ch1 = m.pHi_ch1;
    pLo_ch3 = m.pLo_ch3;  pHi_ch3 = m.pHi_ch3;
    if isfield(m,'smooth_ch1'), smooth_ch1 = m.smooth_ch1; end
    gamma_ch1 = m.gamma_ch1;
    gamma_ch3_perz = m.gamma_ch3_perz;
    gamma_ch3_proj = m.gamma_ch3_proj;
    zoom = m.zoom;  um_per_px_base = m.um_per_px_base;  um_per_px = m.um_per_px;
    range_ch1 = m.range_ch1;  range_ch3 = m.range_ch3;
    scaleBar_um = m.scaleBar_um;  scaleBar_side = m.scaleBar_side;
    rotateCW = m.rotateCW;  flipLR = m.flipLR;
    crop_x_um = m.crop_x_um;  crop_y_um = m.crop_y_um;
    fprintf('[repeat] Loaded params from %s\n', metaFile);
end

%% ======== Auto-detect acquisition params from raw TIFF in parent ========
if ~repeatFromMeta
    parentFolder = fileparts(stitchFolder);
    [~, scan_meta] = detect_session_fps(parentFolder);
    if isfield(scan_meta, 'source') ...
            && any(strcmp(scan_meta.source, {'tif_imfinfo', 'qc'}))
        if isfield(scan_meta, 'zoomFactor') && ~isnan(scan_meta.zoomFactor)
            zoom = scan_meta.zoomFactor;
            fprintf('[auto] zoom = %.2f  (from TIFF)\n', zoom);
        end
        if isfield(scan_meta, 'channelInputRanges') ...
                && numel(scan_meta.channelInputRanges) >= 6
            ranges = scan_meta.channelInputRanges;
            range_ch1 = ranges(1:2);   % Ch1 = first pair
            range_ch3 = ranges(5:6);   % Ch3 = third pair
            fprintf('[auto] range_ch1=[%.2f %.2f]  range_ch3=[%.2f %.2f]  (from TIFF)\n', ...
                range_ch1(1), range_ch1(2), range_ch3(1), range_ch3(2));
        end
    else
        fprintf('[auto] No raw TIFF metadata in parent folder — using manual defaults.\n');
    end
end

um_per_px = um_per_px_base / zoom;
fprintf('zoom=%.2f  um_per_px=%.4f\n', zoom, um_per_px);

%% ======== Read stitched stacks ========
ch1_file = fullfile(stitchFolder, 'stitched_ch1_avgz.tif');
ch3_file = fullfile(stitchFolder, 'stitched_ch3_avgz.tif');

info1 = imfinfo(ch1_file);
nFrames = numel(info1);
H = info1(1).Height;  W = info1(1).Width;

I0 = imread(ch1_file, 1);
inClass = class(I0);
stack_ch1 = zeros(H, W, nFrames, inClass);
stack_ch3 = zeros(H, W, nFrames, inClass);
stack_ch1(:,:,1) = I0;
for k = 2:nFrames
    stack_ch1(:,:,k) = imread(ch1_file, k);
end
for k = 1:nFrames
    stack_ch3(:,:,k) = imread(ch3_file, k);
end
fprintf('Loaded %d frames (%d x %d) class=%s\n', nFrames, W, H, inClass);

%% ======== Apply z_use ========
if isempty(z_use_ch1), z_use_ch1 = 1:nFrames; end
if isempty(z_use_ch3), z_use_ch3 = 1:nFrames; end
stack_ch1 = stack_ch1(:,:, z_use_ch1);
stack_ch3 = stack_ch3(:,:, z_use_ch3);
fprintf('z_use_ch1 = %s  (%d slices)\n', mat2str(z_use_ch1), numel(z_use_ch1));
fprintf('z_use_ch3 = %s  (%d slices)\n', mat2str(z_use_ch3), numel(z_use_ch3));

%% ======== Compute projections ========
doGamma_ch1 = ~isempty(gamma_ch1);
doGamma_ch3_perz = ~isempty(gamma_ch3_perz);
doGamma_ch3_proj = ~isempty(gamma_ch3_proj);
doGamma_ch3 = doGamma_ch3_perz || doGamma_ch3_proj;

% Ch1: max projection, optional smooth, proj gamma
max_ch1 = max(stack_ch1, [], 3);
if ~isempty(smooth_ch1) && smooth_ch1 > 0
    max_ch1 = cast(imgaussfilt(double(max_ch1), smooth_ch1), inClass);
end
max_ch1_lin = clip_rescale_u16(max_ch1, pLo_ch1, pHi_ch1);
if doGamma_ch1
    max_ch1_gam = apply_gamma_u16(max_ch1_lin, gamma_ch1);
end

% Ch3: avg projection, per-Z gamma then proj gamma
if doGamma_ch3_perz
    stack_ch3_g = apply_gamma_stack(stack_ch3, gamma_ch3_perz);
    avg_ch3_src = cast(mean(single(stack_ch3_g), 3), inClass);
else
    avg_ch3_src = cast(mean(single(stack_ch3), 3), inClass);
end
avg_ch3_lin = clip_rescale_u16(avg_ch3_src, pLo_ch3, pHi_ch3);
if doGamma_ch3_proj
    avg_ch3_gam = apply_gamma_u16(avg_ch3_lin, gamma_ch3_proj);
else
    avg_ch3_gam = avg_ch3_lin;
end

%% ======== Save TIFFs (after rotate) ========
% Save is deferred until after rotate so saved images match display

%% ======== Pick display images ========
if doGamma_ch1
    disp_ch1 = max_ch1_gam;
    lbl_ch1 = sprintf('Ch1 Max  \\gamma=%.1f   clip[%.0f%%, %.1f%%]   z=%s', ...
        gamma_ch1, pLo_ch1, pHi_ch1, mat2str(z_use_ch1));
else
    disp_ch1 = max_ch1_lin;
    lbl_ch1 = sprintf('Ch1 Max   clip[%.0f%%, %.1f%%]   z=%s', ...
        pLo_ch1, pHi_ch1, mat2str(z_use_ch1));
end

if doGamma_ch3
    disp_ch3 = avg_ch3_gam;
    gStr = '';
    if doGamma_ch3_perz, gStr = sprintf('perZ=%.1f', gamma_ch3_perz); end
    if doGamma_ch3_proj
        if ~isempty(gStr), gStr = [gStr ' ']; end
        gStr = [gStr sprintf('proj=%.1f', gamma_ch3_proj)];
    end
    lbl_ch3 = sprintf('Ch3 Avg  \\gamma(%s)   clip[%.0f%%, %.1f%%]   z=%s', ...
        gStr, pLo_ch3, pHi_ch3, mat2str(z_use_ch3));
else
    disp_ch3 = avg_ch3_lin;
    lbl_ch3 = sprintf('Ch3 Avg   clip[%.0f%%, %.1f%%]   z=%s', ...
        pLo_ch3, pHi_ch3, mat2str(z_use_ch3));
end

%% ======== Rotate ========
if rotateCW ~= 0
    nRot90 = round(mod(rotateCW, 360) / 90);  % snap to nearest 90
    residual = mod(rotateCW, 360) - nRot90 * 90;
    if nRot90 > 0
        disp_ch1 = rot90(disp_ch1, -nRot90);
        disp_ch3 = rot90(disp_ch3, -nRot90);
    end
    if abs(residual) > 0.01
        disp_ch1 = uint16(imrotate(double(disp_ch1), -residual, 'bilinear', 'crop'));
        disp_ch3 = uint16(imrotate(double(disp_ch3), -residual, 'bilinear', 'crop'));
    end
end
if flipLR
    disp_ch1 = fliplr(disp_ch1);
    disp_ch3 = fliplr(disp_ch3);
end

%% ======== Crop (origin = bottom-left) ========
[H_d, W_d] = size(disp_ch1);
if ~isempty(crop_x_um)
    x1 = round(crop_x_um(1) / um_per_px) + 1;
    x2 = x1 + round(crop_x_um(2) / um_per_px) - 1;
    x1 = max(1, x1);  x2 = min(W_d, x2);
    disp_ch1 = disp_ch1(:, x1:x2);
    disp_ch3 = disp_ch3(:, x1:x2);
    fprintf('[crop_x] x=[%d:%d] px  (%.0f um)\n', x1, x2, (x2-x1+1)*um_per_px);
end
if ~isempty(crop_y_um)
    [H_d, ~] = size(disp_ch1);
    y_bot = round(crop_y_um(1) / um_per_px);
    y_size = round(crop_y_um(2) / um_per_px);
    y2 = H_d - y_bot;
    y1 = y2 - y_size + 1;
    y1 = max(1, y1);  y2 = min(H_d, y2);
    disp_ch1 = disp_ch1(y1:y2, :);
    disp_ch3 = disp_ch3(y1:y2, :);
    fprintf('[crop_y] y=[%d:%d] px  (%.0f um)\n', y1, y2, (y2-y1+1)*um_per_px);
end

%% ======== Merge (before scale bar) ========
g1 = double(disp_ch1) / 65535;
g3 = double(disp_ch3) / 65535;
merged_gm = cat(3, g3, g1, g3);  % R=ch3, G=ch1, B=ch3 → green + magenta
merged_gr = cat(3, g3, g1, zeros(size(g1)));  % R=ch3, G=ch1, B=0 → green + red

%% ======== Burn scale bar ========
disp_ch1 = burn_scalebar(disp_ch1, scaleBar_um, um_per_px, scaleBar_side);
disp_ch3 = burn_scalebar(disp_ch3, scaleBar_um, um_per_px, scaleBar_side);
merged_gm = burn_scalebar_rgb(merged_gm, scaleBar_um, um_per_px, scaleBar_side);
merged_gr = burn_scalebar_rgb(merged_gr, scaleBar_um, um_per_px, scaleBar_side);

%% ======== Save TIFFs (with pixel calibration for ImageJ) ========
if doSave
    write_calibrated_tiff(fullfile(stitchFolder, 'stitched_ch1_maxproj.tif'), ...
        disp_ch1, um_per_px);
    write_calibrated_tiff(fullfile(stitchFolder, 'stitched_ch3_avgproj.tif'), ...
        disp_ch3, um_per_px);
    % Color composites as PNG (uint16 RGB)
    imwrite(uint16(merged_gm * 65535), fullfile(stitchFolder, 'stitched_green_magenta.png'));
    imwrite(uint16(merged_gr * 65535), fullfile(stitchFolder, 'stitched_green_red.png'));
    % Save all parameters
    proj_meta = struct();
    proj_meta.stitchFolder = stitchFolder;
    proj_meta.z_use_ch1 = z_use_ch1;
    proj_meta.z_use_ch3 = z_use_ch3;
    proj_meta.pLo_ch1 = pLo_ch1;  proj_meta.pHi_ch1 = pHi_ch1;
    proj_meta.pLo_ch3 = pLo_ch3;  proj_meta.pHi_ch3 = pHi_ch3;
    proj_meta.smooth_ch1 = smooth_ch1;
    proj_meta.gamma_ch1 = gamma_ch1;
    proj_meta.gamma_ch3_perz = gamma_ch3_perz;
    proj_meta.gamma_ch3_proj = gamma_ch3_proj;
    proj_meta.zoom = zoom;
    proj_meta.um_per_px_base = um_per_px_base;
    proj_meta.um_per_px = um_per_px;
    proj_meta.range_ch1 = range_ch1;
    proj_meta.range_ch3 = range_ch3;
    proj_meta.scaleBar_um = scaleBar_um;
    proj_meta.scaleBar_side = scaleBar_side;
    proj_meta.rotateCW = rotateCW;
    proj_meta.flipLR = flipLR;
    proj_meta.crop_x_um = crop_x_um;
    proj_meta.crop_y_um = crop_y_um;
    save(fullfile(stitchFolder, 'proj_meta.mat'), '-struct', 'proj_meta');

    fprintf('Saved 4 files + proj_meta.mat (%.4f um/px) to:\n  %s\n', um_per_px, stitchFolder);
end

%% ======== Figure: 2x2 ========
figure('Name','Projections','Color','k','Position',[50 50 1200 1000]);
tl = tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');

nexttile; imagesc(disp_ch1); axis image off; colormap gray;
title(lbl_ch1, 'Color','w','FontSize',10);

nexttile; imagesc(disp_ch3); axis image off; colormap gray;
title(lbl_ch3, 'Color','w','FontSize',10);

nexttile; image(merged_gm); axis image off;
title('Green=Ch1  Magenta=Ch3', 'Color','w','FontSize',10);

nexttile; image(merged_gr); axis image off;
title('Green=Ch1  Red=Ch3', 'Color','w','FontSize',10);

%% ======== Local functions ========
function out = clip_rescale_u16(img, pLo, pHi)
    vals = double(img(:));
    lo = prctile(vals, pLo);
    hi = prctile(vals, pHi);
    d = (double(img) - lo) / (hi - lo);
    d = max(0, min(1, d));
    out = uint16(d * 65535);
end

function out = apply_gamma_u16(img, gamma)
    d = double(img) / 65535;
    d = d .^ gamma;
    out = uint16(d * 65535);
end

function out = apply_gamma_stack(stack, gamma)
% Apply gamma per slice, normalizing by each slice's own max.
    cls = class(stack);
    out = zeros(size(stack), cls);
    for k = 1:size(stack,3)
        sl = double(stack(:,:,k));
        mx = max(sl(:));
        if mx <= 0, continue; end
        d = sl / mx;
        d = d .^ gamma;
        out(:,:,k) = cast(d * mx, cls);
    end
end

function patch_tiff_resolution(outFile, um_per_px)
% Reopen TIFF and patch resolution tags + ImageJ description.
% Uses r+ mode + rewriteDirectory for reliability.
    pixels_per_cm = 10000 / um_per_px;
    t = Tiff(outFile, 'r+');
    cln = onCleanup(@() t.close()); %#ok<NASGU>
    nDir = 1;
    while ~t.lastDirectory(), t.nextDirectory(); nDir = nDir+1; end
    % Patch every directory
    t.setDirectory(1);
    for d = 1:nDir
        t.setDirectory(d);
        t.setTag('XResolution', pixels_per_cm);
        t.setTag('YResolution', pixels_per_cm);
        t.setTag('ResolutionUnit', 3);  % Centimeter → ImageJ shows um
        if d == 1
            t.setTag('ImageDescription', sprintf('ImageJ=1.53\nunit=um\n'));
        end
        t.rewriteDirectory();
    end
end

function write_calibrated_tiff(outFile, img, um_per_px)
% Write single-page uint16 grayscale TIFF with pixel calibration.
    imwrite(uint16(img), outFile);
    patch_tiff_resolution(outFile, um_per_px);
end

function write_ij_composite(outFile, pages, nCh, um_per_px)
% Write N-channel ImageJ composite hyperstack.
% pages: cell array of uint16 images, one per channel.
% ImageJ default LUTs: Ch1=Red, Ch2=Green, Ch3=Blue, ...
    [H, W] = size(pages{1});

    desc = sprintf(['ImageJ=1.53\nimages=%d\nchannels=%d\n' ...
        'slices=1\nframes=1\nhyperstack=true\nmode=composite\nunit=um\n'], nCh, nCh);

    if isfile(outFile), delete(outFile); end
    t = Tiff(outFile, 'w');
    cln = onCleanup(@() t.close()); %#ok<NASGU>

    tag.ImageLength  = H;
    tag.ImageWidth   = W;
    tag.Photometric  = Tiff.Photometric.MinIsBlack;
    tag.SamplesPerPixel = 1;
    tag.BitsPerSample = 16;
    tag.SampleFormat  = Tiff.SampleFormat.UInt;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tag.Compression  = Tiff.Compression.None;
    tag.RowsPerStrip = H;

    for p = 1:nCh
        if p == 1
            tag.ImageDescription = desc;
        elseif isfield(tag, 'ImageDescription')
            tag = rmfield(tag, 'ImageDescription');
        end
        t.setTag(tag);
        t.write(uint16(pages{p}));
        if p < nCh
            t.writeDirectory();
        end
    end
    clear cln;

    patch_tiff_resolution(outFile, um_per_px);
    fprintf('[save] %d-ch composite: %s\n', nCh, outFile);
end

function img = burn_scalebar(img, bar_um, um_per_px, side)
% Burn a white scale bar into the image pixel data (bottom corner).
    [H, W] = size(img);
    bar_px = round(bar_um / um_per_px);
    margin = round(0.03 * min(H, W));
    barH   = max(4, round(H * 0.008));

    if side == "left"
        x1 = margin + 1;
        x2 = x1 + bar_px - 1;
    else
        x2 = W - margin;
        x1 = x2 - bar_px + 1;
    end
    y2 = H - margin;
    y1 = y2 - barH + 1;

    x1 = max(1, x1);  x2 = min(W, x2);
    y1 = max(1, y1);  y2 = min(H, y2);

    if isa(img, 'uint16')
        img(y1:y2, x1:x2) = 65535;
    else
        img(y1:y2, x1:x2) = max(img(:));
    end
end

function img = burn_scalebar_rgb(img, bar_um, um_per_px, side)
% Burn a white scale bar into an RGB image (HxWx3 double [0,1]).
    [H, W, ~] = size(img);
    bar_px = round(bar_um / um_per_px);
    margin = round(0.03 * min(H, W));
    barH   = max(4, round(H * 0.008));

    if side == "left"
        x1 = margin + 1;  x2 = x1 + bar_px - 1;
    else
        x2 = W - margin;  x1 = x2 - bar_px + 1;
    end
    y2 = H - margin;  y1 = y2 - barH + 1;
    x1 = max(1,x1); x2 = min(W,x2); y1 = max(1,y1); y2 = min(H,y2);

    img(y1:y2, x1:x2, :) = 1;  % white
end

function draw_scalebar(ax, bar_um, um_per_px, side)
    hold(ax, 'on');
    yl = ylim(ax);  xl = xlim(ax);
    imgH = yl(2) - yl(1);  imgW = xl(2) - xl(1);

    bar_px = bar_um / um_per_px;
    margin = 0.03;
    barH   = imgH * 0.008;

    if side == "left"
        x1 = xl(1) + margin * imgW;
        x2 = x1 + bar_px;
    else
        x2 = xl(2) - margin * imgW;
        x1 = x2 - bar_px;
    end
    y2 = yl(2) - margin * imgH;
    y1 = y2 - barH;

    fill(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], 'w', 'EdgeColor', 'none');
    text(ax, (x1+x2)/2, y1 - imgH*0.01, sprintf('%d \\mum', bar_um), ...
        'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    hold(ax, 'off');
end
