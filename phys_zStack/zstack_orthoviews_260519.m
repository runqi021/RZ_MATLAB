%% zstack_orthoviews_260519.m
% Sagittal + coronal depth-shaded MIP from a ScanImage z-stack TIFF.
%
% Auto-detects channel layout, zoom, and Z step from TIFF metadata.
% Deinterleaves the requested channel, builds a 3D volume (averaging
% framesPerSlice repeats per Z), and renders ImageJ-3D-Project-style
% depth-shaded MIPs along X (sagittal, Y vs Z) and Y (coronal, X vs Z).
%
% Outputs (next to input TIFF):
%   <stem>_ch<N>_sagittal_raw.tif   16-bit, no contrast stretch
%   <stem>_ch<N>_coronal_raw.tif    16-bit, no contrast stretch
%   <stem>_ch<N>_sagittal.png       8-bit, contrast-stretched + aspect-corrected
%   <stem>_ch<N>_coronal.png        8-bit, contrast-stretched + aspect-corrected
%
% Requires on path: loadtiff, saveastiff (in NoRMCorre-master/).

clear; clc; close all;

% Force white-on-black for all figures regardless of MATLAB dark theme.
set(groot, 'defaultFigureColor',         'w');
set(groot, 'defaultAxesColor',           'w');
set(groot, 'defaultAxesXColor',          'k');
set(groot, 'defaultAxesYColor',          'k');
set(groot, 'defaultAxesZColor',          'k');
set(groot, 'defaultTextColor',           'k');
set(groot, 'defaultAxesGridColor',       [0.15 0.15 0.15]);
set(groot, 'defaultAxesTickDir',         'out');
set(groot, 'defaultAxesTickDirMode',     'manual');
set(groot, 'defaultColorbarColor',       'k');   % colorbar tick labels & outline
% Note: no groot defaults exist for axes title/xlabel/ylabel text colors —
% those are set explicitly per-call in the figure code below.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'NoRMCorre-master'));

%% ============================== USER ==============================
tifPath = "C:\Users\Admin\Desktop\260517_shi_cal590\R_1.2x_x1100y850z20-620_00001.tif";

% -------- channel handling --------
AutoDetectChannels  = true;
KeepChannelID       = 2;
NumChannelsFallback = 1;

% -------- voxel scaling --------
PixelSizeBase = 1.7778;   % um/px at 1x zoom (lateral)
ZStepFallback = 5;        % um/slice if metadata missing

% -------- depth-shaded MIP (matches ImageJ 'Brightest Point' + depth-cueing) --------
Alpha       = 0;        % depth-shading strength (0..1)
                          %   0    -> pure MIP, no depth cue, all voxels full brightness
                          %   0.5  -> back of stack is 50% dimmer than front (recommended)
                          %   1.0  -> back of stack fades to black (cells deep in stack disappear)

% -------- preview --------
PreviewPctLo =0;         % low percentile for PNG contrast stretch
PreviewPctHi = 100;      % high percentile

%% ============================== LOAD + METADATA ==============================
assert(isfile(tifPath), "TIFF not found: %s", tifPath);
[recFolder, baseName, ~] = fileparts(char(tifPath));
baseName  = string(baseName);
outFolder = fullfile(string(recFolder), baseName);
if ~isfolder(outFolder), mkdir(char(outFolder)); end
fprintf("[setup] output folder: %s\n", outFolder);

fprintf("[load] %s\n", tifPath);
info   = imfinfo(tifPath);
nPages = numel(info);

% -- channel detection (same pattern as Batch_dffQC_260325.m)
savedCh = [];
if AutoDetectChannels
    savedCh = scanimage_saved_channels_from_info(info(1));
end
if isempty(savedCh)
    savedCh = 1:NumChannelsFallback;
end
savedCh = savedCh(:).';
nSaved  = numel(savedCh);

if mod(nPages, nSaved) ~= 0
    fprintf(2, "[chanDetect] WARNING: nPages=%d not divisible by nSaved=%d. Treating as single-channel.\n", nPages, nSaved);
    savedCh = KeepChannelID;
    nSaved  = 1;
end

kWithin     = find(savedCh == KeepChannelID, 1, "first");
keepID_used = KeepChannelID;
if isempty(kWithin)
    kWithin     = 1;
    keepID_used = savedCh(1);
    fprintf(2, "[chanDetect] KeepChannelID=%d not in savedCh=[%s]. Using %d.\n", ...
        KeepChannelID, num2str(savedCh), keepID_used);
end

fprintf("[chanDetect] savedCh=[%s] (nSaved=%d), keepID=%d (kWithin=%d)\n", ...
    num2str(savedCh), nSaved, keepID_used, kWithin);

% -- voxel sizes from metadata
meta = "";
if isfield(info(1),"Software") && ~isempty(info(1).Software)
    meta = meta + string(info(1).Software) + newline;
end
if isfield(info(1),"ImageDescription") && ~isempty(info(1).ImageDescription)
    meta = meta + string(info(1).ImageDescription) + newline;
end

zoomFactor = parse_scanimage_scalar(meta, "SI.hRoiManager.scanZoomFactor");
if ~isempty(zoomFactor) && isfinite(zoomFactor) && zoomFactor > 0
    pixel_size_um = PixelSizeBase / zoomFactor;
else
    pixel_size_um = PixelSizeBase;
    zoomFactor    = NaN;
    fprintf(2, "[meta] WARNING: could not read scanZoomFactor. Using PixelSizeBase=%.3f um/px.\n", pixel_size_um);
end

zStep_um = parse_scanimage_scalar(meta, "SI.hStackManager.stackZStepSize");
if isempty(zStep_um) || ~isfinite(zStep_um) || zStep_um == 0
    zStep_um = ZStepFallback;
    fprintf(2, "[meta] WARNING: could not read stackZStepSize. Using ZStepFallback=%.3f um.\n", zStep_um);
else
    % ScanImage uses negative stackZStepSize for stacks going down in Z.
    % For aspect-ratio purposes we only care about the magnitude.
    zStep_um = abs(zStep_um);
end

numSlices_meta      = parse_scanimage_scalar(meta, "SI.hStackManager.numSlices");
framesPerSlice_meta = parse_scanimage_scalar(meta, "SI.hStackManager.framesPerSlice");
if isempty(framesPerSlice_meta) || ~isfinite(framesPerSlice_meta) || framesPerSlice_meta < 1
    framesPerSlice = 1;
else
    framesPerSlice = round(framesPerSlice_meta);
end

% Per-slice absolute Z positions in microns (from user's pre-calibrated zero).
% SI.hStackManager.zs is a vector of length numSlices with the true Z value
% acquired at each slice (e.g. [20, 10, 0, -10, ..., -610]).
zs_um = parse_scanimage_vector(meta, "SI.hStackManager.zs");

fprintf("[meta] zoom=%.2fx, pixel=%.4f um/px, zStep=%.3f um, numSlices=%g, framesPerSlice=%d\n", ...
    zoomFactor, pixel_size_um, zStep_um, numSlices_meta, framesPerSlice);

% -- load full TIFF stack
fprintf("[load] reading %d pages...\n", nPages);
F = loadtiff(char(tifPath));

% -- deinterleave the requested channel
if nSaved > 1
    assert(mod(size(F,3), nSaved) == 0, "Pages=%d not divisible by nSaved=%d.", size(F,3), nSaved);
    Fkeep = F(:,:,kWithin:nSaved:end);
else
    Fkeep = F;
end
clear F

% -- collapse framesPerSlice
nZ_raw = size(Fkeep, 3);
if framesPerSlice > 1
    assert(mod(nZ_raw, framesPerSlice) == 0, ...
        "After deinterleave, %d frames not divisible by framesPerSlice=%d.", nZ_raw, framesPerSlice);
    nZ = nZ_raw / framesPerSlice;
    fprintf("[vol] averaging %d frames per Z slice -> Z=%d\n", framesPerSlice, nZ);
    vol = zeros(size(Fkeep,1), size(Fkeep,2), nZ, 'uint16');
    for z = 1:nZ
        idx0 = (z-1)*framesPerSlice + 1;
        idx1 = z*framesPerSlice;
        vol(:,:,z) = uint16(mean(single(Fkeep(:,:,idx0:idx1)), 3));
    end
else
    vol = Fkeep;
end
clear Fkeep

[H, W, Z] = size(vol);
fprintf("[vol] H=%d, W=%d, Z=%d, dtype=%s\n", H, W, Z, class(vol));

% Validate/fallback for per-slice Z positions.
if numel(zs_um) ~= Z
    if ~isempty(zs_um)
        fprintf(2, "[meta] WARNING: zs length %d != Z slices %d. Falling back to slice 1 = 0 um, step = %.2f um.\n", ...
            numel(zs_um), Z, zStep_um);
    else
        fprintf(2, "[meta] WARNING: could not read SI.hStackManager.zs. Falling back to slice 1 = 0 um, step = %.2f um.\n", zStep_um);
    end
    zs_um = (0:Z-1) * zStep_um;   % synthesized
end
zs_um = round(zs_um(:).');   % clean up float noise, force row vector
fprintf("[meta] zs = [%+.0f .. %+.0f] um (%d slices, span = %.0f um)\n", ...
    zs_um(1), zs_um(end), Z, abs(zs_um(end) - zs_um(1)));

%% ============================== DEPTH-SHADED PROJECTIONS ==============================
fprintf("[proj] Alpha=%.2f (0 = pure MIP, 1 = back fades to black)\n", Alpha);

fprintf("[proj] sagittal (project along X, dim 2)...\n");
proj_sag = depth_shaded_proj(vol, 2, Alpha);   % [H, Z]

fprintf("[proj] coronal (project along Y, dim 1)...\n");
proj_cor = depth_shaded_proj(vol, 1, Alpha);   % [W, Z]

%% ============================== TRANSPOSE: Z ON VERTICAL AXIS ==============================
% proj_sag is [H, Z] and proj_cor is [W, Z]. Transpose so Z is rows (vertical).
proj_sag = proj_sag.';   % [Z, H]
proj_cor = proj_cor.';   % [Z, W]

%% ============================== ASPECT-CORRECT (interpolate Z to isotropic) ==============================
% Resize the Z dim (rows) so 1 output px corresponds to 1 lateral px in microns.
% Both raw 16-bit and PNG outputs use this same isotropic resolution.
z_scale = zStep_um / pixel_size_um;
new_Z   = max(1, round(Z * z_scale));
fprintf("[isotropic] z_scale = %.3f lateral-px per Z-step -> resample Z %d -> %d rows\n", z_scale, Z, new_Z);

sag_iso = imresize(proj_sag, [new_Z, H], 'bilinear');
cor_iso = imresize(proj_cor, [new_Z, W], 'bilinear');

%% ============================== SAVE OUTPUTS ==============================
outStem = baseName + "_ch" + string(keepID_used);

out_sag_raw   = char(fullfile(outFolder, outStem + "_sagittal_raw.tif"));
out_cor_raw   = char(fullfile(outFolder, outStem + "_coronal_raw.tif"));
out_sag_png   = char(fullfile(outFolder, outStem + "_sagittal.png"));
out_cor_png   = char(fullfile(outFolder, outStem + "_coronal.png"));
out_avg_stack = char(fullfile(outFolder, outStem + "_avg_zstack.tif"));

opts = struct('overwrite', true, 'message', false, 'compress', 'no');

% Raw 16-bit projections, aspect-corrected.
saveastiff(uint16(round(max(0, min(65535, sag_iso)))), out_sag_raw, opts);
saveastiff(uint16(round(max(0, min(65535, cor_iso)))), out_cor_raw, opts);

% Preview PNGs, percentile-stretched.
sag_disp = stretch_to_uint8(sag_iso, PreviewPctLo, PreviewPctHi);
cor_disp = stretch_to_uint8(cor_iso, PreviewPctLo, PreviewPctHi);
imwrite(sag_disp, out_sag_png);
imwrite(cor_disp, out_cor_png);

% Per-slice averaged Z-stack (multi-page TIFF: one page per Z, each page is
% the framesPerSlice-averaged frame for that depth). Same as the in-memory
% vol — saved so you can scrub through slices in Fiji.
saveastiff(vol, out_avg_stack, opts);

% Per-slice PNGs with a 100 um scale bar burned into the bottom-right corner.
% Contrast-stretched using volume-wide percentiles so slices are comparable.
slicesDir = fullfile(outFolder, "slices");
if ~isfolder(slicesDir), mkdir(char(slicesDir)); end

vol_lims = prctile(double(vol(:)), [PreviewPctLo, PreviewPctHi]);
if vol_lims(2) <= vol_lims(1), vol_lims = [min(vol(:)), max(vol(:))]; end

bar_len_px = max(1, round(100 / pixel_size_um));
bar_thick  = 5;
margin     = 15;

for z = 1:Z
    sl       = double(vol(:,:,z));
    sl_scaled = (sl - vol_lims(1)) / (vol_lims(2) - vol_lims(1));
    sl_scaled = max(0, min(1, sl_scaled));
    sl_u8    = uint8(round(sl_scaled * 255));

    [hh, ww] = size(sl_u8);
    y2 = hh - margin;
    y1 = max(1, y2 - bar_thick + 1);
    x2 = ww - margin;
    x1 = max(1, x2 - bar_len_px + 1);
    sl_u8(y1:y2, x1:x2) = 255;

    out_png_z = char(fullfile(slicesDir, sprintf("%s_slice_%03d_z%+05dum.png", outStem, z, zs_um(z))));
    imwrite(sl_u8, out_png_z);
end

fprintf("Saved sagittal raw:     %s\n", out_sag_raw);
fprintf("Saved coronal  raw:     %s\n", out_cor_raw);
fprintf("Saved sagittal preview: %s\n", out_sag_png);
fprintf("Saved coronal  preview: %s\n", out_cor_png);
fprintf("Saved avg Z-stack:      %s (%d slices)\n", out_avg_stack, Z);
fprintf("Saved %d slice PNGs:    %s (100um bar = %d px)\n", Z, slicesDir, bar_len_px);

%% ============================== 4-PANEL PREVIEW (two-direction views of the full volume) ==============================
% Two views of the FULL volume from opposite sides:
%   Sagittal: viewed from medial side  vs  viewed from lateral side (mirror image)
%   Coronal:  viewed from rostral side vs  viewed from caudal  side (mirror image)
%
% With Alpha=0 (pure MIP) the two views are mirror images of identical content.
% With Alpha>0 the "near side" (front) appears bright, far side dim — so each
% view shows which cells are closer to that viewer.
%
% Convention: rostral = larger physical Y, medial = |X| closer to 0 (per user rule).

motorPos = parse_scanimage_vector(meta, "SI.hMotors.motorPosition");
if numel(motorPos) >= 2
    motor_X = motorPos(1);
    motor_Y = motorPos(2);
else
    motor_X = 0; motor_Y = 0;
    fprintf(2, "[meta] WARNING: could not read motorPosition. Direction assignments may be wrong.\n");
end

% physical coord at the two ends of each lateral axis
phys_X_low  = motor_X + (1   - W/2 - 0.5) * pixel_size_um;
phys_X_high = motor_X + (W   - W/2 - 0.5) * pixel_size_um;
phys_Y_low  = motor_Y + (1   - H/2 - 0.5) * pixel_size_um;
phys_Y_high = motor_Y + (H   - H/2 - 0.5) * pixel_size_um;

% Medial side = the end whose |X| is smaller. If that's the high-col end,
% then "view from medial" means front-of-projection is at high-col → reverse=true.
medial_at_high_col = abs(phys_X_high) < abs(phys_X_low);
% Rostral side = the end whose Y is larger. If that's the high-row end,
% then "view from rostral" means front-of-projection is at high-row → reverse=true.
rostral_at_high_row = phys_Y_high > phys_Y_low;

% -- Sagittal (project along X, dim 2). Same volume, opposite depth directions.
proj_sag_med_view = depth_shaded_proj(vol, 2, Alpha,  medial_at_high_col);   % front = medial side
proj_sag_lat_view = depth_shaded_proj(vol, 2, Alpha, ~medial_at_high_col);   % front = lateral side

% -- Coronal (project along Y, dim 1). Same volume, opposite depth directions.
proj_cor_ros_view = depth_shaded_proj(vol, 1, Alpha,  rostral_at_high_row);  % front = rostral side
proj_cor_cau_view = depth_shaded_proj(vol, 1, Alpha, ~rostral_at_high_row);  % front = caudal side

% Transpose for Z-on-vertical, then aspect-correct
sag_med_iso = imresize(proj_sag_med_view.', [new_Z, H], 'bilinear');
sag_lat_iso = imresize(proj_sag_lat_view.', [new_Z, H], 'bilinear');
cor_ros_iso = imresize(proj_cor_ros_view.', [new_Z, W], 'bilinear');
cor_cau_iso = imresize(proj_cor_cau_view.', [new_Z, W], 'bilinear');

% Mirror the "from opposite side" views so they look anatomically correct
% (looking at the brain from the lateral side flips left-right vs medial).
sag_lat_iso = fliplr(sag_lat_iso);
cor_cau_iso = fliplr(cor_cau_iso);

sag_med_label = sprintf("from medial side (X→0)");
sag_lat_label = sprintf("from lateral side (|X| large, mirrored)");
cor_ros_label = sprintf("from rostral side (Y>0)");
cor_cau_label = sprintf("from caudal side (Y<0, mirrored)");

% display 2x2 — Y-axis uses true ScanImage Z positions (zs_um) rather than
% assuming slice 1 = 0. e.g. for a stack going from 20 down to -610 um,
% top of image shows Z = 20, bottom shows Z = -610.
z_first        = zs_um(1);
z_last         = zs_um(end);
lateral_sag_um = H * pixel_size_um;
lateral_cor_um = W * pixel_size_um;

fig = figure('Color','w', 'Position', [100 100 1400 800]);
tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');

panels = {sag_med_iso, sag_lat_iso, cor_ros_iso, cor_cau_iso};
labels = {"Sagittal — " + sag_med_label, ...
          "Sagittal — " + sag_lat_label, ...
          "Coronal — "  + cor_ros_label, ...
          "Coronal — "  + cor_cau_label};
xls    = {lateral_sag_um, lateral_sag_um, lateral_cor_um, lateral_cor_um};
xlabs  = {"Y (\mum)", "Y (\mum)", "X (\mum)", "X (\mum)"};

for ii = 1:4
    nexttile;
    imagesc([0, xls{ii}], [z_first, z_last], panels{ii});
    axis image ij; colormap gray;
    cb = colorbar; cb.Color = 'k';
    clim_p = prctile(double(panels{ii}(:)), [PreviewPctLo, PreviewPctHi]);
    if clim_p(2) > clim_p(1), clim(clim_p); end
    title(labels{ii}, 'Interpreter','none', 'Color', 'k');
    xlabel(xlabs{ii}, 'Color', 'k');
    ylabel("Z (\mum, ScanImage abs)", 'Color', 'k');
    set(gca, 'XColor', 'k', 'YColor', 'k');
end

st = sgtitle(sprintf("%s  |  motorXY=[%.0f, %.0f]  |  Z range=[%+.0f, %+.0f] um  |  Alpha=%.2f", ...
    baseName, motor_X, motor_Y, z_first, z_last, Alpha), 'Interpreter','none');
st.Color = 'k';

fprintf("\n[done]\n");

%% ====================================================================
%% ====================== LOCAL FUNCTIONS =============================
%% ====================================================================

function out = depth_shaded_proj(vol, projDim, alpha, reverse)
% Max Intensity Projection with linear depth attenuation.
%   vol     : uint16 [H, W, Z]
%   projDim : 1 (project along Y -> coronal [W, Z]) or 2 (project along X -> sagittal [H, Z])
%   alpha   : depth-shading strength
%               0   -> pure MIP (no depth cue, all voxels at full brightness)
%               0.5 -> back of stack is 50% dimmer than front
%               1   -> back of stack fades to black
%   reverse : (optional) if true, "front" is at the highest index along projDim
%             (default false → front is at index 1)
%
% Output is single-precision, same intensity scale as input.

if nargin < 4, reverse = false; end

V = single(vol);
[mip, idx] = max(V, [], projDim);
mip = squeeze(mip);
idx = squeeze(idx);

N = size(V, projDim);
if N < 2
    out = mip;
    return
end

if reverse
    idx = N - idx + 1;
end

% depth_norm: 0 at front, 1 at back
depth_norm  = (single(idx) - 1) / single(N - 1);
attenuation = 1 - alpha * depth_norm;
out         = mip .* attenuation;
end

function out8 = stretch_to_uint8(img, pctLo, pctHi)
% Percentile-clip and rescale to uint8.
v    = double(img(:));
lims = prctile(v, [pctLo, pctHi]);
if lims(2) <= lims(1)
    out8 = zeros(size(img), 'uint8');
    return
end
scaled = (double(img) - lims(1)) / (lims(2) - lims(1));
scaled(scaled < 0) = 0;
scaled(scaled > 1) = 1;
out8 = uint8(round(scaled * 255));
end

%% ====================== ScanImage helpers (copied from Batch_dffQC_260325.m) =========

function savedCh = scanimage_saved_channels_from_info(info1)
meta = "";
if isfield(info1,"Software") && ~isempty(info1.Software)
    meta = meta + string(info1.Software) + newline;
end
if isfield(info1,"ImageDescription") && ~isempty(info1.ImageDescription)
    meta = meta + string(info1.ImageDescription) + newline;
end
meta = string(meta);

if strlength(strtrim(meta)) == 0
    savedCh = [];
    return
end

nAvail = parse_scanimage_scalar(meta, "SI.hChannels.channelsAvailable");
if isempty(nAvail) || ~isfinite(nAvail) || nAvail < 1
    nAvail = 4;
end
nAvail = min(max(round(nAvail),1), 16);

savedCh = parse_scanimage_channels(meta, "SI.hChannels.channelSave", nAvail);
if isempty(savedCh)
    savedCh = parse_scanimage_channels(meta, "SI.hChannels.channelsActive", nAvail);
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

function ch = parse_scanimage_channels(metaStr, key, nAvail)
ch = [];

pat = key + "\s*=\s*([^\r\n]+)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
rhs = strtrim(string(tok{1}));

if ~isempty(regexp(rhs, "^\d+$", "once"))
    m = str2double(rhs);
    if isfinite(m)
        ch = find(bitget(uint32(m), 1:nAvail));
        return;
    end
end

nums = regexp(rhs, "[-+]?\d+\.?\d*", "match");
if ~isempty(nums)
    v = str2double(string(nums));
    v = v(isfinite(v));
    if ~isempty(v) && all(mod(v,1)==0) && all(v>=1) && all(v<=nAvail)
        ch = unique(v(:).', "stable");
        return;
    end
end

tf = regexp(lower(rhs), "(true|false)", "match");
if ~isempty(tf)
    mask = strcmp(tf, "true");
    mask = mask(1:min(numel(mask), nAvail));
    ch = find(mask);
    return;
end
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
