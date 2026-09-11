%% vessel_XZ_projection_260529.m
% -----------------------------------------------------------------------
%  Single-channel blood-vessel volume -> X-Z plane projection (depth view).
%
%  Standalone. Does NOT stitch -- you stitch the vessel volume yourself
%  (manual stitch GUI / NChan_vol_stitch_zUse) and point this at the result,
%  OR point it at any single multi-page z-stack TIFF (one FOV).
%
%  Input  : a single-channel multi-page TIFF (each page = one Z slice).
%  Output : ONE depth projection in the X-Z plane (X across, Z = depth down).
%           MIP along the rostro-caudal axis (image rows / Y).
%
%  Calibration (PixelSize, Zstep) is entered MANUALLY below -- no metadata.
%  Display: percentile clip + gamma 0.5.   NO CLAHE.   NO invert.
%  Origin: one click sets 0,0,0; midline = vertical line through it.
%
%  Outputs (next to input TIFF, in a folder named after the stem):
%    <stem>_XZ_raw.tif   16-bit, aspect-corrected, no contrast stretch
%    <stem>_XZ.png       8-bit, clip + gamma display (no axes)
%    <stem>_XZ_fig.png   figure with calibrated um axes + midline
%    <stem>_XZ_data.mat  projection + calibration (origin px, um axes)
%
%  Requires on path: loadtiff, saveastiff (in NoRMCorre-master/).
% -----------------------------------------------------------------------

clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'NoRMCorre-master'));

%% ============================== USER ==============================
% --- input volume (single channel = blood vessel) ---
tifPath = "C:\Users\Admin\Desktop\live_251104_definitely_not_sst_tdtomato_fitc\Map\matlab_stitch\stitched_ch1_volume.tif";

% --- manual calibration (you enter these by hand) ---
pixel_size_um = 1.4815;   % lateral um per pixel   (e.g. 1.7778/zoom)
zStep_um      = 10;        % um per Z slice

% --- frames-per-slice averaging (set 1 if each page is already one slice) ---
FramesPerSlice = 1;

% --- Z slice subset ([] = all) ---
zUse = [];                % e.g. 2:28 to drop edge slices

% --- depth shading (0 = pure MIP, recommended for vessels) ---
Alpha = 0;                % 0 = no depth cue; 0.5 = back of stack 50% dimmer

% --- display (clip + gamma, NO CLAHE, NO invert) ---
clip_pct  = [0.5 99.5];   % percentile clip for normalization
gamma_val = 0.5;          % gamma correction

% --- midline calibration (clicked on the XY avg-proj, NOT the X-Z image) ---
CalibrateMidline = true;  % true -> click the midline on the XY avg-proj (sets X=0)
                          % false -> X=0 at image center
AvgProjPath = "";         % "" = auto-find *_avgproj.tif in the volume's folder
%% =================================================================

%% ---- load ----
assert(isfile(tifPath), "TIFF not found: %s", tifPath);
[recFolder, baseName, ~] = fileparts(char(tifPath));
baseName  = string(baseName);
outFolder = fullfile(string(recFolder), baseName);
if ~isfolder(outFolder), mkdir(char(outFolder)); end
fprintf("[setup] output folder: %s\n", outFolder);

fprintf("[load] reading %s ...\n", tifPath);
vol = loadtiff(char(tifPath));   % single channel assumed -> all pages are Z
nZ_raw = size(vol, 3);

%% ---- collapse frames-per-slice ----
if FramesPerSlice > 1
    assert(mod(nZ_raw, FramesPerSlice) == 0, ...
        "%d pages not divisible by FramesPerSlice=%d.", nZ_raw, FramesPerSlice);
    nZc = nZ_raw / FramesPerSlice;
    fprintf("[vol] averaging %d frames per Z slice -> Z=%d\n", FramesPerSlice, nZc);
    tmp = zeros(size(vol,1), size(vol,2), nZc, 'uint16');
    for z = 1:nZc
        idx0 = (z-1)*FramesPerSlice + 1;
        idx1 = z*FramesPerSlice;
        tmp(:,:,z) = uint16(mean(single(vol(:,:,idx0:idx1)), 3));
    end
    vol = tmp; clear tmp
end

%% ---- Z subset ----
nZ = size(vol, 3);
if ~isempty(zUse)
    zUse = zUse(:).';
    zUse = zUse(zUse >= 1 & zUse <= nZ);
    assert(~isempty(zUse), "zUse empty after clipping to 1..%d.", nZ);
    vol = vol(:,:,zUse);
end

[H, W, Z] = size(vol);
fprintf("[vol] H=%d, W=%d, Z=%d, dtype=%s\n", H, W, Z, class(vol));
fprintf("[cal] pixel=%.4f um/px, zStep=%.3f um\n", pixel_size_um, zStep_um);

%% ---- X-Z projection (project along Y / rows, dim 1) ----
fprintf("[proj] X-Z plane (project along Y, dim 1), Alpha=%.2f...\n", Alpha);
proj_xz = depth_shaded_proj(vol, 1, Alpha);   % [W, Z]
proj_xz = proj_xz.';                          % [Z, W]  -> Z on vertical

%% ---- aspect-correct (resample Z to isotropic lateral px) ----
z_scale = zStep_um / pixel_size_um;
new_Z   = max(1, round(Z * z_scale));
fprintf("[isotropic] z_scale = %.3f lateral-px per Z-step -> resample Z %d -> %d rows\n", ...
        z_scale, Z, new_Z);
xz_iso = imresize(proj_xz, [new_Z, W], 'bilinear');   % [new_Z, W], 1 row = 1 lateral px

%% ---- display: clip + gamma (no CLAHE, no invert) ----
lims = prctile(double(xz_iso(:)), clip_pct);
if lims(2) <= lims(1), lims = [min(xz_iso(:)), max(xz_iso(:))]; end
disp01 = (double(xz_iso) - lims(1)) / max(lims(2) - lims(1), eps);
disp01 = max(0, min(1, disp01)) .^ gamma_val;     % gamma 0.5

%% ---- midline calibration: click on the XY avg-proj (sets X = 0) ----
% The avg-proj shares the SAME lateral X (columns) as the volume, so a clicked
% column maps 1:1 to the X-Z X-axis. After isotropic resampling 1 px = um.
[Hd, Wd] = size(disp01);
if CalibrateMidline
    % locate the XY avg-proj (same folder as the volume)
    apPath = AvgProjPath;
    if strlength(apPath) == 0
        apGuess = regexprep(char(tifPath), '_volume\.tif$', '_avgproj.tif');
        if isfile(apGuess)
            apPath = string(apGuess);
        else
            d = dir(fullfile(recFolder, '*_avgproj.tif'));
            assert(~isempty(d), "No *_avgproj.tif in %s; set AvgProjPath.", recFolder);
            apPath = string(fullfile(d(1).folder, d(1).name));
        end
    end
    fprintf("[midline] avg-proj: %s\n", apPath);
    AP   = double(imread(char(apPath)));
    apLo = prctile(AP(:), clip_pct(1)); apHi = prctile(AP(:), clip_pct(2));
    APd  = max(0, min(1, (AP - apLo)/max(apHi - apLo, eps))) .^ gamma_val;

    figC = figure('Color','w','Name','Click the midline on the XY avg-proj (sets X=0)', ...
                  'Units','normalized','Position',[0.05 0.06 0.9 0.86]);
    imshow(APd, 'InitialMagnification','fit'); hold on;
    title('Click the MIDLINE (XY plane) \rightarrow sets X = 0', 'FontSize',12);
    [mx, ~] = ginput(1);
    plot([mx mx], [0.5 size(APd,1)+0.5], 'y--', 'LineWidth', 1.0);   % midline
    drawnow;
    % map avg-proj column -> X-Z column (same W axis; rescale if sizes differ)
    ox = mx * (Wd / size(APd,2));
else
    ox = (Wd + 1) / 2;   % image center
end
x_um = ((1:Wd) - ox) * pixel_size_um;   % lateral, 0 at midline
z_um = (0:Hd-1)       * pixel_size_um;  % depth from top slice
fprintf("[midline] x_px(XZ)=%.1f -> X in [%.0f, %.0f] um\n", ox, x_um(1), x_um(end));

%% ---- save raw 16-bit + gamma PNG ----
outStem = baseName + "_XZ";
out_raw = char(fullfile(outFolder, outStem + "_raw.tif"));
out_png = char(fullfile(outFolder, outStem + ".png"));
opts = struct('overwrite', true, 'message', false, 'compress', 'no');
saveastiff(uint16(round(max(0, min(65535, xz_iso)))), out_raw, opts);
imwrite(uint8(round(disp01 * 255)), out_png);

%% ---- calibrated figure ----
fig = figure('Color','w', 'Position', [100 100 1100 600]);
imagesc(x_um, z_um, disp01); hold on;
axis image ij; colormap gray; clim([0 1]);
plot([0 0], z_um([1 end]), 'y--', 'LineWidth', 1.0);   % midline at X=0
xlabel('X (\mum, 0 = midline)', 'Color','k');
ylabel('Z depth (\mum)', 'Color','k');
title(sprintf('%s  |  X-Z projection  |  px=%.3f um, zStep=%.2f um  |  gamma=%.2f', ...
      baseName, pixel_size_um, zStep_um, gamma_val), 'Interpreter','none', 'Color','k');
set(gca, 'XColor','k', 'YColor','k', 'TickDir','out');
exportgraphics(fig, char(fullfile(outFolder, outStem + "_fig.png")), ...
               'Resolution', 200, 'BackgroundColor','white');

%% ---- save data ----
save(char(fullfile(outFolder, outStem + "_data.mat")), ...
     'xz_iso', 'disp01', 'x_um', 'z_um', 'ox', ...
     'pixel_size_um', 'zStep_um', 'Alpha', 'clip_pct', 'gamma_val', ...
     'zUse', 'FramesPerSlice', 'tifPath');

fprintf("Saved raw:  %s\n", out_raw);
fprintf("Saved png:  %s\n", out_png);
fprintf("Saved fig:  %s\n", char(fullfile(outFolder, outStem + "_fig.png")));
fprintf("Saved mat:  %s\n", char(fullfile(outFolder, outStem + "_data.mat")));
fprintf("\n[done]\n");

%% ====================== LOCAL FUNCTIONS =============================
function out = depth_shaded_proj(vol, projDim, alpha)
% Max Intensity Projection with optional linear depth attenuation.
%   projDim : 1 (project along Y -> X-Z plane [W, Z])
%   alpha   : 0 = pure MIP; 1 = back of stack fades to black.
V = single(vol);
[mip, idx] = max(V, [], projDim);
mip = squeeze(mip);
idx = squeeze(idx);
N = size(V, projDim);
if N < 2 || alpha == 0
    out = mip;
    return
end
depth_norm = (single(idx) - 1) / single(N - 1);   % 0 front, 1 back
out        = mip .* (1 - alpha * depth_norm);
end
