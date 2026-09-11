% avg_proj_scalebar_260527.m
% -----------------------------------------------------------------------
%  Average-projection of a *_MC_MC.tif with display polish:
%    - pixel size (um/px) from the RAW TIFF metadata (via detect_session_fps)
%    - crop crop_um from each side
%    - intensity clip to clip_pct percentiles, then gamma
%    - 50 um (scaleBar_um) white scale bar burned in
%
%  Outputs (into the MC_MC folder), image raster + vector scale bar/label:
%    <stem>_avgproj_gamma0.6.pdf   : vector PDF (best for Illustrator)
%    <stem>_avgproj_gamma0.6.png   : 600-dpi raster
%
%  Dependency: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

%% ===================== USER-EDITABLE =====================
mcPath      = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi3_z20_512x512_3x_2000f_00001\pFN_roi3_z20_512x512_3x_2000f_00001_ch1_preproc_MC_MC.tif";
scaleBar_um = 50;          % scale bar length (um)
gamma_val   = 0.6;         % display gamma (<1 brightens midtones)
clip_pct    = [0.5 99.9];  % intensity clip percentiles
crop_um     = 5;          % crop from EACH side (um)
fallback_fps = 30;
doSave      = true;
% =========================================================

%% ---- resolve folder (for raw-TIFF metadata) from the MC_MC path ----
assert(isfile(mcPath), 'MC_MC TIFF not found: %s', mcPath);
folderPath = fileparts(mcPath);
[~, stem]  = fileparts(mcPath);

%% ---- pixel size from RAW tiff metadata ----
[~, scan_meta] = detect_session_fps(folderPath, fallback_fps);
assert(isfield(scan_meta,'pixelSize_um') && isfinite(scan_meta.pixelSize_um), ...
    'detect_session_fps did not return pixelSize_um (need raw TIFF metadata).');
px_um = scan_meta.pixelSize_um;
fprintf('%s : %.4f um/px\n', stem, px_um);

%% ---- average projection ----
V   = tiffreadVolume(mcPath);          % [H x W x nFrames]
avg = mean(double(V), 3);

%% ---- crop crop_um from each side ----
cpx = round(crop_um / px_um);
if 2*cpx < min(size(avg))
    avg = avg(cpx+1:end-cpx, cpx+1:end-cpx);
else
    warning('crop_um too large for image; skipping crop.');
end

%% ---- clip percentiles -> [0,1] -> gamma ----
lo  = prctile(avg(:), clip_pct(1));
hi  = prctile(avg(:), clip_pct(2));
img = (avg - lo) / max(hi - lo, eps);
img = min(max(img, 0), 1);
img = img .^ gamma_val;

%% ---- display image (raster) + VECTOR scale bar (bottom-left) ----
[H, W]  = size(img);
barLen  = min(max(1, round(scaleBar_um / px_um)), W - 2);   % px
margin  = round(0.04 * H);
barThk  = max(3, round(0.012 * H));
barX    = margin;                       % left edge (px)
barY    = H - margin - barThk;          % top edge of bar (px, YDir reverse)

f = figure('Color','w','Name',stem);
imshow(img, [0 1], 'Border','tight'); colormap(gray(256)); hold on;
rectangle('Position',[barX barY barLen barThk], ...   % vector scale bar
          'FaceColor','w', 'EdgeColor','none');
hold off;

if doSave
    base = string(fullfile(folderPath, stem)) + "_avgproj_gamma0.6";
    exportgraphics(f, base + ".pdf", 'ContentType','vector', 'BackgroundColor','white');  % for Illustrator
    exportgraphics(f, base + ".png", 'Resolution',600,        'BackgroundColor','white');  % 600-dpi raster
    fprintf('Saved %s.pdf (vector) + .png (600 dpi)\n', base);
end
