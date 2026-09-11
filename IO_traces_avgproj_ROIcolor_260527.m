% IO_traces_avgproj_ROIcolor_260527.m
% -----------------------------------------------------------------------
%  ONE script, ONE IO recording, TWO matched plots + a combined summary.
%
%    Plot A  IO_traces_ROIcolor          : breath (top) + stacked dF/F
%                                           (bottom), every ROI drawn in its
%                                           OWN color (no text ROI labels).
%    Plot B  IO_avgproj_ROIcolor         : average projection of the
%                                           *_MC_MC.tif (grayscale) with each
%                                           ROI outline drawn in the SAME
%                                           color as its trace in Plot A.
%    Plot C  IO_combined_GIANT_SUMMARY   : A (left) + B (right) in one figure.
%
%  Color is the ONLY ROI identifier and it is shared between A and B, so a
%  trace can be located in the FOV purely by matching hue.
%
%  Raw pixel quality: the avg projection is written to PNG with imwrite at
%  NATIVE resolution (no figure rendering / resampling) -> *_avgproj_raw.png.
%  The on-screen / vector panels add the outlines + scale bars on top.
%
%  Dependencies: detect_session_fps.m, image processing toolbox (bwboundaries)
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

%% ===================== USER-EDITABLE PARAMETERS ======================
io_path  = 'C:\Users\Admin\Desktop\ChAT_analysis\0124\IO\roi3_R_-1000_140_2x_34lp_512x256_00001';

% ROIs to show (these are maskL label values == dFF column indices).
% [] = all ROIs. Order here sets the stacking order (bottom -> top).
io_sel   = [2, 3, 7, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 23, 22, 25, 26, 30, 29, 28, 31, 32, 33, 34];

nDrop        = 30;          % breath frames to toss (match calcium pipeline)
fallback_fps = 30;

% ---- trace panel (Plot A) ----
io_xlim   = [22 82];      % s, time window shown
io_tbar_s = 5;             % time scale bar (s)
io_dffbar = 1;             % dF/F scale bar (dF/F units)

% ---- avg projection panel (Plot B) ----
clip_pct    = [0.5 99.9];   % intensity clip percentiles for display
gamma_val   = 0.6;          % display gamma (<1 brightens midtones)
crop_um     = 5;            % crop from EACH side (um); 0 keeps full FOV
scaleBar_um = 50;           % um scale bar burned onto the projection
outlineLW   = 1.0;          % ROI outline line width
PixelSizeBase = 1.7778;     % um/px at zoom=1 (fallback if no pixelSize_um)

doSave   = true;
% =====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
[~, io_nm] = fileparts(io_path);

%% ===================== LOAD DATA =====================================
% --- locate files ---
dff_hit  = dir(fullfile(io_path, '*_ch1_dFF.mat'));
bp_hit   = dir(fullfile(io_path, '*DLC*breath_peak_data.mat'));
sam_hit  = dir(fullfile(io_path, '*_cpSAM_output.mat'));
mc_hit   = dir(fullfile(io_path, '*_MC_MC.tif'));
assert(~isempty(dff_hit), 'No *_ch1_dFF.mat in %s', io_path);
assert(~isempty(bp_hit),  'No *DLC*breath_peak_data.mat in %s', io_path);
assert(~isempty(sam_hit), 'No *_cpSAM_output.mat in %s', io_path);
assert(~isempty(mc_hit),  'No *_MC_MC.tif in %s', io_path);

fps = detect_session_fps(io_path, fallback_fps);
[~, sm] = detect_session_fps(io_path, fallback_fps);

% pixel size (um/px) from RAW tiff metadata
px_um = NaN;
if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um > 0
    px_um = sm.pixelSize_um;
elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor > 0
    px_um = PixelSizeBase / sm.zoomFactor;
end
assert(isfinite(px_um), 'Could not determine pixel size (um/px) from metadata.');

D    = load(fullfile(dff_hit(1).folder, dff_hit(1).name), 'dFF');
BP   = load(fullfile(bp_hit(1).folder,  bp_hit(1).name));
SAM  = load(fullfile(sam_hit(1).folder, sam_hit(1).name), 'maskL');
maskL = SAM.maskL;

dff_all = double(D.dFF);
nROI_total = size(dff_all, 2);

% --- breath trace (toss nDrop, detrend, demean) ---
bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);

% --- common length ---
T = min(size(dff_all,1), numel(bw));
dff_all = dff_all(1:T, :); bw = bw(1:T);
t = (0:T-1)'/fps;

% --- select ROIs + assign one color per ROI (shared A<->B) ---
if isempty(io_sel)
    sel = 1:nROI_total;
else
    sel = io_sel(io_sel >= 1 & io_sel <= nROI_total);
end
N = numel(sel);
colors = roi_colormap(N);          % [N x 3], max-spread hues
fprintf('IO %s : fps=%.2f, %.4f um/px, %d/%d ROIs selected\n', io_nm, fps, px_um, N, nROI_total);

%% ===================== AVG PROJECTION (raw PNG) =====================
V   = tiffreadVolume(fullfile(mc_hit(1).folder, mc_hit(1).name));
avg = mean(double(V), 3);

% crop crop_um from each side (also shifts ROI-outline coords)
cpx = round(crop_um / px_um);
if cpx > 0 && 2*cpx < min(size(avg))
    avg = avg(cpx+1:end-cpx, cpx+1:end-cpx);
else
    cpx = 0;
end

% clip -> [0,1] -> gamma  (display image; pixel grid untouched)
lo  = prctile(avg(:), clip_pct(1));
hi  = prctile(avg(:), clip_pct(2));
img = min(max((avg - lo) / max(hi - lo, eps), 0), 1) .^ gamma_val;
[H, W] = size(img);

% --- write the projection straight to PNG at NATIVE resolution ---
rawPng = fullfile(io_path, [io_nm '_avgproj_raw.png']);
if doSave
    imwrite(uint16(round(img * 65535)), rawPng);   % 16-bit, no resampling
    fprintf('Saved raw avg projection (native %dx%d px): %s\n', H, W, rawPng);
end

% --- ROI outlines (in cropped coords) ---
B = cell(N,1);
for i = 1:N
    b = bwboundaries(maskL == sel(i), 'noholes');
    for k = 1:numel(b), b{k} = b{k} - cpx; end       % shift for crop
    B{i} = b;
end

% scale bar geometry (px, YDir reverse / image coords)
barLen = min(max(1, round(scaleBar_um / px_um)), W - 2);
margin = round(0.04 * H);
barThk = max(3, round(0.012 * H));
barX   = margin;  barY = H - margin - barThk;

%% ===================== PLOT A : breath + stacked dF/F ===============
mio  = t >= io_xlim(1) & t <= io_xlim(2);
tw   = t(mio); bwW = bw(mio); dffW = dff_all(mio, sel);
norm01 = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);

rng_sel = max(dffW,[],1) - min(dffW,[],1);
sp = max(prctile(rng_sel, 80), 0.3);

fA = figure('Color','w','Name',['IO traces ' io_nm], ...
            'Units','normalized','Position',[0.05 0.08 0.55 0.84]);
axA1 = subplot(4,1,1); hold(axA1,'on');
plot(axA1, tw, norm01(bwW), 'Color',[0 0.35 1], 'LineWidth',0.7);
ylabel(axA1,'breath (norm)'); ylim(axA1,[-0.05 1.05]); xlim(axA1,io_xlim);
title(axA1, sprintf('IO  %s  (fps=%.2f, %d ROIs, color = ROI ID)', io_nm, fps, N), 'Interpreter','none');

axA2 = subplot(4,1,2:4); hold(axA2,'on');
for i = 1:N
    plot(axA2, tw, dffW(:,i) + (N-i)*sp, 'Color', colors(i,:), 'LineWidth', 0.6);  % inverted stack order
end
ylim(axA2,[-sp N*sp]); xlim(axA2,io_xlim);
set(axA2,'YTick',[]);                                % NO text ROI labels
xlabel(axA2,'Time (s)'); ylabel(axA2,'\DeltaF/F (stacked, color = ROI ID)');

% scale bars (vector): time (horizontal) + dF/F (vertical), bottom-left
xb = io_xlim(1) + 0.03*diff(io_xlim);  yb = -0.7*sp;
plot(axA2, [xb xb+io_tbar_s], [yb yb], 'k-', 'LineWidth',2);
plot(axA2, [xb xb], [yb yb+io_dffbar], 'k-', 'LineWidth',2);
text(axA2, xb+io_tbar_s/2, yb, sprintf('%g s', io_tbar_s), ...
     'Horizontal','center','Vertical','top','FontSize',9);
text(axA2, xb, yb+io_dffbar/2, sprintf(' %g \\DeltaF/F', io_dffbar), ...
     'Horizontal','left','Vertical','middle','FontSize',9);
linkaxes([axA1 axA2],'x');

%% ===================== PLOT B : avg proj + colored outlines =========
fB = figure('Color','w','Name',['IO avgproj ' io_nm]);
axB = axes(fB);
imshow(img, [0 1], 'Parent', axB, 'Border','tight'); colormap(axB, gray(256));
hold(axB,'on');
for i = 1:N
    for k = 1:numel(B{i})
        plot(axB, B{i}{k}(:,2), B{i}{k}(:,1), '-', 'Color', colors(i,:), 'LineWidth', outlineLW);
    end
end
rectangle(axB, 'Position',[barX barY barLen barThk], 'FaceColor','w', 'EdgeColor','none');
hold(axB,'off');
title(axB, sprintf('%s  avg proj + ROI outlines (color = ROI ID)', io_nm), 'Interpreter','none');

%% ===================== PLOT C : combined GIANT summary ==============
fC = figure('Color','w','Name',['IO combined ' io_nm], ...
            'Units','normalized','Position',[0.03 0.06 0.94 0.86]);
tl = tiledlayout(fC, 4, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('IO  %s   |   color = ROI ID (shared between panels)', io_nm), ...
      'Interpreter','none', 'FontWeight','bold');

% left column: breath (row 1) + stacked dF/F (rows 2-4)
axC1 = nexttile(tl, 1); hold(axC1,'on');
plot(axC1, tw, norm01(bwW), 'Color',[0 0.35 1], 'LineWidth',0.7);
ylabel(axC1,'breath'); ylim(axC1,[-0.05 1.05]); xlim(axC1,io_xlim);
set(axC1,'XTickLabel',[]);

axC2 = nexttile(tl, 3, [3 1]); hold(axC2,'on');
for i = 1:N
    plot(axC2, tw, dffW(:,i) + (N-i)*sp, 'Color', colors(i,:), 'LineWidth', 0.6);  % inverted stack order
end
ylim(axC2,[-sp N*sp]); xlim(axC2,io_xlim); set(axC2,'YTick',[]);
xlabel(axC2,'Time (s)'); ylabel(axC2,'\DeltaF/F (stacked)');
plot(axC2, [xb xb+io_tbar_s], [yb yb], 'k-', 'LineWidth',2);
plot(axC2, [xb xb], [yb yb+io_dffbar], 'k-', 'LineWidth',2);
text(axC2, xb+io_tbar_s/2, yb, sprintf('%g s', io_tbar_s), 'Horizontal','center','Vertical','top','FontSize',9);
text(axC2, xb, yb+io_dffbar/2, sprintf(' %g \\DeltaF/F', io_dffbar), 'Horizontal','left','Vertical','middle','FontSize',9);

% right column: avg proj + colored outlines (spans all 4 rows)
axC3 = nexttile(tl, 2, [4 1]);
imshow(img, [0 1], 'Parent', axC3, 'Border','tight'); colormap(axC3, gray(256));
hold(axC3,'on');
for i = 1:N
    for k = 1:numel(B{i})
        plot(axC3, B{i}{k}(:,2), B{i}{k}(:,1), '-', 'Color', colors(i,:), 'LineWidth', outlineLW);
    end
end
rectangle(axC3, 'Position',[barX barY barLen barThk], 'FaceColor','w', 'EdgeColor','none');
hold(axC3,'off');
title(axC3, 'avg proj + ROI outlines', 'Interpreter','none');

%% ===================== SAVE =========================================
if doSave
    % traces -> vector PDF (editable lines/text in Illustrator)
    exportgraphics(fA, fullfile(io_path,'IO_traces_ROIcolor.pdf'),  'ContentType','vector', 'BackgroundColor','white');
    % avg proj + outlines -> high-res PNG (raster image, crisp pixels)
    exportgraphics(fB, fullfile(io_path,'IO_avgproj_ROIcolor.png'), 'Resolution',600, 'BackgroundColor','white');
    % combined summary -> keep both as a reference layout
    exportgraphics(fC, fullfile(io_path,'IO_combined_GIANT_SUMMARY.png'), 'Resolution',300, 'BackgroundColor','white');
    exportgraphics(fC, fullfile(io_path,'IO_combined_GIANT_SUMMARY.pdf'), 'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(io_path,'IO_traces_avgproj_ROIcolor.mat'), 'sel','colors','io_xlim','px_um','fps','io_nm');
    fprintf('Saved IO_traces_ROIcolor.pdf + IO_avgproj_ROIcolor.png + IO_combined_GIANT_SUMMARY(.png/.pdf) + .mat to %s\n', io_path);
end

%% ========================= LOCAL FUNCTIONS ==========================
function cmap = roi_colormap(N)
% N maximally-spread, saturated hues. Golden-ratio hue stepping keeps
% stacking neighbors far apart in color so adjacent traces are distinct.
    if N <= 0, cmap = zeros(0,3); return; end
    gr  = 0.618033988749895;
    h   = mod((0:N-1) * gr, 1);
    s   = 0.85 * ones(1, N);
    v   = 0.90 * ones(1, N);
    cmap = hsv2rgb([h(:) s(:) v(:)]);
end
