% roi_hand_dff_260812.m
% -----------------------------------------------------------------------
%  dF/F from hand-drawn ImageJ/Fiji ROIs.
%
%  Draw ROIs in Fiji on the movie (or its average projection), save them as
%  .roi files (or one RoiSet.zip) into the FOV folder, point this at the
%  folder, run. No Cellpose, no segmentation -- just the pixels you circled.
%
%  Parsing is done by imagej_roi_read.m (repo root), which handles rect, oval,
%  polygon, freehand, traced and point ROIs.
%
%  DEFAULTS TO THE RAW TIFF, because that is what a hand ROI is normally drawn
%  on. Note that without motion correction the trace carries whatever motion
%  the FOV had -- a drifting ROI edge shows up as slow dF/F. Set TIF_PICK='mc'
%  if you drew on the registered movie instead.
%
%  Output (into fovPath):
%    roi_hand_dFF.mat   F, dFF, masks, roiInfo, fps, params
%    roi_hand_dFF.png / .pdf
%
%  Runqi Zhang / 2026-08-12

clear; clc;

%% ===================== USER-EDITABLE =====================
fovPath  = 'E:\260730_vgat-g8m_shiverer\phys\roi1_2.4x_x1000y1000_z320_3000f_32lp_00001';

roiSpec  = '';        % '' = every *.roi in fovPath; or a specific .roi / RoiSet.zip
TIF_PICK = 'raw';     % 'raw' | 'preproc' | 'mc'   (raw = the unprocessed acquisition)

nDropFrames    = 30;  % stabilisation frames to drop (raw only; MC already dropped them)
BaselineWinSec = 20;  % dF/F sliding-median baseline window (s)

doSave   = true;
% =========================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ---- 1. fps + the TIFF the ROIs belong to ----
[fps, ~] = detect_session_fps(fovPath);

% NOT fileparts(): FOV folder names contain dots ('roi1_2.4x_...'), and
% fileparts would treat '.4x_...' as an extension and return 'roi1_2'.
parts = strsplit(regexprep(fovPath, '[\\/]+$', ''), {'\', '/'});
stem  = parts{end};

switch lower(TIF_PICK)
    case 'raw',     tifPath = fullfile(fovPath, [stem '.tif']);
    case 'preproc', tifPath = fullfile(fovPath, [stem '_ch1_preproc.tif']);
    case 'mc',      tifPath = fullfile(fovPath, [stem '_ch1_preproc_MC_MC.tif']);
    otherwise, error('TIF_PICK must be raw|preproc|mc');
end
assert(isfile(tifPath), 'No %s TIFF: %s', TIF_PICK, tifPath);

% MC output already had the stabilisation frames removed upstream
if ~strcmpi(TIF_PICK,'raw') && contains(tifPath,'_MC'), nDropFrames = 0; end

tw = Tiff(tifPath,'r');
H  = getTag(tw,'ImageLength');
W  = getTag(tw,'ImageWidth');
T  = 1; while ~lastDirectory(tw); nextDirectory(tw); T = T+1; end

fprintf('%s TIFF: %s\n  %d x %d, %d frames @ %g fps\n', ...
        upper(TIF_PICK), tifPath, H, W, T, fps);

%% ---- 2. ROIs ----
if isempty(roiSpec); roiSpec = fovPath; end
R = imagej_roi_read(roiSpec, H, W);

keep = arrayfun(@(a) ~isempty(a.mask) && any(a.mask(:)), R);
if ~all(keep)
    fprintf('dropping %d ROI(s) with no enclosed area\n', nnz(~keep));
end
R = R(keep);
nR = numel(R);
assert(nR > 0, 'No usable ROIs.');

idx = cell(nR,1);
fprintf('\n%d ROI(s):\n', nR);
for k = 1:nR
    idx{k} = find(R(k).mask);
    fprintf('  %-22s %-9s %5d px   bbox [x %d y %d w %d h %d]\n', ...
            R(k).name, R(k).typeName, numel(idx{k}), R(k).bounds);
end

%% ---- 3. stream the stack, mean inside each ROI ----
F      = zeros(T, nR);
sumImg = zeros(H, W);

fprintf('\nreading %d frames...\n', T);
tic;
for t = 1:T
    setDirectory(tw, t);
    fr = single(read(tw));
    sumImg = sumImg + double(fr);
    for k = 1:nR
        F(t,k) = mean(fr(idx{k}));
    end
    if mod(t,500)==0; fprintf('  %d/%d (%.1fs)\n', t, T, toc); end
end
close(tw);
avgImg = sumImg / T;
fprintf('  done in %.1fs\n', toc);

%% ---- 4. dF/F ----
dout = helper.dFF_RZ(F, 'FPS', fps, ...
                     'DropFirstSec', nDropFrames/fps, ...
                     'BaselineWinSec', BaselineWinSec);
dFF  = dout.dFF;
tt   = (0:size(dFF,1)-1).'/fps;

fprintf('\ndF/F: [%d x %d], %.1f s\n', size(dFF,1), size(dFF,2), tt(end));
for k = 1:nR
    fprintf('  %-22s  mean F %8.1f   dF/F range %.3f .. %.3f\n', ...
            R(k).name, mean(F(:,k)), min(dFF(:,k)), max(dFF(:,k)));
end

%% ---- 5. figure ----
fh = figure('Color','w','Visible','off','Position',[60 60 1500 620]);
tl = tiledlayout(fh, 1, 3, 'TileSpacing','compact','Padding','compact');
cols = lines(max(nR,3));

ax1 = nexttile(tl,1);
lo = prctile(avgImg(:),1); hi = prctile(avgImg(:),99.7);
imagesc(ax1, avgImg, [lo hi]); axis(ax1,'image'); axis(ax1,'off');
colormap(ax1, gray(256)); hold(ax1,'on');
for k = 1:nR
    B = bwboundaries(R(k).mask);
    for b = 1:numel(B)
        plot(ax1, B{b}(:,2), B{b}(:,1), '-', 'Color', cols(k,:), 'LineWidth', 1.4);
    end
    text(ax1, R(k).bounds(1), R(k).bounds(2)-4, sprintf('%d', k), ...
         'Color', cols(k,:), 'FontWeight','bold', 'FontSize', 11);
end
title(ax1, sprintf('mean projection (%s) + hand ROIs', TIF_PICK), 'Interpreter','none');

ax2 = nexttile(tl,2,[1 2]); hold(ax2,'on');
step = 0;
if nR > 1; step = 1.1 * max(range(dFF,1)); end
for k = 1:nR
    plot(ax2, tt, dFF(:,k) + (k-1)*step, '-', 'Color', cols(k,:), 'LineWidth', 0.8);
    text(ax2, tt(1), (k-1)*step, sprintf(' %d: %s', k, R(k).name), ...
         'Color', cols(k,:), 'FontSize', 8, 'VerticalAlignment','bottom', ...
         'Interpreter','none');
end
xlim(ax2, [tt(1) tt(end)]); grid(ax2,'on');
xlabel(ax2,'time (s)'); ylabel(ax2,'\DeltaF/F');
title(ax2, sprintf('hand-ROI dF/F  |  %d ROI(s)  |  %g fps  |  baseline %g s', ...
                   nR, fps, BaselineWinSec));

title(tl, stem, 'Interpreter','none', 'FontWeight','bold');

%% ---- 6. save ----
if doSave
    roiInfo = rmfield(R, 'mask');
    masks   = cat(3, R.mask);
    outMat  = fullfile(fovPath, 'roi_hand_dFF.mat');
    params  = struct('tifPath',tifPath,'TIF_PICK',TIF_PICK,'fps',fps, ...
                     'nDropFrames',nDropFrames,'BaselineWinSec',BaselineWinSec);
    save(outMat, 'F','dFF','tt','masks','roiInfo','avgImg','params');

    base = fullfile(fovPath, 'roi_hand_dFF');
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','none');
    fprintf('\nsaved %s and %s.{png,pdf}\n', outMat, base);
end
close(fh);
