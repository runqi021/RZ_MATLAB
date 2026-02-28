
%%
close all;
clearvars;      % <- safer than clear all; does not mess with pyenv
clc;

%% --- Set Python env for OASIS (Python) ---
PY = 'C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe';  % env with oasis-deconv

pe = pyenv;
if pe.Status == "NotLoaded"
    try
        pyenv('Version', PY, 'ExecutionMode', 'OutOfProcess');
    catch
        pyenv('Version', PY);
    end
else
    if ~contains(string(pe.Version), PY)
        warning("MATLAB Python already set to: %s\nTo change, restart MATLAB and rerun.", pe.Version);
    end
end

try
    py.importlib.import_module('oasis');
    py.importlib.import_module('oasis.functions');
catch ME
    error(['Python OASIS not available in this env.\n' ...
           'Make sure you ran:\n' ...
           '  conda activate cellpose-gpu\n' ...
           '  conda install -c conda-forge oasis-deconv\n\n' ...
           'Original error:\n%s'], ME.message);
end

%% ------------------------------------------------------------
% 0) LOAD DATA (keep your original logic)
% -------------------------------------------------------------
% if exist('out','var') && isstruct(out) && isfield(out,'F') && isfield(out,'t')
%     F = out.F;    t = out.t;
% else
%     [f,p] = uigetfile('*_cpSAM_output.mat','Select cpSAM output MAT');
%     S = load(fullfile(p,f));   % contains F, t, ...
%     F = S.F;    t = S.t;
% end
%
% % Optional downsample for huge recordings
% if numel(t) > 200000
%     stride = ceil(numel(t)/200000);
%     t = t(1:stride:end);
%     F = F(1:stride:end, :);
% end
% Path to TIFF
% Path to TIFF (you can hard-code or use uigetfile)
tiffPath = "C:\Users\zhang\Desktop\251124_live_vglut2_soma_g8s+cy5\phys\chan1\roi3_-650-380-120_46lp_930_x1.4_256x512_3000f_00001_ch1\roi3_-650-380-120_46lp_930_x1.4_256x512_3000f_00001_ch1_MC.tif";
% [fname, fpath] = uigetfile('*.tif','Select MC TIFF');   % optional
% tiffPath = fullfile(fpath, fname);

% base folder of this TIFF
[tdir, tname, ~] = fileparts(tiffPath);

% automatically create output "event" subfolder next to the TIFF
outDir = fullfile(tdir, 'event');
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output folder:\n  %s\n', outDir);
else
    fprintf('Using existing output folder:\n  %s\n', outDir);
end

% file to store / reload ROI selection
roiSelFile = fullfile(tdir, sprintf('%s_roiSelection.mat', tname));

%%
% base folder of the TIFF
[tdir, ~, ~] = fileparts(tiffPath);

% find ANY mat file and blindly load ALL variables into workspace
matFiles = dir(fullfile(tdir, '*_cpSAM_output.mat'));

if isempty(matFiles)
    error('No *_cpSAM_output.mat found in folder:\n  %s', tdir);
end

matPath = fullfile(tdir, matFiles(1).name);

fprintf('Loading EVERYTHING from:\n  %s\n', matPath);

load(matPath);   % <--- THIS loads EVERY variable into your workspace exactly as if you double-clicked

dFF_raw = dFF;

%%
% Read TIFF info
info = imfinfo(tiffPath);
nFrames = numel(info);

H = info(1).Height;
W = info(1).Width;

% Preallocate accumulator
avgImg = zeros(H, W);

% Accumulate frames
for k = 1:nFrames
    frame = double(imread(tiffPath, k));
    avgImg = avgImg + frame;
end

% Divide to get average
avgImg = avgImg / nFrames;

%% ------------------------------------------------------------
% 1) Toss away the first frames
nDrop = 30;   % number of frames to toss

F(1:nDrop, :) = [];
t(1:nDrop)    = [];

%%
%dFF = dFFout.dFF;
%%
[T, N] = size(dFF);
fs = 30;                % Hz
dt = 1/fs;

% %% Plot: stacked ΔF/F (not denoised)
% amp = max(dFF(:)) - min(dFF(:));  % global amplitude
% if amp == 0
%     amp = 1;  % safety
% end
% gap = 0.5 * amp;                   % dynamic spacing
% Y = dFF + gap*(0:N-1);
% 
% figure('Name','All ROI traces (ΔF/F, stacked)');
% plot(t, Y, 'LineWidth', 0.1);
% hold on;   % <<< IMPORTANT so we can draw highlights without erasing traces
% 
% xlim([t(1) t(end)]);
% ylim([-gap, gap*(N-1)+gap]);
% yt  = 0:gap:gap*(N-1);
% ytk = 1:N;
% keepTicks = mod(ytk,5)==1 | ytk==N;
% set(gca,'YTick', yt(keepTicks), 'YTickLabel', string(ytk(keepTicks)));
% xlabel('Time (s)');
% ylabel('ROI (offset)');
% title('All ROI traces (ΔF/F, stacked)');

% define where to store selection for this movie
[tdir, tname, ~] = fileparts(tiffPath);
roiSelFile = fullfile(tdir, sprintf('%s_roiSelection.mat', tname));
%% ------------------------------------------------------------
% 2) --- OASIS deconvolution on F (Python) ---
% -------------------------------------------------------------
fprintf('Running OASIS (default AR) on %d ROIs...\n', N);
g=0.95;
[F_denoised, dFF_denoised, spikes_oasis, baseline_oasis] = ...
    oasis_denoise_and_dff_AR1(F, g);

dFF = dFF_denoised;   % use this for everything downstream
% --- enforce mean ~0 per ROI (or use median if you prefer robustness) ---
dFF = dFF_denoised - mean(dFF_denoised, 1);   % T x N
%%
%figure; plot(dFF(:,9));
%figure; plot(dFF(:,33));

% 3) Plot: stacked ΔF/F (denoised)
amp = max(dFF(:)) - min(dFF(:));  % global amplitude
if amp == 0
    amp = 1;  % safety
end
gap = 0.5 * amp;                   % dynamic spacing
Y = dFF + gap*(0:N-1);
Y_raw = dFF_raw + gap*(0:N-1);

figure('Name','All ROI traces (ΔF/F, stacked)');
plot(t, Y, 'LineWidth', 0.1, 'Color', 'r');
hold on;
plot(t, Y_raw, 'LineWidth', 0.1, 'Color', 'k');

hold on;

% ---- Choose scale lengths ----
dFF_scale = 2;   % vertical bar = 0.5 ΔF/F (adjust as you like)
t_scale   = 2;     % horizontal bar = 2 s

% ---- Position near upper-left corner ----
xRange = t(end) - t(1);
yTop   = gap*(N-3) + 0.8*gap;      % just below top trace

x0 = xRange;         % a bit to the right of left edge
y0 = yTop - dFF_scale/2;           % center vertical bar around yTop

% ---- Vertical scale bar: ΔF/F ----
plot([x0 x0], [y0 y0 + dFF_scale], 'k', 'LineWidth', 2);

text(x0, y0 + dFF_scale + 0.05*dFF_scale, ...
     sprintf('%.2f \\DeltaF/F', dFF_scale), ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize', 9);


hold on;   % <<< IMPORTANT so we can draw highlights without erasing traces

xlim([t(1) t(end)]);
ylim([-gap, gap*(N-1)+gap]);
yt  = 0:gap:gap*(N-1);
ytk = 1:N;
keepTicks = mod(ytk,5)==1 | ytk==N;
set(gca,'YTick', yt(keepTicks), 'YTickLabel', string(ytk(keepTicks)));
xlabel('Time (s)');
ylabel('ROI (offset)');
title('All ROI traces (ΔF/F, stacked)');

%% 3) Plot: stacked ΔF/F (denoised)
amp = max(dFF(:)) - min(dFF(:));  % global amplitude
if amp == 0
    amp = 1;  % safety
end
gap = 0.5 * amp;                   % dynamic spacing
Y = dFF + gap*(0:N-1);
figure('Name','All ROI traces (ΔF/F, stacked)');
plot(t, Y, 'LineWidth', 0.1);

hold on;

% ---- Choose scale lengths ----
dFF_scale = 2;   % vertical bar = 0.5 ΔF/F (adjust as you like)
t_scale   = 2;     % horizontal bar = 2 s

% ---- Position near upper-left corner ----
xRange = t(end) - t(1);
yTop   = gap*(N-3) + 0.8*gap;      % just below top trace

x0 = xRange;         % a bit to the right of left edge
y0 = yTop - dFF_scale/2;           % center vertical bar around yTop

% ---- Vertical scale bar: ΔF/F ----
plot([x0 x0], [y0 y0 + dFF_scale], 'k', 'LineWidth', 2);

text(x0, y0 + dFF_scale + 0.05*dFF_scale, ...
     sprintf('%.2f \\DeltaF/F', dFF_scale), ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize', 9);


hold on;   % <<< IMPORTANT so we can draw highlights without erasing traces

xlim([t(1) t(end)]);
ylim([-gap, gap*(N-1)+gap]);
yt  = 0:gap:gap*(N-1);
ytk = 1:N;
keepTicks = mod(ytk,5)==1 | ytk==N;
set(gca,'YTick', yt(keepTicks), 'YTickLabel', string(ytk(keepTicks)));
xlabel('Time (s)');
ylabel('ROI (offset)');
title('All ROI traces (ΔF/F, stacked)');


% define where to store selection for this movie
[tdir, tname, ~] = fileparts(tiffPath);
roiSelFile = fullfile(tdir, sprintf('%s_roiSelection.mat', tname));
% --- Scale bar for stacked ΔF/F plot ---
%% 3b) Interactively exclude bad ROIs based on the stacked traces
% Click in the y-direction on a trace band to toggle exclusion for that ROI.
% Left-click to toggle, right-click / any key to finish.

centers = gap * (0:N-1);    % vertical centers of each ROI "lane"

% ---- start with everything included ----
excludeMask = false(1, N);

% ---- If we have a saved selection, preload it ----
if exist(roiSelFile, 'file')
    Ssel = load(roiSelFile);

    if isfield(Ssel, 'goodMask') && numel(Ssel.goodMask) == N
        goodMask_prev = logical(Ssel.goodMask);
        excludeMask   = ~goodMask_prev;
        fprintf('Preloaded ROI selection from %s (goodMask)\n', roiSelFile);

    elseif isfield(Ssel, 'goodMask_final') && numel(Ssel.goodMask_final) == N
        goodMask_prev = logical(Ssel.goodMask_final);
        excludeMask   = ~goodMask_prev;
        fprintf('Preloaded ROI selection from %s (goodMask_final)\n', roiSelFile);

    elseif isfield(Ssel, 'excludeIdx')
        prevExcl = Ssel.excludeIdx(:);
        prevExcl = prevExcl(prevExcl >= 1 & prevExcl <= N);
        excludeMask(prevExcl) = true;
        fprintf('Preloaded %d excluded ROIs from %s (excludeIdx)\n', ...
                numel(prevExcl), roiSelFile);

    elseif isfield(Ssel, 'excludeIdx_final')
        prevExcl = Ssel.excludeIdx_final(:);
        prevExcl = prevExcl(prevExcl >= 1 & prevExcl <= N);
        excludeMask(prevExcl) = true;
        fprintf('Preloaded %d excluded ROIs from %s (excludeIdx_final)\n', ...
                numel(prevExcl), roiSelFile);
    else
        fprintf('ROI selection file found but no usable fields. Starting fresh.\n');
    end
end

% ---- draw existing highlights for pre-excluded ROIs (if any) ----
hHighlight = gobjects(1, N);

for roi = find(excludeMask)
    yCenter = centers(roi);   % band center
    hHighlight(roi) = plot([t(1) t(end)], [yCenter yCenter], ...
                           'k-', 'LineWidth', 1.5);
end

disp('Left-click a trace band to toggle exclusion. Right-click or press any key to finish.');

while true
    try
        [~, yClick, button] = ginput(1);
    catch
        % figure closed
        break;
    end

    if isempty(button)
        % key pressed
        break;
    end

    % right or middle click → finish
    if button ~= 1
        break;
    end

    % find nearest ROI lane by y-position
    [~, roi] = min(abs(yClick - centers));

    if roi < 1 || roi > N
        continue;
    end

    % toggle
    excludeMask(roi) = ~excludeMask(roi);

    if excludeMask(roi)
        % draw a bold black line along the center of that band
        if ~isgraphics(hHighlight(roi))
            yCenter = centers(roi);
            hHighlight(roi) = plot([t(1) t(end)], [yCenter yCenter], ...
                                   'k-', 'LineWidth', 1.5);
        end
        fprintf('ROI %d -> EXCLUDE\n', roi);
    else
        % remove highlight
        if isgraphics(hHighlight(roi))
            delete(hHighlight(roi));
            hHighlight(roi) = gobjects(1);
        end
        fprintf('ROI %d -> KEEP\n', roi);
    end
end

hold off;

% ---- build goodMask / keep/exclude lists ----
goodMask = ~excludeMask;
keepIdx  = find(goodMask);
excludeIdx = find(~goodMask);

fprintf('Good ROIs (after interactive step): %d / %d\n', sum(goodMask), N);

%% 3c) Manual extra exclusions (simple)
manExclude = [];   % <---- edit by hand if needed, e.g. [29 72 74]

if ~isempty(manExclude)
    manExclude = manExclude(:);
    manExclude = manExclude(manExclude >= 1 & manExclude <= N);  % valid only
    goodMask(manExclude) = false;
end

% update keep/exclude lists after manual step
keepIdx    = find(goodMask);
excludeIdx = find(~goodMask);

fprintf('After manual exclusion: %d good / %d total.\n', sum(goodMask), N);

% ---- SAVE final selection to disk for next run ----
save(roiSelFile, 'manExclude', 'excludeIdx', 'goodMask', 'keepIdx');
fprintf('Saved ROI selection to %s\n', roiSelFile);

%% 4) Heatmap of ΔF/F (good ROIs only), correct y-axis
dFF_good = dFF(:, goodMask);   % T × Ngood
Ngood    = numel(keepIdx);

% compute percentiles on good ROIs only
pLow  = prctile(dFF_good(:), 0);
pHigh = prctile(dFF_good(:), 99);

figure('Name','ΔF/F heatmap (good ROIs only)');
imagesc(t, 1:Ngood, dFF_good');   % rows = 1:Ngood (even spacing)
set(gca, 'YDir', 'normal');

xlabel('Time (s)');
ylabel('ROI (original indices)');
title('ΔF/F heatmap (good ROIs only)');
colormap(parula);

% percentile scaling
caxis([pLow, pHigh]);

% y-tick labels → original ROI indices
set(gca, 'YTick', 1:Ngood, 'YTickLabel', string(keepIdx));

% colorbar
cb = colorbar;
ylabel(cb, '\DeltaF/F');


%% ------------------------------------------------------------
% 6) Correlation and covariance matrices
% -------------------------------------------------------------
R = corrcoef(dFF(:, goodMask));   % correlation between good ROIs only
      % N x N
%C = cov(dFF);           % N x N

figure('Name','Correlation matrix');
imagesc(R);
axis square;
colorbar;
caxis([0 1]);           % IO can be quite correlated
title('Pairwise correlation (ΔF/F)');
xlabel('ROI'); ylabel('ROI');

% figure('Name','Covariance matrix');
% imagesc(C);
% axis square;
% colorbar;
% title('Covariance (ΔF/F)');
% xlabel('ROI'); ylabel('ROI');

%% event detection
% assume: dFF is T x N, fs = 30 Hz
[T,N] = size(dFF);

% --- peak detection (per ROI) ---
pks  = cell(N,1);
locs = cell(N,1);

for i = 1:N
    if ~goodMask(i)
        % mark bad ROI: no peaks
        pks{i}  = [];
        locs{i} = [];
        continue;
    end

    [pks{i}, locs{i}] = findpeaks(dFF(:, i), ...
        'MinPeakProminence', 0.25, 'MinPeakHeight', 0.16, 'MinPeakDistance', 3);
end

% --- binary spike/event trains S(t, ROI) ---
S = zeros(T, N);   % 0/1 events

for i = 1:N
    if ~isempty(locs{i})
        S(locs{i}, i) = 1;
    end
end

figure('Name','Peak detection (good ROIs only)','Color','w');
tiledlayout('flow');  % automatically fills panels

roiList = find(goodMask);   % list of good ROI indices
M = numel(roiList);

for k = 1:M
    i = roiList(k);         % actual ROI index in original numbering

    nexttile;
    plot(t, dFF(:, i), 'k'); hold on;

    if ~isempty(locs{i})
        plot(t(locs{i}), dFF(locs{i}, i), 'ro', 'MarkerFaceColor','r', 'MarkerSize',3);
    end

    title(sprintf('ROI %d', i));   % keep original ROI index in title
    xlabel('Time (s)');
    ylabel('\DeltaF/F');
end

% Save peak-detection figure
outFig = fullfile(tdir, sprintf('%s_peakDetection.pdf', tname));
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved peak-detection figure:\n  %s\n', outFig);


%%
%figure
%plot(dFF(:,75))
%%
R_peaks = corrcoef(S);   % N x N
figure;
imagesc(R_peaks(keepIdx,keepIdx));
axis square;
colorbar;
caxis([0 1]);            % or [-1 1] if you detrend differently
xlabel('ROI'); ylabel('ROI');
title('Pairwise correlation of peak trains (binary)');

%% 2) Max cross-correlation within ±maxLag for all pairs
activeROIs = find(sum(S,1) > 0);   % row vector of ROI indices
M = numel(activeROIs);             % number of valid ROIs

[T, N] = size(S);
maxLagSec = 0.05;
maxLag = round(maxLagSec * fs);

Rpeak_active   = nan(M);
LagPeak_active = nan(M);


for ii = 1:M
    i = activeROIs(ii);

    for jj = 1:M
        j = activeROIs(jj);

        % compute cross-correlation only for active ROIs
        [cij, lags] = xcorr_nom(S(:,i), S(:,j), maxLag);

        [rmax, idx] = max(cij);
        Rpeak_active(ii,jj)   = rmax;
        LagPeak_active(ii,jj) = lags(idx) / fs;

        
    end
end

Rpeak   = nan(N);
LagPeak     = nan(N);

Rpeak(activeROIs, activeROIs) = Rpeak_active;
LagPeak(activeROIs, activeROIs)   = LagPeak_active;

%%
% 2a) matrix of peak cross-correlation values
figure('Name','Peak cross-corr (binary peaks)');
imagesc(Rpeak_active);
axis square;
colormap(hot);
colorbar;
caxis([0 max(max(Rpeak_active))]);   % adjust if needed
xlabel('ROI'); ylabel('ROI');
title(sprintf('Max cross-correlation within \\pm %.2f s', maxLagSec));

%%
LagAbs = abs(LagPeak_active);   % use abs value for color intensity

figure('Name','Lag of peak cross-corr');

% use ABS value for coloring (symmetry)
LagAbs = abs(LagPeak_active);
imagesc(LagAbs);
axis square;

% --- symmetric blue↔yellow colormap ---
cmap = flipud(parula);     % parula(1)=yellow, parula(end)=blue → flipped
colormap(cmap);

% symmetric color scaling
caxis([0 maxLagSec]);

colorbar;
xlabel('ROI'); 
ylabel('ROI');
title('Lag magnitude (|lag|), yellow = 0-lag (strong synchrony)');

%% Raster (vertical bar style) for good ROIs only
Ngood    = numel(keepIdx);

figRaster = figure('Name','ROI Event Raster (Good ROIs, Bars)','Color','w'); 
hold on;

barHeight = 0.8;   % length of each vertical bar

for k = 1:Ngood
    r = keepIdx(k);               % actual ROI index
    tt = t(S(:, r) == 1);          % event times for that ROI
    
    % draw a vertical line for each detected event
    for i = 1:numel(tt)
        plot([tt(i) tt(i)], [k - barHeight/2, k + barHeight/2], ...
            'k-', 'LineWidth', 2);
    end
end

xlabel('Time (s)');
ylabel('ROI (original index)');
title('ROI Event Raster (Vertical-Bar Style, Good ROIs Only)');

ylim([0 Ngood+1]);
set(gca,'YDir','normal');

% y-tick labels = original ROI IDs
set(gca,'YTick', 1:Ngood, 'YTickLabel', string(keepIdx));

grid off;

% ---- SAVE RASTER FIGURE ----
outRaster = fullfile(tdir, sprintf('%s_raster.pdf', tname));
exportgraphics(figRaster, outRaster, 'ContentType','vector');

fprintf('Saved raster plot:\n  %s\n', outRaster);
%%
binSec = 0.25;              % 0.5 s bin
bin    = round(binSec * fs);   % = 15 samples at fs = 30 Hz

binWin = 0;                % only the central bin → total window = 0.5 s

[T, N] = size(S);

t = (0:T-1)'/fs;

pop_raw = sum(S, 2);             % instantaneous population spike count

nBins = floor(T / bin);
PSTH = zeros(nBins, 1);

for k = 1:nBins
    idx = ( (k-1)*bin + 1 ) : (k*bin);
    PSTH(k) = sum(pop_raw(idx));
end

tPSTH = ((0:nBins-1) * binSec)';   % time of each bin center

PSTH_rate = PSTH / binSec;

figure('Color','w');
bar(tPSTH, PSTH_rate, 'histc');
xlabel('Time (s)');
ylabel('Population rate (spikes/sec)');
title(sprintf('Population PSTH firing rate (bin = %.2f s)', binSec));
grid on;

%%
g = gausswin(5);
g = g/sum(g);

PSTH_sm = conv(PSTH_rate, g, 'same');

figure('Color','w');
plot(tPSTH, PSTH_sm, 'LineWidth',1.2, 'Color', 'k');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
title('PSTH');
grid on;

%%
baseline = median(PSTH_sm);
noise    = mad(PSTH_sm - baseline, 1);   % robust estimate

minHeight = baseline + 2.5*noise;          % event threshold
minProm   = 2*noise;                     % ensures clear rise

[pks_pop, locs_pop] = findpeaks(PSTH_sm, ...
    "MinPeakProminence", minProm, ...
    "MinPeakHeight",     minHeight, ...
    'MinPeakDistance',   binWin*binSec)  % avoid detecting same event twice;

evtTimes = tPSTH(locs_pop);    % event time in seconds
evtFrames = round(evtTimes * fs);   % map to nearest frame in S

figure('Color','w');
plot(tPSTH, PSTH_sm, 'k', 'LineWidth', 1.5); 
hold on;

% Plot peaks
if ~isempty(locs_pop)
    plot(tPSTH(locs_pop), pks_pop, 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 6);
end

xlabel('Time (s)');
ylabel('Firing rate (Hz)');
title('Population activity');
grid on;

legend({'PSTH', 'Detected peaks'}, 'Location', 'best');

events = struct([]);

for k = 1:numel(locs_pop)
    b0 = locs_pop(k);
    binsEvt = (b0 - binWin) : (b0 + binWin);
    binsEvt = binsEvt(binsEvt >= 1 & binsEvt <= nBins);

    % convert bins to frames
    framesEvt = [];
    for b = binsEvt
        fStart = (b-1)*bin + 1;
        fEnd   = min(b*bin, T);
        framesEvt = [framesEvt, fStart:fEnd];
    end

    framesEvt = unique(framesEvt);

    events(k).loc    = locs_pop(k);
    events(k).pks    = pks_pop(k);
    events(k).frames = framesEvt;
    events(k).frame_center = round(mean(framesEvt));
    events(k).t_center     = events(k).frame_center / fs;
end

fprintf("Detected %d global events via findpeaks.\n", numel(events));

%% Collect pairwise |Δt| and distance across all events
allLag  = [];   % will hold |Δt| (s) for each ROI pair, all events
allDist = [];   % will hold pairwise distances (µm or pixels)

for kEvent = 1:numel(events)
    %kEvent = 6;        % choose event
    ev     = events(kEvent);
    
    framesEvt  = ev.frames;
    frameCenter = ev.frame_center;
    
    dt  = nan(1, N);     % Δt per ROI
    act = false(1, N);
    
    for r = 1:N
        idxLocal = find(S(framesEvt, r) == 1, 1, 'first');
        if ~isempty(idxLocal)
            act(r) = true;
            frame_r = framesEvt(idxLocal);
            dt(r) = (frame_r - frameCenter) / fs;   % in seconds
        end
    end
    
    fprintf('Event %d: %d active ROIs.\n', kEvent, sum(act));
    
    % event train map
    % --- 1) ROI centroids from maskL ---
    labels = unique(maskL);
    labels(labels==0) = [];
    N = numel(labels);
    
    props = regionprops(maskL,'Centroid');
    cent = nan(N,2);
    for k = 1:N
        cent(k,:) = props(labels(k)).Centroid;
    end
    xROI = cent(:,1);
    yROI = cent(:,2);
    
    % ---- Accumulate pairwise |Δt| and distance for this event ----
    % active ROIs with a valid Δt in this event
    idxAct = find(act & ~isnan(dt));   % indices of ROIs that fired in this event
    Mact   = numel(idxAct);
    
    if Mact >= 2
        % pairwise distance for these ROIs
        % use precomputed D_um if available and matches size
        if exist('D_um', 'var') && isequal(size(D_um,1), N)
            D_evt = D_um(idxAct, idxAct);   % µm
        else
            % compute from centroids (pixels)
            [X1, X2] = meshgrid(xROI(idxAct), xROI(idxAct));
            [Y1, Y2] = meshgrid(yROI(idxAct), yROI(idxAct));
            D_evt = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);   % pixels
            
            % convert to µm if pixel size is known
            if exist('pxSize_um', 'var')
                D_evt = D_evt * pxSize_um;
            end
        end

        % pairwise |Δt| for these ROIs (s)
        tAct = dt(idxAct);                % 1 x Mact
        [T1, T2] = meshgrid(tAct, tAct);
        lagMat = abs(T2 - T1);            % Mact x Mact, |Δt| in seconds

        % keep only unique pairs i<j
        maskTri = triu(true(Mact, Mact), 1);
        
        allLag  = [allLag;  lagMat(maskTri)];
        allDist = [allDist; D_evt(maskTri)];
    end

    % --- 2) Sort active ROIs by firing time ---
    idxAct = find(act);
    [~, ord] = sort(dt(idxAct));
    seq = idxAct(ord);
    
    xSeq = xROI(seq);
    ySeq = yROI(seq);
    tSeq = dt(seq);   % times used for BOTH dots and line
    
    % --- 3) Background grayscale image (TRUECOLOR RGB) ---
    g = mat2gray(avgImg);
    avgRGB = repmat(g, [1 1 3]);
    
    fig = figure('Color','w');
    ax = axes;
    image(ax, avgRGB);
    axis(ax, 'image', 'ij');
    hold(ax,'on');
    set(ax,'XTick',[],'YTick',[]);
    
    title('Wave path with time-colored line (grayscale background)');
    
    % --- 4) Draw ALL ROI boundaries (light gray) ---
    for k = 1:N
        BW = (maskL == labels(k));
        B  = bwboundaries(BW);
        for b = 1:numel(B)
            plot(ax, B{b}(:,2), B{b}(:,1), '-', ...
                 'Color',[0.6 0.6 0.6], 'LineWidth',0.5);
        end
    end
    
    % --- 5) Set colormap + time range ONCE ---
    cmap = jet(256);
    colormap(ax, cmap);
    
    tMin = min(tSeq);
    tMax = max(tSeq);
    if tMin == tMax
        tMin = tMin - 0.5;
        tMax = tMax + 0.5;
    end
    
    halfWinSec = binWin * binSec;   % = 1 * 0.5 = 0.5 sec
    totalWinSec = (2*binWin + 1) * binSec;   % = 1.5 sec
    
    caxis(ax, [-totalWinSec/2, totalWinSec/2]);   % i.e. [-0.75, 0.75]
    %caxis(ax, [tMin tMax]);   % this range defines colors for both dots and line
    
    % --- 6) TIME-COLORED PATH (using surface trick) ---
    X = [xSeq(:)'; xSeq(:)'];
    Y = [ySeq(:)'; ySeq(:)'];
    C = [tSeq(:)'; tSeq(:)'];   % same tSeq drives the color
    
    surface(ax, X, Y, zeros(size(X)), C, ...
        'FaceColor','none', ...
        'EdgeColor','interp', ...
        'LineWidth', 2);
    
    % --- 7) Active ROI dots, colored by tSeq via CData ---
    scatter(ax, xSeq, ySeq, 40, tSeq, ...  % CData = tSeq
            'filled', ...
            'MarkerEdgeColor','k', ...
            'LineWidth',0.5);
    
    % --- 8) Colorbar ---
    cb = colorbar(ax,'eastoutside');
    cb.Label.String = '\Delta t from event center (s)';
    
    % save as PDF
    outPDF = fullfile(outDir, sprintf('Event_%02d.pdf', kEvent));
    exportgraphics(fig, outPDF, 'ContentType','vector');
    fprintf('Saved %s\n', outPDF);

    % ===== Raw event movie with ROI circles based on event time =====
    persistFrames = 5;            % highlight lasts 5 frames
    framesEvt     = ev.frames;    % indices in S timebase (1..T)
    nEvtFrames    = numel(framesEvt);
    
    highlightTTL = zeros(1, N);   % TTL buffer per ROI
    
    % Map from S rows (1..T) to raw TIFF frames (1..nFrames)
    frameIndexRaw = (1:T) + nDrop;   % nDrop = 30 earlier
    
    % Output path
    %[tdir, tname, ~] = fileparts(outDir);
    outTiff = fullfile(outDir, sprintf('%s_event%02d_rawOverlay.tif', tname, kEvent));
    
    fprintf('Saving raw+overlay event frames to:\n  %s\n', outTiff);
    
    % precompute image size & a blank coordinate grid
    H = info(1).Height;
    W = info(1).Width;
    [Xg, Yg] = meshgrid(1:W, 1:H);   % X = cols, Y = rows
    
    radius   = 14;
    lineHalf = 0.7;                  % line thickness ~1 px
    maxVal   = intmax('uint16');     % 65535 for uint16
    
    for ii = 1:nEvtFrames
        fRow = framesEvt(ii);        % index into S / dFF / t
        fIdx = frameIndexRaw(fRow);  % corresponding raw TIFF frame index
    
        % ---- Read raw frame (uint16), KEEP AS GRAYSCALE ----
        frame = imread(tiffPath, fIdx, 'Info', info);   % H x W, uint16
        frameOut = frame;                               % copy we will draw on
    
        % ---- ROIs firing at THIS analysis frame (event-time masking) ----
        rFired = find(S(fRow, :) > 0);   % use fRow, not fIdx
    
        % Set TTL for newly fired ROIs
        highlightTTL(rFired) = persistFrames;
    
        % ROIs to show (TTL > 0)
        rShow = find(highlightTTL > 0);
    
        % ---- Draw circles directly on frameOut (no RGB, no scaling) ----
        for r = rShow
            cx = round(xROI(r));
            cy = round(yROI(r));
    
            % distance from this ROI center
            dist = sqrt((Xg - cx).^2 + (Yg - cy).^2);
    
            % ring mask: pixels within ~1 px of desired radius
            ringMask = abs(dist - radius) <= lineHalf;
    
            % set ring pixels to max intensity (white)
            frameOut(ringMask) = maxVal;
        end
    
        % ---- Write to multipage TIFF (no scaling, no compression), uint16 GRAY ----
        if ii == 1
            imwrite(frameOut, outTiff, 'tif', ...
                    'WriteMode','overwrite', 'Compression','none');
        else
            imwrite(frameOut, outTiff, 'tif', ...
                    'WriteMode','append', 'Compression','none');
        end
    
        % Decay TTL after writing
        highlightTTL = max(0, highlightTTL - 1);
    end
    
    fprintf('Done. Wrote %d frames for event %d with time-based ROI masking.\n', ...
            nEvtFrames, kEvent);
end

%% ===== Merge event PDFs into a single AllEvents.pdf =====

pdfList   = dir(fullfile(outDir, 'Event_*.pdf'));
fileNames = fullfile(outDir, {pdfList.name});

mergedPDF = fullfile(outDir, 'AllEvents.pdf');

if ~isempty(fileNames)
    append_pdfs_RZ(mergedPDF, fileNames{:});
    %fprintf('Merged %d event PDFs into:\n  %s\n', numel(fileNames), mergedPDF);
    fprintf('Merged %d event PDFs into:\n  %s\n', numel(fileNames), mergedPDF);
else
    fprintf('No Event_*.pdf files found in %s. Nothing to merge.\n', outDir);
end

%%
%% ============================================================
%  Distance vs. time-lag analysis across all events
%  (Are nearby neurons more synchronized?)
% =============================================================

if ~exist('maskL','var')
    warning('maskL not found in workspace. Distance analysis skipped.');
else
    fprintf('\n=== Distance vs |Δt| analysis ===\n');

    % ---- 1) ROI centroids from maskL (once) ----
    labels = unique(maskL);
    labels(labels == 0) = [];      % drop background
    nROI = numel(labels);

    props = regionprops(maskL, 'Centroid');
    cent  = nan(nROI, 2);

    for k = 1:nROI
        cent(k,:) = props(labels(k)).Centroid;
    end

    xROI = cent(:,1);   % column = x
    yROI = cent(:,2);   % row    = y

    % ---- 2) Pairwise distance matrix (pixels -> microns if you know px size) ----
    [X1, X2] = meshgrid(xROI, xROI);
    [Y1, Y2] = meshgrid(yROI, yROI);
    D_pix    = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);   % nROI x nROI

    % set your pixel size here (µm per pixel)
    pxSize_um = 0.9876865;    % *** TODO: replace with your actual value ***
    D_um      = D_pix * pxSize_um;

    % ---- 3) Accumulate distance vs |Δt| over all events ----
    allDist = [];
    allLag  = [];

    nEvents = numel(events);

    for kEvent = 1:nEvents
        ev          = events(kEvent);
        framesEvt   = ev.frames;        % indices into S / dFF / t
        frameCenter = ev.frame_center;  % center frame index

        dt_evt  = nan(1, nROI);   % Δt per ROI for THIS event
        act_evt = false(1, nROI);

        % map ROI indices in maskL space (labels 1..nROI) to S columns
        % here we assume ROI index r corresponds to S(:, r)
        for r = 1:min(nROI, size(S,2))
            idxLocal = find(S(framesEvt, r) == 1, 1, 'first');
            if ~isempty(idxLocal)
                act_evt(r) = true;
                frame_r    = framesEvt(idxLocal);
                dt_evt(r)  = (frame_r - frameCenter) / fs;   % seconds
            end
        end

        idxAct = find(act_evt & ~isnan(dt_evt));
        M      = numel(idxAct);

        if M < 2
            fprintf('Event %d: less than 2 active ROIs, skipping.\n', kEvent);
            continue;
        end

        % subselect distances and times for active ROIs
        D_sub   = D_um(idxAct, idxAct);   % M x M
        t_sub   = dt_evt(idxAct);         % 1 x M

        % pairwise |Δt|
        [T1, T2] = meshgrid(t_sub, t_sub);
        lagMat   = abs(T2 - T1);          % M x M

        % use upper triangle only (i < j)
        maskTri = triu(true(M, M), 1);

        distVec = D_sub(maskTri);
        lagVec  = lagMat(maskTri);

        allDist = [allDist; distVec(:)];
        allLag  = [allLag;  lagVec(:)];

        fprintf('Event %d/%d: used %d active ROIs -> %d pairs.\n', ...
                kEvent, nEvents, M, numel(distVec));
    end

    % ---- 4) Scatter: |Δt| vs distance (all events pooled) ----
    figure('Color','w');
    scatter(allDist, allLag, 8, 'k', 'filled'); 
    hold on;
    
    xlabel('Pairwise distance (µm)');
    ylabel('|Δt| between ROIs (s)');
    title('|Δt| vs distance across all events');
    grid on;
    
    % ---- Clean finite values ----
    valid = isfinite(allDist) & isfinite(allLag);
    x = allDist(valid);
    y = allLag(valid);
    
    % ---- Linear regression (y = a*x + b) ----
    p = polyfit(x, y, 1);         % p(1)=slope, p(2)=intercept
    xfit = linspace(min(x), max(x), 200);
    yfit = polyval(p, xfit);
    
    plot(xfit, yfit, 'r-', 'LineWidth', 2);
    
    % ---- R² (coefficient of determination) ----
    y_pred = polyval(p, x);
    SS_res = sum((y - y_pred).^2);
    SS_tot = sum((y - mean(y)).^2);
    Rsq   = 1 - SS_res/SS_tot;
    
    % ---- Spearman correlation ----
    [rhoSpearman, pSpearman] = corr(x, y, 'Type','Spearman');
    
    % ---- Format significance text ----
    if pSpearman < 1e-3
        pText = sprintf('p < 10^{-3}');
    else
        pText = sprintf('p = %.3f', pSpearman);
    end
    
    % ---- Annotation box on figure ----
    txt = {
        sprintf('Slope = %.3f s/µm', p(1))
        sprintf('R^2 = %.3f', Rsq)
        sprintf('%s', pText)
        };
    
    xpos = max(x) - 0.20*(max(x)-min(x));
    ypos = max(y) - 0.05*(max(y)-min(y));
    
    text(xpos, ypos, txt, ...
        'FontSize', 8, ...
        'BackgroundColor', 'w', ...
        'Margin', 4, ...
        'EdgeColor', 'k');
    
    %legend({'Data', 'Linear fit'}, 'Location','best');


    % % ---- 5) Binned median |Δt| vs distance ----
    % if ~isempty(allDist)
    %     % choose distance bins (in µm)
    %     maxDist = prctile(allDist, 99);   % robust upper bound
    %     nBins   = 15;
    %     edges   = linspace(0, maxDist, nBins+1);
    %     centers = 0.5 * (edges(1:end-1) + edges(2:end));
    % 
    %     medLag  = nan(nBins,1);
    %     nPerBin = nan(nBins,1);
    % 
    %     for b = 1:nBins
    %         inBin = allDist >= edges(b) & allDist < edges(b+1);
    %         nPerBin(b) = sum(inBin);
    %         if nPerBin(b) > 0
    %             medLag(b) = median(allLag(inBin));
    %         end
    %     end
    % 
    %     figure('Color','w');
    %     plot(centers, medLag, '-o', 'LineWidth',1.5);
    %     xlabel('Pairwise distance (µm)');
    %     ylabel('Median |Δt| (s)');
    %     title('Median |Δt| vs distance (binned)');
    %     grid on;

        % ---- 6) Spearman correlation ----
        [rhoSpearman, pSpearman] = corr(allDist, allLag, 'Type','Spearman');
        fprintf('Spearman corr(distance, |Δt|) = %.3f (p = %.3g)\n', ...
                rhoSpearman, pSpearman);
end

% %% ---- 4) Scatter: distance vs |Δt| (event-based pairs) ----
% figDist = figure('Color','w');
% 
% % 1) Keep only finite and within 0.5 s (your event width)
% maxDtPair = 0.15;                 % seconds
% valid = isfinite(allLag) & isfinite(allDist) & (allLag <= maxDtPair);
% 
% x = allLag(valid);     % |Δt| in seconds
% y = allDist(valid);    % distance (µm or pixels)
% 
% % quick sanity check
% fprintf('Using %d ROI pairs after filtering (|Δt| ≤ %.2f s).\n', numel(x), maxDtPair);
% 
% % 2) Scatter (raw data)
% scatter(x, y, 4, 'k', 'filled'); 
% hold on;
% 
% xlabel('|Δt| between ROIs (s)');
% ylabel('Pairwise distance (µm)');
% title('Distance vs |Δt| across all events');
% grid on;
% 
% % 3) Linear regression: distance = m * |Δt| + b
% p = polyfit(x, y, 1);
% xfit = linspace(min(x), max(x), 200);
% yfit = polyval(p, xfit);
% 
% plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
% % R²
% y_pred = polyval(p, x);
% SS_res = sum((y - y_pred).^2);
% SS_tot = sum((y - mean(y)).^2);
% Rsq   = 1 - SS_res/SS_tot;
% 
% % Spearman correlation
% [rhoSpearman, pSpearman] = corr(x, y, 'Type','Spearman');
% 
% if pSpearman < 1e-3
%     pText = 'p < 10^{-3}';
% else
%     pText = sprintf('p = %.3f', pSpearman);
% end
% 
% % 4) Bin |Δt| a bit to avoid near-duplicate x’s from float noise
% %    Here I round to nearest millisecond; tweak 3 -> 2 or 4 if you like
% x_round = round(x, 3);   % seconds, rounded
% 
% [uniqueLag, ~, binIdx] = unique(x_round);   % unique |Δt| bins
% nU = numel(uniqueLag);
% 
% meanDist_u   = accumarray(binIdx, y, [nU 1], @mean);
% medianDist_u = accumarray(binIdx, y, [nU 1], @median);
% 
% % 5) Plot mean/median on top
% plot(uniqueLag, meanDist_u,   'b.-', 'LineWidth', 2, 'MarkerSize', 12);
% plot(uniqueLag, medianDist_u, 'm.-', 'LineWidth', 2, 'MarkerSize', 12);
% 
% legend({'Data', 'Linear fit', 'Mean per |Δt|', 'Median per |Δt|'}, ...
%        'Location','northwest');
% 
% % 6) Annotation
% txt = {
%     sprintf('Slope = %.3f µm/s', p(1))
%     sprintf('R^2 = %.3f', Rsq)
%     pText
%     };
% 
% xpos = min(x) + 0.8*(max(x)-min(x));
% ypos = min(y) + 0.95*(max(y)-min(y));
% 
% text(xpos, ypos, txt, ...
%     'FontSize', 8, ...
%     'BackgroundColor', 'w', ...
%     'Margin', 4, ...
%     'EdgeColor', 'k');
% 
% % 7) SAVE
% outDist = fullfile(tdir, sprintf('%s_distVsLag.pdf', tname));
% exportgraphics(figDist, outDist, 'ContentType','vector');
% fprintf('Saved distance vs |Δt| figure:\n  %s\n', outDist);
% 
% %% ---- Scatter: event-train correlation vs distance (using R_peaks) ----
% 
% Ngood = numel(keepIdx);
% 
% % restrict to good ROIs
% Dg = D_um(keepIdx, keepIdx);      % Ngood x Ngood distances
% Rg = R_peaks(keepIdx, keepIdx);   % Ngood x Ngood correlations
% 
% % upper triangle (i<j) to avoid double counting and diagonals
% maskUT = triu(true(Ngood), 1);
% 
% distVec = Dg(maskUT);
% corrVec = Rg(maskUT);
% 
% % keep only finite
% valid = isfinite(distVec) & isfinite(corrVec);
% distVec = distVec(valid);
% corrVec = corrVec(valid);
% 
% figCorrDist = figure('Color','w');
% scatter(distVec, corrVec, 6, 'k', 'filled');
% hold on;
% xlabel('Pairwise distance (µm)');
% ylabel('Correlation of peak trains');
% title('Peak-train correlation vs distance');
% grid on;
% 
% % ---- optional: linear fit + Spearman ----
% [p_lin, ~] = polyfit(distVec, corrVec, 1);
% xfit = linspace(min(distVec), max(distVec), 200);
% yfit = polyval(p_lin, xfit);
% plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
% [rhoSpearman, pSpearman] = corr(distVec, corrVec, 'Type','Spearman');
% 
% if pSpearman < 1e-3
%     pText = 'p < 10^{-3}';
% else
%     pText = sprintf('p = %.3f', pSpearman);
% end
% 
% txt = {
%     sprintf('Slope = %.4f per µm', p_lin(1))
%     sprintf('Spearman \\rho = %.3f', rhoSpearman)
%     pText
%     };
% 
% xpos = min(distVec) + 0.70*(max(distVec) - min(distVec));
% ypos = min(corrVec) + 0.90*(max(corrVec) - min(corrVec));
% 
% text(xpos, ypos, txt, ...
%     'FontSize', 8, ...
%     'BackgroundColor', 'w', ...
%     'Margin', 4, ...
%     'EdgeColor', 'k');
% 
% % save
% outCorrDist = fullfile(tdir, sprintf('%s_peakCorrVsDist.pdf', tname));
% exportgraphics(figCorrDist, outCorrDist, 'ContentType','vector');
% fprintf('Saved peak-train correlation vs distance figure:\n  %s\n', outCorrDist);
% 
% %%
% %% ---- Scatter: event-train correlation vs distance (using R_peaks), thresholded ----
% 
% corrThresh = 0.15;   % <---- threshold for using a pair
% 
% Ngood = numel(keepIdx);
% 
% % restrict to good ROIs
% Dg = D_um(keepIdx, keepIdx);      % Ngood x Ngood distances
% Rg = R_peaks(keepIdx, keepIdx);   % Ngood x Ngood correlations
% 
% % --- Upper triangle mask (avoid duplicates) ---
% maskUT = triu(true(Ngood), 1);
% 
% distVec = Dg(maskUT);
% corrVec = Rg(maskUT);
% 
% % --- keep only pairs with correlation >= threshold ---
% valid = isfinite(distVec) & isfinite(corrVec) & (corrVec >= corrThresh);
% 
% distVec = distVec(valid);
% corrVec = corrVec(valid);
% 
% fprintf('Using %d ROI pairs with corr ≥ %.2f\n', numel(corrVec), corrThresh);
% 
% % =====================================================
% %         Scatter Plot (thresholded pairs only)
% % =====================================================
% 
% figCorrDist = figure('Color','w');
% scatter(distVec, corrVec, 10, 'k', 'filled');
% hold on;
% xlabel('Pairwise distance (µm)');
% ylabel(sprintf('Correlation (≥ %.2f)', corrThresh));
% title(sprintf('Peak-train correlation vs distance (corr ≥ %.2f)', corrThresh));
% grid on;
% 
% % ---- Optional linear fit ----
% if numel(corrVec) > 5
%     [p_lin, ~] = polyfit(distVec, corrVec, 1);
%     xfit = linspace(min(distVec), max(distVec), 200);
%     yfit = polyval(p_lin, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Spearman
%     [rhoSpearman, pSpearman] = corr(distVec, corrVec, 'Type','Spearman');
% 
%     if pSpearman < 1e-3
%         pText = 'p < 10^{-3}';
%     else
%         pText = sprintf('p = %.3f', pSpearman);
%     end
% 
%     txt = {
%         sprintf('Slope = %.4f per µm', p_lin(1))
%         sprintf('Spearman \\rho = %.3f', rhoSpearman)
%         pText
%     };
% 
%     xpos = min(distVec) + 0.70*(max(distVec) - min(distVec));
%     ypos = min(corrVec) + 0.90*(max(corrVec) - min(corrVec));
% 
%     text(xpos, ypos, txt, ...
%         'FontSize', 8, ...
%         'BackgroundColor', 'w', ...
%         'Margin', 4, ...
%         'EdgeColor', 'k');
% end
% 
% % save
% outCorrDist = fullfile(tdir, sprintf('%s_peakCorrVsDist_thresh%.2f.pdf', ...
%                                       tname, corrThresh));
% exportgraphics(figCorrDist, outCorrDist, 'ContentType','vector');
% 
% fprintf('Saved thresholded peak-train correlation vs distance figure:\n  %s\n', outCorrDist);
% 
% %%
% %% =====================================================
% %     Scatter: ALL pairs shown, correlation USING ONLY
% %     pairs with corr >= threshold
% % =====================================================
% 
% corrThresh = 0.1;   % threshold
% 
% Ngood = numel(keepIdx);
% 
% % restrict to good ROIs
% Dg = D_um(keepIdx, keepIdx);
% Rg = R_peaks(keepIdx, keepIdx);
% 
% % upper triangle only
% maskUT = triu(true(Ngood), 1);
% 
% distVec_all = Dg(maskUT);
% corrVec_all = Rg(maskUT);
% 
% % remove NaNs
% valid_all = isfinite(distVec_all) & isfinite(corrVec_all);
% distVec_all = distVec_all(valid_all);
% corrVec_all = corrVec_all(valid_all);
% 
% % threshold mask
% passMask = corrVec_all >= corrThresh;
% failMask = ~passMask;
% 
% fprintf('Total pairs = %d ; Passing threshold = %d ; Failing = %d\n', ...
%     numel(corrVec_all), sum(passMask), sum(failMask));
% 
% % ===================== PLOT =====================
% figCorrDist = figure('Color','w'); hold on;
% 
% % Plot failing (black)
% scatter(distVec_all(failMask), corrVec_all(failMask), ...
%     10, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
% 
% % Plot passing (blue)
% scatter(distVec_all(passMask), corrVec_all(passMask), ...
%     14, 'b', 'filled', 'MarkerFaceAlpha', 0.9);
% 
% xlabel('Pairwise distance (µm)');
% ylabel('Peak-train correlation');
% title(sprintf('Correlation vs Distance (blue: corr ≥ %.2f)', corrThresh));
% grid on;
% 
% % ===== Regression + Statistics: ONLY threshold-passing pairs =====
% distPass = distVec_all(passMask);
% corrPass = corrVec_all(passMask);
% 
% if numel(corrPass) > 5
%     % linear fit
%     [p_lin, ~] = polyfit(distPass, corrPass, 1);
%     xfit = linspace(min(distPass), max(distPass), 200);
%     yfit = polyval(p_lin, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Spearman ONLY on passing pairs
%     [rhoSpearman, pSpearman] = corr(distPass, corrPass, 'Type','Spearman');
% 
%     if pSpearman < 1e-3
%         pText = 'p < 10^{-3}';
%     else
%         pText = sprintf('p = %.3f', pSpearman);
%     end
% 
%     txt = {
%         sprintf('Slope = %.4f per µm', p_lin(1))
%         sprintf('Spearman \\rho = %.3f', rhoSpearman)
%         pText
%     };
% 
%     xpos = min(distPass) + 0.65*(max(distPass) - min(distPass));
%     ypos = min(corrPass) + 0.90*(max(corrPass) - min(corrPass));
% 
%     text(xpos, ypos, txt, ...
%         'FontSize', 8, ...
%         'BackgroundColor','w', ...
%         'Margin', 4, ...
%         'EdgeColor','k');
% end
% 
% % SAVE
% outCorrDist = fullfile(tdir, sprintf('%s_corrVsDist_ONLYpassThresh%.2f.pdf', ...
%                                      tname, corrThresh));
% exportgraphics(figCorrDist, outCorrDist, 'ContentType','vector');
% 
% fprintf('Saved corr-vs-dist figure using ONLY threshold-passing pairs for statistics:\n  %s\n', outCorrDist);

%%
%% ==== Event-based co-activation vs distance (good ROIs only) ====

Ngood    = numel(keepIdx);
S_good   = S(:, keepIdx);          % T x Ngood
Co       = zeros(Ngood, Ngood);    % # events where i & j co-active
Cnt_pair = zeros(Ngood, Ngood);    % # events where (i OR j) active

N_events = numel(events);
fprintf('Computing event-based co-activation for %d events...\n', N_events);

for e = 1:N_events
    framesEvt = events(e).frames;
    if isempty(framesEvt), continue; end

    % ROIs that fired at least once during this event window
    active = any(S_good(framesEvt, :) > 0, 1);   % 1 x Ngood logical
    idx    = find(active);

    Mact = numel(idx);
    if Mact < 2
        continue;
    end

    % for all unordered pairs within this event
    for a = 1:Mact
        i = idx(a);
        for b = a+1:Mact
            j = idx(b);

            % co-activation count
            Co(i,j) = Co(i,j) + 1;
            Co(j,i) = Co(j,i) + 1;
        end
    end

    % now update "either active" counts (Cnt_pair)
    % any ROI that was active in this event is eligible
    for a = 1:Mact
        i = idx(a);
        for b = a+1:Mact
            j = idx(b);
            Cnt_pair(i,j) = Cnt_pair(i,j) + 1;
            Cnt_pair(j,i) = Cnt_pair(j,i) + 1;
        end
    end
end

% Fraction of events in which the pair co-fired, *conditional* on at least
% one of them being active in that event:
CoFrac = Co %/ max(Cnt_pair, 1);   % Ngood x Ngood, in [0,1]

% Pairwise distances for good ROIs
Dg = D_um(keepIdx, keepIdx);     % Ngood x Ngood

% Upper triangle mask (i<j)
maskUT = triu(true(Ngood), 1);

% Flatten into vectors
distVec = Dg(maskUT);
coVec   = CoFrac(maskUT);

% Keep only finite values
valid = isfinite(distVec) & isfinite(coVec);
distVec = distVec(valid);
coVec   = coVec(valid);

%% ---- Scatter + mean distance per coactivation count ----
figCoDist = figure('Color','w'); hold on;

% raw scatter (each pair)
scatter(coVec, distVec, 6, [0.4 0.4 0.4], 'filled', ...
        'MarkerFaceAlpha', 0.6);

xlabel('Co-activation count');
ylabel('Pairwise distance (µm)');
title('Event-based co-activation vs distance');
grid on;

% ---- Unique coactivation counts present in the data ----
uniqueCounts = unique(coVec);          % e.g. [0 1 2 3 ...]
uniqueCounts = uniqueCounts(:);        % column vector

nC = numel(uniqueCounts);
meanDist_perCount = nan(nC,1);
semDist_perCount  = nan(nC,1);

for ii = 1:nC
    k = uniqueCounts(ii);              % coactivation count value

    idx = (coVec == k);                % ROI pairs with EXACTLY k coactivations
    if ~any(idx), continue; end

    vals = distVec(idx);               % all distances for this count
    meanDist_perCount(ii) = mean(vals);
    semDist_perCount(ii)  = std(vals)/sqrt(numel(vals));
end

% plot mean ± SEM versus count (red)
errorbar(uniqueCounts, meanDist_perCount, semDist_perCount, ...
    'ro-', 'LineWidth', 1, 'MarkerSize', 6, 'CapSize', 0, ...
    'DisplayName', 'Mean±SEM distance');

% ---- Regression of distance vs coactivation count ----
% model: distance = beta0 + beta1 * (coactivation count)
[p_lin, ~] = polyfit(coVec, distVec, 1);
xfit = linspace(min(coVec), max(coVec), 200);
yfit = polyval(p_lin, xfit);
plot(xfit, yfit, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Linear fit');

% Correlation metrics
[rhoSpearman, pSpearman] = corr(coVec, distVec, 'Type','Spearman');

if pSpearman < 1e-3
    pText = 'p < 10^{-3}';
else
    pText = sprintf('p = %.3f', pSpearman);
end

txt = {
    sprintf('Slope = %.4f µm per coactivation', p_lin(1))
    sprintf('Spearman \\rho = %.3f', rhoSpearman)
    pText
};

xpos = min(coVec) + 0.55*(max(coVec) - min(coVec));
ypos = min(distVec) + 0.90*(max(distVec) - min(distVec));

text(xpos, ypos, txt, ...
    'FontSize', 8, ...
    'BackgroundColor','w', ...
    'Margin', 4, ...
    'EdgeColor','k');

legend({'Pairs', 'Mean±SEM distance', 'Linear fit'}, ...
       'Location','northeast');

% save
outCoDist = fullfile(tdir, sprintf('%s_eventCoactivationVsDist.pdf', tname));
exportgraphics(figCoDist, outCoDist, 'ContentType','vector');
fprintf('Saved event-based co-activation vs distance:\n  %s\n', outCoDist);



%% ---- Final export: core analysis variables in ONE file ----
analysisFile = fullfile(outDir, sprintf('%s_analysisCore.mat', tname));

% collect variable names that actually exist in the workspace
varsToSave = {'dFF', 'S', 'F_denoised'};

if exist('pxSize_um', 'var'),  varsToSave{end+1} = 'pxSize_um';  end
if exist('xROI', 'var'),       varsToSave{end+1} = 'xROI';       end
if exist('yROI', 'var'),       varsToSave{end+1} = 'yROI';       end
if exist('D_um', 'var'),       varsToSave{end+1} = 'D_um';       end
if exist('allDist', 'var'),    varsToSave{end+1} = 'allDist';    end
if exist('allLag', 'var'),     varsToSave{end+1} = 'allLag';     end

% save using dynamic var list
save(analysisFile, varsToSave{:}, '-v7.3');

fprintf('Saved core analysis variables to:\n  %s\n', analysisFile);


%%
% c(τ) = (∑x(t)x(t+τ)) / (∑y(t)y(t+τ)∑x(t)y(t+τ)​)

function [cNorm, lags] = xcorr_nom(x, y, maxLag)
%XCORR_NOM  Cross-correlation *without* mean subtraction.
%   [cNorm, lags] = xcorr_nom(x, y, maxLag)
%
%   - x, y: column vectors (will be reshaped if needed)
%   - maxLag: max lag in samples
%   - cNorm: normalized correlation per lag
%
%   For each lag tau, we compute:
%     C_xy(tau) = sum_over_overlap x(t) * y(t+tau)
%     C_xx(tau) = sum_over_overlap x(t) * x(t+tau)
%     C_yy(tau) = sum_over_overlap y(t) * y(t+tau)
%
%   and return:
%     cNorm(tau) = C_xy(tau) / sqrt( C_xx(tau) * C_yy(tau) )

    x = x(:);
    y = y(:);

    if nargin < 2 || isempty(y)
        y = x;
    end

    if numel(x) ~= numel(y)
        error('x and y must have same length');
    end

    N = numel(x);

    if nargin < 3 || isempty(maxLag)
        maxLag = N - 1;
    end

    lags  = -maxLag:maxLag;
    nLags = numel(lags);

    cNorm = zeros(1, nLags);

    for k = 1:nLags
        tau = lags(k);

        if tau >= 0
            ix = 1:(N - tau);
            iy = (1 + tau):N;
        else
            ix = (1 - tau):N;
            iy = 1:(N + tau);
        end

        % overlapping segments
        x0   = x(ix);
        xTau = x(iy);
        y0   = y(ix);
        yTau = y(iy);

        % cross- and auto-correlations at this lag
        Cxy = sum(x0 .* yTau);   % cross
        Cxx = sum(x0 .* xTau);   % auto x
        Cyy = sum(y0 .* yTau);   % auto y

        denom = sqrt(Cxx * Cyy);

        if denom == 0
            cNorm(k) = 0;
        else
            cNorm(k) = Cxy / denom;
        end
    end
end

%%
function [keepIdx, excludeIdx] = selectROIsFromTraces(dFF, t)
% selectROIsFromTraces  Interactively select ROIs to exclude by clicking on traces.
%
%   [keepIdx, excludeIdx] = selectROIsFromTraces(dFF, t)
%   dFF: T x N  (time x ROI)
%   t  : T x 1  (time vector in seconds)
%
%   Usage:
%     - Plots all traces stacked.
%     - Left-click on a trace lane to toggle exclusion of that ROI.
%     - Right-click or press any key to finish.
%
%   Output:
%     keepIdx    - indices of ROIs to keep
%     excludeIdx - indices of ROIs to exclude

[T, N] = size(dFF);

if numel(t) ~= T
    error('Length of time vector t (%d) must match number of rows in dFF (%d).', numel(t), T);
end

% Vertical spacing so traces don't overlap too much
amp = prctile(abs(dFF(:)), 99);     % robust amplitude scale
if amp == 0
    amp = 1;
end
offset = 3 * amp;                   % spacing between traces

% Vertical center of each ROI lane
centers = offset * (0:N-1);         % 1 x N

% Stacked traces
yStack = dFF + centers;             % T x N

hFig = figure('Name','ΔF/F traces (click to exclude ROIs)');
hold on;

hLines = gobjects(1, N);
for i = 1:N
    hLines(i) = plot(t, yStack(:, i), 'Color', [0.3 0.3 0.8]); %#ok<AGROW>
end

set(gca, 'YDir', 'normal');
xlabel('Time (s)');
ylabel('ROI index (stacked)');
title({'ΔF/F traces (stacked)', ...
       'Left-click to toggle exclude; right-click or press key to finish'});

yticks(centers);
yticklabels(string(1:N));

excludeMask = false(1, N);
hHighlight  = gobjects(1, N);

disp('Left-click a trace lane to toggle exclusion. Right-click or press any key to finish.');

while ishandle(hFig)
    try
        [~, yClick, button] = ginput(1);
    catch
        break;  % figure likely closed
    end

    if isempty(button)
        break;  % key press
    end

    % Right or middle click -> finish
    if button ~= 1
        break;
    end

    % Find nearest ROI lane in y
    [~, roi] = min(abs(yClick - centers));
    if roi < 1 || roi > N
        continue;
    end

    % Toggle exclusion
    excludeMask(roi) = ~excludeMask(roi);

    if excludeMask(roi)
        % Add bold black overlay
        if ~isgraphics(hHighlight(roi))
            hHighlight(roi) = plot(t, yStack(:, roi), 'k-', 'LineWidth', 1.5);
        end
        fprintf('ROI %d -> EXCLUDE\n', roi);
    else
        % Remove overlay
        if isgraphics(hHighlight(roi))
            delete(hHighlight(roi));
            hHighlight(roi) = gobjects(1);
        end
        fprintf('ROI %d -> KEEP\n', roi);
    end
end

if ishandle(hFig)
    hold off;
end

excludeIdx = find(excludeMask);
keepMask   = true(1, N);
keepMask(excludeIdx) = false;
keepIdx    = find(keepMask);

fprintf('\nSelected %d ROIs to EXCLUDE.\n', numel(excludeIdx));
fprintf('Keeping %d ROIs.\n', numel(keepIdx));
end

%%
% function [F_denoised, dFF_denoised, spikes, baseline] = oasis_denoise_and_dff(F)
% % oasis_denoise_and_dff  OASIS-based denoising + ΔF/F using Python's default AR model.
% %
% %   [F_denoised, dFF_denoised, spikes, baseline] = oasis_denoise_and_dff(F)
% %
% %   Inputs:
% %     F  : T x N matrix (raw fluorescence per ROI)
% %
% %   Outputs:
% %     F_denoised   : T x N, OASIS denoised calcium trace C for each ROI
% %     dFF_denoised : T x N, ΔF/F computed from C and OASIS baseline b
% %     spikes       : T x N, OASIS spike trains (S)
% %     baseline     : 1 x N, OASIS baseline b for each ROI
% %
% %   Requires:
% %     - pyenv pointing to a Python env with:
% %         conda install -c conda-forge oasis-deconv
% 
%     % --- Make sure Python OASIS is importable ---
%     try
%         np     = py.importlib.import_module('numpy');
%         oasisF = py.importlib.import_module('oasis.functions');
%     catch ME
%         error(['Python OASIS not available.\n' ...
%                'Make sure pyenv points to an environment where you ran:\n' ...
%                '  conda install -c conda-forge oasis-deconv\n\n' ...
%                'Original error:\n%s'], ME.message);
%     end
% 
%     [T, N] = size(F);
% 
%     F_denoised   = zeros(T, N);
%     dFF_denoised = zeros(T, N);
%     spikes       = zeros(T, N);
%     baseline     = zeros(1, N);
% 
%     for i = 1:N
%         y = double(F(:, i));           % T x 1
%         y_np = np.array(y);            % numpy array
% 
%         try
%             % Let OASIS choose its own AR model for this trace
%             out = oasisF.deconvolve(y_np, pyargs('penalty', int32(0)));
% 
%             % out{1} = C, out{2} = S, out{3} = b
%             c_py = out{1};
%             s_py = out{2};
%             b_py = out{3};
% 
%             c = double(py.array.array('d', c_py.tolist()));   % denoised trace
%             s = double(py.array.array('d', s_py.tolist()));   % spikes
%             b = double(b_py);                                 % baseline
%         catch ME
%             warning('OASIS failed on ROI %d: %s\nUsing fallback (no deconvolution).', ...
%                     i, ME.message);
% 
%             % Fallback: no deconvolution, simple baseline from raw F
%             c = y;
%             s = zeros(size(y));
% 
%             finiteY = y(isfinite(y));
%             if isempty(finiteY)
%                 b = 1e-6;
%             else
%                 b = prctile(finiteY, 10);  % 10th percentile as baseline
%             end
%         end
% 
%         % baseline guard
%         if ~isfinite(b) || b <= 0
%             b = max(1e-6, prctile(c, 10));
%         end
% 
%         F_denoised(:, i)   = c(:);
%         dFF_denoised(:, i) = (c - b) / b;
%         spikes(:, i)       = s(:);
%         baseline(1, i)     = b;
%     end
% end
% 
% 
% % function [F_denoised, dFF_denoised, spikes, baseline] = oasis_denoise_and_dff(F, AR)
% % % oasis_denoise_and_dff  OASIS-based denoising + ΔF/F on F (Python OASIS only).
% % %
% % %   [F_denoised, dFF_denoised, spikes, baseline] = oasis_denoise_and_dff(F, AR)
% % %
% % %   Inputs:
% % %     F  : T x N matrix (raw fluorescence per ROI)
% % %     AR : optional AR coefficients for OASIS, e.g. [1.7 -0.712]
% % %          If omitted or empty, defaults to [1.7 -0.712].
% % %
% % %   Outputs:
% % %     F_denoised   : T x N, OASIS denoised calcium trace C for each ROI
% % %     dFF_denoised : T x N, ΔF/F computed from C and OASIS baseline b
% % %     spikes       : T x N, OASIS spike trains (S)
% % %     baseline     : 1 x N, OASIS baseline b for each ROI
% % %
% % %   Requires:
% % %     - pyenv() pointing to Python env with oasis-deconv installed.
% % 
% %     if nargin < 2 || isempty(AR)
% %         AR = [1.7 -0.712];   % default AR(2)
% %     end
% % 
% %     % --- Make sure Python OASIS is importable ---
% %     try
% %         np     = py.importlib.import_module('numpy');
% %         oasisF = py.importlib.import_module('oasis.functions');
% %     catch ME
% %         error(['Python OASIS not available.\n' ...
% %                'Make sure pyenv points to an environment where you ran:\n' ...
% %                '  conda install -c conda-forge oasis-deconv\n\n' ...
% %                'Original error:\n%s'], ME.message);
% %     end
% % 
% %     [T, N] = size(F);
% % 
% %     F_denoised   = zeros(T, N);
% %     dFF_denoised = zeros(T, N);
% %     spikes       = zeros(T, N);
% %     baseline     = zeros(1, N);
% % 
% %     % prepare AR coefficients for Python AR(2) call
% %     g_tuple = py.tuple(num2cell(AR(:).'));
% % 
% %     for i = 1:N
% %         y = double(F(:, i));   % T x 1
% % 
% %         % handle NaNs/Infs quickly
% %         finiteY = y(isfinite(y));
% %         if isempty(finiteY) || std(finiteY) < 1e-8
% %             % trivial trace: nothing to deconvolve
% %             c = y;
% %             s = zeros(size(y));
% %             b = 1e-6;
% %         else
% %             % small pre-normalization to avoid stupid scaling issues
% %             y_clean = y;
% %             y_clean(~isfinite(y_clean)) = median(finiteY);
% %             % shift to start around baseline-ish positive
% %             y_clean = y_clean - min(finiteY) + 1e-3;
% %             y_np = np.array(y_clean);
% % 
% %             try
% %                 % ---- First try: AR(2) with fixed coefficients ----
% %                 out = oasisF.deconvolve(y_np, g_tuple, pyargs('penalty', int32(0)));
% % 
% %             catch ME1
% %                 warning('OASIS AR2 failed on ROI %d: %s\nRetrying with default AR model...', ...
% %                         i, ME1.message);
% %                 try
% %                     % ---- Second try: let OASIS pick AR model by itself ----
% %                     out = oasisF.deconvolve(y_np, pyargs('penalty', int32(0)));
% %                 catch ME2
% %                     warning(['OASIS default model also failed on ROI %d: %s\n' ...
% %                              'Using fallback (no deconvolution).'], ...
% %                              i, ME2.message);
% %                     c = y_clean;
% %                     s = zeros(size(y_clean));
% %                     % robust baseline from low percentile
% %                     b = prctile(finiteY, 10);
% %                     if ~isfinite(b) || b <= 0
% %                         b = 1e-6;
% %                     end
% % 
% %                     F_denoised(:, i)   = c(:);
% %                     dFF_denoised(:, i) = (c - b) / b;
% %                     spikes(:, i)       = s(:);
% %                     baseline(1, i)     = b;
% %                     continue;  % move on to next ROI
% %                 end
% %             end
% % 
% %             % If we get here, at least one OASIS call succeeded
% %             c_py = out{1};
% %             s_py = out{2};
% %             b_py = out{3};
% % 
% %             c = double(py.array.array('d', c_py.tolist()));   % denoised trace
% %             s = double(py.array.array('d', s_py.tolist()));   % spikes
% %             b = double(b_py);                                 % baseline
% % 
% %             % baseline guard (avoid zero or NaN)
% %             if ~isfinite(b) || b <= 0
% %                 b = max(1e-6, prctile(c, 10));
% %             end
% %         end
% % 
% %         F_denoised(:, i)   = c(:);
% %         dFF_denoised(:, i) = (c - b) / b;
% %         spikes(:, i)       = s(:);
% %         baseline(1, i)     = b;
% %     end
% % end

%%
% AR(1), instantaneous rise, not enforcing monotonicity. bound decay time.

% assume exponential decay: c(t)=c(0)e^(−t/τ)
% g = e^(−Δt/τ)
% τ=−Δt/log(g)
% at 30fps, g=0.95 --> τ=650ms​
function [F_denoised, dFF_denoised, spikes, baseline] = oasis_denoise_and_dff_AR1(F, g)
    % F: T x N
    % g: scalar AR(1) coefficient (0 < g < 1)

    try
        np     = py.importlib.import_module('numpy');
        oasisF = py.importlib.import_module('oasis.functions');
    catch ME
        error(['Python OASIS not available.\n' ...
               'Make sure pyenv points to an environment where you ran:\n' ...
               '  conda install -c conda-forge oasis-deconv\n\n' ...
               'Original error:\n%s'], ME.message);
    end

    [T, N] = size(F);
    F_denoised   = zeros(T, N);
    dFF_denoised = zeros(T, N);
    spikes       = zeros(T, N);
    baseline     = zeros(1, N);

    % AR(1) tuple for Python
    g_tuple = py.tuple(num2cell(g(:).'));  % 1-element tuple

    for i = 1:N
        y = double(F(:, i));
        y_np = np.array(y);

        try
            % AR(1) with fixed g
            out = oasisF.deconvolve(y_np, g_tuple, pyargs('penalty', int32(0)));
            c_py = out{1};
            s_py = out{2};
            b_py = out{3};

            c = double(py.array.array('d', c_py.tolist()));
            s = double(py.array.array('d', s_py.tolist()));
            b = double(b_py);

        catch ME
            warning('OASIS AR(1) failed on ROI %d: %s\nUsing fallback.', i, ME.message);
            c = y;
            s = zeros(size(y));

            finiteY = y(isfinite(y));
            if isempty(finiteY)
                b = 1e-6;
            else
                b = prctile(finiteY, 10);
            end
        end

        if ~isfinite(b) || b <= 0
            b = max(1e-6, prctile(c, 10));
        end

        F_denoised(:, i)   = c(:);
        dFF_denoised(:, i) = (c - b) / b;
        spikes(:, i)       = s(:);
        baseline(1, i)     = b;
    end
end
