% IO_population_seq_connectivity_260607.m
% =========================================================================
%  INFERIOR OLIVE (IO) population ENSEMBLE-EVENT viewer (minimal).
%
%  Three things only -- all from GLOBAL dF/F + a calcium threshold (no spike
%  sorting / deconvolution needed):
%    1. POPULATION ACTIVITY TRACE
%         P_pop(t) = mean z-scored dF/F across all ROIs.
%    2. ENSEMBLE EVENTS
%         Detect peaks in P_pop (synchronous co-activation events), then
%         cluster the per-event across-ROI activity vectors into ENSEMBLES.
%    3. PER-ENSEMBLE ACTIVATION SEQUENCE ON THE AVG PROJECTION
%         For each ensemble, every cell that crosses the calcium threshold
%         within a +/- window around the event peak is ordered by its first
%         (rising-edge) crossing and its centroid connected by a line over the
%         average projection, colored by activation time.
%
%  Everything else (PCA, connectivity, modules, sequence stats, split-half,
%  lead-lag, wave fit, time tuning, shuffle controls) has been removed.
%  Spike sorting (ca_spike_data) is NOT used.
%
%  Convention: rows = time frames, columns = ROIs/neurons  (T x N).
%  Dependencies: detect_session_fps.m, helper.dFF_RZ, Image Processing +
%                Statistics toolboxes (regionprops, bwboundaries, kmeans).
% =========================================================================

clear; close all; clc;

%% ========================= USER PARAMETERS =============================
% inputPath: a session FOLDER with *cpSAM_output.mat (+ *ca_spike_data.mat,
% *_AVG_for_CP.tif), or a single .mat with dFF [T x N] (+ spikeProb/maskL).

%inputPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\airpuff_redo\Right_POandDAO_z80_airpuff_00001";

inputPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\sound\DAO_Z55_sound_5khz_00001";

%"D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\sound\Right_POandDAO_z80_sound_5k_00001";

% --- acquisition ---
fps_fallback   = 30;        % Hz, used only if detect_session_fps fails
BaselineWinSec = 20;        % dF/F sliding-median window (s) for helper.dFF_RZ
nDrop          = 30;         % extra leading frames to drop (cpSAM F already trimmed)
PixelSizeBase  = 1.7778;    % um/px at zoom=1 (fallback if no metadata)
pixelSize_override = NaN;   % um/px; set e.g. 0.4444551 to force, else NaN=auto

% --- population signal: GLOBAL dF/F (no spike sorting) ---
seq_smooth_sec  = 0.10;     % light Gaussian smoothing of z dF/F (s)
caThreshSD      = 2.0;      % calcium event threshold: z dF/F crossing (SD). Uniform across
                            % cells (signal is z-scored) -> used for activation timing + raster

% --- ensemble event detection ---
popThresh_kSD   = 2.0;      % event peak threshold = mean(P_pop) + k*SD(P_pop)
popThresh_pct   = [];       % alt: percentile threshold (e.g. 98). [] = use kSD
evMinDist_sec   = .25;      % min separation between events (s)
winSec          = 0.25;      % event window half-width (-0.3..+0.3 s)
minActiveNeurons = 10;       % drop events with fewer active cells than this
nEnsembles      = [];       % # ensembles. [] = auto-select (silhouette); or set a number
maxEns          = 2;        % upper bound on k when auto-selecting

% --- avg-projection display ---
clip_pct        = [0.5 99.9];  % intensity clip percentiles for display
gamma_disp      = 0.6;         % display gamma (<1 brightens midtones)

% --- 2P epoch stacks of ensemble activity (16-bit TIFF, original range) ---
saveEpochTiffs  = true;        % write per-ensemble multipage 16-bit TIFF from the MC movie
vidPreSec       = 0.5;         % epoch: seconds before each event peak
vidPostSec      = 0.5;         % epoch: seconds after each event peak
vidGapFrames    = 5;           % black separator frames between concatenated events

% --- output ---
doSave          = true;
outTag          = 'IO_ensemble_events_260607';
% ========================================================================

%% ========================= PATH SETUP ==================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 8);
rng(1, 'twister');                                % reproducible k-means

%% ============================ LOAD =====================================
DATA = load_io_data(char(inputPath), struct( ...
    'fps_fallback', fps_fallback, 'BaselineWinSec', BaselineWinSec, ...
    'nDrop', nDrop, 'PixelSizeBase', PixelSizeBase, ...
    'pixelSize_override', pixelSize_override));

dFF      = DATA.dFF;
fps      = DATA.fps;
maskL    = DATA.maskL;
roiXY    = DATA.roiXY;
sessName = DATA.name;
folder   = DATA.folder;
[T, N]   = size(dFF);

% Population signal = z-scored dF/F (per ROI), lightly smoothed. No spike
% sorting: events are detected from the global mean and per-cell activation is
% a calcium-threshold crossing. (Detected spikes here are just raw-dF/F
% threshold peaks anyway, so a dF/F threshold loses nothing and captures the
% rising-edge onset better; uniform indicator lag preserves activation order.)
Z = zscore(dFF, 0, 1);                            % [T x N] per-ROI z dF/F
seqSig = smooth_cols_gauss(Z, max(1, round(seq_smooth_sec*fps)));
zThr = caThreshSD;                                % uniform z threshold (signal is z-scored)
sigLabel = sprintf('mean z dF/F (smoothed %.2f s), event thr = %.1f SD', seq_smooth_sec, caThreshSD);
fprintf('\n=== IO session: %s ===\n', sessName);
fprintf('T=%d (%.1f s) | N=%d ROIs | fps=%.2f | %.4f um/px\n', T, T/fps, N, fps, DATA.px_um);
fprintf('Population signal: %s\n', sigLabel);

% output folder
outDir = fullfile(folder, outTag);
if doSave && ~exist(outDir, 'dir'), mkdir(outDir); end

t_img = (0:T-1)' / fps;

% avg projection (backdrop for the ensemble sequence plots)
hasFOV = ~isempty(DATA.avgImg) && ~isempty(roiXY);
if hasFOV
    avgDisp = display_image(DATA.avgImg, clip_pct, gamma_disp);
else
    avgDisp = [];
    warning('No avg projection / centroids -> sequence plots use centroid scatter.');
end

%% ================ 1+2: POPULATION TRACE + ENSEMBLE EVENTS ==============
% Global population activity = mean z dF/F across ROIs; peaks above threshold
% = synchronous (global ensemble) events.
P_pop = mean(seqSig, 2);                          % [T x 1]
if ~isempty(popThresh_pct)
    popThr = prctile(P_pop, popThresh_pct);
    thrLabel = sprintf('%gth pct', popThresh_pct);
else
    popThr = mean(P_pop) + popThresh_kSD * std(P_pop);
    thrLabel = sprintf('mean+%.1fSD', popThresh_kSD);
end
evMinDist = max(1, round(evMinDist_sec * fps));
[~, evFrame] = findpeaks(P_pop, 'MinPeakHeight', popThr, 'MinPeakDistance', evMinDist);

% keep events whose full +/- window fits inside the recording
preF  = round(winSec * fps);
postF = round(winSec * fps);
Lwin  = preF + postF + 1;
tauv  = ((0:Lwin-1) - preF)' / fps;              % time relative to peak (s), [-win..win]
evFrame = evFrame(evFrame - preF >= 1 & evFrame + postF <= T);

% Activation time per event/ROI = first time the cell's z dF/F crosses the
% calcium threshold (zThr) in the window -- the rising-edge onset. If it
% crosses several times, only the first counts. A cell is "active" iff it
% crosses threshold in the window (no spike sorting needed).
nEv0 = numel(evFrame);
Lat   = nan(N, nEv0);                             % first-crossing time (s) relative to event peak
EvVec = nan(nEv0, N);                             % mean z dF/F per event/ROI (clustering feature)
for e = 1:nEv0
    idx = (evFrame(e)-preF) : (evFrame(e)+postF);
    W = seqSig(idx, :)';                          % [N x Lwin] z dF/F
    EvVec(e, :) = mean(W, 2)';
    act = max(W, [], 2) >= zThr;                  % active = crosses calcium threshold
    for i = find(act)'
        f = find(W(i, :) >= zThr, 1, 'first');    % first rising-edge crossing
        Lat(i, e) = tauv(f);
    end
end

% drop events with too few active cells
nAct_ev = sum(isfinite(Lat), 1)';
keep = nAct_ev >= minActiveNeurons;
evFrame = evFrame(keep);  Lat = Lat(:, keep);  EvVec = EvVec(keep, :);
nEv = numel(evFrame);

% Global event strength, tracked SEPARATELY from the (amplitude-blind) cluster
% features so ensembles can be compared by amplitude and #active cells, not
% just spatial pattern.
eventAmp = mean(EvVec, 2);                        % [nEv x 1]
nActive_ev = sum(isfinite(Lat), 1)';             % [nEv x 1] active cells per event
fprintf('Ensemble events: %d (threshold %s = %.3f, window +/-%.2f s)\n', ...
    nEv, thrLabel, popThr, winSec);
assert(nEv >= 2, 'Too few events (%d) to form ensembles.', nEv);

% cluster events into ensembles by their across-ROI activity vector.
% The number of ensembles is data-driven: if nEnsembles is empty we pick k
% by the best mean silhouette over 2..maxEns; set nEnsembles to force a value.
EvVec(~isfinite(EvVec)) = 0;
Xz = zscore(EvVec, 0, 2);                          % [nEv x N], per-event z-scored
if isempty(nEnsembles) || nEnsembles <= 0
    kHi = min(maxEns, nEv - 1);
    if kHi >= 2
        try
            evac = evalclusters(Xz, 'kmeans', 'silhouette', 'KList', 2:kHi);
            kEns = evac.OptimalK;
        catch
            kEns = min(4, kHi);
        end
        if ~isfinite(kEns) || kEns < 1, kEns = min(4, kHi); end
        fprintf('Auto-selected k=%d ensembles (silhouette over 2..%d)\n', kEns, kHi);
    else
        kEns = 1;
    end
else
    kEns = min(nEnsembles, nEv);
end
if kEns > 1
    ensID = kmeans(Xz, kEns, 'Replicates', 20, 'MaxIter', 300, 'EmptyAction', 'singleton');
else
    ensID = ones(nEv, 1);
end
kEns = max(ensID);
cmapE = lines(kEns);
fprintf('Ensembles: k=%d, event counts = [%s]\n', kEns, ...
    num2str(arrayfun(@(c) sum(ensID==c), 1:kEns)));

% --- Figure 1: exact population activity ---
%  TOP   : global population activity P_pop(t) = mean z dF/F + threshold + events
%  BOTTOM: full raster -- every ROI's calcium-threshold crossings over time (so
%          you can see what is actually happening, not just the population mean).
%          Event times shown as faint ensemble-colored vertical lines.
fig1 = figure('Color','w','Name','exact population activity + ensemble events', ...
    'Units','centimeters','Position',[2 2 32 22]);

ax1 = subplot(4,1,1); hold(ax1, 'on');
area(ax1, t_img, P_pop, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none', 'FaceAlpha',0.6);
yline(ax1, popThr, 'k--', 'LineWidth', 1);
yl = ylim(ax1);
for c = 1:kEns
    ii = ensID == c;
    plot(ax1, t_img(evFrame(ii)), repmat(yl(2)*0.97, sum(ii), 1), 'v', ...
        'Color', cmapE(c,:), 'MarkerFaceColor', cmapE(c,:), 'MarkerSize', 6);
end
ylabel(ax1, 'pop. activity (mean z dF/F)'); xlim(ax1, [0 t_img(end)]);
set(ax1, 'XTickLabel', []);
legend(ax1, [{'P_{pop}(t)','threshold'}, arrayfun(@(c) sprintf('ens %d',c), ...
    1:kEns, 'UniformOutput',false)], 'Location','eastoutside','Box','off','FontSize',7);
title(ax1, sprintf('Exact population activity + %d ensemble events | %s', ...
    nEv, strrep(sessName,'_','\_')), 'Interpreter','tex');

ax2 = subplot(4,1,2:4); hold(ax2, 'on');
% faint ensemble-colored event lines behind the raster
for c = 1:kEns
    for e = find(ensID == c)'
        xline(ax2, t_img(evFrame(e)), '-', 'Color', [cmapE(c,:) 0.25], 'LineWidth', 0.5);
    end
end
% raster of rising-edge calcium-threshold crossings (dF/F events, no spike sort)
crossBin = [false(1,N); (seqSig(2:end,:) >= zThr) & (seqSig(1:end-1,:) < zThr)];
[fr, roi] = find(crossBin);
plot(ax2, (fr-1)/fps, roi, '.', 'Color', [0 0 0], 'MarkerSize', 3);
ylabel(ax2, 'ROI (\DeltaF/F threshold crossings)');
xlabel(ax2, 'time (s)'); xlim(ax2, [0 t_img(end)]); ylim(ax2, [0 N+1]);
linkaxes([ax1 ax2], 'x');
savefig_local(fig1, outDir, '01_population_activity_raster', doSave);

%% ============ 3: PER-ENSEMBLE ACTIVATION SEQUENCE ON AVG PROJECTION =====
% For each ensemble, order the ACTIVE cells by their mean first-activation
% time and draw a line connecting their centroids in that temporal order,
% over the avg projection. Nodes are colored by activation time; the open
% circle marks the first-activated cell. A line that sweeps smoothly across
% space = a propagating sequence; a tight tangle = near-synchronous.
latClim = [-winSec winSec];
ensLat  = nan(N, kEns);                           % mean first-activation time per ROI/ensemble
nCols = min(kEns, 3);  nRows = ceil(kEns / nCols);
fig2 = figure('Color','w','Name','per-ensemble activation sequence on avg projection', ...
    'Units','centimeters','Position',[2 2 9*nCols 8*nRows]);
for c = 1:kEns
    ii = find(ensID == c);
    mLat = mean(Lat(:, ii), 2, 'omitnan');        % [N x 1], NaN = never active
    ensLat(:, c) = mLat;
    nActive = nnz(isfinite(mLat));

    ax = subplot(nRows, nCols, c);
    draw_sequence_on_fov(ax, hasFOV, avgDisp, roiXY, mLat, parula(256), latClim);
    title(ax, sprintf('ens %d  (n_{ev}=%d, active=%d, amp=%.3g)', ...
        c, numel(ii), nActive, mean(eventAmp(ii))), ...
        'Interpreter','tex', 'Color', cmapE(c,:));
end
sgtitle(fig2, sprintf('Per-ensemble activation sequence (first-spike order) | %s', ...
    strrep(sessName,'_','\_')), 'Interpreter','tex');
savefig_local(fig2, outDir, '02_ensemble_activation_sequence', doSave);

%% ============ 2P TIFF EPOCH STACKS OF ENSEMBLE ACTIVITY ================
% One multipage 16-bit TIFF per ensemble: the motion-corrected 2P movie cut
% to a [-vidPreSec .. +vidPostSec] window around each of that ensemble's
% event peaks, concatenated (black gap between events). Frames are the
% ORIGINAL 16-bit data, untouched -- no scaling, no clip, no gamma -- so the
% full dynamic range is preserved (LZW lossless). Open in ImageJ/Fiji and set
% playback speed there. Saved in a subfolder of the data folder.
if saveEpochTiffs
    mc = dir(fullfile(folder, '*_MC_MC.tif'));
    if isempty(mc)
        mc = dir(fullfile(folder, '*MC*.tif'));
        mc = mc(~contains({mc.name}, 'AVG', 'IgnoreCase', true));
    end
    if isempty(mc)
        fprintf('No MC movie (*_MC_MC.tif) found -> skipping TIFF epochs.\n');
    else
        fprintf('Loading 2P movie for TIFF epochs: %s\n', mc(1).name);
        V = tiffreadVolume(fullfile(mc(1).folder, mc(1).name));   % [H x W x nMovie]
        nMovie = size(V, 3);
        offset = max(0, nMovie - T);          % movie may retain leading TossFrames
        vpre  = round(vidPreSec  * fps);
        vpost = round(vidPostSec * fps);

        vidDir = fullfile(outDir, 'ensemble_2p_epochs');
        if ~exist(vidDir, 'dir'), mkdir(vidDir); end

        for c = 1:kEns
            evc = find(ensID == c);
            frameIdxList = cell(numel(evc), 1);
            for j = 1:numel(evc)
                cF = evFrame(evc(j)) + offset;
                fi = (cF - vpre) : (cF + vpost);
                frameIdxList{j} = fi(fi >= 1 & fi <= nMovie);
            end
            tpath = fullfile(vidDir, sprintf('ensemble%d_2p_epochs.tif', c));
            nFr = write_epoch_tiff(V, frameIdxList, tpath, vidGapFrames);
            fprintf('  ensemble %d: %d events, %d frames -> %s\n', c, numel(evc), nFr, tpath);
        end
        clear V
        fprintf('16-bit TIFF epoch stacks (original dynamic range) saved in:\n  %s\n', vidDir);
    end
end

%% ============================ SAVE =====================================
if doSave
    results = struct();
    results.sessName = sessName;  results.folder = folder;
    results.fps = fps;  results.px_um = DATA.px_um;  results.T = T;  results.N = N;
    results.params = struct('seq_smooth_sec',seq_smooth_sec,'caThreshSD',caThreshSD, ...
        'popThresh_kSD',popThresh_kSD,'popThresh_pct',popThresh_pct, ...
        'evMinDist_sec',evMinDist_sec,'winSec',winSec, ...
        'minActiveNeurons',minActiveNeurons, ...
        'nEnsembles',nEnsembles,'BaselineWinSec',BaselineWinSec);
    results.P_pop = P_pop;  results.popThr = popThr;
    results.evFrame = evFrame;  results.tauv = tauv;
    results.Lat = Lat;  results.EvVec = EvVec;
    results.eventAmp = eventAmp;  results.nActive_ev = nActive_ev;
    results.ensID = ensID;  results.ensLat = ensLat;
    % per-ensemble summary: count, mean amplitude, mean #active cells
    results.ens_summary = arrayfun(@(c) struct('ensemble',c, ...
        'nEvents',sum(ensID==c), 'meanAmp',mean(eventAmp(ensID==c)), ...
        'meanActive',mean(nActive_ev(ensID==c))), 1:kEns);
    results.roiXY = roiXY;
    save(fullfile(outDir, [outTag '_results.mat']), 'results', '-v7.3');
    fprintf('Saved 2 figures + results to:\n  %s\n', outDir);
end

%% ========================= LOCAL FUNCTIONS ============================
function draw_sequence_on_fov(ax, hasFOV, avgDisp, XY, lat, cmap, clim)
% Draw the activation SEQUENCE over the avg projection: active cells (finite
% lat + finite centroid) ordered by first-activation time and connected by a
% line through their centroids; each line segment & node colored by activation
% time; an open circle marks the first-activated cell. XY = ROI centroids (px).
if hasFOV
    imshow(repmat(avgDisp, [1 1 3]), 'Parent', ax, 'InitialMagnification','fit');
    hold(ax, 'on');
else
    hold(ax, 'on'); set(ax, 'YDir', 'reverse');
end
act = find(isfinite(lat) & all(isfinite(XY), 2));
if numel(act) >= 1
    [~, ord] = sort(lat(act), 'ascend');
    seq = act(ord);
    xy  = XY(seq, :);
    rgb = val2rgb(lat(seq), cmap, clim);
    for k = 1:size(xy,1)-1                         % connect in temporal order
        plot(ax, xy(k:k+1,1), xy(k:k+1,2), '-', 'Color', rgb(k,:), 'LineWidth', 1.2);
    end
    scatter(ax, xy(:,1), xy(:,2), 40, lat(seq), 'filled', 'MarkerEdgeColor', 'k');
    plot(ax, xy(1,1), xy(1,2), 'o', 'MarkerSize', 13, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5);                         % start marker
end
colormap(ax, cmap); caxis(ax, clim);
axis(ax, 'image');
if hasFOV, axis(ax, 'off'); else, xlabel(ax,'x (px)'); ylabel(ax,'y (px)'); end
cb = colorbar(ax); cb.Label.String = 'first-activation time (s)';
end

function nFr = write_epoch_tiff(V, frameIdxList, tpath, gapFrames)
% Write concatenated movie epochs to a multipage 16-bit TIFF, ORIGINAL
% dynamic range untouched (no scaling, no clip, no gamma; LZW lossless).
% V is [H x W x nMovie]; frameIdxList{e} = frame indices for event e.
% gapFrames black frames separate consecutive events. Returns frame count.
[H, W, ~] = size(V);
blk = zeros(H, W, 'like', V);
if isfile(tpath), delete(tpath); end          % fresh stack (no stale append)
nFr = 0;  first = true;
for ev = 1:numel(frameIdxList)
    fi = frameIdxList{ev};
    for k = 1:numel(fi)
        wmode = 'append'; if first, wmode = 'overwrite'; first = false; end
        imwrite(V(:, :, fi(k)), tpath, 'tif', 'WriteMode', wmode, 'Compression', 'lzw');
        nFr = nFr + 1;
    end
    if gapFrames > 0 && ev < numel(frameIdxList)
        for b = 1:gapFrames
            imwrite(blk, tpath, 'tif', 'WriteMode', 'append', 'Compression', 'lzw');
            nFr = nFr + 1;
        end
    end
end
end

function savefig_local(fig, outDir, name, doSave)
% Save a figure as PNG (raster, 200 dpi) into outDir if doSave.
if ~doSave || ~isgraphics(fig), return; end
try
    exportgraphics(fig, fullfile(outDir, [name '.png']), ...
        'Resolution', 200, 'BackgroundColor', 'white');
catch ME
    fprintf(2, 'savefig_local: could not save %s (%s)\n', name, ME.message);
end
end

function img = display_image(raw, clip_pct, g)
% Percentile-clip to [0,1] then gamma -- standard avg-projection display.
clim = prctile(raw(:), clip_pct);
img = (raw - clim(1)) / max(clim(2) - clim(1), eps);
img = min(max(img, 0), 1) .^ g;
end

function rgb = val2rgb(vals, cmap, clim)
% Map values to colormap rows; NaN/out-of-range -> mid gray.
m = size(cmap, 1);
v = (vals(:) - clim(1)) / max(clim(2) - clim(1), eps);
v = min(max(v, 0), 1);
idx = round(v * (m-1)) + 1;
rgb = repmat([0.6 0.6 0.6], numel(vals), 1);
ok  = isfinite(idx);
rgb(ok, :) = cmap(idx(ok), :);
end

function DATA = load_io_data(inputPath, P)
% Load dFF [T x N] + optional spike prob, ROI centroids, avg proj, fps, px.
% Supports a session FOLDER (cpSAM_output.mat) or a single .mat FILE.
DATA = struct('dFF',[],'spikeProb',[],'roiXY',[],'maskL',[],'avgImg',[], ...
    'fps',P.fps_fallback,'px_um',NaN,'name','','folder','');

isFile = (numel(inputPath) > 4) && strcmpi(inputPath(end-3:end), '.mat') ...
    && exist(inputPath, 'file') == 2;

if isFile
    % ---- direct .mat mode ----
    [folder, nm] = fileparts(inputPath);
    DATA.folder = folder;  DATA.name = nm;
    S = load(inputPath);
    assert(isfield(S,'dFF'), 'File has no variable dFF: %s', inputPath);
    DATA.dFF = double(S.dFF);
    if P.nDrop > 0 && size(DATA.dFF,1) > P.nDrop
        DATA.dFF(1:P.nDrop,:) = [];
    end
    if isfield(S,'spikeProb')
        DATA.spikeProb = double(S.spikeProb);
    elseif isfield(S,'spike_train')
        DATA.spikeProb = double(S.spike_train);
    end
    if ~isempty(DATA.spikeProb) && P.nDrop>0 && size(DATA.spikeProb,1)>P.nDrop
        DATA.spikeProb(1:P.nDrop,:) = [];
    end
    if isfield(S,'roiXY')
        DATA.roiXY = double(S.roiXY);
    elseif isfield(S,'maskL')
        DATA.maskL = S.maskL;
    elseif isfield(S,'roiMasks')
        DATA.maskL = S.roiMasks;
    end
    if isempty(DATA.roiXY) && ~isempty(DATA.maskL)
        DATA.roiXY = centroids_from_mask(DATA.maskL, size(DATA.dFF,2));
    end
    if isfield(S,'fps') && S.fps>0
        DATA.fps = S.fps;
    elseif isfield(S,'params') && isfield(S.params,'fps')
        DATA.fps = S.params.fps;
    end
    DATA.px_um = resolve_px(P, folder);
    return
end

% ---- session FOLDER mode (cpSAM pipeline) ----
assert(isfolder(inputPath), 'Not a .mat file or folder: %s', inputPath);
[~, nm] = fileparts(inputPath);
DATA.folder = inputPath;  DATA.name = nm;

sam = dir(fullfile(inputPath, '*cpSAM_output.mat'));
assert(~isempty(sam), 'No *cpSAM_output.mat in %s', inputPath);
SAM = load(fullfile(sam(1).folder, sam(1).name));
assert(isfield(SAM,'F'), 'cpSAM_output has no F: %s', sam(1).name);

% fps + pixel size from metadata
DATA.fps = P.fps_fallback;
sm = struct();
if exist('detect_session_fps','file')==2
    try
        [DATA.fps, sm] = detect_session_fps(inputPath, P.fps_fallback);
    catch
    end
end
DATA.px_um = resolve_px(P, inputPath, sm);

% F -> dFF (cpSAM F already TossFrames-trimmed; nDrop is extra/optional)
F = double(SAM.F);
if P.nDrop > 0 && size(F,1) > P.nDrop, F(1:P.nDrop,:) = []; end
dFFout = helper.dFF_RZ(F, 'FPS', DATA.fps, 'BaselineWinSec', P.BaselineWinSec);
DATA.dFF = double(dFFout.dFF);

% masks -> centroids
if isfield(SAM,'maskL'), DATA.maskL = SAM.maskL;
elseif isfield(SAM,'masks'), DATA.maskL = SAM.masks; end
if ~isempty(DATA.maskL)
    DATA.roiXY = centroids_from_mask(DATA.maskL, size(DATA.dFF,2));
end

% average projection (cellpose input -> aligns perfectly with maskL)
avg_hits = dir(fullfile(inputPath, '*_AVG_for_CP.tif'));
if isempty(avg_hits), avg_hits = dir(fullfile(inputPath, '*AVG*.tif')); end
if ~isempty(avg_hits)
    try
        DATA.avgImg = double(imread(fullfile(avg_hits(1).folder, avg_hits(1).name)));
    catch
    end
end

% optional deconvolved spike probability
sp = dir(fullfile(inputPath, '*ca_spike_data.mat'));
if ~isempty(sp)
    SP = load(fullfile(sp(1).folder, sp(1).name));
    if isfield(SP,'spike_train')
        st = double(SP.spike_train);
    elseif isfield(SP,'spikeProb')
        st = double(SP.spikeProb);
    elseif isfield(SP,'roi_spikes') && isfield(SP.roi_spikes,'spike_train')
        % per-ROI struct array (this repo's ca_spike_data format)
        st = double([SP.roi_spikes.spike_train]);
    else
        st = [];
    end
    if ~isempty(st)
        if P.nDrop>0 && size(st,1)>P.nDrop, st(1:P.nDrop,:) = []; end
        Tc = min(size(st,1), size(DATA.dFF,1));
        DATA.spikeProb = st(1:Tc, :);
        DATA.dFF       = DATA.dFF(1:Tc, :);
        if ~isempty(DATA.roiXY) && size(DATA.spikeProb,2)~=size(DATA.dFF,2)
            DATA.spikeProb = [];   % column mismatch -> ignore spikes
        end
    end
end
end

function px = resolve_px(P, folder, sm)
% pixel size (um/px): override > metadata > zoom-derived > fallback 1.
if isfinite(P.pixelSize_override) && P.pixelSize_override>0
    px = P.pixelSize_override; return
end
px = NaN;
if nargin>=3 && ~isempty(sm)
    if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um>0
        px = sm.pixelSize_um;
    elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor>0
        px = P.PixelSizeBase / sm.zoomFactor;
    end
end
if ~isfinite(px)
    mh = dir(fullfile(folder, '*_meta.mat'));
    if ~isempty(mh)
        m = load(fullfile(mh(1).folder, mh(1).name));
        if isfield(m,'pixelSize_um') && isfinite(m.pixelSize_um) && m.pixelSize_um>0
            px = m.pixelSize_um;
        elseif isfield(m,'zoomFactor') && isfinite(m.zoomFactor) && m.zoomFactor>0
            px = P.PixelSizeBase / m.zoomFactor;
        end
    end
end
if ~isfinite(px), px = 1; warning('Pixel size unknown -> using 1 um/px.'); end
end

function XY = centroids_from_mask(maskL, N)
% Centroids [N x 2] (x,y px) for label values 1..N in a label image.
maskL = double(maskL);
XY = nan(N, 2);
st = regionprops(maskL, 'Centroid');
for i = 1:min(N, numel(st))
    if ~isempty(st(i).Centroid), XY(i,:) = st(i).Centroid; end
end
end

function Y = smooth_cols_gauss(X, w)
% Gaussian smoothing along time (dim 1), per column. w = window in frames.
if w <= 1, Y = X; return; end
g = gausswin(w); g = g/sum(g);
Y = conv2(X, g(:), 'same');
end
