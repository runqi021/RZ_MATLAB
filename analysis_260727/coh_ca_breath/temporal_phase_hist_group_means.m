% temporal_phase_hist_group_means.m
% -----------------------------------------------------------------------
%  POPULATION-MEAN spike-phase summary: ONE compact figure, one ROW per group,
%  each row = that group's per-cycle-normalized MEAN spike distribution:
%     [ time hist | phase hist (square) | polar ]
%
%  Groups = the two binary factors (IO or not) x (significant or not):
%       sig         : non-IO ROIs, r >= confC
%       nonsig      : non-IO ROIs, r <  confC
%       IO_sig      : IO ROIs,     r >= confC
%       IO_nonsig   : IO ROIs,     r <  confC
%  Empty groups are skipped, so a single-genotype experiment (no IO) shows
%  sig + nonsig only.
%
%  Means-only companion to temporal_phase_population.m (which shows every
%  individual cell, paginated). Same underlying math/drawing.
%   - time hist : spikes (left, spk/cycle %) + breath-PEAK distribution
%                 (right, % of cycles, sky), 0 = inspiration ONSET, 0->1.5 median period
%   - phase hist: spike phase 0..2pi (square)
%   - polar     : phase-distribution connecting line (no coherence dot on a mean)
%   RED = inspiration onset, SKY = breath-peak distribution.
%
%  Source: coh_cfg_260727 -> coherence_polar_data.mat. Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot = fileparts(fileparts(scriptDir)); addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence')); addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
addpath(scriptDir);                                   % so coh_cfg_260727 resolves

%% ===================== USER-EDITABLE =====================
cfg      = coh_cfg_260727();
rootPath = cfg.rootPath;
cohData  = cfg.cohData;
outDir   = fullfile(cfg.outRoot, 'temporal_phase_group_means');

nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.1;    nPhaseBins   = 24;   % lead-comp: spikes 3 frames earlier @30Hz (0.1 s)
nTimeBins    = 24;     % time-hist bins over one median breath period
% Rmax (ylim/RLim) auto-snaps per row to 10/20/40/60/80/100 below

onsetCol = [0.90 0.10 0.10];  peakCol = [0.35 0.75 1.00];
meanCol  = [0.00 0.00 0.00];  cohDotCol = [0 0 0];

doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== LOAD PP + folder map =====================
assert(isfile(cohData), 'coherence_polar_data.mat not found: %s', cohData);
S = load(cohData, 'PP','confC','labels','groups');
ioGroup = find(strcmp(S.groups,'IO')); if isempty(ioGroup), ioGroup = -1; end
fprintf('confC=%.4f | %d ROIs total\n', S.confC, numel(S.PP.r));

% one-pass recName -> folder map
allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ============ FOUR GROUPS: (IO?) x (sig?)  (empty ones auto-skipped) ============
isIO  = (S.PP.colorIdx(:) == ioGroup);
isSig = (S.PP.r(:) >= S.confC);
groups = struct('name', {'sig','nonsig','IO_sig','IO_nonsig'}, ...
                'sel',  { find(~isIO &  isSig), find(~isIO & ~isSig), ...
                          find( isIO &  isSig), find( isIO & ~isSig) });

Fcache = containers.Map();   % shared across groups (a folder can appear in several)
G = struct('name',{},'nC',{},'meanC',{},'meanT',{},'meanP',{},'ctrsT',{});
for g = 1:numel(groups)
    [meanC, meanT, meanP, ctrsT, nC] = group_mean_dist( ...
        groups(g).sel, S, folderMap, Fcache, nDrop, fallback_fps, ca_lag_sec, edgesP, nTimeBins);
    fprintf('  %-14s : %d cells\n', groups(g).name, nC);
    if nC == 0
        fprintf('    (skipped -- no ROIs in this group for this experiment)\n');
        continue;   % e.g. IO group empty on a single-genotype experiment
    end
    gi = numel(G) + 1;
    G(gi).name  = groups(g).name; G(gi).nC = nC;
    G(gi).meanC = meanC; G(gi).meanT = meanT; G(gi).meanP = meanP; G(gi).ctrsT = ctrsT;
end
if isempty(G), error('No groups with ROIs to plot (0 significant, 0 IO, 0 nonIO_nonsig).'); end

%% ===================== ONE FIGURE, N rows x [time | phase | polar] =====================
nCol = 3; nRows = numel(G);
fig = figure('Color','w','Name','population-mean spike distributions (per group)', ...
             'Units','normalized','Position',[0.15 0.05 0.55 0.9]);
tl = tiledlayout(fig, nRows, nCol, 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('population-mean spike distributions (spikes/cycle %%)   |   RED=onset  SKY=peak'), ...
      'Interpreter','tex','FontWeight','bold');

for g = 1:numel(G)
    [tt, lt, pt] = triple_tiles(g, 1, nCol);   % perRow=1 -> one group per row
    rmaxE = pick_rmax(max([G(g).meanC(:); G(g).meanT(:)]));
    ttl = sprintf('%s (N=%d)', G(g).name, G(g).nC);
    draw_time (tl, tt, G(g).ctrsT, G(g).meanT, G(g).meanP, meanCol, onsetCol, peakCol, rmaxE, ttl);
    draw_phase(tl, lt, ctrsP, G(g).meanC, meanCol, onsetCol, peakCol, rmaxE, '');
    draw_polar(tl, pt, ctrsP, G(g).meanC, meanCol, 2.0, rmaxE, NaN, NaN, cohDotCol, '');
end

%% ===================== SAVE =====================
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir, 'pop_group_means_spike_phase_polar');
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'G','ctrsP','nPhaseBins','nTimeBins','ca_lag_sec','nDrop');
    fprintf('Saved %s.png/.pdf/.mat\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function [meanC, meanT, meanP, ctrsT, nC] = group_mean_dist(sel, S, folderMap, Fcache, nDrop, fallback_fps, ca_lag_sec, edgesP, nTimeBins)
% Collect per-ROI distributions for a selection and return the population mean.
M = []; dtSpkAll = {}; dtPeakAll = {}; periodAll = []; nCycAll = [];
for k = 1:numel(sel)
    lab = S.labels{sel(k)}; parts = regexp(lab,'/','split');
    if numel(parts)<4, continue; end
    recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
    if ~isKey(folderMap, recName), warning('no folder: %s', lab); continue; end
    folderPath = folderMap(recName);
    try
        if ~isKey(Fcache, folderPath)
            Fcache(folderPath) = load_folder(folderPath, nDrop, fallback_fps);
        end
        R = roi_dist_from_cache(Fcache(folderPath), roi, ca_lag_sec, edgesP);
        M(end+1,:)       = R.pctC;      %#ok<AGROW>
        dtSpkAll{end+1}  = R.dtSpk;     %#ok<AGROW>
        dtPeakAll{end+1} = R.dtPeak;    %#ok<AGROW>
        periodAll(end+1) = R.medPeriod; %#ok<AGROW>
        nCycAll(end+1)   = R.nCyc;      %#ok<AGROW>
    catch ME
        warning('  ERROR %s: %s', lab, ME.message);
    end
end
nC = size(M,1);
if nC == 0        % empty group (e.g. no IO recordings in a single-genotype experiment)
    meanC = []; meanT = []; meanP = []; ctrsT = [];
    return;
end
meanC = mean(M,1);
timeWin = median(periodAll,'omitnan'); if ~isfinite(timeWin)||timeWin<=0, timeWin = 0.5; end
edgesT = linspace(0, timeWin, nTimeBins+1); ctrsT = (edgesT(1:end-1)+edgesT(2:end))/2;
MT = zeros(nC, nTimeBins); MP = zeros(nC, nTimeBins);
for i = 1:nC
    MT(i,:) = 100*histcounts(dtSpkAll{i}, edgesT)/max(nCycAll(i),eps);
    cP = histcounts(dtPeakAll{i}, edgesT);
    MP(i,:) = 100 * cP / max(sum(cP),1);
end
meanT = mean(MT,1); meanP = mean(MP,1);
end

function R = pick_rmax(dataMax)
% snap up to 10/20/40/60/80/100 (or next multiple of 20 if >100).
ladder = [10 20 40 60 80 100];
if ~isfinite(dataMax) || dataMax <= 0, R = 10; return; end
idx = find(ladder >= dataMax, 1);
if isempty(idx), R = ceil(dataMax/20)*20; else, R = ladder(idx); end
end

function F = load_folder(folderPath, nDrop, fallback_fps)
% Load + derive the per-folder pieces shared by all ROIs in a recording.
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
assert(~isempty(bp),'no breath_peak_data'); assert(~isempty(ip),'no breath_insp_start');
assert(isfile(fullfile(folderPath,'ca_spike_data.mat')),'no ca_spike_data');

F.fps = detect_session_fps(folderPath, fallback_fps);
BP = load(fullfile(bp(1).folder, bp(1).name));
IP = load(fullfile(ip(1).folder, ip(1).name));
CA = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes');
nB = numel(BP.breath);
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;  % PEAK
ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;   % ONSET
ef(1:min(nDrop,numel(ef))) = [];
% Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay events.
if contains(folderPath, fullfile('Vglut2','1124'))
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
F.ev = ev; F.ef = ef; F.roi_spikes = CA.roi_spikes;
F.nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); F.nDff = size(D.dFF,1); end
end

function R = roi_dist_from_cache(F, roi, ca_lag_sec, edgesP)
spk = double(F.roi_spikes(roi).spike_train(:));
T = min([F.nDff, numel(F.ev), numel(F.ef), numel(spk)]);
ev = F.ev(1:T); ef = F.ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);
onF = find(ef>0); pkF = find(ev>0);
phi = piecewise_phase_local(pkF, onF, T);
lag = round(ca_lag_sec*F.fps);
sp = find(spk>0) - lag; sp = sp(sp>=1 & sp<=T);
% phase hist
spk_phi = mod(phi(sp), 2*pi); spk_phi = spk_phi(~isnan(spk_phi));
cntW = histcounts(spk_phi, edgesP);
phv = phi(~isnan(phi)); nCyc = (max(phv)-min(phv))/(2*pi);
occ = histcounts(mod(phv,2*pi), edgesP); Ob = mean(occ(occ>0));      % phase occupancy (dwell)
wOcc = ones(1,numel(occ)); wOcc(occ>0) = Ob ./ occ(occ>0);
R.pctC = (100 * cntW / max(nCyc,eps)) .* wOcc;                        % occupancy-compensated spk/cyc %
R.spc  = sum(cntW) / max(nCyc,eps);
R.nCyc = nCyc;
% ONSET-aligned: spike + breath-PEAK times since the PRECEDING inspiration onset (s), 0..period
R.dtSpk  = time_from_prev(sp,  onF, T) / F.fps;      % 0 = insp onset
R.dtPeak = time_from_prev(pkF, onF, T) / F.fps;      % breath peak after onset (distribution)
if numel(onF) >= 2, R.medPeriod = median(diff(onF))/F.fps; else, R.medPeriod = NaN; end
end

function dt = time_from_prev(frames, refIdx, T)
% time (frames) of each event since its PRECEDING reference event (0..period);
% NaN if before the first reference event.
dt = nan(numel(frames),1);
if isempty(refIdx) || isempty(frames), return; end
b = discretize(frames(:), [refIdx(:); T+1]);
ok = ~isnan(b); dt(ok) = frames(ok) - refIdx(b(ok));
end

function [tTile, linTile, polTile] = triple_tiles(e, perRow, nCol)
pidx = e - 1; r = floor(pidx/perRow); c3 = mod(pidx, perRow);
tTile = r*nCol + c3*3 + 1; linTile = tTile + 1; polTile = tTile + 2;
end

function draw_time(tl, tile, ctrsT, spkVals, pkVals, col, onsetCol, peakCol, Rmax, ttl)
% ONSET-aligned time panel: spikes (left axis, spk/cycle %) + breath-PEAK distribution
% (right axis, % of cycles, sky).  0 = inspiration onset; tiled out to 1.5 median periods.
ax = nexttile(tl, tile); hold(ax,'on');
db = (ctrsT(2)-ctrsT(1))/2; period = ctrsT(end)+db;          % one median breath period
nTile = round(0.5*numel(ctrsT));                             % tile single cycle out to 1.5 period
ctrsX = [ctrsT, ctrsT(1:nTile)+period];
spkX  = [spkVals, spkVals(1:nTile)];
pkX   = [pkVals,  pkVals(1:nTile)];
yyaxis(ax,'left');                                           % spikes
bar(ax, ctrsX, spkX, 1, 'FaceColor', col, 'FaceAlpha', 0.75, 'EdgeColor','none');
set(ax,'YColor','k'); ylim(ax,[0 Rmax]); ylabel(ax,'spk/cyc %');
yyaxis(ax,'right');                                          % breath-peak distribution
bar(ax, ctrsX, pkX, 1, 'FaceColor', peakCol, 'FaceAlpha', 0.45, 'EdgeColor','none');
pkMax = max([pkX 1]); set(ax,'YColor', peakCol*0.7); ylim(ax,[0 pkMax*1.10]); ylabel(ax,'peak %');
xline(ax, 0,      'Color', onsetCol, 'LineWidth', 1.0);      % insp ONSET at 0 (red)
xline(ax, period, 'Color', onsetCol, 'LineWidth', 1.0);      % next onset at 1 period (red)
xlim(ax,[0 1.5*period]);
xlabel(ax,'t from insp onset (s)'); set(ax,'FontSize',7);
box(ax,'on'); pbaspect(ax,[1 1 1]);
title(ax, ttl, 'Interpreter','none', 'FontSize',8);
end

function draw_phase(tl, tile, ctrsP, vals, col, onsetCol, peakCol, Rmax, ttl)
ax = nexttile(tl, tile); hold(ax,'on');
bar(ax, ctrsP, vals, 1, 'FaceColor', col, 'FaceAlpha', 0.75, 'EdgeColor','none');
xline(ax, 0,  'Color', onsetCol, 'LineWidth', 1.0);
xline(ax, pi, 'Color', peakCol,  'LineWidth', 1.0);
xlim(ax,[0 2*pi]); ylim(ax,[0 Rmax]);
set(ax,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'},'FontSize',7);
ylabel(ax,'spk/cyc %'); box(ax,'on'); pbaspect(ax,[1 1 1]);
title(ax, ttl, 'Interpreter','none', 'FontSize',8);
end

function draw_polar(tl, tile, ctrsP, vals, col, lw, Rmax, dotTh, dotR, dotCol, ttl)
pax = polaraxes(tl); pax.Layout.Tile = tile; hold(pax,'on');
polarplot(pax,[ctrsP ctrsP(1)],[vals vals(1)],'-','Color',col,'LineWidth',lw);
if isfinite(dotTh) && isfinite(dotR)
    polarplot(pax, dotTh, dotR*Rmax, 'o', ...
        'MarkerFaceColor',dotCol, 'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.5);
end
pax.RLim=[0 Rmax]; pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise';
pax.RAxisLocation=180; pax.FontSize=7; thetaticks(pax,0:90:270);
title(pax, ttl, 'Interpreter','none', 'FontSize',8);
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];
[events, ord] = sort(events); types = types(ord);
keep = true(size(events));
for i = 2:numel(events), if types(i) == types(i-1), keep(i) = false; end, end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases = nan(size(events)); phi_cur = types(1) * pi;
for i = 1:numel(events), phases(i) = phi_cur; phi_cur = phi_cur + pi; end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
