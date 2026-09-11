% Ventral_surface_temporal_phase_population_svd_260606.m
% -----------------------------------------------------------------------
%  POPULATION spike-phase distribution across all SIGNIFICANT ROIs.
%
%  For every significant ROI (PP.r >= confC from the done coherence-polar
%  analysis) compute the per-cycle-normalized spike distributions
%  (spikes/cycle %): each spike gets the exact piecewise breath phase phi
%  (onset=0, peak=pi).  Each entity is drawn as THREE panels:
%     [ time hist | phase hist (square) | polar ]
%   - time hist  : spikes (left axis) + breath-PEAK distribution (right axis),
%                  0 = inspiration ONSET, 0 -> 1.5 median period (tiled)
%   - phase hist : spike phase 0..2pi (square)
%   - polar      : phase distribution connecting line + black coherence
%                  (r,phase) dot read from PP (no dot on the population panel)
%  Entity 1 = population mean across ROIs; then one entity per ROI.
%  RED = inspiration onset (0 / t=0), SKY = breath-peak distribution.
%
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot=fileparts(scriptDir); addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence')); addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_svd_260606', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'temporal_phase_batch_svd_260606');

nDrop        = 30;        fallback_fps = 30;
ca_lag_sec   = 0.1;       nPhaseBins   = 24;   % lead-comp: spikes 3 frames earlier @30Hz (0.1 s)
nTimeBins    = 24;        % time-hist bins over one median breath period
perRow       = 3;         % entities per row (each entity = 3 panels)
% Rmax (ylim/RLim, and coherence-r scaling) auto-snaps to 20/40/60/80/100 below

onsetCol = [0.90 0.10 0.10];  peakCol = [0.35 0.75 1.00];
meanCol  = [0.00 0.00 0.00];  cohDotCol = [0 0 0];

doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== SIGNIFICANT-ROI LIST =====================
assert(isfile(cohData), 'coherence_polar_data.mat not found: %s', cohData);
S = load(cohData, 'PP', 'confC', 'labels');
sigIdx = find(S.PP.r >= S.confC);
fprintf('confC=%.4f | %d significant ROIs\n', S.confC, numel(sigIdx));

%% ===================== COLLECT PER-ROI DISTRIBUTIONS =====================
M = []; shortLab = {}; spcAll = []; cohR = []; cohTh = [];
dtSpkAll = {}; dtPeakAll = {}; periodAll = []; nCycAll = [];
for k = 1:numel(sigIdx)
    lab = S.labels{sigIdx(k)}; parts = regexp(lab,'/','split');
    if numel(parts)<4, warning('skip malformed: %s', lab); continue; end
    group=parts{1}; recDate=parts{2}; recName=strjoin(parts(3:end-1),'/'); roi=str2double(parts{end});
    folderPath = resolve_folder(rootPath, group, recName);
    if isempty(folderPath), warning('skip (no folder): %s', lab); continue; end
    try
        R = roi_phase_dist(folderPath, roi, nDrop, fallback_fps, ca_lag_sec, edgesP);
        M(end+1,:)    = R.pctC;            %#ok<SAGROW>
        shortLab{end+1} = sprintf('%s/%s #%d', group, recDate, roi); %#ok<SAGROW>
        spcAll(end+1) = R.spc;             %#ok<SAGROW>
        cohR(end+1)   = S.PP.r(sigIdx(k)); %#ok<SAGROW>
        cohTh(end+1)  = S.PP.th(sigIdx(k));%#ok<SAGROW>
        dtSpkAll{end+1} = R.dtSpk;         %#ok<SAGROW>
        dtPeakAll{end+1}= R.dtPeak;        %#ok<SAGROW>
        periodAll(end+1)= R.medPeriod;     %#ok<SAGROW>
        nCycAll(end+1)  = R.nCyc;          %#ok<SAGROW>
        fprintf('  [%2d] %-55s  %.2f spk/cycle\n', k, lab, R.spc);
    catch ME
        warning('  ERROR %s: %s', lab, ME.message);
    end
end
N = size(M,1); assert(N>0, 'No ROI distributions collected.');
meanC = mean(M,1);

% --- common time axis = set median breath period; build time-hist matrix ---
timeWin = median(periodAll,'omitnan'); if ~isfinite(timeWin)||timeWin<=0, timeWin = 0.5; end
% ONSET-aligned time hist: 0 = inspiration ONSET, single-cycle bins (tiled to 1.5 period in draw_time)
edgesT = linspace(0, timeWin, nTimeBins+1); ctrsT = (edgesT(1:end-1)+edgesT(2:end))/2;
MT = zeros(N, nTimeBins);   % spikes since onset      (spk/cycle %)
MP = zeros(N, nTimeBins);   % breath PEAK times since onset (% of cycles -> distribution)
for i = 1:N
    MT(i,:) = 100 * histcounts(dtSpkAll{i}, edgesT) / max(nCycAll(i),eps);
    cP = histcounts(dtPeakAll{i}, edgesT);
    MP(i,:) = 100 * cP / max(sum(cP),1);
end
meanT = mean(MT,1); meanP = mean(MP,1);

%% ===================== N+1 ENTITIES, each = [time | phase | polar] =====================
% Rmax (ylim/RLim + coherence-r scaling) is computed PER ENTITY below.
nPan = N + 1; nCol = perRow*3; nRow = ceil(nPan/perRow);
fig = figure('Color','w','Name','population + per-ROI spike distributions', ...
             'Units','normalized','Position',[0.02 0.04 0.96 0.9]);
tl = tiledlayout(fig, nRow, nCol, 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('spike distributions (spikes/cycle %%)   |   N=%d sig ROIs   |   RED=onset  SKY=peak   |   black dot = coherence (r,\\theta)', N), ...
      'Interpreter','tex','FontWeight','bold');

% entity 1: population mean (no coherence dot)
[tt, lt, pt] = triple_tiles(1, perRow, nCol);
rmaxP = pick_rmax(max([meanC(:); meanT(:)]));
draw_time  (tl, tt, ctrsT, meanT, meanP, meanCol, onsetCol, peakCol, rmaxP, sprintf('POP mean (N=%d)', N));
draw_phase (tl, lt, ctrsP, meanC, meanCol, onsetCol, peakCol, rmaxP, '');
draw_polar (tl, pt, ctrsP, meanC, meanCol, 2.0, rmaxP, NaN, NaN, cohDotCol, '');

% entities 2..N+1: each ROI (Rmax per ROI)
for i = 1:N
    [tt, lt, pt] = triple_tiles(i+1, perRow, nCol);
    rmaxE = pick_rmax(max([M(i,:), MT(i,:)]));
    ttl = sprintf('%s (%.2f spk/cyc)', shortLab{i}, spcAll(i));
    draw_time  (tl, tt, ctrsT, MT(i,:), MP(i,:), [0.30 0.30 0.30], onsetCol, peakCol, rmaxE, ttl);
    draw_phase (tl, lt, ctrsP, M(i,:),  [0.30 0.30 0.30], onsetCol, peakCol, rmaxE, '');
    draw_polar (tl, pt, ctrsP, M(i,:),  [0.20 0.20 0.20], 1.4, rmaxE, cohTh(i), cohR(i), cohDotCol, '');
end

%% ===================== SAVE =====================
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir, 'population_perROI_spike_phase_polar');
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'M','MT','MP','ctrsP','ctrsT','meanC','meanT','meanP','shortLab','timeWin','nPhaseBins','ca_lag_sec','nDrop');
    fprintf('Saved %s.png/.pdf/.mat\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function R = pick_rmax(dataMax)
% snap up to 10/20/40/60/80/100 (or next multiple of 20 if >100).
ladder = [10 20 40 60 80 100];
if ~isfinite(dataMax) || dataMax <= 0, R = 10; return; end
idx = find(ladder >= dataMax, 1);
if isempty(idx), R = ceil(dataMax/20)*20; else, R = ladder(idx); end
end

function [tTile, linTile, polTile] = triple_tiles(e, perRow, nCol)
pidx = e-1; r = floor(pidx/perRow); c3 = mod(pidx, perRow);
tTile = r*nCol + c3*3 + 1; linTile = tTile+1; polTile = tTile+2;
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
    polarplot(pax, dotTh, dotR*Rmax, 'o', 'MarkerFaceColor',dotCol, 'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.5);
end
pax.RLim=[0 Rmax]; pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise';
pax.RAxisLocation=180; pax.FontSize=7; thetaticks(pax,0:90:270);
title(pax, ttl, 'Interpreter','none', 'FontSize',8);
end

function R = roi_phase_dist(folderPath, roi, nDrop, fallback_fps, ca_lag_sec, edgesP)
% Per-cycle-normalized phase distribution + spike/peak times-from-onset (s).
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
if isempty(bp), bp = dir(fullfile(folderPath,'breath_peak_pc1.mat')); end
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat')); end
assert(~isempty(bp),'no breath_peak_data'); assert(~isempty(ip),'no breath_insp_start');
assert(isfile(fullfile(folderPath,'ca_spike_data.mat')),'no ca_spike_data');

fps = detect_session_fps(folderPath, fallback_fps);
BP  = load(fullfile(bp(1).folder, bp(1).name));
IP  = load(fullfile(ip(1).folder, ip(1).name));
CA  = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes');
nB  = numel(BP.breath);

ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;  % PEAK
ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;   % ONSET
ef(1:min(nDrop,numel(ef))) = [];
% Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay events.
if contains(folderPath, fullfile('Vglut2','1124'))
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
spk = double(CA.roi_spikes(roi).spike_train(:));

nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); nDff = size(D.dFF,1); end
T = min([nDff, numel(ev), numel(ef), numel(spk)]);
ev = ev(1:T); ef = ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);

onF = find(ef>0); pkF = find(ev>0);
phi = piecewise_phase_local(pkF, onF, T);   % onset=0, peak=pi
lag = round(ca_lag_sec*fps);
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
R.dtSpk  = time_from_prev(sp,  onF, T) / fps;      % 0 = insp onset
R.dtPeak = time_from_prev(pkF, onF, T) / fps;      % breath peak after onset (distribution)
if numel(onF) >= 2, R.medPeriod = median(diff(onF))/fps; else, R.medPeriod = NaN; end
end

function dt = time_from_prev(frames, refIdx, T)
% time (frames) of each event since its PRECEDING reference event (0..period);
% NaN if before the first reference event.
dt = nan(numel(frames),1);
if isempty(refIdx) || isempty(frames), return; end
b = discretize(frames(:), [refIdx(:); T+1]);
ok = ~isnan(b); dt(ok) = frames(ok) - refIdx(b(ok));
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

function folderPath = resolve_folder(rootPath, group, recName)
folderPath = '';
d = dir(fullfile(rootPath, group, '**', recName)); d = d([d.isdir]);
if isempty(d), d = dir(fullfile(rootPath, '**', recName)); d = d([d.isdir]); end
if ~isempty(d), folderPath = fullfile(d(1).folder, d(1).name); end
end
