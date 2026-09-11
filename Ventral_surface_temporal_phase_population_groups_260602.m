% Ventral_surface_temporal_phase_population_groups_260602.m
% -----------------------------------------------------------------------
%  Same per-ROI [linear | polar] spike-phase panels as
%  Ventral_surface_temporal_phase_population_260602.m, but for OTHER cell
%  sets, paginated:
%     'IO'            : all inferior-olive ROIs            (272 -> 20/figure)
%     'nonIO_nonsig'  : non-IO ROIs with r < confC         (14  -> 1 figure)
%
%  Each figure: entity 1 = that set's POPULATION mean, then up to perFigCells
%  cells. Each entity = linear phase histogram (left, spikes/cycle %) + polar
%  (right, connecting line + black coherence (r,phase) dot read from PP).
%
%  Source: coherence_polar_260528/coherence_polar_data.mat (PP.colorIdx names
%  the group; PP.r/PP.th = coherence magnitude/preferred phase).
%
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); addpath(scriptDir);

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_260528', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'temporal_phase_batch_260602');

doSets   = {'IO','nonIO_nonsig'};   % which sets to render
perFigCells = 20;                   % cell panels per figure (IO paginates)
perRow   = 3;                       % entities per row (each = time+phase+polar)

nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.015;  nPhaseBins   = 24;
nTimeBins    = 24;     % time-hist bins over one median breath period
% Rmax (ylim/RLim, coherence-r scaling) auto-snaps per set to 20/40/60/80/100

onsetCol = [0.90 0.10 0.10];  peakCol = [0.35 0.75 1.00];
meanCol  = [0.00 0.00 0.00];  cohDotCol = [0 0 0];

doSave   = true;  savePDF = false;   % big paginated figs: PNG always, PDF optional
closeAfter = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== LOAD PP + folder map =====================
S = load(cohData, 'PP','confC','labels','groups');
ioGroup = find(strcmp(S.groups,'IO'));
fprintf('confC=%.4f | %d ROIs total\n', S.confC, numel(S.PP.r));

% one-pass recName -> folder map (IO recordings live under the scanned tree)
allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ===================== RENDER EACH SET =====================
for s = 1:numel(doSets)
    setName = doSets{s};
    switch setName
        case 'IO',           sel = find(S.PP.colorIdx==ioGroup);
        case 'nonIO_nonsig', sel = find(S.PP.colorIdx~=ioGroup & S.PP.r < S.confC);
        otherwise, warning('unknown set %s', setName); continue;
    end
    fprintf('\n=== %s : %d ROIs ===\n', setName, numel(sel));

    % ---- collect per-ROI distributions ----
    M = []; shortLab = {}; spcAll = []; cohR = []; cohTh = []; Fcache = containers.Map();
    dtSpkAll = {}; dtPeakAll = []; periodAll = []; nCycAll = [];
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
            F = Fcache(folderPath);
            R = roi_dist_from_cache(F, roi, ca_lag_sec, edgesP);
            M(end+1,:)      = R.pctC;                         %#ok<SAGROW>
            shortLab{end+1} = sprintf('%s/%s #%d', parts{1}, parts{2}, roi); %#ok<SAGROW>
            spcAll(end+1)   = R.spc;                          %#ok<SAGROW>
            cohR(end+1)     = S.PP.r(sel(k));                 %#ok<SAGROW>
            cohTh(end+1)    = S.PP.th(sel(k));                %#ok<SAGROW>
            dtSpkAll{end+1} = R.dtSpk;                        %#ok<SAGROW>
            dtPeakAll(end+1)= R.dtPeak;                       %#ok<SAGROW>
            periodAll(end+1)= R.medPeriod;                    %#ok<SAGROW>
            nCycAll(end+1)  = R.nCyc;                         %#ok<SAGROW>
        catch ME
            warning('  ERROR %s: %s', lab, ME.message);
        end
    end
    nC = size(M,1);
    if nC==0, warning('no cells for %s', setName); continue; end
    meanC = mean(M,1);
    % common time axis = set median breath period
    timeWin = median(periodAll,'omitnan'); if ~isfinite(timeWin)||timeWin<=0, timeWin = 0.5; end
    edgesT = linspace(0, timeWin, nTimeBins+1); ctrsT = (edgesT(1:end-1)+edgesT(2:end))/2;
    MT = zeros(nC, nTimeBins);
    for i = 1:nC, MT(i,:) = 100*histcounts(dtSpkAll{i}, edgesT)/max(nCycAll(i),eps); end
    meanT = mean(MT,1); tPeakPop = median(dtPeakAll,'omitnan');

    % ---- paginate ----  (Rmax computed PER ENTITY below)
    nPages = ceil(nC/perFigCells); nCol = perRow*3;
    for pg = 1:nPages
        c0 = (pg-1)*perFigCells + 1; c1 = min(pg*perFigCells, nC);
        cellIdx = c0:c1;
        nEnt = numel(cellIdx) + 1;                  % +1 population
        nRow = ceil(nEnt/perRow);
        fig = figure('Color','w','Name',sprintf('%s page %d/%d',setName,pg,nPages), ...
                     'Units','normalized','Position',[0.02 0.04 0.96 0.9]);
        tl = tiledlayout(fig, nRow, nCol, 'TileSpacing','compact','Padding','compact');
        title(tl, sprintf('%s   page %d/%d   (cells %d-%d of %d)   |   spikes/cycle %%   |   RED=onset SKY=peak   |   black dot = coherence (r,\\theta)', ...
              setName, pg, nPages, c0, c1, nC), 'Interpreter','tex','FontWeight','bold');

        % entity 1: population mean (no coh dot)
        [tt, lt, pt] = triple_tiles(1, perRow, nCol);
        rmaxP = pick_rmax(max([meanC(:); meanT(:)]));
        draw_time (tl, tt, ctrsT, meanT, meanCol, onsetCol, peakCol, rmaxP, tPeakPop, sprintf('%s POP (N=%d)', setName, nC));
        draw_phase(tl, lt, ctrsP, meanC, meanCol, onsetCol, peakCol, rmaxP, '');
        draw_polar(tl, pt, ctrsP, meanC, meanCol, 2.0, rmaxP, NaN, NaN, cohDotCol, '');

        % cells on this page (Rmax per cell)
        for e = 1:numel(cellIdx)
            ci = cellIdx(e); [tt, lt, pt] = triple_tiles(e+1, perRow, nCol);
            rmaxE = pick_rmax(max([M(ci,:), MT(ci,:)]));
            ttl = sprintf('%s (%.2f spk/cyc)', shortLab{ci}, spcAll(ci));
            draw_time (tl, tt, ctrsT, MT(ci,:), [0.30 0.30 0.30], onsetCol, peakCol, rmaxE, dtPeakAll(ci), ttl);
            draw_phase(tl, lt, ctrsP, M(ci,:),  [0.30 0.30 0.30], onsetCol, peakCol, rmaxE, '');
            draw_polar(tl, pt, ctrsP, M(ci,:),  [0.20 0.20 0.20], 1.4, rmaxE, cohTh(ci), cohR(ci), cohDotCol, '');
        end

        if doSave
            base = fullfile(outDir, sprintf('%s_page%02d', setName, pg));
            exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
            if savePDF, exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white'); end
            fprintf('  saved %s (cells %d-%d)\n', [base '.png'], c0, c1);
        end
        if closeAfter, close(fig); end
    end
end
fprintf('\nDone. Output: %s\n', outDir);

%% ===================== LOCAL FUNCTIONS =====================
function R = pick_rmax(dataMax)
% snap up to 10/20/40/60/80/100 (or next multiple of 20 if >100).
ladder = [10 20 40 60 80 100];
if ~isfinite(dataMax) || dataMax <= 0, R = 10; return; end
idx = find(ladder >= dataMax, 1);
if isempty(idx), R = ceil(dataMax/20)*20; else, R = ladder(idx); end
end

function F = load_folder(folderPath, nDrop, fallback_fps)
% Load + derive the per-folder pieces shared by all ROIs in a recording.
bp = dir(fullfile(folderPath,'*DLC*breath_peak_data.mat'));
if isempty(bp), bp = dir(fullfile(folderPath,'*breath_peak_data.mat')); end
ip = dir(fullfile(folderPath,'*DLC*breath_insp_start_data.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'*breath_insp_start_data.mat')); end
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
R.pctC = 100 * cntW / max(nCyc,eps);
R.spc  = sum(cntW) / max(nCyc,eps);
R.nCyc = nCyc;
% time from onset (s)
R.dtSpk = time_from_onset(sp,  onF, T) / F.fps;
dtp     = time_from_onset(pkF, onF, T) / F.fps;
R.dtPeak = median(dtp,'omitnan');
if numel(onF) >= 2, R.medPeriod = median(diff(onF))/F.fps; else, R.medPeriod = NaN; end
end

function dt = time_from_onset(frames, onF, T)
dt = nan(numel(frames),1);
if isempty(onF) || isempty(frames), return; end
b = discretize(frames(:), [onF(:); T+1]);
ok = ~isnan(b); dt(ok) = frames(ok) - onF(b(ok));
end

function [tTile, linTile, polTile] = triple_tiles(e, perRow, nCol)
pidx = e - 1; r = floor(pidx/perRow); c3 = mod(pidx, perRow);
tTile = r*nCol + c3*3 + 1; linTile = tTile + 1; polTile = tTile + 2;
end

function draw_time(tl, tile, ctrsT, vals, col, onsetCol, peakCol, Rmax, tPeak, ttl)
ax = nexttile(tl, tile); hold(ax,'on');
bar(ax, ctrsT, vals, 1, 'FaceColor', col, 'FaceAlpha', 0.75, 'EdgeColor','none');
xline(ax, 0, 'Color', onsetCol, 'LineWidth', 1.0);
if isfinite(tPeak), xline(ax, tPeak, 'Color', peakCol, 'LineWidth', 1.0); end
db = (ctrsT(2)-ctrsT(1))/2; xlim(ax,[0 ctrsT(end)+db]); ylim(ax,[0 Rmax]);
xlabel(ax,'t from onset (s)'); ylabel(ax,'spk/cyc %'); set(ax,'FontSize',7);
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
