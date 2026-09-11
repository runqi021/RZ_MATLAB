% ibi_histogram_260814.m
% -----------------------------------------------------------------------
%  Distribution of inter-breath interval across the archive, at three levels:
%
%    every CYCLE      one entry per inspiration-to-inspiration interval
%    every RECORDING  one entry per recording (its median IBI)
%    every CELL       one entry per cell (median of all its cycles, pooled over
%                     its recordings -- the number printed on the per-cell figures
%                     and the one a fast/slow split would use)
%
%  Drawn to answer whether an IBI = 1 s split is cutting a gap or cutting through
%  a mode. A cut through a dense region assigns near-identical cells to opposite
%  bins; a cut in a trough is defensible.
%
%  Runqi Zhang / 2026-08-14
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath,'pop_analysis_260813');
L   = load(fullfile(outDir,'pop_features.mat'),'T');  T = L.T;
Din = load(fullfile(rootPath,'event_latency_260811','event_latency_data.mat'), ...
           'CELL','OBS','REC');
% T is post-merge, event_latency_data.mat is pre-merge -- pool through this.
% Here it decides which RECORDINGS are visited at all: a recording reachable only
% through the merged-away partner would otherwise never be counted.
obsOf = pooled_obs_260814(Din.CELL, Din.OBS);

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
nDrop = 30;  cut = 1.0;

%% ---- every recording, once ----
recSeen = containers.Map('KeyType','char','ValueType','logical');
allCyc = [];  recIBI = [];  recSess = {};
for ci = 1:height(T)
    c = T.cell(ci);
    if isempty(obsOf{c}), continue; end
    for o = obsOf{c}(:)'
        p  = regexp(Din.OBS(o).label,'/','split');
        fp = Din.REC(Din.OBS(o).rec).folder;
        if ~isfolder(fp) || isKey(recSeen, fp), continue; end
        recSeen(fp) = true;
        [fps,~] = detect_session_fps(fp,30);
        BP = load(fullfile(fp,'breath_peak_pc1.mat'));
        IP = load(fullfile(fp,'breath_insp_start_pc1.mat'));
        nB = numel(BP.breath);
        ef = zeros(nB,1); fi = round(IP.insp_start_idx(:));
        ef(fi(fi>=1 & fi<=nB)) = 1; ef(1:min(nDrop,numel(ef))) = [];
        f = find(ef>0);
        if numel(f) < 3, continue; end
        d = diff(f)/fps;
        allCyc  = [allCyc; d];                       %#ok<AGROW>
        recIBI(end+1) = median(d);                   %#ok<SAGROW>
        recSess{end+1} = [p{1} '_' p{2}];            %#ok<SAGROW>
    end
end
cellIBI = T.IBI;

fprintf('cycles %d | recordings %d | cells %d\n', numel(allCyc), numel(recIBI), numel(cellIBI));
fprintf('cycle  IBI: median %.2f  IQR %.2f-%.2f  range %.2f-%.2f s\n', ...
        median(allCyc), prctile(allCyc,25), prctile(allCyc,75), min(allCyc), max(allCyc));
fprintf('rec    IBI: median %.2f  range %.2f-%.2f s\n', median(recIBI), min(recIBI), max(recIBI));
fprintf('cell   IBI: median %.2f  range %.2f-%.2f s\n', median(cellIBI), min(cellIBI), max(cellIBI));
fprintf('at a %.1f s cut:  cells fast %d, slow %d | recordings fast %d, slow %d\n', ...
        cut, nnz(cellIBI<=cut), nnz(cellIBI>cut), nnz(recIBI<=cut), nnz(recIBI>cut));

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','normalized','Position',[0.06 0.10 0.84 0.72]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf(['inter-breath interval across the archive   |   %d cycles, ' ...
      '%d recordings, %d cells   |   red = the 1 s split'], ...
      numel(allCyc), numel(recIBI), numel(cellIBI)), 'FontWeight','bold');

ax = nexttile(tl); hold(ax,'on');
histogram(ax, allCyc, 0:0.05:6, 'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
xline(ax, cut,'r-','LineWidth',1.6);
xlabel(ax,'IBI (s)'); ylabel(ax,'cycles'); xlim(ax,[0 6]);
title(ax, sprintf('every CYCLE (n=%d)', numel(allCyc)));
box(ax,'on'); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
histogram(ax, allCyc, 0:0.05:6, 'FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
xline(ax, cut,'r-','LineWidth',1.6);
set(ax,'YScale','log'); xlabel(ax,'IBI (s)'); ylabel(ax,'cycles (log)'); xlim(ax,[0 6]);
title(ax,'every CYCLE, log count -- shows the tail');
box(ax,'on'); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
histogram(ax, recIBI, 0:0.1:3.2, 'FaceColor',[0.20 0.45 0.75],'EdgeColor','none');
xline(ax, cut,'r-','LineWidth',1.6);
xlabel(ax,'median IBI of the recording (s)'); ylabel(ax,'recordings'); xlim(ax,[0 3.2]);
title(ax, sprintf('every RECORDING (n=%d)', numel(recIBI)));
box(ax,'on'); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
histogram(ax, cellIBI, 0:0.1:3.2, 'FaceColor',[0.85 0.45 0.10],'EdgeColor','none');
xline(ax, cut,'r-','LineWidth',1.6);
xlabel(ax,'pooled IBI of the cell (s)'); ylabel(ax,'cells'); xlim(ax,[0 3.2]);
title(ax, sprintf('every CELL (n=%d)  fast %d / slow %d', numel(cellIBI), ...
      nnz(cellIBI<=cut), nnz(cellIBI>cut)));
box(ax,'on'); grid(ax,'on');

exportgraphics(fig, fullfile(outDir,'ibi_histogram.png'),'Resolution',200,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'ibi_histogram.pdf'),'ContentType','vector','BackgroundColor','white');
fprintf('\nsaved ibi_histogram.png / .pdf\n');
