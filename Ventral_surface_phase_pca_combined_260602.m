% Ventral_surface_phase_pca_combined_260602.m
% -----------------------------------------------------------------------
%  ONE PCA manifold of breath-phase tuning for the WHOLE population.
%  PCA is fit on ALL neurons together ([N x 24] z-scored phase histograms);
%  each neuron is plotted at its PC1/PC2/PC3 score, colored by preferred
%  phase (circular mean of its OWN spike-phase histogram, onset=0/peak=pi).
%
%  Marker shape by category:
%     significant (r >= confC)   -> circle
%     non-IO non-significant     -> square
%     IO                         -> square (smaller, faint)
%
%  Source: coherence_polar_260528/coherence_polar_data.mat
%  Dependencies: detect_session_fps.m, Statistics TB (pca)
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); addpath(scriptDir);

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_260528', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'temporal_phase_batch_260602');

nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.015;  nPhaseBins   = 24;
minSpikes    = 1;
normMode     = 'zscore';  % 'zscore' | 'sum' | 'peak' | 'raw'
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== LOAD PP + folder map =====================
S = load(cohData, 'PP','confC','labels','groups');
ioGroup = find(strcmp(S.groups,'IO'));
allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ===================== COLLECT ALL NEURONS =====================
% category: 1 = sig, 2 = nonIO non-sig, 3 = IO
isIO  = (S.PP.colorIdx==ioGroup);
isSig = (S.PP.r >= S.confC);
cat_all_PP = zeros(numel(S.PP.r),1);
cat_all_PP(isSig)              = 1;
cat_all_PP(~isIO & ~isSig)     = 2;
cat_all_PP(isIO)               = 3;   % (IO has 0 sig, so no overlap)

M = []; mu = []; cat = []; Fcache = containers.Map();
for k = 1:numel(S.PP.r)
    lab = S.labels{k}; parts = regexp(lab,'/','split');
    if numel(parts)<4, continue; end
    recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
    if ~isKey(folderMap, recName), continue; end
    folderPath = folderMap(recName);
    try
        if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_folder(folderPath, nDrop, fallback_fps); end
        [pctC, nspk] = roi_phase_hist(Fcache(folderPath), roi, ca_lag_sec, edgesP);
        if nspk < minSpikes, continue; end
        M(end+1,:) = pctC;                                       %#ok<SAGROW>
        mu(end+1)  = mod(angle(sum(pctC(:).*exp(1i*ctrsP(:)))),2*pi); %#ok<SAGROW>
        cat(end+1) = cat_all_PP(k);                              %#ok<SAGROW>
    catch ME
        warning('  %s: %s', lab, ME.message);
    end
end
n = size(M,1);
fprintf('Collected %d neurons (sig=%d, nonIO-nonsig=%d, IO=%d)\n', ...
        n, sum(cat==1), sum(cat==2), sum(cat==3));

%% ===================== PCA (all neurons together) =====================
X = M;
switch normMode
    case 'zscore', sd = std(X,0,2); sd(sd==0)=1; X = (X - mean(X,2)) ./ sd;
    case 'sum',    rs = sum(X,2);   rs(rs==0)=1; X = X ./ rs;
    case 'peak',   pk = max(X,[],2);pk(pk==0)=1; X = X ./ pk;
end
ws = warning('off','stats:pca:ColRankDefX');
[~, score, ~, ~, explained] = pca(X);
warning(ws);
Y = score(:,1:3); e3 = explained; e3(end+1:3)=0;
fprintf('PCA top3 = %.1f%% variance\n', sum(explained(1:min(3,end))));

%% ===================== COMBINED FIGURE =====================
fig = figure('Color','w','Name','population phase PCA (all neurons)', ...
             'Units','centimeters','Position',[3 3 18 15]);
ax = axes(fig); hold(ax,'on');
cmap = hsv(256); cidx = max(1, ceil(mu(:)/(2*pi)*256));

markers = {'o','s','s'};            % sig=circle, nonsig=square, IO=square
szs     = [60 60 22];
alphas  = [1.0 1.0 0.45];
edges   = {'k','k','none'};
labs    = {sprintf('sig (n=%d)',sum(cat==1)), ...
           sprintf('nonIO non-sig (n=%d)',sum(cat==2)), ...
           sprintf('IO (n=%d)',sum(cat==3))};
for c = 1:3
    idx = find(cat==c); if isempty(idx), continue; end
    scatter3(ax, Y(idx,1), Y(idx,2), Y(idx,3), szs(c), cmap(cidx(idx),:), 'filled', ...
        markers{c}, 'MarkerEdgeColor',edges{c}, 'LineWidth',0.3, 'MarkerFaceAlpha',alphas(c));
end
grid(ax,'on'); box(ax,'on'); view(ax,30,25); axis(ax,'tight');
xlabel(ax,sprintf('PC1 (%.1f%%)',e3(1))); ylabel(ax,sprintf('PC2 (%.1f%%)',e3(2))); zlabel(ax,sprintf('PC3 (%.1f%%)',e3(3)));
colormap(ax, hsv(256)); caxis(ax,[0 2*pi]);
cb = colorbar(ax); cb.Label.String='preferred phase (rad)'; cb.Ticks=[0 pi 2*pi]; cb.TickLabels={'0','\pi','2\pi'};
title(ax, sprintf('population phase-tuning PCA  |  N=%d neurons  |  norm=%s', n, normMode), 'Interpreter','none');

% shape legend (gray dummies; color = preferred phase)
hL = gobjects(3,1);
for c = 1:3
    hL(c) = scatter3(ax, nan, nan, nan, szs(c), [0.5 0.5 0.5], 'filled', ...
        markers{c}, 'MarkerEdgeColor',edges{c});
end
legend(ax, hL, labs, 'Location','northeast','Box','off','FontSize',8);

if doSave
    base = fullfile(outDir, 'phase_pca_combined');
    exportgraphics(fig, [base '.png'], 'Resolution',220, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'Y','mu','cat','explained','normMode','ctrsP');
    fprintf('Saved %s.png/.pdf/.mat\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function F = load_folder(folderPath, nDrop, fallback_fps)
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
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;
ef(1:min(nDrop,numel(ef))) = [];
F.ev = ev; F.ef = ef; F.roi_spikes = CA.roi_spikes;
F.nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); F.nDff = size(D.dFF,1); end
end

function [pctC, nspk] = roi_phase_hist(F, roi, ca_lag_sec, edgesP)
spk = double(F.roi_spikes(roi).spike_train(:));
T = min([F.nDff, numel(F.ev), numel(F.ef), numel(spk)]);
ev = F.ev(1:T); ef = F.ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);
phi = piecewise_phase_local(find(ev>0), find(ef>0), T);
lag = round(ca_lag_sec*F.fps);
sp = find(spk>0) - lag; sp = sp(sp>=1 & sp<=T);
spk_phi = mod(phi(sp), 2*pi); spk_phi = spk_phi(~isnan(spk_phi));
cntW = histcounts(spk_phi, edgesP);
phv = phi(~isnan(phi)); nCyc = (max(phv)-min(phv))/(2*pi);
pctC = 100 * cntW / max(nCyc,eps);
nspk = numel(spk_phi);
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
