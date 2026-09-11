% Ventral_surface_phase_popstate_traj_svd_260606.m
% -----------------------------------------------------------------------
%  POPULATION TRAJECTORY (state-space) over the breath cycle, ONE figure
%  per set:  'sig' | 'nonIO_nonsig' | 'IO'.
%
%  Build A = [nBins x N neurons] (each column = a neuron's per-cycle tuning
%  curve), z-score each neuron, embed the BINS (nBins points) two ways:
%     PCA    : bins as observations -> score = [nBins x PCs]
%     Isomap : kNN graph on bin-points -> geodesics -> classical MDS
%  Each plotted point is one bin (a population state); connecting them in
%  order draws the trajectory/loop over the cycle.  nBins (=24) dots.
%
%  PHASE trajectory : bins = breath phase (onset=0, peak=pi)
%  TIME  trajectory : bins = within-cycle time fraction (spike-onset)/cycle
%
%  Figure (2x3):
%     row 1  phase PCA | time PCA | PCA scree (cumulative var)
%     row 2  phase Isomap | time Isomap | Isomap residual variance
%
%  Source: coherence_polar_svd_260606/coherence_polar_data.mat
%  Dependencies: detect_session_fps.m, Statistics TB (pca, cmdscale, graph)
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot=fileparts(scriptDir); addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence')); addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_svd_260606', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'temporal_phase_batch_svd_260606');

doSets   = {'sig','nonIO_nonsig','IO'};
nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.1; nPhaseBins   = 24;   % lead-comp: spikes 3 frames earlier @30Hz (0.1 s)
nTimeBins    = 24;     % within-cycle time bins (fraction of cycle, 0..1)
minSpikes    = 1;
kIso         = 4;      % Isomap neighbours among the bin-points (auto-raised)
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;
edgesT = linspace(0, 1, nTimeBins+1);     ctrsT = (edgesT(1:end-1)+edgesT(2:end))/2;

%% ===================== LOAD PP + folder map =====================
S = load(cohData, 'PP','confC','labels','groups');
ioGroup = find(strcmp(S.groups,'IO'));
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
        case 'sig',          sel = find(S.PP.r >= S.confC);
        case 'nonIO_nonsig', sel = find(S.PP.colorIdx~=ioGroup & S.PP.r < S.confC);
        case 'IO',           sel = find(S.PP.colorIdx==ioGroup);
        otherwise, warning('unknown set %s', setName); continue;
    end
    fprintf('\n=== %s : %d ROIs ===\n', setName, numel(sel));

    % ---- collect per-neuron phase + time tuning curves ----
    Mp = []; Mt = []; Fcache = containers.Map();
    for k = 1:numel(sel)
        lab = S.labels{sel(k)}; parts = regexp(lab,'/','split');
        if numel(parts)<4, continue; end
        recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
        if ~isKey(folderMap, recName), continue; end
        folderPath = folderMap(recName);
        try
            if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_folder(folderPath, nDrop, fallback_fps); end
            [pctP, pctT, nspk] = roi_hist(Fcache(folderPath), roi, ca_lag_sec, edgesP, edgesT);
            if nspk < minSpikes, continue; end
            Mp(end+1,:) = pctP;  Mt(end+1,:) = pctT;             %#ok<SAGROW>
        catch ME
            warning('  %s: %s', lab, ME.message);
        end
    end
    N = size(Mp,1);
    fprintf('  %d neurons; trajectory = %d bin-points\n', N, nPhaseBins);
    if N < 3, warning('  too few neurons (%d); skipping', N); continue; end

    % ---- PCA (bins = observations) ----
    [Yp, ep] = pca_traj(Mp);   [Yt, et] = pca_traj(Mt);
    % ---- Isomap (bin-points) ----
    [YpI, kpI, rpI, keepP] = isomap_traj(Mp, kIso, 3);
    [YtI, ktI, rtI, keepT] = isomap_traj(Mt, kIso, 3);
    fprintf('  PCA top3: phase %.1f%% time %.1f%% | Isomap k: phase %d time %d\n', ...
            sum(ep(1:3)), sum(et(1:3)), kpI, ktI);

    %% ---- figure 2x3 ----
    fig = figure('Color','w','Name',sprintf('%s trajectory PCA/Isomap',setName), ...
                 'Units','centimeters','Position',[2 2 36 22]);
    tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('%s  |  population breath-cycle trajectory (state space)  |  N=%d neurons, %d bins', ...
          setName, N, nPhaseBins), 'Interpreter','none','FontWeight','bold');

    pcLbl = @(e) {sprintf('PC1 (%.0f%%)',e(1)), sprintf('PC2 (%.0f%%)',e(2)), sprintf('PC3 (%.0f%%)',e(3))};
    isoLbl = {'Isomap 1','Isomap 2','Isomap 3'};

    % row 1: PCA
    draw_loop(nexttile(tl,1), Yp, pcLbl(ep), ctrsP, hsv(256), [0 2*pi], {'0','\pi','2\pi'}, [0 pi 2*pi], 'phase', 'PHASE  PCA');
    draw_loop(nexttile(tl,2), Yt, pcLbl(et), ctrsT, parula(256), [0 1], {'0','0.5','1'}, [0 .5 1], 'cycle frac', 'TIME  PCA');
    ax3 = nexttile(tl,3); hold(ax3,'on');
    np = min(10,numel(ep));
    plot(ax3,1:np,cumsum(ep(1:np)),'o-','Color',[0.2 0.4 0.8],'LineWidth',1.2,'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(ax3,1:np,cumsum(et(1:np)),'s-','Color',[0.85 0.45 0.1],'LineWidth',1.2,'MarkerFaceColor',[0.85 0.45 0.1]);
    xlabel(ax3,'PC #'); ylabel(ax3,'cumulative % var'); ylim(ax3,[0 100]); grid(ax3,'on'); box(ax3,'on');
    legend(ax3,{sprintf('phase (top3=%.0f%%)',sum(ep(1:3))),sprintf('time (top3=%.0f%%)',sum(et(1:3)))},'Location','southeast','Box','off','FontSize',8);
    title(ax3,'PCA scree');

    % row 2: Isomap
    draw_loop(nexttile(tl,4), YpI, isoLbl, ctrsP(keepP), hsv(256), [0 2*pi], {'0','\pi','2\pi'}, [0 pi 2*pi], 'phase', sprintf('PHASE  Isomap (k=%d)',kpI));
    draw_loop(nexttile(tl,5), YtI, isoLbl, ctrsT(keepT), parula(256), [0 1], {'0','0.5','1'}, [0 .5 1], 'cycle frac', sprintf('TIME  Isomap (k=%d)',ktI));
    ax6 = nexttile(tl,6); hold(ax6,'on');
    plot(ax6,1:numel(rpI),rpI,'o-','Color',[0.2 0.4 0.8],'LineWidth',1.2,'MarkerFaceColor',[0.2 0.4 0.8]);
    plot(ax6,1:numel(rtI),rtI,'s-','Color',[0.85 0.45 0.1],'LineWidth',1.2,'MarkerFaceColor',[0.85 0.45 0.1]);
    xlabel(ax6,'Isomap dim'); ylabel(ax6,'residual variance'); ylim(ax6,[0 1]); grid(ax6,'on'); box(ax6,'on');
    legend(ax6,{'phase','time'},'Location','northeast','Box','off','FontSize',8);
    title(ax6,'Isomap residual variance');

    if doSave
        base = fullfile(outDir, sprintf('phase_time_trajectory_%s', setName));
        exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
        exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        save([base '.mat'], 'Yp','Yt','ep','et','YpI','YtI','rpI','rtI','ctrsP','ctrsT','N','setName');
        fprintf('  saved %s.png/.pdf/.mat\n', base);
    end
end
fprintf('\nDone. Output: %s\n', outDir);

%% ===================== LOCAL FUNCTIONS =====================
function [Y, explained] = pca_traj(M)
A = M';  sd = std(A,0,1); sd(sd==0)=1; A = (A - mean(A,1)) ./ sd;   % [bins x neurons], z-score per neuron
[~, score, ~, ~, explained] = pca(A);
nd = min(3, size(score,2)); Y = score(:,1:nd); if nd<3, Y(:,nd+1:3)=0; end
explained(end+1:3) = 0;
end

function [Y, kUsed, resid, keep] = isomap_traj(M, k0, nDims)
A = M';  sd = std(A,0,1); sd(sd==0)=1; A = (A - mean(A,1)) ./ sd;   % [bins x neurons]
n = size(A,1); D = squareform(pdist(A));
kUsed = min(k0, n-1);
while true
    if max(conncomp(graph(knn_adj(D,kUsed))))==1 || kUsed>=n-1, break; end
    kUsed = kUsed + 1;
end
G = graph(knn_adj(D,kUsed)); comp = conncomp(G);
keep = find(comp == mode(comp));
Dg = distances(subgraph(G, keep));
[Yf,~] = cmdscale(Dg); nd = min(nDims, size(Yf,2)); Y = Yf(:,1:nd); if nd<3, Y(:,nd+1:3)=0; end
mask = tril(true(numel(keep)),-1); dgv = Dg(mask); resid = zeros(1,nd);
for d = 1:nd, Dl = squareform(pdist(Yf(:,1:d))); resid(d) = 1 - corr(dgv, Dl(mask))^2; end
end

function A = knn_adj(D, k)
n = size(D,1); A = zeros(n); [~,ord] = sort(D,2);
for i = 1:n, nb = ord(i,2:min(k+1,n)); A(i,nb) = D(i,nb); end
A = max(A, A');
end

function draw_loop(ax, Y, axLbls, binCtrs, cmap, caxlim, cbLbls, cbTicks, cbName, ttl)
hold(ax,'on');
plot3(ax, [Y(:,1);Y(1,1)], [Y(:,2);Y(1,2)], [Y(:,3);Y(1,3)], '-', 'Color',[0.6 0.6 0.6],'LineWidth',1.0);
ncm = size(cmap,1); cidx = max(1, ceil((binCtrs(:)-caxlim(1))/(caxlim(2)-caxlim(1))*ncm));
scatter3(ax, Y(:,1), Y(:,2), Y(:,3), 55, cmap(min(cidx,ncm),:), 'filled','MarkerEdgeColor','k','LineWidth',0.4);
grid(ax,'on'); box(ax,'on'); view(ax,30,25); axis(ax,'tight');
xlabel(ax,axLbls{1}); ylabel(ax,axLbls{2}); zlabel(ax,axLbls{3});
colormap(ax, cmap); caxis(ax, caxlim);
cb = colorbar(ax); cb.Label.String = cbName; cb.Ticks = cbTicks; cb.TickLabels = cbLbls;
title(ax, ttl);
end

function F = load_folder(folderPath, nDrop, fallback_fps)
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
if isempty(bp), bp = dir(fullfile(folderPath,'breath_peak_pc1.mat')); end
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat')); end
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
% Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay events.
if contains(folderPath, fullfile('Vglut2','1124'))
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
F.ev = ev; F.ef = ef; F.roi_spikes = CA.roi_spikes;
F.nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); F.nDff = size(D.dFF,1); end
end

function [pctP, pctT, nspk] = roi_hist(F, roi, ca_lag_sec, edgesP, edgesT)
spk = double(F.roi_spikes(roi).spike_train(:));
T = min([F.nDff, numel(F.ev), numel(F.ef), numel(spk)]);
ev = F.ev(1:T); ef = F.ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);
onF = find(ef>0);
phi = piecewise_phase_local(find(ev>0), onF, T);
lag = round(ca_lag_sec*F.fps);
sp  = find(spk>0) - lag; sp = sp(sp>=1 & sp<=T);
phv = phi(~isnan(phi)); nCyc = (max(phv)-min(phv))/(2*pi);

spk_phi = mod(phi(sp), 2*pi); spk_phi = spk_phi(~isnan(spk_phi));
occ = histcounts(mod(phv,2*pi), edgesP); Ob = mean(occ(occ>0));      % phase occupancy (dwell)
wOcc = ones(1,numel(occ)); wOcc(occ>0) = Ob ./ occ(occ>0);
pctP = (100 * histcounts(spk_phi, edgesP) / max(nCyc,eps)) .* wOcc;  % occupancy-compensated spk/cyc %
nspk = numel(spk_phi);

pctT = zeros(1, numel(edgesT)-1);
if numel(onF) >= 2
    b = discretize(sp, onF); ok = ~isnan(b);
    dur  = onF(b(ok)+1) - onF(b(ok));
    frac = (sp(ok) - onF(b(ok))) ./ dur;
    pctT = 100 * histcounts(frac, edgesT) / max(nCyc,eps);
end
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
