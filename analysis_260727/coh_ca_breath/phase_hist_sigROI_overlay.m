% phase_hist_sigROI_overlay.m
% -----------------------------------------------------------------------
%  LINEAR (non-polar) overlay of the temporal phase histogram for every
%  SIGNIFICANT ROI on ONE axes.
%
%  Same per-ROI quantity as Ventral_surface_phase_polar_allROI_svd_260606:
%  occupancy-compensated spike-phase histogram (spikes/cycle %),
%  phase = onset 0, peak pi (piecewise), 24 bins, spikes lead-shifted
%  ca_lag=0.1 s, Vglut2/1124 +1-frame fix, nDrop=30.  Only difference:
%  plotted on a linear 0..2pi axis, restricted to sig ROIs (PP.r >= confC),
%  one line per ROI colored by group.
%
%  Source: coherence_polar_svd_260606/coherence_polar_data.mat
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(scriptDir);                                   % so coh_cfg_260727 resolves

%% ===================== USER-EDITABLE =====================
cfg      = coh_cfg_260727();
rootPath = cfg.rootPath;
cohData  = cfg.cohData;
outDir   = fullfile(cfg.outRoot, 'phase_hist_sigROI');

nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.1;    nPhaseBins   = 24;     % spikes 3 frames earlier @30Hz
minSpikes    = 1;
lineW        = 1.0;    % per-ROI line width
showMean     = true;   % overlay thick population mean per group
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== LOAD PP + folder map =====================
S = load(cohData, 'PP','confC','labels','groups','group_colors');
sigIdx = find(S.PP.r >= S.confC);
fprintf('%d significant ROIs (PP.r >= confC=%.4f)\n', numel(sigIdx), S.confC);

allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ===================== COLLECT SIG ROIs =====================
H = []; O = []; gi = []; labs = {}; Fcache = containers.Map();
for kk = 1:numel(sigIdx)
    k = sigIdx(kk);
    lab = S.labels{k}; parts = regexp(lab,'/','split');
    if numel(parts)<4, continue; end
    recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
    if ~isKey(folderMap, recName), warning('no folder: %s', lab); continue; end
    folderPath = folderMap(recName);
    try
        if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_folder(folderPath, nDrop, fallback_fps); end
        [pctP, nspk, occFrac] = roi_phase_hist(Fcache(folderPath), roi, ca_lag_sec, edgesP);
        if nspk < minSpikes, continue; end
        H(end+1,:) = pctP;             %#ok<SAGROW>
        O(end+1,:) = occFrac;          %#ok<SAGROW>
        gi(end+1)  = S.PP.colorIdx(k); %#ok<SAGROW>
        labs{end+1}= lab;              %#ok<SAGROW>
    catch ME
        warning('  %s: %s', lab, ME.message);
    end
end
n = size(H,1);
fprintf('Overlaying %d sig ROIs\n', n);
if n==0, error('no sig ROIs collected'); end

%% ===================== LINEAR OVERLAY =====================
fig = figure('Color','w','Name','sig-ROI phase histogram overlay (linear)', ...
             'Units','centimeters','Position',[3 3 20 13]);
ax = axes(fig); hold(ax,'on');
for i = 1:n
    plot(ax, ctrsP, H(i,:), '-', 'Color',[S.group_colors(gi(i),:) 0.55], 'LineWidth',lineW);
end
if showMean
    ug = unique(gi);
    for j = 1:numel(ug)
        m = mean(H(gi==ug(j),:), 1);
        plot(ax, ctrsP, m, '-', 'Color', S.group_colors(ug(j),:), 'LineWidth', 3);
    end
end
xline(ax, pi, 'k--', 'LineWidth', 0.8);          % breath peak
xlim(ax,[0 2*pi]); set(ax,'XTick',[0 pi 2*pi],'XTickLabel',{'0 (onset)','\pi (peak)','2\pi'});
xlabel(ax,'breath phase'); ylabel(ax,'spikes/cycle (%)'); box(ax,'on');

ug = unique(gi); hL = gobjects(numel(ug),1); leg = cell(numel(ug),1);
for j = 1:numel(ug)
    hL(j) = plot(ax, nan, nan, '-', 'Color', S.group_colors(ug(j),:), 'LineWidth', 2);
    leg{j} = sprintf('%s (n=%d)', S.groups{ug(j)}, sum(gi==ug(j)));
end
legend(ax, hL, leg, 'Location','northeastoutside','Box','off','Interpreter','none','FontSize',9);
title(ax, sprintf('significant-ROI spike-phase histograms (occ-comp spk/cycle %%)  |  N=%d  |  onset=0, peak=\\pi', n), ...
      'Interpreter','tex','FontWeight','bold');

if doSave
    base = fullfile(outDir, 'phase_hist_sigROI_overlay_linear');
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    groups = S.groups; group_colors = S.group_colors;
    save([base '.mat'], 'H','gi','ctrsP','labs','groups','group_colors');
    fprintf('Saved %s.png/.pdf\n', base);
end

%% ===================== OCCUPANCY (DWELL) OVERLAY =====================
% The phase-dwell profile each FOV's compensation divides out: fraction of the
% cycle spent in each phase bin. Two-level step -- low in inspiration (0..pi),
% high in expiration (pi..2pi). (Per-FOV; same-FOV sig ROIs overlap.)
figO = figure('Color','w','Name','sig-ROI phase occupancy (dwell)', ...
              'Units','centimeters','Position',[3 3 20 13]);
axO = axes(figO); hold(axO,'on');
for i = 1:n
    plot(axO, ctrsP, 100*O(i,:), '-', 'Color',[S.group_colors(gi(i),:) 0.55], 'LineWidth',lineW);
end
if showMean
    ug = unique(gi);
    for j = 1:numel(ug)
        plot(axO, ctrsP, 100*mean(O(gi==ug(j),:),1), '-', 'Color', S.group_colors(ug(j),:), 'LineWidth', 3);
    end
end
yline(axO, 100/nPhaseBins, 'k:', 'LineWidth', 0.8);   % flat-dwell reference (uniform)
xline(axO, pi, 'k--', 'LineWidth', 0.8);              % breath peak
xlim(axO,[0 2*pi]); set(axO,'XTick',[0 pi 2*pi],'XTickLabel',{'0 (onset)','\pi (peak)','2\pi'});
xlabel(axO,'breath phase'); ylabel(axO,'dwell (% of cycle per bin)'); box(axO,'on');
ug = unique(gi); hLO = gobjects(numel(ug),1); legO = cell(numel(ug),1);
for j = 1:numel(ug)
    hLO(j) = plot(axO, nan, nan, '-', 'Color', S.group_colors(ug(j),:), 'LineWidth', 2);
    legO{j} = sprintf('%s (n=%d)', S.groups{ug(j)}, sum(gi==ug(j)));
end
legend(axO, hLO, legO, 'Location','northeastoutside','Box','off','Interpreter','none','FontSize',9);
title(axO, sprintf('phase dwell profile (the occ\\_b divided out)  |  N=%d  |  dotted = uniform (%.1f%%)', ...
      n, 100/nPhaseBins), 'Interpreter','tex','FontWeight','bold');
if doSave
    base = fullfile(outDir, 'phase_occupancy_sigROI_overlay_linear');
    exportgraphics(figO, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(figO, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'O','gi','ctrsP','labs');
    fprintf('Saved %s.png/.pdf\n', base);
end

%% ===================== LOCAL FUNCTIONS (verbatim from phase_polar_allROI) =====================
function F = load_folder(folderPath, nDrop, fallback_fps)
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
assert(~isempty(bp),'no breath_peak_pc1'); assert(~isempty(ip),'no breath_insp_start_pc1');
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
if contains(folderPath, fullfile('Vglut2','1124'))   % rising-edge: breath leads 1 frame
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
F.ev = ev; F.ef = ef; F.roi_spikes = CA.roi_spikes;
F.nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); F.nDff = size(D.dFF,1); end
end

function [pctP, nspk, occFrac] = roi_phase_hist(F, roi, ca_lag_sec, edgesP)
spk = double(F.roi_spikes(roi).spike_train(:));
T = min([F.nDff, numel(F.ev), numel(F.ef), numel(spk)]);
ev = F.ev(1:T); ef = F.ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);
phi = piecewise_phase_local(find(ev>0), find(ef>0), T);
lag = round(ca_lag_sec*F.fps);
sp = find(spk>0) - lag; sp = sp(sp>=1 & sp<=T);
spk_phi = mod(phi(sp), 2*pi); spk_phi = spk_phi(~isnan(spk_phi));
phv = phi(~isnan(phi)); nCyc = (max(phv)-min(phv))/(2*pi);
occ = histcounts(mod(phv,2*pi), edgesP); Ob = mean(occ(occ>0));      % phase occupancy (dwell)
wOcc = ones(1,numel(occ)); wOcc(occ>0) = Ob ./ occ(occ>0);
pctP = (100 * histcounts(spk_phi, edgesP) / max(nCyc,eps)) .* wOcc;  % occupancy-compensated spk/cyc %
occFrac = occ / max(sum(occ), eps);                                  % dwell as fraction of cycle
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
