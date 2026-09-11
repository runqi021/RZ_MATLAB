% phase_polar_allROI.m
% -----------------------------------------------------------------------
%  ALL-ROI polar overlay, in the spirit of Ventral_surface_coherence_polar_svd_260606,
%  but using MY phase (the spike-phase histogram), NOT the coherence phase.
%
%  Every ROI from every FOV is drawn on ONE polar axes as a CONNECTING-LINE
%  histogram of its per-cycle spike-phase distribution (spikes/cycle %),
%  colored by GROUP.  No coherence, no dots -- just stack them all and look.
%
%  Convention: onset = 0 (theta right), peak = pi, counter-clockwise.
%
%  Source: coherence_polar_svd_260606/coherence_polar_data.mat
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot = fileparts(fileparts(scriptDir)); addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence')); addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
addpath(scriptDir);                                   % so coh_cfg_260727 resolves

%% ===================== USER-EDITABLE =====================
cfg      = coh_cfg_260727();
rootPath = cfg.rootPath;
cohData  = cfg.cohData;
outDir   = fullfile(cfg.outRoot, 'phase_polar_allROI');

nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.1; nPhaseBins   = 24;   % lead-comp: spikes 3 frames earlier @30Hz (0.1 s)
minSpikes    = 1;
lineW        = 0.5;    % per-ROI line width
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== LOAD PP + folder map =====================
S = load(cohData, 'PP','confC','labels','groups','group_colors');
allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ===================== COLLECT ALL ROIs =====================
P = []; gi = []; Fcache = containers.Map();
for k = 1:numel(S.PP.r)
    lab = S.labels{k}; parts = regexp(lab,'/','split');
    if numel(parts)<4, continue; end
    recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
    if ~isKey(folderMap, recName), continue; end
    folderPath = folderMap(recName);
    try
        if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_folder(folderPath, nDrop, fallback_fps); end
        [pctP, nspk] = roi_phase_hist(Fcache(folderPath), roi, ca_lag_sec, edgesP);
        if nspk < minSpikes, continue; end
        P(end+1,:) = pctP;              %#ok<SAGROW>
        gi(end+1)  = S.PP.colorIdx(k);  %#ok<SAGROW>
    catch ME
        warning('  %s: %s', lab, ME.message);
    end
end
n = size(P,1);
fprintf('Stacking %d ROIs\n', n);

%% ===================== POLAR OVERLAY (full scale + RLim 20) =====================
fig = figure('Color','w','Name','all-ROI phase polar overlay', ...
             'Units','centimeters','Position',[3 3 30 14]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
rmax = max(P(:));

% (1) full scale
pax1 = polaraxes(tl); pax1.Layout.Tile = 1;
draw_overlay(pax1, P, gi, ctrsP, S.group_colors, lineW);
if isfinite(rmax)&&rmax>0, pax1.RLim=[0 rmax]; end
title(pax1, sprintf('full scale (max=%.0f)', rmax),'Interpreter','none');

% (2) RLim capped to 20
pax2 = polaraxes(tl); pax2.Layout.Tile = 2;
draw_overlay(pax2, P, gi, ctrsP, S.group_colors, lineW);
pax2.RLim = [0 10];
title(pax2, 'RLim = [0 10]','Interpreter','none');

% group legend with counts (on the right axes)
ug = unique(gi); hL = gobjects(numel(ug),1); labs = cell(numel(ug),1);
for j = 1:numel(ug)
    hL(j) = polarplot(pax2, nan, nan, '-', 'Color', S.group_colors(ug(j),:), 'LineWidth', 2);
    labs{j} = sprintf('%s (n=%d)', S.groups{ug(j)}, sum(gi==ug(j)));
end
legend(pax2, hL, labs, 'Location','northeastoutside','Box','off','Interpreter','none','FontSize',8);
title(tl, sprintf('all-ROI spike-phase histograms (spikes/cycle %%)  |  N=%d  |  onset=0, peak=\\pi', n), ...
      'Interpreter','tex','FontWeight','bold');

if doSave
    base = fullfile(outDir, 'phase_polar_allROI');
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    groups = S.groups; group_colors = S.group_colors; %#ok<NASGU>
    save([base '.mat'], 'P','gi','ctrsP','groups','group_colors');
    fprintf('Saved %s.png/.pdf\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function draw_overlay(pax, P, gi, ctrsP, group_colors, lineW)
hold(pax,'on'); thc = [ctrsP ctrsP(1)];
for i = 1:size(P,1)
    polarplot(pax, thc, [P(i,:) P(i,1)], '-', 'Color', group_colors(gi(i),:), 'LineWidth', lineW);
end
pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise'; pax.RAxisLocation=180; pax.FontSize=9;
thetaticks(pax,0:45:315);
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

function [pctP, nspk] = roi_phase_hist(F, roi, ca_lag_sec, edgesP)
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
