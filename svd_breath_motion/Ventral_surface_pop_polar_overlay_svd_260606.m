% Ventral_surface_pop_polar_overlay_svd_260606.m
% -----------------------------------------------------------------------
%  Population spike-phase polar overlay for three ROI sets:
%     sig            : PP.r >= confC                      (significant)
%     IO             : PP.colorIdx == IO group            (inferior olive)
%     nonIO_nonsig   : PP.colorIdx ~= IO  &  PP.r < confC (the rest)
%
%  TWO polar panels (both in spk/cyc % units):
%   (1) raw spike-count distribution      : 100 * N_b / nCyc   ("spk/cyc %"; what
%       every other phase histogram in this folder plots).
%   (2) occupancy-COMPENSATED distribution: same units, each bin reweighted by
%       relative dwell time:  sc(b) * (Obar / O_b), where O_b = frames the breath
%       spends in phase bin b and Obar = mean occupancy. This removes the unequal
%       dwell of inspiration vs expiration under the piecewise onset=0/peak=pi
%       phase -- a uniform-firing neuron becomes FLAT. STILL spk/cyc %, not a rate.
%
%  N_b = spikes in phase bin (lead-shifted by ca_lag_sec), binned on the EXACT
%  piecewise phase (piecewise_phase_local). Vglut2/1124 +1-frame fix in the loader.
%
%  Source: coherence_polar_svd_260606/coherence_polar_data.mat
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot = fileparts(scriptDir);
addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence'));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_svd_260606', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'temporal_phase_batch_svd_260606');

nDrop        = 30;     fallback_fps = 30;
ca_lag_sec   = 0.1;    nPhaseBins   = 24;   % lead-comp: spikes 3 frames earlier @30Hz (0.1 s)
normShape    = false;  % true = scale each curve to its own max (pure shape compare, both panels)

% colors / line widths per set
sigCol = [0 0 0];           sigLW = 2.6;     % significant   = black
ioCol  = [0.10 0.45 0.85];  ioLW  = 2.2;     % IO            = blue
nsCol  = [0.60 0.60 0.60];  nsLW  = 2.2;     % nonIO-nonsig  = gray
onsetCol = [0.90 0.10 0.10]; peakCol = [0.35 0.75 1.00];   % 0 = onset (red), pi = peak (sky)

doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
edgesP = linspace(0, 2*pi, nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ===================== LOAD PP + folder map =====================
S = load(cohData, 'PP','confC','labels','groups');
ioGroup = find(strcmp(S.groups,'IO'));
assert(~isempty(ioGroup), 'no IO group in coherence data');
fprintf('confC=%.4f | %d ROIs total\n', S.confC, numel(S.PP.r));

allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ===================== THREE SETS =====================
sets(1) = struct('name','sig',          'sel', find(S.PP.r >= S.confC),                        'col',sigCol,'lw',sigLW);
sets(2) = struct('name','IO',           'sel', find(S.PP.colorIdx==ioGroup),                   'col',ioCol, 'lw',ioLW);
sets(3) = struct('name','nonIO_nonsig', 'sel', find(S.PP.colorIdx~=ioGroup & S.PP.r<S.confC), 'col',nsCol, 'lw',nsLW);

Fcache = containers.Map();
meanRaw = cell(numel(sets),1);   % raw spk/cyc %
meanCmp = cell(numel(sets),1);   % occupancy-compensated spk/cyc % (same units)
nSets   = zeros(numel(sets),1);
for s = 1:numel(sets)
    sel = sets(s).sel; RAW = []; CMP = [];
    for k = 1:numel(sel)
        lab = S.labels{sel(k)}; parts = regexp(lab,'/','split');
        if numel(parts)<4, continue; end
        recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
        if ~isKey(folderMap, recName), warning('no folder: %s', lab); continue; end
        folderPath = folderMap(recName);
        try
            if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_folder(folderPath, nDrop, fallback_fps); end
            R = roi_phasehist(Fcache(folderPath), roi, ca_lag_sec, edgesP);
            if ~R.ok, continue; end
            sc = 100 * R.cntW / max(R.nCyc, eps);          % spk/cyc % (raw)
            nz = R.occ > 0; Obar = mean(R.occ(nz));        % mean dwell (frames/bin)
            w  = ones(1, numel(R.occ)); w(nz) = Obar ./ R.occ(nz);   % dwell reweight (dimensionless)
            cc = sc .* w;                                  % occupancy-compensated spk/cyc % (same units)
            RAW(end+1,:) = sc; CMP(end+1,:) = cc;          %#ok<SAGROW>
        catch ME
            warning('  ERROR %s: %s', lab, ME.message);
        end
    end
    nSets(s) = size(RAW,1);
    meanRaw{s} = mean(RAW,1); meanCmp{s} = mean(CMP,1);
    if normShape
        if max(meanRaw{s})>0, meanRaw{s} = meanRaw{s}/max(meanRaw{s}); end
        if max(meanCmp{s})>0, meanCmp{s} = meanCmp{s}/max(meanCmp{s}); end
    end
    fprintf('  %-14s : %d ROIs\n', sets(s).name, nSets(s));
end

%% ===================== OVERLAY (2 panels, both spk/cyc %) =====================
fig = figure('Color','w','Name','population spike-phase: raw vs occupancy-compensated', ...
             'Units','normalized','Position',[0.06 0.12 0.84 0.72]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');
unit = ternary(normShape, 'norm.', 'spk/cyc %');
draw_overlay(tl, 1, sets, meanRaw, nSets, ctrsP, onsetCol, peakCol, ...
    sprintf('(1) raw distribution  [%s]', unit));
draw_overlay(tl, 2, sets, meanCmp, nSets, ctrsP, onsetCol, peakCol, ...
    sprintf('(2) occupancy-compensated  [%s]', unit));
title(tl, 'population spike-phase  |  onset=0 (red), peak=\pi (sky)  |  spk/cyc %:  N_b/nCyc  vs  (N_b/nCyc)\cdot(Obar/O_b)', ...
      'Interpreter','tex','FontWeight','bold');

%% ===================== SAVE =====================
if doSave
    setNames = {sets.name};
    base = fullfile(outDir, 'pop_polar_overlay');
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'meanRaw','meanCmp','nSets','setNames','ctrsP','edgesP','normShape','ca_lag_sec');
    fprintf('Saved %s.png/.pdf/.mat\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end

function draw_overlay(tl, tile, sets, means, nSets, ctrsP, onsetCol, peakCol, ttl)
Rmax = 0;
for s = 1:numel(sets), if ~isempty(means{s}), Rmax = max(Rmax, max(means{s})); end, end
Rmax = pick_rmax(Rmax);
pax = polaraxes(tl); pax.Layout.Tile = tile; hold(pax,'on');
polarplot(pax, [0 0],   [0 Rmax], '--', 'Color', onsetCol, 'LineWidth', 1.2);
polarplot(pax, [pi pi], [0 Rmax], '--', 'Color', peakCol,  'LineWidth', 1.2);
hLeg = gobjects(numel(sets),1); legTxt = cell(numel(sets),1);
for s = 1:numel(sets)
    mc = means{s}; if isempty(mc), continue; end
    hLeg(s) = polarplot(pax, [ctrsP ctrsP(1)], [mc mc(1)], '-', 'Color', sets(s).col, 'LineWidth', sets(s).lw);
    legTxt{s} = sprintf('%s (n=%d)', sets(s).name, nSets(s));
end
pax.RLim=[0 Rmax]; pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise';
pax.RAxisLocation=180; pax.FontSize=9; thetaticks(pax,0:90:270);
title(pax, ttl, 'Interpreter','none', 'FontWeight','bold');
ok = isgraphics(hLeg);
legend(pax, hLeg(ok), legTxt(ok), 'Location','southoutside', 'Orientation','horizontal', 'Box','off');
end

function R = pick_rmax(dataMax)
ladder = [5 10 20 40 60 80 100 150 200];
if ~isfinite(dataMax) || dataMax <= 0, R = 10; return; end
idx = find(ladder >= dataMax, 1);
if isempty(idx), R = ceil(dataMax/50)*50; else, R = ladder(idx); end
end

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
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;  % PEAK
ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;   % ONSET
ef(1:min(nDrop,numel(ef))) = [];
if contains(folderPath, fullfile('Vglut2','1124'))   % rising-edge +1 fix
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
F.ev = ev; F.ef = ef; F.roi_spikes = CA.roi_spikes; F.nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); F.nDff = size(D.dFF,1); end
end

function R = roi_phasehist(F, roi, ca_lag_sec, edgesP)
R = struct('cntW',[],'occ',[],'nCyc',0,'fps',F.fps,'ok',false);
if roi < 1 || roi > numel(F.roi_spikes), return; end
spk = double(F.roi_spikes(roi).spike_train(:));
T = min([F.nDff, numel(F.ev), numel(F.ef), numel(spk)]);
ev = F.ev(1:T); ef = F.ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);
phi = piecewise_phase_local(find(ev>0), find(ef>0), T);
lag = round(ca_lag_sec*F.fps);
sp  = find(spk>0) - lag; sp = sp(sp>=1 & sp<=T);
spk_phi = mod(phi(sp), 2*pi); spk_phi = spk_phi(~isnan(spk_phi));
phv = phi(~isnan(phi));
R.cntW = histcounts(spk_phi, edgesP);         % spikes per phase bin
R.occ  = histcounts(mod(phv,2*pi), edgesP);   % FRAMES per phase bin (occupancy)
R.nCyc = (max(phv)-min(phv))/(2*pi);
R.ok   = true;
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
