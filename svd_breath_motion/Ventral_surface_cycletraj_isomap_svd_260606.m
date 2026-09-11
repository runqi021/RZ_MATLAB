% Ventral_surface_cycletraj_isomap_svd_260606.m
% -----------------------------------------------------------------------
%  TEMPORAL-TRACE trajectory via ISOMAP on a delay/Hankel embedding.
%  NOT a neuron state-space -- this embeds the 1-D temporal TRACE itself so
%  the breath cycle appears as a nonlinear attractor LOOP (kNN graph ->
%  geodesics -> classical MDS).
%
%  Only the CALCIUM trace is dimension-reduced. Breathing motion is NOT embedded
%  -- its cycle is already defined by the 0/pi breath phase, so embedding it would
%  be circular.
%     calcium = mean dF/F over this FOV's SIGNIFICANT ROIs (PP.r >= confC),
%               shifted earlier by the data-driven GCaMP lead (ca_lag_sec -> 7 fr @30Hz)
%  The trace is z-scored, delay-embedded (window = one median breath cycle), the
%  Hankel rows are uniformly downsampled to maxIsoPts, Isomap-embedded to 3-D, and
%  colored by piecewise breath phase (onset=0 red, peak=pi sky).
%
%  Parallel PCA version: Ventral_surface_cycletraj_pca_svd_260606.m
%  Source FOV list: coherence_polar_svd_260606/coherence_polar_data.mat (labels, PP, confC).
%  Dependencies: detect_session_fps.m, Statistics TB (cmdscale, graph).
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot = fileparts(scriptDir);
addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence'));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_svd_260606', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'traj_isomap_svd_260606');

doFolders    = {};        % {} = every FOV with >=1 significant ROI; else list recNames
nDrop        = 30;        fallback_fps = 30;
ca_lag_sec   = 0.1;       % calcium lead: shift trace EARLIER 3 frames @30Hz (0.1 s)
useSpikes    = true;      % trace = deconvolved spike probability (Gaussian-smoothed spike train), NOT dF/F -> avoids breath-motion artifact in dF/F
spk_smooth_sec = 0.1;     % Gaussian sigma (s) for the spike-probability trace
win_cycles   = 1.0;       % delay-embedding window = win_cycles x median breath period
popCycles    = 50;        % population trajectory: phase-aligned breath cycles taken per FOV
kIso         = 6;         % Isomap neighbours (auto-raised until graph connected)
maxIsoPts    = 1500;      % uniform-downsample Hankel rows to this many points
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
nPhaseBins = 24; phiBin = (0:nPhaseBins-1)/nPhaseBins*2*pi;   % phase of each resampled bin
if useSpikes, traceLbl = 'spike prob'; else, traceLbl = 'dF/F'; end %#ok<UNRCH>

%% ===================== SIG-ROI LIST + FOLDER MAP =====================
assert(isfile(cohData),'coherence_polar_data.mat not found: %s', cohData);
S = load(cohData,'PP','confC','labels');
sigMask = S.PP.r >= S.confC;
sigList = cell(0,2);   % each row = {recName, roi}  (one significant ROI)
for k = 1:numel(S.labels)
    if ~sigMask(k), continue; end
    parts = regexp(S.labels{k},'/','split'); if numel(parts)<4, continue; end
    sigList(end+1,:) = {strjoin(parts(3:end-1),'/'), str2double(parts{end})}; %#ok<SAGROW>
end
if ~isempty(doFolders), sigList = sigList(ismember(sigList(:,1), doFolders), :); end
allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end
fprintf('%d significant ROIs\n', size(sigList,1));

%% ===================== LOOP SIG ROIs =====================
Fcache = containers.Map(); nOK = 0; M = []; labLong = {};
for r = 1:size(sigList,1)
    recName = sigList{r,1}; roi = sigList{r,2};
    if ~isKey(folderMap, recName), warning('no folder for %s', recName); continue; end
    folderPath = folderMap(recName);
    try
        if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_fov(folderPath, nDrop, fallback_fps); end
        F = Fcache(folderPath);
        if roi<1 || roi>size(F.dFF,2) || roi>numel(F.roi_spikes), warning('roi %d out of range in %s', roi, recName); continue; end

        % ---- per-ROI trace (z-scored): spike probability avoids dF/F breath-motion artifact ----
        T = min([size(F.dFF,1), numel(F.phi)]);
        if useSpikes
            st = double(F.roi_spikes(roi).spike_train(:));
            if numel(st) < T, st(end+1:T) = 0; end
            cal0 = zsc(smooth_gauss(st(1:T), round(spk_smooth_sec*F.fps)));
        else
            cal0 = zsc(F.dFF(1:T, roi)); %#ok<UNRCH>
        end
        phi  = F.phi(1:T);
        lag  = round(ca_lag_sec * F.fps);

        % ---- this ROI's popCycles best-breath-correlated cycles, phase-resampled [popCycles x nBins] ----
        [cycMat, ok] = roi_cycle_matrix(cal0, F.breath, phi, T, lag, nPhaseBins, popCycles);
        if ~ok, warning('  %s ROI%d: too few cycles', recName, roi); continue; end

        % ---- per-ROI trajectory = Isomap over its cycles (embed the phase bins) ----
        [Yc, kc, rc, keep] = isomap_traj(cycMat, kIso, 3);
        fig = figure('Color','w','Name',sprintf('%s ROI%d traj Isomap',recName,roi), ...
                     'Units','centimeters','Position',[3 3 17 14]);
        ax = axes('Parent',fig);
        draw_loop(ax, Yc, phiBin(keep), sprintf('%s ROI%d  |  Isomap traj over %d breath-corr cycles  |  %s (+%dfr)  k=%d res=%.2f', ...
                  recName, roi, popCycles, traceLbl, lag, kc, rc(min(3,numel(rc)))));
        if doSave
            base = fullfile(outDir, sprintf('traj_isomap_%s_roi%02d', matlab.lang.makeValidName(recName), roi));
            exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
            exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            save([base '.mat'], 'Yc','rc','kc','keep','cycMat','popCycles','roi','lag','recName');
        end

        % ---- accumulate into the POPULATION matrix (flatten the cycles) ----
        M(end+1,:) = reshape(cycMat.', 1, []); labLong{end+1} = sprintf('%s/roi%d',recName,roi); %#ok<SAGROW>

        nOK = nOK + 1;
        fprintf('[%d/%d] %s ROI%d  %d cycles  k=%d res3=%.2f\n', r, size(sigList,1), recName, roi, popCycles, kc, rc(min(3,numel(rc))));
    catch ME
        warning('  %s ROI%d: %s', recName, roi, ME.message);
    end
end
fprintf('\nDone. %d ROI figures. Output: %s\n', nOK, outDir);

%% ===================== POPULATION TRAJECTORY (all sig ROIs, phase-binned dF/F) =====================
if size(M,1) >= 3
    phiObs = repmat((0:nPhaseBins-1)/nPhaseBins*2*pi, 1, popCycles);   % within-cycle phase per obs
    [Yi, ki, ri, keepB] = isomap_traj(M, kIso, 3);   % M = [Nsig x popCycles*nBins]; embed the obs
    figp = figure('Color','w','Name','population traj Isomap','Units','centimeters','Position',[3 3 18 15]);
    axp = axes('Parent',figp);
    draw_loop(axp, Yi, phiObs(keepB), sprintf('POPULATION Isomap traj  |  %d sig ROIs x %d top breath-corr cycles  k=%d res=%.2f', ...
              size(M,1), popCycles, ki, ri(min(3,numel(ri)))));
    if doSave
        base = fullfile(outDir,'traj_isomap_POPULATION');
        exportgraphics(figp, [base '.png'], 'Resolution',220, 'BackgroundColor','white');
        exportgraphics(figp, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        save([base '.mat'], 'Yi','ri','ki','keepB','M','labLong','phiObs','popCycles','nPhaseBins');
    end
    fprintf('POPULATION Isomap: %d sig ROIs x %d cycles (%d obs), k=%d res3=%.2f\n', size(M,1), popCycles, size(M,2), ki, ri(min(3,numel(ri))));
end

%% ===================== LOCAL FUNCTIONS =====================
function [Y, kUsed, resid, keep] = isomap_traj(M, k0, nDims)
% Population trajectory: M = [neurons x bins]; embed the BINS via Isomap.
A = M'; sd = std(A,0,1); sd(sd==0)=1; A = (A - mean(A,1)) ./ sd;   % [bins x neurons]
n = size(A,1); D = squareform(pdist(A));
kUsed = min(k0, n-1);
while true
    if max(conncomp(graph(knn_adj(D,kUsed))))==1 || kUsed>=n-1, break; end
    kUsed = kUsed + 1;
end
G = graph(knn_adj(D,kUsed)); comp = conncomp(G); keep = find(comp == mode(comp));
Dg = distances(subgraph(G, keep));
[Yf,~] = cmdscale(Dg); nd = min(nDims, size(Yf,2)); Y = Yf(:,1:nd); if nd<3, Y(:,nd+1:3) = 0; end
mask = tril(true(numel(keep),numel(keep)),-1); dgv = Dg(mask); resid = zeros(1,nd);
for d = 1:nd, Dl = squareform(pdist(Yf(:,1:d))); resid(d) = 1 - corr(dgv, Dl(mask))^2; end
end

function [cycMat, ok] = roi_cycle_matrix(cal0, breath, phi, T, lag, nBins, popCycles)
% [popCycles x nBins] phase-resampled trace for the popCycles breath cycles whose
% (lead-shifted) trace is best correlated with the breathing motion. NaN-safe.
ok = false; cycMat = [];
calU = cal0(1+lag:T); brU = breath(1:T-lag);
brU = (brU - mean(brU,'omitnan')) / max(std(brU,'omitnan'),eps);
phiU = phi(1:T-lag);
g = ~isnan(phiU) & ~isnan(calU) & ~isnan(brU); phiU=phiU(g); calU=calU(g); brU=brU(g);
[phiU, iu] = unique(phiU); calU = calU(iu); brU = brU(iu);
if numel(phiU) <= nBins, return; end
firstOn = ceil(phiU(1)/(2*pi))*2*pi; lastOn = floor(phiU(end)/(2*pi))*2*pi;
nCyc = max(1, round((lastOn-firstOn)/(2*pi)));
ccorr = nan(nCyc,1);                               % per-cycle trace<->breath correlation
for c = 1:nCyc
    fr = phiU>=firstOn+(c-1)*2*pi & phiU<firstOn+c*2*pi;
    if nnz(fr) >= 4, ccorr(c) = corr(calU(fr), brU(fr)); end
end
if nCyc < popCycles, return; end                   % need >= popCycles cycles for equal-length stacking
[~,ord] = sort(ccorr,'descend','MissingPlacement','last');  % no-spike (NaN-corr) cycles rank last
pick = sort(ord(1:popCycles));                     % top breath-correlated cycles
ph1 = (0:nBins-1)/nBins*2*pi;
cycMat = nan(popCycles, nBins);
for j = 1:popCycles
    cycMat(j,:) = interp1(phiU, calU, firstOn+(pick(j)-1)*2*pi+ph1, 'linear', NaN);
end
cycMat = fillmissing(cycMat, 'linear', 2, 'EndValues','nearest');
ok = true;
end

function A = knn_adj(D, k)
n = size(D,1); A = zeros(n); [~,ord] = sort(D,2);
for i = 1:n, nb = ord(i,2:min(k+1,n)); A(i,nb) = D(i,nb); end
A = max(A, A');
end

function draw_loop(ax, Y, phaseRows, ttl)
hold(ax,'on'); cmap = hsv(256);
plot3(ax, Y(:,1), Y(:,2), Y(:,3), '-', 'Color',[0.8 0.8 0.8], 'LineWidth',0.3);
cidx = max(1, ceil(phaseRows(:)/(2*pi)*256));
scatter3(ax, Y(:,1), Y(:,2), Y(:,3), 14, cmap(min(cidx,256),:), 'filled', 'MarkerFaceAlpha',0.65);
grid(ax,'on'); box(ax,'on'); view(ax,30,25); axis(ax,'tight');
xlabel(ax,'Isomap 1'); ylabel(ax,'Isomap 2'); zlabel(ax,'Isomap 3');
colormap(ax, hsv(256)); caxis(ax,[0 2*pi]);
cb = colorbar(ax); cb.Label.String='breath phase (rad)'; cb.Ticks=[0 pi 2*pi]; cb.TickLabels={'0','\pi','2\pi'};
title(ax, ttl, 'Interpreter','none');
end

function y = zsc(x)
x = x(:); s = std(x,'omitnan'); if s==0, s = 1; end
y = (x - mean(x,'omitnan')) / s;
end

function y = smooth_gauss(x, sig)
% Gaussian smoothing of a (sparse) spike train -> continuous spike-probability rate.
sig = max(1, sig); w = ceil(3*sig); k = exp(-(-w:w).^2/(2*sig^2)); k = k/sum(k);
y = conv(x(:), k(:), 'same');
end

function F = load_fov(folderPath, nDrop, fallback_fps)
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
assert(~isempty(bp),'no breath_peak_pc1'); assert(~isempty(ip),'no breath_insp_start_pc1');
assert(~isempty(df),'no *_ch1_dFF.mat');
F.fps = detect_session_fps(folderPath, fallback_fps);
BP = load(fullfile(bp(1).folder, bp(1).name));
IP = load(fullfile(ip(1).folder, ip(1).name));
D  = load(fullfile(df(1).folder, df(1).name),'dFF');
F.dFF = double(D.dFF);
assert(isfile(fullfile(folderPath,'ca_spike_data.mat')),'no ca_spike_data');
CA = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes'); F.roi_spikes = CA.roi_spikes;
nB = numel(BP.breath);
br = detrend(double(BP.breath(:))); br(1:min(nDrop,numel(br))) = []; F.breath = br;
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1; ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1; ef(1:min(nDrop,numel(ef))) = [];
% Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay breath.
if contains(folderPath, fullfile('Vglut2','1124'))
    F.breath = [F.breath(1); F.breath(1:end-1)];
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
Tp = min([numel(F.breath), numel(ev), numel(ef), size(F.dFF,1)]);
peakIdx = find(ev(1:Tp)>0); footIdx = find(ef(1:Tp)>0);
F.phi = piecewise_phase_local(peakIdx, footIdx, Tp);
if numel(footIdx) >= 2, F.medPeriod = median(diff(footIdx)); else, F.medPeriod = round(F.fps); end
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
