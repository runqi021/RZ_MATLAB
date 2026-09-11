% spike_trigger_dFF.m
% -----------------------------------------------------------------------
%  Spike-triggered average of dF/F for ALL ACTIVE ROIs, to characterize the
%  calcium transient RISE TIME (so the GCaMP rise can be measured empirically
%  instead of assumed).  Pure calcium: dFF [T x N] + ca-spike train only --
%  no breath, no ca_lag, no rising-edge shift involved.
%
%  Active ROI = more than minEvents spike events (nnz(spike_train>0) > minEvents),
%  matching Ventral_surface_temporal_phase_activeROI_svd_260606.m.
%
%  Per ROI: extract dF/F in a +/- win_sec window around each spike, average,
%  then baseline the dF/F (subtract the pre-spike window mean so it -> 0).
%  Plot: baselined dF/F STA per ROI (gray) + population mean +/- SD, plus a
%  HEATMAP of per-ROI STAs (rows = active ROIs, gray cmap, sorted by amplitude).
%  10-90% rise time is computed from a peak-normalized copy (internal only).
%
%  Source ROI list: coherence_polar_svd_260606/coherence_polar_data.mat (labels).
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(scriptDir);                                   % so coh_cfg_260727 resolves

%% ===================== USER-EDITABLE =====================
cfg      = coh_cfg_260727();
rootPath = cfg.rootPath;
cohData  = cfg.cohData;
outDir   = fullfile(cfg.outRoot, 'spike_triggered_dFF');

ca_lag = 0.10;
fallback_fps = 30;
minEvents    = 5;          % ACTIVE ROI = more than this many spike events (raw nnz)
trainISI_frames = 2;       % spikes with ISI <= this (frames) are one TRAIN; keep only the first
win_sec      = 2;          % +/- spike-triggered window (s)  -> clip to [-1 +1]
gridFps      = 30;         % common time grid for cross-session averaging (Hz)
baseWin      = [-1.0 -0.5];% baseline window for subtraction (s, pre-spike)
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
tg = -win_sec : 1/gridFps : win_sec;        % common time grid (s)

%% ===================== LIST + FOLDER MAP =====================
assert(isfile(cohData), 'coherence_polar_data.mat not found: %s', cohData);
S = load(cohData, 'labels','PP','confC');
allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

%% ===================== COLLECT SPIKE-TRIGGERED dF/F =====================
Fcache = containers.Map();
STAn = []; STAr = []; ttp = []; rise1090 = []; ampl = []; nspkAll = []; labKept = {};
spkSnips = {};   % per-ROI [nspike x numel(tg)] baselined individual-spike snippets
sigKept  = [];   % per-kept-ROI flag: is this ROI significant (PP.r >= confC)?
for k = 1:numel(S.labels)
    lab = S.labels{k}; parts = regexp(lab,'/','split');
    if numel(parts) < 4, continue; end
    recName = strjoin(parts(3:end-1),'/'); roi = str2double(parts{end});
    if ~isKey(folderMap, recName), continue; end
    folderPath = folderMap(recName);
    try
        if ~isKey(Fcache, folderPath), Fcache(folderPath) = load_ca(folderPath, fallback_fps); end
        F = Fcache(folderPath);
        if roi > size(F.dFF,2) || roi > numel(F.roi_spikes), continue; end
        spk = double(F.roi_spikes(roi).spike_train(:));
        Tn  = min(size(F.dFF,1), numel(spk));
        dff = F.dFF(1:Tn, roi); spk = spk(1:Tn);
        if nnz(spk>0) <= minEvents, continue; end          % ACTIVE gate

        win = round(win_sec * F.fps); tau = (-win:win)/F.fps;
        sf  = find(spk>0);
        if ~isempty(sf), sf = sf([true; diff(sf) > trainISI_frames]); end  % collapse trains -> first spike
        sf = sf(sf-win>=1 & sf+win<=Tn);                    % edge-valid spikes
        if numel(sf) < 2, continue; end
        E = zeros(numel(sf), 2*win+1);
        for j = 1:numel(sf), E(j,:) = dff(sf(j)-win : sf(j)+win); end
        sta = mean(E,1);

        bsel = tau>=baseWin(1) & tau<baseWin(2); if ~any(bsel), bsel = tau<0; end
        sta = sta - mean(sta(bsel));                        % baseline the dF/F (pre-spike window -> 0)
        sr  = interp1(tau, sta, tg, 'linear', NaN);         % baselined dF/F on common grid  -> PANEL 2
        pk  = max(sr(tg>=0));
        if ~isfinite(pk) || pk <= 0, continue; end
        sn  = sr / pk;                                      % peak-normalized: baseline 0 -> peak 1  -> PANEL 1

        post = tg>=0; tgp = tg(post); srp = sr(post);
        [~,im] = max(srp); tpk = tgp(im);                  % time-to-peak (s)
        r10 = rise_crossing(tg, sn, 0.1, tpk);
        r90 = rise_crossing(tg, sn, 0.9, tpk);

        STAr(end+1,:) = sr;  STAn(end+1,:) = sn;            %#ok<SAGROW>
        ttp(end+1)    = tpk; rise1090(end+1) = r90 - r10;  %#ok<SAGROW>
        ampl(end+1)   = pk;  nspkAll(end+1) = numel(sf);   %#ok<SAGROW>
        labKept{end+1}= lab;                               %#ok<SAGROW>
        sigKept(end+1) = (S.PP.r(k) >= S.confC);           %#ok<SAGROW>  significant?

        % individual-spike snippets (each baselined to its own pre-spike window)
        Ei = nan(numel(sf), numel(tg));
        for j = 1:numel(sf), Ei(j,:) = interp1(tau, E(j,:) - mean(E(j,bsel)), tg, 'linear', NaN); end
        spkSnips{end+1} = Ei;                               %#ok<SAGROW>
    catch ME
        warning('  %s: %s', lab, ME.message);
    end
end
nROI = size(STAn,1);
assert(nROI>0, 'No active ROIs collected.');
fprintf('%d active ROIs (> %d events) | median t-to-peak=%.0f ms | median 10-90%% rise=%.0f ms\n', ...
        nROI, minEvents, 1000*median(ttp), 1000*median(rise1090,'omitnan'));
sigMask   = logical(sigKept(:));
nSig      = nnz(sigMask);
nSigTotal = nnz(S.PP.r >= S.confC);
fprintf('significant ROIs in STA set: %d / %d total (confC=%.3f)\n', nSig, nSigTotal, S.confC);
if nSig < nSigTotal
    warning('%d significant ROI(s) below the active gate (>%d events) -- not in the sig figure', ...
            nSigTotal - nSig, minEvents);
end

%% ===================== FIGURE =====================
green = [0.20 0.70 0.20];
mR = mean(STAr,1,'omitnan'); sR = std(STAr,0,1,'omitnan');   % baselined
mN = mean(STAn,1,'omitnan'); sN = std(STAn,0,1,'omitnan');   % peak-normalized
medRise = median(rise1090,'omitnan');

fig = figure('Color','w','Name','spike-triggered dF/F (active ROIs)','Position',[40 40 1000 820]);

% --- panel 1: baselined, all ROIs overlaid ---
ax1 = subplot(2,2,1); hold(ax1,'on');
plot(ax1, tg, STAr', '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.3);
xline(ax1, 0-ca_lag, 'r--','LineWidth',1); xlim(ax1,[-win_sec win_sec]);
xlabel(ax1,'time from spike (s)'); ylabel(ax1,'\DeltaF/F (baselined)');
title(ax1, sprintf('baselined: all ROIs overlaid (N=%d)', nROI)); box(ax1,'on');

% --- panel 2: baselined, mean +/- SD ---
ax2 = subplot(2,2,2); hold(ax2,'on');
fill(ax2, [tg fliplr(tg)], [mR+sR fliplr(mR-sR)], green, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(ax2, tg, mR, 'Color',green, 'LineWidth',2);
xline(ax2, 0-ca_lag, 'r--','LineWidth',1); xlim(ax2,[-win_sec win_sec]);
xlabel(ax2,'time from spike (s)'); ylabel(ax2,'\DeltaF/F (baselined)');
title(ax2, sprintf('baselined: mean \\pm SD  |  10-90%% rise=%.0f ms', 1000*medRise),'Interpreter','tex'); box(ax2,'on');

% --- panel 3: peak-normalized, all ROIs overlaid ---
ax3 = subplot(2,2,3); hold(ax3,'on');
plot(ax3, tg, STAn', '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.3);
xline(ax3, 0-ca_lag, 'r--','LineWidth',1); xlim(ax3,[-win_sec win_sec]); ylim(ax3,[-0.5 1.3]);
xlabel(ax3,'time from spike (s)'); ylabel(ax3,'dF/F (peak-normalized)');
title(ax3, 'peak-normalized: all ROIs overlaid'); box(ax3,'on');

% --- panel 4: peak-normalized, mean +/- SD ---
ax4 = subplot(2,2,4); hold(ax4,'on');
fill(ax4, [tg fliplr(tg)], [mN+sN fliplr(mN-sN)], [0.2 0.2 0.2], 'FaceAlpha',0.25, 'EdgeColor','none');
plot(ax4, tg, mN, 'k-','LineWidth',2);
xline(ax4, 0-ca_lag, 'r--','LineWidth',1); xlim(ax4,[-win_sec win_sec]); ylim(ax4,[-0.3 1.1]);
xlabel(ax4,'time from spike (s)'); ylabel(ax4,'dF/F (peak-normalized)');
title(ax4, 'peak-normalized: mean \pm SD','Interpreter','tex'); box(ax4,'on');

sgtitle(sprintf('Spike-triggered dF/F, active ROIs (>%d events), \\pm%g s window', ...
        minEvents, win_sec),'Interpreter','tex','FontWeight','bold');

%% ===================== FIGURE - SIGNIFICANT ROIs ONLY =====================
% Same 2x2 layout as figure 1, but only spikes from significant ROIs (PP.r >= confC).
if nSig > 0
    figSig = draw_sta_2x2(tg, STAr(sigMask,:), STAn(sigMask,:), rise1090(sigMask), win_sec, ca_lag, ...
        sprintf('Spike-triggered dF/F, SIGNIFICANT ROIs (N=%d), \\pm%g s window', nSig, win_sec));
else
    figSig = gobjects(0); warning('no significant ROIs in the STA set -- sig figure skipped');
end

%% ===================== SPIKE-TRIGGERED dF/F HEATMAP =====================
% rows = active ROIs (each = its baselined STA), sorted by peak amplitude; gray cmap.
[~, ordH] = sort(ampl, 'descend');
Hmat = STAr(ordH, :);
hv   = Hmat(isfinite(Hmat));
clH  = prctile(hv, [0.5 99.5]); if ~(clH(2)>clH(1)), clH = [min(hv) max(hv)]; end
figH = figure('Color','w','Name','spike-triggered dF/F heatmap','Position',[60 60 560 700]);
axH  = axes('Parent',figH); hold(axH,'on');
imagesc(axH, tg, 1:nROI, Hmat);
colormap(axH, flipud(gray(256))); caxis(axH, clH);
set(axH,'YDir','reverse'); axis(axH,'tight');
xline(axH, 0-ca_lag, 'r--','LineWidth',1); xlim(axH,[-win_sec win_sec]);
xlabel(axH,'time from spike (s)'); ylabel(axH,'active ROI (sorted by STA peak amplitude)');
title(axH, sprintf('spike-triggered dF/F (baselined)  |  %d active ROIs', nROI));
cb = colorbar(axH); cb.Label.String = '\DeltaF/F (baselined)';

%% ===================== INDIVIDUAL-SPIKE HEATMAP =====================
% rows = every individual spike (pooled across active ROIs), sorted by peak amplitude.
Eall  = vertcat(spkSnips{:});
pkAll = max(Eall(:, tg>=0), [], 2);
[~, ordS] = sort(pkAll, 'descend'); Eall = Eall(ordS, :);
hv2 = Eall(isfinite(Eall));
clS = prctile(hv2, [0.5 99.5]); if ~(clS(2)>clS(1)), clS = [min(hv2) max(hv2)]; end
figS = figure('Color','w','Name','spike-triggered dF/F heatmap (individual spikes)','Position',[640 60 560 700]);
axS  = axes('Parent',figS); hold(axS,'on');
imagesc(axS, tg, 1:size(Eall,1), Eall);
colormap(axS, flipud(gray(256))); caxis(axS, clS);
set(axS,'YDir','reverse'); axis(axS,'tight');
xline(axS, 0-ca_lag, 'r--','LineWidth',1); xlim(axS,[-win_sec win_sec]);
xlabel(axS,'time from spike (s)'); ylabel(axS,'individual spike (sorted by peak amplitude)');
title(axS, sprintf('spike-triggered dF/F (baselined)  |  %d individual spikes, %d ROIs', size(Eall,1), nROI));
cbS = colorbar(axS); cbS.Label.String = '\DeltaF/F (baselined)';

%% ===================== SAVE =====================
if doSave
    base = fullfile(outDir, 'spike_triggered_dFF_active');
    exportgraphics(fig,  [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig,  [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    exportgraphics(figH, [base '_heatmap.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(figH, [base '_heatmap.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    exportgraphics(figS, [base '_heatmap_individual.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(figS, [base '_heatmap_individual.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    if ~isempty(figSig) && isgraphics(figSig)
        exportgraphics(figSig, [base '_sigROI.png'], 'Resolution',200, 'BackgroundColor','white');
        exportgraphics(figSig, [base '_sigROI.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    end
    save([base '.mat'], 'tg','STAn','STAr','ttp','rise1090','ampl','nspkAll', ...
         'labKept','sigKept','win_sec','gridFps','baseWin','minEvents','ordH');
    fprintf('Saved %s.png/.pdf (+_heatmap +_sigROI)/.mat\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function fg = draw_sta_2x2(tg, STAr, STAn, rise1090, win_sec, ca_lag, ttl)
% 2x2 spike-triggered dF/F figure (same layout as the active-ROI figure 1).
green = [0.20 0.70 0.20]; n = size(STAr,1);
mR = mean(STAr,1,'omitnan'); sR = std(STAr,0,1,'omitnan');
mN = mean(STAn,1,'omitnan'); sN = std(STAn,0,1,'omitnan');
medRise = median(rise1090,'omitnan');
fg = figure('Color','w','Name',ttl,'Position',[120 60 1000 820]);
ax1 = subplot(2,2,1); hold(ax1,'on');
plot(ax1, tg, STAr', '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.3);
xline(ax1, 0-ca_lag, 'r--','LineWidth',1); xlim(ax1,[-win_sec win_sec]);
xlabel(ax1,'time from spike (s)'); ylabel(ax1,'\DeltaF/F (baselined)');
title(ax1, sprintf('baselined: ROIs overlaid (N=%d)', n)); box(ax1,'on');
ax2 = subplot(2,2,2); hold(ax2,'on');
fill(ax2, [tg fliplr(tg)], [mR+sR fliplr(mR-sR)], green, 'FaceAlpha',0.25, 'EdgeColor','none');
plot(ax2, tg, mR, 'Color',green, 'LineWidth',2);
xline(ax2, 0-ca_lag, 'r--','LineWidth',1); xlim(ax2,[-win_sec win_sec]);
xlabel(ax2,'time from spike (s)'); ylabel(ax2,'\DeltaF/F (baselined)');
title(ax2, sprintf('baselined: mean \\pm SD  |  10-90%% rise=%.0f ms', 1000*medRise),'Interpreter','tex'); box(ax2,'on');
ax3 = subplot(2,2,3); hold(ax3,'on');
plot(ax3, tg, STAn', '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.3);
xline(ax3, 0-ca_lag, 'r--','LineWidth',1); xlim(ax3,[-win_sec win_sec]); ylim(ax3,[-0.5 1.3]);
xlabel(ax3,'time from spike (s)'); ylabel(ax3,'dF/F (peak-normalized)');
title(ax3, 'peak-normalized: ROIs overlaid'); box(ax3,'on');
ax4 = subplot(2,2,4); hold(ax4,'on');
fill(ax4, [tg fliplr(tg)], [mN+sN fliplr(mN-sN)], [0.2 0.2 0.2], 'FaceAlpha',0.25, 'EdgeColor','none');
plot(ax4, tg, mN, 'k-','LineWidth',2);
xline(ax4, 0-ca_lag, 'r--','LineWidth',1); xlim(ax4,[-win_sec win_sec]); ylim(ax4,[-0.3 1.1]);
xlabel(ax4,'time from spike (s)'); ylabel(ax4,'dF/F (peak-normalized)');
title(ax4, 'peak-normalized: mean \pm SD','Interpreter','tex'); box(ax4,'on');
sgtitle(fg, ttl,'Interpreter','tex','FontWeight','bold');
end

function F = load_ca(folderPath, fallback_fps)
% dF/F + ca-spike trains for a folder (shared across its ROIs).
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
assert(~isempty(df), 'no *_ch1_dFF.mat');
assert(isfile(fullfile(folderPath,'ca_spike_data.mat')), 'no ca_spike_data');
F.fps = detect_session_fps(folderPath, fallback_fps);
D = load(fullfile(df(1).folder, df(1).name), 'dFF'); F.dFF = double(D.dFF);
CA = load(fullfile(folderPath,'ca_spike_data.mat'), 'roi_spikes'); F.roi_spikes = CA.roi_spikes;
end

function tc = rise_crossing(t, y, L, tpk)
% Time (interp) at which y reaches level L on the rising edge up to tpk.
tc = NaN; idx = find(t <= tpk); if isempty(idx), return; end
tt = t(idx); yy = y(idx);
below = find(yy < L);
if isempty(below), tc = tt(1); return; end
b = below(end);
if b >= numel(tt), tc = tt(end); return; end
y1 = yy(b); y2 = yy(b+1); t1 = tt(b); t2 = tt(b+1);
if y2 == y1, tc = t2; else, tc = t1 + (L - y1)*(t2 - t1)/(y2 - y1); end
end
