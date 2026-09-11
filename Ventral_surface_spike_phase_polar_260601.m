% Ventral_surface_spike_phase_polar_260601.m
% -----------------------------------------------------------------------
%  ONE polar plot for ONE ROI: the distribution of Ca-spike phases using
%  YOUR exact piecewise breath phase (inspiration ONSET = 0, breath PEAK =
%  pi, linear ramp in time between events) -- no cosine, no Hilbert.
%
%  Each spike is given the phase phi(t) at its frame; the polar histogram
%  shows how those spike phases are distributed around the breath cycle.
%
%  Overlaid for comparison: the coherence-method preferred phase (dashed
%  radial line) + its phase confidence interval (arc = +/-1.96*phistd),
%  reproduced from Ventral_surface_coherence_polar_260528.
%
%  Convention / colors: ONSET = 0 = RED ray, PEAK = pi = SKY-BLUE ray
%  (matches the HSV 0/pi key). Theta zero at right, counter-clockwise.
%
%  Inputs (in folderPath):
%     *DLC*breath_peak_data.mat       (breath waveform + PEAK idx)
%     *breath_insp_start_data.mat     (insp-start / foot idx)
%     ca_spike_data.mat               (roi_spikes(roi).spike_train)
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
folderPath = 'D:\Ventral_surface_summary\ChAT\0521\cell1\roi5_7x_x-1200y200z-30_3000f_23lp_00001';
roi        = 1;            % ROI index (== dFF column == spike-train index)

nDrop      = 30;           % breath frames to toss (match calcium)
fallback_fps = 30;
nBins      = 24;           % number of phase bins for the spike histogram
trace_xlim = [40 60];      % breath-trace window (s): [] = full trace

% ---- coherence-method overlay (match Ventral_surface_coherence_polar) ----
overlayCoherence = true;
TW_coh     = 4;            % multitaper TW
alpha_coh  = 0.001;        % significance level (confC + jackknife err)
ca_lag_sec = 0.015;        % GCaMP-rise lead compensation (shift spikes earlier)
f_breath_search = [0.2 4]; % Hz, breath PSD peak search band
fwhm_factor = 0.6;         % coherence band = fwhm_factor x FWHM
min_bw      = 0.05;        % Hz, min band width
fmin        = 0.05; fmax = 15;

% ---- colors ----
spikeCol = [0.30 0.30 0.30];      % spike histogram wedges
onsetCol = [0.90 0.10 0.10];      % RED  : onset (phase 0)
peakCol  = [0.35 0.75 1.00];      % SKY  : peak  (phase pi)
cohCol   = [0.00 0.00 0.00];      % coherence overlay (black, dashed)

doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
[~, recName] = fileparts(folderPath);

%% ===================== LOAD + ALIGN =====================
bp = dir(fullfile(folderPath,'*DLC*breath_peak_data.mat'));
ip = dir(fullfile(folderPath,'*DLC*breath_insp_start_data.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'*breath_insp_start_data.mat')); end
assert(~isempty(bp),'No *DLC*breath_peak_data.mat in %s', folderPath);
assert(~isempty(ip),'No *breath_insp_start_data.mat in %s', folderPath);
assert(isfile(fullfile(folderPath,'ca_spike_data.mat')),'No ca_spike_data.mat in %s', folderPath);

fps = detect_session_fps(folderPath, fallback_fps);
BP  = load(fullfile(bp(1).folder, bp(1).name));
IP  = load(fullfile(ip(1).folder, ip(1).name));
CA  = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes');
assert(roi<=numel(CA.roi_spikes),'ROI %d out of range (1..%d)',roi,numel(CA.roi_spikes));

% breath waveform (toss nDrop, detrend, demean) -- only for the breath band
bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);
nB = numel(BP.breath);

% PEAK + FOOT(onset) event trains, same toss
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;   % PEAK
ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;    % ONSET
ef(1:min(nDrop,numel(ef))) = [];

% spike train
spk = double(CA.roi_spikes(roi).spike_train(:));

% common length
T = min([numel(bw), numel(ev), numel(ef), numel(spk)]);
bw = bw(1:T); ev = ev(1:T); ef = ef(1:T); spk = spk(1:T);

peakAll = find(ev>0);
footAll = find(ef>0);
fprintf('%s ROI%d | T=%d @%.3g Hz | %d peaks | %d onsets | %d spikes\n', ...
        recName, roi, T, fps, numel(peakAll), numel(footAll), sum(spk>0));

%% ===================== EXACT BREATH PHASE -> SPIKE PHASES =====================
phi = piecewise_phase_local(peakAll, footAll, T);   % 0 at onset, pi at peak, NaN outside

lagFrames = round(ca_lag_sec*fps);
spk_fr  = find(spk>0);
spk_src = spk_fr - lagFrames; spk_src = spk_src(spk_src>=1 & spk_src<=T);
spk_phi = phi(spk_src); spk_phi = spk_phi(~isnan(spk_phi));
spk_phi = mod(spk_phi, 2*pi);   % WRAP to one cycle (phi ramps cumulatively 0,pi,2pi,...)

edges = linspace(0, 2*pi, nBins+1);
ctrs  = (edges(1:end-1)+edges(2:end))/2;
cnt   = histcounts(spk_phi, edges);

% number of breath cycles spanned by the phase (phi ramps 0,pi,2pi,...)
phi_valid = phi(~isnan(phi));
nCycles   = (max(phi_valid) - min(phi_valid)) / (2*pi);

% normalize per breath cycle: bin value = avg spikes/cycle in that phase (%)
% e.g. 2 spikes/cycle all in one bin -> 200%
pct = 100 * cnt / max(nCycles, eps);
fprintf('nCycles=%.1f | total spikes binned=%d | mean spikes/cycle=%.2f\n', ...
        nCycles, sum(cnt), sum(cnt)/max(nCycles,eps));

% circular mean of the spike phases + 95% CI half-width (for the phase error bar)
if ~isempty(spk_phi)
    mu_spk = mod(angle(mean(exp(1i*spk_phi))), 2*pi);
    ci_spk = circ_conf_local(spk_phi);
else
    mu_spk = NaN; ci_spk = NaN;
end

%% ===================== COHERENCE METHOD (for the overlay) =====================
th_coh = NaN; r_coh = NaN; dphi = NaN; confC = NaN; f_pk = NaN; rlo = NaN; rhi = NaN;
if overlayCoherence && numel(spk_phi) >= 2
    % breath PSD peak + FWHM band
    pB.Fs=fps; pB.tapers=[TW_coh,2*TW_coh-1]; pB.pad=0; pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
    [Sb,fb] = mtspectrumc(bw, pB); Sb=Sb(:); fb=fb(:);
    mm=fb>=f_breath_search(1)&fb<=f_breath_search(2);
    [~,rl]=max(Sb(mm)); ipk=find(mm,1)+rl-1; f_pk=fb(ipk);
    hh=Sb(ipk)/2; lo=ipk; while lo>1&&Sb(lo)>hh, lo=lo-1; end
    hi=ipk;       while hi<numel(fb)&&Sb(hi)>hh, hi=hi+1; end
    f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
    bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
    band=[max(f_pk-bwd/2,fmin), min(f_pk+bwd/2,fmax)];

    ref = cos(phi); ref(isnan(ref)) = 0; ref = ref - mean(ref);
    pcoh.Fs=fps; pcoh.tapers=[TW_coh,2*TW_coh-1]; pcoh.pad=0; pcoh.fpass=band; pcoh.err=[2,alpha_coh];
    [~, Cmag, cphi, ~,~,~, fC, confC, phistd, Cerr] = coherencyc(ref, spk-mean(spk), pcoh);
    fC=fC(:); mbc = fC>=band(1)&fC<=band(2); if ~any(mbc), mbc=true(size(fC)); end
    r_coh  = mean(Cmag(mbc));
    th_coh = wrapToPi(angle(mean(exp(1i*(-cphi(mbc))))) - 2*pi*f_pk*ca_lag_sec);
    dphi   = 1.96*mean(phistd(mbc));
    rlo    = max(0, mean(Cerr(1,mbc)));
    rhi    = min(1, mean(Cerr(2,mbc)));
end

%% ===================== FIGURE: breath trace (top) + polar + linear =====================
fig = figure('Color','w','Name',sprintf('%s ROI%d spike-phase',recName,roi), ...
             'Units','centimeters','Position',[2 2 38 20]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
sgtitle(fig, sprintf('%s  ROI%d  --  spike-phase distribution (n=%d spikes, %d bins)   |   RED=onset(0)  SKY=peak(\\pi)', ...
        recName, roi, numel(spk_phi), nBins), 'Interpreter','tex','FontWeight','bold');
rmax = max([pct, 1]);

% ---- (0) BREATH TRACE + events + wrapped phase sawtooth (top, full width) ----
t = (0:T-1)'/fps;
phiW = mod(phi, 2*pi);                 % wrapped phase 0..2pi for display
if isempty(trace_xlim), w = [t(1) t(end)]; else, w = [min(trace_xlim) max(trace_xlim)]; end
mw = t>=w(1) & t<=w(2);
bax = nexttile(tl,1,[1 3]); hold(bax,'on');
% phase (shows 0->pi inspiration vs pi->2pi expiration timing)
plot(bax, t(mw), phiW(mw), 'k-', 'LineWidth',0.8);
set(bax,'YTick',[0 pi 2*pi],'YTickLabel',{'0','\pi','2\pi'});
ylabel(bax,'phase'); ylim(bax,[0 2*pi]);
% onset (red) + peak (sky) event markers
on_t = footAll/fps; pk_t = peakAll/fps;
on_t = on_t(on_t>=w(1) & on_t<=w(2)); pk_t = pk_t(pk_t>=w(1) & pk_t<=w(2));
yl = ylim(bax);
for x = on_t(:)', plot(bax,[x x],yl,'-','Color',onsetCol,'LineWidth',0.8); end
for x = pk_t(:)', plot(bax,[x x],yl,'-','Color',peakCol,'LineWidth',0.8); end
% spikes as black ticks along the bottom
spk_t = spk_fr/fps; spk_t = spk_t(spk_t>=w(1) & spk_t<=w(2));
plot(bax, spk_t, (yl(1)+0.03*diff(yl))*ones(size(spk_t)), '|','Color','k','MarkerSize',6,'LineWidth',1);
xlim(bax,w); xlabel(bax,'Time (s)'); box(bax,'off');
title(bax,'breath (gray) + onset (red) + peak (sky) + phase sawtooth (purple) + spikes (k ticks)');

% ---- (1) EVENT-PHASE polar: spike-phase distribution ----
pax = polaraxes(tl); pax.Layout.Tile = 5; hold(pax,'on');
polarhistogram(pax,'BinEdges',edges,'BinCounts',pct, ...
        'FaceColor',spikeCol,'FaceAlpha',0.75,'EdgeColor',[0.2 0.2 0.2]);
pax.RLim = [0 rmax];
pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise';
pax.RAxisLocation=180; pax.FontSize=9; thetaticks(pax,0:45:315);
title(pax, sprintf('event-phase polar  |  spikes/cycle (%%)  \\mu=%.0f\\circ', rad2deg(mu_spk)),'Interpreter','tex');

% ---- (2) COHERENCE polar: preferred phase + magnitude + jackknife CI ----
cax = polaraxes(tl); cax.Layout.Tile = 4; hold(cax,'on');
if overlayCoherence && isfinite(th_coh)
    thc = linspace(0,2*pi,361);
    if isfinite(confC), polarplot(cax, thc, confC*ones(size(thc)),'k--','LineWidth',1); end
    polarplot(cax,[th_coh th_coh],[rlo rhi],'-','Color','k','LineWidth',1.8);     % magnitude CI
    if isfinite(dphi) && dphi>0
        a = linspace(th_coh-dphi, th_coh+dphi, 60);
        polarplot(cax, a, r_coh*ones(size(a)),'-','Color','k','LineWidth',1.8);    % phase CI arc
    end
    polarplot(cax, th_coh, r_coh,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
end
cax.RLim=[0 1]; cax.ThetaZeroLocation='right'; cax.ThetaDir='counterclockwise';
cax.RAxisLocation=180; cax.FontSize=9; thetaticks(cax,0:45:315);
title(cax, sprintf('coherence polar  |  %.0f\\circ  r=%.2f  conf=%.2f', ...
      rad2deg(mod(th_coh,2*pi)), r_coh, confC),'Interpreter','tex');

% ---- (3) LINEAR histogram, tiled 0..4pi (single cycle duplicated) ----
lax = nexttile(tl,6); hold(lax,'on');
ctrs_tile = [ctrs, ctrs+2*pi];
bar(lax, ctrs_tile, [pct pct], 1, 'FaceColor',spikeCol,'FaceAlpha',0.75,'EdgeColor','none');
xlim(lax,[0 4*pi]); ylim(lax,[0 rmax*1.05]);
set(lax,'XTick',[0 pi 2*pi 3*pi 4*pi],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
xlabel(lax,'breath phase (onset=0, peak=\pi)'); ylabel(lax,'spikes/cycle (%)');
box(lax,'on');
title(lax, sprintf('linear (0..4\\pi, duplicated)  |  %.1f cycles, %.2f spikes/cycle', ...
      nCycles, sum(cnt)/max(nCycles,eps)),'Interpreter','tex');

%% ===================== SAVE =====================
if doSave
    base = fullfile(folderPath, sprintf('spike_phase_polar_ROI%02d', roi));
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('Saved %s.png/.pdf\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function t = circ_conf_local(alpha)
% 95% CI half-width (rad) for the circular mean direction (Fisher 1993,
% as in CircStat circ_confmean). NaN when undefined (too dispersed / few).
alpha = alpha(:);
n = numel(alpha);
if n < 2, t = NaN; return; end
r  = abs(mean(exp(1i*alpha)));      % mean resultant length
R  = n*r;
c2 = 3.8414588;                     % chi2inv(0.95,1)
if r < 0.9 && r > sqrt(c2/2/n)
    tval = sqrt((2*n*(2*R^2 - n*c2))/(4*n - c2));
elseif r >= 0.9
    tval = sqrt(n^2 - (n^2 - R^2)*exp(c2/n));
else
    t = NaN; return;
end
t = acos(tval/R);
if ~isreal(t) || ~isfinite(t), t = NaN; end
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase reference: FEET at 0/2pi/..., PEAKS at pi/3pi/...
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
[events, ord] = sort(events);
types = types(ord);
keep = true(size(events));
for i = 2:numel(events)
    if types(i) == types(i-1), keep(i) = false; end
end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases  = nan(size(events));
phi_cur = types(1) * pi;      % type=1 (peak) -> pi; type=0 (foot) -> 0
for i = 1:numel(events)
    phases(i) = phi_cur; phi_cur = phi_cur + pi;
end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
