% polar_recheck_coherence_260901.m
% -----------------------------------------------------------------------
%  The re-check polar plot again, with COHERENCE on the radial axis instead of
%  Rayleigh log Z. One dot per cell, nothing else.
%
%      angle  = weighted circular MEAN of the event phases -- THE SAME ANGLE AS
%               polar_recheck_260901. These two figures differ in the RADIUS
%               ONLY: a cell sits at the same bearing on both, and moving
%               between them shows how the two measures of strength -- Rayleigh
%               log Z and coherence -- rank the same set of directions. The
%               coherence phase is computed too, and written to the CSV, but it
%               is not what the dot is drawn at.
%      r      = band-averaged coherence magnitude |C|
%      colour = GENOTYPE  (vagotomised is not distinguished)
%      filled = r >= confC (Chronux confidence level, alpha = 0.01)
%      open   = below it;  IO is always open grey
%
%  THE COHERENCE CONVENTION IS THE PROJECT'S, unchanged -- copied step for step
%  out of Ventral_surface_coherence_polar_svd_260729.m so an r here is the same
%  number as an r there:
%    * reference signal is cos(phi), phi = the piecewise landmark phase
%      (onset = 0, peak = pi), NaN -> 0, mean removed -- NOT the raw breath
%      waveform. The waveform is used only to find the band.
%    * TW = 4, tapers [TW, 2*TW-1] = [4 7], pad 0, err [2, 0.001]
%    * band = breath-PSD peak in [0.2 4] Hz, width fwhm_factor(0.6) x FWHM,
%      floor min_bw 0.05 Hz, clamped back into [0.2 4]
%    * spikes are LEAD-SHIFTED by ca_lag_sec = 0.1 s (3 frames at 30 Hz) BEFORE
%      coherencyc -- a time shift, not a post-hoc rotation of theta. It leaves
%      |C| and confC untouched and moves only the phase.
%    * r  = mean(Cxy) in band,  th = angle(mean(exp(-i*phi_C))) in band
%    * confC = sqrt(1 - alpha^(1/(K-1))), K = 2*TW-1 = 7  ->  0.732 at 0.01
%  The Vglut2/1124 one-frame trigger fix is already applied in the bundle, by
%  the same rule the archive uses.
%
%  THE PHASE IS THE LINEAR TWO-LANDMARK PHASE, as everywhere else: inspiration
%  onset = 0, peak = pi, linear in time between landmarks; event phases pooled
%  with each recording's own occupancy weights (36 bins, w = 1/frames-in-bin) and
%  averaged with Kish n_eff. Identical code to polar_recheck_260901, and NO
%  calcium lag -- the 0.1 s lead-shift applies to the coherence estimate only.
%
%  ORIENTATION matches polar_recheck_260901: 0 to the RIGHT, counterclockwise,
%  so onset is at 3 o'clock and peak at 9 o'clock. (The archive coherence panel
%  drew theta-zero at the TOP going clockwise; only the drawing changed.)
%
%  A CELL, NOT AN ROI. Coherence is computed per recording -- it has to be, the
%  band is set by that recording's own breath rate -- and a cell imaged more than
%  once is pooled as the frame-weighted mean of r and the frame-weighted circular
%  mean of th. Concatenating recordings before coherencyc would put a
%  discontinuity in the middle of the segment.
%
%  Population = the re-check curation: cells tossed whole are dropped, and within
%  a kept cell only the kept recordings contribute, with that GUI's curated
%  events rebuilt into the binary train calcium_spike_gui stores.
%
%  Runqi Zhang / 2026-09-01
% -----------------------------------------------------------------------
clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
bundleDir = fullfile(sumRoot,'per-cell-summary_active_260812','spike_recheck_260901');
outDir    = fullfile(sumRoot,'polar_recheck_260901');

TW              = 4;         % multitaper TW for coherence
alpha_sig       = 0.01;      % primary: the filled/open cut, black dashed circle
alpha_sig2      = 0.001;     % printed in the log only; not drawn
ca_lag_sec      = 0.1;       % spikes lead-shifted 3 frames @30 Hz before coherencyc
minSpikes       = 2;
nBins           = 36;        % phase bins for the occupancy weights (as in the logZ figure)
f_breath_search = [0.2 4];   % Hz, search band for the breath PSD peak
fwhm_factor     = 0.6;
min_bw          = 0.05;      % Hz
fmin            = 0.05;      % Hz, PSD bounds
fmax            = 15;
% =========================================================

K_tap  = 2*TW - 1;
confC  = sqrt(1 - alpha_sig ^(1/(K_tap - 1)));    % 0.732 at alpha = 0.01
confC2 = sqrt(1 - alpha_sig2^(1/(K_tap - 1)));

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

%% ===================== LOAD THE CURATED POPULATION =====================
assert(isfolder(bundleDir), 'No bundle at %s', bundleDir);
CU = ca_recheck_load_curation_260901(bundleDir);
fprintf('curation: %d cells, %d dropped observations, %d cells tossed whole\n', ...
        CU.nCells, CU.nTossed, CU.nTossedCells);

d = dir(fullfile(bundleDir,'cells','*.mat'));
GC = genotype_colors_260817();
GRPS = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};

R = struct('cell',{},'stem',{},'group',{},'th',{},'thCoh',{},'r',{},'nEv',{}, ...
           'nRecKept',{},'f_pk',{},'band',{});
nTossCell = 0; nTossRec = 0; nThin = 0;
tic;
for k = 1:numel(d)
    L = load(fullfile(d(k).folder, d(k).name),'cellInfo','RECc');
    C = L.cellInfo;  REC = L.RECc;

    rAll = []; thAll = []; wAll = []; fAll = []; bAll = []; nEvC = 0; nKept = 0;
    aPool = []; wPool = [];          % event phases + occupancy weights, for the ANGLE
    tossedWhole = false;
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && g.cell_toss, tossedWhole = true; break; end
        if g.found && g.toss, nTossRec = nTossRec + 1; continue; end
        ev = r.spike_idx(:);                      % archive fallback
        if g.found, ev = g.spike_idx(:); end      % curated events win
        nKept = nKept + 1;

        T   = r.T;  fps = r.fps;
        ev  = ev(ev >= 1 & ev <= T);
        nEvC = nEvC + numel(ev);

        % --- the ANGLE: event phases with this recording's occupancy weights ---
        % THIS BLOCK RUNS FIRST AND UNGATED, line for line as in
        % polar_recheck_260901 -- same landmarks, same 36 bins, same guard, raw
        % (unshifted) events. Any extra condition here, such as the per-recording
        % minSpikes the coherence estimate needs, would drop a recording from the
        % angle that the log Z figure keeps, and the same cell would then sit at
        % two different bearings on two figures meant to be compared. That cost
        % up to 3.8 deg before it was moved above the gate.
        phi   = piecewise_phase_local(r.peak, r.foot, T);
        phiW  = mod(phi, 2*pi);
        valid = ~isnan(phiW);
        if nnz(valid) >= 10
            edgesP = linspace(0, 2*pi, nBins+1);
            bidx = nan(T,1);
            bidx(valid) = min(discretize(phiW(valid), edgesP), nBins);
            occ  = accumarray(bidx(valid), 1, [nBins 1]);
            evp  = ev(~isnan(bidx(ev)));
            if ~isempty(evp)
                aPool = [aPool; phiW(evp)];                 %#ok<AGROW>
                wPool = [wPool; 1./max(occ(bidx(evp)),1)];  %#ok<AGROW>
            end
        end

        % --- from here on: the COHERENCE estimate, which has its own gates ---
        if numel(ev) < minSpikes, continue; end
        pk = r.peak(r.peak>=1 & r.peak<=T);
        ft = r.foot(r.foot>=1 & r.foot<=T);
        if numel(pk) < 2 || numel(ft) < 2, continue; end
        ref = cos(phi); ref(isnan(ref)) = 0; ref = ref - mean(ref);

        % --- breath PSD -> coherence band ---
        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
        pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
        [Sb,fb] = mtspectrumc(r.breath(1:T), pB); Sb=Sb(:); fb=fb(:);
        m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        if ~any(m), continue; end
        [~,rl]=max(Sb(m)); ipk=find(m,1)+rl-1; f_pk=fb(ipk);
        h=Sb(ipk)/2; lo=ipk; while lo>1 && Sb(lo)>h, lo=lo-1; end
        hi=ipk;            while hi<numel(fb) && Sb(hi)>h, hi=hi+1; end
        f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
        band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

        pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0;
        pc.fpass=band; pc.err=[2,alpha_sig];

        % --- binary spike train, lead-shifted, then coherencyc ---
        st = zeros(T,1); st(ev) = 1;
        lag = round(ca_lag_sec*fps);
        stL = [st(1+lag:end); zeros(lag,1)];
        [~, Cxy, phiC, ~,~,~, f] = coherencyc(ref, stL - mean(stL), pc);
        f  = f(:);
        mb = f>=band(1) & f<=band(2);
        if ~any(mb), mb = true(size(f)); end

        rAll(end+1,1)  = mean(Cxy(mb));                        %#ok<SAGROW>
        thAll(end+1,1) = angle(mean(exp(1i*(-phiC(mb)))));     %#ok<SAGROW>
        wAll(end+1,1)  = T;                                    %#ok<SAGROW>
        fAll(end+1,1)  = f_pk;                                 %#ok<SAGROW>
        bAll(end+1,:)  = band;                                 %#ok<SAGROW>
    end
    if tossedWhole, nTossCell = nTossCell + 1; continue; end
    if isempty(rAll) || numel(aPool) < minSpikes, nThin = nThin + 1; continue; end

    % RADIUS: frame-weighted mean of the per-recording coherence.
    % ANGLE:  weighted circular mean of the pooled event phases -- pooled over
    %         recordings exactly as the log Z figure pools them, so the two
    %         figures put this cell at the same bearing.
    w  = wAll / sum(wAll);
    rP = sum(w .* rAll);
    cohP = angle(sum(w .* exp(1i*thAll)));     % coherence phase, CSV only
    tP = wresultant_local(aPool, wPool);
    R(end+1) = struct('cell',C.cell, 'stem',C.stem, 'group',C.group, ...
        'th',tP, 'thCoh',cohP, 'r',rP, 'nEv',nEvC, 'nRecKept',nKept, ...
        'f_pk',sum(w.*fAll), 'band',sum(w.*bAll,1)); %#ok<SAGROW>
    if mod(k,50)==0, fprintf('  %d/%d (%.0f s)\n', k, numel(d), toc); end
end

n = numel(R);
fprintf('%d cells plotted  (%d tossed whole, %d recordings tossed, %d with < %d events)\n', ...
        n, nTossCell, nTossRec, nThin, minSpikes);
assert(n > 0, 'no cells survived the curation');

thPlot = [R.th];  rPlot = [R.r];
gname  = string({R.group});
isSig  = rPlot >= confC;

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 20], ...
             'Name','polar: re-checked cells, coherence');
set(fig,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);
ax = polaraxes(fig,'Position',[0.05 0.10 0.66 0.72]);
hold(ax,'on');
ax.ThetaZeroLocation = 'right';
ax.ThetaDir          = 'counterclockwise';

% ONE threshold circle, at the alpha the filled/open cut uses.
tt = linspace(0,2*pi,361);
polarplot(ax, tt, confC*ones(size(tt)), '--','Color',[0.15 0.15 0.15],'LineWidth',0.9);
text(ax, deg2rad(-18), confC, sprintf('  \\alpha=%.2g', alpha_sig), ...
     'FontSize',7,'Color',[0.15 0.15 0.15]);

present = GRPS(ismember(GRPS, cellstr(unique(gname))));
h = gobjects(numel(present),1);
for gi = 1:numel(present)
    g   = present{gi};
    col = GC.(g);
    m   = find(gname == string(g) & isfinite(thPlot));
    if isempty(m), continue; end
    if strcmp(g,'IO')
        h(gi) = polarplot(ax, thPlot(m), rPlot(m), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',[0.50 0.50 0.50], ...
            'LineWidth',0.7,'MarkerSize',6.5,'LineStyle','none');
        continue
    end
    ms = m(isSig(m));  mn = m(~isSig(m));
    if ~isempty(mn)
        polarplot(ax, thPlot(mn), rPlot(mn), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',col, ...
            'LineWidth',0.9,'MarkerSize',6.5,'LineStyle','none');
    end
    if ~isempty(ms)
        h(gi) = polarplot(ax, thPlot(ms), rPlot(ms), 'o', ...
            'MarkerFaceColor',col,'MarkerEdgeColor','w', ...
            'LineWidth',0.5,'MarkerSize',6.5,'LineStyle','none');
    else
        h(gi) = polarplot(ax, thPlot(mn(1)), rPlot(mn(1)), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',col, ...
            'LineWidth',0.9,'MarkerSize',6.5,'LineStyle','none');
    end
end

rlim(ax,[0 1]);
ax.ThetaTick      = 0:30:330;
ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:30:330, 'uni',0);
ax.RAxisLocation  = 180;
ax.GridAlpha      = 0.15;

annotation(fig,'textbox',[0.02 0.845 0.70 0.045],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',10, ...
    'String',['r = coherence |C| in the breath band      ' ...
              'dot = circular MEAN of event phase (same angle as the log Z figure)']);
annotation(fig,'textbox',[0.02 0.90 0.96 0.09],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','tex','FontSize',10, ...
    'String',{sprintf(['N = %d cells   |   multitaper TW = %g, tapers [%g %g]   |   ' ...
                       'insp onset = 0 (right), insp peak = \\pi (left), counterclockwise'], ...
                      n, TW, TW, 2*TW-1), ...
              sprintf(['filled = r \\geq confC(\\alpha=%.3g) = %.3f : %d      ' ...
                       'open = below : %d      IO always open grey'], ...
                      alpha_sig, confC, nnz(isSig), nnz(~isSig))});
annotation(fig,'textbox',[0.02 0.005 0.96 0.045],'EdgeColor',[0.7 0.7 0.7], ...
    'BackgroundColor',[0.96 0.96 0.96],'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontWeight','bold','FontSize',9, ...
    'String',sprintf(['re-check curation %s : %d cells tossed whole, %d recordings tossed   |   ' ...
                      'ref = cos(landmark phase), spikes lead-shifted %g ms before coherencyc'], ...
                      datestr(now,'yyyy-mm-dd'), nTossCell, nTossRec, ca_lag_sec*1000)); %#ok<TNOW1,DATST>

ok = isgraphics(h);
lg = present(ok);
for q = 1:numel(lg)
    mm = gname == string(lg{q});
    lg{q} = sprintf('%s (n=%d, sig %d)', lg{q}, nnz(mm), nnz(mm & isSig));
end
lgd = legend(ax, h(ok), lg, 'Box','off','FontSize',9);
lgd.Units = 'normalized';  lgd.Position = [0.735 0.55 0.24 0.22];

axT = axes(fig,'Position',[0.735 0.10 0.25 0.40]); axis(axT,'off');
Lt = strings(0,1);
Lt(end+1) = sprintf('%-8s %4s %8s %8s %5s', 'genotype','n','meanDeg','med |C|','sig');
Lt(end+1) = string(repmat('-',1,38));
for gi = 1:numel(present)
    m = find(gname == string(present{gi}) & isfinite(thPlot));
    if isempty(m), continue; end
    mu = angle(mean(exp(1i*thPlot(m))));
    Lt(end+1) = sprintf('%-8s %4d %8.0f %8.3f %5d', present{gi}, numel(m), ...
                        mod(rad2deg(mu),360), median(rPlot(m)), nnz(isSig(m))); %#ok<SAGROW>
end
text(axT, 0, 1, Lt, 'FontName','Consolas','FontSize',8, ...
     'VerticalAlignment','top','Interpreter','none');

%% ===================== SAVE =====================
stem = fullfile(outDir,'polar_recheck_coherence');
exportgraphics(fig, [stem '.png'], 'Resolution',300, 'BackgroundColor','white');
exportgraphics(fig, [stem '.pdf'], 'ContentType','vector','BackgroundColor','white');

T = table([R.cell]', string({R.stem})', string({R.group})', ...
          mod(rad2deg([R.th]),360)', mod(rad2deg([R.thCoh]),360)', ...
          [R.r]', [R.f_pk]', ...
          arrayfun(@(x) x.band(1), R)', arrayfun(@(x) x.band(2), R)', ...
          [R.nEv]', [R.nRecKept]', isSig', ...
    'VariableNames',{'cell_idx','stem','group','mean_phase_deg','coh_phase_deg', ...
                     'coh_r','f_peak_hz','band_lo_hz','band_hi_hz','n_events', ...
                     'n_rec_kept','sig'});
writetable(T, fullfile(outDir,'polar_recheck_coherence_percell.csv'));

fprintf('\nconfC(alpha=%.3g) = %.3f    confC(alpha=%.2g) = %.3f\n', ...
        alpha_sig, confC, alpha_sig2, confC2);
fprintf('%-8s %4s %8s %8s %5s\n','genotype','n','meanDeg','med |C|','sig');
for gi = 1:numel(present)
    m = find(gname == string(present{gi}));
    mu = angle(mean(exp(1i*thPlot(m))));
    fprintf('%-8s %4d %8.0f %8.3f %5d\n', present{gi}, numel(m), ...
            mod(rad2deg(mu),360), median(rPlot(m)), nnz(isSig(m)));
end
fprintf('\nwrote %s.png / .pdf and polar_recheck_coherence_percell.csv\n', stem);

%% ===================== LOCAL FUNCTIONS =====================
function [th, Rbar, nEff] = wresultant_local(a, w)
% Weighted circular resultant with Kish's effective sample size. Verbatim from
% polar_selected_260816.m, and the same call polar_recheck_260901 makes, so the
% angle on the two figures is the same number.
S1 = sum(w);  S2 = sum(w.^2);
if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
nEff = S1^2 / max(S2, eps);
v    = sum(w(:) .* exp(1i*a(:))) / S1;
th   = angle(v);
Rbar = min(abs(v), 1);
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
%PIECEWISE_PHASE_LOCAL  Insp onset = 0, peak = pi, next onset = 2pi, linear in
%  TIME between consecutive landmarks. Verbatim from temporal_phase_cell_fig_260812.
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
