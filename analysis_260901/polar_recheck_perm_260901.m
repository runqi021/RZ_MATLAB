% polar_recheck_perm_260901.m
% -----------------------------------------------------------------------
%  Six-panel polar: three radii x two angle definitions, significance by the
%  PEAK-TRIGGERED PSTH CIRCULAR-SHIFT PERMUTATION -- the same p the per-cell
%  summary figure prints (stats.pPeak of temporal_phase_cell_fig_260812).
%
%     columns   r = Rayleigh log Z  |  coherence |C|  |  -log10(permutation p)
%     rows      angle = circular MEAN (top)  |  circular MEDIAN (bottom)
%
%  Same 321 re-checked cells, same linear two-landmark phase (insp onset = 0 at
%  the right, insp peak = pi at the left, counterclockwise). A cell sits at the
%  same bearing across a row; only the radius changes. ONE significance test
%  decides filled vs open on all six panels.
%
%  THE TEST, copied step for step from the histogram block of
%  temporal_phase_cell_fig_260812.m:
%    * trigger = inspiration PEAK; triggers kept only where a full +/-1 IBI
%      window fits inside the recording
%    * histogram of event times relative to every trigger, bins of
%      histBinFrames = 2 frames, window +/- winSecH = 1 pooled IBI, expressed as
%      % per trigger; recordings pooled by summing counts and triggers
%    * null: each kept recording's event train is CIRCULARLY SHIFTED by its own
%      random offset, minShift = max(window+1, 3 breath cycles, 1 s), and the
%      pooled histogram is rebuilt; 2000 shuffles
%    * statistic = max |observed - null mean| over the inner +/- IBI/2, and
%      p = (1 + #{null >= observed}) / (1 + nShuffle)
%
%  THE SPAN is the weighted circular IQR of the cell's own event phases, drawn
%  only on cells the permutation calls significant. Circularity is handled by
%  working in RESIDUALS: each event's offset from the centre, u = angle(exp(i*(a
%  - centre))), wrapped to (-pi, pi], which puts the events on a line; the 25th
%  and 75th weighted percentiles of u are then ordinary percentiles, and the arc
%  is drawn from centre+q25 to centre+q75. The residuals are taken about THE
%  CENTRE THAT ROW DRAWS -- mean on the top row, median on the bottom -- so the
%  dot and its arc always describe the same centre. (The 260808 archive figure
%  moved the dot to the mean but kept a median-centred arc, and they disagreed by
%  >30 deg for 61 of 445 cells.)
%
%  NOT A JACKKNIFE. The old coherence polar drew a jackknife phase bar, which is
%  the uncertainty of the coherence ESTIMATE across leave-one-taper-out replicates
%  -- how well that one number is pinned down. This arc is the SPREAD OF THE
%  EVENTS themselves. The two answer different questions and are not comparable.
%
%  Arcs are restricted to significant cells because the linearisation only holds
%  while residuals stay well inside +/-pi: for a near-uniform cell the quartiles
%  approach +/-pi/2 and the arc would imply a precision the data do not have.
%
%  COHERENCE IS NOT PERMUTED -- it is only a radius here. Its analytic confC is
%  in the CSV for reference but nothing on this figure is marked with it.
%
%  p CANNOT GO BELOW 1/(nShuffle+1). At 2000 shuffles that is 5.0e-4, so the
%  -log10(p) axis ceilings at 3.30 and the alpha = 0.001 cut IS expressible --
%  but only just: p = (1+k)/2001, so a cell needs k <= 1, i.e. at most ONE of
%  2000 shuffles beating it. The test therefore has exactly two passing outcomes
%  (0 or 1 shuffles), and a pile-up on the rim still means "at the resolution
%  floor", not "infinitely significant".
%
%  THE SHUFFLE IS EXACT AND FAST. Rebuilding the trigger histogram 1200 times per
%  recording the direct way is ~9 billion event-trigger differences. Instead the
%  circular correlogram c[m] = #{pairs with mod(event - trigger, T) = m} is built
%  ONCE, and a circular shift of the train by sh is then just a rotation of c:
%  the count at lag d becomes c[mod(d - sh, T)]. Same numbers, ~60 operations per
%  shuffle instead of ~20,000. The script CHECKS this against the direct
%  histogram at shift 0 before trusting it.
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

nShuffle    = 2000;      % raised from the project's 1200 so fewer cells sit on the
                         % p floor; floor = 1/(nShuffle+1) = 5.0e-4, -log10 = 3.30
shiftMinCyc = 3;         % minimum circular shift, in breath cycles
histBinFr   = 2;         % PSTH bin width, in frames
pCrit       = 0.001;     % filled/open cut, i.e. -log10(p) >= 3
rngSeed     = 260901;

TW              = 4;     % coherence only (radius, not tested here)
alphaCoh        = 0.01;  % confC written to the CSV for reference
ca_lag_sec      = 0.1;
minSpikes       = 2;
nBins           = 36;    % phase bins for the occupancy weights
f_breath_search = [0.2 4];
fwhm_factor     = 0.6;
min_bw          = 0.05;
fmin            = 0.05;
fmax            = 15;
% =========================================================

K_tap = 2*TW - 1;
confC = sqrt(1 - alphaCoh^(1/(K_tap-1)));

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
rng(rngSeed);

%% ===================== LOAD =====================
CU = ca_recheck_load_curation_260901(bundleDir);
fprintf('curation: %d cells, %d dropped observations, %d cells tossed whole\n', ...
        CU.nCells, CU.nTossed, CU.nTossedCells);
d  = dir(fullfile(bundleDir,'cells','*.mat'));
GC = genotype_colors_260817();
GRPS = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};

R = struct('cell',{},'stem',{},'group',{},'thMean',{},'thMed',{}, ...
           'qMean',{},'qMed',{},'logZ',{},'cohR',{},'p',{},'nEv',{},'nRecKept',{});
nTossCell = 0; nTossRec = 0; nThin = 0; nNoNull = 0;
checkedShortcut = false;
tAll = tic;

for k = 1:numel(d)
    L = load(fullfile(d(k).folder, d(k).name),'cellInfo','RECc');
    Ci = L.cellInfo;  REC = L.RECc;

    % ---- gather the kept recordings ----
    Q = struct('T',{},'fps',{},'ev',{},'phiW',{},'bidx',{},'occ',{},'peak',{}, ...
               'foot',{},'f_pk',{},'stL',{},'ref',{},'pc',{},'mb',{},'w',{});
    tossedWhole = false;  nKept = 0;  nEvC = 0;
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && g.cell_toss, tossedWhole = true; break; end
        if g.found && g.toss, nTossRec = nTossRec + 1; continue; end
        ev = r.spike_idx(:);
        if g.found, ev = g.spike_idx(:); end
        nKept = nKept + 1;

        T = r.T; fps = r.fps;
        ev = ev(ev >= 1 & ev <= T);
        nEvC = nEvC + numel(ev);

        phi   = piecewise_phase_local(r.peak, r.foot, T);
        phiW  = mod(phi, 2*pi);
        valid = ~isnan(phiW);
        if nnz(valid) < 10 || isempty(ev), continue; end
        edgesP = linspace(0, 2*pi, nBins+1);
        bidx = nan(T,1);
        bidx(valid) = min(discretize(phiW(valid), edgesP), nBins);
        occ  = accumarray(bidx(valid), 1, [nBins 1]);

        % coherence radius
        ref = cos(phi); ref(isnan(ref)) = 0; ref = ref - mean(ref);
        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
        pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
        [Sb,fb] = mtspectrumc(r.breath(1:T), pB); Sb=Sb(:); fb=fb(:);
        mm = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        if ~any(mm), continue; end
        [~,rl]=max(Sb(mm)); ipk=find(mm,1)+rl-1; f_pk=fb(ipk);
        h=Sb(ipk)/2; lo=ipk; while lo>1 && Sb(lo)>h, lo=lo-1; end
        hi=ipk;            while hi<numel(fb) && Sb(hi)>h, hi=hi+1; end
        f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
        band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];
        pc = struct('Fs',fps, 'tapers',dpsschk([TW,2*TW-1],T,fps), 'pad',0, ...
                    'fpass',band, 'err',0, 'trialave',0);
        nfft  = max(2^(nextpow2(T)+0), T);
        fgrid = getfgrid(fps, nfft, band);  fgrid = fgrid(:);
        mb = fgrid>=band(1) & fgrid<=band(2);
        if ~any(mb), mb = true(size(fgrid)); end
        st  = zeros(T,1); st(ev) = 1;
        lag = round(ca_lag_sec*fps);
        stL = [st(1+lag:end); zeros(lag,1)];

        Q(end+1) = struct('T',T,'fps',fps,'ev',ev,'phiW',phiW,'bidx',bidx, ...
                          'occ',occ,'peak',r.peak,'foot',r.foot,'f_pk',f_pk, ...
                          'stL',stL,'ref',ref,'pc',pc,'mb',mb,'w',T); %#ok<SAGROW>
    end
    if tossedWhole, nTossCell = nTossCell + 1; continue; end
    if isempty(Q),  nThin = nThin + 1; continue; end

    % ---- radii and angles ----
    [logZ_obs, thMean, aPool, wPool] = pooled_rayleigh(Q);
    if numel(aPool) < minSpikes, nThin = nThin + 1; continue; end
    thMed  = circ_median_w_local(aPool, wPool);
    cohObs = pooled_coherence(Q);
    % weighted circular IQR about each centre, in residual space
    qMean = wiqr_local(aPool, wPool, thMean);
    qMed  = wiqr_local(aPool, wPool, thMed);

    % ---- pooled IBI and the PSTH grid, exactly as the figure builds them ----
    allIBI = [];
    for i = 1:numel(Q), allIBI = [allIBI; diff(Q(i).foot)/Q(i).fps]; end %#ok<AGROW>
    IBI = median(allIBI);
    if ~isfinite(IBI) || IBI <= 0, IBI = 1/max(median([Q.f_pk]),eps); end
    fpsRef  = median([Q.fps]);
    winSecH = IBI;  testHalf = IBI/2;
    binW    = histBinFr/fpsRef;
    Mb      = floor(winSecH/binW);
    ctrsC   = (-Mb:Mb)*binW;
    edgesC  = ((-Mb-0.5):(Mb+0.5))*binW;
    nB      = numel(ctrsC);

    cnt = zeros(1,nB); nTrig = 0; nullH = zeros(nShuffle, nB); okNull = true;
    for i = 1:numel(Q)
        T = Q(i).T; fps = Q(i).fps; ev = Q(i).ev;
        wH  = max(1, round(winSecH*fps));
        trg = Q(i).peak(Q(i).peak-wH>=1 & Q(i).peak+wH<=T);
        if isempty(trg), continue; end
        cnt   = cnt + trig_hist_local(ev, trg, wH, edgesC, fps);
        nTrig = nTrig + numel(trg);
        if isempty(ev), continue; end
        minShift = max([wH+1, round(shiftMinCyc/max(Q(i).f_pk,eps)*fps), round(fps)]);
        if T - minShift <= minShift, okNull = false; continue; end

        % circular correlogram, built once; a shift is a rotation of it
        m  = mod(double(ev(:)).' - double(trg(:)), T);
        c  = accumarray(m(:)+1, 1, [T 1]);
        lags   = (-wH:wH).';
        binLag = discretize(lags/fps, edgesC);
        okL    = ~isnan(binLag);
        % Every shift at once. c[mod(l-sh,T)] as a function of sh is a circular
        % shift of the reversed correlogram, so the whole null histogram for ALL
        % T possible offsets is nLags circshifts -- after which any number of
        % shuffles is just row indexing, and nShuffle stops costing anything.
        cRev = c([1, T:-1:2]);                       % u(-m mod T)
        H = zeros(T, nB);                            % row sh+1 = null at shift sh
        Ls = lags(okL);  Bs = binLag(okL);
        for q = 1:numel(Ls)
            H(:,Bs(q)) = H(:,Bs(q)) + circshift(cRev, Ls(q));
        end
        if ~checkedShortcut
            % Row 1 is shift 0, which must reproduce the direct histogram exactly
            % -- otherwise every null is silently off by a rotation.
            assert(isequal(H(1,:), trig_hist_local(ev, trg, wH, edgesC, fps)), ...
                   'the all-shifts table disagrees with the direct histogram at shift 0');
            checkedShortcut = true;
        end
        sh = randi([minShift, T-minShift], nShuffle, 1);
        nullH = nullH + H(sh+1, :);
    end

    p = NaN;
    if nTrig > 0 && okNull
        spkH  = 100*cnt/nTrig;
        nullH = 100*nullH/nTrig;
        nullMu = mean(nullH,1);
        tm = abs(ctrsC) <= testHalf;
        if any(tm)
            sObs  = max(abs(spkH(tm) - nullMu(tm)));
            sNull = max(abs(nullH(:,tm) - nullMu(tm)), [], 2);
            p = (1 + nnz(sNull >= sObs))/(1 + numel(sNull));
        end
    end
    if ~isfinite(p), nNoNull = nNoNull + 1; end

    R(end+1) = struct('cell',Ci.cell, 'stem',Ci.stem, 'group',Ci.group, ...
        'thMean',thMean, 'thMed',thMed, 'qMean',qMean, 'qMed',qMed, ...
        'logZ',logZ_obs, 'cohR',cohObs, ...
        'p',p, 'nEv',nEvC, 'nRecKept',nKept); %#ok<SAGROW>

    if mod(numel(R),50)==0
        fprintf('  %d cells done (%.1f min)\n', numel(R), toc(tAll)/60);
    end
end

n = numel(R);
fprintf('\n%d cells  (%d tossed whole, %d recordings tossed, %d too thin, %d without a null)  in %.1f min\n', ...
        n, nTossCell, nTossRec, nThin, nNoNull, toc(tAll)/60);

thMean = [R.thMean];  thMed = [R.thMed];
qMean  = reshape([R.qMean],2,[]).';   % [n x 2] q25 q75 about the mean
qMed   = reshape([R.qMed], 2,[]).';   % [n x 2] about the median
rZ = [R.logZ];  rC = [R.cohR];  pv = [R.p];
rP = -log10(pv);
gname = string({R.group});
sig = pv <= pCrit;
pFloor = -log10(1/(nShuffle+1));
rP(~isfinite(rP)) = 0;

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[1 1 40 26], ...
             'Name','polar: peak-triggered permutation significance');
set(fig,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);

rLimZ = ceil(max(rZ)*1.15);
COL = { {0.020, 'r = Rayleigh log Z',                             rZ, [0 rLimZ],        NaN}
        {0.275, sprintf('r = coherence |C|  (TW=%g)',TW),         rC, [0 1],            NaN}
        {0.530, 'r = -log_{10}( permutation p )',                 rP, [0 pFloor*1.05],  -log10(pCrit)} };
ROW = { {0.53, thMean, 'MEAN',   qMean}
        {0.06, thMed,  'MEDIAN', qMed } };

present = GRPS(ismember(GRPS, cellstr(unique(gname))));
hLeg = gobjects(numel(present),1);
for ci = 1:3
    for ri = 1:2
        ax = polaraxes(fig,'Position',[COL{ci}{1} ROW{ri}{1} 0.225 0.36]);
        hh = draw_panel(ax, ROW{ri}{2}, ROW{ri}{4}, COL{ci}{3}, sig, gname, ...
                        present, GC, COL{ci}{4}, COL{ci}{5}, pCrit);
        if ci==1 && ri==1, hLeg = hh; end
        if ri == 1
            annotation(fig,'textbox',[COL{ci}{1} 0.895 0.225 0.03],'EdgeColor','none', ...
                'HorizontalAlignment','center','FontWeight','bold','FontSize',10, ...
                'Interpreter','tex','String',COL{ci}{2});
        end
    end
end
for ri = 1:2
    annotation(fig,'textbox',[0.002 ROW{ri}{1}+0.15 0.018 0.07],'EdgeColor','none', ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',10, ...
        'String',ROW{ri}{3});
end

annotation(fig,'textbox',[0.02 0.935 0.96 0.055],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','tex','FontSize',10, ...
    'String',{sprintf(['N = %d cells   |   insp onset = 0 (right), insp peak = \\pi ' ...
                       '(left), counterclockwise   |   columns share the angle, ' ...
                       'rows share the radius'], n), ...
              sprintf(['filled = peak-triggered PSTH permutation p \\leq %.2g ' ...
                       '(%d shuffles, min %d breath cycles) : %d of %d cells      ' ...
                       'arc = weighted circular IQR, significant cells only      ' ...
                       'IO always open grey'], pCrit, nShuffle, shiftMinCyc, ...
                      nnz(sig), n)});
annotation(fig,'textbox',[0.02 0.005 0.96 0.04],'EdgeColor',[0.7 0.7 0.7], ...
    'BackgroundColor',[0.96 0.96 0.96],'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontWeight','bold','FontSize',9, ...
    'String',sprintf(['ONE test marks all six panels: the same peak-triggered p the per-cell ' ...
                      'summary prints. Coherence is a radius only, not permuted.   |   ' ...
                      'p floor = 1/%d = %.1e, so -log10(p) ceilings at %.2f'], ...
                      nShuffle+1, 1/(nShuffle+1), pFloor));

lg = present;
for q = 1:numel(lg)
    m = gname == string(lg{q});
    lg{q} = sprintf('%s (n=%d, sig %d)', lg{q}, nnz(m), nnz(m & sig));
end
lgd = legend(hLeg(isgraphics(hLeg)), lg(isgraphics(hLeg)), 'Box','off','FontSize',9);
lgd.Units = 'normalized';  lgd.Position = [0.775 0.60 0.20 0.22];

axT = axes(fig,'Position',[0.775 0.10 0.22 0.42]); axis(axT,'off');
Lt = strings(0,1);
Lt(end+1) = sprintf('%-8s %4s %7s %7s %6s %5s', ...
                    'genotype','n','meanDg','medLogZ','med|C|','sig');
Lt(end+1) = string(repmat('-',1,44));
for gi = 1:numel(present)
    m = gname == string(present{gi});
    mu = angle(mean(exp(1i*thMean(m))));
    Lt(end+1) = sprintf('%-8s %4d %7.0f %7.2f %6.3f %5d', present{gi}, nnz(m), ...
                        mod(rad2deg(mu),360), median(rZ(m)), median(rC(m)), ...
                        nnz(m & sig)); %#ok<SAGROW>
end
text(axT, 0, 1, Lt, 'FontName','Consolas','FontSize',8, ...
     'VerticalAlignment','top','Interpreter','none');

stem = fullfile(outDir,'polar_recheck_perm6');
exportgraphics(fig, [stem '.png'], 'Resolution',300, 'BackgroundColor','white');
exportgraphics(fig, [stem '.pdf'], 'ContentType','vector','BackgroundColor','white');

T = table([R.cell]', string({R.stem})', string({R.group})', ...
          mod(rad2deg(thMean),360)', mod(rad2deg(thMed),360)', ...
          rZ', rC', pv', sig', (rC >= confC)', ...
          rad2deg(qMean(:,1)), rad2deg(qMean(:,2)), ...
          rad2deg(qMed(:,1)),  rad2deg(qMed(:,2)), ...
          [R.nEv]', [R.nRecKept]', ...
    'VariableNames',{'cell_idx','stem','group','mean_deg','median_deg','logZ', ...
                     'coh_r','p_perm_peak','sig_perm','coh_above_confC', ...
                     'iqr_lo_about_mean_deg','iqr_hi_about_mean_deg', ...
                     'iqr_lo_about_median_deg','iqr_hi_about_median_deg', ...
                     'n_events','n_rec_kept'});
writetable(T, fullfile(outDir,'polar_recheck_perm_percell.csv'));

fprintf('\n%-8s %4s %7s %7s %6s %5s\n','genotype','n','meanDg','medLogZ','med|C|','sig');
for gi = 1:numel(present)
    m = gname == string(present{gi});
    mu = angle(mean(exp(1i*thMean(m))));
    fprintf('%-8s %4d %7.0f %7.2f %6.3f %5d\n', present{gi}, nnz(m), ...
            mod(rad2deg(mu),360), median(rZ(m)), median(rC(m)), nnz(m & sig));
end
fprintf('\nwrote %s.png / .pdf and polar_recheck_perm_percell.csv\n', stem);

%% ===================== LOCAL FUNCTIONS =====================
function [logZ, th, aPool, wPool] = pooled_rayleigh(Q)
aPool = []; wPool = [];
for i = 1:numel(Q)
    b  = Q(i).bidx(Q(i).ev);
    ok = ~isnan(b);
    if ~any(ok), continue; end
    aPool = [aPool; Q(i).phiW(Q(i).ev(ok))];             %#ok<AGROW>
    wPool = [wPool; 1./max(Q(i).occ(b(ok)),1)];          %#ok<AGROW>
end
if numel(aPool) < 2, logZ = -Inf; th = NaN; return; end
[th, Rbar, nEff] = wresultant_local(aPool, wPool);
logZ = log(max(nEff * Rbar^2, eps));
end

function rBar = pooled_coherence(Q)
num = 0; den = 0;
for i = 1:numel(Q)
    st = Q(i).stL;
    [~, Cxy] = coherencyc(Q(i).ref, st - mean(st), Q(i).pc);
    num = num + Q(i).w * mean(Cxy(Q(i).mb));
    den = den + Q(i).w;
end
rBar = num / max(den, eps);
end

function h = trig_hist_local(evIdx, trigIdx, winH, edgesC, fps)
%TRIG_HIST_LOCAL  Verbatim from temporal_phase_cell_fig_260812.m.
h = zeros(1, numel(edgesC)-1);
if isempty(evIdx) || isempty(trigIdx), return; end
d = double(evIdx(:)).' - double(trigIdx(:));
d = d(abs(d) <= winH);
if isempty(d), return; end
h = histcounts(d(:)/fps, edgesC);
end

function h = draw_panel(ax, th, q, r, sig, gname, present, GC, rl, critR, pCrit)
hold(ax,'on');
ax.ThetaZeroLocation = 'right';
ax.ThetaDir          = 'counterclockwise';
% The threshold circle is drawn only on the -log10(p) column, where the radius
% IS the test. On the log Z and |C| columns significance comes from the
% permutation, not from the radius, so a circle there would mark a cut that is
% not the one deciding filled versus open.
if isfinite(critR)
    tt = linspace(0,2*pi,361);
    polarplot(ax, tt, critR*ones(size(tt)), '--','Color',[0.15 0.15 0.15],'LineWidth',0.9);
    text(ax, deg2rad(-18), critR, sprintf('  p=%.2g', pCrit), ...
         'FontSize',7,'Color',[0.15 0.15 0.15]);
end
% Spans first, so the dots sit on top of them. Significant cells only.
for gi = 1:numel(present)
    g = present{gi};
    col = GC.(g);  if strcmp(g,'IO'), col = [0.50 0.50 0.50]; end
    m = find(gname == string(g) & isfinite(th) & sig);
    for j = m(:).'
        if ~all(isfinite(q(j,:))), continue; end
        arc = linspace(th(j)+q(j,1), th(j)+q(j,2), 40);
        polarplot(ax, arc, max(r(j),0)*ones(size(arc)), '-', ...
                  'Color',[col 0.45], 'LineWidth',1.3);
    end
end

h = gobjects(numel(present),1);
for gi = 1:numel(present)
    g = present{gi};  col = GC.(g);
    m = find(gname == string(g) & isfinite(th));
    if isempty(m), continue; end
    if strcmp(g,'IO')
        h(gi) = polarplot(ax, th(m), max(r(m),0), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',[0.50 0.50 0.50], ...
            'LineWidth',0.7,'MarkerSize',5,'LineStyle','none');
        continue
    end
    ms = m(sig(m));  mn = m(~sig(m));  hn = gobjects(0);
    if ~isempty(mn)
        hn = polarplot(ax, th(mn), max(r(mn),0), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',col, ...
            'LineWidth',0.9,'MarkerSize',5,'LineStyle','none');
    end
    if ~isempty(ms)
        h(gi) = polarplot(ax, th(ms), max(r(ms),0), 'o', ...
            'MarkerFaceColor',col,'MarkerEdgeColor','w', ...
            'LineWidth',0.5,'MarkerSize',5,'LineStyle','none');
    else
        h(gi) = hn;
    end
end
rlim(ax, rl);
ax.ThetaTick      = 0:45:315;
ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:45:315, 'uni',0);
ax.RAxisLocation  = 180;
ax.GridAlpha      = 0.15;
ax.FontSize       = 7;
end

function [th, Rbar, nEff] = wresultant_local(a, w)
S1 = sum(w);  S2 = sum(w.^2);
if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
nEff = S1^2 / max(S2, eps);
v    = sum(w(:) .* exp(1i*a(:))) / S1;
th   = angle(v);
Rbar = min(abs(v), 1);
end

function q = wiqr_local(a, w, centre)
%WIQR_LOCAL  Weighted circular IQR as residual offsets about a given centre.
%  Returns [q25 q75], both signed offsets in radians, so the arc runs from
%  centre+q25 to centre+q75. Circularity is removed by measuring every event as
%  its wrapped offset from the centre before taking ordinary weighted
%  percentiles -- the same construction polar_selected_260816.m uses.
q = [NaN NaN];
if numel(a) < 2 || ~isfinite(centre), return; end
u = angle(exp(1i*(a(:) - centre)));
q = [wprctile_local(u, w, 25), wprctile_local(u, w, 75)];
end

function qq = wprctile_local(v, w, pr)
%WPRCTILE_LOCAL  Weighted percentile, midpoint rule on the cumulative weight.
%  Verbatim from polar_selected_260816.m.
[v, o] = sort(v(:));  w = w(o);  w = w / sum(w);
c = cumsum(w) - 0.5*w;
if numel(v) < 2, qq = v(1); return; end
[c, iu] = unique(c);  v = v(iu);
qq = interp1(c, v, pr/100, 'linear', 'extrap');
end

function m = circ_median_w_local(a, w)
g = linspace(-pi, pi, 361);
cost = arrayfun(@(x) sum(w(:) .* abs(angle(exp(1i*(a(:)-x))))), g);
[~,i] = min(cost);
g2 = linspace(g(max(i-1,1)), g(min(i+1,numel(g))), 41);
cost2 = arrayfun(@(x) sum(w(:) .* abs(angle(exp(1i*(a(:)-x))))), g2);
[~,j] = min(cost2);
m = g2(j);
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
