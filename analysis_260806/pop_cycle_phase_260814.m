% pop_cycle_phase_260814.m
% -----------------------------------------------------------------------
%  Breath phase two ways, on the SAME 20 randomly chosen cells (logZ > 2), so the
%  parameterisation is the only thing that differs between the two figures.
%
%   'cycle'    onset = 0, NEXT ONSET = 2pi, linear in time. The peak falls where
%              it falls (measured: 0.23pi on average). Occupancy nearly flat.
%   'landmark' onset = 0, PEAK = pi, next onset = 2pi, piecewise linear in time
%              between landmarks. Aligns the inspiratory burst across cycles even
%              when cycle duration varies -- but stretches the fast inspiration
%              over half the circle, so occupancy becomes badly uneven.
%
%  Each figure: the breath waveform over 0 to 4pi (every cycle-pair faint, mean
%  on top), then 20 polar event histograms.
%
%  OCCUPANCY CORRECTION IS APPLIED IN BOTH. Each event is weighted by
%  1/(frames its own recording spent in that phase bin), the ventral convention.
%  Without it a cell firing at a constant rate looks tuned to wherever the phase
%  happens to move slowly. The measured occupancy ratio is printed per mode so
%  the size of the correction is visible rather than assumed.
%
%  Runqi Zhang / 2026-08-14
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath,'pop_analysis_260813');
L   = load(fullfile(outDir,'pop_features.mat'),'T');  T = L.T;
Din = load(fullfile(rootPath,'event_latency_260811','event_latency_data.mat'), ...
           'CELL','OBS','REC');
% T is post-merge, event_latency_data.mat is pre-merge -- pool through this.
obsOf = pooled_obs_260814(Din.CELL, Din.OBS);

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(260814);

logZcut  = 2;
nShowCel = 20;
nPolBins = 24;          % 15 deg bins
nCycPts  = 200;         % resample points per cycle
nDrop    = 30;
maxFaint = 400;

sel  = find(T.logZ >= logZcut);
pick = sel(randperm(numel(sel), min(nShowCel,numel(sel))));   % ONE draw, both figures
fprintf('cells with logZ >= %g : %d   showing the same %d in both figures\n', ...
        logZcut, numel(sel), numel(pick));

MODES  = {'cycle','landmark'};
% Both phase frames are drawn twice: occupancy-CORRECTED (each event weighted by
% 1/frames its own recording spent in that bin) and RAW counts. The raw version
% is the honest picture of where the spikes actually fell; the corrected version
% is what you must use to claim tuning, because in the landmark frame expiration
% occupies most of the real time inside the second half-circle and an unmodulated
% cell will look expiratory. Four figures, one cell draw.
CORRMODE = {true, false};
CORRTAG  = {'occCorr','rawCount'};
EVp    = cell(numel(pick),2);   % event phases
EVw    = cell(numel(pick),2);   % occupancy weights for those events
CYC    = {[],[]};               % pooled cycle-pair waveforms
occTot = {zeros(nPolBins,1), zeros(nPolBins,1)};
edgesP = linspace(0, 2*pi, nPolBins+1);
INF    = struct('lab',{},'logZ',{},'grp',{},'nEv',{});

for ii = 1:numel(pick)
    % pick indexes ROWS of T. The cell NUMBER is T.cell(row) -- those were the same
    % thing only while T held all 276 cells in order, which the 260814 merges ended.
    c = T.cell(pick(ii));
    if isempty(obsOf{c}), continue; end
    for o = obsOf{c}(:)'
        p  = regexp(Din.OBS(o).label,'/','split');
        fp = Din.REC(Din.OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        roi = str2double(p{end});
        [fps,~] = detect_session_fps(fp, 30); %#ok<ASGLU>
        BP = load(fullfile(fp,'breath_peak_pc1.mat'));
        IP = load(fullfile(fp,'breath_insp_start_pc1.mat'));
        CA = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
        nB = numel(BP.breath);
        bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = [];
        bw = bw - mean(bw);
        if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
            ep = double(BP.insp_onsets_train(:) ~= 0);
        else
            ep = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ep(oi(oi>=1 & oi<=nB)) = 1;
        end
        ep(1:min(nDrop,numel(ep))) = [];
        ef = zeros(nB,1); fi = round(IP.insp_start_idx(:));
        ef(fi(fi>=1 & fi<=nB)) = 1; ef(1:min(nDrop,numel(ef))) = [];
        if ~isempty(regexpi(fp,'Vglut2','once')) && strcmp(p{2},'1124')
            bw = [bw(1); bw(1:end-1)]; ep = [0; ep(1:end-1)]; ef = [0; ef(1:end-1)];
        end
        stk = double(CA.roi_spikes(roi).spike_train(:));
        Tn  = min([numel(bw), numel(ep), numel(ef), numel(stk)]);
        bw = bw(1:Tn); ep = ep(1:Tn); ef = ef(1:Tn); stk = stk(1:Tn);
        foot = find(ef>0); peak = find(ep>0);
        if numel(foot) < 3, continue; end
        evIdx = find(stk>0);

        for md = 1:2
            if md == 1
                phi = cycle_phase_local(foot, Tn);
            else
                phi = landmark_phase_local(peak, foot, Tn);
            end
            v = ~isnan(phi);
            bidx = nan(Tn,1);
            bidx(v) = min(discretize(mod(phi(v),2*pi), edgesP), nPolBins);
            occ = accumarray(bidx(v), 1, [nPolBins 1]);
            occTot{md} = occTot{md} + occ;
            e = evIdx(~isnan(bidx(evIdx)));
            EVp{ii,md} = [EVp{ii,md}; mod(phi(e),2*pi)];
            EVw{ii,md} = [EVw{ii,md}; 1./max(occ(bidx(e)),1)];   % own-recording weights
        end

        % cycle-pair waveform, resampled, for each mode's own alignment
        bwz = (bw - mean(bw))/max(std(bw),eps);
        for k = 1:numel(foot)-2
            a = foot(k); b = foot(k+1); d = foot(k+2);
            if d <= b || b <= a, continue; end
            % mode 1: each cycle spans 0..2pi by time
            s1 = interp1(0:(b-a), bwz(a:b), linspace(0,b-a,nCycPts),'linear');
            s2 = interp1(0:(d-b), bwz(b:d), linspace(0,d-b,nCycPts),'linear');
            CYC{1}(end+1,:) = [s1 s2]; %#ok<AGROW>
            % mode 2: onset->peak fills 0..pi, peak->onset fills pi..2pi
            pk1 = peak(peak>a & peak<b); pk2 = peak(peak>b & peak<d);
            if isempty(pk1) || isempty(pk2), continue; end
            pk1 = pk1(1); pk2 = pk2(1);
            h = nCycPts/2;
            q1 = [interp1(0:(pk1-a), bwz(a:pk1), linspace(0,pk1-a,h),'linear'), ...
                  interp1(0:(b-pk1), bwz(pk1:b), linspace(0,b-pk1,h),'linear')];
            q2 = [interp1(0:(pk2-b), bwz(b:pk2), linspace(0,pk2-b,h),'linear'), ...
                  interp1(0:(d-pk2), bwz(pk2:d), linspace(0,d-pk2,h),'linear')];
            CYC{2}(end+1,:) = [q1 q2]; %#ok<AGROW>
        end
    end
    INF(ii) = struct('lab',sprintf('cell %d',c), 'logZ',T.logZ(T.cell==c), ...
                     'grp',T.group{T.cell==c}, 'nEv',numel(EVp{ii,1}));
end

%% ===================== FIGURES =====================
ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;
kgrid = linspace(0, 2*pi, 361);
for md = 1:2
for cm = 1:2
    doCorr = CORRMODE{cm};
    occR = max(occTot{md})/max(min(occTot{md}),1);
    fprintf('\n%-9s %-9s occupancy ratio %.2fx   |  %d cycle-pairs\n', ...
            MODES{md}, CORRTAG{cm}, occR, size(CYC{md},1));

    fig = figure('Color','w','Units','normalized','Position',[0.03 0.04 0.92 0.86]);
    tl = tiledlayout(fig,5,5,'TileSpacing','compact','Padding','compact');
    if md == 1
        sub = 'onset = 0, next onset = 2\pi, linear in time (peak NOT pinned)';
    else
        sub = 'onset = 0, PEAK = \pi, next onset = 2\pi (piecewise linear)';
    end
    if doCorr, ctxt = sprintf('occupancy %.2fx, CORRECTED', occR);
    else,      ctxt = sprintf('occupancy %.2fx, RAW COUNTS (uncorrected)', occR); end
    title(tl, sprintf(['%s   |   %d cells logZ > %g, same %d shown   |   %s   |   ' ...
          'grey = histogram, black = von Mises KDE'], ...
          sub, numel(sel), logZcut, numel(pick), ctxt), 'FontWeight','bold');

    ax = nexttile(tl,1,[1 5]); hold(ax,'on');
    xg = linspace(0, 4*pi, 2*nCycPts);
    C  = CYC{md};
    nF = min(maxFaint, size(C,1));
    plot(ax, xg, C(randperm(size(C,1),nF),:)', '-','Color',[0.6 0.6 0.6 0.08],'LineWidth',0.3);
    mu = mean(C,1,'omitnan'); sd = std(C,0,1,'omitnan');
    fill(ax,[xg fliplr(xg)],[mu+sd fliplr(mu-sd)],[0.2 0.4 0.8], ...
         'FaceAlpha',0.20,'EdgeColor','none');
    plot(ax, xg, mu, '-','Color',[0.1 0.2 0.6],'LineWidth',2);
    for xv = [0 2*pi 4*pi], xline(ax,xv,'-','Color',[0.90 0.10 0.10],'LineWidth',1.4); end
    [~,imx] = max(mu(1:nCycPts));
    xline(ax, xg(imx), '--','Color',[0.35 0.75 1.00],'LineWidth',1.2);
    xlim(ax,[0 4*pi]); xticks(ax,0:pi:4*pi);
    xticklabels(ax,{'0','\pi','2\pi','3\pi','4\pi'});
    xlabel(ax,'phase (red = inspiration onset)'); ylabel(ax,'breath (z per recording)');
    title(ax, sprintf('breath waveform, %d cycle-pairs   |   mean peak at %.2f\\pi', ...
          size(C,1), xg(imx)/pi), 'FontSize',10);
    box(ax,'on'); grid(ax,'on');

    for ii = 1:numel(pick)
        pax = polaraxes(tl); pax.Layout.Tile = 5 + ii;
        a = EVp{ii,md};
        if isempty(a), title(pax,'(no events)','FontSize',7); continue; end
        if doCorr, w = EVw{ii,md}; else, w = ones(size(a)); end

        b = min(discretize(a, edgesP), nPolBins);
        rate = accumarray(b, w, [nPolBins 1]);
        pct  = 100*rate/max(sum(rate),eps);

        % von Mises KDE on the same weights, scaled to the histogram's units so
        % the two are directly comparable on one axis
        [fk, kap] = circ_kde_local(a, w, kgrid);
        fk = fk / max(sum(fk)*mean(diff(kgrid)), eps);      % density, integrates to 1
        fkPct = fk * (2*pi/nPolBins) * 100;                 % -> % per histogram bin

        hold(pax,'on');
        polarplot(pax,[ctrsP ctrsP(1)],[pct' pct(1)],'-','Color',[0.62 0.62 0.62], ...
                  'LineWidth',0.9);
        polarplot(pax, kgrid, fkPct, '-','Color',[0.10 0.10 0.10],'LineWidth',1.6);
        vv = sum(w.*exp(1i*a))/sum(w);
        rmax = max([pct(:); fkPct(:)]);
        polarplot(pax,[0 angle(vv)],[0 rmax*abs(vv)/max(abs(vv),eps)],'-', ...
                  'Color',[0.85 0.20 0.20],'LineWidth',1.6);
        pax.ThetaZeroLocation='right'; pax.ThetaDir='counterclockwise';
        thetaticks(pax,0:90:270); pax.RTickLabel=[]; pax.FontSize=6;
        title(pax, sprintf('%s %s\nlogZ %.1f, %d ev, \\kappa %.0f', INF(ii).grp, ...
              INF(ii).lab, INF(ii).logZ, INF(ii).nEv, kap), 'FontSize',7);
    end

    stem = sprintf('cycle_phase_logZ2_%s_%s', MODES{md}, CORRTAG{cm});
    exportgraphics(fig, fullfile(outDir,[stem '.png']),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
    fprintf('  saved %s.png / .pdf\n', stem);
    close(fig);
end
end

%% ===================== LOCAL =====================
function [f, kappa] = circ_kde_local(a, w, grid)
%CIRC_KDE_LOCAL  Weighted von Mises kernel density on the circle.
%  The circular analogue of a Gaussian KDE: each event contributes a von Mises
%  bump of concentration kappa (kappa is the inverse of a bandwidth -- larger =
%  narrower). Bandwidth from Taylor's (2008) rule of thumb, using Kish's
%  effective n so the occupancy weights do not inflate the sample size.
%
%      f(t) = sum_i w_i exp(kappa*cos(t - a_i)) / (2*pi*I0(kappa))
%
%  Computed with the EXPONENTIALLY SCALED Bessel function and cos()-1 in the
%  exponent, because I0(kappa) overflows for kappa > ~700 and sharply tuned cells
%  reach that easily.
a = a(:);  w = w(:);  w = w/sum(w);
n_eff = 1/sum(w.^2);                       % Kish
R = abs(sum(w.*exp(1i*a)));
kh = kappa_mle_local(R);
kappa = ( 3*n_eff*kh^2*besseli(2,2*kh,1) / ...
          (4*sqrt(pi)*besseli(0,kh,1)^2) )^(2/5);
if ~isfinite(kappa) || kappa <= 0, kappa = 1; end
kappa = min(max(kappa, 0.5), 500);
f = zeros(size(grid));
for i = 1:numel(a)
    f = f + w(i)*exp(kappa*(cos(grid - a(i)) - 1));
end
f = f / (2*pi*besseli(0,kappa,1));
end

function k = kappa_mle_local(R)
%KAPPA_MLE_LOCAL  Standard piecewise approximation to the von Mises MLE.
if R < 0.53
    k = 2*R + R^3 + 5*R^5/6;
elseif R < 0.85
    k = -0.4 + 1.39*R + 0.43/(1-R);
else
    k = 1/(R^3 - 4*R^2 + 3*R);
end
if ~isfinite(k) || k <= 0, k = 0.1; end
end

function phi = cycle_phase_local(foot, T)
%CYCLE_PHASE_LOCAL  onset = 0, next onset = 2pi, linear in time. Peak unpinned.
phi = nan(T,1);
for k = 1:numel(foot)-1
    a = foot(k); b = foot(k+1);
    if b <= a, continue; end
    phi(a:b-1) = 2*pi*(0:(b-a-1))'/(b-a);
end
end

function phi = landmark_phase_local(peak, foot, T)
%LANDMARK_PHASE_LOCAL  onset = 0, PEAK = pi, next onset = 2pi, piecewise linear
%  in time between landmarks, with strict alternation enforced.
phi = nan(T,1);
ev = [peak(:); foot(:)];
ty = [ones(numel(peak),1); zeros(numel(foot),1)];
[ev,o] = sort(ev); ty = ty(o);
keep = true(size(ev));
for i = 2:numel(ev), if ty(i)==ty(i-1), keep(i) = false; end, end
ev = ev(keep); ty = ty(keep);
if numel(ev) < 2, return; end
ph = nan(size(ev)); cur = ty(1)*pi;
for i = 1:numel(ev), ph(i) = cur; cur = cur + pi; end
for i = 1:numel(ev)-1
    a = ev(i); b = ev(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(ph(i), ph(i+1), b-a+1);
end
end
