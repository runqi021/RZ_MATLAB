% pop_event_hist_260814.m
% -----------------------------------------------------------------------
%  POPULATION event histogram on breath phase, 4 panels:
%
%                    cycle phase (peak free)      landmark phase (peak = pi)
%    p < 0.01        pooled events                pooled events
%    p >= 0.01       pooled events                pooled events
%
%  p = the PSTH permutation p; a cell counts as significant if EITHER trigger
%  clears the cut, i.e. min(pOnset,pPeak) < 0.01.
%
%  OCCUPANCY CORRECTED. Every event is weighted by 1/(frames its own recording
%  spent in that phase bin). Mandatory in the landmark frame, where occupancy is
%  ~8x uneven and an unmodulated cell would otherwise look expiratory.
%
%  TWO NORMALISATIONS ARE DRAWN, because "population histogram" is ambiguous:
%    bars  = all events pooled. Honest about total spikes, but a cell with 490
%            events outweighs one with 6, so a handful of cells set the shape.
%    line  = each cell's profile normalised to sum 1 first, then averaged across
%            cells. One cell, one vote.
%  If the two disagree, the pooled shape is being driven by a few busy cells.
%
%  No KDE anywhere.
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
% T is post-merge, event_latency_data.mat is pre-merge: index the pooled list, not
% CELL(c).obs, or a merged cell contributes only one of its recordings.
obsOf = pooled_obs_260814(Din.CELL, Din.OBS);

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

nBins  = 36;                 % 10 deg
nDrop  = 30;
pcut   = 0.01;
edgesP = linspace(0, 2*pi, nBins+1);
ctrsP  = (edgesP(1:end-1)+edgesP(2:end))/2;

nCell = height(T);
W = zeros(nCell, nBins, 2);      % occupancy-weighted event mass, per cell, per frame
NEV = zeros(nCell,1);
peakPhase = nan(nCell,2);        % where the breath peak sits, per frame

fprintf('extracting %d cells ...\n', nCell);
tic;
for ci = 1:nCell
    c = T.cell(ci);
    if isempty(obsOf{c}), continue; end
    for o = obsOf{c}(:)'
        p  = regexp(Din.OBS(o).label,'/','split');
        fp = Din.REC(Din.OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        roi = str2double(p{end});
        BP = load(fullfile(fp,'breath_peak_pc1.mat'));
        IP = load(fullfile(fp,'breath_insp_start_pc1.mat'));
        CA = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
        nB = numel(BP.breath);
        if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
            ep = double(BP.insp_onsets_train(:) ~= 0);
        else
            ep = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ep(oi(oi>=1 & oi<=nB)) = 1;
        end
        ep(1:min(nDrop,numel(ep))) = [];
        ef = zeros(nB,1); fi = round(IP.insp_start_idx(:));
        ef(fi(fi>=1 & fi<=nB)) = 1; ef(1:min(nDrop,numel(ef))) = [];
        if ~isempty(regexpi(fp,'Vglut2','once')) && strcmp(p{2},'1124')
            ep = [0; ep(1:end-1)]; ef = [0; ef(1:end-1)];
        end
        stk = double(CA.roi_spikes(roi).spike_train(:));
        Tn  = min([numel(ep), numel(ef), numel(stk)]);
        ep = ep(1:Tn); ef = ef(1:Tn); stk = stk(1:Tn);
        foot = find(ef>0); peak = find(ep>0);
        if numel(foot) < 3, continue; end
        ev = find(stk>0);

        for md = 1:2
            if md == 1, phi = cycle_phase_local(foot, Tn);
            else,       phi = landmark_phase_local(peak, foot, Tn); end
            v = ~isnan(phi);
            bidx = nan(Tn,1);
            bidx(v) = min(discretize(mod(phi(v),2*pi), edgesP), nBins);
            occ = accumarray(bidx(v), 1, [nBins 1]);
            e   = ev(~isnan(bidx(ev)));
            if isempty(e), continue; end
            % 1/occupancy weight from THIS recording
            W(ci,:,md) = W(ci,:,md) + ...
                accumarray(bidx(e), 1./max(occ(bidx(e)),1), [nBins 1])';
            % where does the breath peak land in this frame?
            pk = peak(~isnan(bidx(peak)));
            if ~isempty(pk), peakPhase(ci,md) = circ_mean_local(mod(phi(pk),2*pi)); end
        end
        NEV(ci) = NEV(ci) + numel(ev);
    end
    if mod(ci,50)==0, fprintf('  %d/%d (%.0f s)\n', ci, nCell, toc); end
end
fprintf('done in %.1f min\n', toc/60);

%% ===================== GROUPS =====================
pMin = min([T.pOnset T.pPeak], [], 2);
SETS = { sprintf('p < %g  (either trigger)', pcut),  pMin <  pcut
         sprintf('p >= %g (neither)',        pcut),  pMin >= pcut };
MODEN = {'cycle phase (onset = 0, next onset = 2\pi)', ...
         'landmark phase (onset = 0, peak = \pi)'};

fig = figure('Color','w','Units','normalized','Position',[0.05 0.08 0.80 0.78]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf(['population event histogram, occupancy-corrected   |   %d active cells' ...
      '   |   bars = all events pooled, line = per-cell mean (one cell one vote)'], nCell), ...
      'FontWeight','bold');

for si = 1:2
    for md = 1:2
        ax = nexttile(tl); hold(ax,'on');
        m = SETS{si,2} & any(W(:,:,md)>0,2);
        Wm = W(m,:,md);
        pooled = sum(Wm,1);  pooled = 100*pooled/max(sum(pooled),eps);
        perCell = Wm ./ max(sum(Wm,2),eps);
        mu = 100*mean(perCell,1);
        se = 100*std(perCell,0,1)/sqrt(max(size(perCell,1),1));

        bar(ax, ctrsP, pooled, 1, 'FaceColor',[0.55 0.55 0.55], ...
            'EdgeColor','none','FaceAlpha',0.85);
        fill(ax,[ctrsP fliplr(ctrsP)],[mu+se fliplr(mu-se)],[0.1 0.3 0.7], ...
             'FaceAlpha',0.20,'EdgeColor','none');
        plot(ax, ctrsP, mu, '-','Color',[0.1 0.2 0.6],'LineWidth',1.8);
        yline(ax, 100/nBins, ':','Color',[0.4 0.4 0.4],'LineWidth',1);   % uniform

        xline(ax, 0, '-','Color',[0.90 0.10 0.10],'LineWidth',1.4);
        pk = circ_mean_local(peakPhase(m,md));
        xline(ax, pk, '--','Color',[0.35 0.75 1.00],'LineWidth',1.4);

        xlim(ax,[0 2*pi]); xticks(ax,0:pi/2:2*pi);
        xticklabels(ax,{'0','\pi/2','\pi','3\pi/2','2\pi'});
        xlabel(ax,'breath phase   (red = onset, blue dashed = mean peak)');
        ylabel(ax,'% of events');
        title(ax, sprintf('%s   |   %s   |   n = %d cells, %d events', ...
              SETS{si,1}, MODEN{md}, nnz(m), sum(NEV(m))), 'FontSize',9);
        box(ax,'on'); grid(ax,'on');
        fprintf('%-26s %-12s n=%3d cells  %6d events  peak at %.2f pi\n', ...
                SETS{si,1}, MODEN{md}(1:5), nnz(m), sum(NEV(m)), pk/pi);
    end
end

exportgraphics(fig, fullfile(outDir,'pop_event_hist_p01.png'), ...
               'Resolution',200,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'pop_event_hist_p01.pdf'), ...
               'ContentType','vector','BackgroundColor','white');
fprintf('\nsaved pop_event_hist_p01.png / .pdf\n');

%% ===================== LOCAL =====================
function m = circ_mean_local(a)
a = a(~isnan(a));
if isempty(a), m = NaN; return; end
m = mod(angle(mean(exp(1i*a(:)))), 2*pi);
end

function phi = cycle_phase_local(foot, T)
phi = nan(T,1);
for k = 1:numel(foot)-1
    a = foot(k); b = foot(k+1);
    if b <= a, continue; end
    phi(a:b-1) = 2*pi*(0:(b-a-1))'/(b-a);
end
end

function phi = landmark_phase_local(peak, foot, T)
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
