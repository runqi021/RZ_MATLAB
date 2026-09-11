% polar_hist_cycle_vglut2_260902.m
% -----------------------------------------------------------------------
%  Population spike-phase polar histogram on a SINGLE-LANDMARK cycle phase.
%  Vglut2 only, two panels: vagal-intact and vagotomised, all active cells.
%
%  THE PHASE. One landmark, not two: phase = 2*pi * (t - onset_k) / (onset_k+1 -
%  onset_k), linear in TIME across the whole cycle, from one inspiration onset to
%  the next. The inspiratory peak is NOT pinned anywhere -- it falls where it
%  falls, and its median position is drawn as the sky dashed ray.
%
%  WHY THIS AND NOT THE TWO-LANDMARK PHASE. The piecewise version maps
%  inspiration onto [0, pi] and expiration onto [pi, 2*pi] whatever their real
%  durations, so it stretches the short half and compresses the long one by a
%  factor that changes cycle by cycle. Every bin then covers a different amount
%  of time, and the amount changes with breathing rate -- which is exactly the
%  variable that differs between the two panels here.
%
%  AND IT REMOVES THE OCCUPANCY PROBLEM ENTIRELY. Under this phase every bin
%  spans T_cycle / nBins seconds in EVERY cycle, so the dwell time is uniform by
%  construction and pooling cycles of different lengths keeps it uniform. The
%  occupancy compensation the old two-landmark histogram needed -- the
%  (Obar/O_b) reweighting, made necessary by a ~7x uneven dwell -- is not needed
%  and is not applied. Raw and compensated are the same curve here.
%
%  UNITS: EVENTS / CYCLE / CELL, POOLED -- EVERY EVENT COUNTS EQUALLY. The
%  numerator is every event of every cell in the group, binned together; the
%  denominator is the total number of (cell x cycle) observations those cells
%  contributed. No per-cell normalisation happens before the pooling, so a cell
%  that fires ten times as often carries ten times the weight in the shape of the
%  curve, which is the point: this is the population's event distribution, not
%  the average of its cells' distributions.
%
%  A cycle in which a cell fired nothing still increments the denominator, so the
%  radius is a rate and not a conditional distribution -- a group that fires less
%  draws a smaller curve rather than the same curve rescaled.
%
%  COORDINATES as in every other figure of this set: 0 at east (= inspiration
%  onset), angles increasing counterclockwise.
%
%  ALL ACTIVE CELLS -- the whole bundle for this genotype. Each panel overlays
%  the SIGNIFICANT cells (green) and the NON-SIGNIFICANT ones (grey), each as its
%  own events/cycle/cell mean, so the two curves are on the same scale and the
%  grey one is the honest "what does a cell that failed the test look like".
%  Significance is the circular-shift permutation on the peak-triggered PSTH,
%  p <= 0.01, read from polar_recheck_perm_percell.csv -- the same test as
%  everywhere else in this set. Tossed recordings and tossed cells are dropped,
%  and the events are the re-check curated ones.
%
%  Runqi Zhang / 2026-09-02
% -----------------------------------------------------------------------
clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot,'analysis_260806'));

%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
bundleDir = fullfile(sumRoot,'per-cell-summary_active_260812','spike_recheck_260901');
outDir    = fullfile(sumRoot,'polar_recheck_260901');

GENOTYPE  = 'Vglut2';
VAG_DATES = {'0824'};
TRIGGER   = 'onset';     % 'onset' = onset-to-onset cycle | 'peak' = peak-to-peak
nB        = 36;          % phase bins (10 deg)
maxCycFac = 3;           % drop cycles longer than this x the recording median
% =========================================================

COL   = [0.10 0.65 0.20];      % significant, Vglut2 green
NCOL  = [0.55 0.55 0.55];      % non-significant, grey
TCOL  = [0.85 0.20 0.10];      % the trigger ray, red
PCOL  = [0.35 0.75 1.00];      % the other landmark's median, sky
pCrit = 0.01;

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

edges = linspace(0, 2*pi, nB+1);
ctrs  = (edges(1:end-1) + edges(2:end))/2;

%% ===================== COLLECT =====================
CU = ca_recheck_load_curation_260901(bundleDir);
d  = dir(fullfile(bundleDir,'cells','*.mat'));

sigCsv = fullfile(outDir,'polar_recheck_perm_percell.csv');
assert(isfile(sigCsv), 'Run polar_recheck_perm_260901.m first (missing %s)', sigCsv);
Tsig   = readtable(sigCsv,'TextType','string');
sigMap = containers.Map(cellstr(Tsig.stem), num2cell(logical(Tsig.sig_perm)));

% {group, 1=sig | 2=nonsig}: POOLED counts and pooled (cell x cycle) count
cntG  = zeros(2,2,nB);
cycG  = zeros(2,2);
nCyc  = zeros(2,1);
nCell = zeros(2,2);
nEvG  = zeros(2,1);
othPh = {[], []};         % the OTHER landmark's phase, per cycle
nDrop = 0;

for k = 1:numel(d)
    stem = erase(d(k).name,'.mat');
    if ~startsWith(stem, [GENOTYPE '_']), continue; end
    gi = 1;
    for v = 1:numel(VAG_DATES)
        if contains(stem, ['_' VAG_DATES{v} '_']), gi = 2; end
    end
    L = load(fullfile(d(k).folder, d(k).name),'RECc');  REC = L.RECc;
    used = false;
    cntC = zeros(1,nB);  cycC = 0;      % THIS cell's own counts and cycles
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && (g.toss || g.cell_toss), continue; end
        ev = r.spike_idx(:);
        if g.found, ev = g.spike_idx(:); end
        if strcmpi(TRIGGER,'peak')
            trg = r.peak(:);  oth = r.foot(:);
        else
            trg = r.foot(:);  oth = r.peak(:);
        end
        if numel(trg) < 2, continue; end
        used = true;

        dt   = diff(trg);
        good = dt <= maxCycFac * median(dt);      % a missed landmark makes a
        nDrop = nDrop + nnz(~good);               % cycle several times too long
        for c = 1:numel(dt)
            if ~good(c), continue; end
            a = trg(c); b = trg(c+1);
            e = ev(ev >= a & ev < b);
            if ~isempty(e)
                ph = 2*pi*(e - a)/(b - a);
                cntC     = cntC + histcounts(ph, edges);
                nEvG(gi) = nEvG(gi) + numel(e);
            end
            o = oth(oth > a & oth < b);           % where the other landmark sits
            if ~isempty(o)
                othPh{gi}(end+1,1) = 2*pi*(o(1) - a)/(b - a); %#ok<SAGROW>
            end
            cycC     = cycC + 1;
            nCyc(gi) = nCyc(gi) + 1;
        end
    end
    if used && cycC > 0
        si = 2;                                     % 1 = sig, 2 = non-sig
        if isKey(sigMap, stem) && sigMap(stem), si = 1; end
        nCell(gi,si) = nCell(gi,si) + 1;
        % POOLED: this cell's raw counts and its cycles go into the group total
        % WITHOUT being divided first, so every event carries the same weight.
        cntG(gi,si,:) = squeeze(cntG(gi,si,:)).' + cntC;
        cycG(gi,si)   = cycG(gi,si) + cycC;
    end
end

% pooled events / cycle / cell: total events in the bin over total cell-cycles
H = zeros(2, 2, nB);
for gi = 1:2
    for si = 1:2
        if cycG(gi,si) > 0, H(gi,si,:) = cntG(gi,si,:) / cycG(gi,si); end
    end
end
fprintf('%s, %s-to-%s cycles\n', GENOTYPE, TRIGGER, TRIGGER);
fprintf('  dropped %d cycles longer than %gx the recording median\n', nDrop, maxCycFac);
NAME  = {'vagal intact','vagotomised'};
SNAME = {'sig','non-sig'};
for gi = 1:2
    for si = 1:2
        h = squeeze(H(gi,si,:)).';
        if nCell(gi,si) == 0, continue; end
        fprintf('  %-14s %-7s %2d cells %6d cell-cycles %5d events   peak %.4f at %.0f deg\n', ...
                NAME{gi}, SNAME{si}, nCell(gi,si), cycG(gi,si), ...
                round(sum(squeeze(cntG(gi,si,:)))), max(h), ...
                rad2deg(ctrs(find(h == max(h), 1))));
    end
end

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[2 3 26 14], ...
             'Name','Vglut2 cycle-phase polar histogram');
set(fig,'DefaultAxesFontSize',9);
rMax = max(H(:))*1.10;
hLeg = gobjects(2,1);

for gi = 1:2
    ax = polaraxes(fig,'Position',[0.055+(gi-1)*0.470, 0.07, 0.355, 0.60]);
    hold(ax,'on');
    ax.ThetaZeroLocation = 'right';
    ax.ThetaDir          = 'counterclockwise';

    % the trigger sits at 0 by definition; the other landmark's median is drawn
    polarplot(ax, [0 0], [0 rMax], '--', 'Color',TCOL, 'LineWidth',1.1);
    if ~isempty(othPh{gi})
        mo = angle(mean(exp(1i*othPh{gi})));
        polarplot(ax, [mo mo], [0 rMax], '--', 'Color',PCOL, 'LineWidth',1.1);
    else
        mo = NaN;
    end

    % closed polygons, so each histogram joins up across 0. Grey first, so the
    % significant curve is never hidden behind it.
    for si = [2 1]
        if nCell(gi,si) == 0, continue; end
        h = squeeze(H(gi,si,:)).';
        cl = COL; lw = 2.0;
        if si == 2, cl = NCOL; lw = 1.5; end
        hh = polarplot(ax, [ctrs ctrs(1)], [h h(1)], '-', 'Color',cl, 'LineWidth',lw);
        if gi == 1, hLeg(si) = hh; end
    end

    rlim(ax,[0 rMax]);
    ax.ThetaTick      = 0:30:330;
    ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:30:330, 'uni',0);
    ax.RAxisLocation  = 180;
    ax.GridAlpha      = 0.15;
    ax.FontSize       = 8;

    annotation(fig,'textbox',[0.055+(gi-1)*0.470, 0.735, 0.355, 0.045], ...
        'EdgeColor','none','HorizontalAlignment','center', ...
        'FontWeight','bold','FontSize',10,'Interpreter','none', ...
        'String',sprintf('%s %s  |  %d sig, %d non-sig', GENOTYPE, NAME{gi}, ...
                         nCell(gi,1), nCell(gi,2)));
    annotation(fig,'textbox',[0.055+(gi-1)*0.470, 0.697, 0.355, 0.04], ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',8.5, ...
        'Interpreter','none', ...
        'String',sprintf('%d cycles, %d events   |   median %s at %.0f deg', ...
                         nCyc(gi), nEvG(gi), ternary(strcmpi(TRIGGER,'peak'),'onset','peak'), ...
                         mod(rad2deg(mo),360)));
end

annotation(fig,'textbox',[0.02 0.945 0.96 0.045],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',9.5,'FontWeight','bold', ...
    'String',sprintf(['%s population spike-phase   |   SINGLE-LANDMARK cycle: ' ...
                      '%s to next %s, linear in time'], GENOTYPE, TRIGGER, TRIGGER));
annotation(fig,'textbox',[0.02 0.905 0.96 0.035],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',8.5,'String',sprintf(['0 at east (red = %s), counterclockwise, %d bins   |   ' ...
                      'sky dashed = median %s   |   r = events / cycle / cell (POOLED)'], TRIGGER, nB, ...
                      ternary(strcmpi(TRIGGER,'peak'),'onset','peak')));
annotation(fig,'textbox',[0.02 0.868 0.96 0.035],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',8.5,'String',['every EVENT counts equally (pooled counts / pooled cell-cycles), NOT ' ...
                      'a mean over per-cell curves   |   no occupancy compensation needed on a ' ...
                      'single-landmark phase']);

lg = gobjects(0); lgs = {};
if isgraphics(hLeg(1)), lg(end+1) = hLeg(1); lgs{end+1} = sprintf('significant (p \\leq %.2g)', pCrit); end
if isgraphics(hLeg(2)), lg(end+1) = hLeg(2); lgs{end+1} = 'not significant'; end
if ~isempty(lg)
    lgd = legend(lg, lgs, 'Box','off','FontSize',9,'Interpreter','tex');
    lgd.Units = 'normalized';  lgd.Position = [0.40 0.005 0.22 0.055];
end

stem = fullfile(outDir, sprintf('%s_polar_hist_cycle_%s_vagotomy', GENOTYPE, TRIGGER));
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');

Tout = table(rad2deg(ctrs(:)), squeeze(H(1,1,:)), squeeze(H(1,2,:)), ...
                              squeeze(H(2,1,:)), squeeze(H(2,2,:)), ...
    'VariableNames',{'phase_deg','intact_sig','intact_nonsig', ...
                     'vagotomised_sig','vagotomised_nonsig'});
writetable(Tout, [stem '.csv']);
fprintf('\nwrote %s.png / .pdf / .csv\n', stem);

%% ===================== LOCAL =====================
function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end
