% polar_recheck_260901.m
% -----------------------------------------------------------------------
%  Phase polar plot for the RE-CHECKED population: one dot per cell, nothing
%  else.
%
%      angle  = weighted circular MEAN of the event phases
%      r      = Rayleigh log Z   (occupancy-weighted, Kish n_eff)
%      colour = GENOTYPE  (vagotomised is not distinguished)
%      filled = significant (logZ >= log(-log(0.01)))
%      open   = not significant;  IO is always open grey
%
%  PHASE IS THE LINEAR TWO-LANDMARK PHASE: inspiration onset = 0, inspiration
%  peak = pi, next onset = 2*pi, linear in time between consecutive landmarks.
%  The axis is drawn counterclockwise with 0 to the RIGHT, so onset is at 3
%  o'clock, peak at 9 o'clock, and the inspiratory half of the cycle is the top.
%
%  NO SPREAD ARC. The 260816 figure drew a weighted-IQR arc and therefore put the
%  dot on the circular MEDIAN, so that dot and arc described the same centre.
%  With the arc gone that reason is gone too, and the dot moves to the MEAN --
%  the direction the Rayleigh Z on the radial axis is actually computed about.
%  Set dotStat = 'median' to get the 260816 convention back.
%
%  THE POPULATION IS THE RE-CHECK CURATION, not popsel. Cells thrown away whole
%  in ca_recheck_gui_260901 are dropped; within a kept cell only the kept
%  recordings contribute, and the events are that GUI's curated ones. Everything
%  else -- the occupancy weighting, Kish n_eff, the 36 phase bins, no calcium lag
%  -- is exactly what temporal_phase_cell_fig_260812 does, so the log Z here is
%  the same quantity as the log Z on every per-cell summary figure.
%
%  OCCUPANCY WEIGHTING IS MANDATORY. The breath cycle is ~7x unevenly sampled in
%  time, so an unweighted mean slides into expiration regardless of when the cell
%  fires. Each event carries w = 1/frames-in-its-phase-bin from ITS OWN recording.
%
%  Runqi Zhang / 2026-09-01
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

dotStat   = 'mean';       % 'mean' (matches the radial log Z) | 'median'
alpha1    = 0.01;         % inner dashed circle, and the filled/open cut
alpha2    = 0.001;        % reported in the CSV only; not drawn
nBins     = 36;           % phase bins for the occupancy weights
minSpikes = 2;            % a direction needs at least this many events
% =========================================================

logZ_crit  = log(-log(alpha1));      % 1.53 at alpha = 0.01
logZ_crit2 = log(-log(alpha2));      % 1.93 at alpha = 0.001

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

R = struct('cell',{},'stem',{},'group',{},'th',{},'thMed',{},'Rbar',{}, ...
           'Z',{},'logZ',{},'nEv',{},'nEff',{},'nRecKept',{},'p_an',{});
nTossCell = 0; nTossRec = 0; nThin = 0;

for k = 1:numel(d)
    L = load(fullfile(d(k).folder, d(k).name),'cellInfo','RECc');
    C = L.cellInfo;  REC = L.RECc;

    aPool = []; wPool = []; nKept = 0;
    tossedWhole = false;
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && g.cell_toss, tossedWhole = true; break; end
        if g.found && g.toss, nTossRec = nTossRec + 1; continue; end
        ev = r.spike_idx(:);                       % archive fallback
        if g.found, ev = g.spike_idx(:); end       % curated events win
        nKept = nKept + 1;

        % --- linear two-landmark phase, onset 0 / peak pi ---
        phi  = piecewise_phase_local(r.peak, r.foot, r.T);
        phiW = mod(phi, 2*pi);
        valid = ~isnan(phiW);
        if nnz(valid) < 10, continue; end

        edgesP = linspace(0, 2*pi, nBins+1);
        bidx = nan(r.T,1);
        bidx(valid) = min(discretize(phiW(valid), edgesP), nBins);
        occ = accumarray(bidx(valid), 1, [nBins 1]);   % FRAMES per bin, THIS rec

        ev = ev(ev >= 1 & ev <= r.T);
        ev = ev(~isnan(bidx(ev)));
        if isempty(ev), continue; end
        aPool = [aPool; phiW(ev)];                       %#ok<AGROW>
        wPool = [wPool; 1./max(occ(bidx(ev)),1)];        %#ok<AGROW>
    end
    if tossedWhole, nTossCell = nTossCell + 1; continue; end
    if numel(aPool) < minSpikes, nThin = nThin + 1; continue; end

    [th, Rbar, nEff] = wresultant_local(aPool, wPool);
    Z = nEff * Rbar^2;
    R(end+1) = struct('cell',C.cell, 'stem',C.stem, 'group',C.group, ...
        'th',th, 'thMed',circ_median_w_local(aPool, wPool), 'Rbar',Rbar, ...
        'Z',Z, 'logZ',log(max(Z,eps)), 'nEv',numel(aPool), 'nEff',nEff, ...
        'nRecKept',nKept, 'p_an',exp(-Z)); %#ok<SAGROW>
end

n = numel(R);
fprintf('%d cells plotted  (%d tossed whole, %d recordings tossed, %d with < %d events)\n', ...
        n, nTossCell, nTossRec, nThin, minSpikes);
assert(n > 0, 'no cells survived the curation');

switch lower(dotStat)
    case 'median', thPlot = [R.thMed];
    otherwise,     thPlot = [R.th];
end
logZ  = [R.logZ];
gname = string({R.group});
isSig = logZ >= logZ_crit;

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 20], ...
             'Name','polar: re-checked cells');
set(fig,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);
ax = polaraxes(fig,'Position',[0.05 0.10 0.66 0.72]);
hold(ax,'on');

% 0 to the RIGHT, angles increasing counterclockwise: onset at 3 o'clock, peak
% at 9 o'clock, inspiration across the top.
ax.ThetaZeroLocation = 'right';
ax.ThetaDir          = 'counterclockwise';

rMax = max(logZ);  if ~isfinite(rMax) || rMax <= 0, rMax = 1; end
rLim = ceil(rMax*1.15);

% ONE threshold circle, at the alpha the filled/open cut uses. A second circle
% at a different alpha invites reading two significance levels off one figure.
tt = linspace(0,2*pi,361);
polarplot(ax, tt, logZ_crit*ones(size(tt)), '--','Color',[0.15 0.15 0.15],'LineWidth',0.9);
text(ax, deg2rad(-18), logZ_crit, sprintf('  \\alpha=%.2g', alpha1), ...
     'FontSize',7,'Color',[0.15 0.15 0.15]);

present = GRPS(ismember(GRPS, cellstr(unique(gname))));
h = gobjects(numel(present),1);
for gi = 1:numel(present)
    g   = present{gi};
    col = GC.(g);
    m   = find(gname == string(g) & isfinite(thPlot));
    if isempty(m), continue; end

    % IO is ALWAYS open grey: it is by far the largest group, and filling it
    % buries every other genotype under one colour.
    if strcmp(g,'IO')
        h(gi) = polarplot(ax, thPlot(m), max(logZ(m),0), 'o', ...
            'MarkerFaceColor','none', 'MarkerEdgeColor',[0.50 0.50 0.50], ...
            'LineWidth',0.7,'MarkerSize',6.5,'LineStyle','none');
        continue
    end
    % everything else: filled if significant, open in its own colour if not
    ms = m(isSig(m));  mn = m(~isSig(m));
    if ~isempty(mn)
        polarplot(ax, thPlot(mn), max(logZ(mn),0), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',col, ...
            'LineWidth',0.9,'MarkerSize',6.5,'LineStyle','none');
    end
    if ~isempty(ms)
        h(gi) = polarplot(ax, thPlot(ms), max(logZ(ms),0), 'o', ...
            'MarkerFaceColor',col,'MarkerEdgeColor','w', ...
            'LineWidth',0.5,'MarkerSize',6.5,'LineStyle','none');
    else
        h(gi) = polarplot(ax, thPlot(mn(1)), max(logZ(mn(1)),0), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',col, ...
            'LineWidth',0.9,'MarkerSize',6.5,'LineStyle','none');
    end
end

rlim(ax,[0 rLim]);
ax.ThetaTick      = 0:30:330;
ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:30:330, 'uni',0);
ax.RAxisLocation  = 180;
ax.GridAlpha      = 0.15;
% Title as an annotation, not title(ax,...): on a polaraxes under software
% OpenGL the axes title renders as a filled grey band on this machine.
% The onset/peak callouts are gone with it -- at r = rLim they landed on top of
% the legend and ran off the left edge, and the header already states the
% convention.
annotation(fig,'textbox',[0.02 0.845 0.70 0.045],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',10, ...
    'String',sprintf('r = log Z      dot = circular %s of event phase', upper(dotStat)));

ok = isgraphics(h);
lg = present(ok);
for q = 1:numel(lg)
    mm = gname == string(lg{q});
    lg{q} = sprintf('%s (n=%d, sig %d)', lg{q}, nnz(mm), nnz(mm & isSig));
end
lgd = legend(ax, h(ok), lg, 'Box','off','FontSize',9);
lgd.Units = 'normalized';  lgd.Position = [0.735 0.55 0.24 0.22];

annotation(fig,'textbox',[0.02 0.90 0.96 0.09],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','tex','FontSize',10, ...
    'String',{sprintf(['N = %d cells   |   linear phase: insp onset = 0 (right), ' ...
                       'insp peak = \\pi (left), counterclockwise'], n), ...
              sprintf(['filled = significant (logZ \\geq %.2f, \\alpha = %.2g) : %d' ...
                       '      open = not significant : %d      IO always open grey'], ...
                      logZ_crit, alpha1, nnz(isSig), nnz(~isSig))});
annotation(fig,'textbox',[0.02 0.005 0.96 0.045],'EdgeColor',[0.7 0.7 0.7], ...
    'BackgroundColor',[0.96 0.96 0.96],'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontWeight','bold','FontSize',9, ...
    'String',sprintf(['re-check curation %s : %d cells tossed whole, %d recordings ' ...
                      'tossed, curated events, NO GCaMP lead correction'], ...
                      datestr(now,'yyyy-mm-dd'), nTossCell, nTossRec)); %#ok<TNOW1,DATST>

% ---- side panel ----
axT = axes(fig,'Position',[0.735 0.10 0.25 0.40]); axis(axT,'off');
Lt = strings(0,1);
Lt(end+1) = sprintf('%-8s %4s %8s %8s %5s', 'genotype','n','meanDeg','medLogZ','sig');
Lt(end+1) = string(repmat('-',1,38));
for gi = 1:numel(present)
    m = find(gname == string(present{gi}) & isfinite(thPlot));
    if isempty(m), continue; end
    mu = angle(mean(exp(1i*thPlot(m))));      % each cell counts once
    Lt(end+1) = sprintf('%-8s %4d %8.0f %8.2f %5d', present{gi}, numel(m), ...
                        mod(rad2deg(mu),360), median(logZ(m)), nnz(isSig(m))); %#ok<SAGROW>
end
text(axT, 0, 1, Lt, 'FontName','Consolas','FontSize',8, ...
     'VerticalAlignment','top','Interpreter','none');

%% ===================== SAVE =====================
stem = fullfile(outDir, sprintf('polar_recheck_%s', lower(dotStat)));
exportgraphics(fig, [stem '.png'], 'Resolution',300, 'BackgroundColor','white');
exportgraphics(fig, [stem '.pdf'], 'ContentType','vector','BackgroundColor','white');

T = table([R.cell]', string({R.stem})', string({R.group})', ...
          mod(rad2deg([R.th]),360)', mod(rad2deg([R.thMed]),360)', ...
          [R.Rbar]', [R.Z]', [R.logZ]', [R.p_an]', [R.nEv]', [R.nEff]', ...
          [R.nRecKept]', isSig', ...
    'VariableNames',{'cell_idx','stem','group','mean_deg','median_deg','Rbar', ...
                     'Z','logZ','p_analytic','n_events','n_eff','n_rec_kept','sig'});
writetable(T, fullfile(outDir,'polar_recheck_percell.csv'));

fprintf('\n%-8s %4s %8s %8s %5s\n','genotype','n','meanDeg','medLogZ','sig');
for gi = 1:numel(present)
    m = find(gname == string(present{gi}));
    mu = angle(mean(exp(1i*thPlot(m))));
    fprintf('%-8s %4d %8.0f %8.2f %5d\n', present{gi}, numel(m), ...
            mod(rad2deg(mu),360), median(logZ(m)), nnz(isSig(m)));
end
fprintf('\nwrote %s.png / .pdf and polar_recheck_percell.csv\n', stem);

%% ===================== LOCAL FUNCTIONS =====================
function phi = piecewise_phase_local(peak_idx, foot_idx, T)
%PIECEWISE_PHASE_LOCAL  Insp onset = 0, peak = pi, next onset = 2pi, linear in
%  TIME between consecutive landmarks. Strict alternation enforced; NaN outside
%  the first/last landmark. Phase accumulates past 2pi (wrap with mod).
%  Verbatim from temporal_phase_cell_fig_260812.m, where it is a local function.
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

function [th, Rbar, nEff] = wresultant_local(a, w)
% Weighted circular resultant with Kish's effective sample size. Using raw n with
% a weighted Rbar breaks the Rayleigh: a few events in a sparsely-sampled
% inspiratory bin carry huge leverage. Verbatim from polar_selected_260816.m.
S1 = sum(w);  S2 = sum(w.^2);
if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
nEff = S1^2 / max(S2, eps);
v    = sum(w(:) .* exp(1i*a(:))) / S1;
th   = angle(v);
Rbar = min(abs(v), 1);
end

function m = circ_median_w_local(a, w)
% Weighted circular median: the direction minimising the weighted sum of angular
% distances, searched on a 1 deg grid and refined. Kept so dotStat='median'
% reproduces the 260816 convention.
g = linspace(-pi, pi, 361);
cost = arrayfun(@(x) sum(w(:) .* abs(angle(exp(1i*(a(:)-x))))), g);
[~,i] = min(cost);
g2 = linspace(g(max(i-1,1)), g(min(i+1,numel(g))), 41);
cost2 = arrayfun(@(x) sum(w(:) .* abs(angle(exp(1i*(a(:)-x))))), g2);
[~,j] = min(cost2);
m = g2(j);
end
