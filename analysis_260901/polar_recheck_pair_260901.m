% polar_recheck_pair_260901.m
% -----------------------------------------------------------------------
%  The two re-check polars side by side: same cells, same angles, two radii.
%
%      left   r = Rayleigh log Z     (occupancy-weighted, Kish n_eff)
%      right  r = coherence |C|      (multitaper TW = 4, breath band)
%      angle  = weighted circular MEAN of the event phases, linear two-landmark
%               phase (insp onset = 0 right, insp peak = pi left, counterclockwise)
%
%  IT READS THE TWO CSVs the standalone scripts write, rather than recomputing.
%  That is the point: whatever is on this figure is exactly what is on those two,
%  and the join asserts the angles agree cell for cell instead of hoping they do.
%  Run polar_recheck_260901.m and polar_recheck_coherence_260901.m first.
%
%  Runqi Zhang / 2026-09-01
% -----------------------------------------------------------------------
clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));

%% ===================== USER-EDITABLE =====================
sumRoot = 'D:\Ventral_surface_summary';
outDir  = fullfile(sumRoot,'polar_recheck_260901');
alphaZ  = 0.01;    % log Z panel: the one threshold circle and the filled/open cut
alphaC  = 0.01;    % coherence panel: same
TW      = 4;
% =========================================================

logZ_crit  = log(-log(alphaZ));
K_tap  = 2*TW - 1;
confC  = sqrt(1 - alphaC ^(1/(K_tap-1)));

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

%% ===================== JOIN THE TWO RESULT TABLES =====================
fZ = fullfile(outDir,'polar_recheck_percell.csv');
fC = fullfile(outDir,'polar_recheck_coherence_percell.csv');
assert(isfile(fZ), 'Run polar_recheck_260901.m first (missing %s)', fZ);
assert(isfile(fC), 'Run polar_recheck_coherence_260901.m first (missing %s)', fC);
Z = readtable(fZ,'TextType','string');
C = readtable(fC,'TextType','string');

[tf, loc] = ismember(C.cell_idx, Z.cell_idx);
assert(all(tf), '%d coherence cells are not in the log Z table', nnz(~tf));
Z = Z(loc, :);

dAng = rad2deg(abs(angle(exp(1i*deg2rad(C.mean_phase_deg - Z.mean_deg)))));
fprintf('joined %d cells; angle agreement max |diff| = %.3g deg\n', height(C), max(dAng));
assert(max(dAng) < 1e-6, ...
   ['the two figures disagree on the angle by up to %.3g deg -- they are supposed ' ...
    'to differ in the radius only, so re-run both scripts before plotting'], max(dAng));

th    = deg2rad(C.mean_phase_deg);
rZ    = Z.logZ;
rC    = C.coh_r;
gname = C.group;
sigZ  = rZ >= logZ_crit;
sigC  = rC >= confC;
n     = height(C);

GC   = genotype_colors_260817();
GRPS = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};
present = GRPS(ismember(GRPS, cellstr(unique(gname))));

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[1 2 38 19], ...
             'Name','polar pair: log Z and coherence');
set(fig,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);

rLimZ = ceil(max(rZ)*1.15);
axZ = polaraxes(fig,'Position',[0.015 0.08 0.365 0.70]);
axC = polaraxes(fig,'Position',[0.400 0.08 0.365 0.70]);

hLeg = draw_panel(axZ, th, rZ, sigZ, gname, present, GC, [0 rLimZ], logZ_crit, alphaZ);
       draw_panel(axC, th, rC, sigC, gname, present, GC, [0 1],     confC,     alphaC);

% Panel titles as figure annotations, not text() inside the polaraxes: at
% r = 1.3*rLim they sat on top of the header block and, on the right panel,
% on the second header line.
annotation(fig,'textbox',[0.015 0.815 0.365 0.04],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',11,'String','r = Rayleigh log Z');
annotation(fig,'textbox',[0.400 0.815 0.365 0.04],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',11, ...
    'String',sprintf('r = coherence |C|   (TW = %g)', TW));

annotation(fig,'textbox',[0.02 0.885 0.96 0.10],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','tex','FontSize',10, ...
    'String',{sprintf(['N = %d cells, SAME angle on both panels: circular mean of ' ...
                       'event phase   |   insp onset = 0 (right), insp peak = \\pi ' ...
                       '(left), counterclockwise'], n), ...
              sprintf(['filled = significant at \\alpha = %.2g   (log Z \\geq %.2f : %d ' ...
                       'cells     |C| \\geq %.3f : %d cells)      IO always open grey'], ...
                      alphaZ, logZ_crit, nnz(sigZ), confC, nnz(sigC))});

lg = present;
for q = 1:numel(lg)
    m = gname == string(lg{q});
    lg{q} = sprintf('%s (n=%d, sig %d | %d)', lg{q}, nnz(m), nnz(m & sigZ), nnz(m & sigC));
end
lgd = legend(axC, hLeg(isgraphics(hLeg)), lg(isgraphics(hLeg)), 'Box','off','FontSize',9);
lgd.Units = 'normalized';  lgd.Position = [0.795 0.55 0.19 0.24];

axT = axes(fig,'Position',[0.782 0.10 0.215 0.40]); axis(axT,'off');
Lt = strings(0,1);
Lt(end+1) = sprintf('%-8s %4s %7s %7s %6s %5s %5s', ...
                    'genotype','n','meanDg','medLogZ','med|C|','sigZ','sigC');
Lt(end+1) = string(repmat('-',1,50));
for gi = 1:numel(present)
    m = gname == string(present{gi});
    mu = angle(mean(exp(1i*th(m))));
    Lt(end+1) = sprintf('%-8s %4d %7.0f %7.2f %6.3f %5d %5d', present{gi}, nnz(m), ...
                        mod(rad2deg(mu),360), median(rZ(m)), median(rC(m)), ...
                        nnz(m & sigZ), nnz(m & sigC)); %#ok<SAGROW>
end
text(axT, 0, 1, Lt, 'FontName','Consolas','FontSize',7.5, ...
     'VerticalAlignment','top','Interpreter','none');

annotation(fig,'textbox',[0.02 0.005 0.96 0.045],'EdgeColor',[0.7 0.7 0.7], ...
    'BackgroundColor',[0.96 0.96 0.96],'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontWeight','bold','FontSize',9, ...
    'String','re-check curation: 39 cells tossed whole, 31 recordings tossed, curated events, NO GCaMP lead correction on the phase');

stem = fullfile(outDir,'polar_recheck_pair');
exportgraphics(fig, [stem '.png'], 'Resolution',300, 'BackgroundColor','white');
exportgraphics(fig, [stem '.pdf'], 'ContentType','vector','BackgroundColor','white');
fprintf('wrote %s.png / .pdf\n', stem);

%% ===================== LOCAL =====================
function h = draw_panel(ax, th, r, sig, gname, present, GC, rl, crit, a1)
hold(ax,'on');
ax.ThetaZeroLocation = 'right';
ax.ThetaDir          = 'counterclockwise';
% ONE threshold circle per panel, at the alpha the filled/open cut uses.
tt = linspace(0,2*pi,361);
polarplot(ax, tt, crit*ones(size(tt)), '--','Color',[0.15 0.15 0.15],'LineWidth',0.9);
text(ax, deg2rad(-18), crit, sprintf('  \\alpha=%.2g', a1), 'FontSize',7,'Color',[0.15 0.15 0.15]);

h = gobjects(numel(present),1);
for gi = 1:numel(present)
    g = present{gi};  col = GC.(g);
    m = find(gname == string(g));
    if isempty(m), continue; end
    if strcmp(g,'IO')
        h(gi) = polarplot(ax, th(m), max(r(m),0), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',[0.50 0.50 0.50], ...
            'LineWidth',0.7,'MarkerSize',6,'LineStyle','none');
        continue
    end
    ms = m(sig(m));  mn = m(~sig(m));
    if ~isempty(mn)
        hn = polarplot(ax, th(mn), max(r(mn),0), 'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',col, ...
            'LineWidth',0.9,'MarkerSize',6,'LineStyle','none');
    end
    if ~isempty(ms)
        h(gi) = polarplot(ax, th(ms), max(r(ms),0), 'o', ...
            'MarkerFaceColor',col,'MarkerEdgeColor','w', ...
            'LineWidth',0.5,'MarkerSize',6,'LineStyle','none');
    else
        h(gi) = hn;
    end
end
rlim(ax, rl);
ax.ThetaTick      = 0:30:330;
ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:30:330, 'uni',0);
ax.RAxisLocation  = 180;
ax.GridAlpha      = 0.15;
end
