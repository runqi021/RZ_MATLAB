% polar_vglut2_vagotomy_260901.m
% -----------------------------------------------------------------------
%  Four polar plots, Vglut2 only: two radii x vagal-intact / vagotomised.
%
%      columns  vagal intact  |  vagotomised
%      row 1    r = coherence |C| in the breath band (multitaper TW = 4)
%      row 2    r = Rayleigh log Z (occupancy-weighted, Kish n_eff)
%      angle    weighted circular MEAN of the event phases, both rows
%      filled   OUTSIDE that row's dashed circle;  open = inside it
%
%  IO IS NOT HERE. IO is its own group label in the registry -- the IO sites of
%  the Vglut2 sessions are recorded as 'IO', not as 'Vglut2' -- so filtering on
%  the group name excludes them by construction, including the ones that live
%  under Vglut2\1124\IO\.
%
%  Linear two-landmark phase throughout: insp onset = 0 at the right, insp peak =
%  pi at the left, counterclockwise.
%
%  ONE TEST PER PANEL, AND THE CIRCLE IS IT. Each row carries its own threshold
%  at the same alpha = 0.01, and a dot is filled exactly when it lies outside
%  that row's circle:
%     coherence  confC   = sqrt(1 - alpha^(1/(K-1))), K = 2*TW-1 = 7  ->  0.732
%     log Z      logZcrit = log(-log(alpha))                          ->  1.53
%
%  The first version of this figure drew the confC circle but filled the dots by
%  the PSTH permutation instead. Those are different tests, so dots inside the
%  circle came out filled and dots outside came out open -- a figure that
%  contradicts itself. Whatever threshold is drawn has to be the threshold that
%  is marked.
%
%  CAVEAT ON confC ITSELF, measured on these data: the circular-shift null for
%  the coherence of a sparse spike train sits at |C| ~ 0.38 and its 99th
%  percentile ranges 0.59-0.97 ACROSS CELLS, so confC is not a per-cell-valid
%  threshold here -- it is an analytic level for two Gaussian processes. The
%  permutation p for every cell is in polar_recheck_perm_percell.csv if a
%  per-cell test is wanted instead.
%
%  Reads polar_recheck_perm_percell.csv rather than recomputing, so every dot is
%  the same number as on the six-panel figure.
%
%  Runqi Zhang / 2026-09-01
% -----------------------------------------------------------------------
clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));

%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
outDir    = fullfile(sumRoot,'polar_recheck_260901');
sigCsv    = fullfile(outDir,'polar_recheck_perm_percell.csv');

GENOTYPE  = 'Vglut2';
VAG_DATES = {'0824'};
TW        = 4;
alphaCoh  = 0.01;      % the dashed confC reference circle
pCrit     = 0.01;      % filled/open, permutation
% =========================================================

confC     = sqrt(1 - alphaCoh^(1/(2*TW-1-1)));
logZ_crit = log(-log(alphaCoh));
COL       = [0.10 0.65 0.20];      % Vglut2 green

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

%% ===================== LOAD =====================
assert(isfile(sigCsv), 'Run polar_recheck_perm_260901.m first (missing %s)', sigCsv);
T = readtable(sigCsv,'TextType','string');
T = T(T.group == GENOTYPE, :);          % IO carries its own label, so it is gone
isVag = false(height(T),1);
for v = 1:numel(VAG_DATES)
    isVag = isVag | contains(T.stem, "_" + VAG_DATES{v} + "_");
end
th  = deg2rad(T.mean_deg);
rC  = T.coh_r;
rZ  = T.logZ;
% THE CIRCLE IS THE RULE, per row. Filling by the PSTH permutation while drawing
% a confC circle put filled dots inside the circle and open dots outside it -- a
% figure that contradicts itself. Whatever threshold is drawn has to be the one
% marked.
sigC = rC >= confC;
sigZ = rZ >= logZ_crit;
fprintf('%s: %d cells (%d intact, %d vagotomised)\n', GENOTYPE, height(T), ...
        nnz(~isVag), nnz(isVag));
fprintf('  |C|  >= %.3f : %d intact / %d vagotomised\n', confC, ...
        nnz(sigC & ~isVag), nnz(sigC & isVag));
fprintf('  logZ >= %.3f : %d intact / %d vagotomised\n', logZ_crit, ...
        nnz(sigZ & ~isVag), nnz(sigZ & isVag));
fprintf('  (permutation p <= %.2g would give %d / %d)\n', pCrit, ...
        nnz(logical(T.sig_perm) & ~isVag), nnz(logical(T.sig_perm) & isVag));

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[1 1 25 24], ...
             'Name','Vglut2 polar: coherence and log Z, intact vs vagotomised');
set(fig,'DefaultAxesFontSize',9);

rLimZ = ceil(max(rZ)*1.15);
ROW = { rC, sigC, [0 1],     confC,     'r = coherence |C|  (TW = 4)'
        rZ, sigZ, [0 rLimZ], logZ_crit, 'r = Rayleigh log Z' };
GRP = { 'vagal intact', ~isVag; 'vagotomised', isVag };
y0  = [0.505 0.055];

for ri = 1:2
    rr   = ROW{ri,1};  sg = ROW{ri,2};  rl = ROW{ri,3};  crit = ROW{ri,4};
    for gi = 1:2
        m  = GRP{gi,2};
        ax = polaraxes(fig,'Position',[0.055+(gi-1)*0.470, y0(ri), 0.355, 0.335]);
        hold(ax,'on');
        ax.ThetaZeroLocation = 'right';
        ax.ThetaDir          = 'counterclockwise';

        tt = linspace(0,2*pi,361);
        polarplot(ax, tt, crit*ones(size(tt)), '--','Color',[0.15 0.15 0.15],'LineWidth',0.9);
        text(ax, deg2rad(-20), crit, sprintf('  \\alpha=%.2g', alphaCoh), ...
             'FontSize',7,'Color',[0.15 0.15 0.15]);

        mn = m & ~sg;  ms = m & sg;
        if any(mn)
            polarplot(ax, th(mn), max(rr(mn),0), 'o', 'MarkerFaceColor','none', ...
                'MarkerEdgeColor',COL, 'LineWidth',1.0,'MarkerSize',6.5,'LineStyle','none');
        end
        if any(ms)
            polarplot(ax, th(ms), max(rr(ms),0), 'o', 'MarkerFaceColor',COL, ...
                'MarkerEdgeColor','w', 'LineWidth',0.5,'MarkerSize',6.5,'LineStyle','none');
        end
        rlim(ax, rl);
        ax.ThetaTick      = 0:30:330;
        ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:30:330, 'uni',0);
        ax.RAxisLocation  = 180;
        ax.GridAlpha      = 0.15;
        ax.FontSize       = 7.5;

        annotation(fig,'textbox',[0.055+(gi-1)*0.470, y0(ri)+0.345, 0.355, 0.030], ...
            'EdgeColor','none','HorizontalAlignment','center', ...
            'FontWeight','bold','FontSize',9.5,'Interpreter','none', ...
            'String',sprintf('%s %s  |  %d cells, %d outside', GENOTYPE, GRP{gi,1}, ...
                             nnz(m), nnz(ms)));
        mu = angle(mean(exp(1i*th(m))));
        annotation(fig,'textbox',[0.055+(gi-1)*0.470, y0(ri)+0.315, 0.355, 0.028], ...
            'EdgeColor','none','HorizontalAlignment','center','FontSize',8, ...
            'Interpreter','none', ...
            'String',sprintf('mean angle %.0f deg   median r %.3f', ...
                             mod(rad2deg(mu),360), median(rr(m))));
    end
    annotation(fig,'textbox',[0.055, y0(ri)+0.378, 0.825, 0.028],'EdgeColor','none', ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',10.5, ...
        'Interpreter','none','String',ROW{ri,5});
end

% Three short lines: a single long one runs off both edges at this width.
annotation(fig,'textbox',[0.02 0.955 0.96 0.040],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',9,'String',sprintf(['%s only, IO excluded   |   angle = circular MEAN of ' ...
                                   'event phase, both rows'], GENOTYPE));
annotation(fig,'textbox',[0.02 0.928 0.96 0.030],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',8.5,'String', ...
    'linear two-landmark phase: insp onset = 0 (right), insp peak = \pi (left), counterclockwise');
annotation(fig,'textbox',[0.02 0.902 0.96 0.030],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',8.5,'String',sprintf(['filled = OUTSIDE that row''s dashed circle   |   ' ...
    'confC(\\alpha=%.2g) = %.3f,  log Z crit = %.2f'], alphaCoh, confC, logZ_crit));

stem = fullfile(outDir, sprintf('%s_polar_coh_logZ_vagotomy', GENOTYPE));
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');

fprintf('\n%-14s %6s %8s %8s %8s %8s\n', ...
        'group','cells','med|C|','>confC','medLogZ','>crit');
for gi = 1:2
    m = GRP{gi,2};
    fprintf('%-14s %6d %8.3f %8d %8.2f %8d\n', GRP{gi,1}, nnz(m), ...
            median(rC(m)), nnz(sigC & m), median(rZ(m)), nnz(sigZ & m));
end
fprintf('\nwrote %s.png / .pdf\n', stem);
