function breath_time_modulation_tail_260727()
%% breath_time_modulation_tail_260727  Sirota-style modulation-distribution figure.
%  Descending complementary CDF of the observed score, with VERTICAL null-derived
%  significance thresholds. One panel per trigger alignment.
% -----------------------------------------------------------------------
%   x  observed modulation score M_z  (shuffle-normalised PETH modulation)
%   y  % of cells with modulation >= x        (descending, 100% -> 0%)
%   vertical dashed lines at the score thresholds where the NULL tail
%   probability equals 0.05 / 0.01 / 0.005
%
% ------------------- RELATION TO SIROTA / KARALIS ------------------------
% Same figure structure, different statistic and different null:
%
%   Sirota                              here
%   ------                              ----
%   logZ = log(n * rbar^2), Rayleigh    M_z, shuffle-normalised PETH modulation
%   analytic uniform-phase null         circular-shift PETH null
%   vertical thresholds from            vertical thresholds from the shuffled
%     the Rayleigh P_null                 P_null
%   observed complementary CDF          observed complementary CDF
%
% The x-axis is deliberately NOT called logZ. Sirota's Z is a strictly positive
% Rayleigh statistic, so log(Z) is well defined; M_z is a z-standardised quantity
% that is legitimately NEGATIVE for unmodulated cells, and log(M_z) would be
% undefined for them. It is the log of a positive test statistic, not the log of a
% z-score -- the two are not the same object.
%
% ------------------- WHY THE THRESHOLDS ARE VERTICAL ---------------------
% A horizontal line at y = 5% only says where the top 5% of the cells I happened to
% record fall -- a sample quantile of this dataset. The significance threshold is a
% SCORE, c_alpha, defined by P_null(M_z >= c_alpha) = alpha, and a score lives on
% the x-axis. With N = 34 cells the observed curve also floors at 1/34 = 2.9%, so
% horizontal 1% and 0.5% lines could never be reached at all. Faint horizontal
% guides can be switched on, but they are guides, not thresholds.
%
% ------------------- THE POOLED NULL, WITHOUT CIRCULARITY ----------------
% Each cell's shuffles are standardised by that cell's own null moments, then
% pooled across cells. Done naively every shuffled value contributes to the mean
% and SD used to standardise itself, which shrinks the tail and makes the
% thresholds too lenient. LEAVE-ONE-OUT removes that exactly and in closed form
% from the running sums:
%
%       mu_(-s)  = (S1 - T_s) / (n-1)
%       var_(-s) = [ (S2 - T_s^2) - (n-1)*mu_(-s)^2 ] / (n-2)
%       M_null_s = (T_s - mu_(-s)) / sqrt(var_(-s))
%
% This also matches how the OBSERVED score is formed: T_obs is not a shuffle, so it
% is standardised by moments that exclude it. Now every null value is too.
%
% Input : breath_time\{onset,peak}\breath_time_peth_data.mat
% Output: breath_time\modulation_tail.png / .pdf / .csv
%
% Runqi Zhang / 2026-07-28
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
TRIGGERS        = {'onset','peak'};
alphas          = [0.05 0.01 0.005];   % null tail probabilities -> vertical lines
SHOW_H_GUIDES   = true;    % faint horizontal guides at 5 / 1 / 0.5 % of cells.
                           %   Guides only -- the thresholds are the verticals.
SHOW_NULL_CURVE = false;   % also draw the null survival curve in grey. Off by
                           %   default: Sirota's panel shows only the observed CDF.
Y_LOG           = false;   % linear % is the published look; log helps the tail
doSave          = true;
% =====================================================================

base = fullfile(cfg.outRoot,'breath_time');
fig = figure('Color','w','Name','modulation distribution', ...
             'Units','centimeters','Position',[2 2 30 12]);
set(fig,'DefaultAxesFontSize',8);
col = cfg.genotype_color;
aCol = [0.15 0.15 0.15; 0.45 0.10 0.55; 0.75 0.10 0.10];
rows = {};

for t = 1:numel(TRIGGERS)
    tg = TRIGGERS{t};
    f = fullfile(base, lower(tg), 'breath_time_peth_data.mat');
    assert(isfile(f), ['Missing %s\nRun breath_time_peth_260727(''%s'') first.'], f, tg);
    S = load(f,'R','sig_exc','params');
    ti = find([S.R.tested]);
    x  = [S.R(ti).mod_exc_z]';  x = x(isfinite(x));  N = numel(x);
    assert(N > 0, 'No tested cells for %s.', tg);

    % ---- pooled null, leave-one-out standardised ----
    z0 = [];
    for c = ti
        T = double(S.R(c).null_exc);  T = T(isfinite(T));
        n = numel(T);
        if n < 12, continue; end
        S1 = sum(T);  S2 = sum(T.^2);
        mu  = (S1 - T) / (n-1);
        ss  = (S2 - T.^2) - (n-1)*mu.^2;
        sd  = sqrt(max(ss,0) / (n-2));
        ok  = sd > 0;
        z0  = [z0; (T(ok) - mu(ok)) ./ sd(ok)]; %#ok<AGROW>
    end
    z0 = sort(z0,'ascend');  M = numel(z0);
    crit = arrayfun(@(a) prctile(z0, 100*(1-a)), alphas);

    ax = axes(fig,'Position',[0.08+0.50*(t-1) 0.15 0.38 0.68]); hold(ax,'on');
    box(ax,'on'); grid(ax,'on');

    % ---- observed descending complementary CDF, in % of cells ----
    xs = sort(x,'ascend');  yo = 100*(N:-1:1)'/N;
    stairs(ax, xs, yo, '-', 'Color',col, 'LineWidth',2);
    plot(ax, xs, yo, '.', 'Color',col, 'MarkerSize',10);

    if SHOW_NULL_CURVE
        plot(ax, z0, 100*(M:-1:1)'/M, '-', 'Color',[.6 .6 .6], 'LineWidth',1.2);
    end
    if SHOW_H_GUIDES
        for a = alphas
            yline(ax, 100*a, ':', 'Color',[.75 .75 .75], 'LineWidth',0.8);
        end
    end

    % ---- VERTICAL null-derived thresholds ----
    yl = ternary(Y_LOG, [100/N*0.5 130], [0 105]);
    for k = 1:numel(alphas)
        cA = crit(k);  pct = 100*mean(x >= cA);
        plot(ax, [cA cA], yl, '--', 'Color',aCol(k,:), 'LineWidth',1.4);
        text(ax, cA, 100 - 7*(k-1), sprintf(' P_{null}<%.3g', alphas(k)), ...
             'Color',aCol(k,:), 'FontSize',7, 'Rotation',0, ...
             'HorizontalAlignment','left', 'VerticalAlignment','top');
        rows(end+1,:) = {string(tg), alphas(k), cA, pct, nnz(x>=cA), N, M}; %#ok<AGROW>
    end

    if Y_LOG, set(ax,'YScale','log'); end
    ylim(ax, yl);
    xlim(ax, [min(xs)-0.5, max(xs)+0.8]);
    xline(ax, 0, 'k-', 'LineWidth',0.8);
    xlabel(ax,'shuffle-normalised PETH modulation,  M_z');
    ylabel(ax,'cells with modulation \geq x  (%)');
    title(ax, {sprintf('%s-triggered', upper(tg)), ...
               sprintf('N = %d cells   |   %d significant (per-cell BH)', N, nnz(S.sig_exc(ti)))}, ...
          'FontSize',8);
end

sgtitle(sprintf(['%s  |  modulation distribution across cells  |  ' ...
                 'vertical dashed = score thresholds at shuffle-null P = %s'], ...
        cfg.genotype, strjoin(compose('%.3g',alphas),', ')), 'FontSize',9);

%% ---- report ----
fprintf('\n===== modulation distribution (vertical null thresholds) =====\n');
for r = 1:size(rows,1)
    fprintf('  %-6s P_null<%5.3g -> M_z >= %5.2f   %4.1f%% of cells (%d of %d)\n', ...
            rows{r,1}, rows{r,2}, rows{r,3}, rows{r,4}, rows{r,5}, rows{r,6});
end
fprintf('  pooled null: %d leave-one-out standardised shuffle values per alignment\n', rows{1,7});
fprintf(['  NOTE: these are population-level thresholds. Per-cell significance is\n' ...
         '  p_exc / q_exc in breath_time_peth_cells.csv (own null + BH across cells).\n']);

if doSave
    exportgraphics(fig, fullfile(base,'modulation_tail.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(base,'modulation_tail.pdf'),'ContentType','vector','BackgroundColor','white');
    T = cell2table(rows,'VariableNames', ...
        {'trigger','alpha_null','critical_Mz','pct_cells_above','n_cells_above','n_cells','n_null_samples'});
    writetable(T, fullfile(base,'modulation_tail.csv'));
    fprintf('\nSaved modulation_tail.png/.pdf/.csv to %s\n', base);
end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
