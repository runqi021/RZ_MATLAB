% pop_pc3d_260814.m
% -----------------------------------------------------------------------
%  PC1-PC2-PC3 in 3D, coloured by genotype/group, with the ONSET-triggered and
%  PEAK-triggered spaces side by side, for three logZ cuts (>=1, >=2, >=3).
%
%  Recomputed here rather than read from pc_*.mat so the two triggers are
%  guaranteed to be built from the SAME rows in the SAME order -- each PCA drops
%  all-NaN rows independently, and if those ever differed the group labels would
%  silently misalign with one of the score matrices.
%
%  PCA signs are arbitrary, so each component's sign is fixed to make its
%  largest-magnitude sample positive; otherwise an axis can flip between the two
%  panels and the comparison is meaningless.
%
%  MOVABLE: each figure is saved as .fig with rotate3d already enabled, so
%  opening it lets you drag the cloud around. The two panels rotate independently
%  (a linkprop tie does not survive being written to a .fig).
%
%  Runqi Zhang / 2026-08-14
% -----------------------------------------------------------------------

clear; close all; clc;
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath,'pop_analysis_260813');
L = load(fullfile(outDir,'pop_features.mat'));
T = L.T; featTau = L.featTau;

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
wid = warning('off','all'); cleanupObj = onCleanup(@() warning(wid));

CUTS      = [1 2 3];
keepOpen  = true;      % leave the figures open so they can be rotated now
TRIGNAME  = {'insp ONSET','insp PEAK'};
TRIGCOL   = {[0.90 0.10 0.10], [0.20 0.45 0.85]};

for ci = 1:numel(CUTS)
    logZcut = CUTS(ci);
    sel = T.logZ >= logZcut;
    Ts  = T(sel,:);
    TRIGX = { L.featShapeOn(sel,:), L.featShape(sel,:) };
    fprintf('\n===== logZ >= %g : %d cells =====\n', logZcut, height(Ts));

    % rows usable in BOTH triggers, so the panels share one cell set
    keepAll = true(height(Ts),1);
    for tg = 1:2, keepAll = keepAll & ~all(isnan(TRIGX{tg}),2); end

    SC = cell(1,2); EX = cell(1,2);
    for tg = 1:2
        X = TRIGX{tg}(keepAll,:);
        for i = 1:size(X,1)
            b = isnan(X(i,:));
            if any(b) && ~all(b)
                idx = 1:size(X,2);
                X(i,b) = interp1(idx(~b), X(i,~b), idx(b), 'linear','extrap');
            end
        end
        X(isnan(X)) = 0;
        [coeff, score, ~, ~, explained] = pca(X);
        for k = 1:size(coeff,2)
            [~,im] = max(abs(coeff(:,k)));
            if coeff(im,k) < 0, coeff(:,k) = -coeff(:,k); score(:,k) = -score(:,k); end
        end
        SC{tg} = score; EX{tg} = explained;
        fprintf('  %-11s PC1-3 = %.1f / %.1f / %.1f %%  (sum %.1f)\n', TRIGNAME{tg}, ...
                explained(1), explained(2), explained(3), sum(explained(1:3)));
    end
    Ts  = Ts(keepAll,:);
    grp = categorical(Ts.group);
    uG  = categories(grp);
    cmap = lines(numel(uG));
    cnts = countcats(grp);
    fprintf('  cells plotted: %d\n', height(Ts));
    for i = 1:numel(uG), fprintf('     %-8s %3d\n', uG{i}, cnts(i)); end

    %% ---- figure ----
    fig = figure('Color','w','Units','normalized','Position',[0.05 0.08 0.86 0.72], ...
                 'Name',sprintf('PC1-3, logZ >= %g',logZcut));
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf(['functional space PC1-PC3   |   logZ >= %g, n = %d cells   |   ' ...
          'dF/F %.0f to %.0f ms around the trigger   |   coloured by group'], ...
          logZcut, height(Ts), 1000*featTau(1), 1000*featTau(end)), ...
          'FontWeight','bold','Interpreter','none');

    axAll = gobjects(1,2);
    for tg = 1:2
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        axAll(tg) = ax;
        sc = SC{tg}; ex = EX{tg};
        for i = 1:numel(uG)
            m = grp == uG{i};
            if ~any(m), continue; end
            scatter3(ax, sc(m,1), sc(m,2), sc(m,3), 46, cmap(i,:), 'filled', ...
                     'MarkerFaceAlpha',0.85,'MarkerEdgeColor',[0.25 0.25 0.25], ...
                     'LineWidth',0.3,'DisplayName',sprintf('%s (%d)',uG{i},nnz(m)));
        end
        % drop lines to the PC1-PC2 floor so depth is readable on a flat page
        zl = min(sc(:,3)) - 0.1*max(range(sc(:,3)),eps);
        for j = 1:size(sc,1)
            plot3(ax,[sc(j,1) sc(j,1)],[sc(j,2) sc(j,2)],[zl sc(j,3)],'-', ...
                  'Color',[0.75 0.75 0.75 0.5],'LineWidth',0.4,'HandleVisibility','off');
        end
        xlabel(ax, sprintf('PC1  (%.1f%%)', ex(1)));
        ylabel(ax, sprintf('PC2  (%.1f%%)', ex(2)));
        zlabel(ax, sprintf('PC3  (%.1f%%)', ex(3)));
        title(ax, sprintf('%s-triggered', TRIGNAME{tg}), 'Color', TRIGCOL{tg});
        view(ax, -37.5, 22); axis(ax,'vis3d');
        if tg == 1, legend(ax,'Location','northeastoutside','Box','off'); end
    end
    rotate3d(fig,'on');            % draggable straight away, and saved into the .fig

    stem = sprintf('pc3d_logZ%g_by_group', logZcut);
    exportgraphics(fig, fullfile(outDir,[stem '.png']),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
    savefig(fig, fullfile(outDir,[stem '.fig']));
    fprintf('  saved %s .png/.pdf/.fig\n', stem);
    if ~keepOpen, close(fig); end
end

fprintf(['\nAll saved to %s\n' ...
         'Open a .fig and drag to rotate (rotate3d is already on).\n'], outDir);
