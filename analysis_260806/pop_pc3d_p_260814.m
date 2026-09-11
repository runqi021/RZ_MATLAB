% pop_pc3d_p_260814.m
% -----------------------------------------------------------------------
%  PC1-PC2-PC3 in 3D coloured by group, selecting cells on the PSTH PERMUTATION
%  p rather than on logZ, at p < 0.05, 0.01 and 0.001.
%
%  EACH TRIGGER USES ITS OWN p AND ITS OWN WAVEFORM:
%     left  panel : cells with pOnset < cut, PCA of the ONSET-triggered dF/F
%     right panel : cells with pPeak  < cut, PCA of the PEAK-triggered dF/F
%  so the two panels are generally DIFFERENT cell sets. That is the point -- a
%  cell can be locked to the onset and not to the peak, or the reverse, and
%  forcing one shared set would hide exactly that.
%
%  The p is the max-statistic permutation p over the central 1 IBI, computed with
%  P.nShuffle = 1200. Its floor is 1/1201 = 8.3e-4, so p < 0.001 is reachable but
%  only by cells sitting exactly at the floor -- that cut is "as significant as
%  1200 shuffles can resolve", not a finer measurement.
%
%  Saved as .png/.pdf/.fig with rotate3d already on.
%
%  Runqi Zhang / 2026-08-14
% -----------------------------------------------------------------------

clear; close all; clc;
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath,'pop_analysis_260813');
L = load(fullfile(outDir,'pop_features.mat'));
T = L.T; featTau = L.featTau;
assert(ismember('pOnset', T.Properties.VariableNames), ...
    'pop_features.mat has no pOnset -- rerun pop_features_260813.m with P.nShuffle > 0');

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
wid = warning('off','all'); cleanupObj = onCleanup(@() warning(wid));

PCUTS    = [0.05 0.01 0.001];
keepOpen = true;
TRIGNAME = {'insp ONSET','insp PEAK'};
TRIGCOL  = {[0.90 0.10 0.10],[0.20 0.45 0.85]};
TRIGP    = {T.pOnset, T.pPeak};
TRIGX    = {L.featShapeOn, L.featShape};

fprintf('permutation p available for %d cells (%d NaN)\n', ...
        nnz(isfinite(T.pOnset)), nnz(~isfinite(T.pOnset)));

for ci = 1:numel(PCUTS)
    pcut = PCUTS(ci);
    fprintf('\n===== p < %g =====\n', pcut);

    fig = figure('Color','w','Units','normalized','Position',[0.05 0.08 0.86 0.72], ...
                 'Name',sprintf('PC1-3, permutation p < %g',pcut));
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    nUsed = zeros(1,2);
    for tg = 1:2
        sel = TRIGP{tg} < pcut;
        X   = TRIGX{tg}(sel,:);
        Ts  = T(sel,:);
        for i = 1:size(X,1)
            b = isnan(X(i,:));
            if any(b) && ~all(b)
                idx = 1:size(X,2);
                X(i,b) = interp1(idx(~b), X(i,~b), idx(b), 'linear','extrap');
            end
        end
        ok = ~all(isnan(X),2);  X = X(ok,:);  Ts = Ts(ok,:);
        X(isnan(X)) = 0;
        nUsed(tg) = size(X,1);

        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        if size(X,1) < 4
            axis(ax,'off');
            text(ax,0.5,0.5,sprintf('%s: only %d cells at p < %g', ...
                 TRIGNAME{tg}, size(X,1), pcut),'Horizontal','center');
            fprintf('  %-11s n=%3d  -- too few to factor\n', TRIGNAME{tg}, size(X,1));
            continue;
        end

        [coeff, score, ~, ~, explained] = pca(X);
        for k = 1:size(coeff,2)
            [~,im] = max(abs(coeff(:,k)));
            if coeff(im,k) < 0, coeff(:,k) = -coeff(:,k); score(:,k) = -score(:,k); end
        end
        fprintf('  %-11s n=%3d  PC1-3 = %.1f / %.1f / %.1f %%  (sum %.1f)\n', ...
                TRIGNAME{tg}, size(X,1), explained(1), explained(2), explained(3), ...
                sum(explained(1:3)));

        grp = categorical(Ts.group);
        uG  = categories(grp);  cnts = countcats(grp);
        cmap = lines(numel(uG));
        for i = 1:numel(uG)
            m = grp == uG{i};
            if ~any(m), continue; end
            scatter3(ax, score(m,1), score(m,2), score(m,3), 46, cmap(i,:), 'filled', ...
                'MarkerFaceAlpha',0.85,'MarkerEdgeColor',[0.25 0.25 0.25], ...
                'LineWidth',0.3,'DisplayName',sprintf('%s (%d)',uG{i},cnts(i)));
        end
        zl = min(score(:,3)) - 0.1*max(range(score(:,3)),eps);
        for j = 1:size(score,1)
            plot3(ax,[score(j,1) score(j,1)],[score(j,2) score(j,2)],[zl score(j,3)], ...
                  '-','Color',[0.75 0.75 0.75 0.5],'LineWidth',0.4,'HandleVisibility','off');
        end
        xlabel(ax, sprintf('PC1  (%.1f%%)', explained(1)));
        ylabel(ax, sprintf('PC2  (%.1f%%)', explained(2)));
        zlabel(ax, sprintf('PC3  (%.1f%%)', explained(3)));
        title(ax, sprintf('%s-triggered   n=%d', TRIGNAME{tg}, size(X,1)), ...
              'Color', TRIGCOL{tg});
        view(ax,-37.5,22); axis(ax,'vis3d');
        legend(ax,'Location','northeastoutside','Box','off','FontSize',8);
        for i = 1:numel(uG), fprintf('       %-8s %3d\n', uG{i}, cnts(i)); end
    end

    title(tl, sprintf(['functional space PC1-PC3   |   permutation p < %g   |   ' ...
          'onset n=%d, peak n=%d   |   dF/F %.0f to %.0f ms   |   coloured by group'], ...
          pcut, nUsed(1), nUsed(2), 1000*featTau(1), 1000*featTau(end)), ...
          'FontWeight','bold','Interpreter','none');
    rotate3d(fig,'on');

    stem = sprintf('pc3d_p%g_by_group', pcut);
    exportgraphics(fig, fullfile(outDir,[stem '.png']),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
    savefig(fig, fullfile(outDir,[stem '.fig']));
    fprintf('  saved %s .png/.pdf/.fig\n', stem);
    if ~keepOpen, close(fig); end
end
fprintf('\nAll saved to %s (open a .fig and drag to rotate)\n', outDir);
