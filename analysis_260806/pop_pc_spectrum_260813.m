% pop_pc_spectrum_260813.m
% -----------------------------------------------------------------------
%  PC spectrum + the first 10 PC loadings, for two cell sets:
%     ALL ACTIVE          every cell in the archive registry (active by
%                         definition: pooled events > 5)
%     logZ >= 1           the breath-modulated subset
%
%  Feature = peak-triggered mean dF/F, +/-1 s, fixed 30 Hz grid, shape-normalised
%  per cell. Same matrix the figures are built from.
%
%  The scree is shown against a TIME-SHUFFLED null: each cell's samples are
%  permuted, which preserves its amplitude distribution and destroys only the
%  temporal ordering, so the null answers "is there a waveform?" rather than
%  "is there variance?".
%
%  Runqi Zhang / 2026-08-13
% -----------------------------------------------------------------------

clear; close all; clc;
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath,'pop_analysis_260813');
L = load(fullfile(outDir,'pop_features.mat'));
T = L.T; featTau = L.featTau;
% Two triggers, same fixed +/-1 s axis, same cells. RED = inspiration onset,
% BLUE = inspiration peak, matching every other figure in this project.
TRIG = { 'onset', L.featShapeOn, [0.90 0.10 0.10]
         'peak',  L.featShape,   [0.20 0.45 0.85] };

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(260813);
wid = warning('off','all');   % pca prints a rank warning on row-z-scored data
cleanupObj = onCleanup(@() warning(wid));

nPCshow = 10;
nPerm   = 500;

SETS = { 'all_active', 'ALL ACTIVE cells',            true(height(T),1)
         'logZ1',      'cells with pooled logZ >= 1', T.logZ >= 1
         'logZ2',      'cells with pooled logZ >= 2', T.logZ >= 2
         'logZ3',      'cells with pooled logZ >= 3', T.logZ >= 3 };

for s = 1:size(SETS,1)
    tag = SETS{s,1}; ttl = SETS{s,2}; sel = SETS{s,3};

    fprintf('\n===== %s =====\n', ttl);
    RES = struct('name',{},'coeff',{},'score',{},'explained',{},'nullHi',{}, ...
                 'nSig',{},'col',{},'n',{});
    for tg = 1:size(TRIG,1)
        X  = TRIG{tg,2}(sel,:);
        Ts = T(sel,:);
        for i = 1:size(X,1)
            b = isnan(X(i,:));
            if any(b) && ~all(b)
                idx = 1:size(X,2);
                X(i,b) = interp1(idx(~b), X(i,~b), idx(b), 'linear','extrap');
            end
        end
        ok = ~all(isnan(X),2);  X = X(ok,:);  Ts = Ts(ok,:);
        X(isnan(X)) = 0;

        [coeff, score, ~, ~, explained] = pca(X);
        % PCA signs are arbitrary; fix them so the two triggers can be compared
        % without one appearing mirrored for no reason.
        for k = 1:size(coeff,2)
            [~,im] = max(abs(coeff(:,k)));
            if coeff(im,k) < 0, coeff(:,k) = -coeff(:,k); score(:,k) = -score(:,k); end
        end
        nK = min(15, size(coeff,2));
        nullExp = zeros(nPerm, nK);
        for p = 1:nPerm
            Xs = X;
            for i = 1:size(Xs,1), Xs(i,:) = Xs(i, randperm(size(Xs,2))); end
            [~,~,~,~,e] = pca(Xs);
            nullExp(p,:) = e(1:nK)';
        end
        nullHi = prctile(nullExp, 95, 1);
        nSig = find(explained(1:nK)' <= nullHi, 1) - 1;
        if isempty(nSig), nSig = nK; end
        RES(tg) = struct('name',TRIG{tg,1},'coeff',coeff,'score',score, ...
            'explained',explained,'nullHi',nullHi,'nSig',nSig,'col',TRIG{tg,3}, ...
            'n',size(X,1));
        fprintf('  %-5s n=%3d  PCs>null %d  var: %s\n', TRIG{tg,1}, size(X,1), nSig, ...
                strjoin(compose('%.1f', explained(1:nPCshow)'), ' '));
    end
    nSig = max([RES.nSig]);
    nK   = min(15, numel(RES(1).explained));

    %% ---- figure ----
    fig = figure('Color','w','Units','normalized','Position',[0.04 0.06 0.86 0.78]);
    % 4x4: spectrum and cumulative take the top row (2 tiles each), PC1-10 fill
    % tiles 5-14. A 3x4 grid has only 12 tiles and cannot hold them.
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf(['%s   |   n = %d cells   |   dF/F %.0f to %.0f ms around the trigger' ...
          '   |   RED = insp onset, BLUE = insp peak'], ttl, RES(1).n, ...
          1000*featTau(1), 1000*featTau(end)), 'FontWeight','bold','Interpreter','none');

    % PC spectrum, both triggers
    ax = nexttile(tl,1,[1 2]); hold(ax,'on');
    for tg = 1:numel(RES)
        plot(ax, 1:nK, RES(tg).explained(1:nK), '-o','LineWidth',1.5,'MarkerSize',4, ...
             'Color',RES(tg).col);
        plot(ax, 1:nK, RES(tg).nullHi, '--','LineWidth',0.9,'Color',[RES(tg).col 0.5]);
    end
    xlabel(ax,'PC'); ylabel(ax,'variance explained (%)');
    legend(ax, {sprintf('onset (%d PCs > null)',RES(1).nSig),'onset null', ...
                sprintf('peak (%d PCs > null)',RES(2).nSig),'peak null'}, 'Box','off');
    title(ax,'PC spectrum'); box(ax,'on'); xlim(ax,[0.5 nK+0.5]); grid(ax,'on');

    % cumulative, both triggers
    ax = nexttile(tl,3,[1 2]); hold(ax,'on');
    for tg = 1:numel(RES)
        plot(ax, 1:nK, cumsum(RES(tg).explained(1:nK)), '-o','LineWidth',1.5, ...
             'MarkerSize',4,'Color',RES(tg).col);
    end
    yline(ax,80,'k:','80%','LineWidth',1);
    xlabel(ax,'PC'); ylabel(ax,'cumulative variance (%)');
    ylim(ax,[0 100]); xlim(ax,[0.5 nK+0.5]);
    title(ax,'cumulative'); box(ax,'on'); grid(ax,'on');

    % first 10 PC loadings, both triggers overlaid
    for k = 1:nPCshow
        ax = nexttile(tl, 4+k); hold(ax,'on');
        for tg = 1:numel(RES)
            sig = k <= RES(tg).nSig;
            plot(ax, featTau, RES(tg).coeff(:,k), '-', ...
                 'LineWidth', ternary_num(sig,1.5,0.8), ...
                 'Color', [RES(tg).col ternary_num(sig,1,0.35)]);
        end
        xline(ax,0,'k:','LineWidth',0.8); yline(ax,0,'-','Color',[.85 .85 .85]);
        xlim(ax,[featTau(1) featTau(end)]);
        title(ax, sprintf('PC%d   %.0f%% / %.0f%%', k, ...
              RES(1).explained(k), RES(2).explained(k)), 'FontSize',9);
        if k >= 7, xlabel(ax,'time from trigger (s)'); end
        if mod(k,4)==1, ylabel(ax,'loading'); end
        box(ax,'on'); grid(ax,'on');
    end

    exportgraphics(fig, fullfile(outDir, sprintf('pc_spectrum_%s.png',tag)), ...
                   'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir, sprintf('pc_spectrum_%s.pdf',tag)), ...
                   'ContentType','vector','BackgroundColor','white');
    close(fig);
    save(fullfile(outDir, sprintf('pc_%s.mat',tag)), 'RES','featTau','Ts','nSig');
    fprintf('  saved                 : pc_spectrum_%s.png\n', tag);
end

function v = ternary_num(tf,a,b)
if tf, v = a; else, v = b; end
end
