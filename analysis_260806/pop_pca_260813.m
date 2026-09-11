% pop_pca_260813.m
% -----------------------------------------------------------------------
%  Step 2: PCA of the functional space, with the diagnostics that decide whether
%  clustering is even worth doing.
%
%  ACTIVE BLOCK = the dF/F phase profile, shape-normalised per cell.
%  The event-rate profile is NOT used as an active feature: with 6 events spread
%  over 36 bins the rate profile of a sparse cell is almost all zeros, so a
%  distance on it would rank cells by event count rather than by pattern. It is
%  carried as a supplementary variable instead.
%
%  THREE QUESTIONS, IN ORDER
%    1. How many PCs are real?      -> scree against a PHASE-SHUFFLED null, which
%                                      keeps each cell's marginal distribution and
%                                      destroys only the phase ordering.
%    2. What do the PCs mean?       -> loadings plotted back in phase.
%    3. Is this just session?       -> eta^2 of session and of group on each PC.
%       With IO_1124 = 44% of cells this is the question that decides whether any
%       clustering is biology or batch.
%
%  Runqi Zhang / 2026-08-13
% -----------------------------------------------------------------------

clear; close all; clc;
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath,'pop_analysis_260813');
L = load(fullfile(outDir,'pop_features.mat'));
T = L.T; featShape = L.featShape; featTau = L.featTau;

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(260813);

%% ---- clean: fill any gap in the waveform (edges of a short recording) ----
% NOT circular: this is a time window around the peak, so the ends do not wrap.
X = featShape;
nBad = sum(isnan(X(:)));
for i = 1:size(X,1)
    b = isnan(X(i,:));
    if any(b) && ~all(b)
        idx = 1:size(X,2);
        X(i,b) = interp1(idx(~b), X(i,~b), idx(b), 'linear', 'extrap');
    end
end
keepCell = ~all(isnan(X),2);
X = X(keepCell,:); T = T(keepCell,:);
X(isnan(X)) = 0;
fprintf('samples imputed: %d of %d (%.2f%%) | cells dropped: %d\n', ...
        nBad, numel(X), 100*nBad/max(numel(X),1), nnz(~keepCell));

%% ---- PCA ----
[coeff, score, latent, ~, explained] = pca(X);
nK = min(10, size(coeff,2));

%% ---- permutation null on the scree ----
% The null shuffles TIME SAMPLES within each cell, keeping every cell's amplitude
% distribution and destroying only the temporal ordering -- so it asks "is there
% a waveform here?", not "is there variance here?".
nPerm = 500;
nullExp = zeros(nPerm, nK);
for p = 1:nPerm
    Xs = X;
    for i = 1:size(Xs,1), Xs(i,:) = Xs(i, randperm(size(Xs,2))); end
    [~,~,~,~,e] = pca(Xs);
    nullExp(p,:) = e(1:nK)';
end
nullHi = prctile(nullExp, 95, 1);
nSig   = find(explained(1:nK)' <= nullHi, 1) - 1;
if isempty(nSig), nSig = nK; end
fprintf('PCs above the 95th pct of the phase-shuffled null: %d\n', nSig);
fprintf('  variance explained: %s\n', ...
        strjoin(compose('%.1f%%', explained(1:min(5,nK))'), ', '));

%% ---- how much of each PC is session / group? ----
sess = categorical(T.session); grp = categorical(T.group);
eta2 = @(y,g) 1 - sum(splitapply(@(v) sum((v-mean(v)).^2), y, findgroups(g))) / ...
                  sum((y-mean(y)).^2);
fprintf('\n  PC   var%%    eta2(session)  eta2(group)\n');
e2s = zeros(1,nSig); e2g = zeros(1,nSig);
for k = 1:nSig
    e2s(k) = eta2(score(:,k), sess);
    e2g(k) = eta2(score(:,k), grp);
    fprintf('  %2d  %5.1f%%      %5.3f         %5.3f\n', k, explained(k), e2s(k), e2g(k));
end

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.85 0.8]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf(['functional space: PCA of the peak-triggered dF/F (%.0f-%.0f ms) ' ...
    '  |   %d cells, %d sessions   |   %d PCs above a time-shuffled null'], ...
    1000*featTau(1), 1000*featTau(end), height(T), numel(categories(sess)), nSig), ...
    'FontWeight','bold','Interpreter','none');

% 1 scree vs null
ax = nexttile(tl); hold(ax,'on');
bar(ax, 1:nK, explained(1:nK), 'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
plot(ax, 1:nK, nullHi, 'r--o','LineWidth',1.2,'MarkerSize',4);
xlabel(ax,'PC'); ylabel(ax,'variance explained (%)');
legend(ax,{'data','shuffled null (95th pct)'},'Box','off','Location','northeast');
title(ax,'scree vs phase-shuffled null'); box(ax,'on');

% 2 loadings in phase
ax = nexttile(tl); hold(ax,'on');
cols = lines(max(nSig,1));
for k = 1:min(nSig,4)
    plot(ax, featTau, coeff(:,k), '-','LineWidth',1.6,'Color',cols(k,:));
end
xline(ax,0,'k:','LineWidth',1);
xlim(ax,[featTau(1) featTau(end)]);
xlabel(ax,'time from insp peak (s)'); ylabel(ax,'loading');
legend(ax, compose('PC%d (%.0f%%)',(1:min(nSig,4))',explained(1:min(nSig,4))), ...
       'Box','off','Location','best');
title(ax,'what the PCs mean'); box(ax,'on'); grid(ax,'on');

% 3 PC1-2 by session
ax = nexttile(tl); hold(ax,'on');
uS = categories(sess); cmapS = turbo(numel(uS));
for i = 1:numel(uS)
    m = sess==uS{i};
    scatter(ax, score(m,1), score(m,2), 14, cmapS(i,:), 'filled', ...
            'MarkerFaceAlpha',0.75);
end
xlabel(ax,'PC1'); ylabel(ax,'PC2'); title(ax,'coloured by SESSION'); box(ax,'on');
[~,iBig] = max(countcats(sess));
text(ax,0.02,0.98,sprintf('%s = %.0f%% of cells', uS{iBig}, ...
     100*max(countcats(sess))/height(T)), 'Units','normalized', ...
     'VerticalAlignment','top','FontSize',8,'Interpreter','none');

% 4 PC1-2 by group
ax = nexttile(tl); hold(ax,'on');
uG = categories(grp); cmapG = lines(numel(uG));
for i = 1:numel(uG)
    m = grp==uG{i};
    scatter(ax, score(m,1), score(m,2), 14, cmapG(i,:), 'filled','MarkerFaceAlpha',0.75);
end
legend(ax, uG, 'Box','off','Location','best');
xlabel(ax,'PC1'); ylabel(ax,'PC2'); title(ax,'coloured by GROUP'); box(ax,'on');

% 5 PC1-2 by logZ (supplementary, not an input)
ax = nexttile(tl);
scatter(ax, score(:,1), score(:,2), 16, T.logZ, 'filled','MarkerFaceAlpha',0.85);
cb = colorbar(ax); cb.Label.String = 'pooled logZ';
caxis(ax, prctile(T.logZ,[2 98]));
xlabel(ax,'PC1'); ylabel(ax,'PC2');
title(ax,'coloured by logZ (supplementary)'); box(ax,'on');

% 6 how much of each PC is session vs group
ax = nexttile(tl); hold(ax,'on');
b = bar(ax, 1:nSig, [e2s; e2g]', 'EdgeColor','none');
b(1).FaceColor = [0.85 0.33 0.10]; b(2).FaceColor = [0.00 0.45 0.74];
legend(ax,{'session','group'},'Box','off','Location','northeast');
xlabel(ax,'PC'); ylabel(ax,'\eta^2 (fraction of PC variance)');
ylim(ax,[0 1]); title(ax,'is the axis just batch?'); box(ax,'on'); grid(ax,'on');

exportgraphics(fig, fullfile(outDir,'pop_pca_diagnostic.png'), ...
               'Resolution',200,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'pop_pca_diagnostic.pdf'), ...
               'ContentType','vector','BackgroundColor','white');
save(fullfile(outDir,'pop_pca.mat'),'coeff','score','latent','explained', ...
     'nSig','e2s','e2g','X','featTau','T');
fprintf('\nsaved %s\n', fullfile(outDir,'pop_pca_diagnostic.png'));
