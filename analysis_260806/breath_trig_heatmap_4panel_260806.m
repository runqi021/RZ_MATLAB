% breath_trig_heatmap_4panel_260806.m
% -----------------------------------------------------------------------
%  Heatmaps only, four panels, ACTIVE split from INACTIVE:
%
%      | active dF/F | active z | inactive dF/F | inactive z |
%
%  The Rayleigh logZ bands (white lines at logZ = 3, 2, 1) are drawn on the
%  ACTIVE panels only -- an inactive ROI has no events, so it has no phase
%  distribution and no logZ, and a band there would be meaningless.
%  Inactive ROIs keep the random order they were given upstream.
%
%  COLOUR SCALE.  The two dF/F panels share one clim and the two z panels share
%  another, both from the [0.5 99.5] percentile of the ACTIVE+INACTIVE data
%  together.  Per-panel percentiles would rescale each block to its own spread
%  and make active and inactive look equally modulated, which is the one
%  comparison this figure exists to make.
%
%  Input : breath_trig_heatmap_GENOTYPES_peak.mat  (written by
%          breath_trig_heatmap_genotypes_260806.m -- run that first)
%  Output: breath_trig_heatmap_4panel_<TRIGGER>.png / .pdf
%
%  Runqi Zhang / 2026-08-06

clear; clc;

%% ===================== USER-EDITABLE =====================
dataFile = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806\breath_trig_heatmap_GENOTYPES_peak.mat';
prcLim   = [0.5 99.5];     % shared percentile clamp
doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
assert(isfile(dataFile), ['Data file not found: %s\n' ...
       'Run breath_trig_heatmap_genotypes_260806.m first.'], dataFile);
S = load(dataFile, 'Ccyc','Zcyc','cycTier','cycRoi','ordC','tierName','tt','TRIGGER','zCuts');

Ccyc = S.Ccyc;  Zcyc = S.Zcyc;  tier = S.cycTier(:);  tt = S.tt;
roiOf   = S.cycRoi(S.ordC);           % ROI index per heatmap row, same order as Ccyc
nTier   = numel(S.tierName);
isAct   = tier <  nTier;              % every tier except the last = active
isOff   = tier == nTier;              % last tier = inactive
nActROI = numel(unique(roiOf(isAct)));
nOffROI = numel(unique(roiOf(isOff)));
maxCyc  = max(accumarray(roiOf(:), 1));

% shared colour limits, computed once over BOTH blocks
cD = prctile(Ccyc(:), prcLim);
cZ = prctile(Zcyc(:), prcLim);

% band edges inside the ACTIVE block (logZ tiers only)
tA   = tier(isAct);
bndA = find(diff(tA) ~= 0) + 0.5;
ctrA = arrayfun(@(t) mean(find(tA==t)), 1:nTier-1);
lblA = arrayfun(@(t) sprintf('%s (%d)', S.tierName{t}, nnz(tA==t)), 1:nTier-1, 'uni',0);
keep = ~isnan(ctrA);  ctrA = ctrA(keep);  lblA = lblA(keep);
[ctrA, iS] = sort(ctrA);  lblA = lblA(iS);

fprintf('active %d cycles (%d ROIs) | inactive %d cycles (%d ROIs) | up to %d cycles/ROI\n', ...
        nnz(isAct), nActROI, nnz(isOff), nOffROI, maxCyc);
fprintf('shared clim: dF/F [%.3f %.3f] | z [%.2f %.2f]\n', cD(1), cD(2), cZ(1), cZ(2));

%% ---- figure ----
fh = figure('Color','w','Position',[30 40 1620 900]);
tl = tiledlayout(fh, 1, 4, 'TileSpacing','compact','Padding','compact');

% dF/F pair first, then the z pair -- the two panels sharing a colour scale sit
% side by side so active vs inactive can be compared without crossing the figure.
draw(nexttile(tl,1), tt, Ccyc(isAct,:), cD, bndA, ctrA, lblA, gray(256),   'dF/F', ...
     sprintf('ACTIVE  dF/F  (%d ROIs)', nActROI), S.TRIGGER);
draw(nexttile(tl,2), tt, Ccyc(isOff,:), cD, [], [], {}, gray(256),   'dF/F', ...
     sprintf('INACTIVE  dF/F  (%d ROIs)', nOffROI), S.TRIGGER);
draw(nexttile(tl,3), tt, Zcyc(isAct,:), cZ, bndA, ctrA, lblA, parula(256), 'z (SD)', ...
     sprintf('ACTIVE  z-scored  (%d ROIs)', nActROI), S.TRIGGER);
draw(nexttile(tl,4), tt, Zcyc(isOff,:), cZ, [], [], {}, parula(256), 'z (SD)', ...
     sprintf('INACTIVE  z-scored  (%d ROIs)', nOffROI), S.TRIGGER);

title(tl, sprintf(['single-cycle breath-triggered dF/F   |   %g ROIs active vs %g inactive, ' ...
                   'up to %d random cycles each   |   white lines = Rayleigh logZ %s'], ...
                   nActROI, nOffROI, maxCyc, mat2str(S.zCuts)), ...
      'FontWeight','bold','Interpreter','none');

%% ---- save ----
if doSave
    base = fullfile(fileparts(dataFile), sprintf('breath_trig_heatmap_4panel_%s', lower(S.TRIGGER)));
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('saved %s.{png,pdf}\n', base);
end

%% ---- local ----
function draw(ax, tt, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, TRIGGER)
imagesc(ax, tt, 1:size(A,1), A);
set(ax, 'CLim', clim, 'YDir','normal', 'TickDir','out', 'Layer','top');
colormap(ax, cmap); hold(ax,'on');
xline(ax, 0, 'r-', 'LineWidth', 1);
for b = bnd(:).', plot(ax, [tt(1) tt(end)], [b b], 'w-', 'LineWidth', 1.2); end
if isempty(ctr)
    set(ax, 'YTick', []);
else
    set(ax, 'YTick', ctr, 'YTickLabel', lbl, 'TickLabelInterpreter','none', 'FontSize',8);
end
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
title(ax, ttl);
cb = colorbar(ax); cb.Label.String = cbl;
end
