% breath_trig_active_timephase_260806.m
% -----------------------------------------------------------------------
%  ACTIVE ROIs only, single cycles, one figure per block of ROIs.
%  Each figure has four panels:
%
%      dF/F on TIME | dF/F on PHASE | z on TIME | z on PHASE
%
%  TIME  x = seconds from the inspiratory peak (0 = peak).
%  PHASE x = piecewise breath phase: 0 = inspiration onset, pi = inspiratory peak,
%            2pi = next onset, linear inside each half.  The two halves are scaled
%            separately because inspiration is much shorter than expiration; a
%            uniform ramp between onsets would put the peak at the wrong phase.
%  So the peak sits at x = 0 in the time panels and at x = pi in the phase panels.
%
%  Row order and the white logZ lines are inherited from the upstream ordering
%  (Rayleigh logZ descending, then active).  Inactive ROIs are excluded entirely.
%  Every panel is drawn on the SAME fixed row count (roiPerFig x cycles per ROI),
%  so a cycle is the same height in every panel and a short block leaves the bottom
%  blank rather than being stretched.
%
%  Input : breath_trig_heatmap_GENOTYPES_peak.mat  (needs Pcyc/PZcyc -- rerun
%          breath_trig_heatmap_genotypes_260806.m if they are missing)
%  Output: breath_trig_active_timephase_<TRIGGER>_blockN.png / .pdf
%
%  Runqi Zhang / 2026-08-06

clear; clc;

%% ===================== USER-EDITABLE =====================
dataFile  = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806\breath_trig_heatmap_GENOTYPES_peak.mat';
roiPerFig = 200;           % ACTIVE ROIs per figure (501 active -> 3 figures,
                           % the last holding 101 ROIs and left blank below them)
prcLim    = [0.5 99.5];    % percentile clamp, computed over the active set
doSave    = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
assert(isfile(dataFile), 'Run breath_trig_heatmap_genotypes_260806.m first: %s', dataFile);
S = load(dataFile, 'Ccyc','Zcyc','Pcyc','PZcyc','phGrid','cycTier','cycRoi','ordC', ...
                   'tierName','tt','TRIGGER','zCuts');
assert(isfield(S,'Pcyc') && ~isempty(S.Pcyc), ...
       'No phase data in the .mat -- rerun breath_trig_heatmap_genotypes_260806.m.');

% top-down order, then keep ACTIVE only (every tier except the last)
nTier = numel(S.tierName);
C = flipud(S.Ccyc);  Z = flipud(S.Zcyc);  P = flipud(S.Pcyc);  PZ = flipud(S.PZcyc);
tier = flipud(S.cycTier(:));  roiOf = flipud(S.cycRoi(S.ordC));
act  = tier < nTier;
C = C(act,:); Z = Z(act,:); P = P(act,:); PZ = PZ(act,:);
tier = tier(act); roiOf = roiOf(act);

roiIdx = cumsum([1; diff(roiOf) ~= 0]);
nROI   = roiIdx(end);
nFig   = ceil(nROI / roiPerFig);
maxCyc = max(accumarray(roiOf(:), 1));
yMax   = roiPerFig * maxCyc;

cD = prctile(C(:),  prcLim);   cZ  = prctile(Z(:),  prcLim);
cP = prctile(P(:),  prcLim);   cPZ = prctile(PZ(:), prcLim);
fprintf('%d ACTIVE ROIs, %d rows -> %d figures of %d ROIs (%d rows each)\n', ...
        nROI, size(C,1), nFig, roiPerFig, yMax);
fprintf('phase coverage: %.0f%% of active cycles\n', 100*nnz(~all(isnan(P),2))/size(P,1));

% Phase panels are TILED TWICE, 0..4pi, so a feature straddling the cycle boundary
% is readable instead of being cut in half at the edge.  The second copy is the
% same data, not new cycles.
ph2 = [S.phGrid(:).', S.phGrid(:).' + 2*pi];
for f = 1:nFig
    lo = (f-1)*roiPerFig + 1;  hi = min(f*roiPerFig, nROI);
    m  = roiIdx >= lo & roiIdx <= hi;
    tc = tier(m);
    bnd = find(diff(tc) ~= 0) + 0.5;
    ctr = arrayfun(@(t) mean(find(tc==t)), 1:nTier-1);
    lbl = arrayfun(@(t) S.tierName{t}, 1:nTier-1, 'uni',0);
    keep = ~isnan(ctr); ctr = ctr(keep); lbl = lbl(keep);
    [ctr, iS] = sort(ctr); lbl = lbl(iS);

    fh = figure('Color','w','Position',[20 20 1280 1020]);
    tl = tiledlayout(fh, 1, 4, 'TileSpacing','compact','Padding','compact');

    drawT(nexttile(tl,1), S.tt, C(m,:),  cD,  bnd, ctr, lbl, gray(256),   'dF/F',   'dF/F  TIME',  S.TRIGGER, yMax);
    Pm = P(m,:);  PZm = PZ(m,:);
    drawP(nexttile(tl,2), ph2, [Pm  Pm ], cP,  bnd, [], {}, gray(256),   'dF/F',   'dF/F  PHASE (2 cycles)', yMax);
    drawT(nexttile(tl,3), S.tt, Z(m,:),  cZ,  bnd, [], {}, parula(256), 'z (SD)', 'z  TIME',                S.TRIGGER, yMax);
    drawP(nexttile(tl,4), ph2, [PZm PZm], cPZ, bnd, [], {}, parula(256), 'z (SD)', 'z  PHASE (2 cycles)',   yMax);

    title(tl, sprintf(['ACTIVE ROIs %d-%d of %d   |   single cycles, %d per ROI   |   ' ...
                       'white lines = Rayleigh logZ %s'], lo, hi, nROI, maxCyc, mat2str(S.zCuts)), ...
          'FontWeight','bold','Interpreter','none');

    fprintf('  figure %d: ROI %d-%d, %d of %d rows (%.0f%% filled)\n', f, lo, hi, nnz(m), yMax, 100*nnz(m)/yMax);
    if doSave
        base = fullfile(fileparts(dataFile), ...
                 sprintf('breath_trig_active_timephase_%s_block%d', lower(S.TRIGGER), f));
        exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
        exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        fprintf('  saved %s.{png,pdf}\n', base);
    end
end

%% ---- local ----
function drawT(ax, tt, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, TRIGGER, yMax)
common(ax, tt, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, yMax);
xline(ax, 0, 'r-', 'LineWidth', 1);
set(ax, 'XTick', [-5 0 5]);
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
end

function drawP(ax, ph, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, yMax)
common(ax, ph, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, yMax);
xline(ax, pi,   'r-', 'LineWidth', 1);        % inspiratory peak, cycle 1
xline(ax, 3*pi, 'r-', 'LineWidth', 1);        %                   cycle 2
xline(ax, 2*pi, 'w:', 'LineWidth', 1);        % where the tile repeats
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], 'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'}, ...
        'TickLabelInterpreter','tex');
xlabel(ax, 'breath phase (0 = onset, \pi = peak)');
end

function common(ax, x, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, yMax)
imagesc(ax, x, 1:size(A,1), A, 'AlphaData', ~isnan(A));
set(ax, 'CLim', clim, 'YDir','reverse', 'TickDir','out', 'Layer','top', 'Color',[1 1 1]);
xlim(ax, [x(1) x(end)]); ylim(ax, [0.5, yMax+0.5]);
colormap(ax, cmap); hold(ax,'on');
for b = bnd(:).', plot(ax, [x(1) x(end)], [b b], 'w-', 'LineWidth', 1.2); end
if isempty(ctr), set(ax, 'YTick', []);
else, set(ax, 'YTick', ctr, 'YTickLabel', lbl, 'TickLabelInterpreter','none', 'FontSize',8);
end
title(ax, ttl, 'FontSize', 10);
cb = colorbar(ax); cb.Label.String = cbl;
end
