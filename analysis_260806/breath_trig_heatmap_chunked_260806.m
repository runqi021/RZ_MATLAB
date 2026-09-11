% breath_trig_heatmap_chunked_260806.m
% -----------------------------------------------------------------------
%  Single-cycle heatmaps, TWO separate figures so each panel is full height:
%
%      figure 1   dF/F   block 1 | block 2 | block 3
%      figure 2   z      block 1 | block 2 | block 3
%
%  ROI order is the one built upstream -- Rayleigh logZ descending, then active,
%  then inactive in random order -- read TOP DOWN, so block 1 starts with the most
%  strongly phase-locked ROIs and the last block is all inactive.  White lines mark
%  logZ = 3, 2, 1 and the active|inactive boundary wherever they fall in a block.
%
%  EQUAL ROW HEIGHT.  Every panel is given the SAME y range,
%      roiPerCol x maxCycPerROI rows,
%  regardless of how many rows it actually holds, so one breath cycle occupies the
%  same number of pixels in every panel and blocks are visually comparable.  A block
%  with fewer rows (short recordings contribute fewer cycles, or the last block is a
%  remainder) simply leaves the bottom BLANK instead of being stretched to fit.
%
%  Both dF/F panels share one clim, all z panels share another, taken over the full
%  dataset.
%
%  Input : breath_trig_heatmap_GENOTYPES_peak.mat
%  Output: breath_trig_heatmap_chunked_dff_<TRIGGER>.png/.pdf
%          breath_trig_heatmap_chunked_z_<TRIGGER>.png/.pdf
%
%  Runqi Zhang / 2026-08-06

clear; clc;

%% ===================== USER-EDITABLE =====================
dataFile  = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806\breath_trig_heatmap_GENOTYPES_peak.mat';
roiPerCol = 300;           % ROIs per panel.  2601 / 300 = 8 full panels + 1 holding
                           % 201 ROIs, which is drawn at the same row height as the
                           % rest and simply left blank at the bottom.
panelWidth = 200;          % px per panel; every panel spans the FULL figure height
                           % (one row, no wrapping) so blocks stay side-by-side
                           % comparable and each is as tall as the figure allows
prcLim    = [0.5 99.5];    % shared percentile clamp
doSave    = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
assert(isfile(dataFile), 'Run breath_trig_heatmap_genotypes_260806.m first: %s', dataFile);
S = load(dataFile, 'Ccyc','Zcyc','cycTier','cycRoi','ordC','tierName','tt','TRIGGER','zCuts');

% flip to TOP-DOWN: upstream rows are bottom-up (row 1 = bottom of the image)
C = flipud(S.Ccyc);  Z = flipud(S.Zcyc);
tier = flipud(S.cycTier(:));  roiOf = flipud(S.cycRoi(S.ordC));
roiIdx = cumsum([1; diff(roiOf) ~= 0]);        % 1..nROI counted from the top
nROI   = roiIdx(end);
nCol   = ceil(nROI / roiPerCol);
maxCyc = max(accumarray(roiOf(:), 1));         % cycles per ROI (nCycShow upstream)
yMax   = roiPerCol * maxCyc;                   % FIXED rows per panel -> equal spacing

cD = prctile(C(:), prcLim);
cZ = prctile(Z(:), prcLim);
fprintf('%d ROIs, %d rows -> %d blocks of %d ROIs | %d cycles/ROI -> %d rows per panel\n', ...
        nROI, size(C,1), nCol, roiPerCol, maxCyc, yMax);
fprintf('shared clim: dF/F [%.3f %.3f] | z [%.2f %.2f]\n', cD(1), cD(2), cZ(1), cZ(2));

%% ---- two figures ----
MEAS = { 'dff', C, cD, gray(256),   'dF/F',   'single-cycle dF/F'
         'z',   Z, cZ, parula(256), 'z (SD)', 'single-cycle z-scored dF/F' };

for k = 1:size(MEAS,1)
    [tag, A, clim, cmap, cbl, ttl] = MEAS{k,:};
    fh = figure('Color','w','Position',[10 20 min(3200, panelWidth*nCol + 220), 1040]);
    tl = tiledlayout(fh, 1, nCol, 'TileSpacing','compact','Padding','compact');

    for c = 1:nCol
        lo = (c-1)*roiPerCol + 1;  hi = min(c*roiPerCol, nROI);
        m  = roiIdx >= lo & roiIdx <= hi;
        tc = tier(m);
        bnd = find(diff(tc) ~= 0) + 0.5;
        ctr = arrayfun(@(t) mean(find(tc==t)), 1:numel(S.tierName));
        lbl = arrayfun(@(t) S.tierName{t}, 1:numel(S.tierName), 'uni',0);
        keep = ~isnan(ctr); ctr = ctr(keep); lbl = lbl(keep);
        [ctr, iS] = sort(ctr); lbl = lbl(iS);

        draw(nexttile(tl, c), S.tt, A(m,:), clim, bnd, ctr, lbl, cmap, cbl, ...
             sprintf('ROI %d-%d', lo, hi), S.TRIGGER, yMax, c == nCol);
        if k == 1
            fprintf('  block %d: ROI %4d-%4d, %5d of %d rows (%.0f%% filled)\n', ...
                    c, lo, hi, nnz(m), yMax, 100*nnz(m)/yMax);
        end
    end

    title(tl, sprintf(['%s   |   %d ROIs in blocks of %d, ordered by Rayleigh logZ ' ...
                       'then active then inactive   |   white lines = logZ %s and active|inactive'], ...
                       ttl, nROI, roiPerCol, mat2str(S.zCuts)), ...
          'FontWeight','bold','Interpreter','none');

    if doSave
        base = fullfile(fileparts(dataFile), ...
                        sprintf('breath_trig_heatmap_chunked_%s_%s', tag, lower(S.TRIGGER)));
        exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
        exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        fprintf('saved %s.{png,pdf}\n', base);
    end
end

%% ---- local ----
function draw(ax, tt, A, clim, bnd, ctr, lbl, cmap, cbl, ttl, TRIGGER, yMax, showCB)
imagesc(ax, tt, 1:size(A,1), A);
set(ax, 'CLim', clim, 'YDir','reverse', 'TickDir','out', 'Layer','top');   % row 1 on top
ylim(ax, [0.5, yMax+0.5]);      % identical in every panel; short blocks leave blank
colormap(ax, cmap); hold(ax,'on');
xline(ax, 0, 'r-', 'LineWidth', 1);
for b = bnd(:).', plot(ax, [tt(1) tt(end)], [b b], 'w-', 'LineWidth', 1.2); end
if isempty(ctr), set(ax, 'YTick', []);
else, set(ax, 'YTick', ctr, 'YTickLabel', lbl, 'TickLabelInterpreter','none', 'FontSize',8);
end
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
title(ax, ttl, 'FontSize', 9);
set(ax, 'XTick', [-5 0 5]);
if showCB       % one colorbar for the figure -- they all share a clim anyway
    cb = colorbar(ax); cb.Label.String = cbl;
end
end
