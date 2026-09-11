% breath_trig_phase_sorted_260806.m
% -----------------------------------------------------------------------
%  Phase heatmaps of the PHASE-LOCKED ROIs only (Rayleigh logZ >= zMin), sorted by
%  each ROI's PREFERRED PHASE.  Two panels:
%
%      dF/F on phase | z-scored dF/F on phase
%
%  Phase is piecewise (0 = inspiration onset, pi = inspiratory peak, 2pi = next
%  onset) and the cycle is tiled twice, 0..4pi, so the sorted diagonal is readable
%  across the wrap point.
%
%  Preferred phase is RY.th from the polar analysis -- the Rayleigh mean direction
%  of that cell's events -- recycled, not recomputed.  Rows are ordered by it, so a
%  cell that fires early in the cycle sits at the top.  No logZ bands are drawn:
%  every ROI here is already above the threshold, and phase is the only sort.
%
%  Input : breath_trig_heatmap_GENOTYPES_peak.mat  (Pcyc/PZcyc/roiLbl/rayZ)
%          polar_coh_vs_rayleigh_data.mat          (RY.th)
%  Output: breath_trig_phase_sorted_<TRIGGER>.png / .pdf
%
%  Runqi Zhang / 2026-08-06

clear; clc;

%% ===================== USER-EDITABLE =====================
dataFile = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806\breath_trig_heatmap_GENOTYPES_peak.mat';
rayFile  = 'D:\Ventral_surface_summary\polar_coh_vs_rayleigh_260808\polar_coh_vs_rayleigh_data.mat';
% REPOINTED 2026-08-11: the 260729 run predates Sst 0806/0807, Vglut2 0810 and the
% Vgat depth bar, so its RY.thMed covers only 389 of the ROIs the heatmap .mat holds.
zMin     = 2;              % keep ROIs with Rayleigh logZ >= this.  The output name
                           % carries the value, so both levels can coexist:
                           % breath_trig_phase_sorted_logZ1_* and _logZ2_*.
prcLim   = [0.5 99.5];     % percentile clamp, over the kept ROIs only
doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
S = load(dataFile, 'Pcyc','PZcyc','phGrid','cycRoi','ordC','roiLbl','rayZ','SpkPh','TRIGGER');
R = load(rayFile,  'RY','labels','obsOfCell');

%% ---- preferred phase per ROI, mapped through the same key as logZ ----
% Preferred phase = RY.th, the MEAN direction of the weighted event phases --
% the same quantity that feeds Rbar and therefore Z, so the sort order and the
% radius in the polar figure come from one statistic.
% NOTE the polar figure plots its DOTS at RY.thMed (weighted circular median)
% instead, and the two differ by >30 deg for 61 of the 445 cells. So rows here
% are ordered by th while the dots there sit at thMed; that is deliberate, not
% a mismatch to be "fixed" by eye.
thObs = nan(numel(R.labels),1);
ok    = R.obsOfCell >= 1 & R.obsOfCell <= numel(R.RY.th);
thObs(ok) = R.RY.th(R.obsOfCell(ok));
key = cellfun(@(s) strjoin(subsref(regexp(s,'/','split'), substruct('()',{2:4})),'/'), ...
              R.labels(:), 'uni', 0);              % MMDD/FOV/roi, no group prefix
[uk, ia] = unique(key);
M = containers.Map(uk, num2cell(thObs(ia)));

nROIall = numel(S.roiLbl);
rayTh   = nan(nROIall,1);
for i = 1:nROIall
    if isKey(M, S.roiLbl{i}), rayTh(i) = M(S.roiLbl{i}); end
end
rayTh = mod(rayTh, 2*pi);                          % into [0, 2pi)

%% ---- keep the locked ROIs, order them by preferred phase ----
keepROI = find(S.rayZ >= zMin & ~isnan(rayTh));
[thSorted, ord] = sort(rayTh(keepROI));            % smallest phase at the TOP
roiOrder = keepROI(ord);
fprintf('%d ROIs with logZ >= %g and a preferred phase\n', numel(roiOrder), zMin);
fprintf('preferred phase: %.0f to %.0f deg (median %.0f)\n', ...
        rad2deg(min(thSorted)), rad2deg(max(thSorted)), rad2deg(median(thSorted)));

roiOfRow = S.cycRoi(S.ordC);                       % ROI index per stored row
rank = nan(nROIall,1); rank(roiOrder) = 1:numel(roiOrder);
sel  = ~isnan(rank(roiOfRow));
[~, rowOrd] = sortrows([rank(roiOfRow(sel)), (1:nnz(sel)).']);   % cycles stay together
idxSel = find(sel);  idxSel = idxSel(rowOrd);

P  = double(S.Pcyc(idxSel,:));
PZ = double(S.PZcyc(idxSel,:));
ph2 = [S.phGrid(:).', S.phGrid(:).' + 2*pi];
fprintf('%d cycle rows\n', size(P,1));

cP  = prctile(P(:),  prcLim);
cPZ = prctile(PZ(:), prcLim);

%% ---- panel 3: spike rate per phase bin, z-scored WITHIN each ROI ----
% SpkPh is spikes per breath cycle in each phase bin, in percent.  z-scoring it
% across bins asks only WHERE in the cycle an ROI fires, not how much -- otherwise
% a couple of high-rate ROIs would set the colour scale for everyone.
% One row per ROI, repeated to the same height as the cycle panels so the three
% panels line up ROI for ROI.
cycPerROI = round(size(P,1) / numel(roiOrder));
Sp  = S.SpkPh(roiOrder, :);
Spz = (Sp - mean(Sp,2,'omitnan')) ./ std(Sp,0,2,'omitnan');
Spz = repelem(Spz, cycPerROI, 1);
cS  = prctile(Spz(:), prcLim);
fprintf('spike-phase rows: %d (%d cycles per ROI) | %d ROIs with no spike histogram\n', ...
        size(Spz,1), cycPerROI, nnz(all(isnan(Sp),2)));

%% ---- figure ----
fh = figure('Color','w','Position',[30 40 1560 1000]);
tl = tiledlayout(fh, 1, 3, 'TileSpacing','compact','Padding','compact');

draw(nexttile(tl,1), ph2, [P  P ], cP,  gray(256),   'dF/F',   'dF/F  (phase)');
draw(nexttile(tl,2), ph2, [PZ PZ], cPZ, parula(256), 'z (SD)', 'z-scored dF/F  (phase)');
draw(nexttile(tl,3), ph2, [Spz Spz], cS, parula(256), 'z (SD)', ...
     '%spk / cycle, z-scored per ROI');

title(tl, sprintf(['%d phase-locked ROIs (logZ >= %g), sorted by preferred phase   |   ' ...
                   'single cycles, 2 cycles shown'], numel(roiOrder), zMin), ...
      'FontWeight','bold','Interpreter','none');

%% ---- save ----
if doSave
    base = fullfile(fileparts(dataFile), ...
             sprintf('breath_trig_phase_sorted_logZ%g_%s', zMin, lower(S.TRIGGER)));
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('saved %s.{png,pdf}\n', base);
end

%% ---- local ----
function draw(ax, ph, A, clim, cmap, cbl, ttl)
imagesc(ax, ph, 1:size(A,1), A, 'AlphaData', ~isnan(A));
set(ax, 'CLim', clim, 'YDir','reverse', 'TickDir','out', 'Layer','top', 'Color',[1 1 1]);
xlim(ax, [ph(1) ph(end)]); colormap(ax, cmap); hold(ax,'on');
xline(ax, pi,   'r-', 'LineWidth', 1);
xline(ax, 3*pi, 'r-', 'LineWidth', 1);
xline(ax, 2*pi, 'w:', 'LineWidth', 1);
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], 'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'}, ...
        'TickLabelInterpreter','tex');
xlabel(ax, 'breath phase (0 = onset, \pi = peak)');
ylabel(ax, 'cycles, ROIs sorted by preferred phase');
title(ax, ttl);
cb = colorbar(ax); cb.Label.String = cbl;
end
