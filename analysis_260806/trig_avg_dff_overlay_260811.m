% trig_avg_dff_overlay_260811.m
% -----------------------------------------------------------------------
%  The SIGNIFICANT cells' triggered-average dF/F, OVERLAID, no population mean.
%
%      row 1   dF/F        as computed in trig_avg_dff_260811
%      row 2   z-scored    each cell's own trace, (x - mean)/std over the
%                          plotted window
%
%      col 1  ONSET-triggered   col 2  PEAK-triggered   col 3  SELF-triggered
%
%  WHY BOTH ROWS. In dF/F a few large-amplitude cells dominate the ink, so the
%  picture is partly about which cells are bright. Z-scoring each trace to its
%  own window removes amplitude entirely and leaves only SHAPE and TIMING, so
%  every cell contributes equally regardless of how strongly it fires. A feature
%  that survives in both rows is a feature of the population; one that appears
%  only in dF/F belongs to a handful of loud cells; one that appears only in z
%  is shared but weak.
%
%  NO MEAN IS DRAWN, deliberately -- the population average and its across-animal
%  SEM live in trig_avg_dff_260811. This figure is for reading the SPREAD: how
%  consistent the cells are with each other, and whether the average describes
%  them or hides a split.
%
%  Traces are loaded from trig_avg_dff.mat, NOT recomputed, so these are exactly
%  the curves the population figure averaged. Baseline conventions therefore carry
%  over unchanged: breath columns have the window mean removed and are in CYCLES,
%  the self column has a pre-event baseline removed and is in SECONDS.
%
%  Input : event_latency_260811\trig_avg_dff.mat
%  Output: trig_avg_dff_sig_overlay.png / .pdf
%
%  Runqi Zhang / 2026-08-11

clear; clc; close all;

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath,'event_latency_260811','trig_avg_dff.mat');
outDir   = fullfile(rootPath,'event_latency_260811');
lineCol  = [0.35 0.35 0.35];
lineAlph = 0.30;           % low enough that overlap reads as density
lineW    = 0.7;
prcYLim  = [0.5 99.5];     % y-limits from the traces themselves, not from outliers
shareY   = true;           % one y-scale per column across ALL SIX figures, so the
                           % criteria can be compared by eye instead of each
                           % figure silently rescaling to its own selection
doSave   = true;

% SIX SELECTIONS, from TWO INDEPENDENT TESTS.
%   p      = the time-domain permutation test (this analysis): is the latency
%            distribution concentrated?
%   logZ   = the occupancy-weighted Rayleigh from the phase summary: is the
%            event PHASE distribution concentrated?
% They are not nested and they disagree on real cells -- the Rayleigh reaches
% cells the permutation never tested (logZ > 1 selects 108 cells, only 78 of
% which cleared the 20-event gate). Running both is the point: a feature that
% survives every panel does not depend on which test defined "locked".
CRIT = { 'p001',  'permutation p < 0.001', @(a) a.p_min <  0.001
         'p01',   'permutation p < 0.01',  @(a) a.p_min <  0.01
         'p05',   'permutation p < 0.05',  @(a) a.p_min <  0.05
         'logZ3', 'Rayleigh logZ > 3',     @(a) a.rayLogZ > 3
         'logZ2', 'Rayleigh logZ > 2',     @(a) a.rayLogZ > 2
         'logZ1', 'Rayleigh logZ > 1',     @(a) a.rayLogZ > 1 };
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
L = load(inFile,'A','xc','xs','SETS','COLS','prm');
A = L.A; xc = L.xc; xs = L.xs; COLS = L.COLS; prm = L.prm;

%% ---- select once per criterion, and pre-compute so the y-scales can be shared ----
SELS = cell(size(CRIT,1),1);  MM = cell(size(CRIT,1),3);  ZZ = cell(size(CRIT,1),3);
fprintf('%-24s %6s   %s\n','criterion','nCells','animals');
for ci = 1:size(CRIT,1)
    keep = arrayfun(CRIT{ci,3}, A);
    % A.tr is ITSELF a 1x3 cell, so {A.tr} is a cell of cells -- check inside it.
    hasTrace = arrayfun(@(a) any(isfinite(a.tr{2})), A);
    keep = keep(:).' & hasTrace(:).';
    SELS{ci} = find(keep);
    for q = 1:3
        Y = cell2mat(arrayfun(@(a) a.tr{q}, A(SELS{ci}), 'uni',0)');
        MM{ci,q} = Y;
        mu = mean(Y,2,'omitnan');  sd = std(Y,0,2,'omitnan');
        ZZ{ci,q} = (Y - mu) ./ max(sd, eps);
    end
    fprintf('%-24s %6d   %d\n', CRIT{ci,2}, numel(SELS{ci}), ...
            numel(unique({A(SELS{ci}).animal})));
end

% ONE y-range per column per row, across every criterion
YL = nan(2,3,2);
if shareY
    for q = 1:3
        for row = 1:2
            v = [];
            for ci = 1:size(CRIT,1)
                if row==1, v = [v; MM{ci,q}(:)]; else, v = [v; ZZ{ci,q}(:)]; end %#ok<AGROW>
            end
            yl = prctile(v(isfinite(v)), prcYLim);  pad = 0.08*diff(yl);
            YL(row,q,:) = [yl(1)-pad, yl(2)+pad];
        end
    end
end

%% ---- one figure per criterion ----
for ci = 1:size(CRIT,1)
    sel = SELS{ci};
    if isempty(sel), fprintf('  %s: no cells, skipped\n', CRIT{ci,2}); continue; end
    fh = figure('Color','w','Position',[30 30 1500 830]);
    tl = tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');
    for row = 1:2
        for q = 1:3
            ax = nexttile(tl,(row-1)*3+q); hold(ax,'on');
            if q < 3, x = xc; else, x = xs; end
            if row == 1, Y = MM{ci,q}; ylab = 'dF/F'; else, Y = ZZ{ci,q}; ylab = 'z (per trace)'; end
            for k = 1:size(Y,1)
                plot(ax, x, Y(k,:), '-', 'Color',[lineCol lineAlph], 'LineWidth',lineW);
            end
            xline(ax, 0, 'r-', 'LineWidth',1);
            yline(ax, 0, ':', 'Color',[0.55 0.55 0.55]);
            xlim(ax, x([1 end]));
            if shareY, ylim(ax, squeeze(YL(row,q,:)).');
            else
                yl = prctile(Y(:), prcYLim);  pad = 0.08*diff(yl);
                ylim(ax, [yl(1)-pad, yl(2)+pad]);
            end
            set(ax,'TickDir','out'); box(ax,'off');
            if q < 3, xlabel(ax, sprintf('time from %s (cycles)', COLS{q}));
            else,     xlabel(ax, 'time from own event (s)');
            end
            if q == 1, ylabel(ax, ylab); end
            title(ax, sprintf('%s  |  %s   (n=%d cells)', ylab, COLS{q}, size(Y,1)), ...
                  'FontSize',9,'Interpreter','none');
        end
    end
    nA = numel(unique({A(sel).animal}));
    title(tl, {sprintf('%s   |   every cell overlaid, NO average   |   n = %d cells, N = %d animals', ...
                       CRIT{ci,2}, numel(sel), nA), ...
               ['top = dF/F (a few loud cells dominate the ink);  bottom = each trace z-scored to its own window ' ...
                '(amplitude removed, shape and timing only).   breath cols in CYCLES, self col in SECONDS.'], ...
               'y-scales are SHARED across all six criterion figures, so they can be compared directly.'}, ...
          'FontWeight','bold','FontSize',10,'Interpreter','none');

    if doSave
        if ~isfolder(outDir), mkdir(outDir); end
        base = fullfile(outDir, sprintf('trig_avg_dff_overlay_%s', CRIT{ci,1}));
        exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
        exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
        fprintf('  saved %s.png/.pdf\n', CRIT{ci,1});
    end
    close(fh);
end
