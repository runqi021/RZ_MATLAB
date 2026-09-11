% cycle_dff_heatmap_groups_260812.m
% -----------------------------------------------------------------------
%  SINGLE-CYCLE dF/F HEATMAPS for the six selection criteria.
%
%      one row  = ONE BREATH CYCLE OF ONE CELL   (not a cell, not an average)
%      left     dF/F around that cycle's INSPIRATION ONSET
%      right    dF/F around that same cycle's INSPIRATORY PEAK
%      colour   raw dF/F, window mean removed per row
%
%  SORTED BY TIME FROM THE ONSET TO THAT CELL'S NEAREST EVENT, and BOTH panels
%  carry that ONE order.
%
%  THE ROWS ARE THE SAME CYCLES IN BOTH PANELS. That is the point, and it is why
%  the cycles are built from the ONSET train and the peak is then looked up
%  INSIDE each cycle rather than the two panels being built from their own
%  trigger lists. Triggering each panel separately gives two different row sets
%  (8840 vs 8849 cycles last time), which cannot share a sort order -- row 500
%  would be a different breath on the left than on the right, and comparing them
%  would be meaningless.
%
%  Cycles where the cell had NO event are kept, sorted last, below a cyan line.
%  They are the majority (a cell at 0.3 events/cycle is silent through 70% of
%  breaths) and dropping them would turn a sparse code into a dense one.
%
%  The orange line is the sort key itself. If the detected events are real
%  calcium transients the bright band must FOLLOW it. A vertical band instead --
%  fluorescence rising at a fixed point in the breath no matter when the event
%  was called -- is what a breathing-locked artifact looks like.
%
%  SIX CRITERIA, from two independent tests: the time-domain permutation p
%  (p < 0.001 / 0.01 / 0.05) and the occupancy-weighted Rayleigh logZ (> 3 / 2 / 1).
%  They are not nested: logZ > 1 selects cells the permutation never tested.
%
%  Input : event_latency_260811\event_latency_data.mat  (+ each recording's *_dFF.mat,
%          and breath_time_summary_260808 for the Rayleigh logZ)
%  Output: cycle_dff_heatmap_<crit>.png / .pdf   x6
%
%  Runqi Zhang / 2026-08-12

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath,'event_latency_260811','event_latency_data.mat');
rayCsv   = fullfile(rootPath,'breath_time_summary_260808','breath_time_summary_cells.csv');
outDir   = fullfile(rootPath,'event_latency_260811');
nCol     = 61;             % samples across the +/- half-cycle window (odd)
cycWin   = 0.5;            % window = +/- this many cycles
prcLim   = [0.5 99.5];     % dF/F colour clamp
maxRows  = 12000;          % even subsample WITHIN the sort order if exceeded
doSave   = true;

CRIT = { 'p001',  'permutation p < 0.001', @(pm,lz) pm <  0.001
         'p01',   'permutation p < 0.01',  @(pm,lz) pm <  0.01
         'p05',   'permutation p < 0.05',  @(pm,lz) pm <  0.05
         'logZ3', 'Rayleigh logZ > 3',     @(pm,lz) lz >  3
         'logZ2', 'Rayleigh logZ > 2',     @(pm,lz) lz >  2
         'logZ1', 'Rayleigh logZ > 1',     @(pm,lz) lz >  1 };
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
D = load(inFile,'S','CELL','OBS','REC','prm');
S = D.S; CELL = D.CELL; OBS = D.OBS; REC = D.REC; prm = D.prm;
assert(isfield(REC,'folder'), 'rerun Ventral_surface_event_latency_260811.m first');

p_min = arrayfun(@(s) min(s.p,[],'omitnan'), S).';
rayLZ = nan(numel(CELL),1);
if isfile(rayCsv)
    Bt = readtable(rayCsv);
    [tf,loc] = ismember(string({CELL.key}), string(Bt.cell_key));
    rayLZ(tf) = Bt.ray_logZ(loc(tf));
end
xg = linspace(-cycWin, cycWin, nCol);

%% ---- extract every cycle of every active cell, ONCE ----
% Both alignments and the sort key come out of the same pass, keyed to the same
% cycle, so the row sets can never diverge.
dffCache = containers.Map('KeyType','char','ValueType','any');
Aon = [];  Apk = [];  key = [];  cidx = [];
for c = 1:numel(CELL)
    for j = CELL(c).obs
        r  = REC(OBS(j).rec);
        on = r.trig{1};  pk = r.trig{2};
        if numel(on) < 2 || isempty(pk), continue; end
        if ~isKey(dffCache, r.folder)
            dh = dir(fullfile(r.folder,'*_dFF.mat'));
            if isempty(dh), dffCache(r.folder) = [];
            else, L = load(fullfile(dh(1).folder,dh(1).name),'dFF'); dffCache(r.folder) = L.dFF;
            end
        end
        dFF = dffCache(r.folder);
        if isempty(dFF), continue; end
        pp  = regexp(OBS(j).label,'/','split');
        rid = str2double(pp{end});
        if ~isfinite(rid) || rid < 1 || rid > size(dFF,2), continue; end
        y  = double(dFF(:,rid));  T = numel(y);
        tv = (1:T)/r.fps;
        ev = OBS(j).ev(:)/r.fps;
        ibi = mean(diff(on));  half = ibi/2;

        for k = 1:numel(on)-1
            t0 = on(k);
            % the PEAK BELONGING TO THIS CYCLE, not the nearest peak overall
            q = pk(pk > on(k) & pk < on(k+1));
            if isempty(q), continue; end
            t1 = q(1);
            qa = t0 + xg*ibi;   qb = t1 + xg*ibi;
            if qa(1) < tv(1) || qa(end) > tv(end), continue; end
            if qb(1) < tv(1) || qb(end) > tv(end), continue; end
            sa = interp1(tv,y,qa,'linear');  sb = interp1(tv,y,qb,'linear');
            Aon(end+1,:) = sa - mean(sa,'omitnan'); %#ok<AGROW>
            Apk(end+1,:) = sb - mean(sb,'omitnan'); %#ok<AGROW>
            % SORT KEY: onset -> nearest event of this cell, in cycle units
            if isempty(ev)
                key(end+1,1) = NaN; %#ok<AGROW>
            else
                [d,im] = min(abs(ev - t0));
                if d <= half, key(end+1,1) = (ev(im)-t0)/ibi; %#ok<AGROW>
                else,         key(end+1,1) = NaN;             %#ok<AGROW>
                end
            end
            cidx(end+1,1) = c; %#ok<AGROW>
        end
    end
    if mod(c,50)==0, fprintf('  %d/%d cells scanned, %d cycles\n', c, numel(CELL), size(Aon,1)); end
end
fprintf('%d cycles from %d cells\n', size(Aon,1), numel(unique(cidx)));

cl = prctile([Aon(:); Apk(:)], prcLim);
fprintf('dF/F colour scale %.3f to %.3f (shared by all six figures)\n', cl(1), cl(2));

%% ---- one figure per criterion ----
for ci = 1:size(CRIT,1)
    keepCell = find(CRIT{ci,3}(p_min, rayLZ));
    m = ismember(cidx, keepCell);
    if ~any(m), fprintf('  %s: no cycles, skipped\n', CRIT{ci,1}); continue; end
    A1 = Aon(m,:);  A2 = Apk(m,:);  kk = key(m);
    [ks,ord] = sort(kk);                       % NaN (no event) sorts LAST
    A1 = A1(ord,:);  A2 = A2(ord,:);
    nEvRow = nnz(~isnan(ks));
    if maxRows > 0 && size(A1,1) > maxRows
        pick = round(linspace(1,size(A1,1),maxRows));
        A1 = A1(pick,:); A2 = A2(pick,:); ks = ks(pick); nEvRow = nnz(~isnan(ks));
    end

    fh = figure('Color','w','Position',[30 30 1400 900]);
    tl = tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
    for q = 1:2
        if q==1, A = A1; nm = 'ONSET'; else, A = A2; nm = 'PEAK'; end
        ax = nexttile(tl,q);
        imagesc(ax, xg, 1:size(A,1), A, 'AlphaData', ~isnan(A));
        set(ax,'CLim',cl,'YDir','reverse','TickDir','out','Layer','top','Color',[1 1 1]);
        colormap(ax, gray(256)); hold(ax,'on');
        xline(ax, 0, 'r-', 'LineWidth',1.2);
        if nEvRow > 0 && nEvRow < size(A,1)
            yline(ax, nEvRow+0.5, 'c-', 'LineWidth',1.2);
        end
        if q == 1     % the sort key belongs on the panel it was measured from
            plot(ax, ks(1:nEvRow), 1:nEvRow, '-', 'Color',[1 0.55 0 0.9], 'LineWidth',1.4);
        end
        xlim(ax,[xg(1) xg(end)]);
        xlabel(ax, sprintf('time from %s (cycles)', nm));
        if q == 1, ylabel(ax,'breath cycles, sorted by time from ONSET to nearest event'); end
        title(ax, sprintf('%s-aligned', nm), 'FontWeight','bold');
        cb = colorbar(ax); cb.Label.String = 'dF/F';
    end
    title(tl, {sprintf('%s   |   %d cells, %d cycles (%d with an event)   |   same rows in both panels', ...
                       CRIT{ci,2}, numel(keepCell), size(A1,1), nEvRow), ...
               ['orange = the sort key (onset -> nearest event).  red = alignment landmark.  ' ...
                'cyan = below it the cell had NO event in that cycle.  window mean removed per row; colour scale shared across all six.']}, ...
          'FontWeight','bold','FontSize',10,'Interpreter','none');

    if doSave
        if ~isfolder(outDir), mkdir(outDir); end
        base = fullfile(outDir, sprintf('cycle_dff_heatmap_%s', CRIT{ci,1}));
        exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
        exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
        fprintf('  %-22s %5d cells %6d cycles -> saved %s\n', CRIT{ci,2}, ...
                numel(keepCell), size(A1,1), CRIT{ci,1});
    end
    close(fh);
end
