% cycle_dff_heatmap_260811.m
% -----------------------------------------------------------------------
%  dF/F HEATMAP, ALL CYCLES OF ALL SIGNIFICANT CELLS STACKED.
%
%      one row   = ONE BREATH CYCLE OF ONE CELL   (not one cell, not an average)
%      x         = time around the trigger, +/- half a cycle, in units of the
%                  cycle so recordings breathing at 0.5 s and 2.9 s can share
%                  an axis
%      colour    = dF/F, raw. Nothing averaged, nothing normalised per row.
%
%      LEFT  panel  ONSET-triggered
%      RIGHT panel  PEAK-triggered
%
%  SORTED BY TIME TO THE NEAREST DETECTED EVENT of that cell in that cycle.
%  Cycles are sorted independently within each panel, because a cycle's nearest
%  event sits at a different latency depending on which landmark you measure
%  from -- sharing one order between the panels would smear one of them.
%
%  WHAT THIS FIGURE IS FOR. The latency result rests entirely on detected EVENT
%  TIMES. This shows the underlying dF/F those events were detected from, with
%  no detector in the loop except for the sort. If the events are real calcium
%  transients, the bright band must run diagonally, tracking the sort key --
%  each cycle's fluorescence rising exactly where its event was called. A
%  diagonal that follows the sort is the detector agreeing with the raw signal.
%  A vertical band instead of a diagonal would mean the dF/F rises at a fixed
%  point in the breath regardless of where the event was called, which is what a
%  breathing-locked artifact looks like.
%
%  CYCLES WITH NO EVENT ARE KEPT, sorted last, and separated by a line. They are
%  the majority (a cell firing 0.3 events/cycle is silent in 70% of breaths) and
%  dropping them would turn a sparse code into an apparently dense one.
%
%  Input : event_latency_260811\event_latency_data.mat  + each recording's
%          *_ch1_dFF.mat (the only thing not already in that file)
%  Output: cycle_dff_heatmap_sig.png / .pdf
%
%  Runqi Zhang / 2026-08-11

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath,'event_latency_260811','event_latency_data.mat');
outDir   = fullfile(rootPath,'event_latency_260811');
nCol     = 61;             % columns across the +/- half-cycle window (odd -> a
                           % column centred exactly on the trigger)
prcLim   = [0.5 99.5];     % dF/F colour clamp, over the drawn rows only
maxRows  = 12000;          % if there are more cycles than this, subsample rows
                           % EVENLY WITHIN the sort order (never randomly, which
                           % would break the diagonal) -- 0 = keep everything
rngSeed  = 260811;
doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(rngSeed);
D = load(inFile,'S','CELL','OBS','REC','groups','prm');
S = D.S; CELL = D.CELL; OBS = D.OBS; REC = D.REC; prm = D.prm;
assert(isfield(REC,'folder'), ['this .mat predates the OBS/REC save -- rerun ' ...
       'Ventral_surface_event_latency_260811.m first']);

sig = find(any(vertcat(S.sig),2));
fprintf('%d significant cells (p<%.4g)\n', numel(sig), prm.alphaPerm);

TRIG = {'onset','peak'};
xg   = linspace(-0.5, 0.5, nCol);      % cycle units, 0 = trigger

%% ---- gather every cycle of every significant cell ----
dffCache = containers.Map('KeyType','char','ValueType','any');
P = struct('A',{[],[]},'key',{[],[]},'cell',{[],[]});
for q = 1:2
    A = [];  key = [];  cid = [];
    for ii = 1:numel(sig)
        c = sig(ii);
        for j = CELL(c).obs
            r  = REC(OBS(j).rec);
            ts = r.trig{q};
            if numel(ts) < 2, continue; end

            % dF/F for this recording, loaded once and cached
            if ~isKey(dffCache, r.folder)
                dh = dir(fullfile(r.folder,'*_dFF.mat'));
                if isempty(dh), dffCache(r.folder) = []; else
                    L = load(fullfile(dh(1).folder, dh(1).name),'dFF');
                    dffCache(r.folder) = L.dFF;
                end
            end
            dFF = dffCache(r.folder);
            if isempty(dFF), continue; end

            % ROI index is the trailing field of the label, and equals the
            % maskL label, which is the dFF column.
            pp  = regexp(OBS(j).label,'/','split');
            rid = str2double(pp{end});
            if ~isfinite(rid) || rid < 1 || rid > size(dFF,2), continue; end
            y   = double(dFF(:,rid));
            T   = numel(y);
            tv  = (1:T)/r.fps;                       % seconds
            ev  = OBS(j).ev(:)/r.fps;                % this cell's event times

            ibi  = mean(diff(ts));
            half = ibi/2;
            for k = 1:numel(ts)
                t0 = ts(k);
                tq = t0 + xg*ibi;                    % the +/- half cycle window
                if tq(1) < tv(1) || tq(end) > tv(end), continue; end
                A(end+1,:) = interp1(tv, y, tq, 'linear'); %#ok<AGROW>
                % sort key: signed latency to this cell's NEAREST event, in
                % cycle units. NaN when the cell was silent through this breath.
                if isempty(ev)
                    key(end+1,1) = NaN; %#ok<AGROW>
                else
                    [d, im] = min(abs(ev - t0));
                    if d <= half, key(end+1,1) = (ev(im)-t0)/ibi; %#ok<AGROW>
                    else,         key(end+1,1) = NaN;             %#ok<AGROW>
                    end
                end
                cid(end+1,1) = c; %#ok<AGROW>
            end
        end
    end
    P(q).A = A;  P(q).key = key;  P(q).cell = cid;
    fprintf('%-5s: %d cycles from %d cells | %d with an event (%.0f%%)\n', ...
        TRIG{q}, size(A,1), numel(unique(cid)), nnz(~isnan(key)), 100*mean(~isnan(key)));
end

%% ---- one colour scale across both panels ----
allv = [P(1).A(:); P(2).A(:)];  allv = allv(isfinite(allv));
cl = prctile(allv, prcLim);
fprintf('dF/F colour scale %.3f to %.3f\n', cl(1), cl(2));

%% ---- figure ----
fh = figure('Color','w','Position',[30 30 1420 900]);
tl = tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
for q = 1:2
    A = P(q).A;  key = P(q).key;
    [ks, ord] = sort(key);                 % NaN (no event) sorts LAST
    A = A(ord,:);
    nEvRow = nnz(~isnan(ks));

    if maxRows > 0 && size(A,1) > maxRows
        % EVEN subsample within the sort order -- never random, which would
        % scramble the diagonal this figure exists to show.
        pick = round(linspace(1, size(A,1), maxRows));
        A = A(pick,:);  ks = ks(pick);  nEvRow = nnz(~isnan(ks));
    end

    ax = nexttile(tl,q);
    imagesc(ax, xg, 1:size(A,1), A, 'AlphaData', ~isnan(A));
    set(ax,'CLim',cl,'YDir','reverse','TickDir','out','Layer','top','Color',[1 1 1]);
    colormap(ax, gray(256)); hold(ax,'on');
    xline(ax, 0, 'r-', 'LineWidth',1.2);                       % the trigger
    if nEvRow > 0 && nEvRow < size(A,1)
        yline(ax, nEvRow+0.5, 'c-', 'LineWidth',1.2);          % event / no-event
    end
    % the sort key itself, drawn over the map: the diagonal the dF/F should follow
    plot(ax, ks(1:nEvRow), 1:nEvRow, '-', 'Color',[1 0.55 0 0.9], 'LineWidth',1.4);
    xlim(ax,[xg(1) xg(end)]);
    xlabel(ax, sprintf('time from %s  (cycles)', upper(TRIG{q})));
    if q == 1, ylabel(ax,'breath cycles of all significant cells, sorted by time to nearest event'); end
    title(ax, sprintf('%s-triggered  |  %d cycles, %d with an event', ...
          upper(TRIG{q}), size(A,1), nEvRow), 'FontWeight','bold');
    cb = colorbar(ax); cb.Label.String = 'dF/F';
end

title(tl, {sprintf(['dF/F, every breath cycle of every significant cell, stacked   |   %d cells, p < %.4g   |   ' ...
                    'orange = the sort key (time to nearest detected event)'], numel(sig), prm.alphaPerm), ...
           ['red = trigger, cyan = boundary below which the cell had NO event in that cycle.   ' ...
            'Nothing averaged or normalised.  If the events are real transients the bright band follows the orange line.']}, ...
      'FontWeight','bold','FontSize',10,'Interpreter','none');

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir,'cycle_dff_heatmap_sig');
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    fprintf('saved %s.{png,pdf}\n', base);
end
