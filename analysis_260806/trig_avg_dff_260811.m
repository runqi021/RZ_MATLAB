% trig_avg_dff_260811.m
% -----------------------------------------------------------------------
%  TRIGGERED-AVERAGE dF/F, three triggers side by side, for the SIGNIFICANT
%  cells and for the tested-but-not-significant cells.
%
%      col 1  ONSET-triggered      breath, inspiration onset
%      col 2  PEAK-triggered       breath, inspiratory peak
%      col 3  SELF-triggered       this cell's OWN detected events
%
%  The third column is the control the other two are read against. It is the
%  same cell, the same trace, the same averaging -- triggered on its own events
%  instead of on the breath. It shows what ONE real calcium transient looks like
%  in this data: the amplitude, the rise, the decay. A breath-triggered average
%  can then be read as a fraction of it. If a cell's peak-triggered average is a
%  tenth of its self-triggered average, the cell fires on ~a tenth of breaths,
%  or fires at scattered latencies -- either way the breath explains little of
%  its activity, however significant the locking test was.
%
%  ---------------------------- ERROR BARS -------------------------------
%  SEM IS ACROSS ANIMALS, NOT ACROSS CELLS. Cells from one mouse share a
%  breathing pattern, a depth, a laser power and a day; they are not independent
%  samples of "a genotype". Averaging within an animal first and taking the SEM
%  over animal means answers the question actually being asked, and it is
%  brutally honest about n -- Sert contributes ONE animal, so it gets no error
%  bar at all. A cell-wise SEM would have looked ~5x tighter and meant nothing.
%  n is printed as BOTH counts on every panel: (N animals, n cells).
%
%  animal = one session date (Genotype/MMDD), the archive's own convention.
%
%  ---------------------------- BASELINES --------------------------------
%  Deliberately different, because the two questions are different:
%    breath-triggered  the WINDOW MEAN is subtracted. A cyclic average has no
%                      meaningful DC -- the trace wraps onto itself -- so what
%                      is plotted is modulation about the cycle mean.
%    self-triggered    a PRE-EVENT baseline is subtracted, so the height is the
%                      transient's real amplitude in dF/F.
%  Consequence: heights compare WITHIN a column freely, and BETWEEN columns only
%  as "modulation depth vs transient amplitude". Stated on the figure.
%
%  x-axis: breath columns in CYCLES (recordings breathe 0.5-2.9 s, so seconds
%  would not superimpose); the self column in SECONDS (a calcium transient has a
%  fixed duration in time, not in breaths).
%
%  Input : event_latency_260811\event_latency_data.mat + each recording's *_dFF.mat
%  Output: trig_avg_dff.png / .pdf / .csv  (+ per-cell traces in a .mat)
%
%  Runqi Zhang / 2026-08-11

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath,'event_latency_260811','event_latency_data.mat');
outDir   = fullfile(rootPath,'event_latency_260811');
nCycPts  = 61;                 % samples across the breath window (odd -> a
                               % sample exactly on the trigger)
cycWin   = 0.5;                % breath window = +/- this many cycles
selfWin  = [-1.0 2.0];         % self-triggered window, SECONDS
nSelfPts = 94;
selfBase = [-1.0 -0.5];        % pre-event baseline for the self column, s
showCells = true;              % draw the per-cell traces behind the mean
doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
D = load(inFile,'S','CELL','OBS','REC','groups','prm');
S = D.S; CELL = D.CELL; OBS = D.OBS; REC = D.REC; prm = D.prm;
assert(isfield(REC,'folder'), 'rerun Ventral_surface_event_latency_260811.m first');

tested = vertcat(S.tested);  sigM = vertcat(S.sig);
anySig = any(sigM,2);  anyTst = any(tested,2);
% Traces are computed for EVERY ACTIVE CELL, not just the tested ones, so that
% downstream figures can select on ANY criterion. The Rayleigh in particular
% reaches cells the permutation never tested: logZ > 1 selects 108 cells of
% which only 78 cleared the 20-event gate. Recomputing per criterion would mean
% re-reading every dF/F file each time.
SETS = { 'SIGNIFICANT',     find(anyTst &  anySig)
         'TESTED, NOT SIG', find(anyTst & ~anySig)
         'UNTESTED',        find(~anyTst) };
COLS = {'ONSET','PEAK','SELF'};
fprintf('%d significant | %d tested-not-sig | %d untested (all get traces)\n', ...
        nnz(anyTst & anySig), nnz(anyTst & ~anySig), nnz(~anyTst));

% Rayleigh logZ, joined by cell_key from the phase summary. Verified 276/276
% match -- both are built from the same cell_link identity, so a mismatch would
% mean one of the two summaries had drifted.
rayCsv = fullfile(rootPath,'breath_time_summary_260808','breath_time_summary_cells.csv');
rayLogZ = nan(numel(CELL),1);
if isfile(rayCsv)
    Bt = readtable(rayCsv);
    [tf,loc] = ismember(string({CELL.key}), string(Bt.cell_key));
    rayLogZ(tf) = Bt.ray_logZ(loc(tf));
    fprintf('Rayleigh logZ joined for %d of %d cells\n', nnz(tf), numel(CELL));
else
    warning('Rayleigh CSV not found -- logZ criteria will be unavailable');
end

xc = linspace(-cycWin, cycWin, nCycPts);          % cycle units
xs = linspace(selfWin(1), selfWin(2), nSelfPts);  % seconds

%% ---- per-cell triggered averages ----
dffCache = containers.Map('KeyType','char','ValueType','any');
A = struct('set',{},'cell',{},'animal',{},'group',{},'tr',{},'nTrig',{}, ...
           'p_onset',{},'p_peak',{},'p_min',{},'rayLogZ',{});
for si = 1:size(SETS,1)
    idx = SETS{si,2};
    for ii = 1:numel(idx)
        c = idx(ii);
        acc = {zeros(0,nCycPts), zeros(0,nCycPts), zeros(0,nSelfPts)};
        for j = CELL(c).obs
            r = REC(OBS(j).rec);
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

            for q = 1:3
                if q < 3
                    ts  = r.trig{q};
                    if numel(ts) < 2, continue; end
                    ibi = mean(diff(ts));
                    for k = 1:numel(ts)
                        tq = ts(k) + xc*ibi;
                        if tq(1) < tv(1) || tq(end) > tv(end), continue; end
                        seg = interp1(tv, y, tq, 'linear');
                        acc{q}(end+1,:) = seg - mean(seg,'omitnan');   % window mean
                    end
                else
                    for k = 1:numel(ev)
                        tq = ev(k) + xs;
                        if tq(1) < tv(1) || tq(end) > tv(end), continue; end
                        seg = interp1(tv, y, tq, 'linear');
                        bm  = mean(seg(xs >= selfBase(1) & xs <= selfBase(2)),'omitnan');
                        acc{q}(end+1,:) = seg - bm;                    % pre-event baseline
                    end
                end
            end
        end
        pp = regexp(CELL(c).label,'/','split');
        A(end+1) = struct('set',SETS{si,1},'cell',c, ...
            'animal',sprintf('%s/%s',pp{1},pp{2}), 'group',pp{1}, ...
            'tr',{cellfun(@(x) mean(x,1,'omitnan'), acc, 'uni',0)}, ...
            'nTrig',cellfun(@(x) size(x,1), acc), ...
            'p_onset',S(c).p(1), 'p_peak',S(c).p(2), ...
            'p_min',min(S(c).p,[],'omitnan'), 'rayLogZ',rayLogZ(c)); %#ok<AGROW>
    end
    fprintf('  %-16s %d cells done\n', SETS{si,1}, numel(idx));
end

%% ---- population: average WITHIN animal first, then SEM ACROSS animals ----
% The population figure keeps its original two sets; UNTESTED cells have traces
% saved but are not drawn there (a mean over 6-event cells is noise).
POP = struct('set',{},'q',{},'m',{},'sem',{},'nA',{},'nC',{},'animals',{});
for si = 1:2
    sel = strcmp({A.set}, SETS{si,1});
    for q = 1:3
        M = cell2mat(arrayfun(@(a) a.tr{q}, A(sel), 'uni',0)');
        an = {A(sel).animal};
        ok = ~all(isnan(M),2);
        M = M(ok,:);  an = an(ok);
        [ua,~,ia] = unique(an);
        AM = zeros(numel(ua), size(M,2));
        for k = 1:numel(ua), AM(k,:) = mean(M(ia==k,:),1,'omitnan'); end
        POP(end+1) = struct('set',SETS{si,1},'q',q, ...
            'm',mean(AM,1,'omitnan'), ...
            'sem',std(AM,0,1,'omitnan')/sqrt(max(size(AM,1),1)), ...
            'nA',size(AM,1),'nC',size(M,1),'animals',{ua}); %#ok<AGROW>
    end
end

fprintf('\n%-16s %-6s %8s %8s %12s\n','set','trig','animals','cells','peak dF/F');
for i = 1:numel(POP)
    fprintf('%-16s %-6s %8d %8d %12.4f\n', POP(i).set, COLS{POP(i).q}, ...
        POP(i).nA, POP(i).nC, max(POP(i).m));
end

%% ---- figure ----
fh = figure('Color','w','Position',[30 30 1500 850]);
tl = tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');
CL = [0.15 0.45 0.85; 0.55 0.55 0.55];       % sig blue, non-sig grey
for si = 1:2
    sel = find(strcmp({A.set}, SETS{si,1}));
    for q = 1:3
        ax = nexttile(tl, (si-1)*3 + q); hold(ax,'on');
        if q < 3, x = xc; else, x = xs; end
        if showCells
            for k = sel
                plot(ax, x, A(k).tr{q}, '-', 'Color',[0.6 0.6 0.6 0.18], 'LineWidth',0.4);
            end
        end
        P = POP(strcmp({POP.set},SETS{si,1}) & [POP.q]==q);
        if P.nA > 1
            fill(ax, [x fliplr(x)], [P.m+P.sem fliplr(P.m-P.sem)], CL(si,:), ...
                 'FaceAlpha',0.30, 'EdgeColor','none');
        end
        plot(ax, x, P.m, '-', 'Color',CL(si,:), 'LineWidth',2);
        xline(ax, 0, 'r-', 'LineWidth',1);
        yline(ax, 0, ':', 'Color',[0.5 0.5 0.5]);
        xlim(ax, x([1 end])); set(ax,'TickDir','out'); box(ax,'off');
        if q < 3, xlabel(ax, sprintf('time from %s (cycles)', COLS{q}));
        else,     xlabel(ax, 'time from own event (s)');
        end
        if q == 1, ylabel(ax, sprintf('%s\ndF/F (window mean removed)', SETS{si,1})); end
        ttl = sprintf('%s  |  %s   (N=%d animals, n=%d cells)', SETS{si,1}, COLS{q}, P.nA, P.nC);
        if P.nA == 1, ttl = [ttl '  -- ONE ANIMAL, no SEM']; end %#ok<AGROW>
        title(ax, ttl, 'FontSize',9, 'Interpreter','none');
    end
end
title(tl, {sprintf(['Triggered-average dF/F   |   thin grey = one cell, thick = mean over ANIMALS, ' ...
                    'band = SEM ACROSS ANIMALS (not cells)   |   sig: p < %.4g'], prm.alphaPerm), ...
           ['breath columns: window mean removed, x in CYCLES.   self column: pre-event baseline removed, x in SECONDS.  ' ...
            'So heights compare within a column, and between columns only as modulation depth vs transient amplitude.']}, ...
      'FontWeight','bold','FontSize',10,'Interpreter','none');

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir,'trig_avg_dff');
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    save([base '.mat'],'A','POP','xc','xs','SETS','COLS','prm','-v7.3');
    T = table({A.set}', {CELL([A.cell]).key}', {A.animal}', {A.group}', ...
        [A.p_onset]', [A.p_peak]', [A.p_min]', [A.rayLogZ]', ...
        arrayfun(@(a) max(a.tr{1}), A)', arrayfun(@(a) max(a.tr{2}), A)', ...
        arrayfun(@(a) max(a.tr{3}), A)', vertcat(A.nTrig), ...
        'VariableNames',{'set','cell_key','animal','group', ...
        'p_onset','p_peak','p_min','ray_logZ', ...
        'peak_onset_dff','peak_peak_dff','peak_self_dff','n_trig_onset_peak_self'});
    writetable(T,[base '.csv']);
    fprintf('\nsaved %s.{png,pdf,mat,csv}\n', base);
end
