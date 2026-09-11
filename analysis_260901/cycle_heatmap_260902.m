% cycle_heatmap_260902.m
% -----------------------------------------------------------------------
%  Per-cycle dF/F heatmap, one row per breath cycle, sorted by the latency
%  between the landmark and the nearest event ON ONE SIDE of it.
%
%      rows      inspiration ONSET triggered  |  inspiratory PEAK triggered
%      columns   sorted by time AFTER  |  sorted by time BEFORE
%      4 panels per figure
%
%      one figure per (set x class): significant and non-significant cells are
%      separate figures, so nothing is pooled across the significance split.
%
%  THE TWO COLUMNS ARE MIRROR RULES, both searched within half the pooled median
%  cycle and both sorted DESCENDING:
%      after  : dt = min(ev - trig) over ev >= trig   (first event after)
%      before : dt = max(ev - trig) over ev <= trig   (last event before)
%  A cycle appears in the after column only if it has an event after the
%  landmark, and in the before column only if it has one before; most cycles
%  appear in neither, some in both. The counts differ per column and are printed
%  in each panel title.
%
%  WHY NOT NEAREST-EVENT. A nearest search lets an event that PRECEDED the
%  landmark claim the row, so a panel meant to read as a response latency ends up
%  containing negative dt. Splitting the two directions into their own panels
%  keeps each one a latency in a defined direction.
%
%  BREATH STRIP above every panel, as in the per-cell summary figure. Each
%  recording's breath is SIGN-CORRECTED first -- the trace is SVD PC1, whose sign
%  is arbitrary per video (project_breath_sign_is_analysis_immune), so averaging
%  raw traces across recordings partly cancels the waveform. Each recording is
%  flipped, if needed, so its own peak-triggered average is positive at t = 0,
%  and only then averaged, per RECORDING rather than per cell or per cycle. The
%  strip is min-max scaled and carries no y axis: only its timing means anything.
%
%  RAW dF/F IN GRAYSCALE with clim = prctile([0.5 99.5]) and NOTHING normalised,
%  the archive convention (project_breath_trig_heatmap_260806); parula is
%  reserved for z-scored panels. useZ = true switches both together. Colour
%  limits are computed PER FIGURE and printed on it.
%
%  TIME, NOT PHASE, on the x axis: a fixed WIN_S-wide window, so a cycle that ran
%  long extends past the edge rather than being stretched to fit. The window and
%  the dt search are separate -- the window shows the neighbouring cycles, the
%  search stays inside this one.
%
%  clear and clc but NOT close all: figures accumulate so sets can be compared
%  side by side (project_breath_trig_heatmap_260806).
%
%  Runqi Zhang / 2026-09-02
% -----------------------------------------------------------------------
clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot,'analysis_260806'));

%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
bundleDir = fullfile(sumRoot,'per-cell-summary_active_260812','spike_recheck_260901');
outDir    = fullfile(sumRoot,'polar_recheck_260901');
sigCsv    = fullfile(outDir,'polar_recheck_perm_percell.csv');

% {display name, file tag, genotype prefix ('' = all), dates, mode}
% Vglut2 intact is NOT one breathing state. Per-session median cycle, measured
% 2026-09-02: 0224 = 0.467 s, 0728 = 0.650 s, 0810 = 1.083 s, 1124 = 2.717 s,
% against 0824 (vagotomised) = 1.267 s. So the fast session 0224 is split out,
% and note that 1124 is SLOWER than the vagotomised animals -- pooling all four
% intact sessions averages a 2.1 Hz animal with a 0.37 Hz one.
SETS = { 'Sert',                     'Sert',              'Sert',   {},               'all'
         'Vglut2 intact FAST (0224)','Vglut2_intact_fast','Vglut2', {'0224'},         'include'
         'Vglut2 intact rest',       'Vglut2_intact_rest','Vglut2', {'0224','0824'},  'exclude'
         'Vglut2 vagal intact',      'Vglut2_intact',     'Vglut2', {'0824'},         'exclude'
         'Vglut2 vagotomised',       'Vglut2_vagotom',    'Vglut2', {'0824'},         'include'
         'Vglut2',                   'Vglut2',            'Vglut2', {},               'all'
         'ALL GENOTYPES',            'ALL',               '',       {},               'all' };
useZ      = false;            % false = raw dF/F in gray (archive convention)
WIN_S     = 2;                % TOTAL width of the plotted window, seconds
nTau      = 201;
climPct   = [0.5 99.5];
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

CU = ca_recheck_load_curation_260901(bundleDir);
Ts = readtable(sigCsv,'TextType','string');
sigMap = containers.Map(cellstr(Ts.stem), num2cell(logical(Ts.sig_perm)));

TRIG = {'onset','peak'};
CLS  = {'significant','non-significant'};
CTAG = {'sig','nonsig'};
DIRN = {'after','before'};
halfWin = WIN_S/2;
tau = linspace(-halfWin, halfWin, nTau);

for gsel = 1:size(SETS,1)
    SETNAME  = SETS{gsel,1};   tag   = SETS{gsel,2};
    GENOTYPE = SETS{gsel,3};   DATES = SETS{gsel,4};   DMODE = SETS{gsel,5};

    d = dir(fullfile(bundleDir,'cells','*.mat'));
    if ~isempty(GENOTYPE)
        d = d(startsWith({d.name}, [GENOTYPE '_']));
    end
    if ~strcmpi(DMODE,'all') && ~isempty(DATES)
        hit = false(1,numel(d));
        for v = 1:numel(DATES), hit = hit | contains({d.name}, ['_' DATES{v} '_']); end
        if strcmpi(DMODE,'include'), d = d(hit); else, d = d(~hit); end
    end
    assert(~isempty(d), 'no cells matched %s', SETNAME);

    % dt search half-width = half the pooled median cycle of THIS set
    allDt = [];
    for k = 1:numel(d)
        L = load(fullfile(d(k).folder,d(k).name),'RECc');
        for i = 1:numel(L.RECc)
            r = L.RECc(i);  g = CU.get(r.folder, r.roi);
            if g.found && (g.toss || g.cell_toss), continue; end
            if numel(r.foot) > 1, allDt = [allDt; diff(r.foot)/r.fps]; end %#ok<AGROW>
        end
    end
    halfCyc = median(allDt)/2;
    fprintf('\n=== %s ===\npooled median cycle %.3f s | display +/- %.2f s | dt search %.3f s each side\n', ...
            SETNAME, 2*halfCyc, halfWin, halfCyc);

    % E{trigger, class, direction}
    E = cell(2,2,2);  DT = cell(2,2,2);  B = cell(2,2);
    for q=1:2, for c=1:2, B{q,c}=[]; for v=1:2, E{q,c,v}=[]; DT{q,c,v}=[]; end, end, end
    nCellCls = zeros(1,2);
    seenRec  = containers.Map('KeyType','char','ValueType','logical');

    for k = 1:numel(d)
        stem  = erase(d(k).name,'.mat');
        isSig = isKey(sigMap,stem) && sigMap(stem);
        cls   = 1 + double(~isSig);
        L = load(fullfile(d(k).folder,d(k).name),'RECc');
        used = false;
        for i = 1:numel(L.RECc)
            r = L.RECc(i);  g = CU.get(r.folder, r.roi);
            if g.found && (g.toss || g.cell_toss), continue; end
            ev = r.spike_idx(:);
            if g.found, ev = g.spike_idx(:); end
            x = r.dff;
            if useZ, x = (x - mean(x,'omitnan')) / max(std(x,'omitnan'), eps); end
            w  = max(2, round(halfWin*r.fps));
            tk = (-w:w)/r.fps;

            % breath, sign fixed on this recording's OWN peak-triggered average
            bz  = (r.breath - mean(r.breath)) / max(std(r.breath), eps);
            pkT = r.peak(r.peak-w >= 1 & r.peak+w <= r.T);
            if ~isempty(pkT)
                Ep = zeros(numel(pkT), 2*w+1);
                for c = 1:numel(pkT), Ep(c,:) = bz(pkT(c)-w : pkT(c)+w); end
                v0 = mean(Ep(:, w+1));
                if isfinite(v0) && v0 < 0, bz = -bz; end
            end

            for q = 1:2
                if q == 1, trg = r.foot(:); else, trg = r.peak(:); end
                trg = trg(trg-w >= 1 & trg+w <= r.T);
                if isempty(trg), continue; end
                used = true;
                M   = zeros(numel(trg), 2*w+1);
                dtA = nan(numel(trg),1);  dtB = nan(numel(trg),1);
                for c = 1:numel(trg)
                    M(c,:) = x(trg(c)-w : trg(c)+w);
                    if isempty(ev), continue; end
                    dd  = (ev - trg(c))/r.fps;
                    aft = dd(dd >= 0 & dd <=  halfCyc);   % first event AFTER
                    bef = dd(dd <= 0 & dd >= -halfCyc);   % last  event BEFORE
                    if ~isempty(aft), dtA(c) = min(aft); end
                    if ~isempty(bef), dtB(c) = max(bef); end
                end
                Mi = single(interp1(tk, M.', tau, 'linear', NaN).');
                E{q,cls,1}  = [E{q,cls,1};  Mi];   DT{q,cls,1} = [DT{q,cls,1}; dtA];
                E{q,cls,2}  = [E{q,cls,2};  Mi];   DT{q,cls,2} = [DT{q,cls,2}; dtB];

                rkey = sprintf('%s|%d|%d', r.folder, q, cls);
                if ~isKey(seenRec, rkey)
                    seenRec(rkey) = true;
                    Eb = zeros(numel(trg), 2*w+1);
                    for c = 1:numel(trg), Eb(c,:) = bz(trg(c)-w : trg(c)+w); end
                    B{q,cls} = [B{q,cls}; interp1(tk, mean(Eb,1,'omitnan'), tau, 'linear', NaN)];
                end
            end
        end
        if used, nCellCls(cls) = nCellCls(cls) + 1; end
    end

    % keep only cycles that HAVE an event on that side, then sort descending
    NC = zeros(2,2,2);
    for q = 1:2
        for c = 1:2
            for v = 1:2
                dt = DT{q,c,v};  keep = isfinite(dt);
                NC(q,c,v) = nnz(~keep);
                E{q,c,v}  = E{q,c,v}(keep,:);  dt = dt(keep);
                [~, ord] = sort(dt, 'descend');
                E{q,c,v} = E{q,c,v}(ord,:);  DT{q,c,v} = dt(ord);
            end
        end
    end

    %% ---------------- one figure per class ----------------
    for c = 1:2
        if nCellCls(c) == 0, continue; end
        sub  = E(:,c,:);
        allv = cell2mat(cellfun(@(m) double(m(:)), sub(:), 'uni',0));
        cl   = prctile(allv(isfinite(allv)), climPct);

        fig = figure('Color','w','Units','centimeters','Position',[1 1 22 20], ...
                     'Name',sprintf('%s / %s per-cycle dF/F', SETNAME, CLS{c}));
        set(fig,'DefaultAxesFontSize',9);
        x0 = [0.095 0.545];  y0 = [0.545 0.085];
        hHm = 0.295;  hStrip = 0.048;

        for q = 1:2
            for v = 1:2
                ax = axes(fig,'Position',[x0(v), y0(q), 0.375, hHm]); %#ok<LAXES>
                M = E{q,c,v};
                if isempty(M), axis(ax,'off'); continue; end
                imagesc(ax, tau, 1:size(M,1), M);
                set(ax,'YDir','normal','CLim',cl,'TickDir','out','Layer','top');
                colormap(ax, ternary(useZ, parula(256), gray(256)));
                hold(ax,'on');
                plot(ax, [0 0], [0.5 size(M,1)+0.5], '-', 'Color',[0.90 0.15 0.10], 'LineWidth',1.0);
                hold(ax,'off');
                xlim(ax,[tau(1) tau(end)]);  ylim(ax,[0.5 size(M,1)+0.5]);
                xlabel(ax, sprintf('time from insp %s (s)', TRIG{q}));
                ylabel(ax, sprintf('cycle (sorted by time %s)', DIRN{v}));

                axb = axes(fig,'Position',[x0(v), y0(q)+hHm+0.004, 0.375, hStrip]); %#ok<LAXES>
                hold(axb,'on');
                if ~isempty(B{q,c})
                    b = mean(B{q,c},1,'omitnan');
                    b = (b - min(b)) / max(max(b)-min(b), eps);
                    plot(axb, tau, b, '-', 'Color',[0.25 0.25 0.25], 'LineWidth',1.2);
                end
                plot(axb, [0 0], [-0.05 1.05], '-', 'Color',[0.90 0.15 0.10], 'LineWidth',1.0);
                hold(axb,'off');
                xlim(axb,[tau(1) tau(end)]);  ylim(axb,[-0.05 1.05]);
                set(axb,'XTick',[],'YTick',[],'Box','on');
                if v == 1, ylabel(axb,'breath','FontSize',7.5); end
                title(axb, {sprintf('%s  |  event %s the landmark', upper(TRIG{q}), upper(DIRN{v})), ...
                            sprintf('%d cycles (%d without), %d recordings', ...
                                    size(M,1), NC(q,c,v), size(B{q,c},1))}, ...
                      'FontSize',8.5,'FontWeight','normal','Interpreter','none');
            end
        end

        cb = colorbar(ax,'Position',[0.945 0.085 0.014 0.295]);
        cb.Label.String = ternary(useZ,'dF/F (z-scored)','dF/F');

        sgtitle({sprintf('%s  |  %s cells (n = %d)  |  per-cycle dF/F', ...
                         SETNAME, CLS{c}, nCellCls(c)), ...
                 sprintf(['left = sorted by time to the FIRST event AFTER, right = to the ' ...
                          'LAST event BEFORE; both descending, window \\pm %.1f s'], halfWin), ...
                 sprintf(['%s  |  clim %.3g to %.3g (%g-%g pct)  |  dt searched %.0f ms ' ...
                          'each side  |  strip = trigger-averaged breath, sign-corrected'], ...
                          ternary(useZ,'z-scored, parula','raw dF/F, grayscale'), ...
                          cl(1), cl(2), climPct(1), climPct(2), 1000*halfCyc)}, ...
                'Interpreter','tex','FontSize',9);

        stem = fullfile(outDir, sprintf('%s_%s_cycle_dff_heatmap', tag, CTAG{c}));
        exportgraphics(fig,[stem '.png'],'Resolution',250,'BackgroundColor','white');
        exportgraphics(fig,[stem '.pdf'],'ContentType','image','BackgroundColor','white');
        fprintf('  %-16s n=%3d cells', CLS{c}, nCellCls(c));
        for q = 1:2
            for v = 1:2
                fprintf('  | %s-%s %d (med %.0f ms)', TRIG{q}(1:2), DIRN{v}(1:3), ...
                        size(E{q,c,v},1), 1000*median(DT{q,c,v}));
            end
        end
        fprintf('\n    -> %s.png\n', stem);
    end
end

%% ===================== LOCAL =====================
function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end
