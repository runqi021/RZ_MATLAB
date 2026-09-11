% cycle_heatmap_phase_260902.m
% -----------------------------------------------------------------------
%  Per-cycle dF/F heatmap on a SINGLE-LANDMARK PHASE axis, three cycles wide.
%
%      rows      inspiration ONSET cycles  |  inspiratory PEAK cycles
%      columns   sorted by EARLIEST event  |  sorted by LATEST event
%      4 panels per figure, one figure per (set x class)
%
%  THE AXIS IS PHASE, 0 to nCycShow. Cycle k runs from its own landmark to the
%  next and is mapped to 0..1, the one after to 1..2, the one after that to 2..3.
%  EACH CYCLE IS NORMALISED BY ITS OWN DURATION, so phases 1.0 and 2.0 are always
%  landmarks no matter how long those particular breaths were -- which is the
%  whole point of plotting in phase rather than in seconds, where the later
%  landmarks smear out across the panel. Every extra cycle costs rows: the last
%  nCycShow cycles of every recording have no successor and drop out.
%
%  ONE LANDMARK, NOT TWO. Phase is 2*pi*(t - t_k)/(t_k+1 - t_k), linear in time
%  across the whole cycle. The other landmark is not pinned anywhere; its median
%  position is drawn as the sky dashed line, and it is a measured quantity.
%
%  SORTING NEEDS NO DIRECTION RULE HERE. On a cycle that starts at its own
%  landmark every event is after it by construction, so there is no
%  first-after / last-before split and no nearest-event ambiguity -- the key is
%  just an event phase in [0,1). The two columns use the EARLIEST and the LATEST
%  event of cycle 1, both descending. THE ROWS ARE THE SAME CYCLES in both
%  columns, only reordered: a cycle qualifies by having at least one event, and
%  for a cycle with exactly one event the two keys are identical. Where the two
%  columns differ, that cycle fired more than once.
%
%  ONLY CYCLES WITH AN EVENT IN CYCLE 1 are drawn; the count without is in each
%  title. A cycle also needs all nCycShow-1 SUCCESSORS to exist inside the
%  recording. Only cycle 1 is selected on and sorted by; cycles 2 and 3 are
%  whatever followed, unsorted and unconditioned, which is why they wash out.
%
%  ROWS DO NOT OVERLAP: strideCyc = nCycShow, so the recording is tiled into
%  CONSECUTIVE, NON-OVERLAPPING blocks of nCycShow cycles and every breath is
%  plotted exactly once, in exactly one slot.
%
%  This was not true before. With a stride of 1 each cycle appeared in up to
%  nCycShow rows -- once as its own anchor, and again as the successor of the
%  previous anchors -- verified numerically: row i's cycle-2 block equalled row
%  i+1's cycle-1 block to 1e-15. Phase 0..1 was still one-cycle-one-row, but the
%  later slots were the same breaths redrawn at a shift, so any structure there
%  was not independent evidence. Set strideCyc = 1 to get the sliding version
%  back, and read the later cycles accordingly if you do.
%
%  THE PRICE OF THE TILING is that a cycle only gets to be cycle 1 if it lands on
%  an anchor. A cycle with an event sitting at slot 2 or 3 of its block is drawn,
%  but not selected on and not sorted by -- so the row count falls by rather more
%  than nCycShow, and which cycles anchor depends on where each recording starts.
%
%  CYCLES MORE THAN maxCycFac x THE RECORDING MEDIAN ARE DROPPED, for cycle 1 and
%  its successor both: a missed landmark produces a "cycle" several times too
%  long, and normalising it to 0..1 would compress a genuine breath into a
%  fraction of the axis and pollute every row it touches.
%
%  BREATH STRIP above each panel, resampled the same way, averaged per RECORDING,
%  with each recording's PC1 sign fixed on its own peak-triggered average first
%  (the sign is arbitrary per video -- project_breath_sign_is_analysis_immune).
%
%  RAW dF/F IN GRAYSCALE, clim = prctile([0.5 99.5]), nothing normalised
%  (project_breath_trig_heatmap_260806).
%
%  clear and clc but NOT close all.
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

% Vglut2 intact is NOT one breathing state. Per-session median cycle, measured
% 2026-09-02: 0224 = 0.467 s, 0728 = 0.650 s, 0810 = 1.083 s, 1124 = 2.717 s,
% against 0824 (vagotomised) = 1.267 s. So the fast session 0224 is split out,
% and note that 1124 is SLOWER than the vagotomised animals -- pooling all four
% intact sessions averages a 2.1 Hz animal with a 0.37 Hz one.
% ChAT/0826 is the vagotomised ChAT session (source folder
% C:\Users\Admin\Desktop\260826_ChAT-soma-g8s_vagotomized); Sert 0828+0829 are
% the vagotomised Sert sessions; Vglut2/0824 the vagotomised Vglut2 one.
SETS = { 'Sert intact',              'Sert_intact',       'Sert',   {'0828','0829'},  'exclude'
         'Sert vagotomised',         'Sert_vagotom',      'Sert',   {'0828','0829'},  'include'
         'Sert',                     'Sert',              'Sert',   {},               'all'
         'ChAT intact',              'ChAT_intact',       'ChAT',   {'0826'},         'exclude'
         'ChAT vagotomised',         'ChAT_vagotom',      'ChAT',   {'0826'},         'include'
         'Vglut2 intact FAST (0224)','Vglut2_intact_fast','Vglut2', {'0224'},         'include'
         'Vglut2 intact rest',       'Vglut2_intact_rest','Vglut2', {'0224','0824'},  'exclude'
         'Vglut2 vagal intact',      'Vglut2_intact',     'Vglut2', {'0824'},         'exclude'
         'Vglut2 vagotomised',       'Vglut2_vagotom',    'Vglut2', {'0824'},         'include'
         'Vglut2',                   'Vglut2',            'Vglut2', {},               'all'
         'ALL GENOTYPES',            'ALL',               '',       {},               'all' };
useZ      = false;
nCycShow  = 3;           % how many consecutive cycles to display
strideCyc = 3;           % anchor step: nCycShow = non-overlapping, 1 = sliding
nPh       = 301;         % samples across the 0..nCycShow phase axis
maxCycFac = 3;           % drop cycles longer than this x the recording median
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
ph   = linspace(0, nCycShow, nPh);

for gsel = 1:size(SETS,1)
    SETNAME  = SETS{gsel,1};   tag   = SETS{gsel,2};
    GENOTYPE = SETS{gsel,3};   DATES = SETS{gsel,4};   DMODE = SETS{gsel,5};

    d = dir(fullfile(bundleDir,'cells','*.mat'));
    if ~isempty(GENOTYPE), d = d(startsWith({d.name}, [GENOTYPE '_'])); end
    if ~strcmpi(DMODE,'all') && ~isempty(DATES)
        hit = false(1,numel(d));
        for v = 1:numel(DATES), hit = hit | contains({d.name}, ['_' DATES{v} '_']); end
        if strcmpi(DMODE,'include'), d = d(hit); else, d = d(~hit); end
    end
    assert(~isempty(d), 'no cells matched %s', SETNAME);

    E = cell(2,2);  PH1 = cell(2,2);  PH2 = cell(2,2);
    B = cell(2,2);  OTH = cell(2,2);
    for q=1:2, for c=1:2
        E{q,c}=[]; PH1{q,c}=[]; PH2{q,c}=[]; B{q,c}=[]; OTH{q,c}=[];
    end, end
    NC = zeros(2,2);  nCellCls = zeros(1,2);  nDrop = 0;
    seenRec = containers.Map('KeyType','char','ValueType','logical');

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

            % breath, sign fixed on this recording's own peak-triggered average
            bz  = (r.breath - mean(r.breath)) / max(std(r.breath), eps);
            wS  = max(2, round(0.25*r.fps));
            pkT = r.peak(r.peak-wS >= 1 & r.peak+wS <= r.T);
            if ~isempty(pkT)
                v0 = mean(arrayfun(@(t) bz(t), pkT));
                if isfinite(v0) && v0 < 0, bz = -bz; end
            end

            for q = 1:2
                if q == 1, trg = r.foot(:); oth = r.peak(:);
                else,      trg = r.peak(:); oth = r.foot(:); end
                if numel(trg) < 3, continue; end
                used = true;
                dt  = diff(trg);
                med = median(dt);
                Mi = [];  p1 = [];  p2 = [];  Bi = [];  op = [];
                for cix = 1:strideCyc:numel(dt)-(nCycShow-1)
                    lens = dt(cix : cix+nCycShow-1);
                    if any(lens > maxCycFac*med)
                        nDrop = nDrop + 1;  continue
                    end
                    a = trg(cix);  b = trg(cix+1);
                    if a < 1 || trg(cix+nCycShow) > r.T, continue; end
                    % every cycle normalised by ITS OWN length, laid end to end
                    uu = [];  vv = [];  bb = [];
                    for cj = 0:nCycShow-1
                        s0 = trg(cix+cj);  s1 = trg(cix+cj+1);
                        idx = (s0:s1)';
                        qj  = cj + (idx - s0)/(s1 - s0);
                        if cj == 0
                            uu = qj;         vv = x(idx);         bb = bz(idx);
                        else
                            uu = [uu; qj(2:end)];   %#ok<AGROW>
                            vv = [vv; x(idx(2:end))];   %#ok<AGROW>
                            bb = [bb; bz(idx(2:end))];  %#ok<AGROW>
                        end
                    end
                    Mi(end+1,:) = interp1(uu, vv, ph, 'linear', NaN); %#ok<AGROW>
                    Bi(end+1,:) = interp1(uu, bb, ph, 'linear', NaN); %#ok<AGROW>

                    e = ev(ev >= a & ev < b);
                    if isempty(e)
                        p1(end+1,1) = NaN;  p2(end+1,1) = NaN; %#ok<AGROW>
                    else
                        p1(end+1,1) = (min(e)-a)/(b-a); %#ok<AGROW>
                        p2(end+1,1) = (max(e)-a)/(b-a); %#ok<AGROW>
                    end
                    o = oth(oth > a & oth < b);
                    if ~isempty(o), op(end+1,1) = (o(1)-a)/(b-a); end %#ok<AGROW>
                end
                if isempty(Mi), continue; end
                E{q,cls}   = [E{q,cls};   single(Mi)];
                PH1{q,cls} = [PH1{q,cls}; p1];
                PH2{q,cls} = [PH2{q,cls}; p2];
                OTH{q,cls} = [OTH{q,cls}; op];
                rkey = sprintf('%s|%d|%d', r.folder, q, cls);
                if ~isKey(seenRec, rkey)
                    seenRec(rkey) = true;
                    B{q,cls} = [B{q,cls}; mean(Bi,1,'omitnan')];
                end
            end
        end
        if used, nCellCls(cls) = nCellCls(cls) + 1; end
    end

    % keep cycles with an event in cycle 1; the SAME rows are then ordered two
    % ways, by the earliest and by the latest event phase, both descending
    ORD = cell(2,2,2);
    for q = 1:2
        for c = 1:2
            p1 = PH1{q,c};  keep = isfinite(p1);
            NC(q,c) = nnz(~keep);
            E{q,c}   = E{q,c}(keep,:);
            PH1{q,c} = p1(keep);  PH2{q,c} = PH2{q,c}(keep);
            [~, ORD{q,c,1}] = sort(PH1{q,c}, 'descend');
            [~, ORD{q,c,2}] = sort(PH2{q,c}, 'descend');
        end
    end

    allv = cell2mat(cellfun(@(m) double(m(:)), E(:), 'uni',0));
    cl   = prctile(allv(isfinite(allv)), climPct);
    fprintf('\n=== %s ===  (%d long cycles dropped)\n', SETNAME, nDrop);

    %% ---------------- one figure per class ----------------
    SKEY = {'EARLIEST','LATEST'};
    CTAG = {'sig','nonsig'};
    for c = 1:2
    if nCellCls(c) == 0, continue; end
    fig = figure('Color','w','Units','centimeters','Position',[1 1 22 20], ...
                 'Name',sprintf('%s / %s per-cycle dF/F, phase', SETNAME, CLS{c}));
    set(fig,'DefaultAxesFontSize',9);
    x0 = [0.095 0.545];  y0 = [0.545 0.085];
    hHm = 0.295;  hStrip = 0.048;

    for q = 1:2
        for sk = 1:2
            ax = axes(fig,'Position',[x0(sk), y0(q), 0.375, hHm]); %#ok<LAXES>
            M = E{q,c};
            if isempty(M), axis(ax,'off'); continue; end
            M = M(ORD{q,c,sk}, :);
            imagesc(ax, ph, 1:size(M,1), M);
            set(ax,'YDir','normal','CLim',cl,'TickDir','out','Layer','top');
            colormap(ax, ternary(useZ, parula(256), gray(256)));
            hold(ax,'on');
            for cb0 = 1:nCycShow-1
                plot(ax, [cb0 cb0], [0.5 size(M,1)+0.5], '-', ...
                     'Color',[0.90 0.15 0.10], 'LineWidth',1.0);
            end
            hold(ax,'off');
            xlim(ax,[0 nCycShow]);  ylim(ax,[0.5 size(M,1)+0.5]);
            set(ax,'XTick',0:0.5:nCycShow);
            xlabel(ax, sprintf('phase (cycles from insp %s)', TRIG{q}));
            ylabel(ax, sprintf('cycle (sorted by phase of %s event)', lower(SKEY{sk})));

            axb = axes(fig,'Position',[x0(sk), y0(q)+hHm+0.004, 0.375, hStrip]); %#ok<LAXES>
            hold(axb,'on');
            mo = NaN;
            if ~isempty(B{q,c})
                b = mean(B{q,c},1,'omitnan');
                b = (b - min(b)) / max(max(b)-min(b), eps);
                plot(axb, ph, b, '-', 'Color',[0.25 0.25 0.25], 'LineWidth',1.2);
            end
            if ~isempty(OTH{q,c})
                mo = median(OTH{q,c});
                for cb0 = 0:nCycShow-1
                    plot(axb, [mo mo]+cb0, [-0.05 1.05], '--', ...
                         'Color',[0.35 0.75 1.00], 'LineWidth',1.0);
                end
            end
            for cb0 = 1:nCycShow-1
                plot(axb, [cb0 cb0], [-0.05 1.05], '-', ...
                     'Color',[0.90 0.15 0.10], 'LineWidth',1.0);
            end
            hold(axb,'off');
            xlim(axb,[0 nCycShow]);  ylim(axb,[-0.05 1.05]);
            set(axb,'XTick',[],'YTick',[],'Box','on');
            if sk == 1, ylabel(axb,'breath','FontSize',7.5); end
            title(axb, {sprintf('%s cycles  |  sorted by %s event', upper(TRIG{q}), SKEY{sk}), ...
                        sprintf('%d cycles with an event (%d without), median %s at phase %.2f', ...
                                size(M,1), NC(q,c), ...
                                ternary(q==1,'peak','onset'), mo)}, ...
                  'FontSize',8.5,'FontWeight','normal','Interpreter','none');
        end
    end

    cb = colorbar(ax,'Position',[0.945 0.085 0.014 0.295]);
    cb.Label.String = ternary(useZ,'dF/F (z-scored)','dF/F');

    sgtitle({sprintf('%s  |  %s cells (n = %d)  |  per-cycle dF/F on PHASE', ...
                     SETNAME, CLS{c}, nCellCls(c)), ...
             ['each cycle normalised by ITS OWN duration: 0 = the landmark, 1 = the next, ' ...
              '2 = the one after   |   red = cycle boundary'], ...
             sprintf(['SAME cycles in both columns, reordered: left by the EARLIEST event ' ...
                      'of cycle 1, right by the LATEST, both descending  |  %s  |  ' ...
                      'clim %.3g to %.3g  |  sky dashed = median other landmark'], ...
                      ternary(useZ,'z-scored, parula','raw dF/F, grayscale'), cl(1), cl(2))}, ...
            'Interpreter','tex','FontSize',9);

    stem = fullfile(outDir, sprintf('%s_%s_cycle_dff_heatmap_phase', tag, CTAG{c}));
    exportgraphics(fig,[stem '.png'],'Resolution',250,'BackgroundColor','white');
    exportgraphics(fig,[stem '.pdf'],'ContentType','image','BackgroundColor','white');

    for q = 1:2
        fprintf('  %-6s %-16s %6d cycles (%5d without)  earliest %.3f  latest %.3f\n', ...
                TRIG{q}, CLS{c}, size(E{q,c},1), NC(q,c), ...
                median(PH1{q,c}), median(PH2{q,c}));
    end
    fprintf('  -> %s.png\n', stem);
    end
end

%% ===================== LOCAL =====================
function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end
