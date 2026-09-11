% cycle_mean_phase_chat_260902.m
% -----------------------------------------------------------------------
%  Mean +/- SD of the per-cycle dF/F, on the same phase axis and from the same
%  cycles as the heatmaps, for two ChAT groups, SIGNIFICANT cells only:
%      vagal-intact (0521, 0522, 0523)
%      vagotomised (0826)
%
%  READ THE n BEFORE READING ANYTHING ELSE. ChAT contributes 5 active cells in
%  total and 3 of them are significant: TWO intact (0521 cell018, 0522 cell020)
%  and ONE vagotomised (0826 cell023). The vagotomised "population average" is a
%  single cell, and the intact one is two. The SD band is across CYCLES, not
%  across cells, so it says nothing about whether either curve generalises to
%  other ChAT neurons -- with n = 1 and n = 2 it cannot. Treat this figure as two
%  worked examples, not as a group comparison.
%
%      rows      the three groups, one per row
%      columns   inspiration ONSET cycles  |  inspiratory PEAK cycles
%      6 panels, ONE GROUP PER PANEL -- the groups are not overlaid
%
%  Every panel shares one y scale, so the panels can still be read against each
%  other; that is what the overlay was for and it is the only thing lost.
%
%  SD, NOT SEM. The band is the spread of the individual cycles, which is the
%  quantity the heatmap shows row by row. An SEM band over thousands of cycles
%  would be a hairline and would say only that the mean is well estimated, which
%  is not in doubt and is not the question.
%
%  THE CYCLES ARE EXACTLY THE HEATMAPS' CYCLES: significant cells only,
%  non-overlapping blocks of nCycShow cycles (so every breath is counted once),
%  cycle 1 must contain an event, and blocks with a cycle longer than
%  maxCycFac x the recording median are dropped. Change any of those here and
%  the two figures stop describing the same data.
%
%  PHASE, each cycle normalised by ITS OWN duration: 0 = the landmark, 1 and 2 =
%  the following ones. Averaging in phase rather than in seconds is what lets
%  three groups whose breathing rates differ by 3x be drawn on one axis at all.
%
%  SOLID = INTACT, DASHED = VAGOTOMISED, the manipulation convention from
%  popsel_population_figure_260831. The two intact groups are separated by
%  lightness, not hue: they are the same genotype.
%
%  THE BREATH IS OVERLAID ON THE SAME AXES, on a right-hand scale, one grey curve
%  per group with that group's linestyle. It is min-max scaled and its axis is
%  unlabelled on purpose: the trace is SVD PC1, whose units and sign are
%  arbitrary per video, so only its TIMING carries meaning. Each recording's sign
%  is fixed on its own peak-triggered average before averaging, or the pooled
%  waveform partly cancels (project_breath_sign_is_analysis_immune).
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

% {display name, dates, mode, class ('sig'|'all'), colour, linestyle}
GRP = { 'ChAT intact, SIGNIFICANT',     {'0826'}, 'exclude', 'sig', [0.85 0.10 0.10], '-'
        'ChAT vagotomised, SIGNIFICANT', {'0826'}, 'include', 'sig', [0.85 0.10 0.10], '--' };
GENOTYPE  = 'ChAT';
nCycShow  = 3;
strideCyc = 3;           % non-overlapping, as in the heatmap
nPh       = 301;
maxCycFac = 3;
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

CU = ca_recheck_load_curation_260901(bundleDir);
Ts = readtable(sigCsv,'TextType','string');
sigMap = containers.Map(cellstr(Ts.stem), num2cell(logical(Ts.sig_perm)));

TRIG = {'onset','peak'};
ph   = linspace(0, nCycShow, nPh);
nG   = size(GRP,1);
MU = cell(2,nG);  SD = cell(2,nG);  NN = zeros(2,nG);  NCELL = zeros(1,nG);
BR = cell(2,nG);

for gi = 1:nG
    DATES = GRP{gi,2};  DMODE = GRP{gi,3};  CLASS = GRP{gi,4};
    d = dir(fullfile(bundleDir,'cells',[GENOTYPE '_*.mat']));
    hit = false(1,numel(d));
    for v = 1:numel(DATES), hit = hit | contains({d.name}, ['_' DATES{v} '_']); end
    if strcmpi(DMODE,'include'), d = d(hit); else, d = d(~hit); end

    E = {[],[]};  B = {[],[]};
    seenRec = containers.Map('KeyType','char','ValueType','logical');
    for k = 1:numel(d)
        stem = erase(d(k).name,'.mat');
        if strcmpi(CLASS,'sig') && ~(isKey(sigMap,stem) && sigMap(stem))
            continue                                   % this group takes sig cells only
        end
        L = load(fullfile(d(k).folder,d(k).name),'RECc');
        used = false;
        for i = 1:numel(L.RECc)
            r = L.RECc(i);  g = CU.get(r.folder, r.roi);
            if g.found && (g.toss || g.cell_toss), continue; end
            ev = r.spike_idx(:);
            if g.found, ev = g.spike_idx(:); end
            x  = r.dff;
            bz = (r.breath - mean(r.breath)) / max(std(r.breath), eps);
            pkT = r.peak(r.peak >= 1 & r.peak <= r.T);
            if ~isempty(pkT) && mean(bz(pkT)) < 0, bz = -bz; end

            for q = 1:2
                if q == 1, trg = r.foot(:); else, trg = r.peak(:); end
                if numel(trg) < nCycShow+1, continue; end
                used = true;
                dt = diff(trg);  med = median(dt);
                Mi = [];  Bi = [];
                for cix = 1:strideCyc:numel(dt)-(nCycShow-1)
                    if any(dt(cix:cix+nCycShow-1) > maxCycFac*med), continue; end
                    a = trg(cix);  b = trg(cix+1);
                    if a < 1 || trg(cix+nCycShow) > r.T, continue; end
                    if ~any(ev >= a & ev < b), continue; end     % cycle 1 must fire
                    uu = [];  vv = [];  bb = [];
                    for cj = 0:nCycShow-1
                        s0 = trg(cix+cj);  s1 = trg(cix+cj+1);
                        idx = (s0:s1)';  qj = cj + (idx - s0)/(s1 - s0);
                        if cj == 0
                            uu = qj;  vv = x(idx);  bb = bz(idx);
                        else
                            uu = [uu; qj(2:end)];      %#ok<AGROW>
                            vv = [vv; x(idx(2:end))];  %#ok<AGROW>
                            bb = [bb; bz(idx(2:end))]; %#ok<AGROW>
                        end
                    end
                    Mi(end+1,:) = interp1(uu, vv, ph, 'linear', NaN); %#ok<AGROW>
                    Bi(end+1,:) = interp1(uu, bb, ph, 'linear', NaN); %#ok<AGROW>
                end
                if isempty(Mi), continue; end
                E{q} = [E{q}; single(Mi)];
                rkey = sprintf('%s|%d', r.folder, q);
                if ~isKey(seenRec, rkey)
                    seenRec(rkey) = true;
                    B{q} = [B{q}; mean(Bi,1,'omitnan')];
                end
            end
        end
        if used, NCELL(gi) = NCELL(gi) + 1; end
    end
    for q = 1:2
        MU{q,gi} = mean(double(E{q}),1,'omitnan');
        SD{q,gi} = std(double(E{q}),0,1,'omitnan');
        NN(q,gi) = size(E{q},1);
        BR{q,gi} = mean(B{q},1,'omitnan');
    end
    fprintf('%-30s %3d cells (%s) | onset %4d cycles | peak %4d cycles\n', ...
            GRP{gi,1}, NCELL(gi), CLASS, NN(1,gi), NN(2,gi));
end

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 17], ...
             'Name','ChAT: mean per-cycle dF/F on phase');
set(fig,'DefaultAxesFontSize',9);

% one y scale for all six panels
yAll = [];
for q=1:2, for gi=1:nG
    yAll = [yAll, MU{q,gi}+SD{q,gi}, MU{q,gi}-SD{q,gi}]; %#ok<AGROW>
end, end
yl = [min(yAll) max(yAll)];  yl = yl + [-1 1]*0.05*diff(yl);

x0 = [0.095 0.545];
y0 = [0.505 0.085];                % one row per group, top to bottom
for gi = 1:nG
    for q = 1:2
        ax = axes(fig,'Position',[x0(q), y0(gi), 0.375, 0.345]); %#ok<LAXES>

        % this group's own breath, right-hand scale, drawn under the calcium
        yyaxis(ax,'right');
        hold(ax,'on');
        b = BR{q,gi};
        b = (b - min(b)) / max(max(b)-min(b), eps);
        plot(ax, ph, b, '-', 'Color',[0.60 0.60 0.60], 'LineWidth',1.1);
        ylim(ax,[-0.05 1.45]);  ax.YColor = [0.60 0.60 0.60];
        if q == 2, ylabel(ax,'breath'); else, set(ax,'YTickLabel',[]); end

        yyaxis(ax,'left');
        hold(ax,'on'); box(ax,'on');
        m = MU{q,gi};  sd = SD{q,gi};  ok = isfinite(m) & isfinite(sd);
        fill(ax, [ph(ok) fliplr(ph(ok))], [m(ok)+sd(ok) fliplr(m(ok)-sd(ok))], ...
             GRP{gi,5}, 'FaceAlpha',0.20, 'EdgeColor','none');
        plot(ax, ph, m, GRP{gi,6}, 'Color',GRP{gi,5}, 'LineWidth',2.0);
        for cb0 = 1:nCycShow-1
            plot(ax, [cb0 cb0], yl, '-', 'Color',[0.90 0.15 0.10], 'LineWidth',0.8);
        end
        hold(ax,'off');
        xlim(ax,[0 nCycShow]); ylim(ax, yl); set(ax,'XTick',0:0.5:nCycShow);
        ax.YColor = [0 0 0];
        if gi == nG
            xlabel(ax, sprintf('phase (cycles from insp %s)', TRIG{q}));
        else
            set(ax,'XTickLabel',[]);
        end
        if q == 1, ylabel(ax,'dF/F'); else, set(ax,'YTickLabel',[]); end
        title(ax, sprintf('%s  |  %s  |  %d cells, %d cycles', ...
                          GRP{gi,1}, upper(TRIG{q}), NCELL(gi), NN(q,gi)), ...
              'FontSize',8.5,'FontWeight','normal','Interpreter','none');
    end
end

annotation(fig,'textbox',[0.02 0.955 0.96 0.035],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',9.5,'FontWeight','bold', ...
    'String','ChAT SIGNIFICANT cells  |  mean per-cycle dF/F \pm 1 SD across cycles');
annotation(fig,'textbox',[0.02 0.925 0.96 0.030],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',8.5,'String',['TOP = intact (n = 2 cells).  BOTTOM = vagotomised (n = 1 cell).  ' ...
    'The SD band is across CYCLES, not cells -- these are worked examples, not a group comparison']);
annotation(fig,'textbox',[0.02 0.898 0.96 0.030],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','tex', ...
    'FontSize',8.5,'String',sprintf(['%d cycles each normalised by its own duration, shared y ' ...
    'scale  |  red = cycle boundary  |  grey = that group''s mean breath, right axis, ' ...
    'timing only'], nCycShow));

stem = fullfile(outDir,'ChAT_cycle_mean_phase');
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');

T = table(ph(:));  T.Properties.VariableNames = {'phase'};
for q = 1:2
    for gi = 1:nG
        nm = matlab.lang.makeValidName(sprintf('%s_%s', TRIG{q}, GRP{gi,1}));
        T.([char(nm) '_mean']) = MU{q,gi}(:);
        T.([char(nm) '_sd'])   = SD{q,gi}(:);
    end
end
writetable(T, [stem '.csv']);

fprintf('\n%-28s %8s %10s %10s %10s\n','group','trigger','peak dF/F','at phase','mean SD');
for q = 1:2
    for gi = 1:nG
        [pk,ip] = max(MU{q,gi});
        fprintf('%-28s %8s %10.4f %10.3f %10.4f\n', GRP{gi,1}, TRIG{q}, pk, ph(ip), ...
                mean(SD{q,gi},'omitnan'));
    end
end
fprintf('\nwrote %s.png / .pdf / .csv\n', stem);
