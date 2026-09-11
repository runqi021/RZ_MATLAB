% popsel_population_figure_260831.m
% -----------------------------------------------------------------------
%  POPULATION AVERAGE for the popsel groups, six panels.
%
%      row 1  INSPIRATION-ONSET triggered
%      row 2  BREATH-PEAK       triggered
%      col 1  raw dF/F          what was recorded; set by the brightest cells
%      col 2  z-scored dF/F     every cell weighted equally; shows the SHAPE
%      col 3  event histogram   events per cycle
%
%  Mean +/- SEM ACROSS CELLS, each cell counting once no matter how many events
%  it fired or how many recordings it appeared in -- otherwise one 278-event cell
%  outweighs twenty 10-event ones.
%
%  X AXIS IS ABSOLUTE TIME IN SECONDS, +/-1.5 s. Note this is NOT the tau/IBI
%  cycle axis the popsel GUI averages on: IBI runs 0.40-2.83 s across these
%  groups, so on a seconds axis a fast breather's next breath falls inside the
%  window while a slow breather's does not, and the flanks mix different phases.
%  That is the price of real time units and it is what was asked for; the
%  cycle-axis version is the one to use for comparing transient SHAPE.
%
%  THE HISTOGRAM CANNOT FILL THE WINDOW. temporal_phase_cell_fig_260812:358
%  hardcodes the histogram half-window to the cell's own IBI and ignores
%  trigWin_sec, so a cell with IBI 0.4 s contributes histogram only out to
%  +/-0.4 s and is NaN beyond. The mean is NaN-aware, so the number of
%  contributing cells FALLS with |t| in column 3 -- the panel prints the |t| at
%  which every cell still contributes, and shades beyond it.
%
%  THE INPUT IS THE CURATION. Cells come from popsel_decisions_<GROUP>.csv,
%  decision == "include", joined to popsel_cache_<GROUP>.mat by cell number.
%  Both are written by popsel_run_260831.m and are always on the same numbering.
%
%  Runqi Zhang / 2026-08-31
% -----------------------------------------------------------------------

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));

% Software OpenGL: the NVIDIA path intermittently writes 9,223-byte blank PNGs
% out of exportgraphics, or silently leaves a stale file. See
% reference_matlab_opengl_export_stubs.
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ===================== USER-EDITABLE =====================
sumRoot = 'D:\Ventral_surface_summary';
inDir   = fullfile(sumRoot,'popsel_260816');
outDir  = fullfile(sumRoot,'popsel_population_260831');

% {group, legend label, colour, line style}.  Colour = genotype, so the same
% genotype reads the same way as on the ventral map.  Line style carries the
% manipulation: solid = intact, dashed = vagotomised.
GRP = { ...
  'IO',                  'IO',                 [0    0    0   ], '-'
  'Vglut2',              'Vglut2',             [0.10 0.65 0.20], '-'
  'Vglut2_vagotomized',  'Vglut2 vagotomised', [0.10 0.65 0.20], '--'
  'Sert',                'Sert',               [0.90 0.45 0.10], '-'
  'Sert_vagotomized',    'Sert vagotomised',   [0.90 0.45 0.10], '--'
  'Sst',                 'Sst',                [0.55 0.20 0.75], '-'
  'Vgat',                'Vgat',               [0.10 0.30 0.85], '-' };

XLIM    = [-1.5 1.5];   % seconds
showIO  = true;         % false = drop IO, ~120 cells, dominates the raw y-scale
LW      = 1.8;
SEMband = true;
%% =========================================================

if ~isfolder(outDir), mkdir(outDir); end

%% ---- gather ----------------------------------------------------------------
G = struct('name',{},'lab',{},'col',{},'sty',{},'n',{},'cells',{}, ...
           'tauS',{},'ctrS',{},'cov',{}, ...
           'dffOn',{},'dffZOn',{},'histOn',{}, ...
           'dffPk',{},'dffZPk',{},'histPk',{});
for i = 1:size(GRP,1)
    nm = GRP{i,1};
    if ~showIO && strcmp(nm,'IO'), continue; end
    cf = fullfile(inDir, sprintf('popsel_cache_%s.mat', nm));
    df = fullfile(inDir, sprintf('popsel_decisions_%s.csv', nm));
    if ~isfile(cf), fprintf(2,'%-20s no cache, skipped\n', nm); continue; end
    if ~isfile(df), fprintf(2,'%-20s not curated yet, skipped\n', nm); continue; end
    K = load(cf);  C = K.C;
    T = readtable(df,'TextType','string');

    assert(isfield(C,'dffOnsetS'), ...
        ['%s cache has no absolute-time curves -- rebuild it with ' ...
         'popsel_run_260831.m (REBUILD = true)'], nm);

    inc = T.cell(T.decision == "include");
    [tf, loc] = ismember(inc, [C.cell].');
    if any(~tf)
        fprintf(2,'%-20s %d included cell(s) not in the cache (ignored): %s\n', ...
                nm, nnz(~tf), mat2str(inc(~tf).'));
    end
    sel = loc(tf);
    if isempty(sel), fprintf(2,'%-20s no included cells, skipped\n', nm); continue; end

    % |t| out to which EVERY included cell still contributes histogram
    cov = min([C(sel).histWinHalfS]);

    G(end+1) = struct('name',nm,'lab',GRP{i,2},'col',GRP{i,3},'sty',GRP{i,4}, ...
        'n',numel(sel),'cells',[C(sel).cell], ...
        'tauS',K.tauS,'ctrS',K.ctrS,'cov',cov, ...
        'dffOn', cell2mat({C(sel).dffOnsetS}.'), ...
        'dffZOn',cell2mat({C(sel).dffOnsetZS}.'), ...
        'histOn',cell2mat({C(sel).histOnsetS}.'), ...
        'dffPk', cell2mat({C(sel).dffPeakS}.'), ...
        'dffZPk',cell2mat({C(sel).dffPeakZS}.'), ...
        'histPk',cell2mat({C(sel).histPeakS}.')); %#ok<SAGROW>
    fprintf('%-20s %3d included of %3d   histogram covers all cells to +/-%.2f s\n', ...
            nm, numel(sel), numel(C), cov);
end
assert(~isempty(G), 'no group had an included cell -- curate something first');

%% ---- figure ----------------------------------------------------------------
fig = figure('Color','w','Units','pixels','Position',[60 60 1500 860],'Visible','off');
tl  = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% {field, ylabel, title, x-grid field}
PAN = { 'dffOn','mean dF/F','raw dF/F','tauS'
        'dffZOn','mean dF/F (z-scored)','z-scored dF/F','tauS'
        'histOn','events per cycle','event histogram','ctrS'
        'dffPk','mean dF/F','raw dF/F','tauS'
        'dffZPk','mean dF/F (z-scored)','z-scored dF/F','tauS'
        'histPk','events per cycle','event histogram','ctrS' };
ROWLAB = {'inspiration onset','breath peak'};

for p = 1:6
    ax = nexttile(tl,p); hold(ax,'on');
    row = ceil(p/3);
    isHist = strcmp(PAN{p,4},'ctrS');
    for k = 1:numel(G)
        M = G(k).(PAN{p,1});
        x = G(k).(PAN{p,4});
        mu = mean(M,1,'omitnan');
        sd = std(M,0,1,'omitnan');
        nn = sum(~isnan(M),1);
        se = sd ./ max(sqrt(nn),1);
        if SEMband
            good = ~isnan(mu) & nn > 1;
            if any(good)
                fill(ax,[x(good) fliplr(x(good))],[mu(good)+se(good) fliplr(mu(good)-se(good))], ...
                     G(k).col,'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
            end
        end
        plot(ax,x,mu,'Color',G(k).col,'LineStyle',G(k).sty,'LineWidth',LW, ...
             'DisplayName',sprintf('%s (n=%d)',G(k).lab,G(k).n));
    end
    yl = ylim(ax);
    % The trigger sits at 0 by construction. Red for onset, blue for peak, the
    % same colours temporal_phase_cell_fig_260812 uses on the per-cell figures.
    tc = [0.90 0.10 0.10]; if row == 2, tc = [0.20 0.50 0.95]; end
    plot(ax,[0 0],yl,'-','Color',tc,'LineWidth',1,'HandleVisibility','off');

    % Shade where the histogram no longer has every cell behind it.
    if isHist
        cmin = min([G.cov]);
        if isfinite(cmin) && cmin < XLIM(2)
            for s = [-1 1]
                xa = sort([s*cmin, s*XLIM(2)]);
                fill(ax,[xa(1) xa(2) xa(2) xa(1)],[yl(1) yl(1) yl(2) yl(2)], ...
                     [0.5 0.5 0.5],'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
            end
            text(ax,0,yl(2),sprintf('all cells only within \\pm%.2f s',cmin), ...
                 'HorizontalAlignment','center','VerticalAlignment','top', ...
                 'FontSize',8,'Color',[0.35 0.35 0.35]);
        end
    end
    ylim(ax,yl); xlim(ax,XLIM);
    xlabel(ax,sprintf('time from %s (s)',ROWLAB{row}));
    ylabel(ax,PAN{p,2});
    title(ax,sprintf('%s  -  %s triggered',PAN{p,3},ROWLAB{row}),'FontWeight','normal');
    grid(ax,'on'); ax.GridAlpha = 0.12; box(ax,'off');
    if p == 3, legend(ax,'Location','northeast','Box','off','FontSize',8); end
end

title(tl, sprintf(['population average, absolute time   |   %d groups, %d cells   ' ...
    '|   top row = inspiration-onset triggered, bottom row = breath-peak triggered'], ...
    numel(G), sum([G.n])), 'FontWeight','normal');

%% ---- save ------------------------------------------------------------------
stem = fullfile(outDir,'popsel_population_onset_peak');
exportgraphics(fig,[stem '.png'],'Resolution',200,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','none');

T = table(string({G.name}.'), string({G.lab}.'), [G.n].', [G.cov].', ...
    'VariableNames',{'group','label','n_cells','hist_full_coverage_s'});
writetable(T, fullfile(outDir,'popsel_population_groups.csv'));

rows = {};
for k = 1:numel(G)
    for c = G(k).cells, rows(end+1,:) = {G(k).name, c}; end %#ok<SAGROW>
end
writetable(cell2table(rows,'VariableNames',{'group','cell'}), ...
           fullfile(outDir,'popsel_population_cells.csv'));

d = dir([stem '.png']);
fprintf('\nsaved -> %s\n', outDir);
fprintf('  popsel_population_onset_peak.png  (%.0f KB)\n', d.bytes/1024);
if d.bytes < 20000
    fprintf(2,'  PNG is suspiciously small -- check it is not an OpenGL blank stub\n');
end
close(fig);
disp(T);
