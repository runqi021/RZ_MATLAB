% cycle_traces_vglut2_sig_260901.m
% -----------------------------------------------------------------------
%  Every individual breath cycle, overlaid. Onset-triggered dF/F for the
%  SIGNIFICANT Vglut2 cells, vagal-intact beside vagotomised.
%
%  One thin line per CYCLE -- not per cell. Every inspiration onset in every kept
%  recording of every selected cell contributes its own dF/F snippet, so the
%  panel shows the raw variability the cell-averaged version hides: how many
%  cycles actually carry a transient, how much the amplitude scatters, and
%  whether the timing is tight or smeared.
%
%  RAW dF/F by default, because "the cycle dF/F traces" are these. Be aware of
%  what that means: a cell with big transients contributes both more amplitude
%  AND, if it fired more, more lines, so the cloud is weighted towards the loud
%  cells twice over. useZ = true z-scores within each recording first, which
%  removes the amplitude half of that.
%
%  DRAWN AS ONE LINE OBJECT PER PANEL, with NaN separators between cycles.
%  Thousands of separate line objects would take minutes to render and would make
%  the vector PDF enormous; a single NaN-separated polyline draws in one pass.
%
%  NOTHING IS AVERAGED. No mean line, no SEM band, no breath curve -- the panel
%  is the cycles and the t = 0 marker, and any structure you see is structure in
%  the cloud itself. showMean and showBreath at the top put those back.
%
%  Runqi Zhang / 2026-09-01
% -----------------------------------------------------------------------
clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot,'analysis_260806'));

%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
bundleDir = fullfile(sumRoot,'per-cell-summary_active_260812','spike_recheck_260901');
outDir    = fullfile(sumRoot,'polar_recheck_260901');
sigCsv    = fullfile(outDir,'polar_recheck_perm_percell.csv');

GENOTYPE  = 'Vglut2';
VAG_DATES = {'0824'};
WIN_S     = 1.5;         % half-window, seconds
useZ      = false;       % true = z-score within each recording before pooling
lineAlpha = 0.035;       % per-cycle line opacity; density reads as shading
clipPct   = [0.2 99.8];  % y limits, percentile of all plotted samples
showMean  = false;       % thick mean-over-cycles line on top of the cloud
showBreath= false;       % trigger-averaged breath on the right axis
% =========================================================

COL  = [0.10 0.65 0.20];
BCOL = [0.45 0.45 0.45];
tauS = linspace(-WIN_S, WIN_S, 181);

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

%% ===================== SELECT =====================
assert(isfile(sigCsv), 'Run polar_recheck_perm_260901.m first (missing %s)', sigCsv);
Tsig = readtable(sigCsv,'TextType','string');
nGenoTot = nnz(Tsig.group == GENOTYPE);
Tsig = Tsig(Tsig.group == GENOTYPE & logical(Tsig.sig_perm), :);
isVag = false(height(Tsig),1);
for v = 1:numel(VAG_DATES)
    isVag = isVag | contains(Tsig.stem, "_" + VAG_DATES{v} + "_");
end
fprintf('%s significant: %d  (%d intact, %d vagotomised)\n', ...
        GENOTYPE, height(Tsig), nnz(~isVag), nnz(isVag));

CU = ca_recheck_load_curation_260901(bundleDir);

%% ===================== EVERY CYCLE =====================
E = {[], []};      % all cycle traces, group 1 = intact, 2 = vagotomised
Bg = {[], []};     % breath, one row per recording
nCellG = [0 0];
for k = 1:height(Tsig)
    f = fullfile(bundleDir,'cells',[char(Tsig.stem(k)) '.mat']);
    if ~isfile(f), continue; end
    L = load(f,'RECc');  REC = L.RECc;
    gi = 1 + double(isVag(k));
    used = false;
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && (g.toss || g.cell_toss), continue; end
        ft = r.foot(:);
        if numel(ft) < 2, continue; end
        x = r.dff;
        if useZ, x = (x - mean(x,'omitnan')) / max(std(x,'omitnan'), eps); end
        w  = max(1, ceil(WIN_S*r.fps) + 1);
        tk = (-w:w)/r.fps;

        Ei = epochs_local(x, ft, w, r.T);
        if isempty(Ei), continue; end
        E{gi} = [E{gi}; interp1(tk, Ei.', tauS, 'linear', NaN).'];
        used = true;

        Eb = epochs_local((r.breath-mean(r.breath))/max(std(r.breath),eps), ft, w, r.T);
        if ~isempty(Eb)
            Bg{gi} = [Bg{gi}; interp1(tk, mean(Eb,1,'omitnan'), tauS, 'linear', NaN)];
        end
    end
    if used, nCellG(gi) = nCellG(gi) + 1; end
end
fprintf('cycles: %d intact, %d vagotomised\n', size(E{1},1), size(E{2},1));

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[2 3 26 11], ...
             'Name','Vglut2 significant: every onset-triggered cycle');
set(fig,'DefaultAxesFontSize',9);

allv = [E{1}(:); E{2}(:)];  allv = allv(isfinite(allv));
yl = prctile(allv, clipPct);
yl = yl + [-1 1]*0.04*max(diff(yl),eps);

NAME = {'vagal intact','vagotomised'};
for gi = 1:2
    ax = axes(fig,'Position',[0.075+(gi-1)*0.470, 0.155, 0.375, 0.60]); %#ok<LAXES>
    hold(ax,'on'); box(ax,'on');

    if showBreath
        yyaxis(ax,'right');
        if ~isempty(Bg{gi})
            b = mean(Bg{gi},1,'omitnan');
            b = (b - min(b)) / max(max(b)-min(b), eps);
            plot(ax, tauS, b, '-', 'Color',BCOL, 'LineWidth',1.3);
        end
        ylim(ax,[-0.05 1.35]);  ax.YColor = BCOL;
        if gi == 2, ylabel(ax,'breath (norm.)'); else, set(ax,'YTickLabel',[]); end
        yyaxis(ax,'left');
    end
    M = E{gi};
    if ~isempty(M)
        % ONE polyline for every cycle, NaN-separated
        X = repmat([tauS, NaN], size(M,1), 1).';
        Y = [M, nan(size(M,1),1)].';
        plot(ax, X(:), Y(:), '-', 'Color',[COL lineAlpha], 'LineWidth',0.4);
        if showMean
            plot(ax, tauS, mean(M,1,'omitnan'), '-', ...
                 'Color',[0.02 0.30 0.08], 'LineWidth',2.4);
        end
    end
    plot(ax, [0 0], yl, 'k-', 'LineWidth',0.8);
    ylim(ax, yl);  ax.YColor = [0 0 0];
    xlim(ax, [-WIN_S WIN_S]);
    xlabel(ax,'time from inspiration ONSET (s)');
    if gi == 1
        ylabel(ax, ternary(useZ,'dF/F (z-scored)','dF/F'));
    else
        set(ax,'YTickLabel',[]);
    end
    title(ax, sprintf('%s %s  |  %d cells  |  %d cycles', ...
                      GENOTYPE, NAME{gi}, nCellG(gi), size(M,1)), ...
          'Interpreter','none','FontSize',9,'FontWeight','normal');
end

sgtitle({sprintf('%s SIGNIFICANT cells: %d of %d   (permutation p \\leq 0.01, peak-triggered PSTH, 1200 shifts)', ...
                 GENOTYPE, sum(nCellG), nGenoTot), ...
         'inspiration ONSET triggered   |   ONE THIN LINE PER CYCLE, nothing averaged', ...
         sprintf(['%s   |   y limits = %.4g to %.4g (%.3g-%.3g percentile of all plotted ' ...
                  'samples)'], ternary(useZ,'z-scored within each recording','raw dF/F'), ...
                  yl(1), yl(2), clipPct(1), clipPct(2))}, ...
        'Interpreter','tex','FontSize',8.5);

stem = fullfile(outDir, sprintf('%s_sig_onset_cycle_traces', GENOTYPE));
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');
dpng = dir([stem '.png']);  dpdf = dir([stem '.pdf']);
fprintf('png %.1f MB, pdf %.1f MB\n', dpng.bytes/1e6, dpdf.bytes/1e6);

for gi = 1:2
    mu = mean(E{gi},1,'omitnan');
    [pk,ip] = max(mu);
    fprintf('%-14s %2d cells %5d cycles   mean peak %.3f at %+.0f ms\n', ...
            NAME{gi}, nCellG(gi), size(E{gi},1), pk, 1000*tauS(ip));
end
fprintf('\nwrote %s.png / .pdf\n', stem);

%% ===================== LOCAL =====================
function E = epochs_local(x, idx, w, T)
idx = idx(idx-w >= 1 & idx+w <= T);
E = [];
if isempty(idx), return; end
E = zeros(numel(idx), 2*w+1);
for k = 1:numel(idx), E(k,:) = x(idx(k)-w : idx(k)+w); end
end

function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end
