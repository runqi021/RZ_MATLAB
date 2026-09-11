% cycle_overlay_vglut2_sig_260901.m
% -----------------------------------------------------------------------
%  Onset-triggered cycle waveform for every SIGNIFICANT Vglut2 cell, all
%  overlaid: vagal-intact beside vagotomised. Two panels.
%
%  One thin line per CELL = that cell's inspiration-onset-triggered average
%  dF/F, pooled over its kept recordings. One thick line = the mean over cells,
%  with a SEM band. Grey = the trigger-averaged breathing waveform, min-max
%  normalised on its own right-hand axis.
%
%  Z-SCORED WITHIN EACH RECORDING, BEFORE ANY POOLING. dF/F amplitude varies by
%  an order of magnitude between cells, so an overlay of raw traces is a picture
%  of which cells are brightest, and a mean over them is set by the two or three
%  loudest. Set useZ = false for raw dF/F if that is what is wanted.
%
%  X IS ABSOLUTE TIME, +/-1.5 s, the same axis the popsel population figures use,
%  so this figure can be read against them. It is wider than one cycle on
%  purpose: the neighbouring inspirations show up as the flanking humps, and how
%  fast they wash out is itself the read-out of how regular the breathing is.
%
%  DON'T READ AMPLITUDE OFF THE GREY CURVE. The breath trace is SVD PC1, whose
%  units and SIGN are arbitrary per video, so it is min-max scaled over the
%  window and averaged per RECORDING (not per cell) -- a FOV with eight cells
%  must not outvote one with a single cell.
%
%  Alignment is loaded from the re-check bundle, which reproduces
%  load_obs_local's rules exactly (nDrop on breath only, the Vglut2/1124
%  one-frame fix, T = min(dFF, breath, events)), and tossed recordings and cells
%  are dropped.
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
useZ      = true;        % z-score dF/F within each recording before pooling
% =========================================================

COL  = [0.10 0.65 0.20];         % Vglut2 green
BCOL = [0.45 0.45 0.45];         % breath grey
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

%% ===================== PER-CELL WAVEFORMS =====================
W = nan(height(Tsig), numel(tauS));     % one row per cell
B = [];                                 % one row per RECORDING, breath
Bvag = [];
nRecUsed = 0;
for k = 1:height(Tsig)
    f = fullfile(bundleDir,'cells',[char(Tsig.stem(k)) '.mat']);
    if ~isfile(f), continue; end
    L = load(f,'RECc');  REC = L.RECc;
    E = [];                              % epochs for this cell, on tauS
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
        E = [E; interp1(tk, Ei.', tauS, 'linear', NaN).']; %#ok<AGROW>

        % breath, once per RECORDING (not per cell)
        bwz = (r.breath - mean(r.breath)) / max(std(r.breath), eps);
        Eb  = epochs_local(bwz, ft, w, r.T);
        if ~isempty(Eb)
            bmu = interp1(tk, mean(Eb,1,'omitnan'), tauS, 'linear', NaN);
            if isVag(k), Bvag = [Bvag; bmu]; else, B = [B; bmu]; end %#ok<AGROW>
            nRecUsed = nRecUsed + 1;
        end
    end
    if ~isempty(E), W(k,:) = mean(E,1,'omitnan'); end
end
ok = any(isfinite(W),2);
fprintf('waveforms for %d cells over %d recording-contributions\n', nnz(ok), nRecUsed);

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[2 3 26 11], ...
             'Name','Vglut2 significant: onset-triggered cycle waveforms');
set(fig,'DefaultAxesFontSize',9);

GRP = { 'vagal intact', ~isVag & ok, B
        'vagotomised',   isVag & ok, Bvag };
yl = [min(W(ok,:),[],'all','omitnan'), max(W(ok,:),[],'all','omitnan')];
yl = yl + [-1 1]*0.05*max(diff(yl),eps);

for gi = 1:2
    m = GRP{gi,2};
    ax = axes(fig,'Position',[0.075+(gi-1)*0.470, 0.155, 0.375, 0.60]); %#ok<LAXES>
    hold(ax,'on'); box(ax,'on');

    % breath first, on its own right-hand axis, so the calcium sits on top
    yyaxis(ax,'right');
    bb = GRP{gi,3};
    if ~isempty(bb)
        b = mean(bb,1,'omitnan');
        b = (b - min(b)) / max(max(b)-min(b), eps);
        plot(ax, tauS, b, '-', 'Color',BCOL, 'LineWidth',1.3);
    end
    ylim(ax,[-0.05 1.35]);
    ax.YColor = BCOL;
    if gi == 2, ylabel(ax,'breath (norm.)'); else, set(ax,'YTickLabel',[]); end

    yyaxis(ax,'left');
    for k = find(m).'
        plot(ax, tauS, W(k,:), '-', 'Color',[COL 0.35], 'LineWidth',0.7);
    end
    mu = mean(W(m,:),1,'omitnan');
    nn = sum(isfinite(W(m,:)),1);
    se = std(W(m,:),0,1,'omitnan') ./ max(sqrt(nn),1);
    good = isfinite(mu) & nn > 1;
    fill(ax, [tauS(good) fliplr(tauS(good))], [mu(good)+se(good) fliplr(mu(good)-se(good))], ...
         COL, 'FaceAlpha',0.20, 'EdgeColor','none');
    plot(ax, tauS, mu, '-', 'Color',COL*0.6, 'LineWidth',2.2);
    plot(ax, [0 0], yl, 'k-', 'LineWidth',0.8);
    ylim(ax, yl);
    ax.YColor = [0 0 0];
    xlim(ax, [-WIN_S WIN_S]);
    xlabel(ax,'time from inspiration ONSET (s)');
    if gi == 1
        ylabel(ax, ternary(useZ,'dF/F (z-scored)','dF/F'));
    else
        set(ax,'YTickLabel',[]);
    end
    title(ax, sprintf('%s %s  |  %d cells', GENOTYPE, GRP{gi,1}, nnz(m)), ...
          'Interpreter','none','FontSize',9,'FontWeight','normal');
end

sgtitle({sprintf('%s SIGNIFICANT cells: %d of %d   (permutation p \\leq 0.01, peak-triggered PSTH, 1200 shifts)', ...
                 GENOTYPE, nnz(ok), nGenoTot), ...
         sprintf(['inspiration ONSET triggered   |   thin = one cell, thick = mean ' ...
                  '\\pm SEM over cells, grey = trigger-averaged breath (min-max, ' ...
                  'right axis)']), ...
         ternary(useZ, ['dF/F z-scored WITHIN each recording before pooling, so the ' ...
                        'overlay is not a picture of which cells are brightest'], ...
                       'raw dF/F')}, ...
        'Interpreter','tex','FontSize',8.5);

stem = fullfile(outDir, sprintf('%s_sig_onset_cycle_overlay', GENOTYPE));
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');

Tout = array2table([tauS(:), mean(W(GRP{1,2},:),1,'omitnan').', ...
                              mean(W(GRP{2,2},:),1,'omitnan').'], ...
    'VariableNames',{'time_s','intact_mean','vagotomised_mean'});
writetable(Tout, [stem '.csv']);

for gi = 1:2
    m = GRP{gi,2};
    mu = mean(W(m,:),1,'omitnan');
    [pk,ip] = max(mu);
    fprintf('%-14s %2d cells   peak %.3f at %+.0f ms\n', ...
            GRP{gi,1}, nnz(m), pk, 1000*tauS(ip));
end
fprintf('\nwrote %s.png / .pdf / .csv\n', stem);

%% ===================== LOCAL =====================
function E = epochs_local(x, idx, w, T)
%EPOCHS_LOCAL  Rows of x cut around every idx that fits a full window.
idx = idx(idx-w >= 1 & idx+w <= T);
E = [];
if isempty(idx), return; end
E = zeros(numel(idx), 2*w+1);
for k = 1:numel(idx), E(k,:) = x(idx(k)-w : idx(k)+w); end
end

function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end
