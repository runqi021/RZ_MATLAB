% chat_sig_heatmap_260902.m
% -----------------------------------------------------------------------
%  ONE per-cycle dF/F heatmap for ChAT: every SIGNIFICANT cell, intact and
%  vagotomised POOLED INTO THE SAME HEATMAP -- not split into columns.
%
%      ONSET cycles  |  PEAK cycles, two panels of one pooled heatmap each
%
%  READ THE n FIRST. ChAT contributes 5 active cells and 3 are significant:
%  TWO intact (0521 cell018, 0522 cell020) and ONE vagotomised (0826 cell023).
%  The right-hand column is one neuron. Two worked examples against one, not a
%  group comparison.
%
%  ChAT/0826 IS THE VAGOTOMISED SESSION -- source folder
%  C:\Users\Admin\Desktop\260826_ChAT-soma-g8s_vagotomized.
%
%  Everything else follows the other cycle figures exactly, so the panels can be
%  read against them: single-landmark phase with each cycle normalised by its own
%  duration, nCycShow cycles laid end to end, non-overlapping blocks so every
%  breath is counted once, cycle 1 must contain an event, blocks holding a cycle
%  longer than maxCycFac x the recording median dropped, rows sorted by the phase
%  of the earliest event descending, raw dF/F in grayscale at
%  prctile([0.5 99.5]), and the breath sign-corrected per recording before
%  averaging (SVD PC1's sign is arbitrary per video).
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

GENOTYPE  = 'ChAT';
VAG_DATES = {'0826'};
COL       = [0.85 0.10 0.10];      % ChAT red
nCycShow  = 3;
strideCyc = 3;
nPh       = 301;
maxCycFac = 3;
climPct   = [0.5 99.5];
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

CU = ca_recheck_load_curation_260901(bundleDir);
Ts = readtable(sigCsv,'TextType','string');
sigMap = containers.Map(cellstr(Ts.stem), num2cell(logical(Ts.sig_perm)));

GRPN = {'intact','vagotomised'};
TRIG = {'onset','peak'};
ph   = linspace(0, nCycShow, nPh);

E = cell(2,1);  P1 = cell(2,1);  B = cell(2,1);
for q=1:2, E{q,1}=[]; P1{q,1}=[]; B{q,1}=[]; end
NCELL = 0;  CELLS = {};  NVAG = 0;

d = dir(fullfile(bundleDir,'cells',[GENOTYPE '_*.mat']));
for k = 1:numel(d)
    stem = erase(d(k).name,'.mat');
    if ~(isKey(sigMap,stem) && sigMap(stem)), continue; end       % SIG only
    isVag = false;
    for v = 1:numel(VAG_DATES), isVag = isVag || contains(stem, ['_' VAG_DATES{v} '_']); end
    g = 1;                       % pooled: intact and vagotomised share one heatmap
    L = load(fullfile(d(k).folder,d(k).name),'RECc');
    used = false;
    for i = 1:numel(L.RECc)
        r = L.RECc(i);  gg = CU.get(r.folder, r.roi);
        if gg.found && (gg.toss || gg.cell_toss), continue; end
        ev = r.spike_idx(:);
        if gg.found, ev = gg.spike_idx(:); end
        x  = r.dff;
        bz = (r.breath - mean(r.breath)) / max(std(r.breath), eps);
        pkT = r.peak(r.peak >= 1 & r.peak <= r.T);
        if ~isempty(pkT) && mean(bz(pkT)) < 0, bz = -bz; end

        for q = 1:2
            if q == 1, trg = r.foot(:); else, trg = r.peak(:); end
            if numel(trg) < nCycShow+1, continue; end
            used = true;
            dt = diff(trg);  med = median(dt);
            Mi = []; p1 = []; Bi = [];
            for cix = 1:strideCyc:numel(dt)-(nCycShow-1)
                if any(dt(cix:cix+nCycShow-1) > maxCycFac*med), continue; end
                a = trg(cix);  b = trg(cix+1);
                if a < 1 || trg(cix+nCycShow) > r.T, continue; end
                e = ev(ev >= a & ev < b);
                if isempty(e), continue; end            % cycle 1 must fire
                uu=[]; vv=[]; bb=[];
                for cj = 0:nCycShow-1
                    s0 = trg(cix+cj); s1 = trg(cix+cj+1);
                    idx = (s0:s1)';  qj = cj + (idx-s0)/(s1-s0);
                    if cj==0, uu=qj; vv=x(idx); bb=bz(idx);
                    else, uu=[uu;qj(2:end)]; vv=[vv;x(idx(2:end))]; bb=[bb;bz(idx(2:end))]; end %#ok<AGROW>
                end
                Mi(end+1,:) = interp1(uu,vv,ph,'linear',NaN); %#ok<AGROW>
                Bi(end+1,:) = interp1(uu,bb,ph,'linear',NaN); %#ok<AGROW>
                p1(end+1,1) = (min(e)-a)/(b-a);              %#ok<AGROW>
            end
            if isempty(Mi), continue; end
            E{q,g}  = [E{q,g};  Mi];
            P1{q,g} = [P1{q,g}; p1];
            B{q,g}  = [B{q,g};  mean(Bi,1,'omitnan')];
        end
    end
    if used, NCELL = NCELL + 1; CELLS{end+1} = stem; NVAG = NVAG + double(isVag); end
end

for q=1:2
    [~,ord] = sort(P1{q,1},'descend');
    E{q,1} = E{q,1}(ord,:);  P1{q,1} = P1{q,1}(ord);
end

allv = cell2mat(cellfun(@(m) m(:), E(:), 'uni',0));
cl   = prctile(allv(isfinite(allv)), climPct);
fprintf('ChAT significant, POOLED: %d cells (%d intact, %d vagotomised)\n', ...
        NCELL, NCELL-NVAG, NVAG);
fprintf('   %s\n', strjoin(CELLS, ', '));
for q = 1:2
    fprintf('   %-6s %4d cycles  median earliest phase %.3f\n', ...
            TRIG{q}, size(E{q,1},1), median(P1{q,1}));
end

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[1 1 23 15], ...
             'Name','ChAT significant: pooled per-cycle dF/F');
set(fig,'DefaultAxesFontSize',9);
x0 = [0.105 0.565];  hHm = 0.615;  hStrip = 0.055;

for q = 1:2
    ax = axes(fig,'Position',[x0(q), 0.105, 0.375, hHm]); %#ok<LAXES>
    M = E{q,1};
    imagesc(ax, ph, 1:size(M,1), M);
    set(ax,'YDir','normal','CLim',cl,'TickDir','out','Layer','top');
    colormap(ax, gray(256));
    hold(ax,'on');
    for cb0 = 1:nCycShow-1
        plot(ax,[cb0 cb0],[0.5 size(M,1)+0.5],'-','Color',COL,'LineWidth',0.9);
    end
    hold(ax,'off');
    xlim(ax,[0 nCycShow]); ylim(ax,[0.5 size(M,1)+0.5]);
    set(ax,'XTick',0:0.5:nCycShow);
    xlabel(ax,'phase (cycles from the landmark)');
    if q == 1, ylabel(ax,'cycle (sorted by phase of earliest event)'); end

    axb = axes(fig,'Position',[x0(q), 0.105+hHm+0.005, 0.375, hStrip]); %#ok<LAXES>
    hold(axb,'on');
    b = mean(B{q,1},1,'omitnan');  b = (b-min(b))/max(max(b)-min(b),eps);
    plot(axb, ph, b, '-', 'Color',[0.35 0.35 0.35], 'LineWidth',1.1);
    for cb0 = 1:nCycShow-1
        plot(axb,[cb0 cb0],[-0.05 1.05],'-','Color',COL,'LineWidth',0.9);
    end
    hold(axb,'off');
    xlim(axb,[0 nCycShow]); ylim(axb,[-0.05 1.05]);
    set(axb,'XTick',[],'YTick',[],'Box','on');
    if q == 1, ylabel(axb,'breath','FontSize',7); end
    title(axb, sprintf('%s cycles  |  %d cycles', upper(TRIG{q}), size(M,1)), ...
          'FontSize',9,'FontWeight','normal','Interpreter','none');
end

annotation(fig,'textbox',[0.02 0.955 0.96 0.035],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold', ...
    'FontSize',10,'Interpreter','tex', ...
    'String',sprintf('ChAT SIGNIFICANT cells, POOLED  |  %d cells (%d intact, %d vagotomised)', ...
                     NCELL, NCELL-NVAG, NVAG));
annotation(fig,'textbox',[0.02 0.928 0.96 0.028],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8.5, ...
    'Interpreter','tex','String',sprintf(['%d cycles per row, each normalised by its own ' ...
    'duration  |  heatmap rows sorted by the phase of the earliest event  |  red = cycle ' ...
    'boundary'], nCycShow));
annotation(fig,'textbox',[0.02 0.901 0.96 0.028],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8.5, ...
    'Interpreter','tex','String',sprintf(['raw dF/F grayscale, clim %.3g to %.3g  |  intact ' ...
    'and vagotomised cycles are interleaved by their sort key, not separated'], cl(1), cl(2)));

stem = fullfile(outDir,'ChAT_sig_cycle_dff_heatmap');
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','image','BackgroundColor','white');
fprintf('\nwrote %s.png / .pdf\n', stem);
