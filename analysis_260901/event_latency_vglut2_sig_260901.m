% event_latency_vglut2_sig_260901.m
% -----------------------------------------------------------------------
%  Onset-triggered event-latency population histogram for the SIGNIFICANT
%  Vglut2 cells, vagal-intact beside vagotomised. Two panels, nothing else.
%
%  Built to match the population panels of Ventral_surface_event_latency_260811
%  (the *_population.png figures): same latency definition, same epc unit, same
%  frame-snapped bin grid, same bar + outline drawing. What differs is the
%  population (the re-check curation, significant cells only), the split
%  (vagotomy), and that only the ONSET trigger is drawn.
%
%  THE LATENCY, as in the archive: every event goes to its NEAREST inspiration
%  onset with a SIGNED latency, in ms, negative = the event LEADS the trigger.
%  Events further than half a cycle from every trigger are dropped, so each event
%  is counted once and only once. Triggers are wrapped by +/- the recording
%  duration before the nearest-trigger search, so events near the ends are not
%  forced onto a distant trigger.
%
%  epc = EVENTS PER CYCLE PER BIN: the pooled count in a bin divided by the total
%  number of breath cycles the contributing cells were exposed to. A cell
%  contributes in proportion to its events, which is what a population panel is
%  for -- it is NOT a mean over per-cell densities.
%
%  BINS ARE A WHOLE ODD NUMBER OF FRAMES WIDE, centred on zero lag, via the
%  archive's lat_axis rule. A latency can only take values on the frame lattice,
%  so an arbitrary bin width puts a different number of possible values in each
%  bin and the histogram alternates high-low from the grid alone.
%
%  ONE SHARED AXIS FOR BOTH PANELS. Vagotomy changes breathing rate, so a
%  per-group window would compare two different spans of time and the panels
%  could not be read against each other. The window is half the median cycle over
%  all the selected cells; each group's own median IBI is printed in its title.
%
%  SELECTION: the cells the six-panel polar calls significant -- circular-shift
%  permutation on the PEAK-triggered PSTH, p <= 0.01. The test is peak-triggered
%  and the display here is onset-triggered; that is deliberate, the selection is
%  simply "cells that are locked to the breath at all".
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
VAG_DATES = {'0824'};    % Vglut2 vagotomised session; everything else is intact
nBins     = 25;          % display bins across one IBI, as in the archive
% =========================================================

COL_T   = [0.85 0.20 0.10];      % onset red, as in the archive
COL_POP = [0.25 0.25 0.25];

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

%% ===================== SELECT THE CELLS =====================
assert(isfile(sigCsv), 'Run polar_recheck_perm_260901.m first (missing %s)', sigCsv);
Tsig = readtable(sigCsv,'TextType','string');
nGenoTot = nnz(Tsig.group == GENOTYPE);        % before the significance cut
keep = Tsig.group == GENOTYPE & logical(Tsig.sig_perm);
Tsig = Tsig(keep,:);
isVag = false(height(Tsig),1);
for v = 1:numel(VAG_DATES)
    isVag = isVag | contains(Tsig.stem, "_" + VAG_DATES{v} + "_");
end
fprintf('%s significant cells: %d  (%d intact, %d vagotomised)\n', ...
        GENOTYPE, height(Tsig), nnz(~isVag), nnz(isVag));

CU = ca_recheck_load_curation_260901(bundleDir);

%% ===================== PER-CELL LATENCIES =====================
S = struct('stem',{},'vag',{},'L',{},'nCyc',{},'ibi',{},'fps',{},'nEv',{});
for k = 1:height(Tsig)
    f = fullfile(bundleDir,'cells',char(Tsig.stem(k)) + ".mat");
    f = fullfile(bundleDir,'cells',[char(Tsig.stem(k)) '.mat']);
    if ~isfile(f), fprintf(2,'  missing bundle file: %s\n', f); continue; end
    L = load(f,'RECc');  REC = L.RECc;

    Lall = []; nCyc = 0; dtAll = []; fpsAll = [];
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && (g.toss || g.cell_toss), continue; end
        ev = r.spike_idx(:);
        if g.found, ev = g.spike_idx(:); end
        ts = r.foot(:)/r.fps;                    % inspiration ONSET, seconds
        if numel(ts) < 2 || isempty(ev), continue; end
        Trec = r.T / r.fps;
        dtAll  = [dtAll; diff(ts)];              %#ok<AGROW>
        nCyc   = nCyc + numel(ts) - 1;
        fpsAll = [fpsAll; r.fps];                %#ok<AGROW>

        % nearest trigger, with the trigger train wrapped so the ends behave
        te  = [ts(:).'-Trec, ts(:).', ts(:).'+Trec];
        edg = [-inf, (te(1:end-1)+te(2:end))/2, inf];
        tE  = ev/r.fps;  tE = tE(tE >= 0 & tE <= Trec);
        if isempty(tE), continue; end
        L0  = tE - te(discretize(tE, edg)).';
        Lall = [Lall; 1000*L0];                  %#ok<AGROW>
    end
    if isempty(dtAll), continue; end
    ibi  = mean(dtAll);
    Lall = Lall(abs(Lall) <= 1000*ibi/2);        % |latency| <= half a cycle
    S(end+1) = struct('stem',Tsig.stem(k), 'vag',isVag(k), 'L',Lall, ...
                      'nCyc',nCyc, 'ibi',ibi, 'fps',mode(fpsAll), ...
                      'nEv',numel(Lall)); %#ok<SAGROW>
end
fprintf('latencies built for %d cells (%d events total)\n', numel(S), sum([S.nEv]));

%% ===================== POOLED HISTOGRAMS =====================
medIBI = median([S.ibi]);
medFPS = mode([S.fps]);
[edges, ctrs] = lat_axis_local(medIBI, nBins, medFPS);

GRP = { 'vagal intact', ~[S.vag]
        'vagotomised',   [S.vag] };
H = cell(2,1); nCell = zeros(2,1); nEvG = zeros(2,1); ibiG = nan(2,1);
for gi = 1:2
    m = GRP{gi,2};
    cnt = zeros(1,numel(ctrs)); cyc = 0;
    for k = find(m)
        cnt = cnt + histcounts(S(k).L, edges);
        cyc = cyc + S(k).nCyc;
    end
    H{gi}     = cnt / max(cyc,1);
    nCell(gi) = nnz(m);
    nEvG(gi)  = sum([S(m).nEv]);
    ibiG(gi)  = median([S(m).ibi]);
end

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[2 4 26 10], ...
             'Name','Vglut2 significant: onset-triggered latency');
set(fig,'DefaultAxesFontSize',9);
yTop = max([H{1}, H{2}, eps])*1.12;      % one y scale, or the panels lie

for gi = 1:2
    ax = axes(fig,'Position',[0.075+(gi-1)*0.480, 0.165, 0.385, 0.56]); %#ok<LAXES>
    hold(ax,'on'); box(ax,'on');
    bar(ax, ctrs, H{gi}, 1, 'FaceColor',COL_POP, 'EdgeColor','none');
    plot(ax, ctrs, H{gi}, '-', 'Color',COL_T, 'LineWidth',1.2);
    xline(ax, 0, 'k-', 'LineWidth',0.5);
    xlim(ax, [ctrs(1)-diff(ctrs(1:2))/2, ctrs(end)+diff(ctrs(1:2))/2]);
    ylim(ax, [0 yTop]);
    xlabel(ax,'latency from inspiration ONSET (ms)');
    if gi == 1, ylabel(ax,'epc'); end
    title(ax, sprintf('%s %s  |  %d cells  |  %d events  |  median IBI %.0f ms', ...
                      GENOTYPE, GRP{gi,1}, nCell(gi), nEvG(gi), 1000*ibiG(gi)), ...
          'Interpreter','none','FontSize',9,'FontWeight','normal');
end

% Three short lines, not two long ones: at 26 cm the old second line ran off
% both edges of the canvas.
sgtitle({sprintf('%s SIGNIFICANT cells: %d of %d   (permutation p \\leq 0.01, peak-triggered PSTH, 1200 shifts)', ...
                 GENOTYPE, numel(S), nGenoTot), ...
         'inspiration ONSET triggered, negative latency LEADS the trigger   |   y = epc (events per cycle per bin)', ...
         sprintf('window = \\pm half the pooled median cycle (%.0f ms), shared by both panels   |   no GCaMP lag compensation', ...
                 1000*medIBI/2)}, ...
        'Interpreter','tex','FontSize',8.5);

stem = fullfile(outDir, sprintf('%s_sig_onset_latency_vagotomy', GENOTYPE));
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');

Tout = table(repmat(string(GRP{1,1}),numel(ctrs),1), ctrs(:), H{1}(:), ...
    'VariableNames',{'group','latency_ms','epc'});
Tout = [Tout; table(repmat(string(GRP{2,1}),numel(ctrs),1), ctrs(:), H{2}(:), ...
    'VariableNames',{'group','latency_ms','epc'})];
writetable(Tout, [stem '.csv']);

fprintf('\n%-14s %6s %8s %10s %10s\n','group','cells','events','medIBI_ms','peak_ms');
for gi = 1:2
    [~,im] = max(H{gi});
    fprintf('%-14s %6d %8d %10.0f %10.0f\n', GRP{gi,1}, nCell(gi), nEvG(gi), ...
            1000*ibiG(gi), ctrs(im));
end
fprintf('\nwrote %s.png / .pdf / .csv\n', stem);

%% ===================== LOCAL =====================
function [edges, ctrs] = lat_axis_local(ibi, nBins, fps)
%LAT_AXIS_LOCAL  Verbatim from Ventral_surface_event_latency_260811.m.
%  Bin width snapped to an ODD whole number of frames, one bin centred on zero
%  lag, so every bin spans exactly the same number of frame-lattice points and no
%  latency lands on an edge.
if ~isfinite(ibi) || ibi <= 0, ibi = 1; end
if nargin < 3 || ~isfinite(fps) || fps <= 0, fps = 30; end
dt = 1000/fps;
k  = max(1, round(1000*ibi/nBins/dt));
if mod(k,2) == 0
    if 1000*ibi/nBins/dt >= k, k = k + 1; else, k = max(1, k - 1); end
end
nb    = max(1, floor(((1000*ibi/2) - k*dt/2)/(k*dt)));
ctrs  = (-nb:nb)*k*dt;
edges = [ctrs - k*dt/2, ctrs(end) + k*dt/2];
end
