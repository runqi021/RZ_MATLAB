% polar_selected_260816.m
% -----------------------------------------------------------------------
%  polar_coh_vs_rayleigh_260808 redone for the SELECTED cells, with the dot on
%  the circular MEDIAN.
%
%      r      = Rayleigh log Z          (occupancy-weighted, Kish n_eff)
%      angle  = weighted circular MEDIAN of the event phases
%      arc    = weighted circular IQR   (25th-75th, offsets about that median)
%      colour = class from cell_classes_260816
%
%  WHY THE MEDIAN MATTERS HERE. The archive figure moved the dot to the MEAN on
%  2026-08-11 but left the arc as the median-centred IQR, and its own comment
%  records that mean and median disagree by >30 deg for 61 of 445 cells -- so the
%  arc was being translated onto a dot it was not computed about. Putting the dot
%  back on the median makes dot and arc describe the same centre again.
%
%  EVERYTHING IS RECOMPUTED FROM RAW EVENT PHASES, not from the 36-bin profile: a
%  median taken off 10 deg bins is quantised to 10 deg, which is coarse against
%  the between-class differences being read off this figure. The four statistical
%  helpers at the bottom are copied verbatim from
%  Ventral_surface_polar_coh_vs_rayleigh_260808.m so the numbers stay comparable
%  with the archive, and the script CHECKS itself against the archive CSV for
%  every cell that exists in both.
%
%  OCCUPANCY WEIGHTING IS MANDATORY, NOT OPTIONAL. The breath cycle is ~7x
%  unevenly sampled in time (most frames fall in the slow expiratory half), so an
%  unweighted median slides into expiration regardless of when the cell fires.
%  Each event carries w = 1/frames-in-its-phase-bin from ITS OWN recording.
%
%  Runqi Zhang / 2026-08-16
% -----------------------------------------------------------------------
clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));
addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
addpath(fullfile(repoRoot,'falloff-analysis-260805'));
%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
regFile   = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
archiveCsv= fullfile(sumRoot,'polar_coh_vs_rayleigh_260808','polar_scores_percell.csv');
outDir    = fullfile(sumRoot,'polar_selected_260816');
sertMode  = 'included';     % see cell_classes_260816
alpha1    = 0.05;           % inner dashed circle
alpha2    = 0.001;          % outer dashed circle
nShuffle  = 0;              % per-cell Rayleigh shuffle; 0 = analytic p only
% =========================================================
if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
logZ_crit  = log(-log(alpha1));     % 1.10
logZ_crit2 = log(-log(alpha2));     % 1.93
%% ---- the selected cells ----
% The classified cells, PLUS the IO cells kept in the popsel GUI. IO was curated
% but not classified, so it carries the class label 'IO (included)'.
[CL, INFO] = cell_classes_260816(sertMode);
ioDec = fullfile(sumRoot,'popsel_260816','popsel_decisions_IO.csv');
if isfile(ioDec)
    Ti  = readtable(ioDec,'TextType','string');
    ioC = Ti.cell(Ti.decision == "include");
    ioC = setdiff(ioC, CL.cell);                 % never double-count
    if ~isempty(ioC)
        add = table(ioC(:), repmat("IO",numel(ioC),1), strings(numel(ioC),1), ...
                    repmat("IO (included)",numel(ioC),1), ...
                    repmat("popsel include",numel(ioC),1), ...
                    'VariableNames',CL.Properties.VariableNames);
        CL = [CL; add];
    end
    fprintf('added %d IO cells from the popsel curation\n', numel(ioC));
end
fprintf('%d selected cells total (Sert rule: %s)\n', height(CL), INFO.sertMode);
% genotype colours, taken from Ventral_surface_polar_coh_vs_rayleigh_260808 so
% this figure reads the same way as the archive one
GRPS = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};
GCOL = [0    0    0   ;     % IO      black
        0.85 0.10 0.10;     % ChAT    red
        0.10 0.65 0.20;     % Vglut2  green
        0.10 0.30 0.85;     % Vgat    blue
        0.55 0.20 0.75;     % Sst     purple
        0.90 0.45 0.10];    % Sert    orange
D = load(regFile,'CELL','OBS','REC');
obsOf = pooled_obs_260814(D.CELL, D.OBS);
P = struct('doCoh',false,'nDrop',30,'fallback_fps',30,'TW_spec',6,'alpha_sig',0.01, ...
    'minSpikes',2,'ca_lag_sec',0,'f_breath_search',[0.2 4],'fwhm_factor',0.6, ...
    'min_bw',0.05,'fmin',0.05,'fmax',15,'trigWinIBI',2, ...
    'ylim_dff',[],'ylim_epc',[],'histBinFrames',2,'nShuffle',0,'shiftMinCyc',3, ...
    ... % FIXED epoch window in seconds. Without this the epochs are cut at
    ... % +/-1 IBI of each cell, so on a seconds axis a 0.47 s cell would vanish
    ... % from the average past 0.47 s and the mean would be over a different set
    ... % of cells at every lag.
    'trigWin_sec', 2*1.5, ...
    'pad_um',20,'clip_pct',[0.5 99.9],'scalebar_um',50,'rayPhaseBins',36,'featFs',30, ...
    'gamma_val',1,'PixelSizeBase',1.7778,'outlineLW',0.8,'sortMode','none', ...
    'dffColor',[0.2 0.7 0.2],'onsetCol',[0.9 0.1 0.1],'peakCol',[0.35 0.75 1], ...
    'statsOnly',true);
% SECONDS, not IBI. Requested 2026-08-16: the group overlays are plotted on an
% absolute time axis.
%
% This needs a FIXED epoch window, not the default one. By default the figure
% function cuts +/-trigWinIBI/2 of each cell's OWN IBI, so a 0.47 s cell would
% only ever supply +/-0.47 s and would drop out of the average past that lag,
% making the number of contributing cells a function of time. Forcing
% P.trigWin_sec below gives every cell the same +/-WIN_S seconds regardless of
% how fast it breathes, so the mean is over the same cells at every lag.
%
% CAVEAT that cannot be fixed the same way: the event histogram window is tied to
% IBI inside the figure function (winSecH = IBI) and is shared with the
% permutation test, so changing it would move every p-value in the project. The
% histogram is therefore resampled from its native +/-IBI onto the seconds axis,
% and coverage genuinely thins past the fastest cells' IBI -- the panel reports
% how many cells contribute.
WIN_S = 1.5;                      % half-width of the plotted window, seconds
tauN  = linspace(-WIN_S, WIN_S, 201);
edgeN = linspace(-WIN_S, WIN_S, 42);  ctrN = (edgeN(1:end-1)+edgeN(2:end))/2;
n = height(CL);
R = struct('cell',num2cell(CL.cell),'key',[],'group',[],'class',[], ...
           'th',[],'Rbar',[],'Z',[],'logZ',[],'nEv',[],'nEff',[], ...
           'thMed',[],'q25',[],'q75',[],'p_an',[],'dffN',[],'histN',[], ...
           'dOnPk',[]);
tic;
for k = 1:n
    c = CL.cell(k);
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = obsOf{c}(:)'
        p = regexp(D.OBS(o).label,'/','split');
        fp = D.REC(D.OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp,'roi',str2double(p{end}), ...
            'recName',strjoin(p(3:end-1),'/'),'group',p{1},'recDate',p{2}); %#ok<SAGROW>
    end
    [~,~,st] = temporal_phase_cell_fig_260812(O, P);
    a = st.evPhase(:);  w = st.evWeight(:);
    R(k).key   = D.CELL(c).key;
    R(k).group = CL.group(k);
    R(k).class = CL.class(k);
    R(k).nEv   = numel(a);
    % ONSET-triggered curves on the SECONDS axis (tau is already in seconds)
    R(k).dffNo  = interp1(st.tau,      st.muOnset,   tauN, 'linear', NaN);
    R(k).dffNoz = interp1(st.tau,      st.muOnsetZ,  tauN, 'linear', NaN);
    R(k).histNo = interp1(st.histCtrs, st.histOnset, ctrN, 'linear', NaN);
    % Median inspiration onset->peak interval for this cell's recordings, in
    % seconds. On an onset-triggered axis this is where the breath PEAK sits.
    R(k).dOnPk  = st.dOnPk;
    if numel(a) >= 2
        [th, Rbar, nEff] = wresultant(a, w);
        Z = nEff * Rbar^2;
        med = circ_median_w(a, w);
        u   = angle(exp(1i*(a - med)));         % residuals about the MEDIAN
        R(k).th   = th;    R(k).Rbar = Rbar;   R(k).Z = Z;
        R(k).logZ = log(max(Z,eps));           R(k).nEff = nEff;
        R(k).thMed = med;
        R(k).q25 = wprctile(u, w, 25);
        R(k).q75 = wprctile(u, w, 75);
        R(k).p_an = exp(-Z);                    % p ~ exp(-Z), n-independent
    else
        [R(k).th,R(k).Rbar,R(k).Z,R(k).logZ,R(k).nEff] = deal(NaN);
        [R(k).thMed,R(k).q25,R(k).q75,R(k).p_an] = deal(NaN);
    end
    if mod(k,10)==0, fprintf('  %d/%d (%.0f s)\n', k, n, toc); end
end
fprintf('computed %d cells in %.1f min\n\n', n, toc/60);
%% ---- SELF-CHECK ----
% PRIMARY reference is pop_features.mat, because that is the CURRENT pipeline --
% the same temporal_phase_cell_fig_260812 that produced every per-cell figure,
% every star and every logZ used elsewhere today. If this script disagrees with
% it, this script is wrong.
fprintf('=== check 1: against the CURRENT pipeline (pop_features) ===\n');
PF = load(fullfile(sumRoot,'pop_analysis_260813','pop_features.mat'),'T');
dz = []; miss = 0;
for k = 1:n
    j = find(PF.T.cell == R(k).cell, 1);
    if isempty(j), miss = miss + 1; continue; end
    dz(end+1) = R(k).logZ - PF.T.logZ(j); %#ok<SAGROW>
end
fprintf('  compared %d cells, max |dlogZ| = %.3e  (%d not in pop_features)\n', ...
        numel(dz), max(abs(dz)), miss);
assert(isempty(dz) || max(abs(dz)) < 1e-9, ...
       'recomputed Rayleigh disagrees with the current pipeline');
fprintf('  EXACT match.\n');
% SECONDARY: the 2026-08-11 archive. It is a DIFFERENT VINTAGE and is expected to
% differ -- reported, not asserted. Same cells and the same event counts, but the
% phases differ; sweeping P.rayLagSec shows the archive sits nearest a +0.05 s
% calcium-lead correction (mean |dlogZ| 0.15 vs 0.52 at zero lag) without any lag
% reproducing it exactly, so more than one thing changed. The current pipeline
% uses NO lag, which is the deliberate convention here.
fprintf('\n=== check 2: against the 2026-08-11 archive (informational) ===\n');
if isfile(archiveCsv)
    A = readtable(archiveCsv,'TextType','string');
    dZ = []; dM = []; nIn = 0;
    for k = 1:n
        j = find(A.cell_key == string(R(k).key), 1);
        if isempty(j), continue; end
        nIn = nIn + 1;
        dZ(end+1) = R(k).logZ - A.ray_logZ(j); %#ok<SAGROW>
        dM(end+1) = rad2deg(abs(angle(exp(1i*(deg2rad(A.ray_med_deg(j)) - R(k).thMed))))); %#ok<SAGROW>
    end
    if nIn > 0
        fprintf('  %d of %d cells exist in the archive\n', nIn, n);
        fprintf('  |dlogZ|   median %.3f  max %.3f\n', median(abs(dZ)), max(abs(dZ)));
        fprintf('  |dmedian| median %.1f deg  max %.1f deg\n', median(dM), max(dM));
        fprintf('  event counts agree, so this is a phase-convention difference,\n');
        fprintf('  not a different event set. Current pipeline = no lag.\n');
    end
else
    fprintf('  archive CSV not found\n');
end
%% ---- figure ----
% colour by GENOTYPE, as in the archive figure
gname = string([R.group]);
present = GRPS(ismember(GRPS, cellstr(unique(gname))));
colOf = @(g) GCOL(find(strcmp(GRPS, g),1), :);
% Layout copied from Ventral_surface_polar_coh_vs_rayleigh_260808 so the two
% figures can be read side by side: square-ish canvas, polar axes left, legend
% right, two-line header, footer banner naming the lead correction.
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 20], ...
             'Name','polar: selected cells, median phase');
set(fig,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);
ax = polaraxes(fig,'Position',[0.05 0.10 0.66 0.72]);
hold(ax,'on');
rMax = max([R.logZ]);  if ~isfinite(rMax) || rMax <= 0, rMax = 1; end
rLim = ceil(rMax*1.15);
% significance circles, labelled inline as in the archive
tt = linspace(0,2*pi,361);
polarplot(ax, tt, logZ_crit *ones(size(tt)), '--','Color',[0.55 0.55 0.55],'LineWidth',0.9);
polarplot(ax, tt, logZ_crit2*ones(size(tt)), '--','Color',[0.15 0.15 0.15],'LineWidth',0.9);
text(ax, deg2rad(-18), logZ_crit,  '  \alpha=0.05',  'FontSize',7,'Color',[0.35 0.35 0.35]);
text(ax, deg2rad(-30), logZ_crit2, '  \alpha=0.001', 'FontSize',7,'Color',[0.15 0.15 0.15]);
% IQR arcs first, so dots sit on top. Arcs only on cells outside the alpha=0.05
% circle: below it the direction is not distinguishable from uniform, so drawing
% a spread there implies a precision the data does not have.
for k = 1:n
    if ~isfinite(R(k).thMed) || ~isfinite(R(k).logZ), continue; end
    if R(k).logZ < logZ_crit, continue; end
    arc = linspace(R(k).thMed + R(k).q25, R(k).thMed + R(k).q75, 40);
    polarplot(ax, arc, max(R(k).logZ,0)*ones(size(arc)), '-', ...
              'Color',[colOf(R(k).group) 0.5],'LineWidth',1.4);
end
h = gobjects(numel(present),1);
for gi = 1:numel(present)
    m = find(gname == string(present{gi}) & isfinite([R.thMed]));
    if isempty(m), continue; end
    % IO is drawn OPEN, as in the archive: it is by far the largest group and
    % filled markers would bury every other genotype under a black cloud.
    if strcmp(present{gi},'IO')
        fc = 'none';  ec = [0.20 0.20 0.20];  lw = 0.6;
    else
        fc = colOf(present{gi});  ec = 'w';   lw = 0.5;
    end
    h(gi) = polarplot(ax, [R(m).thMed], max([R(m).logZ],0), 'o', ...
        'MarkerFaceColor',fc,'MarkerEdgeColor',ec,'LineWidth',lw, ...
        'MarkerSize',6.5,'LineStyle','none');
end
rlim(ax,[0 rLim]);
ax.ThetaTick = 0:30:330;
ax.ThetaTickLabel = arrayfun(@(t) sprintf('%d',t), 0:30:330, 'uni',0);
ax.RAxisLocation = 180;
ax.GridAlpha = 0.15;
title(ax, {'r = log Z      dot = circular MEDIAN,  arc = IQR', ...
           sprintf('(arcs on cells outside the \\alpha = %.2g circle)', alpha1)}, ...
      'Interpreter','tex','FontWeight','bold');
% legend text matches the archive: n and how many clear logZ >= 2
ok = isgraphics(h);
lg = present(ok);
for q = 1:numel(lg)
    mm = gname == string(lg{q});
    lg{q} = sprintf('%s (n=%d, logZ>=2: %d)', lg{q}, nnz(mm), nnz([R(mm).logZ] >= 2));
end
lgd = legend(ax, h(ok), lg, 'Box','off','FontSize',9);
lgd.Units = 'normalized';  lgd.Position = [0.735 0.52 0.24 0.22];
% two-line header, and a footer banner naming the lead correction -- the archive
% carries one saying "events shifted 30 ms EARLIER", so this states the opposite
% explicitly rather than leaving it to be inferred.
nObsTot = sum(arrayfun(@(x) numel(obsOf{x}), [R.cell]));
annotation(fig,'textbox',[0.02 0.90 0.96 0.09],'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','tex','FontSize',10, ...
    'String',{sprintf('N = %d cells  (%d ROI-observations pooled)   |   insp onset = 0, insp peak = \\pi', ...
                      n, nObsTot), ...
              sprintf('logZ crit(\\alpha=%.2g) = %.2f \\rightarrow %d sig      logZ crit(\\alpha=%.3g) = %.2f \\rightarrow %d sig', ...
                      alpha1, logZ_crit, nnz([R.logZ]>=logZ_crit), ...
                      alpha2, logZ_crit2, nnz([R.logZ]>=logZ_crit2))});
annotation(fig,'textbox',[0.02 0.005 0.96 0.045],'EdgeColor',[0.7 0.7 0.7], ...
    'BackgroundColor',[0.96 0.96 0.96],'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontWeight','bold','FontSize',9, ...
    'String','NO GCaMP lead correction: events used as detected (archive figure shifts them 30 ms earlier)');
% ---- side panel: per-genotype, then per-class ----
axT = axes(fig,'Position',[0.735 0.10 0.25 0.38]); axis(axT,'off');
L = strings(0,1);
L(end+1) = sprintf('%-14s %3s %7s %7s %6s', 'genotype','n','medDeg','logZ','n sig');
L(end+1) = string(repmat('-',1,42));
for gi = 1:numel(present)
    m = find(gname == string(present{gi}) & isfinite([R.thMed]));
    if isempty(m), continue; end
    % population direction = circular mean of the per-cell MEDIANS, unweighted:
    % each cell counts once whatever its event count
    mu = angle(mean(exp(1i*[R(m).thMed])));
    L(end+1) = sprintf('%-14s %3d %7.0f %7.2f %6d', present{gi}, numel(m), ...
                       mod(rad2deg(mu),360), median([R(m).logZ]), ...
                       nnz([R(m).logZ] >= logZ_crit));
end
L(end+1) = "";
L(end+1) = sprintf('%-24s %3s %7s', 'class','n','medDeg');
L(end+1) = string(repmat('-',1,42));
for cc = unique([R.class],'stable')
    m = find([R.class] == cc & isfinite([R.thMed]));
    if isempty(m), continue; end
    mu = angle(mean(exp(1i*[R(m).thMed])));
    L(end+1) = sprintf('%-24s %3d %7.0f', cc, numel(m), mod(rad2deg(mu),360));
end
L(end+1) = "";
L(end+1) = sprintf('logZ >= %.2f (a=%.3g): %d of %d', logZ_crit, alpha1, ...
                   nnz([R.logZ] >= logZ_crit), n);
L(end+1) = sprintf('logZ >= %.2f (a=%.3g): %d of %d', logZ_crit2, alpha2, ...
                   nnz([R.logZ] >= logZ_crit2), n);
L(end+1) = "";
L(end+1) = "occupancy-weighted throughout;";
L(end+1) = "Kish n_eff, not raw n.";
text(axT, 0, 1, strjoin(L, newline), 'VerticalAlignment','top', ...
     'FontName','Consolas','FontSize',8,'Interpreter','none');
exportgraphics(fig, fullfile(outDir,'polar_selected.png'),'Resolution',600,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'polar_selected.pdf'),'ContentType','vector','BackgroundColor','white');
%% ---- group overlays: mean dF/F and mean event histogram ----
% Every genotype on one pair of axes, so the phase relationships are read against
% each other rather than across figures. Mean +/- SEM ACROSS CELLS: each cell
% counts once regardless of how many events it fired, otherwise one 315-event
% Vglut2 cell would outweigh twenty 8-event IO cells.
% The overlay splits genotypes that contain more than one class -- a single
% "Vglut2" line would average pre-I and I cells whose peaks sit ~60 deg apart and
% show neither. The polar panel above stays coloured by genotype.
dispG = strings(1,n);
for k = 1:n
    if R(k).group == "IO", dispG(k) = "IO";
    else,                  dispG(k) = R(k).group + " " + R(k).class; end
end
% Colours come from ONE definition shared with the cartoon map, so a green trace
% here and a green dot there are the same population. Edit them there, not here.
DISP = group_display_colors_260816();
% Genotypes dropped from the OVERLAY ONLY (2026-08-17). The polar panel above and
% the selection itself are untouched -- this removes lines from these three axes
% and nothing else. IO is 107 of ~156 cells and its flat black trace dominates a
% figure that is about the tuned populations; Sst is entirely 'null', so it adds
% a flat line with no tuning to read.
OVL_EXCLUDE_GENO = ["IO","Sst"];
dispGeno = string(DISP(:,1));
for xg = OVL_EXCLUDE_GENO
    dispGeno(startsWith(dispGeno, xg)) = "";      % blanked -> filtered out below
end
keep = ismember(string(DISP(:,1)), unique(dispG)) & dispGeno ~= "";
DISP = DISP(keep,:);
% The cells actually drawn. Everything else on this figure -- the breath-peak
% line, the contributing-cell count and the header n -- is computed from THESE
% cells, not from all of R, or the figure would describe a population it does
% not show.
inOvl = ismember(dispG, string(DISP(:,1)));
fprintf('overlay: %d of %d cells (excluded genotypes: %s)\n', ...
        nnz(inOvl), numel(dispG), strjoin(cellstr(OVL_EXCLUDE_GENO), ', '));
lg2 = cell(size(DISP,1),1);
for q = 1:size(DISP,1)
    lg2{q} = sprintf('%s (n=%d)', DISP{q,1}, nnz(dispG==string(DISP{q,1})));
end
% One 3-panel figure per trigger. Raw dF/F and z-scored dF/F are BOTH shown
% because they answer different questions: raw is the signal actually recorded
% but is dominated by the brightest cells, z-scored gives every cell equal weight
% and so shows the SHAPE of the modulation. Presenting only one of them would
% hide either the amplitude or the timing.
TRIGS = { 'onset', 'dffNo', 'dffNoz', 'histNo', 'insp onset', 'group_overlays_onset' };
for ti = 1:size(TRIGS,1)
    figT = figure('Color','w','Units','centimeters','Position',[2 2 36 12], ...
                  'Name',['group overlays - ' TRIGS{ti,1}]);
    set(figT,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);
    tlT = tiledlayout(figT,1,3,'TileSpacing','compact','Padding','compact');
    ax1t = nexttile(tlT); hold(ax1t,'on');
    ax2t = nexttile(tlT); hold(ax2t,'on');
    ax3t = nexttile(tlT); hold(ax3t,'on');
    hT = gobjects(size(DISP,1),1);
    for gi = 1:size(DISP,1)
        m = find(dispG == string(DISP{gi,1}));
        if isempty(m), continue; end
        col = DISP{gi,2};
        trio = { TRIGS{ti,2}, ax1t, tauN
                 TRIGS{ti,3}, ax2t, tauN
                 TRIGS{ti,4}, ax3t, ctrN };
        for q = 1:3
            M  = cell2mat({R(m).(trio{q,1})}.');
            mu = mean(M,1,'omitnan');
            se = std(M,0,1,'omitnan') ./ max(sqrt(sum(isfinite(M),1)),1);
            xx = trio{q,3};
            fill(trio{q,2},[xx fliplr(xx)],[mu+se fliplr(mu-se)], col, ...
                 'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
            hh = plot(trio{q,2}, xx, mu, '-','Color',col,'LineWidth',1.8);
            if q == 1, hT(gi) = hh; end
        end
    end
    % Breath landmarks. RED = inspiration onset, which is the trigger and so sits
    % at 0 by construction. BLUE = the breath PEAK, at the population average of
    % each cell's median onset->peak interval. The peak time is a property of the
    % breathing, not of the cells: it is averaged over the cells plotted here only
    % because that is the population the panels describe.
    dPk = [R(inOvl).dOnPk];  dPk = dPk(isfinite(dPk));
    tPk = mean(dPk);
    for axq = [ax1t ax2t ax3t]
        xline(axq, 0, '-', 'Color',[0.90 0.10 0.10],'LineWidth',1.2);
        if ~isempty(dPk) && tPk <= WIN_S
            xline(axq, tPk, '-', 'Color',[0.35 0.75 1.00],'LineWidth',1.2);
        end
        xlim(axq,[-WIN_S WIN_S]); grid(axq,'on'); box(axq,'on');
        xlabel(axq, sprintf('time from %s (s)', TRIGS{ti,5}));
    end
    % How many cells actually contribute at each lag. The dF/F panels are complete
    % because the epoch window was forced to a fixed number of SECONDS; the
    % histogram is not, because its window is tied to IBI inside the figure
    % function and shared with the permutation test. Stated rather than hidden.
    Hall = cell2mat({R(inOvl).histNo}.');
    nEdge = min(sum(isfinite(Hall),1));
    text(ax3t, 0.02, 0.02, sprintf('cells contributing: %d at 0 s, %d at the edges', ...
         max(sum(isfinite(Hall),1)), nEdge), 'Units','normalized', ...
         'FontSize',7.5,'Interpreter','none');
    ylabel(ax1t,'\DeltaF/F');
    title(ax1t,'mean \DeltaF/F \pm SEM','FontWeight','normal');
    ylabel(ax2t,'z-scored \DeltaF/F');
    title(ax2t,'mean z-scored \DeltaF/F \pm SEM','FontWeight','normal');
    ylabel(ax3t,'spk/cyc %');
    title(ax3t,'event histogram','FontWeight','normal');
    okT = isgraphics(hT);
    legend(ax1t, hT(okT), lg2(okT), 'Location','northwest','Box','off','FontSize',7.5);
    title(tlT, sprintf(['inspiration-%s-triggered, absolute time  ' ...
          '(red = onset, blue = mean breath peak at %.0f ms, n = %d cells)'], ...
          TRIGS{ti,1}, 1000*tPk, numel(dPk)), 'FontWeight','normal');
    exportgraphics(figT, fullfile(outDir,[TRIGS{ti,6} '.png']), ...
                   'Resolution',600,'BackgroundColor','white');
    exportgraphics(figT, fullfile(outDir,[TRIGS{ti,6} '.pdf']), ...
                   'ContentType','vector','BackgroundColor','white');
end

%% ---- CSV ----
T = table([R.cell].', string({R.key}).', string([R.group]).', string([R.class]).', ...
          [R.nEv].', [R.nEff].', [R.Rbar].', [R.Z].', [R.logZ].', ...
          mod(rad2deg([R.th]).',360), mod(rad2deg([R.thMed]).',360), ...
          rad2deg([R.q25]).', rad2deg([R.q75]).', rad2deg([R.q75]-[R.q25]).', [R.p_an].', ...
          'VariableNames',{'cell','cell_key','group','class','n_events','n_eff', ...
          'Rbar','Z','logZ','mean_deg','median_deg','iqr_lo_off_deg','iqr_hi_off_deg', ...
          'iqr_width_deg','p_analytic'});
T = sortrows(T,{'class','group','cell'});
writetable(T, fullfile(outDir,'polar_selected_scores.csv'));
save(fullfile(outDir,'polar_selected_260816.mat'),'R','CL','INFO','logZ_crit','logZ_crit2');
fprintf('\nsaved -> %s\n', outDir);
disp(T(:,{'cell','group','class','n_events','logZ','median_deg','iqr_width_deg'}));

%% ================= helpers, copied verbatim from =====================
%% Ventral_surface_polar_coh_vs_rayleigh_260808.m so the maths is identical
function m = circ_median_w(a, w)
% Weighted circular median: the direction minimising the total weighted angular
% distance to the data. Solved on a 0.5 deg grid -- a closed form does not exist
% on the circle, and the grid is exact enough next to a 10 deg histogram.
    g = linspace(-pi, pi, 721);  g(end) = [];
    d = abs(angle(exp(1i*(a(:) - g))));      % nEv x nGrid angular distances
    [~, k] = min(sum(w(:) .* d, 1));
    m = g(k);
end

function q = wprctile(v, w, p)
% Weighted percentile, midpoint rule on the cumulative weight.
    [v, o] = sort(v(:));  w = w(o);  w = w / sum(w);
    c = cumsum(w) - 0.5*w;
    if numel(v) < 2, q = v(1); return; end
    [c, iu] = unique(c);  v = v(iu);
    q = interp1(c, v, p/100, 'linear', 'extrap');
end

function [th, Rbar, nEff] = wresultant(a, w)
% Weighted circular resultant with Kish's effective sample size. Using raw n with
% a weighted Rbar breaks the Rayleigh: a few events in a sparsely-sampled
% inspiratory bin carry huge leverage.
    S1 = sum(w);  S2 = sum(w.^2);
    if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
    nEff = S1^2 / max(S2, eps);
    v    = sum(w(:) .* exp(1i*a(:))) / S1;
    th   = angle(v);
    Rbar = min(abs(v), 1);
end
