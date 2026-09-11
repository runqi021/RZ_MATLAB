% polar_cloud_wn_260811.m
% -----------------------------------------------------------------------
%  Polar plot of the PHASE-LOCKED cells (logZ >= zMin) drawn as WRAPPED NORMAL
%  clouds instead of dots.
%
%  Each cell becomes a petal centred on its preferred phase:
%
%      angle    mu    = RY.th          the weighted mean direction
%      spread   sigma = sqrt(-2 ln Rbar)
%      radius   logZ                   same radial axis as the polar figure
%      AREA     proportional to EVENTS PER BREATH CYCLE
%
%  WHY THE WRAPPED NORMAL IS THE RIGHT OBJECT, AND WHY sigma NEEDS NO NEW FIT.
%  For a wrapped normal the first trigonometric moment is Rbar = exp(-sigma^2/2),
%  so sigma = sqrt(-2 ln Rbar) EXACTLY. The Rayleigh statistic already gives Rbar
%  (occupancy-weighted), so the WN sigma falls out of what is already computed --
%  no extra estimation, and the spread inherits the same weighting as the test.
%  The density is evaluated by summing the k = -2..2 wraps, which converges long
%  before the Jacobi theta form is needed at these sigmas.
%
%  sigma DIVERGES as Rbar -> 0. An unlocked cell has no meaningful width, which is
%  precisely why this figure is restricted to logZ >= zMin.
%
%  AREA = EVENTS PER CYCLE. The petal's area is set to (events / breath cycles)
%  times areaScale, so a cell that fires once every cycle is one unit of ink and a
%  cell that fires every tenth cycle is a tenth of it. Breath cycles are counted
%  per recording from the curated onsets, summed over every recording the cell
%  appears in -- the same recordings its events were pooled from.
%
%  Input : polar_coh_vs_rayleigh_260808\polar_coh_vs_rayleigh_data.mat
%          the archive (breath_insp_start_pc1.mat, for cycle counts)
%  Output: polar_cloud_wn.png / .pdf / .csv
%
%  Runqi Zhang / 2026-08-11

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

%% ===================== USER-EDITABLE =====================
rootPath  = 'D:\Ventral_surface_summary';
rayFile   = fullfile(rootPath,'polar_coh_vs_rayleigh_260808','polar_coh_vs_rayleigh_data.mat');
outDir    = fullfile(rootPath,'breath_trig_heatmap_260806');
zMin      = 2;             % keep cells with logZ >= this
areaScale = 0.55;          % petal area per (event / cycle), in polar-axis units
nDrop     = 30;            % breath frames tossed to align with calcium
spanSD    = 1.5;           %% petal spans mu +/- spanSD*sigma. Cutting at a fixed
                           %% FRACTION of the peak instead (e.g. 2%%) reaches 2.8 sigma,
                           %% which at the median sigma of 54 deg wraps a petal most of
                           %% the way round the circle and makes the figure unreadable.
maxSpanDeg = 120;          %% hard cap, so a near-uniform cell cannot encircle the plot
alphaFace = 0.45;
doSave    = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
R = load(rayFile,'RY','cellLabels','labels','obsOfCell','groups','group_colors','logZ_crit','logZ_crit2');
RY = R.RY;

%% ---- index every recording folder in the archive, once ----
recDir = containers.Map('KeyType','char','ValueType','char');
hits = dir(fullfile(rootPath,'**','breath_insp_start_pc1.mat'));
for h = 1:numel(hits)
    pp = strsplit(hits(h).folder, filesep);
    if ~isKey(recDir, pp{end}), recDir(pp{end}) = hits(h).folder; end
end
cycCache = containers.Map('KeyType','char','ValueType','double');

%% ---- per cell: which recordings, how many breath cycles ----
nCell = numel(RY.logZ);
nCyc  = zeros(nCell,1);
for c = 1:nCell
    obs = find(R.obsOfCell == c);
    for o = obs(:).'
        parts = regexp(R.labels{o}, '/', 'split');
        if numel(parts) < 4, continue; end
        rec = parts{3};
        if ~isKey(cycCache, rec)
            n = NaN;
            if isKey(recDir, rec)
                IP = load(fullfile(recDir(rec),'breath_insp_start_pc1.mat'),'insp_start_idx');
                on = round(IP.insp_start_idx(:)) - nDrop;
                n  = max(0, numel(on(on>=1)) - 1);       % complete cycles
            end
            cycCache(rec) = n;
        end
        if ~isnan(cycCache(rec)), nCyc(c) = nCyc(c) + cycCache(rec); end
    end
end

%% ---- select and derive ----
sel = find(RY.logZ >= zMin & ~isnan(RY.th) & RY.Rbar > 0 & nCyc > 0);
mu    = RY.th(sel);
sigma = sqrt(-2*log(min(RY.Rbar(sel), 1-1e-9)));      % WN sigma, radians
rad   = RY.logZ(sel);
evPer = RY.nEv(sel) ./ nCyc(sel);
grpOf = cellfun(@(s) strtok(s,'/'), R.cellLabels(sel), 'uni', 0);

fprintf('%d cells with logZ >= %g\n', numel(sel), zMin);
fprintf('  sigma  %.0f - %.0f deg (median %.0f)\n', rad2deg(min(sigma)), rad2deg(max(sigma)), rad2deg(median(sigma)));
fprintf('  events/cycle  %.3f - %.3f (median %.3f)\n', min(evPer), max(evPer), median(evPer));
for g = 1:numel(R.groups)
    m = strcmp(grpOf, R.groups{g});
    if any(m), fprintf('  %-7s n=%2d | sigma %3.0f deg | ev/cyc %.2f\n', ...
        R.groups{g}, nnz(m), rad2deg(median(sigma(m))), median(evPer(m))); end
end

%% ---- figure ----
% Drawn on CARTESIAN axes with a hand-built polar frame: R2021b refuses to make
% a patch a child of polaraxes, and the clouds must be filled patches.
fh = figure('Color','w','Position',[60 40 950 900]);
ax = axes(fh); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
rMax = max(rad)*1.15;
tt = linspace(0,2*pi,361);

% rings: the two significance circles dashed, plain integers light
for z = 1:floor(rMax)
    plot(ax, z*cos(tt), z*sin(tt), '-', 'Color',[0.88 0.88 0.88], 'LineWidth',0.7);
    text(ax, z*cosd(64), z*sind(64), sprintf('%d',z), 'FontSize',8, 'Color',[0.5 0.5 0.5]);
end
for z = [R.logZ_crit R.logZ_crit2]
    plot(ax, z*cos(tt), z*sin(tt), '--', 'Color',[0.45 0.45 0.45], 'LineWidth',0.9);
end
plot(ax, rMax*cos(tt), rMax*sin(tt), '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1);

% spokes + phase labels
lblT = {'0 (onset)','','\pi/2','','\pi (peak)','','3\pi/2',''};
for k = 0:7
    a = k*pi/4;
    plot(ax, [0 rMax]*cos(a), [0 rMax]*sin(a), '-', 'Color',[0.9 0.9 0.9], 'LineWidth',0.6);
    if ~isempty(lblT{k+1})
        text(ax, 1.07*rMax*cos(a), 1.07*rMax*sin(a), lblT{k+1}, 'FontSize',10, ...
             'HorizontalAlignment','center');
    end
end

for i = 1:numel(sel)
    g = find(strcmp(R.groups, grpOf{i}), 1);
    if isempty(g), continue; end
    col = R.group_colors(g,:);
    hw = min(spanSD*sigma(i), deg2rad(maxSpanDeg));      %% angular half-width
    th = mu(i) + linspace(-hw, hw, 181);
    d  = zeros(size(th));
    for k = -2:2
        d = d + exp(-((th - mu(i) - 2*pi*k).^2) / (2*sigma(i)^2));
    end
    d = d / max(d);
    A = areaScale * evPer(i);
    w = A / max(rad(i) * trapz(th, d), eps);
    rIn = rad(i) - 0.5*w*d;  rOut = rad(i) + 0.5*w*d;
    thP = [th fliplr(th)];  rP = [rOut fliplr(rIn)];
    patch(ax, rP.*cos(thP), rP.*sin(thP), col, 'FaceAlpha',alphaFace, ...
          'EdgeColor',col, 'EdgeAlpha',0.7, 'LineWidth',0.5);
end

lg = gobjects(numel(R.groups),1);
for g = 1:numel(R.groups)
    lg(g) = patch(ax, NaN, NaN, R.group_colors(g,:), 'EdgeColor','none');
end
n_g = arrayfun(@(g) nnz(strcmp(grpOf, R.groups{g})), 1:numel(R.groups));
legend(ax, lg(n_g>0), arrayfun(@(g) sprintf('%s (n=%d)', R.groups{g}, n_g(g)), find(n_g>0), 'uni',0), ...
       'Location','southoutside','Orientation','horizontal','Box','off');
xlim(ax, 1.2*rMax*[-1 1]); ylim(ax, 1.2*rMax*[-1 1]);
title(ax, {sprintf('wrapped-normal clouds, %d cells with logZ \\geq %g', numel(sel), zMin), ...
           'radius = logZ,  \mu = RY.th,  width = \sigma = \surd(-2 ln Rbar),  AREA \propto events per breath cycle'}, ...
      'FontWeight','bold');

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir,'polar_cloud_wn');
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    writetable(table(grpOf, R.cellLabels(sel), rad, rad2deg(mod(mu,2*pi)), rad2deg(sigma), ...
                     RY.nEv(sel), nCyc(sel), evPer, ...
        'VariableNames',{'group','cell','logZ','mu_deg','sigma_deg','n_events','n_cycles','ev_per_cycle'}), ...
        [base '.csv']);
    fprintf('saved %s.{png,pdf,csv}\n', base);
end
