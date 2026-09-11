% polar_event_latency_260811.m
% -----------------------------------------------------------------------
%  POLAR view of the time-domain event latencies. One dot per tested cell.
%
%      ANGLE   = latency / IBI, wrapped to [0, 1) of the breath cycle.
%                0 = the trigger itself, 1 = the next one. This is a LINEAR
%                normalisation by time -- milliseconds divided by that cell's
%                own cycle -- NOT the piecewise breath phase. Piecewise phase
%                stretched inspiration (~19% of the cycle) onto half the circle
%                and made a fixed GCaMP lag into a rate-dependent angle. Here a
%                cell at 0.25 fired a quarter of a cycle after the trigger, in
%                whatever milliseconds that is for that cell.
%
%      RADIUS  = -log10(p) from the permutation test. FURTHER OUT = STRONGER.
%                Rings at p = 0.05, 0.01 and 0.001 (dashed = the alpha in use).
%                This is why the source run needs 10000 shifts: at 1000 the
%                smallest attainable p is 1/1001, so every significant cell
%                would pile onto exactly r = 3 and the radial axis would say
%                nothing. At 10000 the floor is 1e-4 and the axis spans 0..4.
%
%      ARM     = the circular IQR, Q1 to Q3, drawn as an arc at that radius,
%                with the DOT at the circular MEDIAN. Both are anchored on the
%                peak (see cell_latency): the window is a whole cycle and WRAPS,
%                so a cell peaking near the window edge would get a linear
%                median on the opposite side of the circle.
%
%      DOT AREA = events per cycle (epc), clipped to the 0.5-99.5 percentiles
%                so one hyperactive cell cannot set the scale for everyone.
%
%      COLOUR  = genotype, same key as the rest of the archive figures.
%
%  A LONG ARM IS NOT A BUG. A cell can be strongly significant (large radius)
%  and still have an arm most of the way round: that is a sharp peak sitting on
%  a broad background, where half the events are nowhere near the peak. The
%  radius says "there IS a peak", the arm says "and this is how much of the
%  cycle the middle half of the events actually covers". They answer different
%  questions and are meant to disagree.
%
%  Input : event_latency_260811\event_latency_data.mat  (written by
%          Ventral_surface_event_latency_260811.m -- nothing is recomputed here)
%  Output: polar_event_latency_<trigger>.png / .pdf  + a combined figure
%
%  Runqi Zhang / 2026-08-11

clear; clc; close all;

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath,'event_latency_260811','event_latency_data.mat');
outDir   = fullfile(rootPath,'event_latency_260811');
radialVar = 'all';          % 'neglogp' | 'z' | 'logz' | 'all'
zFloor    = 0.5;            % log-radius anchor, see the 'logz' case below
%  'neglogp'  radius = -log10(p).  What was asked for, and it SATURATES: p can
%             never go below 1/(nPerm+1), and on this archive 100% of the
%             onset-significant and 97% of the peak-significant cells sit
%             exactly on that floor, so every one of them lands on the same
%             ring and the radial axis cannot rank them.
%  'z'        radius = (T_obs - mean(T_null)) / std(T_null), i.e. how many null
%             SDs above the shuffles the peak stands. Monotonic with p, so the
%             ordering agrees wherever p can still discriminate, but it keeps
%             separating cells long after p has bottomed out.
sizeLim  = [0.5 99.5];      % percentile clip on epc -> dot area
sizeRng  = [12 320];        % dot area range, points^2
pRings   = [0.05 0.01 0.001];
sigOnly  = false;           % true = drop cells that did not reach alpha
armAlpha = 0.55;
doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
D = load(inFile,'S','CELL','groups','prm');
S = D.S; CELL = D.CELL; groups = D.groups; prm = D.prm;

% genotype colours, same key as ventral_pie_summary / the cartoon map
cmap = containers.Map( ...
    {'ChAT','Sst','Vglut2','Vgat','Sert','IO'}, ...
    {[0.85 0.10 0.10],[0.55 0.20 0.75],[0.10 0.65 0.20], ...
     [0.10 0.30 0.85],[0.90 0.45 0.10],[0.35 0.35 0.35]});
gCol = cell2mat(cellfun(@(g) cmap(g), groups(:), 'uni',0));

TRIG = {'onset','peak'};
tested = vertcat(S.tested);  sigM = vertcat(S.sig);
gi     = [CELL.gi];
fprintf('%d cells | tested %d/%d (onset/peak) | sig %d/%d | alpha %.4g, %d shifts\n', ...
        numel(S), nnz(tested(:,1)), nnz(tested(:,2)), nnz(sigM(:,1)), nnz(sigM(:,2)), ...
        prm.alphaPerm, prm.nPerm);
pFloor = 1/(prm.nPerm+1);

%% ---- assemble the plotted quantities, per trigger ----
P = struct('th',{},'r',{},'a1',{},'a3',{},'epc',{},'g',{},'sig',{});
for q = 1:2
    % NOTE the orientations: tested is N x 2 (column), arrayfun over the 1 x N
    % struct array returns a ROW. Combining them without transposing would
    % implicitly expand to an N x N logical and quietly select the wrong cells.
    ibiQ = arrayfun(@(s) s.ibi(q), S);            % 1 x N
    keep = find(tested(:,q).' & ~isnan(ibiQ));    % 1 x N -> row of indices
    if sigOnly, keep = keep(sigM(keep,q).'); end
    ibi = arrayfun(@(k) S(k).ibi(q)*1000, keep);      % ms
    med = arrayfun(@(k) S(k).med(q),      keep);
    q1  = arrayfun(@(k) S(k).q1(q),       keep);
    q3  = arrayfun(@(k) S(k).q3(q),       keep);
    pv  = arrayfun(@(k) S(k).p(q),        keep);
    nEv = arrayfun(@(k) S(k).nEv(q),      keep);
    nCy = arrayfun(@(k) S(k).nCyc(q),     keep);
    P(q).th  = 2*pi*mod(med./ibi, 1);                 % fraction of a cycle
    P(q).a1  = 2*pi*(med - q1)./ibi;                  % arm, backwards from median
    P(q).a3  = 2*pi*(q3 - med)./ibi;                  % arm, forwards
    P(q).rp  = -log10(max(pv, pFloor));
    P(q).rz  = arrayfun(@(k) S(k).z(q), keep);
    P(q).p   = pv;
    P(q).epc = nEv ./ max(nCy,1);
    P(q).g   = gi(keep);
    P(q).sig = sigM(keep,q).';
    P(q).keep = keep;
end

nAtFloor = arrayfun(@(q) nnz(P(q).p <= pFloor*1.001 & P(q).sig), 1:2);
nSigQ    = arrayfun(@(q) nnz(P(q).sig), 1:2);
fprintf(['p SATURATION: %d/%d onset and %d/%d peak significant cells sit exactly on ' ...
         'the p floor,\n              so -log10(p) puts them all on one ring. ' ...
         'z does not saturate.\n'], nAtFloor(1), nSigQ(1), nAtFloor(2), nSigQ(2));

% ONE size scale across both panels, so a dot means the same thing in each
allE = [P(1).epc, P(2).epc];
sLim = prctile(allE, sizeLim);
fprintf('epc %.3f - %.3f (clip %.3f - %.3f) | radius max %.2f\n', ...
        min(allE), max(allE), sLim(1), sLim(2), max([P(1).r P(2).r]));
sz = @(e) sizeRng(1) + diff(sizeRng) * ...
          (min(max(e, sLim(1)), sLim(2)) - sLim(1)) / max(diff(sLim), eps);

switch lower(radialVar)
    case 'neglogp', VARS = {'neglogp'};
    case 'z',       VARS = {'z'};
    case 'logz',    VARS = {'logz'};
    otherwise,      VARS = {'neglogp','z','logz'};
end
zCrit = min([P(1).rz(P(1).sig), P(2).rz(P(2).sig)]);   % smallest z that passed alpha

for vi = 1:numel(VARS)
V = VARS{vi};
rTick = []; rTickLab = {};
switch V
    case 'neglogp'
        for q=1:2, P(q).r = P(q).rp; end
        rLab = '-log_{10} p';
        rMax = max([P(1).r, P(2).r, -log10(min(pRings))]) * 1.08;
        % rings ARE the p levels
        ringR = -log10(pRings);  ringL = arrayfun(@(x) sprintf('p=%g',x), pRings, 'uni',0);
        ringD = abs(pRings - prm.alphaPerm) < eps(pRings);
    case 'z'
        for q=1:2, P(q).r = P(q).rz; end
        rLab = 'z  (null SDs)';
        rMax = max([P(1).r, P(2).r, 3]) * 1.08;
        % a p level has no single z, so the rings show where the SIGNIFICANT
        % cells actually start: the smallest z among cells that passed alpha.
        ringR = zCrit; ringL = {sprintf('p<%g starts here (z=%.1f)', prm.alphaPerm, zCrit)};
        ringD = true;
    case 'logz'
        % LOG of the same z. z spans ~0 to 19 here, so a linear radius spends
        % most of its area on the handful of very strong cells and crushes the
        % 3-8 band where most significant cells live; the log spreads that band
        % out. Plotted as log(z/zFloor) so the axis still starts at 0, and the
        % ticks are LABELLED WITH z, not with its log -- the spacing is
        % logarithmic, the numbers you read are the real ones.
        %
        % z can be <= 0 (a cell whose peak is below the shuffle mean), and log
        % is undefined there, so anything at or below zFloor is PINNED to the
        % centre. Those are all far from significant; the count is printed.
        for q=1:2, P(q).r = log(max(P(q).rz, zFloor)/zFloor); end
        nPin = sum(arrayfun(@(q) nnz(P(q).rz <= zFloor), 1:2));
        fprintf('logz: %d cell-triggers have z <= %.2f and are pinned at the centre\n', nPin, zFloor);
        rLab = 'z  (log spacing)';
        rMax = max([P(1).r, P(2).r]) * 1.08;
        zt   = [0.5 1 2 3 5 8 12 20 30];
        zt   = zt(zt >= zFloor & log(zt/zFloor) <= rMax);
        rTick = log(zt/zFloor);  rTickLab = arrayfun(@(x) sprintf('%g',x), zt, 'uni',0);
        ringR = log(zCrit/zFloor);
        ringL = {sprintf('p<%g starts here (z=%.1f)', prm.alphaPerm, zCrit)};
        ringD = true;
end

%% ---- figure ----
fh = figure('Color','w','Position',[40 60 1500 780]);
tl = tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
for q = 1:2
    ax = polaraxes('Parent',tl); ax.Layout.Tile = q;
    hold(ax,'on');
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir          = 'clockwise';
    ax.ThetaTick         = 0:45:315;
    ax.ThetaTickLabel    = {'0','.125','.25','.375','.5','.625','.75','.875'};
    ax.RLim              = [0 rMax];
    if isempty(rTick), ax.RTick = 0:1:floor(rMax);
    else,              ax.RTick = rTick;  ax.RTickLabel = rTickLab;
    end
    ax.RAxisLocation     = 22;
    ax.GridColor         = [0.85 0.85 0.85];
    ax.FontSize          = 9;

    % significance rings
    tt = linspace(0,2*pi,361);
    for ri = 1:numel(ringR)
        rr = ringR(ri);
        if ~isfinite(rr) || rr > rMax, continue; end
        ls = ':'; lw = 0.9;
        if ringD(ri), ls = '--'; lw = 1.5; end
        polarplot(ax, tt, rr*ones(size(tt)), ls, 'Color',[0.45 0.45 0.45], 'LineWidth',lw);
        text(ax, deg2rad(-6), rr, ringL{ri}, 'FontSize',7.5, ...
             'Color',[0.35 0.35 0.35], 'HorizontalAlignment','center', ...
             'BackgroundColor','w', 'Margin',0.5);
    end

    % arms first, dots on top
    for i = 1:numel(P(q).th)
        c  = gCol(P(q).g(i),:);
        sp = min(P(q).a1(i) + P(q).a3(i), 2*pi);      % never exceed one turn
        aa = linspace(P(q).th(i) - min(P(q).a1(i), sp), ...
                      P(q).th(i) + min(P(q).a3(i), sp), 60);
        polarplot(ax, aa, P(q).r(i)*ones(size(aa)), '-', ...
                  'Color',[c armAlpha], 'LineWidth',1.1);
    end
    for gg = 1:numel(groups)
        m = P(q).g == gg;
        if ~any(m), continue; end
        polarscatter(ax, P(q).th(m), P(q).r(m), sz(P(q).epc(m)), gCol(gg,:), 'filled', ...
                     'MarkerFaceAlpha',0.75, 'MarkerEdgeColor','w', 'LineWidth',0.4);
    end
    title(ax, sprintf('%s-triggered   |   %d tested, %d sig', ...
          upper(TRIG{q}), numel(P(q).th), nnz(P(q).sig)), 'FontWeight','bold');
end

% legend + key
lg = gobjects(numel(groups),1);
axl = axes('Parent',fh,'Position',[0 0 1 1],'Visible','off'); hold(axl,'on');
for gg = 1:numel(groups)
    n = nnz(P(2).g == gg);
    lg(gg) = scatter(axl, NaN, NaN, 60, gCol(gg,:), 'filled', ...
                     'DisplayName', sprintf('%s (n=%d)', groups{gg}, n));
end
legend(axl, lg, 'Orientation','horizontal', 'Box','off', ...
       'Position',[0.15 0.005 0.7 0.035], 'FontSize',9);

if strcmp(V,'z')
    note = sprintf(['radius = z, the peak height in null SDs -- it does NOT saturate. ' ...
        'z of the peak in epc is the SAME number (a constant divisor cancels in a z-score).']);
elseif strcmp(V,'logz')
    note = sprintf(['radius = log z, LOGARITHMIC SPACING but the tick labels are z itself. ' ...
        'Spreads the crowded z=%.0f-8 band where most significant cells sit; anything with ' ...
        'z <= %.1f is pinned at the centre.'], zCrit, zFloor);
else
    note = sprintf(['%d shifts -> p floor %.1e. SATURATED: %d/%d onset and %d/%d peak ' ...
        'significant cells sit ON that floor, so this radius cannot rank them -- see the z version.'], ...
        prm.nPerm, pFloor, nAtFloor(1), nSigQ(1), nAtFloor(2), nSigQ(2));
end
title(tl, {sprintf(['Event latency on the breath cycle   |   angle = latency / IBI (0 = trigger, 1 = next),   ' ...
                    'radius = %s,   arm = circular IQR,   dot area = events/cycle (%.2f-%.2f)'], ...
                    rLab, sLim(1), sLim(2)), note, ...
           sprintf(['alpha %.4g, RAW p with NO family-wise correction.   NO GCaMP lag compensation.   ' ...
                    'Angle is a LINEAR time normalisation, not piecewise breath phase.'], prm.alphaPerm)}, ...
      'FontWeight','bold','FontSize',10,'Interpreter','tex');

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir, ['polar_event_latency_' V]);
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    T = table();
    for q = 1:2
        k = P(q).keep(:);
        T = [T; table(repmat(TRIG(q),numel(k),1), {CELL(k).key}', {CELL(k).label}', ...
             groups(P(q).g(:))', P(q).th(:)/(2*pi), P(q).r(:), ...
             (P(q).a1(:)+P(q).a3(:))/(2*pi), P(q).epc(:), P(q).sig(:), ...
             'VariableNames',{'trigger','cell_key','label','group','frac_of_cycle', ...
                              'radius','iqr_frac_of_cycle','epc','sig'})]; %#ok<AGROW>
    end
    writetable(T,[base '.csv']);
    fprintf('saved %s.{png,pdf,csv}\n', base);
end

%% ---- what the figure shows, in numbers ----
if vi > 1, continue; end          % the numbers are radius-independent
for q = 1:2
    m = P(q).sig;
    if ~any(m), continue; end
    fprintf('\n%s, significant cells (n=%d):\n', TRIG{q}, nnz(m));
    fprintf('  fraction of cycle : median %.3f (IQR arm median %.2f of a cycle)\n', ...
            median(P(q).th(m))/(2*pi), median(P(q).a1(m)+P(q).a3(m))/(2*pi));
    for gg = 1:numel(groups)
        mm = m & P(q).g == gg;
        if ~any(mm), continue; end
        fprintf('  %-7s n=%2d | frac %.3f | arm %.2f cyc | epc %.2f\n', groups{gg}, nnz(mm), ...
                median(P(q).th(mm))/(2*pi), median(P(q).a1(mm)+P(q).a3(mm))/(2*pi), ...
                median(P(q).epc(mm)));
    end
end
end   % radial-variable loop
