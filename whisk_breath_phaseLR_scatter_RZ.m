function whisk_breath_phaseLR_scatter_RZ()
% whisk_breath_phaseLR_scatter_RZ
%
% Per-FRAME bilateral whisking phase scatter:  LEFT whisker instantaneous phase
% (phiL) vs RIGHT whisker instantaneous phase (phiR), one dot PER SAMPLE.
% If the two whiskers oscillate in phase, dots fall on the diagonal.
%
%   phiL = angle(hilbert(Lbp)),  phiR = angle(hilbert(Rbp))   [BP whisker angle]
%
% Every frame is plotted. By default frames are restricted to whisking epochs
% (USE_EPOCHS=true; set false for literally every frame including quiet periods,
% where phase is just noise). Optionally each frame is colored by breath phase
% (COLOR_BY_BREATH) so you can see whether L/R synchrony depends on the breath
% cycle. Detection (whisk band, epoch gate, inspiration onset) matches
% whisk_breath_raster_RZ / the xcorr script.
%
% Pooled over all included sessions. Nothing is written to disk.

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW    = 400;
LIK_THR = 0.6;

WHISK_BP  = [5 50];   % Hz, whisker band
BREATH_BP = [2 15];   % Hz, breathing band (only for breath-phase coloring)

EXPOSURE_OFFSET_S = 0.36324e-3;

USE_EPOCHS      = true;    % true = only frames inside whisking epochs; false = every frame
COLOR_BY_BREATH = true;    % true = color each frame by breath phase (0 = insp onset)
MARKER_SIZE     = 1.4;       % dot size (points); bigger = bolder dots

% whisk-epoch gating (global L+R envelope)
ENV_THR = 5; ENV_SMOOTH = 0.05; MIN_DUR = 1; MERGE_GAP = 0;

% inspiration onset (trough of inhale-up breath) -- for breath-phase coloring
MIN_INSP_DIST_S = 0.06; BREATH_PROM_MAD = 1.5; MAX_CYCLE_S = 1.0;

EXCLUDE      = "5840027";
EXCLUDE_SESS = {{"5916297",3}, {"5916297",4}, {"5916300",4}};

% --- saving: scatter layer is RASTERIZED (embedded image), axes/text/colorbar
%     stay VECTOR, so the PDF is small and opens in Illustrator with editable
%     labels. The light histogram stays fully vector.
SAVE_FIG  = false;
SAVE_DPI  = 300;                 % raster resolution of the dot layer
SAVE_PATH = "";                  % "" -> <repoRoot>/whisk_phaseLR_scatter.pdf
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
ord = 4;
[b_wb,a_wb] = butter(3,   WHISK_BP/(fpsW/2),  'bandpass');
[b_br,a_br] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

phiL = []; phiR = []; bph = [];     % deg per frame (pooled); bph = breath phase (deg) or NaN
nSess = 0;

fprintf('\n=== SESSION LOOP ===\n');
for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE)>0 && strcmp(animal,char(EXCLUDE)), continue; end
    if is_excluded_sess(animal,kRun,EXCLUDE_SESS), continue; end

    % -------- whisker --------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), LIK_THR);
    catch ME
        warning('whisk load failed %s n%d: %s', animal,kRun,ME.message); continue;
    end
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    tW = (0:numel(La)-1)'/fpsW;
    La0 = fillmissing(La(:),'linear','EndValues','nearest'); La0 = La0 - mean(La0,'omitnan');
    Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest'); Ra0 = Ra0 - mean(Ra0,'omitnan');
    xL = filtfilt(b_wb,a_wb, La0);
    xR = filtfilt(b_wb,a_wb, Ra0);

    % -------- breath: load, align to whisk grid, crop to overlap --------
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kRun)), dataRoot);
        if ~isfile(Pn.breath), warning('no breath file %s n%d',animal,kRun); continue; end
        Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
        tB = (0:numel(br)-1)'/fb + EXPOSURE_OFFSET_S;
        brw = interp1(tB, br, tW, 'linear', NaN);
    catch ME
        warning('breath load failed %s n%d: %s', animal,kRun,ME.message); continue;
    end
    keep = isfinite(brw);
    if nnz(keep) < fpsW, warning('overlap too short %s n%d',animal,kRun); continue; end
    xL = xL(keep); xR = xR(keep); brw = brw(keep);
    N  = numel(brw); tC = (0:N-1)'/fpsW;

    % -------- per-frame whisker phase --------
    pL = mod(rad2deg(angle(hilbert(xL))), 360);   % [0,360)
    pR = mod(rad2deg(angle(hilbert(xR))), 360);

    % -------- frame selection: whisking epochs (or all) --------
    if USE_EPOCHS
        env = movmean(0.5*(abs(hilbert(xL))+abs(hilbert(xR))), max(1,round(ENV_SMOOTH*fpsW)));
        ep  = bool_to_epochs(env > ENV_THR, tC, MIN_DUR, MERGE_GAP);
        sel = false(N,1); for q=1:size(ep,1), sel = sel | (tC>=ep(q,1)&tC<=ep(q,2)); end
    else
        sel = true(N,1);
    end
    if ~any(sel), warning('no frames selected %s n%d',animal,kRun); continue; end

    % -------- breath phase per frame (for coloring) --------
    cph = nan(N,1);
    if COLOR_BY_BREATH
        Bbp = filtfilt(b_br,a_br, brw);
        prom = BREATH_PROM_MAD * mad(Bbp,1); if ~isfinite(prom)||prom<=0, prom = 0.5*std(Bbp,'omitnan'); end
        [~,locs] = findpeaks(-Bbp,'MinPeakDistance',round(MIN_INSP_DIST_S*fpsW),'MinPeakProminence',prom);
        if numel(locs) >= 2, cph = breath_phase(tC, tC(locs), MAX_CYCLE_S); end
    end

    phiL = [phiL; pL(sel)]; phiR = [phiR; pR(sel)]; bph = [bph; cph(sel)]; %#ok<AGROW>
    nSess = nSess + 1;
    fprintf('%s n%d : %d frames (%.1f s)\n', animal, kRun, nnz(sel), nnz(sel)/fpsW);
end

assert(~isempty(phiL), 'no frames found');
fprintf('\nIncluded %d sessions, %d frames total\n', nSess, numel(phiL));

% -------- circular stats on per-frame L vs R phase --------
aL = deg2rad(phiL); aR = deg2rad(phiR);
dphi = wrapToPiLocal(aL - aR);
cmeanDiff = rad2deg(angle(mean(exp(1i*dphi))));
R_conc    = abs(mean(exp(1i*dphi)));
rho       = circ_corr_local(aL, aR);
frac45    = mean(abs(dphi) <= pi/4);     % fraction of frames with |phiL-phiR| <= 45 deg
fprintf('fraction within +/-pi/4 (45 deg): %.3f (%d/%d frames)\n', ...
    frac45, nnz(abs(dphi)<=pi/4), numel(dphi));

% ============================== PLOT ==================================
figure('Color','w','Position',[50 90 1640 500], 'Name','whisk phase L vs R (per frame)');
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Panel 1: plain scatter (no breath coloring)
ax0 = nexttile; hold(ax0,'on'); grid(ax0,'on');
plot(ax0, phiL, phiR, '.', 'Color',[0 0 0], 'MarkerSize',MARKER_SIZE);
plot(ax0, [0 360],[0 360], 'r-', 'LineWidth',1);
axis(ax0,'square'); xlim(ax0,[0 360]); ylim(ax0,[0 360]);
set(ax0,'XTick',0:90:360,'YTick',0:90:360);
xlabel(ax0,'LEFT whisker phase (deg)'); ylabel(ax0,'RIGHT whisker phase (deg)');
title(ax0, sprintf('phiL vs phiR per frame  (n=%d frames, %d sessions)\ncirc-corr \\rho=%.2f', ...
    numel(phiL), nSess, rho));

% Panel 2: same scatter, colored by breath phase
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
if COLOR_BY_BREATH && any(isfinite(bph))
    ok = isfinite(bph);
    scatter(ax1, phiL(ok), phiR(ok), MARKER_SIZE^2, bph(ok), '.');
    colormap(ax1, hsv); cb = colorbar(ax1); cb.Label.String = 'breath phase (deg, 0 = insp onset)';
    caxis(ax1,[0 360]);
    title(ax1, sprintf('colored by breath phase  (%d colored frames)', nnz(ok)));
else
    plot(ax1, phiL, phiR, '.', 'Color',[0.10 0.45 0.70], 'MarkerSize',MARKER_SIZE);
    title(ax1, 'colored by breath phase  (COLOR\_BY\_BREATH = false)');
end
plot(ax1, [0 360],[0 360], 'k-', 'LineWidth',1);
axis(ax1,'square'); xlim(ax1,[0 360]); ylim(ax1,[0 360]);
set(ax1,'XTick',0:90:360,'YTick',0:90:360);
xlabel(ax1,'LEFT whisker phase (deg)'); ylabel(ax1,'RIGHT whisker phase (deg)');

% Panel 3: L-R phase difference histogram
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
histogram(ax2, rad2deg(dphi), linspace(-180,180,37), 'Normalization','probability', ...
    'FaceColor',[0.10 0.45 0.70],'EdgeColor','none');
xline(ax2, 0, 'k--'); xline(ax2, cmeanDiff, 'r-', 'LineWidth',1.5);
xline(ax2, -45, 'k:'); xline(ax2, 45, 'k:');     % +/- pi/4 band
xlim(ax2,[-180 180]); set(ax2,'XTick',-180:90:180);
xlabel(ax2,'L-R whisker phase difference  \phi_L - \phi_R (deg)'); ylabel(ax2,'fraction of frames');
title(ax2, sprintf('L-R phase diff:  circ mean=%+.1f deg,  R=%.2f,  %.0f%% within \\pm45 deg', ...
    cmeanDiff, R_conc, 100*frac45));

title(tl, sprintf('Bilateral whisking phase synchrony (per frame)  (%d sessions)', nSess), ...
    'FontWeight','bold','Interpreter','none');

% -------- save: rasterize the heavy scatter layers, keep everything vector ----
if SAVE_FIG
    if strlength(SAVE_PATH)==0, SAVE_PATH = fullfile(repoRoot,'whisk_phaseLR_scatter.pdf'); end
    rasterize_markers(ax0, SAVE_DPI);     % plain scatter
    rasterize_markers(ax1, SAVE_DPI);     % breath-colored scatter
    exportgraphics(gcf, char(SAVE_PATH), 'ContentType','vector');
    fprintf('saved %s  (dots rasterized @%d dpi, axes/text vector)\n', char(SAVE_PATH), SAVE_DPI);
end
end

% ============================= HELPERS =============================
function rasterize_markers(ax, dpi)
% Replace the dense dot layers (scatter, or line with markers + no line) in ax
% by a single high-res RASTER image, so a vector PDF export stays small. Vector
% decorations (grid, diagonal, ticks, labels, colorbar) are untouched.
    sc = findobj(ax,'Type','scatter');
    ln = findobj(ax,'Type','line');
    ln = ln(arrayfun(@(h) strcmp(h.LineStyle,'none') && ~strcmp(h.Marker,'none'), ln));
    markers = [sc; ln];
    if isempty(markers), return; end

    xl = ax.XLim; yl = ax.YLim;
    % offscreen canvas sized to the ACTUAL axis (pixels), so dots keep their real
    % size; supersample SS for crisp edges. dpi/96 sets the supersample factor.
    ou = ax.Units; ax.Units = 'pixels'; p = ax.Position; ax.Units = ou;
    W = max(200, round(p(3))); H = max(200, round(p(4)));
    ftmp = figure('Visible','off','Color','w','Units','pixels','Position',[1 1 W H]);
    atmp = axes(ftmp,'Position',[0 0 1 1]); hold(atmp,'on');
    copyobj(markers, atmp);
    set(atmp,'XLim',xl,'YLim',yl,'YDir',ax.YDir,'Visible','off');
    colormap(atmp, colormap(ax)); atmp.CLim = ax.CLim;   % preserve color mapping

    tmpPng = [tempname '.png'];
    exportgraphics(atmp, tmpPng, 'Resolution', dpi);     % output ~ W*(dpi/96) px
    close(ftmp);
    img = imread(tmpPng); delete(tmpPng);

    % OPAQUE image (no AlphaData): transparency would force exportgraphics to
    % rasterize the WHOLE figure. Push the axes grid/ticks on top so they (and
    % the diagonal) stay vector, drawn over the opaque white-background raster.
    delete(markers);                                     % drop the vector dots
    ax.Layer = 'top';
    hImg = image(ax,'XData',xl,'YData',yl,'CData',flipud(img));
    uistack(hImg,'bottom');
end

function phi = breath_phase(tEv, tInsp, maxCycleS)
    tEv = tEv(:); phi = nan(size(tEv));
    k = discretize(tEv, tInsp); ok = ~isnan(k) & k < numel(tInsp);
    t0 = tInsp(k(ok)); t1 = tInsp(k(ok)+1); T = t1 - t0;
    p  = 360 .* (tEv(ok) - t0) ./ T; p(T > maxCycleS) = NaN; phi(ok) = p;
end

function w = wrapToPiLocal(a), w = mod(a+pi, 2*pi) - pi; end

function rho = circ_corr_local(a, b)
    abar = angle(mean(exp(1i*a))); bbar = angle(mean(exp(1i*b)));
    sa = sin(a-abar); sb = sin(b-bbar);
    rho = sum(sa.*sb) / sqrt(sum(sa.^2)*sum(sb.^2));
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a=logical(a(:)); d=diff([false;a;false]); s=find(d==1); e=find(d==-1)-1;
    ep=[t(s),t(e)]; if isempty(ep), return; end
    ep2=ep(1,:);
    for i=2:size(ep,1)
        if ep(i,1)-ep2(end,2)<=mergeGap, ep2(end,2)=ep(i,2); else, ep2(end+1,:)=ep(i,:); end %#ok<AGROW>
    end
    ep=ep2; ep=ep(ep(:,2)-ep(:,1)>=minDur,:);
end

function tf = is_excluded_sess(animal, kRun, excl)
    tf=false;
    for q=1:numel(excl)
        if strcmp(char(excl{q}{1}),char(animal)) && double(excl{q}{2})==kRun, tf=true; return; end
    end
end

function S = list_sessions(dirPath)
    d=dir(fullfile(char(dirPath),'*_whisk_n*DLC*.csv')); S={}; key={};
    for i=1:numel(d)
        tok=regexp(d(i).name,'^(\d+)_whisk_n(\d+)','tokens','once'); if isempty(tok), continue; end
        k=sprintf('%s_%s',tok{1},tok{2}); if any(strcmp(key,k)), continue; end
        key{end+1}=k; S{end+1}={tok{1},str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an=cellfun(@(c) str2double(c{1}),S); rn=cellfun(@(c) c{2},S);
        [~,o]=sortrows([an(:) rn(:)]); S=S(o);
    end
end

function csv = pick_csv(dirPath, prefix)
    d=dir(fullfile(char(dirPath),[char(prefix) '*DLC*.csv']));
    assert(~isempty(d),'no DLC csv matching %s* in %s',prefix,dirPath);
    [~,ix]=max(arrayfun(@(x) bestnum(x.name),d)); csv=fullfile(d(ix).folder,d(ix).name);
end
function n = bestnum(name)
    tok=regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
