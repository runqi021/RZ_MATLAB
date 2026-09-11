function falloff_gcamp_260730(dataDir, roiTag, fitDepth, powerLut, zSurface, surfMode)
% FALLOFF_GCAMP_260730  Optical falloff + resliced projections from ONE ramped
% GCaMP stack (260730_shiver_vgat_g8m), with the pia found automatically.
%
% Written for a session that does not fit either of the two earlier acquisition
% styles, so it is a standalone script rather than a row in
% falloff_run_all_260805.m:
%
%   - ONE stack per FOV, power ramped by pzAdjust (10 -> 50 %,
%     lengthConstants = -292 um), 471 slices, 50 frames/slice, 1 um step,
%     motor z +20 -> -450, zoom 1.2, single saved channel (GCaMP, ch 1).
%   - The filename carries no z token, so depth comes from the MOTOR z that was
%     SET (stackZStartPos + k*stackZStepSize), never from per-slice zs.
%   - Z_SURFACE (motor z of the pia) is MEASURED, not assumed.  Two things make
%     the earlier sessions' rule ("pia = where raw (top - bot) peaks") wrong
%     here, and both were checked on this stack:
%       * the power RAMPS 10 -> 50 %, so raw signal climbs with depth by
%         construction and its peak lands at motor z = -256, 250 um into the
%         tissue.  The search therefore runs on (top - bot)/mW^2.
%       * the label is SOMA GCaMP, not vessel dye.  Pial vessels sit exactly at
%         the surface so their peak IS the pia; a soma label is sparse in L1, so
%         its peak sits below the surface (here motor -13).  The default marker
%         is instead the HALF-RISE of the normalised curve -- where the signal
%         has climbed half way from the above-tissue value to the plateau, i.e.
%         the tissue/water interface the beam actually enters.
%     SURF_MODE picks between them, both are always printed, and passing
%     zSurface overrides the search entirely.
%
% The stack is read ONCE and produces both standard deliverables:
%   falloff_gcamp_<tag>.png/.pdf/.mat        4-row intensity-vs-depth figure
%   falloff_gcamp_proj_<tag>.png/.pdf        depth-corrected MAX projections
% plus two caches next to the data (delete either to force a re-read):
%   falloff_gcamp_<tag>_metrics.mat          per-slice top/bot/histograms
%   falloff_gcamp_<tag>_vol.mat              raw 1x1x1 um volume + per-plane %
%
% ANATOMY OF THIS SESSION (given by Runqi, not derived):
%   image ROW    (+y) -> more LATERAL
%   image COLUMN (+x) -> CAUDAL,  -x -> ROSTRAL
% so the ML axis is the image row and the AP axis is the image column:
%   coronal  view = ML vs depth = max along the AP axis (columns), medial -> lateral
%   sagittal view = AP vs depth = max along the ML axis (rows),    rostral -> caudal
%
% USAGE
%   falloff_gcamp_260730
%   falloff_gcamp_260730('C:\260730_shiver_vgat_g8m', 'roi1')

%% --------------------------- USER PARAMETERS -------------------------------
DATA_DIR   = 'C:\260730_shiver_vgat_g8m';
ROI_TAG    = 'roi1';

POWER_LUT  = 'post_260723';  % acquired 2026-07-30
P_EXPONENT = 2;              % two-photon: signal ~ power^2
FIT_DEPTH  = [10 500];       % shiverer window (this mouse is a shiverer)
PCT_FRAC   = 0.05;           % top / bottom 5 % of pixels, as a FRACTION
Z_SURFACE  = [];             % [] = find the pia from the data (see above)
SURF_MODE  = 'halfrise';     % 'halfrise' = tissue/water interface (default for a
                             %              soma label), 'peak' = brightest plane
                             %              (the vessel-session rule)
SURF_SMOOTH = 5;             % um of moving-mean smoothing before the search

PIX_BASE   = 1.7778;         % um/px at zoom 1 for 512 px lines (rig constant)
PIX_TARGET = 1.0;            % isotropic target voxel, um
BIN_UM     = 100;            % depth bin for the x-y projections
DEPTH_MAX  = 600;            % FIXED depth axis in every falloff figure
CLIM_PCT   = [0.5 99.5];     % display window, percentiles per panel
FLOOR_FRAC = 0.05;           % per-plane floor = mean of the dimmest fraction

PLANE_DEPTHS = 0:50:400;     % single-plane figure: one panel per depth, each the
                             % 50-frame average of the NEAREST slice.  Raw ADU,
                             % linear (no gamma), each panel clipped to CLIM_PCT.
% Image pixel -> stage micron.  MEASURED rig property, not a guess: an image
% COLUMN runs along +stage x and an image ROW along -stage y
% ([[reference_stage_axis_convention]], 8-way phase-correlation test on Vgat/0730).
% Anatomy on top of that, as given by Runqi for this session: down the rows
% (+image y) is LATERAL, right along the columns (+image x) is CAUDAL.
ROW_SIGN_Y = -1;             % image row  -> ROW_SIGN_Y * stage y
COL_SIGN_X = +1;             % image col  -> COL_SIGN_X * stage x

DIST_SPEC   = [50 100];      % row 1: F distribution at 50, 150, 250 ... um
DEPTH_CLIM  = [0 600];       % FIXED colour range, depth
POWER_CLIM  = [10 350];      % FIXED colour range, power (mW)
DEPTH_XLIM  = [-50 600];
DEPTH_XTICK = 0:100:600;
ADU_XLIM    = [-50 32767];
SAVE_FIG    = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(dataDir),  DATA_DIR  = dataDir;  end
if nargin >= 2 && ~isempty(roiTag),   ROI_TAG   = roiTag;   end
if nargin >= 3 && ~isempty(fitDepth), FIT_DEPTH = fitDepth; end
if nargin >= 4 && ~isempty(powerLut), POWER_LUT = powerLut; end
if nargin >= 5 && ~isempty(zSurface), Z_SURFACE = zSurface; end
if nargin >= 6 && ~isempty(surfMode), SURF_MODE = surfMode; end

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);              % laser_power_calibration lives here
addpath(fileparts(thisDir));   % repo root, unconditionally

f = dir(fullfile(DATA_DIR, [ROI_TAG '*.tif']));
assert(numel(f) == 1, ['Expected exactly ONE stack matching ''%s'' under %s, ' ...
    'found %d -- this script is for the single-ramped-stack layout.'], ...
    ROI_TAG, DATA_DIR, numel(f));
TIF = fullfile(f.folder, f.name);
TAG = regexprep(ROI_TAG, '[\*\.\s/\\]+', '_');

% Software-quality OpenGL: the NVIDIA path silently writes blank 9 kB PNGs on
% roughly 1 % of exportgraphics calls.  See [[reference_matlab_opengl_export_stubs]].
try, opengl('software'); catch, end  %#ok<CTCH>

%% --- metadata (page 1 only -- imfinfo on 23550 IFDs is minutes and GBs) ------
meta = read_si_meta(TIF);
nSl  = si_num(meta, "SI.hStackManager.numSlices");
fpsl = si_num(meta, "SI.hStackManager.framesPerSlice");
z0   = si_num(meta, "SI.hStackManager.stackZStartPos");
dzS  = si_num(meta, "SI.hStackManager.stackZStepSize");   % SIGNED, -1 here
P0   = si_num(meta, "SI.hBeams.stackStartPower");
P1   = si_num(meta, "SI.hBeams.stackEndPower");
LC   = si_num(meta, "SI.hBeams.lengthConstants");
zRef = si_num(meta, "SI.hStackManager.zPowerReference");
pzOn = startsWith(lower(strtrim(si_line(meta, "SI.hBeams.pzAdjust"))), "true");
zoom = si_num(meta, "SI.hRoiManager.scanZoomFactor");
ppl  = si_num(meta, "SI.hRoiManager.pixelsPerLine");
chSv = si_vec(meta, "SI.hChannels.channelSave");
nCh  = max(numel(chSv), 1);
mPos = si_vec(meta, "SI.hMotors.motorPosition");   % [x y z] of the FOV CENTRE
assert(nCh == 1, ['%s saves %d channels (%s); this script assumes the single ' ...
    'GCaMP channel.'], f.name, nCh, mat2str(chSv(:).'));

% Motor z that was SET, uniform by construction.  Per-slice SI zs is not used:
% its steps come back non-uniform and its zero moves between acquisitions.
zSet = z0 + (0:nSl-1).' * dzS;
if ~isfinite(zRef), zRef = zSet(1); end
ramped = pzOn && ~(isfinite(P0) && isfinite(P1) && abs(P1-P0) < 1e-9);
if ramped
    assert(isfinite(LC) && LC ~= 0, 'ramp on but lengthConstants unusable');
    pctZ = P0 .* exp((zSet - zRef) ./ LC);
else
    pctZ = repmat(P0, nSl, 1);
end
pxUm = PIX_BASE * (512/ppl) / zoom;

fprintf('\n%s\n', TIF);
fprintf('  %d slices x %d frames, 1 channel, zoom %g, %.4f um/px lateral\n', ...
    nSl, fpsl, zoom, pxUm);
fprintf('  motor z %+g -> %+g, step %g um\n', zSet(1), zSet(end), dzS);
if ramped
    fprintf('  power RAMPED %g -> %g %% (Lz = %g um, ref z = %+g)\n', ...
        pctZ(1), pctZ(end), LC, zRef);
else
    fprintf('  power flat at %g %%\n', P0);
end

%% --- read: per-slice metrics AND the isotropic volume, in ONE pass ----------
mCache = fullfile(DATA_DIR, sprintf('falloff_gcamp_%s_metrics.mat', TAG));
vCache = fullfile(DATA_DIR, sprintf('falloff_gcamp_%s_vol.mat', TAG));
if isfile(mCache) && isfile(vCache)
    fprintf('using caches (delete to re-read):\n  %s\n  %s\n', mCache, vCache);
    M = load(mCache);  top = M.top;  bot = M.bot;  Hst = M.Hst;
    hCenters = M.hCenters;  satFrac = M.satFrac;  satLev = M.satLev;
    Q = load(vCache);   V = Q.V;
else
    [top, bot, Hst, hCenters, satFrac, satLev, V] = ...
        read_stack(TIF, nSl, fpsl, PCT_FRAC, pxUm, PIX_TARGET);
    save(mCache, 'top', 'bot', 'Hst', 'hCenters', 'satFrac', 'satLev', ...
        'zSet', 'pctZ', 'TIF');
    save(vCache, 'V', 'zSet', 'pctZ', 'pxUm', 'PIX_TARGET', 'TIF', '-v7.3');
    fprintf('cached -> %s\n         %s\n', mCache, vCache);
end
diffRaw = top - bot;

%% --- power ------------------------------------------------------------------
% Comes BEFORE the pia search: with a 10 -> 50 % ramp the raw signal climbs with
% depth whatever the tissue does, so the surface can only be found on the
% power-normalised curve.
[mW, cal] = laser_power_calibration(pctZ, POWER_LUT);
fprintf('\npower calibration: %s  %g-%g %% -> %.1f-%.1f mW%s\n', cal.name, ...
    pctZ(1), pctZ(end), mW(1), mW(end), tern(any(cal.clamped), '  [CLAMPED]', ''));

diffN = diffRaw ./ (mW(:) .^ P_EXPONENT);

%% --- the pia, measured ------------------------------------------------------
% Both markers are always computed and printed, because they answer different
% questions and on a soma label they do not coincide:
%   'peak'     the brightest plane of (top - bot)/mW^2.  This is the rule used on
%              the vessel sessions, where the pial vessels sit ON the surface.
%              With GCaMP in cell-sparse L1 it lands BELOW the pia.
%   'halfrise' where the normalised curve has climbed half way from its
%              above-tissue value (the shallowest slice) to that peak -- the
%              blurred edge of the tissue itself, which is where the beam enters
%              and therefore where depth 0 belongs for an attenuation measurement.
w  = max(1, round(SURF_SMOOTH / abs(dzS)));
sm = movmean(diffN, w);
[pk, kPeak] = max(sm);
base  = sm(1);                       % shallowest slice = above tissue
half  = (pk + base) / 2;
kHalf = find(sm(1:kPeak) <= half, 1, 'last');
if isempty(kHalf)
    zHalf = zSet(1);  kHalf = 1;
elseif kHalf < kPeak                 % linear interpolation between the two slices
    f = (half - sm(kHalf)) / (sm(kHalf+1) - sm(kHalf));
    zHalf = zSet(kHalf) + f * (zSet(kHalf+1) - zSet(kHalf));
else
    zHalf = zSet(kHalf);
end
fprintf('\npia search on (top-bot)/mW^2, %g um moving mean\n', SURF_SMOOTH);
fprintf('  above-tissue value %.2f at motor z %+g, plateau %.2f at motor z %+g\n', ...
    base, zSet(1), pk, zSet(kPeak));
fprintf('  half-rise  motor z = %+.1f\n', zHalf);
fprintf('  peak       motor z = %+g\n', zSet(kPeak));
lo = max(1, kHalf-8);  hi = min(nSl, kHalf+8);
fprintf('  %8s %12s %12s\n', 'motor z', '(t-b)/mW^2', 'smoothed');
for k = lo:hi
    fprintf('  %+8g %12.2f %12.2f%s\n', zSet(k), diffN(k), sm(k), ...
        tern(k == kHalf, '   <- half-rise', ''));
end
if isempty(Z_SURFACE)
    switch lower(SURF_MODE)
    case 'halfrise', Z_SURFACE = zHalf;
    case 'peak',     Z_SURFACE = zSet(kPeak);
    otherwise, error('SURF_MODE must be ''halfrise'' or ''peak''');
    end
    fprintf('  USING %s: pia at motor z = %+.1f (depth 0)\n', SURF_MODE, Z_SURFACE);
    if kPeak <= 2 || kPeak >= nSl-1
        warning(['The normalised peak is at the very %s of the stack -- the pia ' ...
            'may be outside it.  Check before trusting the depth axis.'], ...
            tern(kPeak <= 2, 'start', 'end'));
    end
else
    fprintf('  pia SUPPLIED: motor z = %+g (search result ignored)\n', Z_SURFACE);
end
depth = Z_SURFACE - zSet;     % motor z is inverted: more negative = deeper

%% --- fit --------------------------------------------------------------------
fitMask = depth >= FIT_DEPTH(1) & depth <= FIT_DEPTH(2) & diffN > 0;
[ell, A, r2] = fit_exp(depth(fitMask), diffN(fitMask));
fprintf('\n(top %g%% - bot %g%%) / mW^%g, fit over %g-%g um (%d slices)\n', ...
    100*PCT_FRAC, 100*PCT_FRAC, P_EXPONENT, FIT_DEPTH(1), FIT_DEPTH(2), nnz(fitMask));
fprintf('  l = %.1f um   (EAL = 2l = %.1f um)   R^2 = %.4f\n', ell, 2*ell, r2);
fprintf('  saturated slices (>0.1%% of pixels at %g): %d of %d\n', ...
    satLev, nnz(satFrac > 1e-3), nSl);

%% --- attenuation F(0)/F(z) ---------------------------------------------------
kRef = find(abs(depth) < 0.5, 1);
assert(~isempty(kRef), 'No slice at depth 0 -- F(0) is not measured.');
F0raw = diffRaw(kRef);  F0n = diffN(kRef);
fprintf('\nF(0) at depth 0 (motor z %+g): raw %.1f ADU, /mW^2 %.4g\n', ...
    zSet(kRef), F0raw, F0n);
good = diffRaw > 0;
att  = nan(nSl,1);  att(good)  = F0raw ./ diffRaw(good);
attP = nan(nSl,1);  attP(good) = F0n   ./ diffN(good);

mA = fitMask & isfinite(attP) & attP > 0;
pA = polyfit(depth(mA), log(attP(mA)), 1);
POOLatt = struct('d', depth(mA), 'A', exp(pA(2)), 'L', 1/pA(1), ...
    'r2', 1 - sum((log(attP(mA)) - polyval(pA, depth(mA))).^2) / ...
              sum((log(attP(mA)) - mean(log(attP(mA)))).^2));
POOL = struct('d', depth(fitMask), 'A', A, 'L', ell, 'r2', r2);

%% ============================ FIGURE 1: falloff =============================
nC = 256;  cmap = jet(nC);
ramp = @(v, lim) cmap(min(nC, max(1, round(1 + (nC-1) * ...
    (v - lim(1)) / diff(lim)))), :);
col = ramp(median(mW), POWER_CLIM);
lblP = sprintf('%.0f-%.0f mW', mW(1), mW(end));

spec = { ...
 top,   sprintf('top %g%% (ADU)',100*PCT_FRAC),  sprintf('raw TOP %g%% of pixels -- faded = above pia', 100*PCT_FRAC), 'linear', [],                                'stacks'
 bot,   sprintf('bot %g%% (ADU)',100*PCT_FRAC),  sprintf('raw BOT %g%% of pixels', 100*PCT_FRAC),                      'linear', [],                                'none'
 diffRaw, sprintf('top %g%% - bot %g%% (ADU)',100*PCT_FRAC,100*PCT_FRAC), 'raw TOP - BOT -- floor removed',            'linear', [],                                'none'
 diffN, sprintf('(top - bot) / mW^%g',P_EXPONENT), '(TOP - BOT) normalised by power',                                  'log',    @(ax) draw_fit(ax, POOL, -1, true), 'fit'
 att,   'attenuation  F(0) / F(z)', 'ATTENUATION (top - bot), raw, single reference F(0) at z = 0',                    'linear', [],                                'none'
 attP,  'attenuation  F(0) / F(z)', 'ATTENUATION (top - bot), power-normalised, F(0) at z = 0',                        'log',    @(ax) draw_fit(ax, POOLatt, +1, false), 'none' };

fg = figure('Color','w','Position',[40 10 1400 1550],'Name','falloff');
tl = tiledlayout(fg, 4, 2, 'TileSpacing','compact','Padding','compact');

ax = nexttile([1 2]); hold(ax,'on');
dWant = DIST_SPEC(1) : DIST_SPEC(2) : max(depth);
for j = 1:numel(dWant)
    [~, k] = min(abs(depth - dWant(j)));
    h = double(Hst(k,:));  m = h > 0;  h(~m) = NaN;   % break the line at empty bins
    plot(ax, hCenters, h, '-', 'Color', ramp(dWant(j), DEPTH_CLIM), 'LineWidth', 1.1);
end
set(ax, 'YScale','log');  xlim(ax, ADU_XLIM);  grid(ax,'on'); box(ax,'on');
colormap(ax, cmap);  caxis(ax, DEPTH_CLIM);
cb = colorbar(ax);  cb.Label.String = 'depth below pia (\mum)';  cb.Direction = 'reverse';
xlabel(ax, 'pixel value F (ADU)');  ylabel(ax, 'count');
title(ax, sprintf('F distribution at %g:%g:%g \\mum', DIST_SPEC(1), DIST_SPEC(2), ...
    max(dWant)), 'FontWeight','normal');

up = depth < 0;    % above the pia: shown, never fitted
for q = 1:size(spec,1)
    [y, ylab, sub, ys, fitFcn, legMode] = spec{q,:};
    ax = nexttile; hold(ax,'on');
    plot(ax, depth(up), y(up), '.', 'Color', (col+2)/3, 'MarkerSize', 6, ...
        'HandleVisibility','off');
    plot(ax, depth(~up), y(~up), '.', 'Color', col, 'MarkerSize', 9, ...
        'DisplayName', lblP, 'HandleVisibility', tern(strcmp(legMode,'stacks'),'on','off'));
    if ~isempty(fitFcn), fitFcn(ax); end
    set(ax, 'YScale', ys);
    xlim(ax, DEPTH_XLIM);  set(ax, 'XTick', DEPTH_XTICK);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'depth (\mum)');  ylabel(ax, ylab);
    title(ax, sub, 'FontWeight','normal');
    if ~strcmp(legMode,'none'), legend(ax, 'Location','southwest', 'Box','off'); end
    if strcmp(legMode,'stacks')
        colormap(ax, cmap);  caxis(ax, POWER_CLIM);
        cbp = colorbar(ax);  cbp.Label.String = 'laser power (mW)';
    end
end
title(tl, {sprintf('%s  --  %s  --  GCaMP optical falloff', DATA_DIR, ROI_TAG), ...
    sprintf('pia at motor z %+.1f (%s of (top-bot)/mW^2)  |  LUT %s, signal ~ mW^%g  |  fit %g-%g um', ...
        Z_SURFACE, SURF_MODE, cal.name, P_EXPONENT, FIT_DEPTH(1), FIT_DEPTH(2))}, ...
    'FontWeight','bold','Interpreter','none');

if SAVE_FIG
    base = fullfile(DATA_DIR, sprintf('falloff_gcamp_%s', TAG));
    exportgraphics(fg, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(fg, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'top', 'bot', 'diffRaw', 'diffN', 'depth', 'zSet', ...
        'pctZ', 'mW', 'ell', 'A', 'r2', 'Z_SURFACE', 'SURF_MODE', 'zHalf', ...
        'kPeak', 'FIT_DEPTH', 'POWER_LUT', ...
        'P_EXPONENT', 'PCT_FRAC', 'satFrac', 'satLev', 'TIF');
    fprintf('\nsaved %s.png/.pdf/.mat\n', base);
end

%% --- stage coordinates of the 1 um grid --------------------------------------
% One value per pixel of the resampled volume, from the FOV centre in
% SI.hMotors.motorPosition.  These are REAL stage microns, so two FOVs of the
% same session can be read against each other directly.
nR = size(V,1);  nCol = size(V,2);
axML = (0:nR-1)   * PIX_TARGET;      % panel coordinate along the rows (ML)
axAP = (0:nCol-1) * PIX_TARGET;      % panel coordinate along the cols (AP)
stgY = mPos(2) + ROW_SIGN_Y * (axML - (nR-1)  /2*PIX_TARGET);   % rows -> stage y
stgX = mPos(1) + COL_SIGN_X * (axAP - (nCol-1)/2*PIX_TARGET);   % cols -> stage x
fprintf('\nFOV centre (SI.hMotors.motorPosition): x %+.1f  y %+.1f\n', mPos(1), mPos(2));
fprintf('  stage x across the columns: %+.1f -> %+.1f  (rostral -> caudal)\n', stgX(1), stgX(end));
fprintf('  stage y down    the rows  : %+.1f -> %+.1f  (medial  -> lateral)\n', stgY(1), stgY(end));

%% ===================== FIGURE 2: single planes, raw ==========================
% One panel per requested depth, each the 50-frame average of THAT ONE slice.
% Raw ADU, linear grey (no gamma, no depth correction), clipped per panel to
% CLIM_PCT -- so the panels show what the microscope actually collected and
% deeper ones are not made to look like shallow ones.
fs = figure('Color','w','Position',[60 40 1450 1500],'Name','planes');
nP  = numel(PLANE_DEPTHS);
nRw = ceil(sqrt(nP));  nCl = ceil(nP/nRw);
ts  = tiledlayout(fs, nRw, nCl, 'TileSpacing','compact','Padding','compact');
fprintf('\nsingle planes (50-frame average, raw ADU)\n');
for j = 1:nP
    [~, kP] = min(abs(depth - PLANE_DEPTHS(j)));
    ax = nexttile;
    show(ax, V(:,:,kP), CLIM_PCT, axAP, axML);
    axis(ax, 'image');
    title(ax, sprintf('z = %g \\mum   (%g %% lp, %.0f mW)', PLANE_DEPTHS(j), ...
        round(pctZ(kP),1), mW(kP)), 'FontWeight','normal');
    if j == 1
        ylabel(ax, 'medial \rightarrow lateral');
        xlabel(ax, 'rostral \rightarrow caudal');
        add_scalebar(ax, axAP(end), axML(end), 100);
    end
    fprintf('  %4g um -> slice at depth %6.1f (motor z %+g), %5.1f %%, %.0f mW\n', ...
        PLANE_DEPTHS(j), depth(kP), zSet(kP), pctZ(kP), mW(kP));
end
title(ts, {sprintf('%s  --  %s  --  single planes, 50-frame average, RAW', DATA_DIR, ROI_TAG), ...
    sprintf('%g um isotropic, bicubic   |   linear grey, no gamma   |   each panel clipped %g-%g%%   |   pia at motor z %+.1f', ...
        PIX_TARGET, CLIM_PCT(1), CLIM_PCT(2), Z_SURFACE), ''}, ...
    'FontWeight','bold','Interpreter','none','FontSize',10);
if SAVE_FIG
    base = fullfile(DATA_DIR, sprintf('falloff_gcamp_planes_%s', TAG));
    exportgraphics(fs, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(fs, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('saved %s.png/.pdf\n', base);
end

%% ========================= FIGURE 3: projections ============================
% (a) floor out FIRST -- l was measured on (top - bot), and the floor itself
%     grows with depth, so amplifying it brightens the deep end for a reason
%     that has nothing to do with cells.
nFloor = max(1, round(FLOOR_FRAC * size(V,1) * size(V,2)));
for k = 1:size(V,3)
    pl = V(:,:,k);
    V(:,:,k) = pl - mean(mink(pl(:), nFloor));
end
% (b) gain = exp(z/l)/mW^2, clamped beyond the fitted range where the signal
%     stops decaying and the exponential would just explode the background.
zGain = min(depth, FIT_DEPTH(2));
if any(depth > FIT_DEPTH(2))
    fprintf('  gain clamped beyond %g um for %d of %d planes\n', ...
        FIT_DEPTH(2), nnz(depth > FIT_DEPTH(2)), numel(depth));
end
V = V .* reshape(single(exp(zGain / ell) ./ (mW(:) .^ P_EXPONENT)), 1, 1, []);

edges = 0 : BIN_UM : DEPTH_MAX;
nBin  = numel(edges) - 1;
XY = cell(nBin,1);  lbl = cell(nBin,1);
fprintf('\nx-y MAX projections\n');
for b = 1:nBin
    m = depth >= edges(b) & depth < edges(b+1);
    lbl{b} = sprintf('%g-%g \\mum', edges(b), edges(b+1));
    if any(m), XY{b} = max(V(:,:,m), [], 3); else, XY{b} = []; end
    fprintf('  bin %-14s %4d planes\n', lbl{b}, nnz(m));
end

% ANATOMY (given for this session):  image ROW +y = more LATERAL,
%                                    image COL +x = CAUDAL.
% So ML is dimension 1 and AP is dimension 2, and:
%   coronal  = ML vs depth = max along AP (dim 2)
%   sagittal = AP vs depth = max along ML (dim 1)
inT = depth >= 0 & depth <= DEPTH_MAX;
COR = squeeze(max(V(:,:,inT), [], 2));   % (ML, z)  medial -> lateral down the rows
SAG = squeeze(max(V(:,:,inT), [], 1));   % (AP, z)  rostral -> caudal along the cols
dT  = depth(inT);

fp = figure('Color','w','Position',[30 30 max(1500, 300*nBin) 1000]);
tl2 = tiledlayout(fp, 2, nBin, 'TileSpacing','compact','Padding','compact');
tl2.Title.FontSize = 10;

for b = 1:nBin
    ax = nexttile(b);
    if isempty(XY{b})
        set(ax, 'XTick',[], 'YTick',[], 'Color',[0.93 0.93 0.93]);
        text(ax, 0.5, 0.5, 'no data', 'Units','normalized', ...
            'HorizontalAlignment','center', 'Color',[0.45 0.45 0.45]);
        box(ax,'on');
    else
        show(ax, XY{b}, CLIM_PCT, axAP, axML);   % x = AP (cols), y = ML (rows)
        axis(ax,'image');
        if b == 1
            ylabel(ax, 'medial \rightarrow lateral (\mum)');
            xlabel(ax, 'rostral \rightarrow caudal (\mum)');
            add_scalebar(ax, axAP(end), axML(end), 100);
        end
    end
    title(ax, lbl{b}, 'FontWeight','normal', 'FontSize', 10);
end

% The lateral axis of each side view carries the REAL STAGE POSITION, not a
% distance from the corner: stage y across the coronal view, stage x across the
% sagittal one.  The ticks therefore run in whichever direction the stage does --
% down the rows is -stage y, which for this prep is the LATERAL direction.
nHalf = max(1, floor(nBin/2));
ax = nexttile(nBin+1, [1 nHalf]);
show(ax, COR.', CLIM_PCT, axML, dT);      % rows = depth, cols = ML
side_axes(ax, 'stage y (\mum),  medial \rightarrow lateral', axML, stgY, DEPTH_MAX, BIN_UM);
title(ax, 'coronal view  (ML vs depth, max along AP)', 'FontWeight','normal');

ax = nexttile(nBin+1+nHalf, [1 nBin-nHalf]);
show(ax, SAG.', CLIM_PCT, axAP, dT);
side_axes(ax, 'stage x (\mum),  rostral \rightarrow caudal', axAP, stgX, DEPTH_MAX, BIN_UM);
title(ax, 'sagittal view  (AP vs depth, max along ML)', 'FontWeight','normal');

title(tl2, {sprintf('%s  --  %s  --  GCaMP MAX projections', DATA_DIR, ROI_TAG), ...
    sprintf(['depth-corrected by exp(z/l), l = %.1f um   |   %g um isotropic, bicubic   |   ' ...
             'display %g-%g%%   |   pia at motor z %+.1f'], ell, PIX_TARGET, ...
             CLIM_PCT(1), CLIM_PCT(2), Z_SURFACE), ''}, ...
    'FontWeight','bold','Interpreter','none','FontSize',10);

if SAVE_FIG
    base = fullfile(DATA_DIR, sprintf('falloff_gcamp_proj_%s', TAG));
    exportgraphics(fp, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(fp, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('saved %s.png/.pdf\n', base);
end
end

% ============================== read =========================================

function [pTop, pBot, Hst, hCenters, satFrac, satLev, V] = ...
    read_stack(TIF, nSl, fpsl, PCT_FRAC, pxUm, PIX_TARGET)
% ONE sequential pass over every page.  Pages are walked with nextDirectory
% rather than setDirectory(t,k): the IFD chain is walked from the start on every
% setDirectory, which on a 23550-page file is quadratic.

t = Tiff(TIF, 'r');  cl = onCleanup(@() t.close());
I0 = t.read();
H = size(I0,1);  W = size(I0,2);
satLev = double(intmax(class(I0)));
nPct = max(1, round(PCT_FRAC * H * W));

hEdges = -1024 : 8 : 32768;      % fixed bins, shared by every slice
Hst = zeros(nSl, numel(hEdges)-1, 'single');
hCenters = (hEdges(1:end-1) + hEdges(2:end)) / 2;
pTop = nan(nSl,1);  pBot = nan(nSl,1);  satFrac = nan(nSl,1);

sc = pxUm / PIX_TARGET;                                    % >1 upsamples
probe = imresize(zeros(H, W, 'single'), sc, 'bicubic');
V = zeros(size(probe,1), size(probe,2), nSl, 'single');
fprintf('\nreading %d slices x %d frames -> volume %dx%dx%d @ %g um/px\n', ...
    nSl, fpsl, size(V,1), size(V,2), nSl, PIX_TARGET);

tS = tic;  first = true;
for k = 1:nSl
    acc = zeros(H, W);  nSat = 0;
    for fr = 1:fpsl
        if first, first = false; else, t.nextDirectory(); end
        F = double(t.read());
        nSat = nSat + nnz(F >= satLev);
        acc = acc + F;
    end
    I = acc / fpsl;
    satFrac(k) = nSat / (fpsl*H*W);
    v = I(:);
    pTop(k) = mean(maxk(v, nPct));
    pBot(k) = mean(mink(v, nPct));
    Hst(k,:) = single(histcounts(v, hEdges));
    V(:,:,k) = imresize(single(I), sc, 'bicubic');          % RAW ADU
    if mod(k, max(1,floor(nSl/10))) == 0 || k == nSl
        fprintf('  %3d/%3d  top %9.1f  bot %8.1f  sat %.3f%%  (%.0f s)\n', ...
            k, nSl, pTop(k), pBot(k), 100*satFrac(k), toc(tS));
    end
end
fprintf('  read in %.1f s\n', toc(tS));
end

% ============================== fitting ======================================

function [L, A, r2] = fit_exp(d, y)
% y ~ A*exp(-d/L), OLS on log y.
m = isfinite(d) & isfinite(y) & y > 0;
d = d(m); y = y(m);
if numel(d) < 5, L = NaN; A = NaN; r2 = NaN; return; end
p  = polyfit(d, log(y), 1);
L  = -1/p(1);
A  = exp(p(2));
res = log(y) - polyval(p, d);
r2 = 1 - sum(res.^2) / sum((log(y) - mean(log(y))).^2);
end

% ============================== plotting =====================================

function draw_fit(ax, P, sgn, showLeg)
dF = linspace(min(P.d), max(P.d), 200);
h = plot(ax, dF, P.A*exp(sgn*dF/P.L), 'k--', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('linear fit:  l = %.1f \\mum,  R^2 = %.3f', P.L, P.r2));
if ~showLeg, set(h, 'HandleVisibility', 'off'); end
end

function show(ax, I, pct, xc, yc)
I = double(I);
v = I(isfinite(I));
lo = prctile(v, pct(1));  hi = prctile(v, pct(2));
if hi <= lo, hi = lo + 1; end
imagesc(ax, xc, yc, I, [lo hi]);
colormap(ax, gray);   % grayscale, NOT jet -- these are images, not measurements
set(ax, 'XTick', [], 'YTick', []);
end

function side_axes(ax, xlab, lat, stage, depthMax, tick)
% Depth on y, fixed 0-depthMax at integer `tick` steps.  The x axis is the panel
% coordinate `lat` (microns from the first pixel, which is what the image is
% drawn on) but LABELLED with the stage position `stage` of that same pixel, at
% round `tick`-micron stage values.  Nothing is flipped: the picture is untouched
% and only the numbers change, so the tick values may run either way depending on
% which way the stage axis points.
ylim(ax, [0 depthMax]);
xlim(ax, [lat(1) lat(end)]);
sLo = min(stage);  sHi = max(stage);
sTick = ceil(sLo/tick)*tick : tick : floor(sHi/tick)*tick;
xTick = interp1(stage, lat, sTick);
[xTick, o] = sort(xTick);  sTick = sTick(o);
set(ax, 'YTick', 0:tick:depthMax, 'YTickMode','manual', ...
        'XTick', xTick, 'XTickMode','manual', ...
        'XTickLabel', compose('%.0f', sTick));
daspect(ax, [1 1 1]);   % 1x1x1 um voxels, so the display must be isometric too
xlabel(ax, xlab);  ylabel(ax, 'depth (\mum)');
box(ax,'on');
end

function add_scalebar(ax, W, H, umLen)
x0 = W*0.06;  y0 = H*0.94;
hold(ax,'on');
plot(ax, [x0 x0+umLen], [y0 y0], 'w-', 'LineWidth', 3);
text(ax, x0, y0-H*0.05, sprintf('%g \\mum', umLen), 'Color','w', ...
    'FontSize', 9, 'VerticalAlignment','bottom');
end

% ========================= ScanImage metadata helpers ========================

function meta = read_si_meta(tp)
% Page 1 only.  imfinfo would parse all 23550 IFDs.
t = Tiff(tp, 'r');  cl = onCleanup(@() t.close());
meta = "";
try, meta = string(t.getTag('Software')); catch, end  %#ok<CTCH>
if strlength(meta) == 0
    try, meta = string(t.getTag('ImageDescription')); catch, end  %#ok<CTCH>
end
assert(strlength(meta) > 0, 'No ScanImage metadata in %s', tp);
meta = replace(meta, char(8629), newline);
meta = replace(meta, char(13), newline);
end

function v = si_num(meta, key)
ln = si_line(meta, key);
if strlength(ln) == 0, v = NaN; return; end
tok = regexp(ln, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function vec = si_vec(meta, key)
vec = [];
ln = si_line(meta, key);
if strlength(ln) == 0, return; end
m = regexp(ln, '\[(.*)\]', 'tokens', 'once');
if isempty(m), src = ln; else, src = m{1}; end
nums = regexp(src, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(nums), return; end
vec = str2double(nums(:));
end

function ln = si_line(meta, key)
ln = "";
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
ln = extractBefore(s, newline);
if isempty(ln) || ismissing(ln), ln = s; end
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
