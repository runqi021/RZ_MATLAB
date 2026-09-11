function falloff_260804(dataDir, roiTag, zSurface, powerLut, fitDepth, depthMode, distSpec, channel)
% FALLOFF_260804  Tissue optical penetration from a vessel depth series.
%
% Reads a set of ScanImage stacks of ONE ROI taken down the same column --
% either one power-ramped stack or several overlapping fixed-power stacks --
% and plots signal vs depth, normalised by power, with an exponential fit.
%
%   top 1 % = mean of the brightest PCT_FRAC of pixels in the frame-averaged
%             slice;  bot 1 % = mean of the dimmest PCT_FRAC.  Fractions, not
%             pixel counts, so FOV size and zoom stay comparable.
%
% Every round is treated on its own.  Nothing is pooled or rescaled across
% rounds, so nothing depends on an assumption about why they differ.
%
% Two things that are not optional:
%   - saturated slices are dropped.  Counted on the RAW frames: a pixel
%     clipped in 10 of 30 frames averages below the ceiling and looks clean.
%   - power comes from mW, not %.  ScanImage % is a Pockels setpoint; mW goes
%     as ~%^1.9 near 10-20 %, so signal goes as ~%^3.8 there, not %^2.
%
% USAGE
%   falloff_260804
%   falloff_260804('C:\fall-off\260728_vglut2_vessel', 'roi1')
%
% OUTPUT (into dataDir)
%   falloff_<roiTag>.png / .pdf / .mat
%   falloff_<roiTag>_metrics.mat   cache; delete to force a re-read

%% --------------------------- USER PARAMETERS -------------------------------
DATA_DIR   = 'C:\fall-off\260728_vglut2_vessel';
ROI_TAG    = 'roi1';

POWER_LUT  = 'pre_260723';  % 'pre_260723' | 'post_260723' | acquisition date
P_EXPONENT = 2;             % two-photon: signal ~ power^2

% SI motor z of the PIA.  ScanImage z is inverted -- negative is deeper -- so
% depth below pia = Z_SURFACE - zs, and with the pia at motor 0 that is just
% -zs.  Confirmed by the filenames on this session: zs runs 10 -> -100 for the
% file named z-10-100, -40 -> -160 for z40-160, etc, so the filename z IS depth
% below pia and the stacks START 10 um ABOVE the pia.
% Do NOT default this to the shallowest slice: that put the pia at depth 10 and
% shifted every depth by 10 um.
Z_SURFACE  = 0;
Z_TOP      = Inf;           % ignore slices above this motor z.  Keep Inf: this
                            % is a hard truncation of the raw data, and the
                            % motor zero is not the pia in every session, so a
                            % fixed cut silently deletes real slices.  Points
                            % above the pia are faded in the plots instead.
PCT_FRAC   = 0.05;          % 0.05 = top 5 % and bottom 5 % of pixels.
                            % Each metric is the MEAN over that fraction, and
                            % it is a fraction not a pixel count so different
                            % FOV sizes and zooms stay comparable.

% SI channel number carrying the VESSEL dye, e.g. 3 on 260806 where
% channelSave = [1;3] and channel 1 is the cell label.  [] is fine only when the
% file has a single saved channel.  Channels are INTERLEAVED page by page, so
% reading with the wrong channel index does not error -- it silently measures
% the other label -- which is why this is required rather than guessed.
CHANNEL      = [];

SAT_LEVEL    = [];          % [] = auto from pixel class (int16 -> 32767)
SAT_FRAC_MAX = 1e-3;        % slice is saturated above this clipped fraction

% Depth range (um) used for every exponential fit.  GROUP-DEPENDENT and passed
% in by the caller: 10-200 for wildtype, 10-500 for shiverer.  The lower bound
% of 10 skips the surface layer; the upper bound is where that group still has
% real signal, so the fit is not dragged by the noise floor.
FIT_DEPTH  = [10 200];

% How the depth axis is built.  NOT the same across sessions:
%   'filename' - depth = filename z token + uniform set step.  Use when the
%                filename z IS depth below pia (260728, 260804).
%   'motor'    - depth = Z_SURFACE - (stackZStartPos + k*stackZStepSize).  Use
%                when the filename z is MOTOR z (260114).  Z_SURFACE is then
%                the motor z of the pia and MUST be given.
DEPTH_MODE = 'filename';
MIN_DZ     = 10;            % panel 4: skip slices closer than this to the
                            % reference slice, where (z-z0)/log(y0/y) is 0/0
% Row 1 depths, [first depth, step].  Each curve is ONE PLANE -- the slice
% nearest that depth -- not a slab average.  Sampled at bin CENTRES so a curve
% sits inside a stack rather than on the handover between two of them.
% GROUP-DEPENDENT, passed in by the caller: wildtype gets a 50 um step
% ([25 50] -> 25, 75, 125, ...) because those stacks stop near 200, shiverer a
% 100 um step ([50 100] -> 50, 150, 250, ...) because they run to 600.
DIST_SPEC   = [50 100];
DEPTH_CLIM  = [0 600];      % FIXED colour range for the row-1 depth colourbar.
POWER_CLIM  = [10 350];     % FIXED colour range for the per-stack power colour.
DEPTH_XLIM  = [-50 600];    % FIXED x range for every depth panel (rows 2-4).
DEPTH_XTICK = 0:100:600;    % integer 100 um ticks, identical in every figure.
                            % The -50 margin only keeps the faded above-pia
                            % points visible; the plotted axis is 0-600.
ADU_XLIM    = [-50 32767];  % FIXED x range for the row-1 F distribution, out to
                            % the int16 ADC ceiling
                            % Both fixed on purpose: the same depth, and the
                            % same laser power, get the same colour in every
                            % ROI and every session, so figures can be compared
                            % side by side without reading the legend.
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(dataDir),  DATA_DIR  = dataDir;  end
if nargin >= 2 && ~isempty(roiTag),   ROI_TAG   = roiTag;   end
if nargin >= 3 && ~isempty(zSurface), Z_SURFACE = zSurface; end
if nargin >= 4 && ~isempty(powerLut), POWER_LUT = powerLut; end
if nargin >= 5 && ~isempty(fitDepth), FIT_DEPTH  = fitDepth;  end
if nargin >= 6 && ~isempty(depthMode), DEPTH_MODE = depthMode; end
if nargin >= 7 && ~isempty(distSpec),  DIST_SPEC  = distSpec;  end
if nargin >= 8 && ~isempty(channel),   CHANNEL    = channel;   end

[~, TAG] = match_stacks(DATA_DIR, ROI_TAG);   % filesystem-safe output name
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);              % laser_power_calibration lives here
addpath(fileparts(thisDir));   % repo root, unconditionally

cacheFile = fullfile(DATA_DIR, sprintf('falloff_%s_metrics.mat', TAG));
useCache = false;
if exist(cacheFile, 'file')
    S = load(cacheFile);
    % A cache written from a different channel holds the WRONG label entirely,
    % and nothing downstream could tell.  Re-read instead of trusting it.
    cCh = [];  if isfield(S, 'CHANNEL'), cCh = S.CHANNEL; end
    useCache = isequal(cCh, CHANNEL);
    if ~useCache
        fprintf('cache was made with channel %s, now asking for %s -- re-reading\n', ...
            mat2str(cCh), mat2str(CHANNEL));
    end
end
if useCache
    fprintf('using cached metrics: %s  (delete to re-read)\n', cacheFile);
    R = S.R;
else
    R = read_metrics(DATA_DIR, ROI_TAG, Z_TOP, PCT_FRAC, SAT_LEVEL, SAT_FRAC_MAX, CHANNEL);
    save(cacheFile, 'R', 'Z_TOP', 'PCT_FRAC', 'DATA_DIR', 'ROI_TAG', 'CHANNEL');
    fprintf('cached -> %s\n', cacheFile);
end

%% --- depth axis --------------------------------------------------------------
% DEPTH COMES FROM THE FILENAME, not from the motor.  The z token is the depth
% below pia that was recorded at acquisition, and it spans the stack linearly.
% ScanImage motor z is NOT a common frame across stacks in a session -- the
% zero moves between acquisitions -- so motor z cannot be used to put stacks on
% one depth axis.  The filename can: on both sessions checked, every stack's
% slice count matches its filename range at ~1 um steps.
% DEPTH = the range and step size that were SET, applied uniformly.  The start
% comes from the filename z token (depth below pia at acquisition) and the step
% from SI.hStackManager.stackZStepSize.  Per-slice zs readings are not used:
% they come back with non-uniform steps (0.995 um on some stacks) and their
% motor zero moves between acquisitions, so they cannot put stacks on one axis.
switch lower(DEPTH_MODE)
case 'filename'
    fprintf('\ndepth below pia = filename z start + uniform set step size\n');
    fprintf('  %-12s %8s %8s %20s %14s\n', 'stack', 'nSlice', 'step', 'depth below pia', 'filename z');
    for i = 1:numel(R)
        t = regexp(R(i).name, '_z(-?\d+)-(\d+)_', 'tokens', 'once');
        assert(~isempty(t), ['No _z<start>-<end>_ token in %s -- if the ' ...
            'filename z is MOTOR z, use DEPTH_MODE = ''motor'''], R(i).name);
        zA = str2double(t{1});  zB = str2double(t{2});
        n  = numel(R(i).z);
        dz = abs(R(i).stepSet);
        if ~isfinite(dz) || dz == 0, dz = (zB - zA)/(n-1); end
        R(i).depth = zA + (0:n-1).' * dz;
        fprintf('  %-12s %8d %8.3f %8.1f -> %-10.1f %6.0f -> %-6.0f\n', R(i).label, n, dz, ...
            R(i).depth(1), R(i).depth(end), zA, zB);
    end
case 'motor'
    assert(~isempty(Z_SURFACE) && isfinite(Z_SURFACE), ...
        'DEPTH_MODE ''motor'' needs Z_SURFACE = motor z of the pia');
    fprintf(['\ndepth below pia = %g - motor z, from the SET start and step\n' ...
             '  (pia at motor z = %g)\n'], Z_SURFACE, Z_SURFACE);
    fprintf('  %-12s %8s %8s %18s %20s\n', 'stack', 'nSlice', 'step', 'motor z', 'depth below pia');
    for i = 1:numel(R)
        n  = numel(R(i).z);
        dz = R(i).stepSet;                        % SIGNED: negative goes deeper
        z0 = R(i).zStartSet;
        if ~isfinite(dz) || dz == 0, dz = median(diff(R(i).z)); end
        if ~isfinite(z0), z0 = R(i).z(1); end
        zSet = z0 + (0:n-1).' * dz;
        R(i).depth = Z_SURFACE - zSet;
        fprintf('  %-12s %8d %8.2f %7.0f -> %-8.0f %8.1f -> %-10.1f\n', R(i).label, ...
            n, dz, zSet(1), zSet(end), R(i).depth(1), R(i).depth(end));
    end
otherwise
    error('DEPTH_MODE must be ''filename'' or ''motor'', got ''%s''', DEPTH_MODE);
end

%% --- power -------------------------------------------------------------------
fprintf('\npower calibration: %s\n', POWER_LUT);
for i = 1:numel(R)
    [R(i).mW, cal] = laser_power_calibration(R(i).pctZ, POWER_LUT);
    R(i).lutName = cal.name;
    R(i).clamped = any(cal.clamped);
end
fprintf('  %-12s %10s %16s %8s\n', 'round', '%', 'mW', 'nSlice');
for i = 1:numel(R)
    if R(i).ramped
        sp = sprintf('%g-%g', R(i).pctZ(1), R(i).pctZ(end));
        sm = sprintf('%.1f-%.1f', R(i).mW(1), R(i).mW(end));
    else
        sp = sprintf('%g', R(i).pct);
        sm = sprintf('%.1f', R(i).mW(1));
    end
    fprintf('  %-12s %10s %16s %8d%s\n', R(i).label, sp, sm, ...
        numel(R(i).z), tern(R(i).clamped, '   [CLAMPED]', ''));
end

%% --- normalise ---------------------------------------------------------------
for i = 1:numel(R)
    pw = R(i).mW(:) .^ P_EXPONENT;
    % top - bot removes the additive floor: bot 1 % is the detector offset plus
    % whatever background sits under the whole frame, and it does not scale
    % with P^2, so subtracting it before dividing by power is what makes the
    % deep end mean anything.
    R(i).diff   = R(i).top - R(i).bot;
    R(i).diff_n = R(i).diff ./ pw;
    R(i).fitMask = R(i).depth >= FIT_DEPTH(1) & ...
                   R(i).depth <= FIT_DEPTH(2) & R(i).diff_n > 0;
end

%% --- overlap ownership --------------------------------------------------------
% Consecutive stacks overlap in depth.  Each overlap is split at its MIDPOINT:
% the shallower stack keeps the top half, the deeper stack takes the bottom half.
% Every depth then belongs to exactly one stack.  Without this the POOLED fit
% counts overlap depths twice, once from each stack, and the two disagree, so
% the fit is pulled toward wherever the overlaps happen to fall.
% Per-stack fits still use ALL of that stack's own slices -- the split only
% decides who contributes to the pooled fit.
[~, ord] = sort(arrayfun(@(r) min(r.depth), R));
R = R(ord);
lo = -inf(numel(R),1);  hi = inf(numel(R),1);
for i = 1:numel(R)-1
    ovA = max(min(R(i).depth), min(R(i+1).depth));
    ovB = min(max(R(i).depth), max(R(i+1).depth));
    if ovB > ovA
        mid = (ovA + ovB) / 2;
        hi(i) = mid;  lo(i+1) = mid;
    else
        hi(i) = (max(R(i).depth) + min(R(i+1).depth)) / 2;
        lo(i+1) = hi(i);
    end
end
fprintf('\noverlap split at the midpoint -- who owns which depths (pooled fit)\n');
fprintf('  %-10s %18s %18s %8s\n', 'stack', 'own range', 'owns for pooling', 'n');
for i = 1:numel(R)
    R(i).poolMask = R(i).fitMask & R(i).depth > lo(i) & R(i).depth <= hi(i);
    fprintf('  %-10s %7.0f -> %-8.0f %7.0f -> %-8.0f %8d\n', R(i).label, ...
        min(R(i).depth), max(R(i).depth), ...
        max(lo(i), min(R(i).depth)), min(hi(i), max(R(i).depth)), sum(R(i).poolMask));
end

%% --- fits, per stack and pooled, for each metric ------------------------------
POOL = struct();
for fld = {'diff_n'}
    f = fld{1};
    fprintf('\n%s: exponential fit (depth >= %g um)\n', ...
        sprintf('(TOP %g%% - BOT %g%%) / mW^2', 100*PCT_FRAC, 100*PCT_FRAC), FIT_DEPTH(1));
    fprintf('  %-10s %6s %14s %10s %10s\n', 'round', 'n', 'depth range', 'l (um)', 'R^2');
    dP2 = []; yP2 = [];
    for i = 1:numel(R)
        m = R(i).fitMask;
        [L, A, rr] = fit_exp(R(i).depth(m), R(i).(f)(m));
        R(i).([f '_L']) = L;  R(i).([f '_A']) = A;
        fprintf('  %-10s %6d %6.0f-%-7.0f %10.1f %10.4f\n', R(i).label, sum(m), ...
            min(R(i).depth(m)), max(R(i).depth(m)), L, rr);
        mp  = R(i).poolMask;   % overlap-owned subset only, no double counting
        dP2 = [dP2; R(i).depth(mp)];  yP2 = [yP2; R(i).(f)(mp)];  %#ok<AGROW>
    end
    [dP2, o] = sort(dP2);  yP2 = yP2(o);
    [L, A, rr] = fit_exp(dP2, yP2);
    POOL.(f) = struct('d', dP2, 'L', L, 'A', A, 'r2', rr);
    fprintf('  %-10s %6d %6.0f-%-7.0f %10.1f %10.4f   <- POOLED\n', ...
        'all', numel(dP2), min(dP2), max(dP2), L, rr);
end

%% --- per-slice attenuation factor --------------------------------------------
% factor(z) = y(z0) / y(z), one value per z slice, straight from the tail.
% Dimensionless: how many times dimmer than the reference slice.
% F(0) is the value at depth 0 -- literally at the pia -- and it is the single
% reference for every stack.  Only the shallowest stack reaches z = 0.
iRef = []; kRef = [];
for i = 1:numel(R)
    k = find(abs(R(i).depth) < 0.5, 1);
    if ~isempty(k), iRef = i; kRef = k; break; end
end
assert(~isempty(iRef), ...
    'No slice at depth 0 -- F(0) is not measured, so F(0)/F(z) cannot be formed.');
dRefAll = R(iRef).depth(kRef);
F0raw   = R(iRef).diff(kRef);
F0n     = R(iRef).diff_n(kRef);

fprintf('\nF(0) taken from %s at depth %.0f um (SI z = %+.0f)\n', ...
    R(iRef).label, dRefAll, R(iRef).z(kRef));
fprintf('  raw (top-bot) = %.1f ADU     /mW^2 = %.4g\n', F0raw, F0n);

fprintf('\nattenuation F(0)/F(z), single reference at z = 0\n');
fprintf('  %-10s %12s %20s %20s\n', 'round', 'n', 'raw factor', 'power-normalised');
for i = 1:numel(R)
    good = R(i).diff > 0;
    R(i).att  = F0raw ./ R(i).diff;      R(i).att(~good)  = NaN;
    R(i).attP = F0n   ./ R(i).diff_n;    R(i).attP(~good) = NaN;
    a = R(i).att(isfinite(R(i).att));  b = R(i).attP(isfinite(R(i).attP));
    fprintf('  %-10s %12d %8.2f - %-8.2f %9.2f - %-8.2f\n', R(i).label, numel(a), ...
        min(a), max(a), min(b), max(b));
end

% Pooled fit across all stacks.  y_ref/y rises as exp(+z/l), so the slope of
% log(attP) vs depth is +1/l.  This is the only fit that spans the full depth
% range; the per-round fits in panel 2 each cover ~80 um.
dP = []; aP = [];
for i = 1:numel(R)
    m  = R(i).poolMask & isfinite(R(i).attP) & R(i).attP > 0;
    dP = [dP; R(i).depth(m)];   %#ok<AGROW>
    aP = [aP; R(i).attP(m)];    %#ok<AGROW>
end
[dP, o] = sort(dP);  aP = aP(o);
pP  = polyfit(dP, log(aP), 1);
Lp  = 1/pP(1);
Ap  = exp(pP(2));
res = log(aP) - polyval(pP, dP);
r2p = 1 - sum(res.^2) / sum((log(aP) - mean(log(aP))).^2);
POOLatt = struct('d', dP, 'A', Ap, 'L', Lp, 'r2', r2p);
fprintf('  pooled exp fit over %.0f-%.0f um (%d points): l = %.1f um, R^2 = %.4f\n', ...
    min(dP), max(dP), numel(dP), Lp, r2p);

%% --- figures ------------------------------------------------------------------
% One jet ramp used for BOTH the depth colourbar (row 1) and the per-stack power
% colour (rows 2-4), each on its own fixed range: blue = shallow / low power,
% red = deep / high power.
nC   = 256;
cmap = jet(nC);
ramp = @(v, lim) cmap(min(nC, max(1, round(1 + (nC-1) * ...
    (v - lim(1)) / diff(lim)))), :);
depth2col = @(d) ramp(d, DEPTH_CLIM);
% A ramped stack has no single power; colour it by its MEDIAN mW, and let the
% legend carry the full range.
cols = cell2mat(arrayfun(@(r) ramp(median(r.mW), POWER_CLIM), R(:), 'uni', 0));

for i = 1:numel(R)
    if R(i).ramped
        R(i).labelmW = sprintf('%.0f-%.0f mW', R(i).mW(1), R(i).mW(end));
    else
        R(i).labelmW = sprintf('%.0f mW', R(i).mW(1));
    end
end

%   field     y label                        subtitle                    yscale  fade  fit   legend
spec = { ...
 'top',   sprintf('top %g%% (ADU)',100*PCT_FRAC),           sprintf('raw TOP %g%% of pixels -- faded = above pia', 100*PCT_FRAC), ...
                                                                     'linear', true,  [],                                    'stacks'
 'bot',   sprintf('bot %g%% (ADU)',100*PCT_FRAC),           sprintf('raw BOT %g%% of pixels', 100*PCT_FRAC), ...
                                                                     'linear', true,  [],                                    'none'
 'diff',  sprintf('top %g%% - bot %g%% (ADU)',100*PCT_FRAC,100*PCT_FRAC),  'raw TOP - BOT -- floor removed', 'linear', true,  [],                                    'none'
 'diff_n',sprintf('(top - bot) / mW^%g',P_EXPONENT), '(TOP - BOT) normalised by power', ...
                                                                     'log',    true,  @(ax) draw_fit(ax, POOL.diff_n, -1, true),  'fit'
 'att',   'attenuation  F(0) / F(z)','ATTENUATION (top - bot), raw, single reference F(0) at z = 0', ...
                                                                     'linear', true,  [],                                    'none'
 'attP',  'attenuation  F(0) / F(z)','ATTENUATION (top - bot), power-normalised, F(0) at z = 0', ...
                                                                     'log',    true,  @(ax) draw_fit(ax, POOLatt, +1, false), 'none' };

f = figure('Color','w','Position',[40 10 1400 1550],'Name','falloff');
tl = tiledlayout(f, 4, 2, 'TileSpacing','compact','Padding','compact');

% --- row 1: the F distribution itself, every DIST_STEP um --------------------
% Each depth is taken from the stack whose own range is centred nearest it, so
% every line is one real slice.  Raw ADU -- the power differs between stacks,
% which is why the legend carries it.
maxD  = max(arrayfun(@(r) max(r.depth), R));
dSpec = DIST_SPEC;
dWant = dSpec(1) : dSpec(2) : maxD;
ax = nexttile([1 2]); hold(ax,'on');
ctr = arrayfun(@(r) mean([min(r.depth) max(r.depth)]), R);
xhi = 0;
for j = 1:numel(dWant)
    ok = find(dWant(j) >= arrayfun(@(r) min(r.depth), R) & ...
              dWant(j) <= arrayfun(@(r) max(r.depth), R));
    if isempty(ok), continue; end
    [~, w] = min(abs(dWant(j) - ctr(ok)));  i = ok(w);
    [~, k] = min(abs(R(i).depth - dWant(j)));
    h = double(R(i).hist(k,:));
    m = h > 0;
    h(~m) = NaN;   % NaN not 0, so the line breaks at empty bins instead of
                   % dropping to the axis floor on a log scale
    plot(ax, R(i).hCenters, h, '-', 'Color', depth2col(dWant(j)), ...
        'LineWidth', 1.1, 'HandleVisibility', 'off');
    xhi = max(xhi, max(R(i).hCenters(m)));
end
set(ax, 'YScale', 'log');
xlim(ax, ADU_XLIM);
grid(ax,'on'); box(ax,'on');
colormap(ax, cmap);  caxis(ax, DEPTH_CLIM);
cb = colorbar(ax);  cb.Label.String = 'depth below pia (\mum)';
cb.Direction = 'reverse';   % 0 at the top, deep at the bottom, like the tissue
xlabel(ax, 'pixel value F (ADU)');  ylabel(ax, 'count');
title(ax, sprintf('F distribution at %g:%g:%g \\mum', ...
    dSpec(1), dSpec(2), max(dWant)), 'FontWeight','normal');
for q = 1:size(spec,1)
    [fld, ylab, sub, ys, fade, fitFcn, legMode] = spec{q,:};
    ax = nexttile; hold(ax,'on');
    for i = 1:numel(R)
        y = R(i).(fld);
        % Above the pia (depth < 0) is shown but never analysed -- faded, and
        % kept out of the legend.  Nothing else is faded.
        up = R(i).depth < 0;
        if fade && any(up)
            plot(ax, R(i).depth(up), y(up), '.', 'Color', (cols(i,:)+2)/3, ...
                'MarkerSize', 6, 'HandleVisibility','off');
        end
        plot(ax, R(i).depth(~up), y(~up), '.', 'Color', cols(i,:), ...
            'MarkerSize', 9, 'DisplayName', R(i).labelmW, ...
            'HandleVisibility', tern(strcmp(legMode,'stacks'), 'on', 'off'));
    end
    if ~isempty(fitFcn), fitFcn(ax); end
    set(ax, 'YScale', ys);
    xlim(ax, DEPTH_XLIM);  set(ax, 'XTick', DEPTH_XTICK);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax, 'depth (\mum)');  ylabel(ax, ylab);
    title(ax, sub, 'FontWeight','normal');
    if ~strcmp(legMode, 'none')
        legend(ax, 'Location', tern(q==1,'northeast','southwest'), 'Box','off');
    end
    if strcmp(legMode, 'stacks')
        % Legend names the discrete stack powers; the colourbar puts them on the
        % same fixed 10-350 mW scale used in every other figure.
        colormap(ax, cmap);  caxis(ax, POWER_CLIM);
        cbp = colorbar(ax);  cbp.Label.String = 'laser power (mW)';
    end
end
title(tl, sprintf('%s  --  %s  --  optical falloff  (LUT %s, signal ~ mW^%g)', ...
    DATA_DIR, ROI_TAG, R(1).lutName, P_EXPONENT), ...
    'FontWeight','bold','Interpreter','none');

if SAVE_FIG
    base = fullfile(DATA_DIR, sprintf('falloff_%s', TAG));
    exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'R', 'POWER_LUT', 'P_EXPONENT', 'PCT_FRAC', ...
        'FIT_DEPTH', 'MIN_DZ', 'Lp', 'Ap', 'r2p', 'dRefAll');
    fprintf('\nsaved %s.png/.pdf/.mat\n', base);
end
end

function draw_fit(ax, P, sgn, showLeg)
% One pooled exponential fit.  sgn = -1 for a decaying quantity, +1 for the
% attenuation factor, which is its reciprocal and therefore rises.  Both are
% the SAME fit, so the constant is quoted once and the other panel just draws
% the line.
dF = linspace(min(P.d), max(P.d), 200);
h = plot(ax, dF, P.A*exp(sgn*dF/P.L), 'k--', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('linear fit:  l = %.1f \\mum,  R^2 = %.3f', P.L, P.r2));
if ~showLeg, set(h, 'HandleVisibility', 'off'); end
end


% ============================== read =========================================

function R = read_metrics(DATA_DIR, ROI_TAG, Z_TOP, PCT_FRAC, SAT_LEVEL, SAT_FRAC_MAX, CHANNEL)
files = match_stacks(DATA_DIR, ROI_TAG);
assert(~isempty(files), 'No stacks matching ''%s'' under %s', ROI_TAG, DATA_DIR);

R = struct('name',{},'label',{},'pct',{},'pctZ',{},'z',{},'top',{},'bot',{}, ...
           'satFrac',{},'sat',{},'satLevel',{},'ramped',{});
P0all = nan(numel(files),1);
M = cell(numel(files),1);
for i = 1:numel(files)
    info = imfinfo(fullfile(files(i).folder, files(i).name));
    M{i} = si_meta(info);
    P0all(i) = si_num(M{i}, "SI.hBeams.stackStartPower");
    M{i} = {M{i}, info};
end
[~, ord] = sort(P0all);
files = files(ord);  M = M(ord);

for i = 1:numel(files)
    meta = M{i}{1};  info = M{i}{2};
    tp = fullfile(files(i).folder, files(i).name);

    zs   = si_vec(meta, "SI.hStackManager.zs");
    nSl  = si_num(meta, "SI.hStackManager.numSlices");
    fpsl = si_num(meta, "SI.hStackManager.framesPerSlice");
    P0   = si_num(meta, "SI.hBeams.stackStartPower");
    P1   = si_num(meta, "SI.hBeams.stackEndPower");
    LC   = si_num(meta, "SI.hBeams.lengthConstants");
    pzOn = si_flag(meta, "SI.hBeams.pzAdjust");
    chSave = si_vec(meta, "SI.hChannels.channelSave");
    nCh    = max(numel(chSave), 1);
    chIdx  = pick_channel(chSave, CHANNEL, files(i).name);

    if isempty(zs)
        z0 = si_num(meta, "SI.hStackManager.stackZStartPos");
        z1 = si_num(meta, "SI.hStackManager.stackZEndPos");
        dz = si_num(meta, "SI.hStackManager.stackZStepSize");
        assert(all(isfinite([z0 z1 dz])) && dz ~= 0, 'No usable z in %s', files(i).name);
        zs = (z0:dz:z1).';
    end
    assert(numel(zs) == nSl, '%s: zs/numSlices mismatch', files(i).name);
    assert(numel(info) == nSl*fpsl*nCh, '%s: page count mismatch', files(i).name);

    % Ramp only if ScanImage says so.  lengthConstants is present even when
    % pzAdjust is off, and applying it then invents a rise across a flat stack.
    ramped = pzOn && ~(isfinite(P0) && isfinite(P1) && abs(P1-P0) < 1e-9);
    if ramped
        assert(isfinite(LC) && LC ~= 0, '%s: ramp on but lengthConstants unusable', files(i).name);
        zRef = si_num(meta, "SI.hStackManager.zPowerReference");
        if ~isfinite(zRef), zRef = zs(1); end
        pctZ = P0 .* exp((zs - zRef) ./ LC);
    else
        pctZ = repmat(P0, numel(zs), 1);
    end

    H = info(1).Height;  W = info(1).Width;
    idxZ = find(zs <= Z_TOP);  nK = numel(idxZ);
    nPct = max(1, round(PCT_FRAC * H * W));

    fprintf('\n[%d/%d] %s  %g%%%s  %d/%d slices x %d frames%s\n', i, numel(files), ...
        files(i).name, P0, tern(ramped,' RAMPED',''), nK, nSl, fpsl, ...
        tern(nCh > 1, sprintf('  [ch %g of %s]', chSave(chIdx), mat2str(chSave(:).')), ''));

    pTop = nan(nK,1);  pBot = nan(nK,1);  satFrac = nan(nK,1);
    % Per-slice pixel-value histogram on FIXED edges shared by every slice and
    % every stack, so distributions at different depths are directly
    % comparable.  ~13 MB for a whole ROI, so it is cached and any depth can be
    % replotted without touching the TIFFs again.
    hEdges = -1024 : 8 : 32768;
    Hst = zeros(nK, numel(hEdges)-1, 'single');
    t = Tiff(tp,'r');  cl = onCleanup(@() t.close());
    if isempty(SAT_LEVEL)
        setDirectory(t,1);  satLev = double(intmax(class(t.read())));
    else
        satLev = SAT_LEVEL;
    end

    tS = tic;
    for k = 1:nK
        acc = zeros(H,W);  nSat = 0;
        for fr = 1:fpsl
            % Channels are interleaved page by page within each frame.
            setDirectory(t, ((idxZ(k)-1)*fpsl + fr-1)*nCh + chIdx);
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

        if mod(k, max(1,floor(nK/5))) == 0 || k == nK
            fprintf('    z %6.1f (%3d/%3d)  top1%% %9.1f  bot1%% %8.1f  sat %.3f%%\n', ...
                zs(idxZ(k)), k, nK, pTop(k), pBot(k), 100*satFrac(k));
        end
    end
    clear cl
    fprintf('    %.1f s\n', toc(tS));

    n = numel(R)+1;
    R(n).name = files(i).name;
    if ramped
        R(n).label = sprintf('%g-%g%% lp', P0, P1);
    else
        R(n).label = sprintf('%g%% lp', P0);
    end
    R(n).pct  = P0;
    R(n).pctZ = pctZ(idxZ);
    R(n).z    = zs(idxZ);
    R(n).top = pTop;
    R(n).bot = pBot;
    R(n).satFrac = satFrac;
    R(n).sat  = satFrac > SAT_FRAC_MAX;
    R(n).satLevel = satLev;
    R(n).ramped = ramped;
    R(n).stepSet   = si_num(meta, "SI.hStackManager.stackZStepSize");
    R(n).zStartSet = si_num(meta, "SI.hStackManager.stackZStartPos");
    R(n).hist    = Hst;
    R(n).hCenters = (hEdges(1:end-1) + hEdges(2:end)) / 2;
    R(n).nCh     = nCh;
    R(n).channel = chSave(chIdx);
end
end

% --------------------------------------------------------------------------
function chIdx = pick_channel(chSave, CHANNEL, name)
% Position within the interleaved pages of the requested SI channel.
% Single-channel files need no choice; multi-channel ones must be told, because
% every channel reads without error and only the values differ.
nCh = max(numel(chSave), 1);
if nCh <= 1
    chIdx = 1;
    if ~isempty(CHANNEL) && ~isempty(chSave) && chSave(1) ~= CHANNEL
        warning('%s saves only channel %g; CHANNEL = %g ignored', name, chSave(1), CHANNEL);
    end
    return
end
assert(~isempty(CHANNEL), ['%s saves %d channels (%s) -- set CHANNEL to the ' ...
    'one carrying the vessel dye'], name, nCh, mat2str(chSave(:).'));
chIdx = find(chSave == CHANNEL, 1);
assert(~isempty(chIdx), '%s: channel %g not saved (saved: %s)', ...
    name, CHANNEL, mat2str(chSave(:).'));
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

function [dc, L, ci, nPts] = local_length(d, y, win)
% One fit per NON-OVERLAPPING depth bin, each on its own raw points, with the
% 95 % CI on the slope propagated to L.  Sliding windows share most of their
% data, so they look smoother and better-determined than the evidence supports.
dc = []; L = []; ci = []; nPts = [];
if numel(d) < 10, return; end
e = (min(d) : win : max(d) + win).';
for i = 1:numel(e)-1
    m = d >= e(i) & d < e(i+1) & y > 0;
    if nnz(m) < 8, continue; end
    [b, bint] = regress(log(y(m)), [ones(nnz(m),1), d(m)]);
    if b(2) >= 0, continue; end
    dc(end+1,1)   = mean(d(m));            %#ok<AGROW>
    L(end+1,1)    = -1/b(2);               %#ok<AGROW>
    nPts(end+1,1) = nnz(m);                %#ok<AGROW>
    hi = -1/min(bint(2,:));  lo = -1/max(bint(2,:));
    ci(end+1,:) = [min(lo,hi), max(lo,hi)]; %#ok<AGROW>
end
end

% ============================== plotting =====================================

function draw_split(ax, x, y, bad, col, lbl)
if any(bad)
    ys = y; ys(~bad) = NaN;
    yc = y; yc(bad)  = NaN;
    plot(ax, x, ys, ':', 'Color', col, 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(ax, x, yc, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', lbl);
else
    plot(ax, x, y, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', lbl);
end
end

function finish(ax, xl, yl, ttl)
grid(ax,'on'); box(ax,'on');
xlabel(ax, xl); ylabel(ax, yl);
title(ax, ttl, 'FontWeight','normal');
end

% ========================= ScanImage metadata helpers ========================

function meta = si_meta(info)
meta = "";
if isfield(info(1),'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1),'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
end
assert(strlength(meta) > 0, 'No ScanImage metadata in Software/ImageDescription');
meta = replace(meta, char(8629), newline);
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

function tf = si_flag(meta, key)
tf = startsWith(lower(strtrim(si_line(meta, key))), "true");
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

function [files, tag] = match_stacks(DATA_DIR, ROI_TAG)
% ROI_TAG may be a plain prefix ('roi1', 'tile15') OR a glob when the part that
% distinguishes FOVs sits in the MIDDLE of the name ('*col00_row02*'), as on
% 251104 where every file starts with the same experiment string.
if contains(ROI_TAG, '*') || endsWith(lower(ROI_TAG), '.tif')
    pat = ROI_TAG;
    if ~endsWith(lower(pat), '.tif'), pat = [pat '.tif']; end
else
    pat = [ROI_TAG '*.tif'];
end
files = dir(fullfile(DATA_DIR, pat));
tag = regexprep(ROI_TAG, '[\*\.\s/\\]+', '_');
tag = regexprep(tag, '^_|_$', '');
end
