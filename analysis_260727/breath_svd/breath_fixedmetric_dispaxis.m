function breath_fixedmetric_dispaxis(rootDir)
% BREATH_FIXEDMETRIC_DISPAXIS  Diagnostic figure for the rigid-displacement metric.
%
% breath_fixedmetric_analyze.m shows the metric AFTER it exists (traces,
% amplitudes, per-breath scatter).  It never shows how the displacement is
% MEASURED, so the two assumptions the whole metric rests on go unchecked:
%
%   1. the chest wall moves along ONE axis, so collapsing (dx, dy) onto a single
%      frozen direction loses nothing;
%   2. that axis is stable across the session, so a direction frozen on the
%      baseline runs is still the right direction an hour later.
%
% Both are testable from what breath_fixedmetric_extract.py already saved --
% D_all holds the raw per-frame (dx, dy) for every run -- so this is a pure
% re-plot and no video is decoded.
%
% Method being visualised (breath_fixedmetric_extract.py:348):
%   each frame is phase-correlated (cv2.phaseCorrelate, Hanning-windowed)
%   against THAT RUN'S OWN mean frame, giving a sub-pixel (dx, dy); the pooled
%   baseline (dx, dy) are then SVD'd and the top right-singular vector is frozen
%   as the motion axis.  Phase correlation whitens the Fourier magnitude and
%   keeps only the phase ramp, which is why brightness drift does not bias it.
%
% USAGE
%   breath_fixedmetric_dispaxis                        % uses ROOT_DIR below
%   breath_fixedmetric_dispaxis('E:\260724_heteroSHI_nRuns\cam1')
%
% OUTPUT
%   <ROOT_DIR>/breath_fixedmetric_dispaxis.png

%% --------------------------- USER PARAMETERS -------------------------------
ROOT_DIR   = 'E:\260723_heteroSHI_nRuns\cam1';

START_SEC  = 10;      % start of the example trace window
WIN_SEC    = 5;       % length of the example trace window
FMAX       = 10;      % spectrum x-limit (Hz)
MAX_CLOUD  = 20000;   % max points drawn in the (dx,dy) cloud (plot speed only)

% A per-run axis direction only means something if that run's motion actually
% HAS a dominant direction.  When the (dx, dy) cloud is isotropic the top
% singular vector is a random direction and its angle is noise, so runs below
% this variance fraction are drawn but excluded from the drift statistics.
% 0.5 = perfectly isotropic (no axis at all), 1.0 = perfectly 1-D.
AXIS_OK    = 0.70;

% ...and being 1-D is not enough either: a slow postural drift is perfectly 1-D
% while carrying no breathing at all.  A run only tests the frozen axis if the
% motion along its OWN axis is actually breathing, so it must also put this
% fraction of its power in the breathing band.  Measured on the run's own axis,
% never the frozen one, or the test would assume its own conclusion.
BF_OK      = 0.50;
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(rootDir), ROOT_DIR = rootDir; end

% Two levels below the repo root (analysis_260727\breath_svd).
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
addpath(here); addpath(repoRoot);

matFile = fullfile(ROOT_DIR, 'breath_fixedmetric.mat');
assert(isfile(matFile), ...
    'breath_fixedmetric.mat not found in %s -- run breath_fixedmetric_extract.py first.', ROOT_DIR);
S = load(matFile);

fps   = double(S.fps);                 % EFFECTIVE rate, already stride-corrected
D_all = double(S.D_all);               % [T x 2 x nR], columns are (dx, dy)
ax    = double(S.motion_axis(:));      % unit vector in the (dx, dy) plane
ax    = ax / norm(ax);
perp  = [-ax(2); ax(1)];               % orthogonal direction, same units
[T, ~, nR] = size(D_all);
isBase = logical(double(S.is_baseline(:)));
runNames = cellstr(string(S.run_names(:)));

fprintf('breath_fixedmetric_dispaxis: %d runs (%d basis), T = %d frames @ %.3f Hz\n', ...
    nR, sum(isBase), T, fps);
fprintf('  frozen axis (dx,dy) = (%+.3f, %+.3f) = %+.1f deg, %.1f%% of pooled variance\n', ...
    ax(1), ax(2), atan2d(ax(2), ax(1)), 100 * double(S.axis_var));

% Wall-clock run axis, same precedence as breath_fixedmetric_analyze.m: the
% folder-name start time is preferred, because mtime is the END of the run.
tRun = (1:nR)'; tUnit = 'run index';
if isfield(S, 'run_tname') && all(isfinite(double(S.run_tname(:))))
    t0 = double(S.run_tname(:));
    tRun = (t0 - t0(1)) / 3600;  tUnit = 'time from first run (h)';
elseif isfield(S, 'run_hours')
    tRun = double(S.run_hours(:)) - (T / fps) / 3600;
    tRun = tRun - tRun(1);       tUnit = 'time from first run (h)';
end
if ~strcmp(tUnit, 'run index') && max(tRun) * 60 < 120
    tRun = tRun * 60;            tUnit = 'time from first run (min)';
end

%% --- per-run displacement geometry ------------------------------------------
% Each run gets its OWN axis by the same SVD used to freeze the shared one.  A
% run whose axis has swung away from the frozen one is a run where the frozen
% projection is measuring the wrong direction -- which shows up as a spuriously
% small amplitude, not as an obvious failure.
runAng   = nan(nR, 1);   % per-run axis direction, deg, folded to [-90, 90]
runVar   = nan(nR, 1);   % fraction of that run's (dx,dy) variance on its own axis
runAlign = nan(nR, 1);   % |cos| between the run axis and the frozen axis
runSD    = nan(nR, 1);   % SD of the frozen projection, px
runSDo   = nan(nR, 1);   % SD of the orthogonal residual, px
runBF    = nan(nR, 1);   % band fraction along the run's OWN axis
runDev   = nan(nR, 1);   % signed angle from the frozen axis, deg, in [-90, 90]

band = [0.5 8];
if isfield(S, 'band'), band = double(S.band(:))'; end

aFroz = fold90(atan2d(ax(2), ax(1)));

for j = 1:nR
    Dj = D_all(:, :, j);
    Dj = Dj - mean(Dj, 1);                       % per-run centring, as in extract
    [~, sv, V] = svd(Dj, 'econ');
    v = V(:, 1);
    sv = diag(sv);
    runVar(j)   = sv(1)^2 / max(sum(sv.^2), eps);
    runAlign(j) = abs(v' * ax);
    % Sign of a singular vector is arbitrary; fold the angle so the plot does
    % not show a fake 180 deg jump every time SVD flips it.
    runAng(j) = fold90(atan2d(v(2), v(1)));
    % Deviation FROM the frozen axis, folded the same way.  Plotting the raw
    % angle instead puts breathing (~+88 deg here) right on the fold boundary,
    % where a 2 deg wobble reads as a 177 deg swing and the SD is meaningless.
    runDev(j) = fold90(runAng(j) - aFroz);
    runSD(j)  = std(Dj * ax);
    runSDo(j) = std(Dj * perp);
    runBF(j)  = local_bandfrac(Dj * v, fps, band);
end

% A run tests the frozen axis only if it has an axis AND that axis carries
% breathing.  Everything else is excluded from the drift statistics rather
% than being allowed to manufacture a failure.
ok = runVar >= AXIS_OK & runBF >= BF_OK;

% Example run: the basis run with the cleanest breathing to show a trace from.
cand = find(isBase);
if isempty(cand), cand = (1:nR)'; end
[~, k] = max(runBF(cand));
jEx = cand(k);
fprintf('  example run: %s (band fraction %.3f)\n', runNames{jEx}, runBF(jEx));

Dex  = D_all(:, :, jEx) - mean(D_all(:, :, jEx), 1);
pEx  = Dex * ax;                                  % the metric itself
oEx  = Dex * perp;                                % what it discards
tEx  = (0:T-1)' / fps;
i0 = max(1, round(START_SEC * fps) + 1);
i1 = min(T, i0 + round(WIN_SEC * fps) - 1);
w  = i0:i1;

%% --- figure ------------------------------------------------------------------
f = figure('Color', 'w', 'Position', [60 40 1560 980], 'Name', 'displacement axis');
tl = tiledlayout(f, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

cBase = [0.20 0.35 0.85];
cTest = [0.85 0.25 0.20];
cPerp = [0.55 0.55 0.55];

% (1) where the displacement is measured ------------------------------------
nexttile;
if isfield(S, 'ref_mean_img')
    img = double(S.ref_mean_img);
    imagesc(img); colormap(gca, gray); axis image off; hold on;
    roiFile = fullfile(ROOT_DIR, 'breath_roi.mat');
    if isfile(roiFile)
        Rr = load(roiFile);
        if isfield(Rr, 'roi_xywh')
            r = double(Rr.roi_xywh(:))';          % [x y w h], 0-based from cv2
            rectangle('Position', [r(1)+1, r(2)+1, r(3), r(4)], ...
                'EdgeColor', [0 0.85 0], 'LineWidth', 1.5);
            text(r(1)+1, r(2)-6, 'phase-correlation ROI', ...
                'Color', [0 0.6 0], 'FontSize', 8, 'VerticalAlignment', 'bottom');
        end
    end
    % Frozen axis drawn on the image, so the direction is anatomical rather
    % than an abstract pair of numbers.
    cx = size(img, 2) / 2; cy = size(img, 1) / 2;
    L = 0.22 * min(size(img));
    plot(cx + [-1 1] * L * ax(1), cy + [-1 1] * L * ax(2), '-', ...
        'Color', cTest, 'LineWidth', 2);
    title('reference mean frame + frozen axis', 'FontWeight', 'normal');
else
    axis off; text(0.5, 0.5, 'no ref\_mean\_img in .mat', 'HorizontalAlignment', 'center');
end

% (2) the (dx, dy) cloud and the axis it was frozen from ---------------------
nexttile; hold on;
P = [];
for j = find(isBase)'
    Dj = D_all(:, :, j) - mean(D_all(:, :, j), 1);
    P = [P; Dj]; %#ok<AGROW>
end
if isempty(P)
    P = D_all(:, :, jEx) - mean(D_all(:, :, jEx), 1);
end
if size(P, 1) > MAX_CLOUD
    P = P(round(linspace(1, size(P, 1), MAX_CLOUD)), :);
end
% A handful of frames jump ~10 px (a swallow, a postural shift).  They are real
% and stay in every statistic; the VIEW is clipped to a robust limit because
% otherwise breathing -- which is under a pixel -- is a dot at the origin.
Lp = max(prctile(abs(P(:)), 99.5), eps);
% Per-run axes first, so the cloud draws on top of them.  Only runs that HAVE
% an axis are shown: an isotropic run contributes a random line and, at 256
% runs, those alone fill the panel and hide the data.
nOk = 0;
for j = 1:nR
    if ~ok(j), continue; end
    nOk = nOk + 1;
    v = [cosd(runAng(j)); sind(runAng(j))];
    plot([-1 1] * 1.5 * Lp * v(1), [-1 1] * 1.5 * Lp * v(2), '-', ...
        'Color', [0.80 0.80 0.82], 'LineWidth', 0.5);
end
plot(P(:, 1), P(:, 2), '.', 'Color', [0.35 0.40 0.60], 'MarkerSize', 2);
plot([-1 1] * 1.5 * Lp * ax(1),   [-1 1] * 1.5 * Lp * ax(2),   '-', ...
    'Color', cTest, 'LineWidth', 2);
plot([-1 1] * 1.5 * Lp * perp(1), [-1 1] * 1.5 * Lp * perp(2), '--', ...
    'Color', cPerp, 'LineWidth', 1);
axis equal; grid on; box on;
xlim([-Lp Lp]); ylim([-Lp Lp]);
xlabel('dx (px)'); ylabel('dy (px)');
title(sprintf('pooled baseline (dx,dy);  axis = %.1f%% of variance  [view \\pm%.2f px]', ...
    100 * double(S.axis_var), Lp), 'FontWeight', 'normal');

% (3) both raw channels, example window --------------------------------------
nexttile; hold on;
plot(tEx(w), Dex(w, 1), '-', 'Color', cBase, 'LineWidth', 1);
plot(tEx(w), Dex(w, 2), '-', 'Color', cTest, 'LineWidth', 1);
grid on; box on; xlim(tEx([w(1) w(end)]));
xlabel('time (s)'); ylabel('displacement (px)');
legend({'dx', 'dy'}, 'Location', 'best', 'Box', 'off');
title(sprintf('raw phase-correlation output  (%s)', ...
    strrep(runNames{jEx}, '_', '\_')), 'FontWeight', 'normal');

% (4) what the projection keeps vs what it throws away ------------------------
nexttile; hold on;
plot(tEx(w), pEx(w), '-', 'Color', cTest,  'LineWidth', 1.2);
plot(tEx(w), oEx(w), '-', 'Color', cPerp, 'LineWidth', 0.9);
grid on; box on; xlim(tEx([w(1) w(end)]));
xlabel('time (s)'); ylabel('displacement (px)');
legend({'on axis (= the metric)', 'orthogonal (discarded)'}, ...
    'Location', 'best', 'Box', 'off');
title(sprintf('SD %.4f px on axis vs %.4f px orthogonal', ...
    runSD(jEx), runSDo(jEx)), 'FontWeight', 'normal');

% (5) spectra of both, same run ----------------------------------------------
% If the axis is right, the breathing peak lives in the projection and the
% orthogonal residual is broadband -- that is the check, not the SD ratio.
nexttile; hold on;
[fA, PA] = local_psd(pEx, fps);
[~,  PO] = local_psd(oEx, fps);
m = fA <= FMAX;
plot(fA(m), PA(m), '-', 'Color', cTest,  'LineWidth', 1.2);
plot(fA(m), PO(m), '-', 'Color', cPerp, 'LineWidth', 0.9);
set(gca, 'YScale', 'log'); grid on; box on;
yl = ylim;
patch([band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.2 0.6 0.2], 'FaceAlpha', 0.07, 'EdgeColor', 'none');
ylim(yl); xlim([0 FMAX]);
xlabel('frequency (Hz)'); ylabel('power (px^2/Hz)');
legend({'on axis', 'orthogonal'}, 'Location', 'best', 'Box', 'off');
title('spectrum -- breathing should be on-axis only', 'FontWeight', 'normal');

% (6) amplitude scale: is sub-pixel estimation actually required? -------------
nexttile;
mag = [];
for j = 1:nR
    Dj = D_all(:, :, j) - mean(D_all(:, :, j), 1);
    mag = [mag; sqrt(sum(Dj.^2, 2))]; %#ok<AGROW>
end
histogram(mag, 80, 'FaceColor', [0.45 0.5 0.7], 'EdgeColor', 'none');
hold on; grid on; box on;
yl = ylim;
plot([1 1], yl, '-', 'Color', cTest, 'LineWidth', 1.5);
text(1, yl(2), ' 1 px', 'Color', cTest, 'VerticalAlignment', 'top', 'FontSize', 9);
ylim(yl);
xlabel('|displacement| from run mean (px)'); ylabel('frames');
title(sprintf('%.1f%% of frames move < 1 px  (median %.3f px)', ...
    100 * mean(mag < 1), median(mag)), 'FontWeight', 'normal');

% (7) axis stability across the session --------------------------------------
nexttile; hold on;
% Runs with no axis, or an axis carrying no breathing, are shown hollow and
% grey: their angle is the direction of a noise cloud.  Reading their scatter
% as "the axis drifted" is the trap this panel exists to prevent.
plot(tRun(~ok), runDev(~ok), 'o', 'Color', [0.72 0.72 0.72], 'MarkerSize', 3);
plot(tRun(ok & ~isBase), runDev(ok & ~isBase), 'o', 'Color', cTest, ...
    'MarkerFaceColor', cTest, 'MarkerSize', 4);
plot(tRun(ok & isBase),  runDev(ok & isBase),  'o', 'Color', cBase, ...
    'MarkerFaceColor', cBase, 'MarkerSize', 4);
plot(xlim, [0 0], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
grid on; box on; ylim([-95 95]);
xlabel(tUnit); ylabel('deviation from frozen axis (deg)');
title(sprintf('axis drift: %+.1f \\pm %.1f deg over %d breathing runs', ...
    mean(runDev(ok)), std(runDev(ok)), sum(ok)), 'FontWeight', 'normal');

% (8) per-run alignment with the frozen axis ---------------------------------
% |cos| = how much of a run's real motion the frozen projection recovers.  It
% is the amplitude error the metric silently absorbs when the axis has moved.
nexttile; hold on;
plot(tRun(~isBase), runAlign(~isBase), 'o', 'Color', cTest, ...
    'MarkerFaceColor', cTest, 'MarkerSize', 4);
plot(tRun(isBase),  runAlign(isBase),  'o', 'Color', cBase, ...
    'MarkerFaceColor', cBase, 'MarkerSize', 4);
plot(tRun, runVar, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
plot(xlim, [AXIS_OK AXIS_OK], ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
grid on; box on; ylim([0 1.02]);
xlabel(tUnit); ylabel('|cos| with frozen axis');
legend({'test', 'basis', 'run 1-D-ness', 'axis-valid cutoff'}, ...
    'Location', 'south', 'Box', 'off');
% Reported over axis-valid runs only.  A low |cos| on an isotropic run is not
% evidence that the frozen axis is wrong -- there is no right axis to compare
% against -- so mixing those in would manufacture a failure that is not there.
title(sprintf('breathing runs: |cos| \\geq %.2f;  %d of %d runs excluded', ...
    min(runAlign(ok)), sum(~ok), nR), 'FontWeight', 'normal');

% (9) method, in words --------------------------------------------------------
nexttile; axis off;
if isfield(S, 'reg_shift')
    rs = double(S.reg_shift);            % [nR x 2], (dy, dx) integer, from phase 2
    drift = sprintf('%+d..%+d dy, %+d..%+d dx px', ...
        min(rs(:,1)), max(rs(:,1)), min(rs(:,2)), max(rs(:,2)));
else
    drift = 'n/a';
end
strideTxt = 'n/a';
if isfield(S, 'stride'), strideTxt = sprintf('%d', round(double(S.stride))); end
fpsCam = fps;
if isfield(S, 'fps_camera'), fpsCam = double(S.fps_camera); end
txt = {
    '\bfHow the pixel displacement is measured\rm'
    ''
    'Per frame: cv2.phaseCorrelate(run mean, frame, Hanning)'
    '  \rightarrow sub-pixel (dx, dy), phase-only, magnitude whitened'
    '  \rightarrow NOT normxcorr2: no local-variance normalisation'
    '  \rightarrow immune to brightness drift, O(N log N)'
    ''
    'Reference = that run''s OWN mean frame,'
    'so hour-scale FOV drift cannot leak into breathing.'
    ''
    'Axis = top right-singular vector of the pooled'
    'baseline (dx, dy); frozen once, applied to every run.'
    'Polarity fixed once on the pooled projection.'
    ''
    sprintf('camera %.2f fps, stride %s \\rightarrow %.2f fps effective', fpsCam, strideTxt, fps)
    sprintf('between-run FOV drift (integer, separate step): %s', drift)
    sprintf('band of interest: %.2f-%.2f Hz', band(1), band(2))
    };
text(0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'Interpreter', 'tex');

% Interpreter 'none': a Windows path is full of backslashes, and tex reads
% \c, \2 etc. as commands and refuses to render the title at all.
title(tl, sprintf('%s  --  rigid displacement diagnostics  (blue = basis runs, red = test)', ...
    ROOT_DIR), 'FontWeight', 'bold', 'Interpreter', 'none');

if SAVE_FIG
    outPng = fullfile(ROOT_DIR, 'breath_fixedmetric_dispaxis.png');
    try
        exportgraphics(f, outPng, 'Resolution', 200);
    catch
        print(f, outPng, '-dpng', '-r200');
    end
    fprintf('saved %s\n', outPng);
end

fprintf('  %d of %d runs test the axis (1-D >= %.2f AND own-axis band frac >= %.2f)\n', ...
    sum(ok), nR, AXIS_OK, BF_OK);
fprintf('  those runs: deviation %+.1f +/- %.1f deg, min |cos| with frozen = %.3f\n', ...
    mean(runDev(ok)), std(runDev(ok)), min(runAlign(ok)));
fprintf('  those runs: on-axis / orthogonal SD ratio, median %.2f\n', ...
    median(runSD(ok) ./ max(runSDo(ok), eps)));
if any(~ok)
    fprintf('  excluded: %d not 1-D, %d 1-D but not breathing\n', ...
        sum(runVar < AXIS_OK), sum(runVar >= AXIS_OK & runBF < BF_OK));
    fprintf(['  NOTE: on those runs the frozen projection is measuring noise, not\n' ...
             '        breathing -- check them before trusting their amplitudes.\n']);
end

end

% ============================== local functions ==============================

function [f, P] = local_psd(x, fps)
% Hanning-windowed periodogram.  Deliberately plain FFT rather than pwelch or
% Chronux: this panel only needs to show where the power sits, and a local
% function keeps the script free of toolbox dependencies.
x = x(:) - mean(x(:));
n = numel(x);
w = hanning(n);
X = fft(x .* w);
P = abs(X(1:floor(n/2)+1)).^2 / (fps * sum(w.^2));
P(2:end-1) = 2 * P(2:end-1);
f = (0:floor(n/2))' * fps / n;
end

function r = local_bandfrac(x, fps, band)
% Fraction of 0.1-Nyquist power inside `band` -- mirrors band_frac() in
% breath_fixedmetric_extract.py so the run choice here matches its QC table.
[f, P] = local_psd(x, fps);
inb = f >= band(1) & f < band(2);
tot = f >= 0.1 & f < 0.98 * fps / 2;
r = sum(P(inb)) / max(sum(P(tot)), eps);
end

function a = fold90(a)
% Fold an AXIAL angle into [-90, 90].  A singular vector and its negative are
% the same axis, so +170 deg and -10 deg must land on the same number.
a = mod(a + 90, 180) - 90;
end

function w = hanning(n)
w = 0.5 - 0.5 * cos(2*pi*(0:n-1)'/(n-1));
end
