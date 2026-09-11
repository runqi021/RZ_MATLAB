function breath_fixedmetric_analyze(rootDir)
% breath_fixedmetric_analyze  Baseline-referenced breathing metrics + QC.
%
% Consumes breath_fixedmetric.mat written by breath_fixedmetric_extract.py and
% answers the question that motivated the whole thing: is the readout
% comparable across runs (and, by extension, across animals)?
%
% Three metrics are carried side by side:
%   disp -- rigid chest-wall translation in PIXELS on a frozen axis (primary)
%   fb   -- frozen-basis SVD, one spatial mode shared by every run
%   pv   -- per-video SVD PC1, the incumbent, for benchmarking only
%
% Amplitude is expressed as a multiple of the pooled-baseline median breath
% peak-to-trough, NOT as a z-score.  Z-scoring divides by the SD of the whole
% trace, so a run containing one large body movement gets an inflated SD and
% its breathing silently shrinks -- exactly backwards for comparability.  A
% physiological reference does not have that failure mode, and "1.0" keeps a
% concrete meaning (one baseline breath) across animals.
%
% The run-to-run scatter of baseline amplitude is reported explicitly: with
% several baseline runs, that scatter IS the measurement noise floor, and it is
% the empirical test of the assumption that baseline breathing is a stable
% standard.  If it drifts monotonically, the assumption is false in your prep
% and you want to know before it reaches the group data.
%
% USAGE
%   breath_fixedmetric_analyze                     % uses ROOT_DIR below
%   breath_fixedmetric_analyze('D:\breath_tracking_test')

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR    = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
BAND        = [0.1 20];    % Hz, bandpass applied before breath detection
MIN_SEP_SEC = 0.12;       % refractory between breaths (s); 0.12 s -> <8.3 Hz
MIN_PROM_F  = 0.5;        % peak prominence, as a fraction of baseline median p2p
NDROP       = 0;         % drop leading frames (camera/illumination settling)
USE_CURATED = true;       % prefer hand-curated peaks/troughs from the GUIs

% ---- run selection: both decided here, not in the extractor ----------------
% Changing either costs seconds, because neither touches the video or the cache.
%
% BASELINE_RUNS defines the amplitude reference -- the runs whose median breath
% peak-to-trough becomes "1.0".  Make it a reference CONDITION, not "everything
% before the manipulation": a block spanning many hours of drifting prep gives a
% yardstick that describes none of the runs it is applied to.
%   [] = use whatever the extractor flagged (its BASELINE_RUNS)
%
% Which runs are INCLUDED is decided by the detection, not by a hand-set cutoff:
% a run needs at least MIN_BREATHS detected events to enter any statistic, plot
% or reference.  Once the peak/trough GUIs have run that is your own curation
% deciding -- a run where the animal has stopped breathing yields no events and
% drops out by itself, with nothing inferred from signal shape.
BASELINE_RUNS  = 1:16;    % e.g. 1:20
MIN_BREATHS    = 1;       % fewer detected breaths than this = run excluded

SAVE_FIGS   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(rootDir), ROOT_DIR = rootDir; end

% Two levels below the repo root (analysis_260727\breath_svd) -- see the note
% in the analysis_260727 README about why one fileparts is not enough.
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
addpath(here); addpath(repoRoot);

matFile = fullfile(ROOT_DIR, 'breath_fixedmetric.mat');
assert(isfile(matFile), ...
    'breath_fixedmetric.mat not found in %s -- run breath_fixedmetric_extract.py first.', ROOT_DIR);
S = load(matFile);

fps      = double(S.fps);
runNames = cellstr(string(S.run_names(:)));
runNum   = double(S.run_num(:));
nR       = numel(runNames);
metrics  = {'disp', 'fb', 'pv'};
labels   = {'rigid displacement (px)', 'frozen-basis PC1', 'per-video PC1'};


% Wall-clock axis, anchored to each run's START.
%
% This matters for the per-breath scatter, which adds a peak's offset WITHIN the
% run to the run's time: added to an end-time it would place every breath one
% run-duration late.  run_hours is derived from the video mtime, which on this
% rig is the END of the recording (start + 60 s), so the folder-name timestamp
% is preferred and mtime is corrected by the run duration when it is not.
tRun = (1:nR)'; tUnit = 'run index'; haveTime = false;
if isfield(S, 'run_tname') && all(isfinite(double(S.run_tname(:))))
    t0 = double(S.run_tname(:));                    % folder-name start, epoch s
    tRun = (t0 - t0(1)) / 3600;
    haveTime = true;
elseif isfield(S, 'run_hours')
    durH = (double(S.T) / fps) / 3600;              % one run, in hours
    tRun = double(S.run_hours(:)) - durH;
    tRun = tRun - tRun(1);
    haveTime = true;
end
if haveTime
    if max(tRun) * 60 < 120
        tRun = tRun * 60; tUnit = 'time from first run (min)';
    else
        tUnit = 'time from first run (h)';
    end
end

% Baseline/valid counts are reported after detection, since the run selection
% now depends on it.
fprintf('breath_fixedmetric_analyze: %d runs, fps = %.4f\n', nR, fps);
if haveTime
    spanH = (tRun(end) - tRun(1)) / tern(contains(tUnit, '(min)'), 60, 1);
    fprintf('  session spans %.2f h wall clock; x-axis in %s\n', spanH, tUnit);
end
fprintf('  motion axis (dx,dy) = (%+.3f, %+.3f), %.1f%% of displacement variance\n', ...
    S.motion_axis(1), S.motion_axis(2), 100 * S.axis_var);

%% --- pass 1: bandpass + detect ----------------------------------------------
% Detection has to finish before run selection, because a run is included only
% if the detection found breaths in it -- see the selection block below.
R = struct();
for m = 1:numel(metrics)
    X = double(S.(metrics{m}));
    X = X(NDROP+1:end, :);
    Xf = bandpass_fft(X, fps, BAND);

    % The DETECTION threshold is per-run, tied to that run's own robust scale,
    % so a run whose amplitude sits far from the pooled reference still gets its
    % breaths found.  The amplitude that is REPORTED is then divided by the
    % shared baseline reference, so cross-run comparability is untouched.
    % Tying detection to the global reference instead makes low-amplitude runs
    % return no peaks and a NaN amplitude -- which is a detection failure
    % masquerading as a measurement.
    prom = MIN_PROM_F * arrayfun(@(j) iqr_p2p(Xf(:, j)), 1:nR);
    [pk, amp] = deal_detect(Xf, fps, MIN_SEP_SEC, prom);

    % Hand-curated events win over the automatic ones.  breathing_peak_gui_pc1
    % writes breath_peak_pc1.mat and breathing_trough_gui_pc1 writes
    % breath_insp_start_pc1.mat; if either is present for a run its indices
    % replace the automatic detection here.  Only the INDICES are taken -- the
    % GUI stores its amplitudes in z-units (it detrends and z-scores what it
    % loads), so amplitude is always re-measured on the physical-unit trace.
    if USE_CURATED
        [pk, amp, nCur, hasCur] = apply_curated(pk, amp, Xf, ROOT_DIR, runNames, NDROP);
        R(m).nCurated = nCur;
        R(m).hasCur = hasCur(:);
    else
        R(m).nCurated = 0;
        R(m).hasCur = false(nR, 1);
    end

    R(m).name = metrics{m};
    R(m).label = labels{m};
    R(m).X = Xf;
    R(m).pk = pk;
    R(m).amp = amp;
    % Band fraction MUST be measured on the raw trace.  Computing it on Xf
    % measures the bandpass filter, not the data, and returns ~1.0 for every
    % run including ones with no breathing at all.
    R(m).bandfrac = arrayfun(@(j) band_frac(X(:, j), fps, BAND), 1:nR);
    R(m).nBreath = cellfun(@numel, pk);
end

%% --- run selection: decided by the detection itself -------------------------
% CURATION IS AUTHORITATIVE.  As soon as any run has been curated in the
% peak/trough GUIs, only curated runs are analysed: an uncurated run is one you
% have not vetted, and filling it in with automatic detection would silently mix
% vetted and unvetted events in the same summary.  A run where the animal has
% stopped breathing simply never gets curated, so it drops out by itself.
%
% Before any curation exists the automatic detection is used, so a first pass
% still works and you can see the traces before committing time to the GUIs.
% Under curation all metrics share one event list, so metric 1 speaks for all.
anyCurated = any(R(1).hasCur);
if anyCurated
    isValid = R(1).hasCur & R(1).nBreath(:) >= MIN_BREATHS;
else
    isValid = R(1).nBreath(:) >= MIN_BREATHS;
end
if isempty(BASELINE_RUNS)
    isBase = logical(S.is_baseline(:));
    baseSrc = 'from extractor';
else
    isBase = ismember(runNum, BASELINE_RUNS(:));
    baseSrc = 'set here';
end
isBase = isBase & isValid;      % an excluded run can never define the reference
assert(any(isBase), ['no baseline runs left: none of BASELINE_RUNS reached ' ...
    'MIN_BREATHS=%d detected breaths.'], MIN_BREATHS);

%% --- pass 2: amplitude reference, now that the baseline set is known --------
for m = 1:numel(metrics)
    ref = median(cell2mat(R(m).amp(isBase)'), 'omitnan');
    R(m).ref = ref;
    R(m).ampNorm = cellfun(@(a) median(a, 'omitnan') / ref, R(m).amp);
    R(m).rate = cellfun(@(p) rate_from_peaks(p, fps), R(m).pk);
end

%% --- run selection report ---------------------------------------------------
if anyCurated
    fprintf('\nevents: CURATED (peak/trough GUI) -- %d of %d runs curated\n', ...
        R(1).nCurated, nR);
    fprintf('included: %d run(s)  %s\n', sum(isValid), compact_ranges(runNum(isValid)));
    unc = ~R(1).hasCur;
    if any(unc)
        fprintf('  excluded, not curated: %d run(s)  %s\n', ...
            sum(unc), compact_ranges(runNum(unc)));
    end
    thin = R(1).hasCur & ~isValid;
    if any(thin)
        fprintf('  excluded, curated but < %d events: %s\n', ...
            MIN_BREATHS, compact_ranges(runNum(thin)));
    end
else
    fprintf('\nevents: AUTOMATIC (no curated files found -- run the peak GUI)\n');
    fprintf('included: %d of %d runs with >= %d detected breaths\n', ...
        sum(isValid), nR, MIN_BREATHS);
    if any(~isValid)
        fprintf('  excluded: %s\n', compact_ranges(runNum(~isValid)));
    end
end
fprintf('baseline (%s): %s  -> %d run(s) define the amplitude reference\n', ...
    baseSrc, compact_ranges(runNum(isBase)), sum(isBase));

% Signal quality is REPORTED, never used to include or exclude anything.  It is
% here only to flag runs worth a second look in the GUI -- a run that passed
% MIN_BREATHS but looks poor may be one where detection latched onto noise.
qBand = R(1).bandfrac(:);
qDom  = arrayfun(@(j) dom_freq(R(1).X(:, j), fps), 1:nR)';
lowq = find(isValid & (qBand < 0.60 | qDom < 0.5 | qDom > 6));
if ~isempty(lowq)
    fprintf('  FYI, included but low quality (bandfrac<0.60 or rate outside 0.5-6 Hz): %s\n', ...
        compact_ranges(runNum(lowq)));
end

% The excluded runs are an empirical noise floor: whatever a metric reports
% there is what it reports with no breathing present.  That sets the smallest
% amplitude worth believing, which no amount of internal statistics can.
if any(~isValid) && any(isValid)
    fprintf('\nnoise floor from the %d excluded run(s):\n', sum(~isValid));
    for m = 1:numel(metrics)
        sdV = std(R(m).X(:, isValid), 0, 1);
        sdN = std(R(m).X(:, ~isValid), 0, 1);
        fprintf('  %-22s noise sd %.4g vs valid sd %.4g   -> floor = %.2f of a baseline breath\n', ...
            R(m).label, median(sdN), median(sdV), median(sdN) / R(m).ref);
    end
end

%% --- comparability report --------------------------------------------------
fprintf('\n%-22s %10s %10s %10s\n', 'metric', 'CV_base', 'ref', 'unit');
for m = 1:numel(metrics)
    a = R(m).ampNorm(isBase); a = a(isfinite(a));
    fprintf('%-22s %9.1f%% %10.4g %10s\n', R(m).label, ...
        100 * std(a) / mean(a), R(m).ref, tern(m == 1, 'px', 'a.u.'));
end

% CV_base above still contains the real session drift.  The successive-
% difference estimate is trend-robust, so that one is the honest measurement
% noise floor -- report both or the drift gets mistaken for imprecision.
fprintf('\nbaseline run-to-run amplitude scatter (CV_base includes real drift):\n');
for m = 1:numel(metrics)
    a = R(m).ampNorm(isBase); a = a(isfinite(a));
    nz = std(diff(a)) / sqrt(2) / mean(a);
    fprintf('  %-22s %.3f +/- %.3f  (n=%d)   noise-only CV = %.1f%%\n', ...
        R(m).label, mean(a), std(a), numel(a), 100 * nz);
end

% Does baseline amplitude drift with run order?  If yes, "baseline is a stable
% standard" is false for this prep.
bi = find(isBase);
for m = 1:numel(metrics)
    % Correlate against real elapsed time, not run index -- the two differ once
    % a run is skipped or the spacing is uneven.
    [rho, pv] = corr_spearman(tRun(bi), R(m).ampNorm(bi)');
    fprintf('  %-22s baseline drift vs wall-clock: rho=%+.3f p=%.3f%s\n', ...
        R(m).label, rho, pv, tern(pv < 0.05, '   <-- DRIFTS', ''));
end

fprintf('\nper-video basis rotation |u1 . pv1| (1 = same mode as frozen basis):\n');
fprintf('  min %.3f  median %.3f  max %.3f\n', ...
    min(S.pv_align), median(S.pv_align), max(S.pv_align));
nflip = sum(S.pv_align_sign < 0);
fprintf('  per-video PC1 sign disagreed with the frozen mode in %d/%d runs\n', nflip, nR);

%% --- figures ---------------------------------------------------------------
t = (0:size(R(1).X, 1) - 1)' / fps;

f1 = figure('Color', 'w', 'Position', [60 60 1500 900], 'Name', 'metric comparison');
tl = tiledlayout(f1, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Example traces: first and last baseline, plus the first included non-baseline
% run.  Excluded runs are never shown -- there is no breathing in them to see.
jShow = [find(isBase, 1), find(isBase, 1, 'last'), find(~isBase & isValid, 1)];
jShow = unique(jShow(~cellfun(@isempty, num2cell(jShow))), 'stable');
for m = 1:3
    nexttile;
    hold on;
    w = t >= 5 & t < 10;
    for k = 1:numel(jShow)
        x = R(m).X(:, jShow(k));
        plot(t(w), x(w) / R(m).ref + (k - 1) * 2.5, 'LineWidth', 0.9);
    end
    hold off; box off;
    title(R(m).label, 'FontWeight', 'normal');
    xlabel('time (s)');
    if m == 1, ylabel('amplitude / baseline p2p'); end
    if m == 1
        legend(runNames(jShow), 'Interpreter', 'none', 'Location', 'northoutside', ...
            'Box', 'off', 'FontSize', 7);
    end
end

for m = 1:3
    nexttile;
    a = R(m).ampNorm(:);
    hold on;
    sel = isBase;
    plot(tRun(sel), a(sel), 'o-', 'MarkerFaceColor', [.2 .4 .8], ...
        'Color', [.2 .4 .8], 'MarkerSize', 4);
    sel = ~isBase & isValid;
    if any(sel)
        plot(tRun(sel), a(sel), 'o-', 'MarkerFaceColor', [.85 .3 .2], ...
            'Color', [.85 .3 .2], 'MarkerSize', 4);
    end
    yline(1, ':', 'Color', [.4 .4 .4]);
    ab = a(isBase); ab = ab(isfinite(ab));
    mu = mean(ab); sd = std(ab);
    if isfinite(mu) && isfinite(sd) && sd > 0
        yline(mu + sd, '--', 'Color', [.7 .7 .7]);
        yline(mu - sd, '--', 'Color', [.7 .7 .7]);
    end
    hold off; box off;
    xlabel(tUnit); title(sprintf('CV_{base} = %.1f%%', 100 * sd / mu), ...
        'FontWeight', 'normal');
    if m == 1, ylabel('breath amplitude / baseline'); end
end

for m = 1:3
    nexttile;
    hold on;
    rr = R(m).rate(:);
    sel = isBase;
    plot(tRun(sel), rr(sel), 'o-', 'MarkerFaceColor', [.2 .4 .8], ...
        'Color', [.2 .4 .8], 'MarkerSize', 4);
    sel = ~isBase & isValid;
    if any(sel)
        plot(tRun(sel), rr(sel), 'o-', 'MarkerFaceColor', [.85 .3 .2], ...
            'Color', [.85 .3 .2], 'MarkerSize', 4);
    end
    hold off; box off;
    xlabel(tUnit);
    title(sprintf('band frac %.3f', mean(R(m).bandfrac)), 'FontWeight', 'normal');
    if m == 1, ylabel('breath rate (Hz)'); end
end
title(tl, sprintf('%s  --  blue = baseline, red = test  (%d of %d runs included)', ...
    ROOT_DIR, sum(isValid), nR), ...
    'Interpreter', 'none');

f2 = figure('Color', 'w', 'Position', [80 80 1200 420], 'Name', 'frozen basis');
subplot(1, 3, 1);
imagesc(S.ref_mean_img); axis image off; colormap(gca, gray);
title('reference mean frame', 'FontWeight', 'normal');
subplot(1, 3, 2);
u = S.u1_img; lim = max(abs(u(:)));
imagesc(u, [-lim lim]); axis image off; colormap(gca, redblue(256));
title(sprintf('frozen spatial mode u_1 (%.1f%% var)', 100 * S.basis_var_frac(1)), ...
    'FontWeight', 'normal');
subplot(1, 3, 3);
plot(S.reg_shift(:, 1), 'o-'); hold on; plot(S.reg_shift(:, 2), 's-'); hold off;
box off; xlabel('run index'); ylabel('drift (px)');
legend({'dy', 'dx'}, 'Box', 'off', 'Location', 'best');
title('FOV drift vs reference', 'FontWeight', 'normal');

%% --- per-breath scatter -----------------------------------------------------
% One dot per BREATH, not per run: instantaneous rate (1/IBI) and per-breath
% amplitude against session time.  The run-level panels above show medians,
% which hide the spread -- this shows every breath, so variability, outliers
% and abrupt transitions are visible rather than averaged away.
%
% Rate comes from the event times, so it is identical for all three metrics.
% Amplitude is read off `disp`, the metric in physical pixels.
mAmp = 1;                                  % 1=disp, 2=fb, 3=pv
secToT = 1 / 3600;                         % tRun in hours
if contains(tUnit, '(min)'), secToT = 1 / 60; end

% Typed empties: growing from [] would promote the logical mask to double, and
% a 0/1 double then indexes instead of masking.
bT = zeros(0, 1); bF = zeros(0, 1); bA = zeros(0, 1); bBase = false(0, 1);
for j = find(isValid(:))'
    p = double(R(mAmp).pk{j}(:));
    a = double(R(mAmp).amp{j}(:));
    if numel(p) < 2, continue; end
    ibi = diff(p) / fps;                   % s between consecutive breaths
    % amp(i) is peak i -> the trough before peak i+1, and ibi(i) is peak i ->
    % peak i+1, so entry i of both belongs to peak i.  Same length by
    % construction; the min() only guards a malformed curated file.
    n = min(numel(ibi), numel(a));
    bT = [bT; tRun(j) + (p(1:n) / fps) * secToT];      %#ok<AGROW>
    bF = [bF; 1 ./ ibi(1:n)];                          %#ok<AGROW>
    bA = [bA; a(1:n) / R(mAmp).ref];                   %#ok<AGROW>
    bBase = [bBase; repmat(logical(isBase(j)), n, 1)]; %#ok<AGROW>
end

ok3 = isfinite(bT) & isfinite(bF) & isfinite(bA);
bT = bT(ok3); bF = bF(ok3); bA = bA(ok3); bBase = bBase(ok3);

f3 = figure('Color', 'w', 'Position', [100 100 1250 760], 'Name', 'per-breath scatter');
tl3 = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
cB = [.20 .40 .80]; cT = [.85 .30 .20];

nexttile; hold on;
scatter(bT(bBase),  bF(bBase),  6, cB, 'filled', 'MarkerFaceAlpha', 0.25);
scatter(bT(~bBase), bF(~bBase), 6, cT, 'filled', 'MarkerFaceAlpha', 0.25);
hold off; box off; grid on;
ylabel('instantaneous rate  1/IBI  (Hz)');
title(sprintf('%d breaths from %d runs', numel(bF), sum(isValid)), 'FontWeight', 'normal');

nexttile; hold on;
scatter(bT(bBase),  bA(bBase),  6, cB, 'filled', 'MarkerFaceAlpha', 0.25);
scatter(bT(~bBase), bA(~bBase), 6, cT, 'filled', 'MarkerFaceAlpha', 0.25);
yline(1, ':', 'Color', [.4 .4 .4]);
hold off; box off; grid on;
xlabel(tUnit); ylabel('breath amplitude / baseline');
title(sprintf('amplitude from %s', R(mAmp).label), 'FontWeight', 'normal');
title(tl3, sprintf('%s  --  per-breath  (blue = baseline, red = test)', ROOT_DIR), ...
    'Interpreter', 'none');

if SAVE_FIGS
    save_fig(f1, fullfile(ROOT_DIR, 'breath_fixedmetric_comparison.png'));
    save_fig(f2, fullfile(ROOT_DIR, 'breath_fixedmetric_basis.png'));
    save_fig(f3, fullfile(ROOT_DIR, 'breath_fixedmetric_perbreath.png'));
    fprintf('\nsaved figures to %s\n', ROOT_DIR);
end

% Per-breath vectors are saved too, so you can re-plot or pool across sessions
% without re-running the detection.
breath = struct('t', bT, 'rate_hz', bF, 'amp_norm', bA, 'is_baseline', bBase, ...
    'amp_metric', R(mAmp).name, 'amp_ref', R(mAmp).ref, 't_unit', tUnit);

out = fullfile(ROOT_DIR, 'breath_fixedmetric_results.mat');
save(out, 'R', 'runNames', 'isBase', 'isValid', 'qBand', 'qDom', 'fps', 'BAND', 'MIN_SEP_SEC', 'MIN_PROM_F', 'NDROP', 'tRun', 'tUnit', 'breath');
fprintf('saved %s\n', out);

end

%% =============================== LOCAL FUNCTIONS ===============================

function Y = bandpass_fft(X, fps, band)
% Zero-phase brick-wall bandpass. No filter design toolbox dependency, and
% zero phase matters because breath timing is the downstream measurement.
T = size(X, 1);
X = X - mean(X, 1);
nfft = 2^nextpow2(2 * T);
F = fft(X, nfft, 1);
f = (0:nfft-1)' * (fps / nfft);
f(f > fps/2) = fps - f(f > fps/2);
keep = f >= band(1) & f <= band(2);
F(~keep, :) = 0;
Y = real(ifft(F, [], 1));
Y = Y(1:T, :);
Y = Y - mean(Y, 1);
end

function d = dom_freq(x, fps)
x = x(:) - mean(x);
n = numel(x);
w = 0.5 - 0.5 * cos(2 * pi * (0:n-1)' / (n - 1));
P = abs(fft(x .* w)).^2;
f = (0:n-1)' * (fps / n);
h = 1:floor(n/2);
P = P(h); f = f(h);
m = f > 0.3;
[~, i] = max(P(m));
fm = f(m);
d = fm(i);
end

function s = compact_ranges(v)
% "37-89" rather than 53 comma-separated numbers.
v = sort(v(:))';
if isempty(v), s = '(none)'; return; end
b = [true, diff(v) ~= 1];
starts = v(b); ends = v([b(2:end), true]);
parts = arrayfun(@(a, z) tern2(a == z, sprintf('%d', a), sprintf('%d-%d', a, z)), ...
    starts, ends, 'UniformOutput', false);
s = strjoin(parts, ', ');
end

function o = tern2(c, a, b)
if c, o = a; else, o = b; end
end

function save_fig(f, path)
% exportgraphics is preferred, but it refuses some figure states; falling back
% to print keeps a long analysis from dying at the very last step.
if ~isgraphics(f, 'figure')
    fprintf(2, '  figure handle invalid, not saved: %s\n', path);
    return;
end
try
    exportgraphics(f, path, 'Resolution', 200);
catch
    try
        print(f, path, '-dpng', '-r200');
    catch ME
        fprintf(2, '  could not save %s (%s)\n', path, ME.message);
    end
end
end

function [pk, amp, nCur, hasCur] = apply_curated(pk, amp, Xf, rootDir, runNames, nDrop)
% Replace automatic peaks with hand-curated ones where the GUIs have saved them.
%
% Index bookkeeping matters here: the GUIs work on the FULL trace written by
% breath_fixedmetric_to_gui, while this script has already dropped nDrop leading
% frames, so every curated index shifts by nDrop before use.  Getting that wrong
% would silently misplace every event by a fixed lag.
% Curated files are never deleted by a re-run, which is right -- they are hours
% of manual work.  But changing the ROI changes every trace, and the old events
% would then be applied to data they were not marked on.  The ROI file's
% timestamp settles it: anything curated before the current ROI was chosen is
% stale and is refused rather than used silently.
nCur = 0; nStale = 0;
hasCur = false(numel(pk), 1);
roiInfo = dir(fullfile(rootDir, 'breath_roi.mat'));
roiTime = 0;
if ~isempty(roiInfo), roiTime = roiInfo(1).datenum; end

for j = 1:numel(pk)
    pkFile = fullfile(rootDir, runNames{j}, 'breath_peak_pc1.mat');
    if ~isfile(pkFile), continue; end
    fi = dir(pkFile);
    if roiTime > 0 && fi(1).datenum < roiTime
        nStale = nStale + 1;
        continue;
    end
    P = load(pkFile);
    if ~isfield(P, 'insp_onset_idx') || isempty(P.insp_onset_idx), continue; end

    x = Xf(:, j);
    locs = round(double(P.insp_onset_idx(:))) - nDrop;
    locs = locs(locs >= 1 & locs <= numel(x));
    if numel(locs) < 2, continue; end
    locs = unique(sort(locs));

    % Troughs from the companion GUI if it has run; otherwise the minimum
    % between consecutive peaks, matching the automatic path.
    trFile = fullfile(rootDir, runNames{j}, 'breath_insp_start_pc1.mat');
    tr = [];
    if isfile(trFile)
        Tr = load(trFile);
        if isfield(Tr, 'insp_start_idx') && ~isempty(Tr.insp_start_idx)
            tr = round(double(Tr.insp_start_idx(:))) - nDrop;
            tr = tr(tr >= 1 & tr <= numel(x));
        end
    end

    a = nan(numel(locs) - 1, 1);
    for i = 1:numel(locs) - 1
        if isempty(tr)
            a(i) = x(locs(i)) - min(x(locs(i):locs(i+1)));
        else
            k = tr(tr > locs(i) & tr < locs(i+1));
            if isempty(k)
                a(i) = x(locs(i)) - min(x(locs(i):locs(i+1)));
            else
                a(i) = x(locs(i)) - min(x(k));
            end
        end
    end
    pk{j} = locs;
    amp{j} = a;
    hasCur(j) = true;
    nCur = nCur + 1;
end

if nStale > 0
    fprintf(2, ['\n  WARNING: ignored curated events for %d run(s) -- they predate\n' ...
        '           the current ROI, so they were marked on different traces.\n' ...
        '           Re-curate those runs, or restore the previous ROI.\n'], nStale);
end
end

function [pk, amp] = deal_detect(X, fps, minSepSec, minProm)
% Breath peaks + per-breath peak-to-trough amplitude, per column.
% minProm may be a scalar or one value per column.
N = size(X, 2);
pk = cell(1, N); amp = cell(1, N);
md = max(1, round(minSepSec * fps));
if isscalar(minProm), minProm = repmat(minProm, 1, N); end
for j = 1:N
    x = X(:, j);
    [~, locs] = findpeaks(x, 'MinPeakDistance', md, 'MinPeakProminence', minProm(j));
    pk{j} = locs;
    if numel(locs) < 2
        amp{j} = NaN;
        continue;
    end
    a = nan(numel(locs) - 1, 1);
    for i = 1:numel(locs) - 1
        a(i) = x(locs(i)) - min(x(locs(i):locs(i+1)));
    end
    amp{j} = a;
end
end

function v = iqr_p2p(x)
% Robust provisional amplitude scale before any peak is detected.
v = diff(prctile(x, [10 90]));
end

function r = rate_from_peaks(locs, fps)
if numel(locs) < 2, r = NaN; return; end
r = fps / median(diff(locs));
end

function v = band_frac(x, fps, band)
x = x(:) - mean(x);
n = numel(x);
w = 0.5 - 0.5 * cos(2 * pi * (0:n-1)' / (n - 1));   % Hann, inlined (no toolbox)
P = abs(fft(x .* w)).^2;
f = (0:n-1)' * (fps / n);
h = 1:floor(n/2);
P = P(h); f = f(h);
inb = f >= band(1) & f < band(2);
tot = f >= 0.1 & f < 0.98 * fps / 2;
v = sum(P(inb)) / max(sum(P(tot)), eps);
end

function [rho, p] = corr_spearman(x, y)
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = tiedrank_local(x(ok)); y = tiedrank_local(y(ok));
n = numel(x);
if n < 4, rho = NaN; p = NaN; return; end
rho = (mean(x .* y) - mean(x) * mean(y)) / (std(x, 1) * std(y, 1));
tstat = rho * sqrt((n - 2) / max(1 - rho^2, eps));
p = 2 * (1 - tcdf_local(abs(tstat), n - 2));
end

function r = tiedrank_local(x)
[~, i] = sort(x);
r = zeros(size(x));
r(i) = 1:numel(x);
end

function p = tcdf_local(t, v)
% Student-t CDF via the incomplete beta function (no Statistics toolbox call).
p = 1 - 0.5 * betainc(v / (v + t.^2), v / 2, 0.5);
end

function c = redblue(n)
h = floor(n / 2);
c = [[linspace(0.15, 1, h)', linspace(0.3, 1, h)', linspace(0.7, 1, h)']; ...
     [linspace(1, 0.8, n - h)', linspace(1, 0.25, n - h)', linspace(1, 0.15, n - h)']];
end

function o = tern(c, a, b)
if c, o = a; else, o = b; end
end
