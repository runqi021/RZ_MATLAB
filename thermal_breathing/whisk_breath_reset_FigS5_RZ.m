% whisk_breath_reset_population_blackdots_ms_RZ_FIXED.m
%
% Figure S5-style population plot:
%
% A) Whisking rhythm response to breathing
%    x = Tr: inspiration onset time within whisk cycle, in ms
%    y = shift in whisking period by breathing, T1 - T0, in ms
%
% B) Breathing rhythm response to whisking
%    x = Tr: whisk protraction time within breathing cycle, in ms
%    y = shift in breathing period by whisking, T1 - T0, in ms
%
% Important fixes:
%   1) Breath event is inspiration ONSET, not breath peak.
%   2) Panel B uses ALL whisk events within each breath cycle, not only first.
%   3) X-axis is raw time in ms, not normalized.
%   4) Default whisk event train uses LEFT whisker only to avoid fake L/R double cycles.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz
BREATH_BP = [1 20];   % Hz

% Breath event mode:
%   "onset" = upward zero-crossing before each positive breath peak
%   "peak"  = positive breath peak itself
BREATH_EVENT_MODE = "onset";

% Whisk protraction event:
%   "left"     = left whisker only, safest for rhythm reset analysis
%   "right"    = right whisker only
%   "mergedLR" = merge L/R events
WHISK_EVENT_SOURCE = "right";

% Whisk crossing direction:
%   "up"   = upward zero crossing of BP angle
%   "down" = downward zero crossing of BP angle
%
WHISK_CROSSING_DIRECTION = "up";

% Whisk epoch mask
ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

% Breath peak detection used to find inspiration cycles
MIN_INSP_DIST_S = 0.06;
PROM_MAD = 0.5;

% Whisk event detection
MIN_WHISK_DIST_S = 0.03;
MERGE_LR_WHISK_S = 0.060;   % larger merge window if using mergedLR

% Cycle sanity limits
MIN_WHISK_PERIOD_S  = 0.025; % 40 Hz max
MAX_WHISK_PERIOD_S  = 0.20; % 4 Hz min

MIN_BREATH_PERIOD_S = 0.060; % 16.7 Hz max
MAX_BREATH_PERIOD_S = 0.600; % 1.67 Hz min

% Reset logic
USE_ALL_BREATHS_WITHIN_WHISK_CYCLE = true;
USE_ALL_WHISKS_WITHIN_BREATH_CYCLE = true;

EXCLUDE = "5840027";

% Plot settings
DOT_SIZE  = 1;
DOT_ALPHA = 1;

XLIM_A = [0 200];
YLIM_A = [-150 150];

XLIM_B = [0 800];
YLIM_B = [-400 300];

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

R_whisk_by_breath = init_reset_struct();
R_breath_by_whisk = init_reset_struct();

nSess = 0;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1};
    kRun   = sess{e}{2};

    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE))
        continue;
    end

    % ============================ WHISKER ============================

    try
        whiskCsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun));
        M = dlc_gate_interp(whiskCsv, 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

    tW = (0:numel(La)-1)' / fpsW;
    N  = numel(tW);

    La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
    Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');

    La0 = La0 - mean(La0, 'omitnan');
    Ra0 = Ra0 - mean(Ra0, 'omitnan');

    xL = filtfilt(b_wbp, a_wbp, La0);
    xR = filtfilt(b_wbp, a_wbp, Ra0);

    % ======================= WHISK EPOCH MASK ========================

    env_raw = abs(hilbert(xL));
    env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));

    is_whisk = env > ENV_THR;
    ep = bool_to_epochs(is_whisk, tW, MIN_DUR, MERGE_GAP);

    m = false(N,1);
    for q = 1:size(ep,1)
        m = m | (tW >= ep(q,1) & tW <= ep(q,2));
    end

    if nnz(m) < 50
        warning('no valid whisk epoch %s n%d', animal, kRun);
        continue;
    end

    % Candidate protraction events
    tOnL = detect_zero_crossings(xL, tW, m, MIN_WHISK_DIST_S, WHISK_CROSSING_DIRECTION);
    tOnR = detect_zero_crossings(xR, tW, m, MIN_WHISK_DIST_S, WHISK_CROSSING_DIRECTION);

    switch lower(char(WHISK_EVENT_SOURCE))
        case 'left'
            tWhisk = tOnL;
        case 'right'
            tWhisk = tOnR;
        case 'mergedlr'
            tWhisk = merge_close_events(sort([tOnL(:); tOnR(:)]), MERGE_LR_WHISK_S);
        otherwise
            error('Unknown WHISK_EVENT_SOURCE: %s', WHISK_EVENT_SOURCE);
    end

    tWhisk = sort(tWhisk(:));
    tWhisk = tWhisk(isfinite(tWhisk));

    if numel(tWhisk) < 5
        warning('too few whisk events %s n%d', animal, kRun);
        continue;
    end

    % ============================ BREATH =============================

    try
        noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun));
        Pn = thermal_resolve_paths(noseCsv, dataRoot);

        if ~isfile(Pn.breath)
            warning('canonical LP-sub breath file not found %s n%d', animal, kRun);
            continue;
        end

        Bs = load(Pn.breath);

        br = Bs.breath(:);      % already LP<1 Hz subtracted
                                % already inverted: inhale = positive
        fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb;

        brw = interp1(tB, br, tW, 'linear', NaN);

        if nnz(isfinite(brw)) < 50
            warning('too few finite breath samples %s n%d', animal, kRun);
            continue;
        end

        brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');

        % IMPORTANT: DO NOT invert.
        brw_bp1_20 = filtfilt(b_br, a_br, brw);

        tInsp = detect_inspiration_events( ...
            brw_bp1_20, tW, ...
            MIN_INSP_DIST_S, PROM_MAD, ...
            BREATH_EVENT_MODE);

    catch ME
        warning('breath load/detect failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    tInsp = sort(tInsp(:));
    tInsp = tInsp(isfinite(tInsp));

    if numel(tInsp) < 5
        warning('too few inspiration events %s n%d', animal, kRun);
        continue;
    end

    % ========================== RESET ANALYSIS ==========================

    % A) Whisk rhythm response to inspiration.
    % Perturbed rhythm = whisk events.
    % Interrupting event = inspiration onset.
    R1 = collect_reset_effects( ...
        tWhisk, tInsp, ...
        MIN_WHISK_PERIOD_S, MAX_WHISK_PERIOD_S, ...
        USE_ALL_BREATHS_WITHIN_WHISK_CYCLE);

    % B) Breathing rhythm response to whisk protraction.
    % Perturbed rhythm = inspiration onsets.
    % Interrupting event = whisk protractions.
    %
    % IMPORTANT:
    % Use ALL whisks within the breath cycle, otherwise x-values collapse
    % near the first whisk after inspiration.
    R2 = collect_reset_effects( ...
        tInsp, tWhisk, ...
        MIN_BREATH_PERIOD_S, MAX_BREATH_PERIOD_S, ...
        USE_ALL_WHISKS_WITHIN_BREATH_CYCLE);

    R1.session = repmat({sprintf('%s_n%d', animal, kRun)}, size(R1.dT));
    R2.session = repmat({sprintf('%s_n%d', animal, kRun)}, size(R2.dT));

    R_whisk_by_breath = append_reset_struct(R_whisk_by_breath, R1);
    R_breath_by_whisk = append_reset_struct(R_breath_by_whisk, R2);

    nSess = nSess + 1;

    fprintf('%s n%d: whisk events=%d, insp events=%d, whisk<-breath resets=%d, breath<-whisk resets=%d\n', ...
        animal, kRun, numel(tWhisk), numel(tInsp), numel(R1.dT), numel(R2.dT));
end

fprintf('\nSessions used: %d\n', nSess);
fprintf('Whisk rhythm response to inspiration: %d reset events\n', numel(R_whisk_by_breath.dT));
fprintf('Breath rhythm response to whisking:   %d reset events\n', numel(R_breath_by_whisk.dT));

print_reset_summary('Panel A whisk<-breath', R_whisk_by_breath);
print_reset_summary('Panel B breath<-whisk', R_breath_by_whisk);

% ============================== PLOT ==============================

figure('Color','w','Position',[200 80 760 960]);

% ----------------------------- PANEL A -----------------------------

ax1 = subplot(2,1,1); hold(ax1,'on');

xA = R_whisk_by_breath.treset * 1000;   % ms, raw time
yA = R_whisk_by_breath.dT     * 1000;   % ms

goodA = isfinite(xA) & isfinite(yA);
xA = xA(goodA);
yA = yA(goodA);

scatter_black_alpha(ax1, xA, yA, DOT_SIZE, DOT_ALPHA);
plot(ax1, XLIM_A, [0 0], 'k-', 'LineWidth', 1.5);

xlim(ax1, XLIM_A);
ylim(ax1, YLIM_A);

xlabel(ax1, 'Tr: Time of breath within whisk cycle (ms)');
ylabel(ax1, sprintf('Shift in whisking period by breathing\nT1 - T0 (ms)'));

set(ax1, 'FontSize', 14, 'LineWidth', 1.5, ...
    'TickDir','out', 'Box','off');

text(ax1, -0.27, 1.13, 'a', ...
    'Units','normalized', ...
    'FontSize', 34, ...
    'FontWeight','bold');

% ----------------------------- PANEL B -----------------------------

ax2 = subplot(2,1,2); hold(ax2,'on');

xB = R_breath_by_whisk.treset * 1000;   % ms, raw time
yB = R_breath_by_whisk.dT     * 1000;   % ms

goodB = isfinite(xB) & isfinite(yB);
xB = xB(goodB);
yB = yB(goodB);

scatter_black_alpha(ax2, xB, yB, DOT_SIZE, DOT_ALPHA);
plot(ax2, XLIM_B, [0 0], 'k-', 'LineWidth', 1.5);

xlim(ax2, XLIM_B);
ylim(ax2, YLIM_B);

xlabel(ax2, 'Tr: Time of whisk within breath cycle (ms)');
ylabel(ax2, sprintf('Shift in breathing period by whisking\nT1 - T0 (ms)'));

set(ax2, 'FontSize', 14, 'LineWidth', 1.5, ...
    'TickDir','out', 'Box','off');

text(ax2, -0.27, 1.13, 'b', ...
    'Units','normalized', ...
    'FontSize', 34, ...
    'FontWeight','bold');

% saveas(gcf, 'FigS5_population_blackdots_ms_FIXED_RZ.png');

% ============================= HELPERS =============================

function R = init_reset_struct()
    R.treset = [];
    R.T0 = [];
    R.T1 = [];
    R.dT = [];
    R.session = {};
end

function Rout = append_reset_struct(Rout, Rin)

    f = {'treset','T0','T1','dT'};

    for i = 1:numel(f)
        Rout.(f{i}) = [Rout.(f{i}); Rin.(f{i})(:)];
    end

    if isfield(Rin,'session') && ~isempty(Rin.session)
        Rout.session = [Rout.session(:); Rin.session(:)];
    end
end

function R = collect_reset_effects(Aevents, Bevents, minPeriod, maxPeriod, useAllInterrupters)
% Aevents = event times defining the perturbed rhythm.
% Bevents = interrupting/resetting event times.
%
% For cycle A(i) -> A(i+1):
%   T0 = A(i)   - A(i-1)
%   T1 = A(i+1) - A(i)
%   B occurs inside [A(i), A(i+1)]
%   dT = T1 - T0
%   treset = B - A(i)
%
% useAllInterrupters:
%   true  = use every B event inside this A cycle
%   false = use only first B event inside this A cycle

    Aevents = sort(Aevents(:));
    Bevents = sort(Bevents(:));

    R = init_reset_struct();

    if numel(Aevents) < 3 || isempty(Bevents)
        return;
    end

    for i = 2:(numel(Aevents)-1)

        aPrev = Aevents(i-1);
        a0    = Aevents(i);
        a1    = Aevents(i+1);

        T0 = a0 - aPrev;
        T1 = a1 - a0;

        if ~isfinite(T0) || ~isfinite(T1)
            continue;
        end

        if T0 < minPeriod || T0 > maxPeriod || T1 < minPeriod || T1 > maxPeriod
            continue;
        end

        ix = find(Bevents > a0 & Bevents < a1);

        if isempty(ix)
            continue;
        end

        if ~useAllInterrupters
            ix = ix(1);
        end

        for k = 1:numel(ix)

            b = Bevents(ix(k));
            treset = b - a0;

            if treset <= 0 || treset >= T1
                continue;
            end

            R.treset(end+1,1) = treset;
            R.T0(end+1,1)     = T0;
            R.T1(end+1,1)     = T1;
            R.dT(end+1,1)     = T1 - T0;
        end
    end
end

function tInsp = detect_inspiration_events(event_trace, t, minDistS, promMad, modeStr)
% Detect inspiration timing.
%
% "peak":
%   positive peak of breath trace.
%
% "onset":
%   upward zero-crossing immediately before each positive peak.
%   This is usually what you want for reset timing.

    event_trace = event_trace(:);
    t = t(:);

    prom = promMad * mad(event_trace, 1);
    if ~isfinite(prom) || prom <= 0
        prom = std(event_trace, 'omitnan') * 0.5;
    end

    [~, peakLocs] = findpeaks(event_trace, ...
        'MinPeakDistance', round(minDistS / median(diff(t))), ...
        'MinPeakProminence', prom);

    if isempty(peakLocs)
        tInsp = zeros(0,1);
        return;
    end

    switch lower(char(modeStr))

        case 'peak'
            tInsp = t(peakLocs);

        case 'onset'

            % upward zero crossings of breath trace
            zc = find(event_trace(1:end-1) < 0 & event_trace(2:end) >= 0);

            if isempty(zc)
                warning('No breath upward zero-crossings found; falling back to peaks.');
                tInsp = t(peakLocs);
                return;
            end

            tZC = nan(numel(zc),1);

            for j = 1:numel(zc)
                ii = zc(j);
                y1 = event_trace(ii);
                y2 = event_trace(ii+1);
                t1 = t(ii);
                t2 = t(ii+1);

                if y2 ~= y1
                    frac = -y1 / (y2 - y1);
                    tZC(j) = t1 + frac * (t2 - t1);
                else
                    tZC(j) = t1;
                end
            end

            tInsp = nan(numel(peakLocs),1);

            for p = 1:numel(peakLocs)
                iiPeak = peakLocs(p);

                jj = find(zc < iiPeak, 1, 'last');

                if isempty(jj)
                    tInsp(p) = t(iiPeak);
                else
                    tInsp(p) = tZC(jj);
                end
            end

            tInsp = tInsp(isfinite(tInsp));
            tInsp = unique(round(tInsp * 1e6) / 1e6);
            tInsp = enforce_min_interval(tInsp, minDistS);

        otherwise
            error('Unknown BREATH_EVENT_MODE: %s', modeStr);
    end
end

function scatter_black_alpha(ax, x, y, dotSize, alphaVal)

    if isempty(x)
        text(ax, 0.5, 0.5, 'no reset events', ...
            'Units','normalized', ...
            'HorizontalAlignment','center');
        return;
    end

    try
        scatter(ax, x, y, dotSize, 'k', 'filled', ...
            'MarkerFaceAlpha', alphaVal, ...
            'MarkerEdgeAlpha', alphaVal);
    catch
        scatter(ax, x, y, dotSize, 'k', 'filled');
    end
end

function tOn = detect_zero_crossings(x, t, mask, minDistS, directionStr)
% Zero-crossings of BP angle, restricted to valid whisking mask.

    x = x(:);
    t = t(:);
    mask = logical(mask(:));

    switch lower(char(directionStr))
        case 'up'
            i0 = find(mask(1:end-1) & mask(2:end) & ...
                      isfinite(x(1:end-1)) & isfinite(x(2:end)) & ...
                      x(1:end-1) < 0 & x(2:end) >= 0);
        case 'down'
            i0 = find(mask(1:end-1) & mask(2:end) & ...
                      isfinite(x(1:end-1)) & isfinite(x(2:end)) & ...
                      x(1:end-1) > 0 & x(2:end) <= 0);
        otherwise
            error('Unknown crossing direction: %s', directionStr);
    end

    if isempty(i0)
        tOn = zeros(0,1);
        return;
    end

    tOn = zeros(numel(i0),1);

    for k = 1:numel(i0)
        i = i0(k);

        x1 = x(i);
        x2 = x(i+1);
        t1 = t(i);
        t2 = t(i+1);

        if x2 ~= x1
            frac = -x1 / (x2 - x1);
            tOn(k) = t1 + frac * (t2 - t1);
        else
            tOn(k) = t1;
        end
    end

    if nargin >= 4 && ~isempty(minDistS) && numel(tOn) > 1
        tOn = enforce_min_interval(tOn, minDistS);
    end
end

function tOut = merge_close_events(tIn, mergeWinS)

    tIn = sort(tIn(:));
    tIn = tIn(isfinite(tIn));

    if isempty(tIn)
        tOut = zeros(0,1);
        return;
    end

    clusters = {};
    cur = tIn(1);

    for i = 2:numel(tIn)

        if tIn(i) - cur(end) <= mergeWinS
            cur(end+1,1) = tIn(i); %#ok<AGROW>
        else
            clusters{end+1} = cur; %#ok<AGROW>
            cur = tIn(i);
        end
    end

    clusters{end+1} = cur;

    tOut = zeros(numel(clusters),1);

    for i = 1:numel(clusters)
        tOut(i) = median(clusters{i});
    end
end

function tKeep = enforce_min_interval(t, minDistS)

    t = sort(t(:));
    t = t(isfinite(t));

    if isempty(t)
        tKeep = zeros(0,1);
        return;
    end

    keep = true(size(t));
    last = t(1);

    for k = 2:numel(t)
        if (t(k) - last) < minDistS
            keep(k) = false;
        else
            last = t(k);
        end
    end

    tKeep = t(keep);
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)

    a = logical(a(:));
    d = diff([false; a; false]);

    s = find(d == 1);
    e = find(d == -1) - 1;

    ep = [t(s), t(e)];

    if isempty(ep)
        return;
    end

    ep2 = ep(1,:);

    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= mergeGap
            ep2(end,2) = ep(i,2);
        else
            ep2(end+1,:) = ep(i,:); %#ok<AGROW>
        end
    end

    ep = ep2;
    ep = ep(ep(:,2) - ep(:,1) >= minDur, :);
end

function print_reset_summary(nameStr, R)

    fprintf('\n%s\n', nameStr);

    if isempty(R.dT)
        fprintf('  no events\n');
        return;
    end

    x = R.treset * 1000;
    T0 = R.T0 * 1000;
    T1 = R.T1 * 1000;
    dT = R.dT * 1000;

    fprintf('  n = %d\n', numel(dT));
    fprintf('  Tr ms: min %.2f, median %.2f, max %.2f\n', ...
        min(x), median(x), max(x));
    fprintf('  T0 ms: min %.2f, median %.2f, max %.2f\n', ...
        min(T0), median(T0), max(T0));
    fprintf('  T1 ms: min %.2f, median %.2f, max %.2f\n', ...
        min(T1), median(T1), max(T1));
    fprintf('  dT ms: min %.2f, median %.2f, max %.2f\n', ...
        min(dT), median(dT), max(dT));
end

function S = list_sessions(dirPath)

    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv'));

    S = {};
    key = {};

    for i = 1:numel(d)

        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');

        if isempty(tok)
            continue;
        end

        k = sprintf('%s_%s', tok{1}, tok{2});

        if any(strcmp(key,k))
            continue;
        end

        key{end+1} = k; %#ok<AGROW>
        S{end+1} = {tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end

    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S);
        rn = cellfun(@(c) c{2}, S);
        [~,o] = sortrows([an(:) rn(:)]);
        S = S(o);
    end
end

function csv = pick_csv(dirPath, prefix)

    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);

    [~,ix] = max(arrayfun(@(x) bestnum(x.name), d));
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)

    tok = regexp(name, 'best-(\d+)', 'tokens');

    if isempty(tok)
        n = 0;
    else
        n = str2double(tok{1}{1});
    end
end