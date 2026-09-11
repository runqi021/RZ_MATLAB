% whisk_insp_triggered_pooled_LR_angle_3panels_RZ.m
%
% Three panels, ALL triggered by inspiration event:
%
%   1) all inspiration events
%   2) inspiration events with a whisk onset within 100 ms
%   3) inspiration events with NO whisk onset within 100 ms
%
% All panels contain pooled L/R BP whisker angle.
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   BP 5-50 Hz, zero-phase filtfilt
%   whisk epoch mask from Hilbert envelope of LEFT BP whisker
%   hard envelope threshold = 5 deg
%
% Breathing:
%   load canonical _breath.mat / Pn.breath
%   Bs.breath is already LP<1 Hz subtracted AND already inverted
%   inhale is already positive/upward
%   DO NOT invert again
%   BP 1-20 Hz, zero-phase filtfilt
%
% Inspiration event:
%   findpeaks(brw_bp1_20)
%
% Classification:
%   inspiration with nearby whisk = any L/R whisk onset within 100 ms
%   inspiration without nearby whisk = no L/R whisk onset within 100 ms
%
% Plot:
%   mean pooled L/R BP angle +/- pointwise 95% CI

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz
BREATH_BP = [1 20];   % Hz

ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

PRE_S  = 0.30;        % seconds before inspiration event
POST_S = 0.30;        % seconds after inspiration event

MIN_INSP_DIST_S = 0.06;   % minimum distance between inspiration events
PROM_MAD = 0.5;           % breath peak prominence = PROM_MAD * MAD(event trace)

INSP_WHISK_WIN_S = 0.05;  % 100 ms rule
MIN_WHISK_DIST_S = 0.03;  % minimum distance between detected whisk onsets

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

tRel = (-round(PRE_S*fpsW):round(POST_S*fpsW))' / fpsW;
nRel = numel(tRel);

POOL_allInsp       = zeros(0,nRel);
POOL_inspWithWhisk = zeros(0,nRel);
POOL_inspNoWhisk   = zeros(0,nRel);

nSess = 0;
nAllInsp = 0;
nWithWhisk = 0;
nNoWhisk = 0;

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

    xL = filtfilt(b_wbp, a_wbp, La0);   % L BP angle
    xR = filtfilt(b_wbp, a_wbp, Ra0);   % R BP angle

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

        event_trace = brw_bp1_20;  % inspiration events are positive peaks

        prom = PROM_MAD * mad(event_trace, 1);
        if ~isfinite(prom) || prom <= 0
            prom = std(event_trace, 'omitnan') * 0.5;
        end

        [~, locs] = findpeaks(-event_trace, ...
            'MinPeakDistance', round(MIN_INSP_DIST_S * fpsW), ...
            'MinPeakProminence', prom);

        insp_t = tW(locs);

    catch ME
        warning('breath load/detect failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    if isempty(insp_t)
        warning('no inspiration events detected %s n%d', animal, kRun);
        continue;
    end

    % ======================= WHISK ONSETS FOR GROUPING ==================

    tOnL = detect_upcrossings(xL, tW, m, MIN_WHISK_DIST_S);
    tOnR = detect_upcrossings(xR, tW, m, MIN_WHISK_DIST_S);

    tOnAll = sort([tOnL(:); tOnR(:)]);

    % ======================= INSPIRATION-TRIGGERED SEGMENTS ============

    nAllSess = 0;
    nWithSess = 0;
    nNoSess = 0;

    for ii = 1:numel(insp_t)

        t0 = insp_t(ii);
        tt = t0 + tRel;

        if tt(1) < tW(1) || tt(end) > tW(end)
            continue;
        end

        segL = interp1(tW, xL, tt, 'linear', NaN);
        segR = interp1(tW, xR, tt, 'linear', NaN);

        if mean(isfinite(segL)) < 0.95 || mean(isfinite(segR)) < 0.95
            continue;
        end

        % classify this inspiration by whether any L/R whisk onset is nearby
        if isempty(tOnAll)
            hasNearbyWhisk = false;
        else
            hasNearbyWhisk = min(abs(tOnAll - t0)) <= INSP_WHISK_WIN_S;
        end

        % panel 1: all inspirations, pool L and R
        POOL_allInsp(end+1,:) = segL(:)'; %#ok<SAGROW>
        POOL_allInsp(end+1,:) = segR(:)'; %#ok<SAGROW>
        nAllSess = nAllSess + 2;

        % panel 2/3: same inspiration trigger, grouped by nearby whisk onset
        if hasNearbyWhisk
            POOL_inspWithWhisk(end+1,:) = segL(:)'; %#ok<SAGROW>
            POOL_inspWithWhisk(end+1,:) = segR(:)'; %#ok<SAGROW>
            nWithSess = nWithSess + 2;
        else
            POOL_inspNoWhisk(end+1,:) = segL(:)'; %#ok<SAGROW>
            POOL_inspNoWhisk(end+1,:) = segR(:)'; %#ok<SAGROW>
            nNoSess = nNoSess + 2;
        end
    end

    if nAllSess > 0
        nSess = nSess + 1;
    end

    nAllInsp   = nAllInsp   + nAllSess;
    nWithWhisk = nWithWhisk + nWithSess;
    nNoWhisk   = nNoWhisk   + nNoSess;

    fprintf('%s n%d: all=%d, insp+whisk=%d, insp-no-whisk=%d pooled L/R segments\n', ...
        animal, kRun, nAllSess, nWithSess, nNoSess);
end

fprintf('\nSessions used: %d\n', nSess);
fprintf('All inspiration-triggered pooled segments:       %d\n', nAllInsp);
fprintf('Inspirations WITH nearby whisk pooled segments:  %d\n', nWithWhisk);
fprintf('Inspirations WITHOUT nearby whisk pooled segments: %d\n', nNoWhisk);

% ============================== PLOT ==============================

figure('Color','w','Position',[100 200 1400 420]);

ax1 = subplot(1,3,1); hold(ax1,'on'); grid(ax1,'on');
plot_triggered_mean_ci95(ax1, tRel, POOL_allInsp, [0.15 0.15 0.15]);
xline(ax1, 0, 'k--', 'inspiration');
yline(ax1, 0, 'k:');
xlabel(ax1, 'time from inspiration (s)');
ylabel(ax1, 'pooled L/R BP angle (deg)');
title(ax1, sprintf('All inspiration-triggered whisk angle\nmean +/- 95%% CI, n=%d pooled L/R segments', ...
    size(POOL_allInsp,1)), 'Interpreter','none');
axis(ax1, 'square');
box(ax1, 'off');

ax2 = subplot(1,3,2); hold(ax2,'on'); grid(ax2,'on');
plot_triggered_mean_ci95(ax2, tRel, POOL_inspWithWhisk, [0.15 0.15 0.15]);
xline(ax2, 0, 'k--', 'inspiration');
yline(ax2, 0, 'k:');
xlabel(ax2, 'time from inspiration (s)');
ylabel(ax2, 'pooled L/R BP angle (deg)');
title(ax2, sprintf('Inspirations with nearby whisk onset\nwhisk onset within %.0f ms, mean +/- 95%% CI, n=%d pooled L/R segments', ...
    1000*INSP_WHISK_WIN_S, size(POOL_inspWithWhisk,1)), 'Interpreter','none');
axis(ax2, 'square');
box(ax2, 'off');

ax3 = subplot(1,3,3); hold(ax3,'on'); grid(ax3,'on');
plot_triggered_mean_ci95(ax3, tRel, POOL_inspNoWhisk, [0.15 0.15 0.15]);
xline(ax3, 0, 'k--', 'inspiration');
yline(ax3, 0, 'k:');
xlabel(ax3, 'time from inspiration (s)');
ylabel(ax3, 'pooled L/R BP angle (deg)');
title(ax3, sprintf('Inspirations without nearby whisk onset\nnearest whisk onset > %.0f ms, mean +/- 95%% CI, n=%d pooled L/R segments', ...
    1000*INSP_WHISK_WIN_S, size(POOL_inspNoWhisk,1)), 'Interpreter','none');
axis(ax3, 'square');
box(ax3, 'off');

sgtitle(sprintf(['All panels triggered by inspiration: whisk BP %g-%g Hz, breath BP %g-%g Hz, ' ...
                 'no breath inversion, shading = pointwise 95%% CI'], ...
                 WHISK_BP(1), WHISK_BP(2), BREATH_BP(1), BREATH_BP(2)));

% ============================= HELPERS =============================

function plot_triggered_mean_ci95(ax, tRel, X, col)
% Plot mean +/- pointwise 95% confidence interval.
%
% X is events/segments x time.
% This treats pooled segments as samples.
% For animal/session-level inference, average within session first.

    if isempty(X)
        title(ax, 'no events');
        return;
    end

    mu = mean(X, 1, 'omitnan');

    n  = sum(isfinite(X), 1);
    sd = std(X, 0, 1, 'omitnan');
    sem = sd ./ sqrt(n);

    ci95 = nan(size(mu));

    ok = n > 1 & isfinite(sem);

    % Use t critical value when available.
    % If Statistics Toolbox is unavailable, fall back to 1.96.
    if exist('tinv', 'file') == 2 || exist('tinv', 'builtin') == 5
        tcrit = nan(size(mu));
        tcrit(ok) = tinv(0.975, n(ok)-1);
        ci95(ok) = tcrit(ok) .* sem(ok);
    else
        ci95(ok) = 1.96 .* sem(ok);
    end

    upper = mu + ci95;
    lower = mu - ci95;

    fill(ax, [tRel; flipud(tRel)], ...
        [upper(:); flipud(lower(:))], ...
        col, 'FaceAlpha', 0.22, 'EdgeColor', 'none');

    plot(ax, tRel, mu, 'Color', col, 'LineWidth', 2);
    xlim([-0.2, 0.2]);
end

function tOn = detect_upcrossings(x, t, mask, minDistS)
% Upward zero-crossings of BP angle, restricted to valid whisking mask.

    x = x(:);
    t = t(:);
    mask = logical(mask(:));

    i0 = find(mask(1:end-1) & mask(2:end) & ...
              isfinite(x(1:end-1)) & isfinite(x(2:end)) & ...
              x(1:end-1) < 0 & x(2:end) >= 0);

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

    % prune crossings that are too close
    if nargin >= 4 && ~isempty(minDistS) && numel(tOn) > 1
        keep = true(size(tOn));
        last = tOn(1);

        for k = 2:numel(tOn)
            if (tOn(k) - last) < minDistS
                keep(k) = false;
            else
                last = tOn(k);
            end
        end

        tOn = tOn(keep);
    end
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