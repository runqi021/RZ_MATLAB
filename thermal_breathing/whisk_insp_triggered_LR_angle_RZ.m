% whisk_insp_triggered_LR_angle_RZ.m
%
% Goal:
%   Detect inspiration events from canonical LP-subtracted breathing trace,
%   then trigger-average L/R whisker ANGLE around those events.
%
% Output:
%   Figure with two plots:
%       1) inspiration-triggered L BP angle
%       2) inspiration-triggered R BP angle
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   BP 5-50 Hz, zero-phase filtfilt
%   whisk epoch mask from Hilbert envelope of LEFT BP whisker
%   hard envelope threshold = 5 deg
%
% Breathing:
%   load _breath.mat / Pn.breath
%   Bs.breath is already LP<1 Hz subtracted AND already inverted
%   inhale is already positive/upward
%   DO NOT invert again
%   BP 1-20 Hz, zero-phase filtfilt
%
% Event:
%   inspiration event = findpeaks(brw_bp1_20)
%
% Event time:
%   t = 0 is detected inspiration peak/event.
%
% Animal 5840027 excluded by default.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz, whisker angle band
BREATH_BP = [1 20];   % Hz, breath event-detection band

ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

PRE_S  = 0.30;        % seconds before inspiration event
POST_S = 0.30;        % seconds after inspiration event

MIN_INSP_DIST_S = 0.06;   % minimum distance between inspiration events
PROM_MAD = 0.5;           % peak prominence = PROM_MAD * MAD(event trace)

MIN_MASK_FRAC = 0.80;     % fraction of trigger window that must be in whisk epoch

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

tRel = (-round(PRE_S*fpsW):round(POST_S*fpsW))' / fpsW;
nRel = numel(tRel);

L_evt = zeros(0,nRel);
R_evt = zeros(0,nRel);

nSess = 0;
nEvtTotal = 0;

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

    % ANGLE, not set-point:
    % BP-filtered whisker angle, not LP<1 Hz setpoint.
    xL = filtfilt(b_wbp, a_wbp, La0);   % L BP 5-50 Hz angle
    xR = filtfilt(b_wbp, a_wbp, Ra0);   % R BP 5-50 Hz angle

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

        br = Bs.breath(:);       % already LP<1 Hz subtracted
                                  % already inverted: inhale = positive
        fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb;

        brw = interp1(tB, br, tW, 'linear', NaN);

        if nnz(isfinite(brw)) < 50
            warning('too few finite breath samples %s n%d', animal, kRun);
            continue;
        end

        brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');

        % IMPORTANT:
        % DO NOT invert.
        % _breath.mat is already inverted, so inhale is already positive.
        brw_bp1_20 = filtfilt(b_br, a_br, brw);

        event_trace = brw_bp1_20;  % inspiration peaks are positive

        prom = PROM_MAD * mad(event_trace, 1);
        if ~isfinite(prom) || prom <= 0
            prom = std(event_trace, 'omitnan') * 0.5;
        end

        [~, locs] = findpeaks(event_trace, ...
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

    % ======================= TRIGGERED AVERAGE =======================

    nEvtSess = 0;

    for ii = 1:numel(insp_t)

        tt = insp_t(ii) + tRel;

        if tt(1) < tW(1) || tt(end) > tW(end)
            continue;
        end

        segL = interp1(tW, xL, tt, 'linear', NaN);
        segR = interp1(tW, xR, tt, 'linear', NaN);

        segM = interp1(tW, double(m), tt, 'nearest', 0) > 0.5;

        if mean(isfinite(segL)) < 0.95 || mean(isfinite(segR)) < 0.95
            continue;
        end

        if mean(segM) < MIN_MASK_FRAC
            continue;
        end

        L_evt(end+1,:) = segL(:)'; %#ok<SAGROW>
        R_evt(end+1,:) = segR(:)'; %#ok<SAGROW>

        nEvtSess = nEvtSess + 1;
    end

    if nEvtSess > 0
        nSess = nSess + 1;
        nEvtTotal = nEvtTotal + nEvtSess;
        fprintf('%s n%d: %d inspiration events used\n', animal, kRun, nEvtSess);
    end
end

fprintf('\nSessions used: %d\n', nSess);
fprintf('Total inspiration events used: %d\n', nEvtTotal);

% ============================== PLOT ==============================

figure('Color','w','Position',[160 220 980 420]);

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
plot_triggered_mean_sem(ax1, tRel, L_evt, [0.15 0.15 0.15]);
xline(ax1, 0, 'k--', 'inspiration');
yline(ax1, 0, 'k:');
xlabel(ax1, 'time from inspiration event (s)');
ylabel(ax1, 'L BP angle (deg)');
title(ax1, sprintf('Inspiration-triggered L angle\nBP %g-%g Hz, n=%d events, %d sessions', ...
    WHISK_BP(1), WHISK_BP(2), size(L_evt,1), nSess), 'Interpreter','none');
axis(ax1, 'square');
box(ax1, 'off');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
plot_triggered_mean_sem(ax2, tRel, R_evt, [0.15 0.15 0.15]);
xline(ax2, 0, 'k--', 'inspiration');
yline(ax2, 0, 'k:');
xlabel(ax2, 'time from inspiration event (s)');
ylabel(ax2, 'R BP angle (deg)');
title(ax2, sprintf('Inspiration-triggered R angle\nBP %g-%g Hz, n=%d events, %d sessions', ...
    WHISK_BP(1), WHISK_BP(2), size(R_evt,1), nSess), 'Interpreter','none');
axis(ax2, 'square');
box(ax2, 'off');

sgtitle(sprintf('Triggered average: whisk BP %g-%g Hz, canonical breath BP %g-%g Hz, no breath inversion', ...
    WHISK_BP(1), WHISK_BP(2), BREATH_BP(1), BREATH_BP(2)));

% ============================= HELPERS =============================

function plot_triggered_mean_sem(ax, tRel, X, col)

    if isempty(X)
        title(ax, 'no events');
        return;
    end

    mu = mean(X, 1, 'omitnan');
    n  = sum(isfinite(X), 1);
    se = std(X, 0, 1, 'omitnan') ./ sqrt(n);

    fill(ax, [tRel; flipud(tRel)], ...
        [(mu+se)'; flipud((mu-se)')], ...
        col, 'FaceAlpha', 0.22, 'EdgeColor', 'none');

    plot(ax, tRel, mu, 'Color', col, 'LineWidth', 2);
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