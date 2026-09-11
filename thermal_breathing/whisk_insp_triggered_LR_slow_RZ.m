% whisk_insp_triggered_LR_slow_RZ.m
%
% Goal:
%   Detect inspiration events from thermal breath trace,
%   then trigger-average slow L/R whisker angle around those events.
%
% Output:
%   Figure with two plots:
%       1) inspiration-triggered L angle
%       2) inspiration-triggered R angle
%
% Whisker:
%   raw angle -> fill -> LP < LP_HZ -> demean
%
% Breathing:
%   load *_breath_raw.mat
%   breath_raw = raw temp, no LP-sub, no invert
%   BP 1-20 Hz for event detection
%   inspiration event = findpeaks(-br_bp)
%
% Event time:
%   t = 0 is detected inspiration event / cooling trough.
%
% Animal 5840027 excluded by default.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

LP_HZ = 1;              % whisker slow set-point cutoff
BREATH_BP = [1 20];     % breath band for inspiration event detection

PRE_S  = 0.30;          % seconds before inspiration event
POST_S = 0.30;          % seconds after inspiration event

MIN_INSP_DIST_S = 0.06; % minimum distance between detected inspiration events
PROM_MAD = 0.5;         % adaptive peak prominence = PROM_MAD * MAD(event trace)

EXCLUDE = "5840027";

% If raw thermal breath gets colder during inspiration, keep -1.
% If your breath trace is already inverted so inspiration is positive, use +1.
INSP_POLARITY = -1;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

[b_wslow, a_wslow] = butter(3, LP_HZ/(fpsW/2), 'low');

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

    La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
    Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');

    % slow set-point angle
    sL = filtfilt(b_wslow, a_wslow, La0);
    sR = filtfilt(b_wslow, a_wslow, Ra0);

    % session demean, same spirit as your scatter script
    sL = sL - mean(sL, 'omitnan');
    sR = sR - mean(sR, 'omitnan');

    % ============================ BREATH =============================

    try
        noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun));
        Pn = thermal_resolve_paths(noseCsv, dataRoot);

        rawMat = resolve_breath_raw_mat(Pn, noseCsv);

        if strlength(rawMat) == 0 || ~isfile(rawMat)
            warning('raw breath file not found %s n%d', animal, kRun);
            continue;
        end

        Bs = load(rawMat);

        if ~isfield(Bs, 'breath_raw') || ~isfield(Bs, 'fps')
            warning('breath_raw or fps missing in %s', rawMat);
            continue;
        end

        br = Bs.breath_raw(:);    % raw temp, no LP-sub, no invert
        fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb;

        br = fillmissing(br, 'linear', 'EndValues', 'nearest');
        br = br - mean(br, 'omitnan');

        [b_br, a_br] = butter(4, BREATH_BP/(fb/2), 'bandpass');
        br_bp = filtfilt(b_br, a_br, br);

        % Inspiration detection trace.
        % Default: raw thermal inspiration is cooling, so -br_bp makes inspiration positive.
        event_trace = INSP_POLARITY * br_bp;

        prom = PROM_MAD * mad(event_trace, 1);
        if ~isfinite(prom) || prom <= 0
            prom = std(event_trace, 'omitnan') * 0.5;
        end

        [~, insp_t] = findpeaks(event_trace, tB, ...
            'MinPeakDistance', MIN_INSP_DIST_S, ...
            'MinPeakProminence', prom);

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

        segL = interp1(tW, sL, tt, 'linear', NaN);
        segR = interp1(tW, sR, tt, 'linear', NaN);

        if mean(isfinite(segL)) < 0.95 || mean(isfinite(segR)) < 0.95
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
plot_triggered_mean_sem(ax1, tRel, L_evt, [0.1 0.1 0.1]);
xline(ax1, 0, 'k--', 'inspiration');
xlabel(ax1, 'time from inspiration event (s)');
ylabel(ax1, 'L slow angle (deg)');
title(ax1, sprintf('Inspiration-triggered L angle\nLP<%g Hz, n=%d events, %d sessions', ...
    LP_HZ, size(L_evt,1), nSess), 'Interpreter','none');
axis(ax1, 'square');
box(ax1, 'off');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
plot_triggered_mean_sem(ax2, tRel, R_evt, [0.1 0.1 0.1]);
xline(ax2, 0, 'k--', 'inspiration');
xlabel(ax2, 'time from inspiration event (s)');
ylabel(ax2, 'R slow angle (deg)');
title(ax2, sprintf('Inspiration-triggered R angle\nLP<%g Hz, n=%d events, %d sessions', ...
    LP_HZ, size(R_evt,1), nSess), 'Interpreter','none');
axis(ax2, 'square');
box(ax2, 'off');

sgtitle(sprintf('Triggered average around inspiration events: breath BP %g-%g Hz, polarity=%+d', ...
    BREATH_BP(1), BREATH_BP(2), INSP_POLARITY));

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

    yline(ax, 0, 'k:');
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

function rawMat = resolve_breath_raw_mat(Pn, noseCsv)
% Find <stem>_breath_raw.mat next to the thermal .ats.

    rawMat = "";

    candidates = strings(0,1);

    directFields = {'breath_raw', 'breathRaw', 'raw_breath', 'rawBreath'};

    for i = 1:numel(directFields)
        fn = directFields{i};
        if isfield(Pn, fn)
            try
                candidates(end+1,1) = string(Pn.(fn)); %#ok<AGROW>
            catch
            end
        end
    end

    if isfield(Pn, 'ats')
        atsPath = string(Pn.ats);

        if isfile(atsPath)
            [p, stem] = fileparts(atsPath);
            candidates(end+1,1) = string(fullfile(p, stem + "_breath_raw.mat")); %#ok<AGROW>
        end
    end

    if isfield(Pn, 'breath')
        oldBreath = string(Pn.breath);
        [p, ~] = fileparts(oldBreath);

        if isfolder(p)
            d = dir(fullfile(p, '*_breath_raw.mat'));
            for k = 1:numel(d)
                candidates(end+1,1) = string(fullfile(d(k).folder, d(k).name)); %#ok<AGROW>
            end
        end
    end

    if nargin >= 2
        noseCsv = string(noseCsv);
        [p, ~] = fileparts(noseCsv);

        if isfolder(p)
            d = dir(fullfile(p, '*_breath_raw.mat'));
            for k = 1:numel(d)
                candidates(end+1,1) = string(fullfile(d(k).folder, d(k).name)); %#ok<AGROW>
            end
        end
    end

    candidates = unique(candidates, 'stable');

    for i = 1:numel(candidates)
        if strlength(candidates(i)) > 0 && isfile(candidates(i))
            rawMat = candidates(i);
            return;
        end
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