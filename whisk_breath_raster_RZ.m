function whisk_breath_raster_RZ()

% whisk_breath_raster_RZ
%
% Breath-cycle-triggered raster of L, R, and pooled L/R whisker PROTRACTION
% events relative to inspiration event.
%
% x = time from inspiration event (s)
% y = breath number, ordered by breath duration, short -> long
%
% teal = large protraction events
% red  = next inspiration event = breath duration
%
% Three panels:
%   1) LEFT whisker protraction events
%   2) RIGHT whisker protraction events
%   3) pooled LEFT + RIGHT protraction events
%
% Filtering/detection:
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   BP 5-50 Hz, zero-phase filtfilt
%   whisk epoch mask from Hilbert envelope of LEFT BP whisker
%   hard envelope threshold = 5 deg
%
% Breathing:
%   _breath.mat / Pn.breath
%   Bs.breath is already LP<1 Hz subtracted AND already inverted
%   inhale is positive/upward
%   DO NOT invert again
%   BP 1-20 Hz, zero-phase filtfilt
%
% Inspiration onset:
%   trough (local min) of the inhale-up BP breath = most-exhaled point just
%   before inhalation begins -> findpeaks(-brw_bp1_20)
%   (breath is already inhale-up in the saved .mat; inversion sets DIRECTION,
%    not onset-vs-peak, so onset = trough, not the inhalation peak)
%
% Protraction onset:
%   trough of BP whisker angle (protraction = increasing angle, so the trough
%   is the most-retracted point where the forward sweep begins), inside valid
%   whisking epoch; prominence = WHISK_PROM * std(BP whisk angle within mask)

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 60];   % Hz
BREATH_BP = [2 15];   % Hz

ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0;     % s

WHISK_PROM       = 1;   % protraction trough prominence = WHISK_PROM * std(BP whisk)
MIN_WHISK_DIST_S = 0.05;  % s, minimum distance between protraction events

MIN_INSP_DIST_S = 0.06;   % s, minimum distance between inspiration events
BREATH_PROM_MAD = 1.5;    % inspiration peak prominence = BREATH_PROM_MAD * MAD(breath)

PRE_S = 0;             % s before t0 (current inspiration onset) shown on raster
XMAX  = 1.2;            % raster x-limit max; [] = auto (through 3rd inspiration onset)

EXCLUDE = "5840027";
EXCLUDE_SESS = {{"5916297",3}, {"5916297",4}, {"5916300",4}};  % {animal,run} drops

% If protraction is a trough in your BP angle, keep -1.
% If protraction is a peak in your BP angle, change to +1.
PROTRACTION_POLARITY = -1;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

% ============================ SESSION LIST ============================

ad = dir(char(dataRoot));
sess = {};

for ai = 1:numel(ad)

    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit'))
        continue;
    end

    if any(strcmp(ad(ai).name, EXCLUDE))
        continue;
    end

    animal = ad(ai).name;

    rr = dir(fullfile(char(dataRoot), animal, 'cam1_*'));
    [~,o] = sort({rr.name});
    rr = rr(o);

    for kk = 1:numel(rr)
        if is_excluded_sess(animal, kk, EXCLUDE_SESS)
            continue;
        end
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', animal, kk))))
            sess{end+1} = {animal, kk}; %#ok<AGROW>
        end
    end
end

assert(~isempty(sess), 'no sessions with whisk csv');

% ============================ ACCUMULATORS ============================

dur  = [];       % time of onset 1 (i+1) relative to t0 = breath duration (sort key)
dur3 = [];       % time of onset 3 (i+3) relative to t0 = whisk-window end

protL_byBreath = {};
protR_byBreath = {};
protP_byBreath = {};  % pooled L/R

nSessUsed = 0;

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1};
    kk     = sess{e}{2};

    % ============================ WHISKER ============================

    try
        whiskCsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk));
        M = dlc_gate_interp(whiskCsv, 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kk, ME.message);
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
        warning('no valid whisk epoch %s n%d', animal, kk);
        continue;
    end

    % ======================= PROTRACTION EVENTS =======================

    tProtL = detect_large_protractions( ...
        xL, tW, m, MIN_WHISK_DIST_S, WHISK_PROM, PROTRACTION_POLARITY);

    tProtR = detect_large_protractions( ...
        xR, tW, m, MIN_WHISK_DIST_S, WHISK_PROM, PROTRACTION_POLARITY);

    tProtPool = sort([tProtL(:); tProtR(:)]);

    if isempty(tProtL) && isempty(tProtR)
        warning('no protraction events detected %s n%d', animal, kk);
        continue;
    end

    % ============================ BREATH =============================

    try
        noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk));
        Pn = thermal_resolve_paths(noseCsv, dataRoot);

        if ~isfile(Pn.breath)
            warning('canonical LP-sub breath file not found %s n%d', animal, kk);
            continue;
        end

        Bs = load(Pn.breath);

        br = Bs.breath(:);      % already LP<1 Hz subtracted
                                % already inverted: inhale = positive
        fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb;

        brw = interp1(tB, br, tW, 'linear', NaN);

        if nnz(isfinite(brw)) < 50
            warning('too few finite breath samples %s n%d', animal, kk);
            continue;
        end

        brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');

        % IMPORTANT: DO NOT re-invert. Saved breath is already inhale-up
        % (inhale = positive); inversion sets DIRECTION, not onset-vs-peak.
        brw_bp1_20 = filtfilt(b_br, a_br, brw);

        % Inspiration ONSET = trough (local min) of the inhale-up breath, i.e.
        % the most-exhaled point just before inhalation begins. Same logic as
        % protraction onset (trough of the rising signal): detect peaks of -x.
        event_trace = -brw_bp1_20;

        prom = BREATH_PROM_MAD * mad(brw_bp1_20, 1);
        if ~isfinite(prom) || prom <= 0
            prom = std(brw_bp1_20, 'omitnan') * 0.5;
        end

        [~, locs] = findpeaks(event_trace, ...
            'MinPeakDistance', round(MIN_INSP_DIST_S * fpsW), ...
            'MinPeakProminence', prom);

        tInsp = tW(locs);

    catch ME
        warning('breath load/detect failed %s n%d: %s', animal, kk, ME.message);
        continue;
    end

    if numel(tInsp) < 4
        warning('too few inspiration events %s n%d', animal, kk);
        continue;
    end

    % ============================ PER BREATH ===========================
    % onset numbering with the trigger = onset 0:
    %   onset 0 = t0 (trigger), onset 1 = i+1 (sort key), ..., onset 3 = i+3

    nBreathSess = 0;

    for i = 1:numel(tInsp)-3

        d  = tInsp(i+1) - tInsp(i);   % onset 1 (sort key = breath duration)
        d3 = tInsp(i+3) - tInsp(i);   % onset 3 (whisk window end)

        lo = tInsp(i) - PRE_S;
        hi = tInsp(i+3);             % collect whisk through onset 3

        relL = tProtL(tProtL >= lo & tProtL <= hi) - tInsp(i);
        relR = tProtR(tProtR >= lo & tProtR <= hi) - tInsp(i);
        relP = tProtPool(tProtPool >= lo & tProtPool <= hi) - tInsp(i);

        dur(end+1,1)  = d;  %#ok<AGROW>
        dur3(end+1,1) = d3; %#ok<AGROW>

        protL_byBreath{end+1,1} = relL(:)'; %#ok<AGROW>
        protR_byBreath{end+1,1} = relR(:)'; %#ok<AGROW>
        protP_byBreath{end+1,1} = relP(:)'; %#ok<AGROW>

        nBreathSess = nBreathSess + 1;
    end

    if nBreathSess > 0
        nSessUsed = nSessUsed + 1;

        fprintf('%s n%d: %d breaths, L=%d, R=%d, pooled=%d protractions\n', ...
            animal, kk, nBreathSess, numel(tProtL), numel(tProtR), numel(tProtPool));
    end
end

assert(~isempty(dur), 'no breaths found');

% ============================ SORT + SCATTER DATA ======================

[durs, ord] = sort(dur);
dur3s = dur3(ord);

protL_byBreath = protL_byBreath(ord);
protR_byBreath = protR_byBreath(ord);
protP_byBreath = protP_byBreath(ord);

N = numel(durs);

% x-axis upper limit: cover through the third inspiration onset.
% XMAX==[] -> auto (99th pct of third-onset time); else use the fixed value.
if isempty(XMAX)
    xhi = prctile(dur3s, 99);
else
    xhi = XMAX;
end

[XL, YL] = breath_cells_to_xy(protL_byBreath);
[XR, YR] = breath_cells_to_xy(protR_byBreath);
[XP, YP] = breath_cells_to_xy(protP_byBreath);

% ============================== PLOT ==================================

figure('Color','w','Position',[80 80 1400 650]);

tl = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile; hold(ax1,'on');
plot_raster_panel(ax1, XL, YL, durs, N, PRE_S, xhi, ...
    'LEFT whisker');

ax2 = nexttile; hold(ax2,'on');
plot_raster_panel(ax2, XR, YR, durs, N, PRE_S, xhi, ...
    'RIGHT whisker');

ax3 = nexttile; hold(ax3,'on');
plot_raster_panel(ax3, XP, YP, durs, N, PRE_S, xhi, ...
    'POOLED L/R whisker');

title(tl, sprintf(['Large protraction raster relative to inspiration onset (t0)\n' ...
                   'teal = protraction onset, red = next inspiration onset (sort key); ' ...
                   'whisker shown through 3rd onset; %d breaths, %d sessions'], ...
                   N, nSessUsed), ...
      'Interpreter','none');

end

% ============================= HELPERS =============================

function [X, Y] = breath_cells_to_xy(C)

    X = [];
    Y = [];

    for r = 1:numel(C)
        X = [X C{r}]; %#ok<AGROW>
        Y = [Y r * ones(1, numel(C{r}))]; %#ok<AGROW>
    end
end

function plot_raster_panel(ax, X, Y, durs, N, PRE_S, XMAX, ttl)

    plot(ax, X, Y, '.', ...
        'Color', [0.10 0.55 0.55], ...
        'MarkerSize', 3);

    % next inspiration onset (i+1 = breath duration, the sort key)
    plot(ax, durs, 1:N, '.', ...
        'Color', [0.85 0.10 0.10], ...
        'MarkerSize', 3);

    xlim(ax, [-PRE_S XMAX]);
    ylim(ax, [1 N]);
    set(ax, 'YDir', 'reverse');

    xlabel(ax, 'Time from inspiration onset (s)');
    ylabel(ax, 'Breath number, ordered by duration');
    title(ax, ttl, 'Interpreter','none');

    box(ax, 'on');
    grid(ax, 'on');
end

function tProt = detect_large_protractions(x, t, mask, minDistS, whiskProm, polarity)
% Detect large protraction events from BP whisker angle.
%
% polarity = -1 means protraction is a trough, so detect peaks of -x.
% polarity = +1 means protraction is a peak, so detect peaks of x.
%
% Detection is restricted to the valid whisking epoch mask.
% Peak prominence = whiskProm * std(x within mask).

    x = x(:);
    t = t(:);
    mask = logical(mask(:));

    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        error('bad time vector');
    end

    minDistSamples = max(1, round(minDistS / dt));

    good = mask & isfinite(x);

    prom = whiskProm * std(x(good), 'omitnan');

    if ~isfinite(prom) || prom <= 0
        prom = 0.5 * std(x(isfinite(x)), 'omitnan');
    end

    if ~isfinite(prom) || prom <= 0
        tProt = zeros(0,1);
        return;
    end

    detectTrace = polarity * x;

    [~, idx] = findpeaks(detectTrace, ...
        'MinPeakProminence', prom, ...
        'MinPeakDistance', minDistSamples);

    idx = idx(mask(idx));

    tProt = t(idx);
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

function tf = is_excluded_sess(animal, kRun, excl)

    tf = false;
    for q = 1:numel(excl)
        if strcmp(char(excl{q}{1}), char(animal)) && double(excl{q}{2}) == kRun
            tf = true;
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