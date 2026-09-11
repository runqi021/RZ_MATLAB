% whisk_LR_integral_vs_breath_freq_RZ.m
%
% Plot int(L-R) versus instantaneous breathing frequency.
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   BP 5-50 Hz, zero-phase filtfilt
%   whisk epoch mask from Hilbert envelope of LEFT BP whisker
%
% Breathing:
%   load Bs.breath
%   interpolate to whisker timebase
%   BP 1-20 Hz
%   instantaneous breathing frequency from Hilbert phase:
%       f_inst = abs(d phase / dt) / (2*pi)
%
% Output:
%   Figure 1:
%       all points: int(L-R) vs instantaneous breath frequency
%       binned mean +/- SEM
%
%   Figure 2:
%       same, but session mean curves overlaid

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz
BREATH_BP = [1 20];   % Hz, for instantaneous phase/frequency

ENV_THR    = 5;       % deg, hard envelope threshold
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

ASYM_INT_WIN = 0.050; % s, local integral window

TAKE_ABS_LR = true;   % true: int(|L-R|), false: int(L-R)

FREQ_SMOOTH = 0.050;  % s, smooth instantaneous breathing frequency
FREQ_LIM    = [0.5 20];   % Hz, discard absurd values

FREQ_BIN_W  = 0.5;    % Hz
FREQ_BIN_EDGES = 0.5:FREQ_BIN_W:20;

MAX_POINTS_PER_SESSION = 3000;   % for scatter downsample only

MIN_SAMPLES = 50;

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

winSamp = max(1, round(ASYM_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

freqSmoothSamp = max(1, round(FREQ_SMOOTH * fpsW));

% ============================ ACCUMULATORS ============================

F_all = [];
Y_all = [];
SessID_all = [];

sessCurve = struct( ...
    'animal', {}, ...
    'trial', {}, ...
    'fBin', {}, ...
    'yMean', {}, ...
    'n', {} );

nsess = 0;

% ============================ SESSION LOOP ============================

ad = dir(char(dataRoot));

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

        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', animal, kk))))
            continue;
        end

        % ============================ WHISKER ============================

        try
            M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
        catch ME
            warning('whisk load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

        t = (0:numel(La)-1)' / fpsW;
        N = numel(t);

        La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
        Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');

        La0 = La0 - mean(La0, 'omitnan');
        Ra0 = Ra0 - mean(Ra0, 'omitnan');

        xL = filtfilt(b_wbp, a_wbp, La0);   % L BP 5-50 Hz
        xR = filtfilt(b_wbp, a_wbp, Ra0);   % R BP 5-50 Hz

        % ======================= WHISK EPOCH MASK ========================

        env_raw = abs(hilbert(xL));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));

        is_whisk = env > ENV_THR;
        ep = bool_to_epochs(is_whisk, t, MIN_DUR, MERGE_GAP);

        m = false(N,1);
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        % ============================ BREATH =============================

        brw_bp1_20 = nan(N,1);
        f_breath   = nan(N,1);

        try
            Pn = thermal_resolve_paths( ...
                pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);

            if isfile(Pn.breath)

                Bs = load(Pn.breath);

                br = Bs.breath(:);
                fb = double(Bs.fps);
                tBr = (0:numel(br)-1)' / fb;

                brw = interp1(tBr, br, t, 'linear', NaN);

                if nnz(isfinite(brw)) < MIN_SAMPLES
                    continue;
                end

                brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');

                % Keep same direction as current script.
                brw_bp1_20 = filtfilt(b_br, a_br, brw);

                % Instantaneous breathing frequency from Hilbert phase.
                ph = unwrap(angle(hilbert(brw_bp1_20)));
                f_breath = abs(gradient(ph) * fpsW / (2*pi));

                % Smooth and clean.
                f_breath = movmedian(f_breath, freqSmoothSamp, 'omitnan');
                f_breath = movmean(f_breath, freqSmoothSamp, 'omitnan');

            end

        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        % ============================ INTEGRAL ===========================

        if TAKE_ABS_LR
            lr_inst = abs(xL - xR);
        else
            lr_inst = xL - xR;
        end

        int_LR = movsum(lr_inst, winSamp, 'omitnan') / fpsW;

        % ======================== COLLECT POINTS =========================

        mask = m & ...
               isfinite(int_LR) & ...
               isfinite(f_breath) & ...
               f_breath >= FREQ_LIM(1) & ...
               f_breath <= FREQ_LIM(2);

        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        f_this = f_breath(mask);
        y_this = int_LR(mask);

        % downsample scatter only, keep time order
        idx = find(mask);
        stride = max(1, ceil(numel(idx) / MAX_POINTS_PER_SESSION));
        idx_ds = idx(1:stride:end);

        nsess = nsess + 1;

        F_all = [F_all; f_breath(idx_ds)]; %#ok<SAGROW>
        Y_all = [Y_all; int_LR(idx_ds)];   %#ok<SAGROW>
        SessID_all = [SessID_all; nsess * ones(numel(idx_ds),1)]; %#ok<SAGROW>

        % session binned curve
        [fb, ym, nn] = bin_mean_curve(f_this, y_this, FREQ_BIN_EDGES);

        sessCurve(nsess).animal = animal;
        sessCurve(nsess).trial  = kk;
        sessCurve(nsess).fBin   = fb;
        sessCurve(nsess).yMean  = ym;
        sessCurve(nsess).n      = nn;

    end
end

fprintf('%d sessions contributed\n', nsess);

% ============================ POPULATION BIN ============================

[fBin, yMean, ySEM, yN] = bin_mean_sem(F_all, Y_all, FREQ_BIN_EDGES);

if TAKE_ABS_LR
    yLabelText = sprintf('local integral of |L-R|, deg*s, win=%.0f ms', 1000*ASYM_INT_WIN);
    ttlText = 'int(|L-R|) vs instantaneous breathing frequency';
else
    yLabelText = sprintf('local integral of L-R, deg*s, win=%.0f ms', 1000*ASYM_INT_WIN);
    ttlText = 'int(L-R) vs instantaneous breathing frequency';
end

% ============================== FIGURE 1 ==============================
% All pooled points plus population binned mean +/- SEM

figure('Color','w','Position',[140 140 720 560]);

ax1 = axes; hold(ax1,'on'); grid(ax1,'on');

scatter(ax1, F_all, Y_all, 6, ...
    'MarkerFaceColor', [0.6 0.6 0.6], ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.18);

good = isfinite(fBin) & isfinite(yMean) & isfinite(ySEM) & yN >= 3;

fill(ax1, [fBin(good); flipud(fBin(good))], ...
    [yMean(good)+ySEM(good); flipud(yMean(good)-ySEM(good))], ...
    [0.2 0.2 0.2], ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');

plot(ax1, fBin(good), yMean(good), 'k-', 'LineWidth', 2);

xlabel(ax1, 'instantaneous breathing frequency (Hz)');
ylabel(ax1, yLabelText);

title(ax1, sprintf('%s\npooled points + binned mean +/- SEM, n=%d sessions', ...
    ttlText, nsess), ...
    'Interpreter','none');

xlim(ax1, FREQ_LIM);
box(ax1,'off');

% ============================== FIGURE 2 ==============================
% Each session binned curve overlaid

figure('Color','w','Position',[180 180 720 560]);

ax2 = axes; hold(ax2,'on'); grid(ax2,'on');

C = lines(max(nsess,1));

for si = 1:nsess

    f = sessCurve(si).fBin;
    y = sessCurve(si).yMean;
    n = sessCurve(si).n;

    good = isfinite(f) & isfinite(y) & n >= 3;

    if ~any(good)
        continue;
    end

    ci = mod(si-1, size(C,1)) + 1;

    plot(ax2, f(good), y(good), '-', ...
        'Color', C(ci,:), ...
        'LineWidth', 1.0);
end

plot(ax2, fBin(good), yMean(good), 'k-', 'LineWidth', 3);

xlabel(ax2, 'instantaneous breathing frequency (Hz)');
ylabel(ax2, yLabelText);

title(ax2, sprintf('%s\nsession curves overlaid, black = population mean', ttlText), ...
    'Interpreter','none');

xlim(ax2, FREQ_LIM);
box(ax2,'off');

% ============================= HELPERS =============================

function [fb, ym, ys, nn] = bin_mean_sem(x, y, edges)

    x = x(:);
    y = y(:);

    nb = numel(edges) - 1;

    fb = edges(1:end-1)' + diff(edges(:))/2;
    ym = nan(nb,1);
    ys = nan(nb,1);
    nn = zeros(nb,1);

    for i = 1:nb
        m = x >= edges(i) & x < edges(i+1) & isfinite(x) & isfinite(y);
        nn(i) = nnz(m);

        if nn(i) >= 1
            ym(i) = mean(y(m), 'omitnan');
            ys(i) = std(y(m), 0, 'omitnan') / sqrt(nn(i));
        end
    end
end

function [fb, ym, nn] = bin_mean_curve(x, y, edges)

    x = x(:);
    y = y(:);

    nb = numel(edges) - 1;

    fb = edges(1:end-1)' + diff(edges(:))/2;
    ym = nan(nb,1);
    nn = zeros(nb,1);

    for i = 1:nb
        m = x >= edges(i) & x < edges(i+1) & isfinite(x) & isfinite(y);
        nn(i) = nnz(m);

        if nn(i) >= 1
            ym(i) = mean(y(m), 'omitnan');
        end
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