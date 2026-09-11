% whisk_LR_integral_xcorr_filtered_RZ.m
%
% Population cross-correlation:
%   local int(L+R) vs inverted BP breathing
%   local int(L-R) vs inverted BP breathing
%
% Whisker:
%   demean angle
%   LP < 5 Hz
%   BP 5-50 Hz
%   detection from hard 5 deg Hilbert envelope of BP whisking
%
% Breathing:
%   input breath is already LP<1 Hz subtracted
%   invert temperature signal
%   BP 1-20 Hz
%
% Lag convention:
%   corr( whisk_signal(t+lag), breath(t) )
%   lag > 0: breath leads whisker
%   lag < 0: whisker leads breath

clear; clc;

take_abs = false;

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_LP = 5;        % Hz
WHISK_BP = [5 50];   % Hz
BREATH_BP = [1 20];  % Hz

ENV_THR    = 5;      % hard envelope threshold, deg
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;      % s
MERGE_GAP  = 0.2;    % s

ASYM_INT_WIN = 0.050;  % s
XC_MAXLAG    = 0.5;    % s
MIN_SAMPLES  = 50;

EXCLUDE = "5840027";
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wlp,a_wlp] = butter(ord, WHISK_LP/(fpsW/2), 'low');
[b_wbp,a_wbp] = butter(ord, WHISK_BP/(fpsW/2), 'bandpass');
[b_br, a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

nLag   = round(XC_MAXLAG * fpsW);
lags_t = (-nLag:nLag)' / fpsW;

winSamp = max(1, round(ASYM_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

xc_sum = zeros(0, numel(lags_t));
xc_lr  = zeros(0, numel(lags_t));
ac_br  = zeros(0, numel(lags_t));
nsess  = 0;

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
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8)))));

        t = (0:numel(La)-1)' / fpsW;
        N = numel(t);

        La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
        Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');

        La0 = La0 - mean(La0, 'omitnan');
        Ra0 = Ra0 - mean(Ra0, 'omitnan');

        La_lp5    = filtfilt(b_wlp, a_wlp, La0);
        Ra_lp5    = filtfilt(b_wlp, a_wlp, Ra0);
        La_bp5_50 = filtfilt(b_wbp, a_wbp, La0);
        Ra_bp5_50 = filtfilt(b_wbp, a_wbp, Ra0);

        % Actual whisker time series used for xcorr.
        % Use BP whisking. If you want slow+fast angle instead, use:
        % xL = La_lp5 + La_bp5_50; xR = Ra_lp5 + Ra_bp5_50;
        xL = La_bp5_50;
        xR = Ra_bp5_50;

        % ======================= WHISK EPOCH MASK ========================
        env_raw = abs(hilbert(La_bp5_50));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));

        is_whisk = env > ENV_THR;
        ep = bool_to_epochs(is_whisk, t, MIN_DUR, MERGE_GAP);

        m = false(N,1);
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        % ============================ BREATH =============================
        brw_bp1_20 = nan(N,1);

        try
            Pn = thermal_resolve_paths( ...
                pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);

            if isfile(Pn.breath)
                Bs = load(Pn.breath);

                br = Bs.breath(:);       % already LP<1 Hz subtracted
                fb = double(Bs.fps);
                tBr = (0:numel(br)-1)' / fb;

                brw = interp1(tBr, br, t, 'linear', NaN);
                brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');

                brw_inv = -brw;          % invert temperature signal
                brw_bp1_20 = filtfilt(b_br, a_br, brw_inv);
            end

        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
        end

        % ============================ INTEGRAL ===========================
        if take_abs
            lr_inst  = abs(xL - xR);
            sum_inst = abs(xL + xR);
        else
            lr_inst  = xL - xR;
            sum_inst = xL + xR;
        end

        int_LR  = movsum(lr_inst,  winSamp, 'omitnan') / fpsW;
        int_SUM = movsum(sum_inst, winSamp, 'omitnan') / fpsW;

        % ============================= XCORR =============================
        mask = m & isfinite(int_LR) & isfinite(int_SUM) & isfinite(brw_bp1_20);

        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        cc_sum = lagged_corr(int_SUM, brw_bp1_20, mask, nLag, MIN_SAMPLES);
        cc_lr  = lagged_corr(int_LR,  brw_bp1_20, mask, nLag, MIN_SAMPLES);
        acc_br = lagged_corr(brw_bp1_20, brw_bp1_20, mask, nLag, MIN_SAMPLES);

        if any(isfinite(cc_lr))
            xc_sum(end+1,:) = cc_sum'; %#ok<SAGROW>
            xc_lr(end+1,:)  = cc_lr';  %#ok<SAGROW>
            ac_br(end+1,:)  = acc_br'; %#ok<SAGROW>
            nsess = nsess + 1;
        end

    end
end

fprintf('%d sessions contributed\n', nsess);

% =============================== FIGURE ===============================
col = [0.2 0.2 0.2];

figure('Color','w','Position',[120 200 1240 380]);

ax0 = subplot(1,3,1); hold on; grid on;
plot_xcorr_panel(ax0, xc_sum, lags_t, col, ...
    'int(L+R) vs breath', ...
    'lag (s), corr[int(L+R)(t+lag), breath(t)]', 'corr', true);

ax1 = subplot(1,3,2); hold on; grid on;
plot_xcorr_panel(ax1, xc_lr, lags_t, col, ...
    'int(L-R) vs breath', ...
    'lag (s), corr[int(L-R)(t+lag), breath(t)]', 'corr', true);

ax2 = subplot(1,3,3); hold on; grid on;
plot_xcorr_panel(ax2, ac_br, lags_t, col, ...
    'breath autocorr', ...
    'lag (s)', 'corr', false);
ylim(ax2, [-1 1]);

% ============================= HELPERS =============================

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

function cc = lagged_corr(x, y, mask, nLag, minSamp)

    x = x(:);
    y = y(:);
    mask = logical(mask(:));

    N = numel(x);
    lags = -nLag:nLag;
    cc = nan(numel(lags),1);

    for li = 1:numel(lags)

        k = lags(li);
        ti = max(1,1-k) : min(N,N-k);

        v = mask(ti) & isfinite(x(ti+k)) & isfinite(y(ti));

        if nnz(v) >= minSamp
            cc(li) = corr(x(ti(v)+k), y(ti(v)));
        end
    end
end

function plot_xcorr_panel(ax, mat, lags_t, col, ttl, xlab, ylab, showPeak)

    if isempty(mat)
        title(ax, sprintf('%s\nn=0', ttl), 'Interpreter','none');
        return;
    end

    n  = size(mat,1);
    mu = mean(mat, 1, 'omitnan');
    se = std(mat, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(mat),1));

    fill(ax, [lags_t; flipud(lags_t)], ...
        [(mu+se)'; flipud((mu-se)')], ...
        col, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    plot(ax, lags_t, mu, 'Color', col, 'LineWidth', 2);
    xline(ax, 0, 'k--');

    if showPeak
        [pkv,pki] = max(abs(mu));
        title(ax, sprintf('%s\npeak |r|=%.3f @ %.0f ms, n=%d', ...
            ttl, pkv, 1000*lags_t(pki), n), 'Interpreter','none');
    else
        title(ax, sprintf('%s\nn=%d', ttl, n), 'Interpreter','none');
    end

    xlabel(ax, xlab);
    ylabel(ax, ylab);
    xlim(ax, [lags_t(1), lags_t(end)]);
    axis(ax, 'square');
    box(ax, 'off');
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