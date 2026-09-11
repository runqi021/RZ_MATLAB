% whisk_LR_sum_breath_xcorr_quick_RZ.m
%
% Quick population xcorr:
%   sum(L+R) vs breath
%
% Same preprocessing as working script:
%   Whisker: La/Ra -> fill NaNs -> demean -> BP 5-50 Hz
%   Epochs: Hilbert envelope of LEFT BP whisker, threshold 5 deg
%   Breath: brw already LP<1 Hz subtracted -> invert -> BP 1-20 Hz
%
% NO integration.
% NO shuffle.
%
% Lag convention:
%   corr( sum(L+R)(t+lag), breath(t) )
%   lag > 0: breath leads sum(L+R)
%   lag < 0: sum(L+R) leads breath

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

XC_MAXLAG    = 0.5;   % s
MIN_SAMPLES  = 50;

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

nLag   = round(XC_MAXLAG * fpsW);
lags_t = (-nLag:nLag)' / fpsW;
nL     = numel(lags_t);

xc_sum = zeros(0,nL);
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

        try
            Pn = thermal_resolve_paths( ...
                pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);

            if isfile(Pn.breath)

                Bs = load(Pn.breath);

                br = Bs.breath(:);       % already LP<1 Hz subtracted
                fb = double(Bs.fps);
                tBr = (0:numel(br)-1)' / fb;

                brw = interp1(tBr, br, t, 'linear', NaN);

                if nnz(isfinite(brw)) < MIN_SAMPLES
                    continue;
                end

                brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');

                brw_inv = -brw;          % invert temperature breathing signal
                brw_bp1_20 = filtfilt(b_br, a_br, brw_inv);  % BP 1-20 Hz
            end

        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        % ========================= SUM, NO INTEGRAL ======================

        sum_inst = xL + xR;

        % ============================= XCORR =============================

        mask = m & isfinite(sum_inst) & isfinite(brw_bp1_20);

        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        cc_sum = lagged_corr_fast(sum_inst, brw_bp1_20, mask, nLag, MIN_SAMPLES);

        if any(isfinite(cc_sum))
            xc_sum(end+1,:) = cc_sum'; %#ok<SAGROW>
            nsess = nsess + 1;
        end
    end
end

fprintf('%d sessions contributed\n', nsess);

% ============================== PLOT ==================================

figure('Color','w','Position',[220 220 520 480]);
ax = axes; hold(ax,'on'); grid(ax,'on');

plot_xcorr_panel(ax, xc_sum, lags_t, [0.2 0.2 0.2], ...
    'sum(L+R) vs breath', ...
    'lag (s), corr[sum(L+R)(t+lag), breath(t)]', ...
    'corr');

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

function cc = lagged_corr_fast(x, y, mask, nLag, minSamp)
% Pearson corr( x(t+lag), y(t) ) for lag = -nLag:nLag.
% lag > 0 means y/ref leads x/signal.

    x = x(:);
    y = y(:);
    mask = logical(mask(:));

    N = numel(x);
    lags = -nLag:nLag;
    cc = nan(numel(lags),1);

    for li = 1:numel(lags)

        k = lags(li);

        if k >= 0
            ix = (1+k):N;      % x(t+k)
            iy = 1:(N-k);      % y(t)
        else
            ix = 1:(N+k);      % x(t+k)
            iy = (1-k):N;      % y(t)
        end

        v = mask(iy) & isfinite(x(ix)) & isfinite(y(iy));

        if nnz(v) < minSamp
            continue;
        end

        xx = x(ix(v));
        yy = y(iy(v));

        xx = xx - mean(xx);
        yy = yy - mean(yy);

        den = sqrt(sum(xx.^2) * sum(yy.^2));

        if den > 0
            cc(li) = sum(xx .* yy) / den;
        end
    end
end

function plot_xcorr_panel(ax, mat, lags_t, col, ttl, xlab, ylab)

    if isempty(mat)
        title(ax, sprintf('%s\nn=0', ttl), 'Interpreter','none');
        xlabel(ax, xlab);
        ylabel(ax, ylab);
        xlim(ax, [lags_t(1), lags_t(end)]);
        xline(ax, 0, 'k--');
        ylim(ax, [-1 1]);
        axis(ax, 'square');
        box(ax, 'off');
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

    tmp = abs(mu);
    tmp(~isfinite(tmp)) = -Inf;
    [~,pki] = max(tmp);
    pkv = mu(pki);

    plot(ax, lags_t(pki), pkv, 'ko', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',6);

    title(ax, sprintf('%s\npeak r=%.3f @ %.0f ms, n=%d', ...
        ttl, pkv, 1000*lags_t(pki), n), ...
        'Interpreter','none');

    fprintf('%s: peak r=%.4f @ %.0f ms, n=%d\n', ...
        ttl, pkv, 1000*lags_t(pki), n);

    xlabel(ax, xlab);
    ylabel(ax, ylab);
    xlim(ax, [lags_t(1), lags_t(end)]);
    ylim(ax, [-1 1]);
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