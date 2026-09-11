% whisk_abs_LR_vs_breath_xcorr_shuffle_RZ.m
%
% Minimal script:
%   Only plot cross-correlation of local integral of abs(L-R) with breathing.
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   BP 5-50 Hz
%   whisk epoch mask from LEFT BP whisker Hilbert envelope
%   hard envelope threshold = 5 deg
%
% Signal:
%   abs_lr_inst = abs(xL - xR)
%   abs_int_LR  = local integral of abs(xL - xR)
%
% Breathing:
%   brw is already LP<1 Hz subtracted
%   interpolate to whisker timebase
%   optional sign flip by BREATH_SIGN
%   BP 1-20 Hz
%
% Shuffle:
%   circular time-shift breathing within each session
%
% Lag convention:
%   corr(signal(t+lag), ref(t))
%   lag > 0: breath leads abs(L-R)
%   lag < 0: abs(L-R) leads breath

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz
BREATH_BP = [1 20];   % Hz

% Match your current script with +1.
% Set to -1 if you want inverted temperature breathing.
BREATH_SIGN = +1;

ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

ABS_INT_WIN = 0.050;  % s, local integral window
XC_MAXLAG   = 5;    % s
MIN_SAMPLES = 50;

DO_SHUFFLE  = true;
N_SHUFF     = 200;
MIN_SHIFT_S = 5.0;

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

winSamp = max(1, round(ABS_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

% ============================ ACCUMULATORS ============================

xc_abs = zeros(0,nL);                 % real abs int(L-R) vs breath
xc_abs_shuf = zeros(0,nL,N_SHUFF);    % shuffled breathing control

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

                brw_signed = BREATH_SIGN * brw;
                brw_bp1_20 = filtfilt(b_br, a_br, brw_signed);
            end

        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        % ======================== ABS L-R INTEGRAL =======================

        abs_lr_inst = abs(xL - xR);
        abs_int_LR  = movsum(abs_lr_inst, winSamp, 'omitnan') / fpsW;

        mask = m & isfinite(abs_int_LR) & isfinite(brw_bp1_20);

        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        % ============================= REAL XCORR ========================

        cc_abs = lagged_corr_fast(abs_int_LR, brw_bp1_20, mask, nLag, MIN_SAMPLES);

        if ~any(isfinite(cc_abs))
            continue;
        end

        % ============================= SHUFFLE ===========================

        if DO_SHUFFLE

            cc_abs_sh = nan(nL, N_SHUFF);

            minShift = round(MIN_SHIFT_S * fpsW);
            maxShift = N - minShift;

            if maxShift > minShift

                for si = 1:N_SHUFF

                    sh = randi([minShift maxShift]);
                    br_sh = circshift(brw_bp1_20, sh);

                    cc_abs_sh(:,si) = lagged_corr_fast(abs_int_LR, br_sh, mask, nLag, MIN_SAMPLES);
                end
            end
        end

        % ============================= SAVE ==============================

        xc_abs(end+1,:) = cc_abs'; %#ok<SAGROW>

        if DO_SHUFFLE
            xc_abs_shuf(end+1,:,:) = reshape(cc_abs_sh, [1 nL N_SHUFF]); %#ok<SAGROW>
        end

        nsess = nsess + 1;
    end
end

fprintf('%d sessions contributed\n', nsess);

% ============================ SHUFFLE BAND ============================

if DO_SHUFFLE && ~isempty(xc_abs_shuf)

    mu_abs_sh = squeeze(mean(xc_abs_shuf, 1, 'omitnan'));  % nLag x nShuff

    abs_sh_lo = prctile(mu_abs_sh,  2.5, 2);
    abs_sh_hi = prctile(mu_abs_sh, 97.5, 2);
end

% =============================== PLOT =================================

figure('Color','w','Position',[180 240 520 460]);

ax = axes; hold(ax,'on'); grid(ax,'on');

if DO_SHUFFLE && exist('abs_sh_lo','var')
    plot_shuffle_band(ax, lags_t, abs_sh_lo, abs_sh_hi);
end

plot_xcorr_panel_abs(ax, xc_abs, lags_t, ...
    'abs ∫(L-R) vs breath', ...
    'lag (s), corr[local ∫|L-R|dt(t+lag), breath(t)]', ...
    'corr');

ylim(ax, [-1 1]);

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
% lag > 0: y/ref leads x/signal.
% lag < 0: x/signal leads y/ref.

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

function plot_shuffle_band(ax, lags_t, lo, hi)

    fill(ax, [lags_t; flipud(lags_t)], ...
        [hi(:); flipud(lo(:))], ...
        [0.75 0.75 0.75], ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', 'none');
end

function plot_xcorr_panel_abs(ax, mat, lags_t, ttl, xlab, ylab)

    if isempty(mat)
        title(ax, sprintf('%s\nn=0', ttl), 'Interpreter','none');
        xlabel(ax, xlab);
        ylabel(ax, ylab);
        xlim(ax, [lags_t(1), lags_t(end)]);
        xline(ax, 0, 'k--');
        axis(ax, 'square');
        box(ax, 'off');
        return;
    end

    n  = size(mat,1);
    mu = mean(mat, 1, 'omitnan');
    se = std(mat, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(mat),1));

    col = [0.2 0.2 0.2];

    fill(ax, [lags_t; flipud(lags_t)], ...
        [(mu+se)'; flipud((mu-se)')], ...
        col, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    plot(ax, lags_t, mu, 'Color', col, 'LineWidth', 2);
    xline(ax, 0, 'k--');

    % report positive peak and negative trough separately
    mu2 = mu;
    mu2(~isfinite(mu2)) = NaN;

    [pkPos, ixPos] = max(mu2);
    [pkNeg, ixNeg] = min(mu2);

    title(ax, sprintf('%s\nmax r=%.3f @ %.0f ms | min r=%.3f @ %.0f ms | n=%d', ...
        ttl, ...
        pkPos, 1000*lags_t(ixPos), ...
        pkNeg, 1000*lags_t(ixNeg), ...
        n), ...
        'Interpreter','none');

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