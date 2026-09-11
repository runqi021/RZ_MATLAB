% whisk_LR_integral_xcorr_raw_RZ.m
%
% Population cross-correlation on RAW traces (NO bandpass) -- a new version of
% whisk_LR_integral_xcorr_filtered_RZ.m. Only xcorr/autocorr on raw traces.
%
% Whisker:
%   La/Ra -> fill NaNs -> demean        (RAW: NO bandpass on the correlated signals)
%   whisk epoch mask from the Hilbert envelope of a BP 5-50 Hz COPY of LEFT
%     whisker -- the band-pass is used ONLY to detect active whisking, never on
%     the signals that get correlated.
%   hard envelope threshold = 5 deg
%
% Breathing:
%   brw is already LP<1 Hz subtracted (detrended) AND inhale-up
%   interpolate to the whisker timebase; use AS-IS (NO bandpass, NO inversion)
%
% Signals correlated (all RAW):
%   int(L+R), int(L-R) = local movsum integrals of (xL+xR),(xL-xR)  [abs if take_abs]
%   xL, xR             = raw demeaned L/R whisker angle
%   breath (brw)       = raw detrended breathing
%
% Shuffle:
%   circular time-shift control for:
%       int(L+R) vs breath
%       int(L-R) vs breath
%       L vs R
%
% Lag convention:
%   corr(signal(t+lag), ref(t))
%   lag > 0: ref leads signal
%   lag < 0: signal leads ref
%
% Sign: breath used as-is (already inhale-up) -- matches the sign-fixed filtered
% version (brw_inv = brw, i.e. NO extra inversion). NEW file; the filtered
% version is left untouched.

clear; clc;

take_abs = true;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

ENV_BP = [5 50];      % Hz, used ONLY to build the whisk-epoch envelope (NOT the signals)

ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

ASYM_INT_WIN = 0.050; % s
XC_MAXLAG    = 0.5;   % s
MIN_SAMPLES  = 50;
XC_NORM      = 'pearson';   % normalization mode for lagged_corr_fast:
                            %   'pearson' = center + normalize -> Pearson r in [-1,1] (current)
                            %   'coeff'   = no mean removal, normalize -> uncentered NCC in [-1,1]
                            %   'cov'     = center, NO normalize -> covariance (units of signal*ref)

DO_SHUFFLE  = false;
N_SHUFF     = 200;
MIN_SHIFT_S = 5.0;    % circular shift must be at least this many seconds

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

% band-pass for the whisk-epoch ENVELOPE ONLY (signals stay raw)
[b_env, a_env] = butter(ord, ENV_BP/(fpsW/2), 'bandpass');

nLag   = round(XC_MAXLAG * fpsW);
lags_t = (-nLag:nLag)' / fpsW;
nL     = numel(lags_t);

winSamp = max(1, round(ASYM_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

% ============================ ACCUMULATORS ============================

% Figure 1
xc_sum = zeros(0,nL);   % int(L+R) vs breath
xc_lr  = zeros(0,nL);   % int(L-R) vs breath
ac_br  = zeros(0,nL);   % breath autocorr

% Figure 2
ac_L   = zeros(0,nL);   % L angle autocorr
ac_R   = zeros(0,nL);   % R angle autocorr
xc_LR  = zeros(0,nL);   % L-R cross-correlation

% Figure 3
xc_sumL = zeros(0,nL);  % int(L+R) vs L
xc_sumR = zeros(0,nL);  % int(L+R) vs R
xc_lrL  = zeros(0,nL);  % int(L-R) vs L
xc_lrR  = zeros(0,nL);  % int(L-R) vs R

% Shuffle controls
xc_sum_shuf = zeros(0,nL,N_SHUFF);  % shuffled int(L+R) vs breath
xc_lr_shuf  = zeros(0,nL,N_SHUFF);  % shuffled int(L-R) vs breath
xc_LR_shuf  = zeros(0,nL,N_SHUFF);  % shuffled L vs R

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

        % RAW demeaned whisker angles (NO bandpass on the correlated signals)
        xL = La0 - mean(La0, 'omitnan');
        xR = Ra0 - mean(Ra0, 'omitnan');

        % ======================= WHISK EPOCH MASK ========================
        % BP copy of LEFT whisker -> Hilbert envelope, used ONLY for the mask.

        xL_env  = filtfilt(b_env, a_env, xL);
        env_raw = abs(hilbert(xL_env));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));

        is_whisk = env > ENV_THR;
        ep = bool_to_epochs(is_whisk, t, MIN_DUR, MERGE_GAP);

        m = false(N,1);
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        % ============================ BREATH =============================

        brw = nan(N,1);

        try
            Pn = thermal_resolve_paths( ...
                pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);

            if isfile(Pn.breath)

                Bs = load(Pn.breath);

                br = Bs.breath(:);       % already LP<1 Hz subtracted (detrended), inhale-up
                fb = double(Bs.fps);
                tBr = (0:numel(br)-1)' / fb;

                brw = interp1(tBr, br, t, 'linear', NaN);

                if nnz(isfinite(brw)) < MIN_SAMPLES
                    continue;
                end

                brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');
                % use AS-IS: no bandpass, no inversion (already detrended + inhale-up)
            end

        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
            continue;
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

        mask = m & isfinite(int_LR) & isfinite(int_SUM) & isfinite(brw);

        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        % Figure 1: integrals vs breathing
        cc_sum = lagged_corr_fast(int_SUM, brw, mask, nLag, MIN_SAMPLES, XC_NORM);
        cc_lr  = lagged_corr_fast(int_LR,  brw, mask, nLag, MIN_SAMPLES, XC_NORM);
        acc_br = lagged_corr_fast(brw,     brw, mask, nLag, MIN_SAMPLES, XC_NORM);

        % Figure 2: angle auto/cross
        acc_L = lagged_corr_fast(xL, xL, mask, nLag, MIN_SAMPLES, XC_NORM);
        acc_R = lagged_corr_fast(xR, xR, mask, nLag, MIN_SAMPLES, XC_NORM);
        cc_LR = lagged_corr_fast(xL, xR, mask, nLag, MIN_SAMPLES, XC_NORM);
        % cc_LR = corr[L(t+lag), R(t)]
        % lag > 0 means R leads L

        % Figure 3: integrals vs angles
        cc_sumL = lagged_corr_fast(int_SUM, xL, mask, nLag, MIN_SAMPLES, XC_NORM);
        cc_sumR = lagged_corr_fast(int_SUM, xR, mask, nLag, MIN_SAMPLES, XC_NORM);
        cc_lrL  = lagged_corr_fast(int_LR,  xL, mask, nLag, MIN_SAMPLES, XC_NORM);
        cc_lrR  = lagged_corr_fast(int_LR,  xR, mask, nLag, MIN_SAMPLES, XC_NORM);

        % ============================= SHUFFLE ============================

        if DO_SHUFFLE

            cc_sum_sh = nan(nL, N_SHUFF);
            cc_lr_sh  = nan(nL, N_SHUFF);
            cc_LR_sh  = nan(nL, N_SHUFF);

            minShift = round(MIN_SHIFT_S * fpsW);
            maxShift = N - minShift;

            if maxShift > minShift

                for si = 1:N_SHUFF

                    sh = randi([minShift maxShift]);

                    br_sh = circshift(brw, sh);
                    xR_sh = circshift(xR, sh);

                    cc_sum_sh(:,si) = lagged_corr_fast(int_SUM, br_sh, mask, nLag, MIN_SAMPLES, XC_NORM);
                    cc_lr_sh(:,si)  = lagged_corr_fast(int_LR,  br_sh, mask, nLag, MIN_SAMPLES, XC_NORM);
                    cc_LR_sh(:,si)  = lagged_corr_fast(xL, xR_sh, mask, nLag, MIN_SAMPLES, XC_NORM);

                end
            end
        end

        % ============================= SAVE ===============================

        if any(isfinite(cc_lr))

            xc_sum(end+1,:) = cc_sum'; %#ok<SAGROW>
            xc_lr(end+1,:)  = cc_lr';  %#ok<SAGROW>
            ac_br(end+1,:)  = acc_br'; %#ok<SAGROW>

            ac_L(end+1,:)  = acc_L'; %#ok<SAGROW>
            ac_R(end+1,:)  = acc_R'; %#ok<SAGROW>
            xc_LR(end+1,:) = cc_LR'; %#ok<SAGROW>

            xc_sumL(end+1,:) = cc_sumL'; %#ok<SAGROW>
            xc_sumR(end+1,:) = cc_sumR'; %#ok<SAGROW>
            xc_lrL(end+1,:)  = cc_lrL';  %#ok<SAGROW>
            xc_lrR(end+1,:)  = cc_lrR';  %#ok<SAGROW>

            if DO_SHUFFLE
                xc_sum_shuf(end+1,:,:) = reshape(cc_sum_sh, [1 nL N_SHUFF]); %#ok<SAGROW>
                xc_lr_shuf(end+1,:,:)  = reshape(cc_lr_sh,  [1 nL N_SHUFF]); %#ok<SAGROW>
                xc_LR_shuf(end+1,:,:)  = reshape(cc_LR_sh,  [1 nL N_SHUFF]); %#ok<SAGROW>
            end

            nsess = nsess + 1;
        end
    end
end

fprintf('%d sessions contributed\n', nsess);

% ============================ SHUFFLE BANDS ============================

if DO_SHUFFLE && ~isempty(xc_lr_shuf)

    mu_sum_sh = squeeze(mean(xc_sum_shuf, 1, 'omitnan'));  % nLag x nShuff
    mu_lr_sh  = squeeze(mean(xc_lr_shuf,  1, 'omitnan'));
    mu_LR_sh  = squeeze(mean(xc_LR_shuf,  1, 'omitnan'));

    sum_sh_lo = prctile(mu_sum_sh, 2.5, 2);
    sum_sh_hi = prctile(mu_sum_sh, 97.5, 2);

    lr_sh_lo = prctile(mu_lr_sh, 2.5, 2);
    lr_sh_hi = prctile(mu_lr_sh, 97.5, 2);

    LR_sh_lo = prctile(mu_LR_sh, 2.5, 2);
    LR_sh_hi = prctile(mu_LR_sh, 97.5, 2);
end

col = [0.2 0.2 0.2];

% ============================== FIGURE 1 ==============================

figure('Color','w','Position',[120 230 1240 380]);

ax0 = subplot(1,3,1); hold(ax0,'on'); grid(ax0,'on');
if DO_SHUFFLE && exist('sum_sh_lo','var')
    plot_shuffle_band(ax0, lags_t, sum_sh_lo, sum_sh_hi);
end
plot_xcorr_panel(ax0, xc_sum, lags_t, col, ...
    'int(L+R) vs breath (raw)', ...
    'lag (s), corr[int(L+R)(t+lag), breath(t)]', ...
    'corr', true);

ax1 = subplot(1,3,2); hold(ax1,'on'); grid(ax1,'on');
if DO_SHUFFLE && exist('lr_sh_lo','var')
    plot_shuffle_band(ax1, lags_t, lr_sh_lo, lr_sh_hi);
end
plot_xcorr_panel(ax1, xc_lr, lags_t, col, ...
    'int(L-R) vs breath (raw)', ...
    'lag (s), corr[int(L-R)(t+lag), breath(t)]', ...
    'corr', true);

ax2 = subplot(1,3,3); hold(ax2,'on'); grid(ax2,'on');
plot_xcorr_panel(ax2, ac_br, lags_t, col, ...
    'breath autocorr (raw)', ...
    'lag (s)', ...
    'corr', false);
ylim(ax2, [-1 1]);

% ============================== FIGURE 2 ==============================

figure('Color','w','Position',[140 180 1240 380]);

ax10 = subplot(1,3,1); hold(ax10,'on'); grid(ax10,'on');
plot_xcorr_panel(ax10, ac_L, lags_t, col, ...
    'L angle autocorr (raw)', ...
    'lag (s)', ...
    'corr', false);
ylim(ax10, [-1 1]);

ax11 = subplot(1,3,2); hold(ax11,'on'); grid(ax11,'on');
plot_xcorr_panel(ax11, ac_R, lags_t, col, ...
    'R angle autocorr (raw)', ...
    'lag (s)', ...
    'corr', false);
ylim(ax11, [-1 1]);

ax12 = subplot(1,3,3); hold(ax12,'on'); grid(ax12,'on');
if DO_SHUFFLE && exist('LR_sh_lo','var')
    plot_shuffle_band(ax12, lags_t, LR_sh_lo, LR_sh_hi);
end
plot_xcorr_panel(ax12, xc_LR, lags_t, col, ...
    'L-R cross-correlation (raw)', ...
    'lag (s), corr[L(t+lag), R(t)]', ...
    'corr', true);

% ============================== FIGURE 3 ==============================

figure('Color','w','Position',[160 120 920 760]);

ax20 = subplot(2,2,1); hold(ax20,'on'); grid(ax20,'on');
plot_xcorr_panel(ax20, xc_sumL, lags_t, col, ...
    'int(L+R) vs L angle (raw)', ...
    'lag (s), corr[int(L+R)(t+lag), L(t)]', ...
    'corr', true);

ax21 = subplot(2,2,2); hold(ax21,'on'); grid(ax21,'on');
plot_xcorr_panel(ax21, xc_sumR, lags_t, col, ...
    'int(L+R) vs R angle (raw)', ...
    'lag (s), corr[int(L+R)(t+lag), R(t)]', ...
    'corr', true);

ax22 = subplot(2,2,3); hold(ax22,'on'); grid(ax22,'on');
plot_xcorr_panel(ax22, xc_lrL, lags_t, col, ...
    'int(L-R) vs L angle (raw)', ...
    'lag (s), corr[int(L-R)(t+lag), L(t)]', ...
    'corr', true);

ax23 = subplot(2,2,4); hold(ax23,'on'); grid(ax23,'on');
plot_xcorr_panel(ax23, xc_lrR, lags_t, col, ...
    'int(L-R) vs R angle (raw)', ...
    'lag (s), corr[int(L-R)(t+lag), R(t)]', ...
    'corr', true);

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

function cc = lagged_corr_fast(x, y, mask, nLag, minSamp, normMode)
% Lagged corr( x(t+lag), y(t) ) for lag = -nLag:nLag, over masked samples.
% lag > 0 means y/ref leads x/signal.
%
% normMode (all PER-LAG, evaluated over the valid/masked overlap at each lag):
%   'pearson' : center both, divide by product of L2 norms -> Pearson r in
%               [-1,1] (mean removed). Autocorr(0)=1.   [default]
%   'coeff'   : do NOT remove the mean; divide by product of L2 norms ->
%               uncentered normalized cross-correlation (cosine) in [-1,1].
%   'cov'     : center both, do NOT normalize -> covariance (units of x*y);
%               unbounded, so the autocorr panels' ylim([-1 1]) will clip it.
% NOTE: normalization is per-lag (re-evaluated over each lag's overlap), the
% right choice for masked/variable-overlap data. MATLAB's xcorr(...,'coeff')
% instead divides every lag by ONE global factor sqrt(Rxx(0)Ryy(0)) (no mask,
% no centering), so its tails decay with the shrinking overlap.

    if nargin < 6 || isempty(normMode), normMode = 'pearson'; end

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
        nv = nnz(v);

        if nv < minSamp
            continue;
        end

        xx = x(ix(v));
        yy = y(iy(v));

        switch normMode
            case 'pearson'                          % centered + normalized (Pearson r)
                xx = xx - mean(xx);
                yy = yy - mean(yy);
                den = sqrt(sum(xx.^2) * sum(yy.^2));
                if den > 0
                    cc(li) = sum(xx .* yy) / den;
                end
            case 'coeff'                            % uncentered normalized cross-corr
                den = sqrt(sum(xx.^2) * sum(yy.^2));
                if den > 0
                    cc(li) = sum(xx .* yy) / den;
                end
            case 'cov'                              % centered covariance (unnormalized)
                xx = xx - mean(xx);
                yy = yy - mean(yy);
                cc(li) = sum(xx .* yy) / nv;
            otherwise
                error('lagged_corr_fast: unknown normMode "%s"', normMode);
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

function plot_xcorr_panel(ax, mat, lags_t, col, ttl, xlab, ylab, showPeak)

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

    fill(ax, [lags_t; flipud(lags_t)], ...
        [(mu+se)'; flipud((mu-se)')], ...
        col, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    plot(ax, lags_t, mu, 'Color', col, 'LineWidth', 2);
    xline(ax, 0, 'k--');

    if showPeak
        tmp = abs(mu);
        tmp(~isfinite(tmp)) = -Inf;
        [~,pki] = max(tmp);
        pkv = mu(pki);

        title(ax, sprintf('%s\npeak r=%.3f @ %.0f ms, n=%d', ...
            ttl, pkv, 1000*lags_t(pki), n), ...
            'Interpreter','none');
    else
        title(ax, sprintf('%s\nn=%d', ttl, n), ...
            'Interpreter','none');
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
