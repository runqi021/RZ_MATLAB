% whisk_LR_setpoint_xcorr_RZ.m
%
% Setpoint analysis:
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   LP < 4 Hz = setpoint
%   BP 5-50 Hz = fast whisking, used only for whisk-epoch mask
%
% Analysis 1:
%   Lset vs Rset cross-correlation
%   dLset/dt vs dRset/dt cross-correlation
%
% Analysis 2:
%   Paired Lset/Rset session summary
%   one gray line per session: median Lset -> median Rset during whisking
%
% Lag convention:
%   corr(signal(t+lag), ref(t))
%   For Lset vs Rset:
%       corr[Lset(t+lag), Rset(t)]
%       lag > 0 means Rset leads Lset

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW = 400;

SET_LP   = 4;        % Hz, setpoint low-pass
WHISK_BP = [5 50];   % Hz, fast whisking band for epoch mask only

ENV_THR    = 5;      % deg, hard envelope threshold for whisking epochs
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;      % s
MERGE_GAP  = 0.2;    % s

XC_MAXLAG   = 1.0;   % s, wider lag for slow setpoint
MIN_SAMPLES = 50;

DO_SHUFFLE  = false;
N_SHUFF     = 200;
MIN_SHIFT_S = 5.0;

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;

[b_lp, a_lp] = butter(ord, SET_LP/(fpsW/2), 'low');
[b_bp, a_bp] = butter(ord, WHISK_BP/(fpsW/2), 'bandpass');

nLag   = round(XC_MAXLAG * fpsW);
lags_t = (-nLag:nLag)' / fpsW;
nL     = numel(lags_t);

% ============================ ACCUMULATORS ============================

% Xcorr
xc_set  = zeros(0,nL);   % Lset vs Rset
xc_dset = zeros(0,nL);   % dLset vs dRset

xc_set_shuf  = zeros(0,nL,N_SHUFF);
xc_dset_shuf = zeros(0,nL,N_SHUFF);

% Paired L/R setpoint summary
Lset_med = [];
Rset_med = [];
r_LR_set = [];

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

        % ============================ LOAD WHISKER ============================

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

        % ============================ FILTER ============================

        Lset = filtfilt(b_lp, a_lp, La0);   % LP4 left setpoint
        Rset = filtfilt(b_lp, a_lp, Ra0);   % LP4 right setpoint

        xL = filtfilt(b_bp, a_bp, La0);     % BP 5-50, mask only

        dLset = gradient(Lset) * fpsW;
        dRset = gradient(Rset) * fpsW;

        % ======================= WHISK EPOCH MASK ========================

        env_raw = abs(hilbert(xL));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));

        is_whisk = env > ENV_THR;
        ep = bool_to_epochs(is_whisk, t, MIN_DUR, MERGE_GAP);

        m = false(N,1);
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        mask = m & isfinite(Lset) & isfinite(Rset) & ...
                 isfinite(dLset) & isfinite(dRset);

        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        % ============================ XCORR ============================

        cc_set  = lagged_corr_fast(Lset,  Rset,  mask, nLag, MIN_SAMPLES);
        cc_dset = lagged_corr_fast(dLset, dRset, mask, nLag, MIN_SAMPLES);

        % ==================== PAIRED SESSION SUMMARY ====================

        Ltmp = Lset(mask);
        Rtmp = Rset(mask);

        r_LR_set(end+1,1) = corr(Ltmp, Rtmp, 'Rows','complete'); %#ok<SAGROW>

        Lset_med(end+1,1) = median(Ltmp, 'omitnan'); %#ok<SAGROW>
        Rset_med(end+1,1) = median(Rtmp, 'omitnan'); %#ok<SAGROW>

        % ============================ SHUFFLE ============================

        if DO_SHUFFLE

            cc_set_sh  = nan(nL, N_SHUFF);
            cc_dset_sh = nan(nL, N_SHUFF);

            minShift = round(MIN_SHIFT_S * fpsW);
            maxShift = N - minShift;

            if maxShift > minShift
                for si = 1:N_SHUFF

                    sh = randi([minShift maxShift]);

                    Rset_sh  = circshift(Rset,  sh);
                    dRset_sh = circshift(dRset, sh);

                    cc_set_sh(:,si)  = lagged_corr_fast(Lset,  Rset_sh,  mask, nLag, MIN_SAMPLES);
                    cc_dset_sh(:,si) = lagged_corr_fast(dLset, dRset_sh, mask, nLag, MIN_SAMPLES);
                end
            end
        end

        % ============================ SAVE ============================

        if any(isfinite(cc_set))

            xc_set(end+1,:)  = cc_set';  %#ok<SAGROW>
            xc_dset(end+1,:) = cc_dset'; %#ok<SAGROW>

            if DO_SHUFFLE
                xc_set_shuf(end+1,:,:)  = reshape(cc_set_sh,  [1 nL N_SHUFF]); %#ok<SAGROW>
                xc_dset_shuf(end+1,:,:) = reshape(cc_dset_sh, [1 nL N_SHUFF]); %#ok<SAGROW>
            end

            nsess = nsess + 1;
        end
    end
end

fprintf('%d sessions contributed\n', nsess);

% ============================ SHUFFLE BANDS ============================

if DO_SHUFFLE && ~isempty(xc_set_shuf)

    mu_set_sh  = squeeze(mean(xc_set_shuf,  1, 'omitnan'));
    mu_dset_sh = squeeze(mean(xc_dset_shuf, 1, 'omitnan'));

    set_sh_lo  = prctile(mu_set_sh,  2.5, 2);
    set_sh_hi  = prctile(mu_set_sh, 97.5, 2);

    dset_sh_lo = prctile(mu_dset_sh,  2.5, 2);
    dset_sh_hi = prctile(mu_dset_sh, 97.5, 2);
end

col = [0.2 0.2 0.2];

% ============================== FIGURE 1 ==============================
% L/R setpoint xcorr

figure('Color','w','Position',[120 230 900 380]);

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');

if DO_SHUFFLE && exist('set_sh_lo','var')
    plot_shuffle_band(ax1, lags_t, set_sh_lo, set_sh_hi);
end

plot_xcorr_panel(ax1, xc_set, lags_t, col, ...
    'LP4 L-R setpoint xcorr', ...
    'lag (s), corr[Lset(t+lag), Rset(t)]', ...
    'corr', true);

ylim(ax1, [-1 1]);

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');

if DO_SHUFFLE && exist('dset_sh_lo','var')
    plot_shuffle_band(ax2, lags_t, dset_sh_lo, dset_sh_hi);
end

plot_xcorr_panel(ax2, xc_dset, lags_t, col, ...
    'LP4 setpoint velocity xcorr', ...
    'lag (s), corr[dLset(t+lag), dRset(t)]', ...
    'corr', true);

ylim(ax2, [-1 1]);

% ============================== FIGURE 2 ==============================
% Paired Lset/Rset session summary

figure('Color','w','Position',[160 180 430 500]);

ax3 = axes; hold(ax3,'on'); grid(ax3,'on');

for i = 1:numel(Lset_med)
    plot(ax3, [1 2], [Lset_med(i) Rset_med(i)], '-', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1);
end

plot(ax3, ones(size(Lset_med)), Lset_med, 'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 5);

plot(ax3, 2*ones(size(Rset_med)), Rset_med, 'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 5);

muL = mean(Lset_med, 'omitnan');
muR = mean(Rset_med, 'omitnan');

plot(ax3, [1 2], [muL muR], 'r-', 'LineWidth', 3);
plot(ax3, [1 2], [muL muR], 'ro', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);

xlim(ax3, [0.5 2.5]);
xticks(ax3, [1 2]);
xticklabels(ax3, {'L setpoint', 'R setpoint'});

ylabel(ax3, 'LP4 setpoint angle, deg');

title(ax3, sprintf('L vs R setpoint during whisking\nmean L-R = %.2f deg, mean session r = %.3f, n = %d', ...
    mean(Lset_med - Rset_med, 'omitnan'), ...
    mean(r_LR_set, 'omitnan'), ...
    sum(isfinite(r_LR_set))), ...
    'Interpreter','none');

box(ax3,'off');

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

    x = x(:);
    y = y(:);
    mask = logical(mask(:));

    N = numel(x);
    lags = -nLag:nLag;
    cc = nan(numel(lags),1);

    for li = 1:numel(lags)

        k = lags(li);

        if k >= 0
            ix = (1+k):N;
            iy = 1:(N-k);
        else
            ix = 1:(N+k);
            iy = (1-k):N;
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
        col, ...
        'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');

    plot(ax, lags_t, mu, ...
        'Color', col, ...
        'LineWidth', 2);

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