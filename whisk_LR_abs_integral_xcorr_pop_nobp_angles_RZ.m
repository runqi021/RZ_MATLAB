% whisk_LR_abs_integral_xcorr_pop_nobp_angles_RZ  (script)
%
% Same population, same masking, same NO-bandpass signals as
% whisk_LR_abs_integral_xcorr_pop_nobp_RZ, but emits THREE figures:
%
%   FIG 1 (the one you liked, reproduced unchanged):
%     int(L+R) vs breath | int(L-R) vs breath | breath autocorr
%
%   FIG 2 (angle auto/cross):
%     L angle autocorr | R angle autocorr | L-R cross-correlation
%
%   FIG 3 (2x2, integrals vs each angle):
%     int(L+R) xcorr L angle | int(L+R) xcorr R angle
%     int(L-R) xcorr L angle | int(L-R) xcorr R angle
%
% Definitions (all RAW, demeaned; no bandpass on the correlated signals):
%   aLR_inst = (xL - xR)            [or abs(...) if take_abs]
%   sLR_inst = (xL + xR)            [or abs(...) if take_abs]
%   aLR/sLR  = local integral of L-/+R over a sliding window (deg*s)
%   xL, xR   = raw demeaned L / R whisker angles (deg)
%   brw      = raw breathing (demeaned) interpolated to the whisk timebase
%
% For each lag k:  corr( signal(t+k), ref(t) )
%   k > 0 : ref leads signal      k < 0 : signal leads ref
%   FIG 1 ref = breath;  FIG 2 L-R xcorr: signal=L, ref=R (k>0 => R leads L);
%   FIG 3 ref = the L or R angle.
% Pearson corr per lag on masked (whisk-epoch) samples. FIG 2/3 use the SAME mask
% and SAME session set as FIG 1 (so they also require breath) -- to relax that,
% drop isfinite(brw) from `mask` below.
%
% This is a NEW file -- whisk_LR_abs_integral_xcorr_pop_nobp_RZ.m is untouched.

%%
take_abs = false;

% ============================ USER-EDITABLE ============================
dataRoot      = "D:\260615_thermalNbasler";
whiskDir      = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir       = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

THR_FRAC      = 0.15;        % whisk-epoch threshold (x 95th-pct envelope)
MIN_DUR       = 1;          % s, min epoch duration
MERGE_GAP     = 0.20;       % s, merge epochs closer than this
ENV_BP        = [5 60];     % Hz, bandpass used ONLY to build the whisk-epoch envelope

ASYM_INT_WIN  = 0.050;      % s, local integration window for (xL-xR)
XC_MAXLAG     = 0.5;        % s, max lag magnitude
MIN_SAMPLES   = 50;         % min masked samples to use a session at a given lag

EXCLUDE       = "5840027";
fpsW          = 400;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

% bandpass used only for envelope-based whisk-epoch detection (not the signals)
[be,ae] = butter(3, ENV_BP/(fpsW/2), 'bandpass');

nLag    = round(XC_MAXLAG * fpsW);
lags_t  = (-nLag:nLag)' / fpsW;     % signed lag axis in seconds
nL      = numel(lags_t);

% Force integration window odd so movsum is centered around current frame.
winSamp = max(1, round(ASYM_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

col = [0.20 0.20 0.20];

% FIG 1 accumulators (rows = sessions)
xc_sum  = zeros(0,nL);   % int(L+R) vs breath
xc      = zeros(0,nL);   % int(L-R) vs breath
ac      = zeros(0,nL);   % breath autocorr
% FIG 2 accumulators
ac_L    = zeros(0,nL);   % L angle autocorr
ac_R    = zeros(0,nL);   % R angle autocorr
xc_LR   = zeros(0,nL);   % L-R cross-correlation: corr(L(t+k),R(t))
% FIG 3 accumulators
xc_sumL = zeros(0,nL);   % int(L+R) vs L angle
xc_sumR = zeros(0,nL);   % int(L+R) vs R angle
xc_aLRL = zeros(0,nL);   % int(L-R) vs L angle
xc_aLRR = zeros(0,nL);   % int(L-R) vs R angle
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
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));

        t  = (0:numel(La)-1)'/fpsW;
        N  = numel(t);

        % RAW (no bandpass) demeaned whisker angles.
        xL = fillmissing(La - mean(La,'omitnan'), 'linear');
        xR = fillmissing(Ra - mean(Ra,'omitnan'), 'linear');

        % Local integral of L-R / L+R imbalance (centered movsum, deg*s).
        if take_abs
            aLR_inst = abs(xL - xR);
            sLR_inst = abs(xL + xR);
        else
            aLR_inst = (xL - xR);
            sLR_inst = (xL + xR);
        end
        aLR = movsum(aLR_inst, winSamp, 'omitnan') / fpsW;   % L-R integral
        sLR = movsum(sLR_inst, winSamp, 'omitnan') / fpsW;   % L+R integral

        % ===================== WHISK EPOCH MASK ==========================
        % Envelope uses a bandpassed copy only to find active whisking; the
        % correlated signals above stay raw.
        xLbp = filtfilt(be,ae, xL);
        xRbp = filtfilt(be,ae, xR);
        env  = abs(hilbert(xLbp)) + abs(hilbert(xRbp));
        ep   = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);

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
                br = Bs.breath(:);
                fb = double(Bs.fps);
                tBr = (0:numel(br)-1)'/fb;

                % RAW (no bandpass) demeaned breathing.
                brr = fillmissing(br - mean(br,'omitnan'), 'linear');
                brw = interp1(tBr, brr, t, 'linear', NaN);
            end
        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
        end

        % ========================= POPULATION MASK =======================
        mask = m & isfinite(aLR) & isfinite(sLR) & isfinite(brw);
        if nnz(mask) < MIN_SAMPLES
            continue;
        end

        % ----- FIG 1: vs breath (breath sign flip: inhale up) -----
        cc_sum  = lagged_corr(sLR, -brw, mask, nLag, MIN_SAMPLES);   % int(L+R) vs breath
        cc      = lagged_corr(aLR, -brw, mask, nLag, MIN_SAMPLES);   % int(L-R) vs breath
        acc     = lagged_corr(brw,  brw, mask, nLag, MIN_SAMPLES);   % breath autocorr

        % ----- FIG 2: angle autocorr + L-R cross-correlation -----
        acc_L   = lagged_corr(xL, xL, mask, nLag, MIN_SAMPLES);      % L autocorr
        acc_R   = lagged_corr(xR, xR, mask, nLag, MIN_SAMPLES);      % R autocorr
        cc_LR   = lagged_corr(xL, xR, mask, nLag, MIN_SAMPLES);      % corr(L(t+k),R(t)); k>0 => R leads L

        % ----- FIG 3: integrals vs each angle -----
        cc_sumL = lagged_corr(sLR, xL, mask, nLag, MIN_SAMPLES);     % int(L+R) vs L
        cc_sumR = lagged_corr(sLR, xR, mask, nLag, MIN_SAMPLES);     % int(L+R) vs R
        cc_aLRL = lagged_corr(aLR, xL, mask, nLag, MIN_SAMPLES);     % int(L-R) vs L
        cc_aLRR = lagged_corr(aLR, xR, mask, nLag, MIN_SAMPLES);     % int(L-R) vs R

        if any(isfinite(cc))
            xc_sum(end+1,:)  = cc_sum';  %#ok<SAGROW>
            xc(end+1,:)      = cc';      %#ok<SAGROW>
            ac(end+1,:)      = acc';     %#ok<SAGROW>
            ac_L(end+1,:)    = acc_L';   %#ok<SAGROW>
            ac_R(end+1,:)    = acc_R';   %#ok<SAGROW>
            xc_LR(end+1,:)   = cc_LR';   %#ok<SAGROW>
            xc_sumL(end+1,:) = cc_sumL'; %#ok<SAGROW>
            xc_sumR(end+1,:) = cc_sumR'; %#ok<SAGROW>
            xc_aLRL(end+1,:) = cc_aLRL'; %#ok<SAGROW>
            xc_aLRR(end+1,:) = cc_aLRR'; %#ok<SAGROW>
            nsess = nsess + 1;
        end
    end
end

fprintf('%d sessions contributed\n', nsess);

% ============================== FIGURE 1 ==============================
% int(L+R) vs breath | int(L-R) vs breath | breath autocorr.
figure('Color','w','Position',[120 230 1240 380]);

ax0 = subplot(1,3,1); hold(ax0,'on'); grid(ax0,'on');
plot_xcorr_panel(ax0, xc_sum, lags_t, col, ...
    'population int(L+R) xcorr (no bandpass)', ...
    'lag (s), corr[local \int(L+R)dt(t+lag), breath(t)]', 'corr', true);

ax1 = subplot(1,3,2); hold(ax1,'on'); grid(ax1,'on');
plot_xcorr_panel(ax1, xc, lags_t, col, ...
    'population int(L-R) xcorr (no bandpass)', ...
    'lag (s), corr[local \int(L-R)dt(t+lag), breath(t)]', 'corr', true);

ax2 = subplot(1,3,3); hold(ax2,'on'); grid(ax2,'on');
plot_xcorr_panel(ax2, ac, lags_t, col, ...
    'population breath autocorr (no bandpass)', ...
    'lag (s)', 'breath autocorr', false);
ylim(ax2, [-1 1]);

% ============================== FIGURE 2 ==============================
% L angle autocorr | R angle autocorr | L-R cross-correlation.
figure('Color','w','Position',[140 180 1240 380]);

ax10 = subplot(1,3,1); hold(ax10,'on'); grid(ax10,'on');
plot_xcorr_panel(ax10, ac_L, lags_t, col, ...
    'population L angle autocorr (no bandpass)', ...
    'lag (s)', 'L angle autocorr', false);
ylim(ax10, [-1 1]);

ax11 = subplot(1,3,2); hold(ax11,'on'); grid(ax11,'on');
plot_xcorr_panel(ax11, ac_R, lags_t, col, ...
    'population R angle autocorr (no bandpass)', ...
    'lag (s)', 'R angle autocorr', false);
ylim(ax11, [-1 1]);

ax12 = subplot(1,3,3); hold(ax12,'on'); grid(ax12,'on');
plot_xcorr_panel(ax12, xc_LR, lags_t, col, ...
    'population L-R cross-correlation (no bandpass)', ...
    'lag (s), corr[L angle(t+lag), R angle(t)]', 'corr', true);

% ============================== FIGURE 3 ==============================
% 2x2: int(L+R) xcorr L | int(L+R) xcorr R | int(L-R) xcorr L | int(L-R) xcorr R
figure('Color','w','Position',[160 120 920 760]);

ax20 = subplot(2,2,1); hold(ax20,'on'); grid(ax20,'on');
plot_xcorr_panel(ax20, xc_sumL, lags_t, col, ...
    'population int(L+R) xcorr L angle (no bandpass)', ...
    'lag (s), corr[local \int(L+R)dt(t+lag), L angle(t)]', 'corr', true);

ax21 = subplot(2,2,2); hold(ax21,'on'); grid(ax21,'on');
plot_xcorr_panel(ax21, xc_sumR, lags_t, col, ...
    'population int(L+R) xcorr R angle (no bandpass)', ...
    'lag (s), corr[local \int(L+R)dt(t+lag), R angle(t)]', 'corr', true);

ax22 = subplot(2,2,3); hold(ax22,'on'); grid(ax22,'on');
plot_xcorr_panel(ax22, xc_aLRL, lags_t, col, ...
    'population int(L-R) xcorr L angle (no bandpass)', ...
    'lag (s), corr[local \int(L-R)dt(t+lag), L angle(t)]', 'corr', true);

ax23 = subplot(2,2,4); hold(ax23,'on'); grid(ax23,'on');
plot_xcorr_panel(ax23, xc_aLRR, lags_t, col, ...
    'population int(L-R) xcorr R angle (no bandpass)', ...
    'lag (s), corr[local \int(L-R)dt(t+lag), R angle(t)]', 'corr', true);

% ============================= LOCAL FUNCTIONS =============================

function plot_xcorr_panel(ax, mat, lags_t, col, ttl, xlab, ylab, showPeak)
% Mean +/- SEM across sessions, with lag=0 marker; optional peak annotation.
    if isempty(mat)
        title(ax, sprintf('%s\nn=0', ttl), 'FontSize', 9, 'Interpreter','none');
    else
        n  = size(mat,1);
        mu = mean(mat, 1, 'omitnan');
        se = std(mat, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(mat),1));

        fill(ax, [lags_t; flipud(lags_t)], [(mu+se)'; flipud((mu-se)')], ...
            col, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        plot(ax, lags_t, mu, '-', 'Color', col, 'LineWidth', 2);

        if showPeak
            [pkv,pki] = max(mu);
            title(ax, sprintf('%s\npeak r=%.3f @ %.0f ms  n=%d', ...
                ttl, pkv, 1000*lags_t(pki), n), 'FontSize', 9, 'Interpreter','none');
        else
            title(ax, sprintf('%s\nn=%d', ttl, n), 'FontSize', 9, 'Interpreter','none');
        end
    end
    xline(ax, 0, 'k--', 'LineWidth', 1);
    xlabel(ax, xlab);
    ylabel(ax, ylab);
    xlim(ax, [lags_t(1) lags_t(end)]);
    axis(ax, 'square'); box(ax, 'off');
end

function cc = lagged_corr(x, y, mask, nLag, minSamp)
% Pearson corr( x(t+k), y(t) ) at each lag k = -nLag..nLag, masked anchors t.
    x = x(:); y = y(:); mask = logical(mask(:));
    N = numel(x);
    lags = -nLag:nLag;
    cc = nan(numel(lags),1);
    for li = 1:numel(lags)
        k  = lags(li);
        ti = max(1,1-k) : min(N,N-k);
        v  = mask(ti) & isfinite(x(ti+k)) & isfinite(y(ti));
        if nnz(v) >= minSamp
            cc(li) = corr(x(ti(v)+k), y(ti(v)));
        end
    end
end

function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac * prctile(env,95);
    d = diff([0; a; 0]);
    s = find(d == 1);
    e = find(d == -1) - 1;
    ep = [t(s) t(e)];
    if isempty(ep), return; end
    mm = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - mm(end,2) <= mergeGap
            mm(end,2) = ep(i,2);
        else
            mm(end+1,:) = ep(i,:); %#ok<AGROW>
        end
    end
    ep = mm(mm(:,2)-mm(:,1) >= minDur,:);
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    [~,ix] = max(arrayfun(@(x) bestnum(x.name), d));
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
