% whisk_LR_abs_integral_xcorr_pop_nobp_RZ  (script)
%
% Population-only cross-correlation between locally integrated L-R whisker
% asymmetry and breathing, pooled over sessions, with NO bandpass on either
% signal (raw demeaned whisker angles, raw demeaned breathing).
%
%   aLR_inst = (xL - xR)            [or abs(...) if take_abs]
%   aLR      = local integral of L-R imbalance over sliding window
%   brw      = raw breathing (demeaned) interpolated to whisk timebase
%
% For each lag k: corr( aLR(t+k), breath(t) )
%   k > 0 : breath leads integrated L-R asymmetry
%   k < 0 : integrated L-R asymmetry leads breath
%
% Pearson corr per lag on masked (whisk-epoch) samples.

%%
take_abs = false;

% ============================ USER-EDITABLE ============================
dataRoot      = "D:\260615_thermalNbasler";
whiskDir      = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir       = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

THR_FRAC      = 0.15;        % whisk-epoch threshold (x 95th-pct envelope)
MIN_DUR       = 1;        % s, min epoch duration
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

% Force integration window odd so movsum is centered around current frame.
winSamp = max(1, round(ASYM_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

col = [0.20 0.20 0.20];

xc_sum = zeros(0, numel(lags_t));   % local int(L+R) vs breath, rows = sessions
xc     = zeros(0, numel(lags_t));   % local int(L-R) vs breath, rows = sessions
ac     = zeros(0, numel(lags_t));   % breath autocorr, rows = sessions
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

        % Local integral of L-R imbalance (centered movsum, deg*s).
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

        % cc(k) = corr(signal(t+k), -breath(t))   (breath sign flip: inhale up)
        cc_sum = lagged_corr(sLR, -brw, mask, nLag, MIN_SAMPLES);   % L+R
        cc     = lagged_corr(aLR, -brw, mask, nLag, MIN_SAMPLES);   % L-R

        % breath autocorrelation over same masked anchors (sign irrelevant)
        acc = lagged_corr(brw, brw, mask, nLag, MIN_SAMPLES);

        if any(isfinite(cc))
            xc_sum(end+1,:) = cc_sum'; %#ok<SAGROW>
            xc(end+1,:)     = cc';     %#ok<SAGROW>
            ac(end+1,:)     = acc';    %#ok<SAGROW>
            nsess = nsess + 1;
        end
    end
end

fprintf('%d sessions contributed\n', nsess);

% =============================== FIGURE ===============================
% 3 panels: int(L+R) vs breath | int(L-R) vs breath | breath autocorr.
figure('Color','w','Position',[120 200 1240 380]);

% ---- panel 1: local int(L+R) vs breath ----
ax0 = subplot(1,3,1); hold(ax0,'on'); grid(ax0,'on');
plot_xcorr_panel(ax0, xc_sum, lags_t, col, ...
    'population int(L+R) xcorr (no bandpass)', ...
    'lag (s), corr[local \int(L+R)dt(t+lag), breath(t)]', 'corr', true);

% ---- panel 2: local int(L-R) vs breath ----
ax1 = subplot(1,3,2); hold(ax1,'on'); grid(ax1,'on');
plot_xcorr_panel(ax1, xc, lags_t, col, ...
    'population int(L-R) xcorr (no bandpass)', ...
    'lag (s), corr[local \int(L-R)dt(t+lag), breath(t)]', 'corr', true);

% ---- panel 3: breath autocorrelation ----
ax2 = subplot(1,3,3); hold(ax2,'on'); grid(ax2,'on');
plot_xcorr_panel(ax2, ac, lags_t, col, ...
    'population breath autocorr (no bandpass)', ...
    'lag (s)', 'breath autocorr', false);
ylim(ax2, [-1 1]);

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
