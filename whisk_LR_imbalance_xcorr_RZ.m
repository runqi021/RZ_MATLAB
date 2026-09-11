% whisk_LR_abs_integral_xcorr_RZ  (script)
%
% Cross-correlation between locally integrated |xL - xR|
% and breathing, across 5 conditions, pooled over sessions.
%   aLR_inst = abs(xL - xR)
%   aLR      = local integral of abs L-R imbalance over sliding window
%   brw      = BP breathing interpolated to whisk timebase
% For each lag k: corr( aLR(t+k), brw(t) )
% Lag convention:
%   k > 0 : breath leads integrated L-R asymmetry
%   k < 0 : integrated L-R asymmetry leads breath
% Pearson corr per lag on masked samples.
% MATLAB xcorr cannot handle masks / NaNs cleanly.

% function cc = lagged_corr(x, y, mask, nLag, minSamp)
%     x = x(:); y = y(:);
%     mask = logical(mask(:)); N = numel(x);
%     lags = -nLag:nLag;
%     cc = nan(numel(lags),1);
%     for li = 1:numel(lags)
%         k = lags(li);
%         ti = max(1,1-k) : min(N,N-k);  % Valid anchors t where both t and t+k are inside vector.
%         v = mask(ti) & isfinite(x(ti+k)) & isfinite(y(ti));
%         if nnz(v) >= minSamp
%             cc(li) = corr(x(ti(v)+k), y(ti(v)));
%         end
%     end
% end

%%
take_abs = false;

% ============================ USER-EDITABLE ============================
dataRoot      = "D:\260615_thermalNbasler";
whiskDir      = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir       = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

BP            = [2 100];     % whisk bandpass (Hz)
PROT_PROM     = 1;          % L-peak prominence (x std of BP whisk)
MINDIST       = 0.025;      % s, min distance between protraction peaks
WHISK_ITI_MAX = 0.3;        % s, max centered ITI for valid local whisk freq

THR_FRAC      = 0.1;        % whisk-epoch threshold (x 95th-pct envelope)
MIN_DUR       = 1;        % s, min epoch duration
MERGE_GAP     = 0.1;       % s, merge epochs closer than this

BR_BP         = [1 30];     % breath bandpass (Hz)
BR_PROM       = 0.5;        % inspiration trough prominence (x std of BP breath)
CYC_MAX       = 0.5;        % s, drop breath cycles longer than this

SLOW_HZ       = 10;         % whisk freq < this  = slow whisking
FAST_HZ       = 20;         % whisk freq > this  = fast whisking
BASAL_HZ      = 5;          % breath rate < this = basal breathing
SNIFF_HZ      = 8;          % breath rate > this = sniffing

ASYM_INT_WIN  = 0.050;      % s, local integration window for |xL-xR|
XC_MAXLAG     = 0.5;        % s, max lag magnitude
MIN_SAMPLES   = 50;         % min masked samples to use a session at a given lag

EXCLUDE       = "5840027";
fpsW          = 400;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

nLag    = round(XC_MAXLAG * fpsW);
lags_t  = (-nLag:nLag)' / fpsW;     % signed lag axis in seconds

% Force integration window to be odd so movsum is centered around current frame.
winSamp = max(1, round(ASYM_INT_WIN * fpsW));
if mod(winSamp,2) == 0
    winSamp = winSamp + 1;
end

condNames = {'population', ...
    sprintf('slow whisk <%g Hz',   SLOW_HZ), ...
    sprintf('fast whisk >%g Hz',   FAST_HZ), ...
    sprintf('basal breath <%g Hz', BASAL_HZ), ...
    sprintf('sniffing >%g Hz',     SNIFF_HZ)};

condCols  = {[0.20 0.20 0.20], ...
             [0.48 0.25 0.70], ...
             [0.00 0.61 0.48], ...
             [0.04 0.17 0.42], ...
             [0.94 0.29 0.14]};

nCond = numel(condNames);

xc = repmat({zeros(0,numel(lags_t))}, 1, nCond);   % cross-corr, rows = sessions
ac = repmat({zeros(0,numel(lags_t))}, 1, nCond);   % breath autocorr, rows = sessions
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

        % Left and right whisker angles.
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));

        t  = (0:numel(La)-1)'/fpsW;
        N  = numel(t);

        % Bandpass whisker angles in whisking band.
        xL = filtfilt(bw,aw, fillmissing(La - mean(La,'omitnan'), 'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra - mean(Ra,'omitnan'), 'linear'));

        % =================================================================
        % LOCAL INTEGRAL OF ABSOLUTE L-R IMBALANCE
        %
        % aLR_inst is instantaneous absolute asymmetry.
        % aLR is the local time integral over ASYM_INT_WIN.
        %
        % Units are approximately degree * second.
        %
        % Example:
        %   ASYM_INT_WIN = 0.050 s
        %   fpsW = 400 Hz
        %   winSamp ≈ 21 samples
        %
        % Because winSamp is forced odd, movsum is centered around t.
        % =================================================================
        if take_abs; aLR_inst = abs(xL - xR); else; aLR_inst = (xL - xR); end
        aLR      = movsum(aLR_inst, winSamp, 'omitnan') / fpsW;

        % ===================== LOCAL WHISK FREQUENCY =====================
        [~,iLpk] = findpeaks(xL, ...
            'MinPeakProminence', PROT_PROM * std(xL), ...
            'MinPeakDistance', round(MINDIST * fpsW));

        fwFrame = nan(N,1);

        if numel(iLpk) >= 3
            tLpk = (subsamp(xL,iLpk)-1)/fpsW;
            fLpk = nan(size(tLpk));

            for k = 2:numel(tLpk)-1
                dt = tLpk(k+1) - tLpk(k-1);
                if dt > 0 && dt <= 2*WHISK_ITI_MAX
                    fLpk(k) = 2/dt;
                end
            end

            ok = isfinite(fLpk);

            if nnz(ok) >= 2
                fwFrame = interp1(tLpk(ok), fLpk(ok), t, 'linear', NaN);
            end
        end

        % ========================= WHISK EPOCH MASK =======================
        env = abs(hilbert(xL)) + abs(hilbert(xR));
        ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);

        m = false(N,1);
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        % ============================ BREATH ==============================
        brw = nan(N,1);
        fbr = nan(N,1);

        try
            Pn = thermal_resolve_paths( ...
                pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), ...
                dataRoot);

            if isfile(Pn.breath)

                Bs = load(Pn.breath);
                br = Bs.breath(:);
                fb = double(Bs.fps);

                tBr = (0:numel(br)-1)'/fb;

                [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
                brf = filtfilt(b2,a2, fillmissing(br,'linear'));

                % Interpolate breathing to whisker timebase.
                brw = interp1(tBr, brf, t, 'linear', NaN);

                % Inspiration troughs.
                [~,iloc] = findpeaks(-brf, ...
                    'MinPeakProminence', BR_PROM * std(brf), ...
                    'MinPeakDistance', round(fb/BR_BP(2)));

                if numel(iloc) >= 2
                    tInsp = (subsamp(-brf,iloc)-1)/fb;
                    tInsp = tInsp(:);

                    cycdur = diff(tInsp);
                    bfreq  = 1 ./ cycdur;

                    bfreq(cycdur > CYC_MAX) = NaN;

                    bn  = discretize(t, tInsp);
                    okb = ~isnan(bn) & bn >= 1 & bn <= numel(bfreq);

                    fbr(okb) = bfreq(bn(okb));
                end
            end

        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
        end

        % ========================== CONDITION MASKS =======================
        base = m & isfinite(aLR) & isfinite(brw);

        masks = {base, ...
                 base & isfinite(fwFrame) & fwFrame < SLOW_HZ, ...
                 base & isfinite(fwFrame) & fwFrame > FAST_HZ, ...
                 base & isfinite(fbr)     & fbr < BASAL_HZ, ...
                 base & isfinite(fbr)     & fbr > SNIFF_HZ};

        % ===================== MASKED LAGGED CORRELATION ==================
        used = false;

        for c = 1:nCond

            if nnz(masks{c}) < MIN_SAMPLES
                continue;
            end

            % cc(k) = corr(aLR(t+k), breath(t))
            %
            % k > 0: breath leads local integrated L-R asymmetry
            % k < 0: local integrated L-R asymmetry leads breath
            cc = lagged_corr(aLR, -brw, masks{c}, nLag, MIN_SAMPLES);

            % breath autocorrelation over the same masked anchors
            % (sign of breath irrelevant for autocorr)
            acc = lagged_corr(brw, brw, masks{c}, nLag, MIN_SAMPLES);

            if any(isfinite(cc))
                xc{c}(end+1,:) = cc';
                ac{c}(end+1,:) = acc';
                used = true;
            end
        end

        if used
            nsess = nsess + 1;
        end
    end
end

fprintf('%d sessions contributed\n', nsess);

% =============================== FIGURE ===============================
% 7 panels: 5 cross-correlation + 2 breath autocorrelation (basal, sniffing).
figure('Color','w','Position',[40 80 1700 320]);

% ---- panels 1-5: cross-correlation ----
for c = 1:nCond

    ax = subplot(1,7,c);
    hold(ax,'on');
    grid(ax,'on');

    col = condCols{c};
    mat = xc{c};

    if isempty(mat)

        title(ax, sprintf('%s\nn=0', condNames{c}), ...
            'Color', col, ...
            'FontSize', 8, ...
            'Interpreter', 'none');

    else

        n = size(mat,1);

        mu = mean(mat, 1, 'omitnan');

        % Valid session count per lag.
        nLagValid = sum(isfinite(mat), 1);
        se = std(mat, 0, 1, 'omitnan') ./ sqrt(nLagValid);

        fill(ax, ...
            [lags_t; flipud(lags_t)], ...
            [(mu+se)'; flipud((mu-se)')], ...
            col, ...
            'FaceAlpha', 0.25, ...
            'EdgeColor', 'none');

        plot(ax, lags_t, mu, '-', ...
            'Color', col, ...
            'LineWidth', 2);

        [pkv,pki] = max(mu);

        title(ax, sprintf('%s\npeak r=%.3f @ %.0f ms  n=%d', ...
            condNames{c}, pkv, 1000*lags_t(pki), n), ...
            'Color', col, ...
            'FontSize', 8, ...
            'Interpreter', 'none');
    end

    xline(ax, 0, 'k--', 'LineWidth', 1);
    xlabel(ax, 'lag (s), corr[local \int|L-R|dt(t+lag), breath(t)]');
    if c == 1, ylabel(ax, 'corr'); end
    xlim(ax, [lags_t(1) lags_t(end)]);
    axis(ax, 'square');
    box(ax, 'off');
end

% ---- panels 6-7: breath autocorrelation (basal=4, sniffing=5) ----
acPanels = [4 5];   % indices into condNames / ac
for p = 1:numel(acPanels)

    c   = acPanels(p);
    ax  = subplot(1,7,nCond+p);
    hold(ax,'on');
    grid(ax,'on');

    col = condCols{c};
    mat = ac{c};

    if isempty(mat)
        title(ax, sprintf('breath ACG: %s\nn=0', condNames{c}), ...
            'Color', col, 'FontSize', 8, 'Interpreter', 'none');
    else
        n  = size(mat,1);
        mu = mean(mat, 1, 'omitnan');
        nLagValid = sum(isfinite(mat), 1);
        se = std(mat, 0, 1, 'omitnan') ./ sqrt(nLagValid);

        fill(ax, [lags_t; flipud(lags_t)], [(mu+se)'; flipud((mu-se)')], ...
            col, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        plot(ax, lags_t, mu, '-', 'Color', col, 'LineWidth', 2);

        title(ax, sprintf('breath ACG: %s\nn=%d', condNames{c}, n), ...
            'Color', col, 'FontSize', 8, 'Interpreter', 'none');
    end

    xline(ax, 0, 'k--', 'LineWidth', 1);
    ylim(ax, [-1 1]);
    xlabel(ax, 'lag (s)');
    ylabel(ax, 'breath autocorr');
    xlim(ax, [lags_t(1) lags_t(end)]);
    axis(ax, 'square');
    box(ax, 'off');
end

sgtitle(sprintf('local \\int|xL-xR|dt vs breath cross-correlation + breath autocorrelation  (%d sessions)', nsess), ...
    'FontSize', 12);

% ============================= LOCAL FUNCTIONS =============================

function cc = lagged_corr(x, y, mask, nLag, minSamp)
% Pearson corr( x(t+k), y(t) ) at each lag k = -nLag..nLag.
%
% Here:
%   x = local integrated absolute L-R whisker asymmetry
%   y = breathing signal
%
% Lag convention:
%   k > 0 : y/breath leads x/whisk-asymmetry
%   k < 0 : x/whisk-asymmetry leads y/breath
%
% We use masked anchor times t. This is why this function is used instead
% of xcorr.

    x = x(:);
    y = y(:);
    mask = logical(mask(:));

    N = numel(x);
    lags = -nLag:nLag;
    cc = nan(numel(lags),1);

    for li = 1:numel(lags)

        k = lags(li);

        % Valid anchors t where both t and t+k are inside vector.
        ti = max(1,1-k) : min(N,N-k);

        v = mask(ti) & isfinite(x(ti+k)) & isfinite(y(ti));

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

    if isempty(ep)
        return;
    end

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

function p = subsamp(x, idx)

    x = x(:);
    idx = double(idx(:));
    p = idx;

    in = idx > 1 & idx < numel(x);
    i = idx(in);

    den = x(i-1) - 2*x(i) + x(i+1);
    delta = 0.5 * (x(i-1) - x(i+1)) ./ den;

    delta(~isfinite(delta) | abs(delta) > 0.5) = 0;

    p(in) = i + delta;
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