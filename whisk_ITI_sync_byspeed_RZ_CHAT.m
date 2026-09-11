% whisk_ITI_sync_byspeed_RZ  (script)
%
% Like whisk_ITI_sync_RZ, but SPLITS the L-R protraction-peak synchrony by whisk
% SPEED: each L protraction peak is classified SLOW (centered local whisk freq <
% SLOW_HZ) or FAST (> FAST_HZ) from its bracketing inter-peak intervals.
%
% IMPORTANT:
% Main synchrony metric is ABSOLUTE TIME:
%       absLag_ms = |L peak time - nearest R peak time| * 1000
%
% Not phase-normalized.
% This avoids unfairly penalizing fast whisking just because the cycle is shorter.
%
% Whisk epochs only, pooled over sessions; animal EXCLUDE dropped.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

BP        = [5 30];     % whisk bandpass (Hz)
MINDIST   = 0.025;      % s, min distance between peaks/onsets
PROT_PROM = 1;          % peak/onset prominence (x std of BP whisk)
MAXLAG    = 0.08;       % s, max |L-R| lag to count as a pair

ITI_MAX   = 0.3;        % s, max bracketing interval for valid centered cycle
SLOW_HZ   = 10;         % centered whisk freq < this = SLOW whisk
FAST_HZ   = 20;         % centered whisk freq > this = FAST whisk

SYNC_THRESH_MS = 7.7;    % synchronous if |L-R peak lag| <= this many ms

LAG_BIN_MS  = 5;        % histogram bin width in ms
LAG_XMAX_MS = 80;       % x-axis max in ms; should match 1000*MAXLAG

THR_FRAC  = 0;
MIN_DUR   = 0.5;
MERGE_GAP = 0.20;
ENV_SMOOTH = 0;         % whisk-epoch gating

EXCLUDE   = "5840027";  % animal ids to drop from the analysis
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

S = struct( ...
    'slow', struct('lag',[],'cyc',[],'freq',[]), ...
    'fast', struct('lag',[],'cyc',[],'freq',[]) );

nsess = 0;

ad = dir(char(dataRoot));

for ai = 1:numel(ad)

    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit'))
        continue;
    end

    if any(strcmp(ad(ai).name, EXCLUDE))
        continue;
    end

    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*'));
    [~,o] = sort({rr.name});
    rr = rr(o);

    for kk = 1:numel(rr)

        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk))))
            continue;
        end

        % ---------------- load DLC and compute protraction-positive angles ----------------
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

        t = (0:numel(La)-1)'/fpsW;

        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));

        % ---------------- whisk epoch mask ----------------
        env = movmean((abs(hilbert(xL)) + abs(hilbert(xR))) / 2, ...
            max(1,round(ENV_SMOOTH*fpsW)));

        ep = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);

        m = false(size(t));
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        % ---------------- protraction PEAKS, inside whisk epochs ----------------
        [~,iL] = findpeaks(xL, ...
            'MinPeakProminence', PROT_PROM*std(xL), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        [~,iR] = findpeaks(xR, ...
            'MinPeakProminence', PROT_PROM*std(xR), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        iL = iL(m(iL));
        iR = iR(m(iR));

        nsess = nsess + 1;

        if numel(iL) < 3 || isempty(iR)
            continue;
        end

        % sub-frame peak times
        tL = (subsamp(xL,iL)-1)/fpsW;
        tR = (subsamp(xR,iR)-1)/fpsW;

        % ---------------- classify each L peak by centered local whisk frequency ----------------
        dL = diff(tL);

        for k = 2:numel(tL)-1

            dprev = dL(k-1);
            dnext = dL(k);

            if dprev <= 0 || dnext <= 0
                continue;
            end

            if dprev > ITI_MAX || dnext > ITI_MAX
                continue;
            end

            per = (dprev + dnext) / 2;   % centered cycle period, s
            f   = 1 / per;               % centered cycle frequency, Hz

            % nearest R peak
            [d,j] = min(abs(tR - tL(k)));

            if d > MAXLAG
                continue;
            end

            lag = 1000 * (tL(k) - tR(j));   % signed L - R lag in ms

            if f < SLOW_HZ
                S.slow.lag(end+1,1)  = lag;
                S.slow.cyc(end+1,1)  = per;
                S.slow.freq(end+1,1) = f;

            elseif f > FAST_HZ
                S.fast.lag(end+1,1)  = lag;
                S.fast.cyc(end+1,1)  = per;
                S.fast.freq(end+1,1) = f;
            end
        end
    end
end

fprintf('\n%d sessions | paired peaks: slow %d, fast %d\n', ...
    nsess, numel(S.slow.lag), numel(S.fast.lag));

assert(~isempty(S.slow.lag) && ~isempty(S.fast.lag), ...
    'need paired peaks in both speed classes');

% ============================ SUMMARY STATS ============================

classes = {'slow','fast'};
labs = {
    sprintf('slow whisk <%g Hz', SLOW_HZ), ...
    sprintf('fast whisk >%g Hz', FAST_HZ)
};

cols = {
    [0.20 0.40 0.80], ...  % slow: blue
    [0.85 0.10 0.10]  ...  % fast: red
};

for c = 1:2
    cls = classes{c};

    lagSigned = S.(cls).lag(:);
    absLag    = abs(lagSigned);
    cyc       = S.(cls).cyc(:);
    freq      = S.(cls).freq(:);

    S.(cls).absLag = absLag;

    S.(cls).medianAbsLag = median(absLag,'omitnan');
    S.(cls).iqrAbsLag    = iqr(absLag);
    S.(cls).fracSync     = mean(absLag <= SYNC_THRESH_MS,'omitnan');

    S.(cls).medianSignedLag = median(lagSigned,'omitnan');
    S.(cls).iqrSignedLag    = iqr(lagSigned);

    S.(cls).meanCycleMs = 1000 * mean(cyc,'omitnan');
    S.(cls).meanFreqHz  = mean(freq,'omitnan');

    fprintf(['%s: n=%d | mean cycle %.1f ms (%.1f Hz) | ', ...
        'median |lag| %.1f ms | IQR %.1f ms | sync %.1f%% | ', ...
        'signed median %.1f ms\n'], ...
        labs{c}, ...
        numel(absLag), ...
        S.(cls).meanCycleMs, ...
        S.(cls).meanFreqHz, ...
        S.(cls).medianAbsLag, ...
        S.(cls).iqrAbsLag, ...
        100*S.(cls).fracSync, ...
        S.(cls).medianSignedLag);
end

% ============================ FIGURE: ABSOLUTE TIME LAG ============================

lagEdges = 0:LAG_BIN_MS:LAG_XMAX_MS;

figure('Color','w','Position',[200 180 680 560]);
ax = axes;
hold(ax,'on');
grid(ax,'on');

hLeg = gobjects(1,2);
txt = cell(1,2);

for c = 1:2

    cls = classes{c};
    col = cols{c};

    absLag = S.(cls).absLag(:);

    hLeg(c) = histogram(ax, absLag, lagEdges, ...
        'Normalization','probability', ...
        'FaceColor',col, ...
        'FaceAlpha',0.50, ...
        'EdgeColor',col, ...
        'LineWidth',0.7);

    txt{c} = sprintf('%s  (n=%d, %.1f Hz): median %.1f ms, IQR %.1f ms, sync %.0f%%', ...
        labs{c}, ...
        numel(absLag), ...
        S.(cls).meanFreqHz, ...
        S.(cls).medianAbsLag, ...
        S.(cls).iqrAbsLag, ...
        100*S.(cls).fracSync);
end

xline(ax, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.5);

xlim(ax,[0 LAG_XMAX_MS]);

xlabel(ax,'|L-R protraction peak lag| (ms)');
ylabel(ax,'fraction of paired whisk cycles');

legend(ax, hLeg, txt, ...
    'Location','northoutside', ...
    'Box','off');

title(ax, sprintf('L-R protraction peak synchrony by whisk speed  (%d sessions, threshold = %.0f ms)', ...
    nsess, SYNC_THRESH_MS));

axis(ax,'square');

% ============================ OPTIONAL: SIGNED LAG SUPPLEMENT ============================
% This keeps direction information: negative = L before R, positive = L after R.
% Uncomment if useful.
%
% signedEdges = -LAG_XMAX_MS:LAG_BIN_MS:LAG_XMAX_MS;
%
% figure('Color','w','Position',[240 220 680 560]);
% ax2 = axes; hold(ax2,'on'); grid(ax2,'on');
%
% for c = 1:2
%     cls = classes{c};
%     col = cols{c};
%     histogram(ax2, S.(cls).lag(:), signedEdges, ...
%         'Normalization','probability', ...
%         'FaceColor',col, ...
%         'FaceAlpha',0.45, ...
%         'EdgeColor',col, ...
%         'LineWidth',0.7);
% end
%
% xline(ax2,0,'k--','LineWidth',1.2);
% xline(ax2,SYNC_THRESH_MS,'k:','LineWidth',1.0);
% xline(ax2,-SYNC_THRESH_MS,'k:','LineWidth',1.0);
%
% xlabel(ax2,'L peak - R peak lag (ms)');
% ylabel(ax2,'fraction of paired whisk cycles');
% legend(ax2,labs,'Location','northoutside','Box','off');
% title(ax2,'Signed L-R protraction peak lag by whisk speed');
% axis(ax2,'square');

% ============================ HELPERS ============================

function p = subsamp(x, idx)
    % Parabolic sub-sample refinement of peak locations.
    % Removes frame quantization.
    x = x(:);
    idx = double(idx(:));
    p = idx;

    in = idx > 1 & idx < numel(x);
    i = idx(in);

    ym = x(i-1);
    y0 = x(i);
    yp = x(i+1);

    den = ym - 2*y0 + yp;
    delta = 0.5*(ym - yp) ./ den;

    delta(~isfinite(delta) | abs(delta) > 0.5) = 0;

    p(in) = i + delta;
end

function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac * prctile(env,95);

    d = diff([0; a; 0]);
    s = find(d == 1);
    e = find(d == -1) - 1;

    ep = [t(s) t(e)];

    if ~isempty(ep)
        m = ep(1,:);

        for i = 2:size(ep,1)
            if ep(i,1) - m(end,2) <= mergeGap
                m(end,2) = ep(i,2);
            else
                m(end+1,:) = ep(i,:); %#ok<AGROW>
            end
        end

        ep = m(m(:,2)-m(:,1) >= minDur, :);
    end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);

    bn = arrayfun(@(x) bestnum(x.name), d);
    [~,ix] = max(bn);

    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens');

    if isempty(tok)
        n = 0;
    else
        n = str2double(tok{1}{1});
    end
end