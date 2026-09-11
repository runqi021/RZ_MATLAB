% whisk_ITI_sync_RZ  (script)
% Combines whisk_ITI_dist_RZ + whisk_LR_sync_RZ in one pass (pooled over all
% sessions, whisk epochs only). Parameters match whisk_LR_sync_RZ.
%
%   ITI   : inter-onset intervals for protraction trough/onset events, L vs R.
%   SYNC  : L-R protraction-PEAK synchrony — each L peak paired to nearest R
%           peak within +/-MAXLAG.
%
% Figure:
%   subplot 1: protraction onset ITI, L and R
%   subplot 2: absolute |L-R protraction peak lag| in ms
%
% Synchrony cutoff:
%   sync_ms = SYNC_FRAC * mean protraction onset ITI from subplot 1

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

BP        = [5 30];     % fast-whisking bandpass (Hz)
MINDIST   = 0.025;      % s, min distance between peaks/onsets
PROT_PROM = 1;          % peak/onset prominence, x std of BP whisk
MAXLAG    = 0.08;       % s, max |L-R| lag to count as a pair

SYNC_FRAC = 0.10;       % synchronous = |lag| < 10% of mean protraction onset ITI

ITI_MAX   = 0.3;        % s, drop ITIs longer than this
NBINS     = 40;         % histogram bins for ITI
LAG_BIN_MS  = 2.5;        % ms, absolute lag histogram bin width
LAG_XMAX_MS = 100;      % ms, absolute lag x-axis max

THR_FRAC  = 0;
MIN_DUR   = 0.5;
MERGE_GAP = 0.20;
ENV_SMOOTH = 0;         % whisk-epoch gating

EXCLUDE   = "5840027";  % animal ids to drop
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
maxLag = round(MAXLAG*fpsW);

lag_ms = [];        % signed L - R protraction peak lag, ms
nsess  = 0;

ITI = struct('protL',[],'protR',[]);

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

        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);

        % protraction-positive angle convention
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

        % ---- ITI: protraction trough/onset events, L and R ----
        ITI.protL = [ITI.protL; iti(-xL, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)]; %#ok<AGROW>
        ITI.protR = [ITI.protR; iti(-xR, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)]; %#ok<AGROW>

        % ---- SYNC: protraction PEAKS, pair L peak to nearest R peak ----
        [~,iL] = findpeaks(xL, ...
            'MinPeakProminence', PROT_PROM*std(xL), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        [~,iR] = findpeaks(xR, ...
            'MinPeakProminence', PROT_PROM*std(xR), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        iL = iL(m(iL));
        iR = iR(m(iR));

        nsess = nsess + 1;

        if isempty(iL) || isempty(iR)
            continue;
        end

        for k = 1:numel(iL)
            [d,j] = min(abs(iR - iL(k)));

            if d <= maxLag
                lag_ms(end+1,1) = 1000 * (iL(k) - iR(j)) / fpsW; %#ok<AGROW>
            end
        end
    end
end

fprintf('\n%d sessions | ITI counts: protL %d protR %d\n', ...
    nsess, numel(ITI.protL), numel(ITI.protR));

assert(~isempty(lag_ms),'no paired peaks');

% ============================ SUMMARY STATS ============================

% Left subplot data: protraction onset/trough ITI
prot = [ITI.protL; ITI.protR];

mP   = mean(prot);
medP = median(prot);

% Right subplot data: absolute protraction peak L/R lag
absLag_ms = abs(lag_ms);

sync_ms = 1000 * SYNC_FRAC * mP;   % 10% of left-subplot mean ITI

medianAbsLag = median(absLag_ms);
iqrAbsLag    = iqr(absLag_ms);
fracSync     = mean(absLag_ms < sync_ms);

medianSignedLag = median(lag_ms);
iqrSignedLag    = iqr(lag_ms);

fprintf('protraction onset ITI: mean %.1f ms %.1f Hz | median %.1f ms %.1f Hz\n', ...
    1000*mP, 1/mP, 1000*medP, 1/medP);

fprintf('%d paired peaks | sync window = |lag| < %.1f ms = %.0f%% of left-subplot mean ITI\n', ...
    numel(absLag_ms), sync_ms, 100*SYNC_FRAC);

fprintf('absolute |L-R| peak lag: median %.1f ms | IQR %.1f ms | sync %.1f%%\n', ...
    medianAbsLag, iqrAbsLag, 100*fracSync);

fprintf('signed L-R peak lag: median %.1f ms | IQR %.1f ms\n', ...
    medianSignedLag, iqrSignedLag);

% ============================ FIGURE 1x2 ============================

edges = linspace(0, ITI_MAX, NBINS+1);

colL = [0.85 0.10 0.10];
colR = [0.00 0.40 0.85];

figure('Color','w','Position',[120 160 1100 460]);

% ---- subplot 1: PROTRACTION onset ITI, L and R ----
ax1 = subplot(1,2,1);
hold(ax1,'on');
grid(ax1,'on');

hh(ax1, ITI.protL, edges, colL);
hh(ax1, ITI.protR, edges, colR);

xline(ax1, mP,   'b--','LineWidth',1.4);
xline(ax1, medP, 'k--','LineWidth',1.4);

xlim(ax1,[0 ITI_MAX]);

xlabel(ax1,'inter-onset interval (s)');
ylabel(ax1,'probability');

legend(ax1, ...
    {'L','R', ...
    sprintf('mean %.0f ms (%.1f Hz)',1000*mP,1/mP), ...
    sprintf('median %.0f ms (%.1f Hz)',1000*medP,1/medP)}, ...
    'Location','northeast', ...
    'Box','off');

title(ax1, sprintf('Protraction onset ITI, L+R combined  (n=%d)  mean %.0f ms / median %.0f ms', ...
    numel(prot), 1000*mP, 1000*medP));

axis(ax1,'square');

% ---- subplot 2: ABSOLUTE L-R protraction peak lag ----
ax4 = subplot(1,2,2);
hold(ax4,'on');
grid(ax4,'on');

lagEdges = 0:LAG_BIN_MS:LAG_XMAX_MS;

histogram(ax4, absLag_ms, lagEdges, ...
    'Normalization','probability', ...
    'FaceColor',[0.40 0.40 0.40], ...
    'FaceAlpha',0.70, ...
    'EdgeColor','none');

xline(ax4, sync_ms, 'r--', 'LineWidth', 1.5);

xlabel(ax4,'|L-R protraction peak lag| (ms)');
ylabel(ax4,'% of pairs');
ytickformat(ax4,'percentage');

xlim(ax4,[0 max(lagEdges)]);
axis(ax4,'square');

title(ax4, sprintf(['L-R protraction peak synchrony  (n=%d)  ', ...
    'median |lag| %.1f ms, IQR %.1f ms, sync %.0f%% ', ...
    '(|lag| < %.1f ms = %.0f%% mean onset ITI)'], ...
    numel(absLag_ms), medianAbsLag, iqrAbsLag, ...
    100*fracSync, sync_ms, 100*SYNC_FRAC));

% ================= OPTIONAL SIGNED LAG SUPPLEMENT =================
% Uncomment if you want direction bias: negative = L before R, positive = L after R.
%
% figure('Color','w','Position',[180 180 540 460]);
% axS = axes; hold(axS,'on'); grid(axS,'on');
%
% signedEdges = -100:5:100;
%
% histogram(axS, lag_ms, signedEdges, ...
%     'Normalization','probability', ...
%     'FaceColor',[0.40 0.40 0.40], ...
%     'FaceAlpha',0.70, ...
%     'EdgeColor','none');
%
% xline(axS,0,'k--','LineWidth',1.2);
% xline(axS, sync_ms,'r:','LineWidth',1.0);
% xline(axS,-sync_ms,'r:','LineWidth',1.0);
%
% xlabel(axS,'L peak - R peak lag (ms)');
% ylabel(axS,'% of pairs');
% ytickformat(axS,'percentage');
% title(axS, sprintf('Signed L-R peak lag  median %.1f ms, IQR %.1f ms', ...
%     medianSignedLag, iqrSignedLag));
% axis(axS,'square');

% ============================ HELPERS ============================

function d = iti(sig, m, t, prom, fpsW, mindist, itimax)
    [~,loc] = findpeaks(sig, ...
        'MinPeakProminence', prom*std(sig), ...
        'MinPeakDistance', round(mindist*fpsW));

    loc = loc(m(loc));                       % keep onsets inside whisk epochs

    d = diff((loc-1)/fpsW);                  % inter-onset intervals, s
    d = d(d > 0 & d <= itimax);
    d = d(:);
end

function hh(ax, x, edges, col)
    histogram(ax, x, edges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'EdgeColor',col, ...
        'LineWidth',1.4);
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
