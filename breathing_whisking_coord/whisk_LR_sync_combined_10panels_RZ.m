% whisk_LR_sync_combined_10panels_RZ.m
% Combined L/R whisker synchrony — 10-panel publication figure.
% Single event-extraction pass; all panels share one detection pipeline.
%
% 5×2 tiled layout:
%   row 1:  ITI L/R         | all L/R lag
%   row 2:  slow (<SLOW_HZ) | fast (>FAST_HZ)
%   row 3:  basal breath    | sniff breath
%   row 4:  basal+slow      | basal+fast
%   row 5:  sniff+slow      | sniff+fast

% ========================= USER-EDITABLE ==========================
dataRoot   = "D:\260615_thermalNbasler";
whiskDir   = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

BP         = [5 30];       % whisk bandpass (Hz)
MINDIST    = 0.025;        % s, min distance between protraction peaks
PROT_PROM  = 1;            % protraction peak prominence (× std)
MAXLAG     = 0.08;         % s, max accepted L/R pairing lag

ITI_MAX    = 0.3;          % s, max ITI kept in histogram
SLOW_HZ    = 10;           % Hz, L-only whisk freq upper bound for "slow"
FAST_HZ    = 20;           % Hz, L-only whisk freq lower bound for "fast"

BR_BP      = [2 15];       % breath bandpass (Hz)
BR_PROM    = 0.5;          % inspiration trough prominence (× std)
BASAL_HZ   = 6;            % Hz, per-breath rate upper bound for "basal"
SNIFF_HZ   = 8;            % Hz, per-breath rate lower bound for "sniff"
CYC_MAX    = 0.5;          % s, max accepted breath cycle duration

SYNC_THRESH_MS = 7.7;      % ms, L/R synchrony threshold (vertical line)

LAG_BIN_MS  = 2.5;         % ms, histogram bin width
LAG_XMAX_MS = 80;          % ms, x-axis maximum for all lag panels

THR_FRAC    = 0.2;            % epoch threshold = THR_FRAC × 95th-pct(envelope)
MIN_DUR     = 0.5;          % s, minimum epoch duration
MERGE_GAP   = 0.20;         % s, merge epochs separated by less than this
ENV_SMOOTH  = 0;            % frames, movmean on envelope; 0 = off

EXCLUDE     = "5840027";
fpsW        = 400;

% Per-panel x/y limits — panels 1-10 in order:
%   1:ITI | 2:all | 3:slow | 4:fast | 5:basal | 6:sniff | 7:B+sl | 8:B+fa | 9:Sn+sl | 10:Sn+fa
% XLIM: controls both histogram window and xlim; same thing.
% YLIM: [] = auto-scale from data.
XLIM = {[0  ITI_MAX],    ...  % 1  ITI
        [0  LAG_XMAX_MS],...  % 2  all
        [0  LAG_XMAX_MS],...  % 3  slow
        [0  LAG_XMAX_MS],...  % 4  fast
        [0  LAG_XMAX_MS],...  % 5  basal
        [0  LAG_XMAX_MS],...  % 6  sniff
        [0  LAG_XMAX_MS],...  % 7  basal+slow
        [0  LAG_XMAX_MS],...  % 8  basal+fast
        [0  LAG_XMAX_MS],...  % 9  sniff+slow
        [0  LAG_XMAX_MS]};    % 10 sniff+fast
YLIM = {[0 0.4], [0 0.4], [0 0.4], [0 0.4], [0 0.4], [0 0.4], [0 0.4], [0 0.4], [0 0.4], [0 0.4]};
SMOOTH_HIST_LINE = false;   % true = light gaussian smooth on overlay line; false = exact bin counts
% ==================================================================

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

% --- session discovery ---
ad = dir(char(dataRoot)); sess = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    if any(strcmp(ad(a).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*'));
    [~,o] = sort({rr.name}); rr = rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), ...
                sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess), 'no sessions with whisk csv in %s', whiskDir);

% --- event containers ---
D.all        = struct('absLag',[], 'signedLag',[], 'wfreq',[], 'tPair',[]);
D.slow       = struct('absLag',[]);
D.fast       = struct('absLag',[]);
D.basal      = struct('absLag',[], 'brate',[]);
D.sniff      = struct('absLag',[], 'brate',[]);
D.basal_slow = struct('absLag',[], 'brate',[]);
D.basal_fast = struct('absLag',[], 'brate',[]);
D.sniff_slow = struct('absLag',[], 'brate',[]);
D.sniff_fast = struct('absLag',[], 'brate',[]);
ITI.protL = []; ITI.protR = [];

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');
lagEdges = 0 : LAG_BIN_MS : LAG_XMAX_MS;

nSess = 0; nBreathValid = 0;

% --- main session loop ---
for e = 1:numel(sess)
    animal = sess{e}{1}; kk = sess{e}{2};

    % whisk DLC (required)
    try
        M = dlc_gate_interp(pick_csv(whiskDir, ...
            sprintf('%s_whisk_n%d', animal, kk)), 0.6);
    catch ME
        warning('whisk DLC load failed %s n%d: %s', animal, kk, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)' / fpsW;

    xL = filtfilt(bw,aw, fillmissing(La - mean(La,'omitnan'), 'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra - mean(Ra,'omitnan'), 'linear'));

    % whisk epoch mask
    env = movmean((abs(hilbert(xL)) + abs(hilbert(xR))) / 2, ...
                  max(1, round(ENV_SMOOTH*fpsW)));
    ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
    m   = false(size(t));
    for q = 1:size(ep,1), m = m | (t >= ep(q,1) & t <= ep(q,2)); end
    if ~any(m), continue; end

    % protraction peaks (within epochs only)
    [~,iL] = findpeaks(xL, 'MinPeakProminence', PROT_PROM*std(xL), ...
                           'MinPeakDistance',   round(MINDIST*fpsW));
    [~,iR] = findpeaks(xR, 'MinPeakProminence', PROT_PROM*std(xR), ...
                           'MinPeakDistance',   round(MINDIST*fpsW));
    iL = iL(m(iL)); iR = iR(m(iR));
    if numel(iL) < 3 || isempty(iR), continue; end

    tL = (subsamp(xL, iL) - 1) / fpsW; tL = tL(:);
    tR = (subsamp(xR, iR) - 1) / fpsW; tR = tR(:);

    % L-only centered local frequency (inner peaks only)
    fL = centeredFreqAtPeaks(tL, ITI_MAX);

    % pair each L peak (k=2:N-1) to nearest R peak within MAXLAG
    sess_abs = []; sess_sgn = []; sess_wfreq = []; sess_tPair = [];
    for k = 2:numel(tL)-1
        wfreq = fL(k);
        if isnan(wfreq), continue; end
        [d,j] = min(abs(tR - tL(k)));
        if d > MAXLAG, continue; end
        sess_abs(end+1,1)   = abs(tL(k)-tR(j))*1000; %#ok<AGROW>
        sess_sgn(end+1,1)   = (tL(k)-tR(j))*1000;    %#ok<AGROW>
        sess_wfreq(end+1,1) = wfreq;                  %#ok<AGROW>
        sess_tPair(end+1,1) = 0.5*(tL(k)+tR(j));      %#ok<AGROW>
    end

    % accumulate WHISK containers
    D.all.absLag    = [D.all.absLag;    sess_abs];
    D.all.signedLag = [D.all.signedLag; sess_sgn];
    D.all.wfreq     = [D.all.wfreq;     sess_wfreq];
    D.all.tPair     = [D.all.tPair;     sess_tPair];
    if ~isempty(sess_wfreq)
        sl = sess_wfreq < SLOW_HZ; fa = sess_wfreq > FAST_HZ;
        D.slow.absLag = [D.slow.absLag; sess_abs(sl)];
        D.fast.absLag = [D.fast.absLag; sess_abs(fa)];
    end

    % ITI from protraction onsets (troughs of angle = peaks of -angle)
    ITI.protL = [ITI.protL; iti(-xL, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)];
    ITI.protR = [ITI.protR; iti(-xR, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)];

    nSess = nSess + 1;

    % breath (non-fatal — WHISK events already committed above)
    if isempty(sess_abs), continue; end
    try
        Pn = thermal_resolve_paths( ...
            pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);
        if ~isfile(Pn.breath), continue; end
        Bs  = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
        tBr = (0:numel(br)-1)' / fb;
        [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
        brf = filtfilt(b2,a2, fillmissing(br,'linear'));
        [~,il] = findpeaks(-brf, ...
            'MinPeakProminence', BR_PROM*std(brf), ...
            'MinPeakDistance',   round(fb/BR_BP(2)));
        if numel(il) < 2, continue; end
        tInsp  = (subsamp(-brf, il) - 1) / fb; tInsp = tInsp(:);
        cycdur = diff(tInsp);

        for q = 1:numel(sess_abs)
            bn = discretize(sess_tPair(q), tInsp);
            if isnan(bn) || bn < 1 || bn > numel(cycdur), continue; end
            bdur = cycdur(bn);
            if bdur <= 0 || bdur > CYC_MAX, continue; end
            brate = 1/bdur;
            nBreathValid = nBreathValid + 1;
            alag  = sess_abs(q);
            wfreq = sess_wfreq(q);
            isB  = brate < BASAL_HZ;
            isSn = brate > SNIFF_HZ;
            isSl = wfreq < SLOW_HZ;
            isFa = wfreq > FAST_HZ;
            if isB
                D.basal.absLag = [D.basal.absLag; alag];
                D.basal.brate  = [D.basal.brate;  brate];
                if isSl
                    D.basal_slow.absLag = [D.basal_slow.absLag; alag];
                    D.basal_slow.brate  = [D.basal_slow.brate;  brate];
                elseif isFa
                    D.basal_fast.absLag = [D.basal_fast.absLag; alag];
                    D.basal_fast.brate  = [D.basal_fast.brate;  brate];
                end
            elseif isSn
                D.sniff.absLag = [D.sniff.absLag; alag];
                D.sniff.brate  = [D.sniff.brate;  brate];
                if isSl
                    D.sniff_slow.absLag = [D.sniff_slow.absLag; alag];
                    D.sniff_slow.brate  = [D.sniff_slow.brate;  brate];
                elseif isFa
                    D.sniff_fast.absLag = [D.sniff_fast.absLag; alag];
                    D.sniff_fast.brate  = [D.sniff_fast.brate;  brate];
                end
            end
        end
    catch ME
        warning('breath failed %s n%d: %s', animal, kk, ME.message);
    end
end

assert(~isempty(D.all.absLag), 'no paired L/R protraction events found');

% --- summary print ---
nW     = numel(D.all.absLag);
nSlow  = numel(D.slow.absLag);  nFast  = numel(D.fast.absLag);
nMidW  = nW - nSlow - nFast;
nBasal = numel(D.basal.absLag); nSniff = numel(D.sniff.absLag);
nMidB  = nBreathValid - nBasal - nSniff;

fprintf('\n=== L/R whisk sync: %d sessions, %d whisk pairs ===\n', nSess, nW);
fprintf('  whisk speed:  slow %d | mid %d | fast %d\n', nSlow, nMidW, nFast);
fprintf('  breath-valid: %d  (basal %d | mid %d | sniff %d)\n', ...
    nBreathValid, nBasal, nMidB, nSniff);
fprintf('\n%-18s %6s  %8s  %6s  %6s  %10s  %10s\n', ...
    'group','n','med|ms|','IQR','sync%','wHz','bHz');
fprintf('%s\n', repmat('-',1,75));
printRow('all',        D.all.absLag,       SYNC_THRESH_MS, mstr(D.all.wfreq),        '-');
printRow('slow',       D.slow.absLag,       SYNC_THRESH_MS, sprintf('<%.0f',SLOW_HZ),  '-');
printRow('fast',       D.fast.absLag,       SYNC_THRESH_MS, sprintf('>%.0f',FAST_HZ),  '-');
printRow('basal',      D.basal.absLag,      SYNC_THRESH_MS, '-',   mstr(D.basal.brate));
printRow('sniff',      D.sniff.absLag,      SYNC_THRESH_MS, '-',   mstr(D.sniff.brate));
printRow('basal+slow', D.basal_slow.absLag, SYNC_THRESH_MS, sprintf('<%.0f',SLOW_HZ), mstr(D.basal_slow.brate));
printRow('basal+fast', D.basal_fast.absLag, SYNC_THRESH_MS, sprintf('>%.0f',FAST_HZ), mstr(D.basal_fast.brate));
printRow('sniff+slow', D.sniff_slow.absLag, SYNC_THRESH_MS, sprintf('<%.0f',SLOW_HZ), mstr(D.sniff_slow.brate));
printRow('sniff+fast', D.sniff_fast.absLag, SYNC_THRESH_MS, sprintf('>%.0f',FAST_HZ), mstr(D.sniff_fast.brate));
fprintf('\n');

% --- colors ---
colL      = [0.20 0.40 0.80];
colR      = [0.85 0.10 0.10];
colAll    = [0.50 0.50 0.50];
colSlow   = [0.20 0.40 0.80];
colFast   = [0.85 0.10 0.10];
colBasal  = [0.00 0.00 0.00];
colSniff  = [0.90 0.40 0.00];
colBSlow  = [0.20 0.40 0.80];
colBFast  = [0.00 0.55 0.55];
colSSlow  = [0.95 0.60 0.00];
colSFast  = [0.55 0.00 0.10];

% --- figure ---
hdr = sprintf('L/R whisk synchrony  |  %d sessions  |  whisk pairs: %d  |  breath-valid: %d', ...
    nSess, nW, nBreathValid);
figure('Color','w', 'Position',[30 30 1380 1020], 'Name','L/R whisk synchrony');
tlo = tiledlayout(5, 2, 'TileSpacing','compact', 'Padding','compact');
title(tlo, hdr, 'FontSize', 9, 'FontWeight','normal');

itiEdges = 0 : 0.01 : ITI_MAX;

ax = nexttile(tlo);
plotITIHist(ax, ITI.protL, ITI.protR, itiEdges, colL, colR, XLIM{1}, YLIM{1});

ax = nexttile(tlo);
plotLagHist(ax, D.all.absLag, lagEdges, colAll, ...
    'all pairs', SYNC_THRESH_MS, XLIM{2}, YLIM{2}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.slow.absLag, lagEdges, colSlow, ...
    sprintf('slow  (<%.0f Hz)', SLOW_HZ), SYNC_THRESH_MS, XLIM{3}, YLIM{3}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.fast.absLag, lagEdges, colFast, ...
    sprintf('fast  (>%.0f Hz)', FAST_HZ), SYNC_THRESH_MS, XLIM{4}, YLIM{4}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.basal.absLag, lagEdges, colBasal, ...
    sprintf('basal breath  (<%.0f Hz)', BASAL_HZ), SYNC_THRESH_MS, XLIM{5}, YLIM{5}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.sniff.absLag, lagEdges, colSniff, ...
    sprintf('sniff breath  (>%.0f Hz)', SNIFF_HZ), SYNC_THRESH_MS, XLIM{6}, YLIM{6}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.basal_slow.absLag, lagEdges, colBSlow, ...
    'basal + slow', SYNC_THRESH_MS, XLIM{7}, YLIM{7}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.basal_fast.absLag, lagEdges, colBFast, ...
    'basal + fast', SYNC_THRESH_MS, XLIM{8}, YLIM{8}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.sniff_slow.absLag, lagEdges, colSSlow, ...
    'sniff + slow', SYNC_THRESH_MS, XLIM{9}, YLIM{9}, SMOOTH_HIST_LINE);

ax = nexttile(tlo);
plotLagHist(ax, D.sniff_fast.absLag, lagEdges, colSFast, ...
    'sniff + fast', SYNC_THRESH_MS, XLIM{10}, YLIM{10}, SMOOTH_HIST_LINE);

% ===================== LOCAL FUNCTIONS =====================

function plotLagHist(ax, absLag, edges, col, label, thresh, xl, yl, smoothLine)
    hold(ax,'on'); grid(ax,'on');
    if isempty(absLag)
        title(ax, sprintf('%s  (n=0)', label));
        xlim(ax, xl); axis(ax,'square'); box(ax,'off'); return;
    end
    n    = numel(absLag);
    med  = median(absLag,'omitnan');
    iqrV = iqr(absLag);
    sync = 100 * mean(absLag <= thresh, 'omitnan');

    % filled bars
    histogram(ax, absLag, edges, 'Normalization','probability', ...
        'FaceColor', col, 'EdgeColor','none');

    % overlay line — exact bin counts or optionally gaussian-smoothed
    [counts, binEdges] = histcounts(absLag, edges, 'Normalization','probability');
    if smoothLine
        yLine = smoothdata(counts, 'gaussian', 3);
    else
        yLine = counts;
    end
    stairs(ax, binEdges(1:end-1), yLine, 'Color', col, 'LineWidth', 1.8);

    xline(ax, thresh, 'k--', 'LineWidth', 1.2);
    xlim(ax, xl);
    if ~isempty(yl), ylim(ax, yl); end
    xlabel(ax,'|L-R lag| (ms)');
    ylabel(ax,'probability');
    title(ax, sprintf('%s\nn=%d   med=%.1f ms   IQR=%.1f ms   sync=%.0f%%', ...
        label, n, med, iqrV, sync));
    axis(ax,'square'); box(ax,'off');
end

function plotITIHist(ax, itiL, itiR, edges, colL, colR, xl, yl)
    hold(ax,'on'); grid(ax,'on');
    if ~isempty(itiL)
        histogram(ax, itiL, edges, 'Normalization','probability', ...
            'DisplayStyle','stairs','EdgeColor',colL,'LineWidth',1.6,'FaceColor','none');
    end
    if ~isempty(itiR)
        histogram(ax, itiR, edges, 'Normalization','probability', ...
            'DisplayStyle','stairs','EdgeColor',colR,'LineWidth',1.6,'FaceColor','none');
    end
    if ~isempty(itiL)
        xline(ax, median(itiL,'omitnan'), '-', 'Color',colL,'LineWidth',1.6, ...
            'Label','L med','LabelVerticalAlignment','top');
        xline(ax, mean(itiL,'omitnan'),  '--','Color',colL,'LineWidth',1.0);
    end
    if ~isempty(itiR)
        xline(ax, median(itiR,'omitnan'), '-', 'Color',colR,'LineWidth',1.6, ...
            'Label','R med','LabelVerticalAlignment','bottom');
        xline(ax, mean(itiR,'omitnan'),  '--','Color',colR,'LineWidth',1.0);
    end
    xlabel(ax,'inter-onset interval (s)');
    ylabel(ax,'probability');
    xlim(ax, xl);
    if ~isempty(yl), ylim(ax, yl); end
    if ~isempty(itiL) && ~isempty(itiR)
        legend(ax, sprintf('L  n=%d',numel(itiL)), sprintf('R  n=%d',numel(itiR)), ...
            'Location','northeast','Box','off');
        title(ax, sprintf('protraction onset ITI\nL med=%.3f s     R med=%.3f s', ...
            median(itiL,'omitnan'), median(itiR,'omitnan')));
    elseif ~isempty(itiL)
        title(ax, sprintf('protraction onset ITI\nL med=%.3f s', median(itiL,'omitnan')));
    else
        title(ax,'protraction onset ITI  (no data)');
    end
    axis(ax,'square'); box(ax,'off');
end

function printRow(name, lag, thresh, wStr, bStr)
    n = numel(lag);
    if n == 0
        fprintf('%-18s %6d  %8s  %6s  %6s  %10s  %10s\n', ...
            name,0,'-','-','-',wStr,bStr);
        return;
    end
    med  = median(lag,'omitnan');
    iqrV = iqr(lag);
    sync = 100 * mean(lag <= thresh,'omitnan');
    fprintf('%-18s %6d  %8.2f  %6.2f  %6.1f  %10s  %10s\n', ...
        name,n,med,iqrV,sync,wStr,bStr);
end

function s = mstr(v)
    if isempty(v), s = '-'; else, s = sprintf('%.1f', median(v,'omitnan')); end
end

function p = subsamp(x, idx)
    x = x(:); idx = double(idx(:)); p = idx;
    in = idx > 1 & idx < numel(x); i = idx(in);
    ym = x(i-1); y0 = x(i); yp = x(i+1);
    den = ym - 2*y0 + yp;
    delta = 0.5*(ym-yp) ./ den;
    delta(~isfinite(delta) | abs(delta) > 0.5) = 0;
    p(in) = i + delta;
end

function f = centeredFreqAtPeaks(t, iti_max)
    n = numel(t); f = nan(n,1);
    for k = 2:n-1
        dt = t(k+1) - t(k-1);
        if dt > 0 && dt <= 2*iti_max
            f(k) = 2/dt;
        end
    end
end

function intervals = iti(sig, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)
    if ~any(m), intervals = []; return; end
    sig = sig(:);
    [~,idx] = findpeaks(sig, ...
        'MinPeakProminence', PROT_PROM * std(sig(m)), ...
        'MinPeakDistance',   round(MINDIST * fpsW));
    idx = idx(m(idx));
    if numel(idx) < 2, intervals = []; return; end
    pk_t = (subsamp(sig, idx) - 1) / fpsW; pk_t = pk_t(:);
    ivals = diff(pk_t);
    intervals = ivals(ivals > 0 & ivals < ITI_MAX);
end

function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac * prctile(env, 95);
    d = diff([0; a; 0]); s = find(d==1); e = find(d==-1)-1;
    ep = [t(s) t(e)];
    if isempty(ep), return; end
    m = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - m(end,2) <= mergeGap
            m(end,2) = ep(i,2);
        else
            m(end+1,:) = ep(i,:); %#ok<AGROW>
        end
    end
    ep = m(m(:,2) - m(:,1) >= minDur, :);
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv for %s in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix] = max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
