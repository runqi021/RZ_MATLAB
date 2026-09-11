% whisk_LR_sync_2x2_breathXwhisk_RZ  (script)
%
% 2x2 stratification of L/R protraction-peak synchrony:
%
%                         basal breath < BASAL_HZ      sniff > SNIFF_HZ
%   slow whisk < SLOW_HZ       |L-R lag| ms              |L-R lag| ms
%   fast whisk > FAST_HZ       |L-R lag| ms              |L-R lag| ms
%
% Main metric:
%   absLag_ms = abs(tL_peak - tR_peak) * 1000;
%
% Not phase-normalized.
%
% For each group:
%   n paired cycles
%   median abs lag
%   IQR abs lag
%   % sync within SYNC_THRESH_MS
%
% Figures:
%   Figure 1: four histograms
%   Figure 2: summary point plot, median |lag| ms with IQR

% ============================ USER-EDITABLE ============================

dataRoot  = "D:\260615_thermalNbasler";

whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

% whisk processing
BP        = [5 30];     % whisk bandpass, Hz
MINDIST   = 0.025;      % s, min distance between whisk peaks
PROT_PROM = 1;          % whisk peak prominence, x std of BP whisk
MAXLAG    = 0.08;       % s, max |L-R| lag to count as a pair

% whisk speed classes
SLOW_HZ   = 10;         % slow whisk if centered whisk freq < this
FAST_HZ   = 20;         % fast whisk if centered whisk freq > this
ITI_MAX   = 0.3;        % s, max bracketing whisk interval for valid local freq

% breath processing
BR_BP     = [2 15];     % breath bandpass, Hz
BR_PROM   = 0.5;        % inspiration trough prominence, x std of BP breath
BASAL_HZ  = 6;          % basal breath if breath rate < this
SNIFF_HZ  = 8;          % sniff if breath rate > this
CYC_MAX   = 0.5;        % s, drop breath cycles longer than this

% synchrony metric
SYNC_THRESH_MS = 7.7;    % sync if |L-R peak lag| <= this many ms
% SYNC_THRESH_MS = 10;  % use this instead if you want stricter threshold

% display
LAG_BIN_MS  = 2.5;
LAG_XMAX_MS = 80;

% whisk epoch gating
THR_FRAC    = 0;
MIN_DUR     = 0.5;
MERGE_GAP   = 0.20;
ENV_SMOOTH  = 0;

EXCLUDE     = "5840027";
fpsW        = 400;

% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

[bbr,abr] = butter(3, BR_BP/(100/2), 'bandpass'); %#ok<NASGU>
% Actual breath filter is rebuilt per session after reading fb.

% ---------------- group storage ----------------

G = struct();

G.basal_slow = emptyGroup();
G.basal_fast = emptyGroup();
G.sniff_slow = emptyGroup();
G.sniff_fast = emptyGroup();

G.midWhisk  = emptyGroup();   % 10-20 Hz whisking, any valid breath rate
G.midBreath = emptyGroup();   % 5-8 Hz breathing, any valid whisk rate

nsessWhisk = 0;
nsessBreath = 0;

ad = dir(char(dataRoot));

for ai = 1:numel(ad)

    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit'))
        continue;
    end

    animal = ad(ai).name;

    if any(strcmp(animal, EXCLUDE))
        continue;
    end

    rr = dir(fullfile(char(dataRoot), animal, 'cam1_*'));
    [~,o] = sort({rr.name});
    rr = rr(o);

    for kk = 1:numel(rr)

        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', animal, kk))))
            continue;
        end

        % ============================ WHISK ============================

        try
            wcsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk));
            M = dlc_gate_interp(wcsv, 0.6);
        catch ME
            warning('Skipping %s n%d: whisk load failed: %s', animal, kk, ME.message);
            continue;
        end

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),   -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));

        t = (0:numel(La)-1)'/fpsW;

        xL = filtfilt(bw,aw, fillmissing(La - mean(La,'omitnan'), 'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra - mean(Ra,'omitnan'), 'linear'));

        env = movmean((abs(hilbert(xL)) + abs(hilbert(xR))) / 2, ...
            max(1, round(ENV_SMOOTH*fpsW)));

        ep = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);

        m = false(size(t));
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end

        [~,iL] = findpeaks(xL, ...
            'MinPeakProminence', PROT_PROM*std(xL), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        [~,iR] = findpeaks(xR, ...
            'MinPeakProminence', PROT_PROM*std(xR), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        iL = iL(m(iL));
        iR = iR(m(iR));

        if numel(iL) < 3 || numel(iR) < 3
            continue;
        end

        tL = (subsamp(xL,iL) - 1) / fpsW;
        tR = (subsamp(xR,iR) - 1) / fpsW;

        fL = centeredFreqAtPeaks(tL, ITI_MAX);
        fR = centeredFreqAtPeaks(tR, ITI_MAX);

        nsessWhisk = nsessWhisk + 1;

        % ============================ BREATH ============================

        try
            ncsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk));
            Pn = thermal_resolve_paths(ncsv, dataRoot);

            if ~isfile(Pn.breath)
                warning('Skipping %s n%d: no breath file', animal, kk);
                continue;
            end

            Bs = load(Pn.breath);
            br = Bs.breath(:);
            fb = double(Bs.fps);

        catch ME
            warning('Skipping %s n%d: breath load failed: %s', animal, kk, ME.message);
            continue;
        end

        [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
        brf = filtfilt(b2,a2, fillmissing(br,'linear'));

        [~,il] = findpeaks(-brf, ...
            'MinPeakProminence', BR_PROM*std(brf), ...
            'MinPeakDistance', round(fb/BR_BP(2)));

        if numel(il) < 2
            continue;
        end

        % sub-frame inspiration times, inspiration = trough of thermal breath
        tInsp = (subsamp(-brf, il) - 1) / fb;
        tInsp = tInsp(:);

        cycdur = diff(tInsp);

        if isempty(cycdur)
            continue;
        end

        nsessBreath = nsessBreath + 1;

        % ============================ PAIR L/R PEAKS AND CLASSIFY ============================

        for k = 1:numel(tL)

            [d,j] = min(abs(tR - tL(k)));

            if d > MAXLAG
                continue;
            end

            tPair = 0.5 * (tL(k) + tR(j));

            absLag_ms = abs(tL(k) - tR(j)) * 1000;
            signedLag_ms = (tL(k) - tR(j)) * 1000;

            % local whisk frequency:
            % use mean of L and R centered frequencies if both available;
            % if only one side available, mean(...,'omitnan') uses the available side.
            wfreq = mean([fL(k), fR(j)], 'omitnan');

            if isnan(wfreq)
                continue;
            end

            % breath cycle containing this paired whisk event
            bn = discretize(tPair, tInsp);

            if isnan(bn) || bn < 1 || bn > numel(cycdur)
                continue;
            end

            bdur = cycdur(bn);

            if bdur <= 0 || bdur > CYC_MAX
                continue;
            end

            brate = 1 / bdur;

            % diagnostic: accumulate middle states before 2x2 filter
            isMidW = wfreq >= SLOW_HZ && wfreq <= FAST_HZ;
            isMidB = brate >= BASAL_HZ && brate <= SNIFF_HZ;

            if isMidW
                G.midWhisk = addPair(G.midWhisk, absLag_ms, signedLag_ms, wfreq, brate, tPair);
            end

            if isMidB
                G.midBreath = addPair(G.midBreath, absLag_ms, signedLag_ms, wfreq, brate, tPair);
            end

            % drop middle whisk and middle breath states
            isSlowW = wfreq < SLOW_HZ;
            isFastW = wfreq > FAST_HZ;

            isBasal = brate < BASAL_HZ;
            isSniff = brate > SNIFF_HZ;

            if isBasal && isSlowW
                G.basal_slow = addPair(G.basal_slow, absLag_ms, signedLag_ms, wfreq, brate, tPair);

            elseif isBasal && isFastW
                G.basal_fast = addPair(G.basal_fast, absLag_ms, signedLag_ms, wfreq, brate, tPair);

            elseif isSniff && isSlowW
                G.sniff_slow = addPair(G.sniff_slow, absLag_ms, signedLag_ms, wfreq, brate, tPair);

            elseif isSniff && isFastW
                G.sniff_fast = addPair(G.sniff_fast, absLag_ms, signedLag_ms, wfreq, brate, tPair);
            end
        end
    end
end

fprintf('\nWhisk sessions used: %d | breath sessions used: %d\n', nsessWhisk, nsessBreath);

% ============================ SUMMARY ============================

keys = {'basal_slow','basal_fast','sniff_slow','sniff_fast'};

labs = { ...
    sprintf('basal <%.1f Hz + slow whisk <%.1f Hz', BASAL_HZ, SLOW_HZ), ...
    sprintf('basal <%.1f Hz + fast whisk >%.1f Hz', BASAL_HZ, FAST_HZ), ...
    sprintf('sniff >%.1f Hz + slow whisk <%.1f Hz', SNIFF_HZ, SLOW_HZ), ...
    sprintf('sniff >%.1f Hz + fast whisk >%.1f Hz', SNIFF_HZ, FAST_HZ)};

shortLabs = { ...
    'basal + slow whisk', ...
    'basal + fast whisk', ...
    'sniff + slow whisk', ...
    'sniff + fast whisk'};

cols = { ...
    [0.10 0.35 0.75], ...  % basal slow
    [0.05 0.65 0.55], ...  % basal fast
    [0.85 0.35 0.10], ...  % sniff slow
    [0.75 0.10 0.10]  ...  % sniff fast
};

n = zeros(4,1);
medAbs = nan(4,1);
iqrAbs = nan(4,1);
q1Abs = nan(4,1);
q3Abs = nan(4,1);
syncPct = nan(4,1);
medSigned = nan(4,1);
medWhiskHz = nan(4,1);
medBreathHz = nan(4,1);

fprintf('\n2x2 L/R protraction peak synchrony\n');
fprintf('Sync threshold: |L-R peak lag| <= %.1f ms\n\n', SYNC_THRESH_MS);

for i = 1:4

    x = G.(keys{i}).absLag(:);

    n(i) = numel(x);

    if ~isempty(x)
        medAbs(i) = median(x,'omitnan');
        iqrAbs(i) = iqr(x);
        q1Abs(i) = prctile(x,25);
        q3Abs(i) = prctile(x,75);
        syncPct(i) = 100 * mean(x <= SYNC_THRESH_MS,'omitnan');

        medSigned(i) = median(G.(keys{i}).signedLag,'omitnan');
        medWhiskHz(i) = median(G.(keys{i}).wfreq,'omitnan');
        medBreathHz(i) = median(G.(keys{i}).brate,'omitnan');
    end

    fprintf('%-28s | n=%5d | median |lag| %6.2f ms | IQR %6.2f ms | sync %6.1f%% | med whisk %.2f Hz | med breath %.2f Hz | signed med %.2f ms\n', ...
        shortLabs{i}, n(i), medAbs(i), iqrAbs(i), syncPct(i), medWhiskHz(i), medBreathHz(i), medSigned(i));
end

Summary = table( ...
    string(shortLabs(:)), n, medAbs, iqrAbs, q1Abs, q3Abs, syncPct, medSigned, medWhiskHz, medBreathHz, ...
    'VariableNames', {'Group','N','MedianAbsLag_ms','IQRAbsLag_ms','Q1AbsLag_ms','Q3AbsLag_ms','SyncPct','MedianSignedLag_ms','MedianWhiskHz','MedianBreathHz'});

disp(Summary);

% ============================ FIGURE 1: FOUR HISTOGRAMS ============================

lagEdges = 0:LAG_BIN_MS:LAG_XMAX_MS;

figure('Color','w','Position',[120 80 900 760]);

for i = 1:4

    ax = subplot(2,2,i);
    hold(ax,'on');
    grid(ax,'on');

    x = G.(keys{i}).absLag(:);

    histogram(ax, x, lagEdges, ...
        'Normalization','probability', ...
        'FaceColor',cols{i}, ...
        'FaceAlpha',0.70, ...
        'EdgeColor','none');

    xline(ax, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.4);

    xlim(ax,[0 LAG_XMAX_MS]);
    ylim(ax,[0 inf]);

    xlabel(ax,'|L-R protraction peak lag| (ms)');
    ylabel(ax,'fraction of paired cycles');

    title(ax, sprintf('%s\nn=%d, median %.1f ms, IQR %.1f, sync %.0f%%', ...
        shortLabs{i}, n(i), medAbs(i), iqrAbs(i), syncPct(i)));

    axis(ax,'square');
end

sgtitle(sprintf('L-R protraction peak synchrony: breath state x whisk speed  | sync <= %.1f ms', ...
    SYNC_THRESH_MS));

% ============================ FIGURE 2: SUMMARY POINT/IQR PLOT ============================

figure('Color','w','Position',[220 140 620 520]);

ax = axes;
hold(ax,'on');
grid(ax,'on');

% x-axis:
% 1 = slow whisk
% 2 = fast whisk
%
% color / offset:
% basal vs sniff

xSlow = 1;
xFast = 2;

off = 0.12;

xpos = [xSlow-off; xFast-off; xSlow+off; xFast+off];

% order here matches keys:
% basal_slow, basal_fast, sniff_slow, sniff_fast

for i = 1:4

    if isnan(medAbs(i))
        continue;
    end

    plot(ax, [xpos(i) xpos(i)], [q1Abs(i) q3Abs(i)], ...
        '-', ...
        'Color', cols{i}, ...
        'LineWidth', 3);

    scatter(ax, xpos(i), medAbs(i), 90, ...
        'filled', ...
        'MarkerFaceColor', cols{i}, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8);

    text(ax, xpos(i), q3Abs(i) + 1.5, ...
        sprintf('n=%d\n%.0f%%', n(i), syncPct(i)), ...
        'HorizontalAlignment','center', ...
        'FontSize',8);
end

xlim(ax,[0.5 2.5]);

xticks(ax,[1 2]);
xticklabels(ax,{sprintf('slow whisk\n<%.0f Hz',SLOW_HZ), sprintf('fast whisk\n>%.0f Hz',FAST_HZ)});

ylabel(ax,'median |L-R protraction peak lag| (ms)');
title(ax, sprintf('2x2 synchrony summary: median with IQR, sync <= %.1f ms', SYNC_THRESH_MS));

legendHandles = gobjects(2,1);

legendHandles(1) = plot(ax, nan, nan, 'o', ...
    'MarkerFaceColor',[0.10 0.35 0.75], ...
    'MarkerEdgeColor','k', ...
    'LineStyle','none');

legendHandles(2) = plot(ax, nan, nan, 'o', ...
    'MarkerFaceColor',[0.85 0.35 0.10], ...
    'MarkerEdgeColor','k', ...
    'LineStyle','none');

legend(ax, legendHandles, ...
    {sprintf('basal breath <%.0f Hz', BASAL_HZ), sprintf('sniff >%.0f Hz', SNIFF_HZ)}, ...
    'Location','northwest', ...
    'Box','off');

axis(ax,'square');

% ============================ SIMPLE INTERPRETATION PRINT ============================

fprintf('\nInterpretation checks:\n');

if all(~isnan(medAbs))
    fprintf('Slow vs fast during basal: slow %.2f ms, fast %.2f ms\n', ...
        medAbs(1), medAbs(2));

    fprintf('Slow vs fast during sniff: slow %.2f ms, fast %.2f ms\n', ...
        medAbs(3), medAbs(4));

    fprintf('Basal vs sniff during slow whisk: basal %.2f ms, sniff %.2f ms\n', ...
        medAbs(1), medAbs(3));

    fprintf('Basal vs sniff during fast whisk: basal %.2f ms, sniff %.2f ms\n', ...
        medAbs(2), medAbs(4));

    [bestVal,bestIdx] = min(medAbs);
    fprintf('Smallest median |lag|: %s = %.2f ms\n', shortLabs{bestIdx}, bestVal);
else
    fprintf('Some groups are empty or missing. Interpret only groups with sufficient n.\n');
end

% ============================ FIGURE 3: MIDDLE-STATE DIAGNOSTIC HISTOGRAMS ============================

midKeys = {'midWhisk','midBreath'};
midLabs = { ...
    sprintf('middle whisk %.0f-%.0f Hz, any valid breathing', SLOW_HZ, FAST_HZ), ...
    sprintf('middle breathing %.0f-%.0f Hz, any valid whisking', BASAL_HZ, SNIFF_HZ)};

midCols = { ...
    [0.45 0.45 0.45], ...
    [0.55 0.25 0.70]};

figure('Color','w','Position',[180 180 900 380]);

for i = 1:2

    ax = subplot(1,2,i);
    hold(ax,'on');
    grid(ax,'on');

    x = G.(midKeys{i}).absLag(:);

    histogram(ax, x, lagEdges, ...
        'Normalization','probability', ...
        'FaceColor',midCols{i}, ...
        'FaceAlpha',0.70, ...
        'EdgeColor','none');

    xline(ax, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.4);

    if isempty(x)
        medX  = NaN;
        iqrX  = NaN;
        syncX = NaN;
        nX    = 0;
    else
        medX  = median(x,'omitnan');
        iqrX  = iqr(x);
        syncX = 100 * mean(x <= SYNC_THRESH_MS,'omitnan');
        nX    = numel(x);
    end

    xlim(ax,[0 LAG_XMAX_MS]);

    xlabel(ax,'|L-R protraction peak lag| (ms)');
    ylabel(ax,'fraction of paired cycles');

    title(ax, sprintf('%s\nn=%d, median %.1f ms, IQR %.1f, sync %.0f%%', ...
        midLabs{i}, nX, medX, iqrX, syncX));

    axis(ax,'square');
end

sgtitle(sprintf('Excluded middle-state diagnostic: L/R protraction peak synchrony | sync <= %.1f ms', ...
    SYNC_THRESH_MS));

% ============================ MIDDLE-STATE DIAGNOSTIC PRINT ============================

fprintf('\nMiddle-state diagnostic groups:\n');

for i = 1:2
    x = G.(midKeys{i}).absLag(:);

    if isempty(x)
        fprintf('%s | n=0\n', midLabs{i});
    else
        fprintf('%s | n=%d | median |lag| %.2f ms | IQR %.2f ms | sync %.1f%% | median whisk %.2f Hz | median breath %.2f Hz\n', ...
            midLabs{i}, ...
            numel(x), ...
            median(x,'omitnan'), ...
            iqr(x), ...
            100 * mean(x <= SYNC_THRESH_MS,'omitnan'), ...
            median(G.(midKeys{i}).wfreq,'omitnan'), ...
            median(G.(midKeys{i}).brate,'omitnan'));
    end
end

% ============================ HELPERS ============================

function G = emptyGroup()
    G.absLag = [];
    G.signedLag = [];
    G.wfreq = [];
    G.brate = [];
    G.tPair = [];
end

function G = addPair(G, absLag_ms, signedLag_ms, wfreq, brate, tPair)
    G.absLag(end+1,1) = absLag_ms;
    G.signedLag(end+1,1) = signedLag_ms;
    G.wfreq(end+1,1) = wfreq;
    G.brate(end+1,1) = brate;
    G.tPair(end+1,1) = tPair;
end

function f = centeredFreqAtPeaks(tpk, itimax)
    f = nan(size(tpk));

    if numel(tpk) < 3
        return;
    end

    d = diff(tpk);

    for k = 2:numel(tpk)-1
        dprev = d(k-1);
        dnext = d(k);

        if dprev > 0 && dnext > 0 && dprev <= itimax && dnext <= itimax
            per = (dprev + dnext) / 2;
            f(k) = 1 / per;
        end
    end
end

function p = subsamp(x, idx)
    % Parabolic sub-sample refinement of peak locations.
    % Use x = -signal for trough timing.
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