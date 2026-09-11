% whisk_LR_lag_vs_breath_whisk_freq_RZ  (script)
%
% Continuous analysis of L/R whisker protraction-peak synchrony.
% All valid paired events are collected (no 2x2 category filtering).
%
% Metric: absLag_ms = |tL_peak - tR_peak| * 1000;  not phase-normalized.
%
% Figures:
%   Figure 1: |L-R lag| vs breathing frequency  (scatter + binned)
%   Figure 2: |L-R lag| vs whisking frequency   (scatter + binned)
%   Figure 3: 2D heatmap — median |lag| (ms) by breath-rate x whisk-freq
%             + event-count map
%
% Detection code and helpers are identical to
% whisk_LR_sync_2x2_breathXwhisk_RZ_CHAT.m

% ============================ USER-EDITABLE ============================

dataRoot  = "D:\260615_thermalNbasler";

whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

% whisk processing
BP        = [5 30];     % whisk bandpass, Hz
MINDIST   = 0.025;      % s, min distance between whisk peaks
PROT_PROM = 1;          % whisk peak prominence, x std of BP whisk
MAXLAG    = 0.08;       % s, max |L-R| lag to count as a pair

% local whisk frequency
ITI_MAX   = 0.3;        % s, max bracketing whisk interval for valid local freq

% breath processing
BR_BP     = [2 15];     % breath bandpass, Hz
BR_PROM   = 0.5;        % inspiration trough prominence, x std of BP breath
CYC_MAX   = 0.5;        % s, drop breath cycles longer than this

% synchrony metric
SYNC_THRESH_MS = 7.7;

% display
LAG_YMAX_MS = 80;
SHOW_FIT    = true;   % show linear fit + R/p on scatter plots

% frequency binning
BREATH_EDGES = 2:0.5:12;   % Hz bin edges
WHISK_EDGES  = 5:1:30;     % Hz bin edges
MIN_BIN_N    = 10;          % min events to show a bin

% reference lines (for plotting only; no event filtering)
BASAL_HZ  = 6;
SNIFF_HZ  = 8;
SLOW_HZ   = 10;
FAST_HZ   = 20;

% whisk epoch gating
THR_FRAC    = 0;
MIN_DUR     = 0.5;
MERGE_GAP   = 0.20;
ENV_SMOOTH  = 0;

EXCLUDE     = "5840027";
fpsW        = 400;

% ======================================================================

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

% breath filter rebuilt per session from actual fps

% -------- WHISK: old-script-compatible, L-only freq, no breath requirement --------
WHISK.absLag       = [];
WHISK.signedLag    = [];
WHISK.wfreq_Lonly  = [];
WHISK.wfreq_meanLR = [];
WHISK.tPair        = [];
WHISK.animal       = {};
WHISK.session      = [];

% -------- BREATH: breath-valid subset, requires valid breath cycle --------
BREATH.absLag       = [];
BREATH.signedLag    = [];
BREATH.wfreq_Lonly  = [];
BREATH.wfreq_meanLR = [];
BREATH.brate        = [];
BREATH.tPair        = [];
BREATH.animal       = {};
BREATH.session      = [];

nsessWhisk  = 0;
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
            M    = dlc_gate_interp(wcsv, 0.6);
        catch ME
            warning('Skipping %s n%d: whisk load failed: %s', animal, kk, ME.message);
            continue;
        end

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),   -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));

        t = (0:numel(La)-1)' / fpsW;

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
            'MinPeakDistance',   round(MINDIST*fpsW));

        [~,iR] = findpeaks(xR, ...
            'MinPeakProminence', PROT_PROM*std(xR), ...
            'MinPeakDistance',   round(MINDIST*fpsW));

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

        % --- WHISK events: L-only freq, nearest R within MAXLAG, no breath ---
        for k = 2:numel(tL)-1

            wfreq_Lonly = fL(k);

            if isnan(wfreq_Lonly)
                continue;
            end

            [d,j] = min(abs(tR - tL(k)));

            if d > MAXLAG
                continue;
            end

            tPair        = 0.5 * (tL(k) + tR(j));
            absLag_ms    = abs(tL(k) - tR(j)) * 1000;
            signedLag_ms = (tL(k) - tR(j)) * 1000;
            wfreq_meanLR = mean([fL(k), fR(j)], 'omitnan');

            WHISK.absLag(end+1,1)       = absLag_ms;
            WHISK.signedLag(end+1,1)    = signedLag_ms;
            WHISK.wfreq_Lonly(end+1,1)  = wfreq_Lonly;
            WHISK.wfreq_meanLR(end+1,1) = wfreq_meanLR;
            WHISK.tPair(end+1,1)        = tPair;
            WHISK.animal{end+1,1}       = animal;
            WHISK.session(end+1,1)      = kk;
        end

        % ============================ BREATH ============================

        hasBreath = false;
        tInsp     = [];
        cycdur    = [];

        try
            ncsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk));
            Pn   = thermal_resolve_paths(ncsv, dataRoot);

            if isfile(Pn.breath)
                Bs  = load(Pn.breath);
                br  = Bs.breath(:);
                fb  = double(Bs.fps);

                [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
                brf = filtfilt(b2,a2, fillmissing(br,'linear'));

                [~,il] = findpeaks(-brf, ...
                    'MinPeakProminence', BR_PROM*std(brf), ...
                    'MinPeakDistance',   round(fb/BR_BP(2)));

                if numel(il) >= 2
                    tInsp  = (subsamp(-brf, il) - 1) / fb;
                    tInsp  = tInsp(:);
                    cycdur = diff(tInsp);

                    if ~isempty(cycdur)
                        hasBreath   = true;
                        nsessBreath = nsessBreath + 1;
                    end
                end
            else
                warning('No breath file: %s n%d', animal, kk);
            end
        catch ME
            warning('Breath load failed: %s n%d: %s', animal, kk, ME.message);
        end

        if ~hasBreath
            continue;
        end

        % --- BREATH events: same pairing as WHISK, but also requires valid breath cycle ---
        for k = 2:numel(tL)-1

            wfreq_Lonly = fL(k);

            if isnan(wfreq_Lonly)
                continue;
            end

            [d,j] = min(abs(tR - tL(k)));

            if d > MAXLAG
                continue;
            end

            tPair        = 0.5 * (tL(k) + tR(j));
            absLag_ms    = abs(tL(k) - tR(j)) * 1000;
            signedLag_ms = (tL(k) - tR(j)) * 1000;
            wfreq_meanLR = mean([fL(k), fR(j)], 'omitnan');

            bn = discretize(tPair, tInsp);

            if isnan(bn) || bn < 1 || bn > numel(cycdur)
                continue;
            end

            bdur = cycdur(bn);

            if bdur <= 0 || bdur > CYC_MAX
                continue;
            end

            brate = 1 / bdur;

            BREATH.absLag(end+1,1)       = absLag_ms;
            BREATH.signedLag(end+1,1)    = signedLag_ms;
            BREATH.wfreq_Lonly(end+1,1)  = wfreq_Lonly;
            BREATH.wfreq_meanLR(end+1,1) = wfreq_meanLR;
            BREATH.brate(end+1,1)        = brate;
            BREATH.tPair(end+1,1)        = tPair;
            BREATH.animal{end+1,1}       = animal;
            BREATH.session(end+1,1)      = kk;
        end
    end
end

fprintf('\nWhisk sessions: %d | breath sessions: %d\n', nsessWhisk, nsessBreath);
fprintf('WHISK events: %d | BREATH events: %d\n', numel(WHISK.absLag), numel(BREATH.absLag));

% ============================ VALIDATION PRINTOUT ============================

fprintf('\nOld-compatible WHISK event set:\n');
fprintf('WHISK total paired events: %d\n', numel(WHISK.absLag));
fprintf('WHISK slow  < %.1f Hz: %d\n',   SLOW_HZ, sum(WHISK.wfreq_Lonly < SLOW_HZ));
fprintf('WHISK middle %.1f-%.1f Hz: %d\n', SLOW_HZ, FAST_HZ, ...
    sum(WHISK.wfreq_Lonly >= SLOW_HZ & WHISK.wfreq_Lonly <= FAST_HZ));
fprintf('WHISK fast  > %.1f Hz: %d\n',   FAST_HZ, sum(WHISK.wfreq_Lonly > FAST_HZ));

fprintf('\nBREATH-valid event set:\n');
fprintf('BREATH total paired events: %d\n', numel(BREATH.absLag));
fprintf('BREATH slow  < %.1f Hz: %d\n',   SLOW_HZ, sum(BREATH.wfreq_Lonly < SLOW_HZ));
fprintf('BREATH middle %.1f-%.1f Hz: %d\n', SLOW_HZ, FAST_HZ, ...
    sum(BREATH.wfreq_Lonly >= SLOW_HZ & BREATH.wfreq_Lonly <= FAST_HZ));
fprintf('BREATH fast  > %.1f Hz: %d\n',   FAST_HZ, sum(BREATH.wfreq_Lonly > FAST_HZ));

slowMask = WHISK.wfreq_Lonly < SLOW_HZ;
fastMask = WHISK.wfreq_Lonly > FAST_HZ;

fprintf('\nOld by-speed comparable lag stats using WHISK:\n');
fprintf('slow: n=%d | median |lag| %.2f ms | IQR %.2f ms | sync %.1f%%\n', ...
    sum(slowMask), ...
    median(WHISK.absLag(slowMask),'omitnan'), ...
    iqr(WHISK.absLag(slowMask)), ...
    100*mean(WHISK.absLag(slowMask) <= SYNC_THRESH_MS,'omitnan'));

fprintf('fast: n=%d | median |lag| %.2f ms | IQR %.2f ms | sync %.1f%%\n', ...
    sum(fastMask), ...
    median(WHISK.absLag(fastMask),'omitnan'), ...
    iqr(WHISK.absLag(fastMask)), ...
    100*mean(WHISK.absLag(fastMask) <= SYNC_THRESH_MS,'omitnan'));

% ============================ BINNED STATS ============================

breathCenters = BREATH_EDGES(1:end-1) + diff(BREATH_EDGES)/2;
whiskCenters  = WHISK_EDGES(1:end-1)  + diff(WHISK_EDGES)/2;

nBb = numel(breathCenters);
nBw = numel(whiskCenters);

% Bb uses BREATH (needs brate)
% Bw uses WHISK.wfreq_Lonly (comparable to old by-speed script)
Bb = binnedStats(BREATH.brate,        BREATH.absLag, BREATH_EDGES, SYNC_THRESH_MS);
Bw = binnedStats(WHISK.wfreq_Lonly,   WHISK.absLag,  WHISK_EDGES,  SYNC_THRESH_MS);

% 2D bins — use BREATH (requires both brate and wfreq_Lonly)
medMap = nan(nBw, nBb);
nMap   = zeros(nBw, nBb);

for ib = 1:nBb
    for iw = 1:nBw
        mask = BREATH.brate        >= BREATH_EDGES(ib) & BREATH.brate        < BREATH_EDGES(ib+1) & ...
               BREATH.wfreq_Lonly  >= WHISK_EDGES(iw)  & BREATH.wfreq_Lonly  < WHISK_EDGES(iw+1);
        cnt = sum(mask);
        nMap(iw,ib) = cnt;
        if cnt >= MIN_BIN_N
            medMap(iw,ib) = median(BREATH.absLag(mask), 'omitnan');
        end
    end
end

% ============================ PRINT TABLES ============================

fprintf('\nBreath-frequency binned table [BREATH set] (sync threshold = %.1f ms):\n', SYNC_THRESH_MS);
fprintf('%-12s %6s %10s %8s %8s %8s %8s\n', ...
    'Center(Hz)', 'n', 'med(ms)', 'Q1(ms)', 'Q3(ms)', 'IQR(ms)', 'sync%');

for i = 1:nBb
    if Bb.n(i) == 0; continue; end
    fprintf('%-12.2f %6d %10.2f %8.2f %8.2f %8.2f %8.1f\n', ...
        breathCenters(i), Bb.n(i), Bb.med(i), Bb.q1(i), Bb.q3(i), Bb.iqr(i), Bb.sync(i));
end

fprintf('\nWhisk-frequency binned table [WHISK set, L-only freq] (sync threshold = %.1f ms):\n', SYNC_THRESH_MS);
fprintf('%-12s %6s %10s %8s %8s %8s %8s\n', ...
    'Center(Hz)', 'n', 'med(ms)', 'Q1(ms)', 'Q3(ms)', 'IQR(ms)', 'sync%');

for i = 1:nBw
    if Bw.n(i) == 0; continue; end
    fprintf('%-12.2f %6d %10.2f %8.2f %8.2f %8.2f %8.1f\n', ...
        whiskCenters(i), Bw.n(i), Bw.med(i), Bw.q1(i), Bw.q3(i), Bw.iqr(i), Bw.sync(i));
end

% ============================ FIGURE 1: LAG vs BREATH FREQ  [BREATH set] ============================

figure('Color','w','Position',[100 200 1000 440]);
sgtitle(sprintf('|L-R protraction peak lag| vs breathing frequency  [BREATH set, n=%d]', ...
    numel(BREATH.absLag)));

% left: raw scatter
axB1 = subplot(1,2,1);
hold(axB1,'on');  grid(axB1,'on');

scatter(axB1, BREATH.brate, BREATH.absLag, 8, ...
    'filled', ...
    'MarkerFaceColor', [0.25 0.25 0.25], ...
    'MarkerFaceAlpha', 0.15, ...
    'MarkerEdgeAlpha', 0.15);

if SHOW_FIT
    pb    = polyfit(BREATH.brate, BREATH.absLag, 1);
    xfitB = linspace(BREATH_EDGES(1), BREATH_EDGES(end), 200);
    plot(axB1, xfitB, polyval(pb, xfitB), 'r-', 'LineWidth', 1.5);
    [rB, pB] = corr(BREATH.brate, BREATH.absLag, 'Type', 'Pearson');
    text(axB1, 0.04, 0.93, sprintf('R = %.3f,  p = %.2e', rB, pB), ...
        'Units','normalized', 'FontSize', 9, 'Color', [0.75 0 0]);
end

yline(axB1, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.2);
xline(axB1, BASAL_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
xline(axB1, SNIFF_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

xlim(axB1, [BREATH_EDGES(1) BREATH_EDGES(end)]);
ylim(axB1, [0 LAG_YMAX_MS]);
xlabel(axB1, 'breathing frequency (Hz)');
ylabel(axB1, '|L-R protraction peak lag| (ms)');
title(axB1, 'Raw paired events');

% right: binned median
axB2 = subplot(1,2,2);
hold(axB2,'on');  grid(axB2,'on');

okB  = Bb.n >= MIN_BIN_N;
xbOk = breathCenters(okB);

errorbar(axB2, xbOk, Bb.med(okB), ...
    Bb.med(okB) - Bb.q1(okB), ...
    Bb.q3(okB)  - Bb.med(okB), ...
    'o-', ...
    'Color',           [0.25 0.55 0.85], ...
    'LineWidth',       1.8, ...
    'CapSize',         5, ...
    'MarkerFaceColor', [0.25 0.55 0.85], ...
    'MarkerSize',      7);

yline(axB2, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.2);
xline(axB2, BASAL_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
xline(axB2, SNIFF_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

xlim(axB2, [BREATH_EDGES(1) BREATH_EDGES(end)]);
ylim(axB2, [0 LAG_YMAX_MS]);
xlabel(axB2, 'breathing frequency (Hz)');
ylabel(axB2, 'median |L-R protraction peak lag| (ms)');
title(axB2, sprintf('Binned (%.1f Hz bins, min n=%d, bars=Q1–Q3)', ...
    diff(BREATH_EDGES(1:2)), MIN_BIN_N));
subtitle(axB2, sprintf('dashed = basal/sniff boundaries | sync \\leq %.1f ms', SYNC_THRESH_MS));

% ============================ FIGURE 2: LAG vs WHISK FREQ  [WHISK set, L-only] ============================

figure('Color','w','Position',[140 160 1000 440]);
sgtitle(sprintf('|L-R protraction peak lag| vs whisking frequency  [WHISK set, L-only freq, n=%d]', ...
    numel(WHISK.absLag)));

% left: raw scatter
axW1 = subplot(1,2,1);
hold(axW1,'on');  grid(axW1,'on');

scatter(axW1, WHISK.wfreq_Lonly, WHISK.absLag, 8, ...
    'filled', ...
    'MarkerFaceColor', [0.25 0.25 0.25], ...
    'MarkerFaceAlpha', 0.15, ...
    'MarkerEdgeAlpha', 0.15);

if SHOW_FIT
    pw    = polyfit(WHISK.wfreq_Lonly, WHISK.absLag, 1);
    xfitW = linspace(WHISK_EDGES(1), WHISK_EDGES(end), 200);
    plot(axW1, xfitW, polyval(pw, xfitW), 'r-', 'LineWidth', 1.5);
    [rW, pW] = corr(WHISK.wfreq_Lonly, WHISK.absLag, 'Type', 'Pearson');
    text(axW1, 0.04, 0.93, sprintf('R = %.3f,  p = %.2e', rW, pW), ...
        'Units','normalized', 'FontSize', 9, 'Color', [0.75 0 0]);
end

yline(axW1, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.2);
xline(axW1, SLOW_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
xline(axW1, FAST_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

xlim(axW1, [WHISK_EDGES(1) WHISK_EDGES(end)]);
ylim(axW1, [0 LAG_YMAX_MS]);
xlabel(axW1, 'whisking frequency (Hz)');
ylabel(axW1, '|L-R protraction peak lag| (ms)');
title(axW1, 'Raw paired events');

% right: binned median
axW2 = subplot(1,2,2);
hold(axW2,'on');  grid(axW2,'on');

okW  = Bw.n >= MIN_BIN_N;
xwOk = whiskCenters(okW);

errorbar(axW2, xwOk, Bw.med(okW), ...
    Bw.med(okW) - Bw.q1(okW), ...
    Bw.q3(okW)  - Bw.med(okW), ...
    'o-', ...
    'Color',           [0.10 0.65 0.30], ...
    'LineWidth',       1.8, ...
    'CapSize',         5, ...
    'MarkerFaceColor', [0.10 0.65 0.30], ...
    'MarkerSize',      7);

yline(axW2, SYNC_THRESH_MS, 'k--', 'LineWidth', 1.2);
xline(axW2, SLOW_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
xline(axW2, FAST_HZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

xlim(axW2, [WHISK_EDGES(1) WHISK_EDGES(end)]);
ylim(axW2, [0 LAG_YMAX_MS]);
xlabel(axW2, 'whisking frequency (Hz)');
ylabel(axW2, 'median |L-R protraction peak lag| (ms)');
title(axW2, sprintf('Binned (1 Hz bins, min n=%d, bars=Q1–Q3)', MIN_BIN_N));
subtitle(axW2, sprintf('dashed = slow/fast boundaries | sync \\leq %.1f ms', SYNC_THRESH_MS));

% ============================ FIGURE 3: 2D HEATMAP  [BREATH set, L-only wfreq] ============================

figure('Color','w','Position',[180 60 1200 520]);
sgtitle(sprintf('Joint: median |L-R lag| (ms)  breath × whisk  [BREATH set, n=%d, bins\\geq%d]', ...
    numel(BREATH.absLag), MIN_BIN_N));

% --- median lag map (left) ---
ax1 = subplot(1,2,1);

h1 = imagesc(ax1, breathCenters, whiskCenters, medMap);
set(h1, 'AlphaData', ~isnan(medMap));
set(ax1, 'Color', [0.88 0.88 0.88]);
set(ax1, 'YDir', 'normal');
axis(ax1, 'tight');

colormap(ax1, flipud(parula));
%caxis(ax1, [0 20]);
cb1 = colorbar(ax1);
cb1.Label.String = 'median |L-R lag| (ms)';

hold(ax1,'on');
xline(ax1, BASAL_HZ, 'w-',  'LineWidth', 1.2);
xline(ax1, SNIFF_HZ, 'w-',  'LineWidth', 1.2);
yline(ax1, SLOW_HZ,  'w--', 'LineWidth', 1.2);
yline(ax1, FAST_HZ,  'w--', 'LineWidth', 1.2);
hold(ax1,'off');

xlabel(ax1, 'breathing frequency (Hz)');
ylabel(ax1, 'whisking frequency (Hz)');
title(ax1, 'Median |L-R protraction peak lag|');

% --- event count map (right) ---
ax2 = subplot(1,2,2);

nMapDisp = double(nMap);
nMapDisp(nMap == 0) = NaN;

h2 = imagesc(ax2, breathCenters, whiskCenters, nMapDisp);
set(h2, 'AlphaData', ~isnan(nMapDisp));
set(ax2, 'Color', [0.88 0.88 0.88]);
set(ax2, 'YDir', 'normal');
axis(ax2, 'tight');

colormap(ax2, hot);
cb2 = colorbar(ax2);
cb2.Label.String = 'n paired cycles';

hold(ax2,'on');
xline(ax2, BASAL_HZ, 'w-',  'LineWidth', 1.2);
xline(ax2, SNIFF_HZ, 'w-',  'LineWidth', 1.2);
yline(ax2, SLOW_HZ,  'w--', 'LineWidth', 1.2);
yline(ax2, FAST_HZ,  'w--', 'LineWidth', 1.2);
hold(ax2,'off');

xlabel(ax2, 'breathing frequency (Hz)');
ylabel(ax2, 'whisking frequency (Hz)');
title(ax2, 'Bin count');

% ============================ HELPERS ============================

function B = binnedStats(x, y, edges, syncThresh)
    nc    = numel(edges) - 1;
    B.n   = zeros(nc,1);
    B.med = nan(nc,1);
    B.q1  = nan(nc,1);
    B.q3  = nan(nc,1);
    B.iqr = nan(nc,1);
    B.sync= nan(nc,1);

    for ii = 1:nc
        mask = x >= edges(ii) & x < edges(ii+1);
        B.n(ii) = sum(mask);
        if B.n(ii) > 0
            v = y(mask);
            B.med(ii)  = median(v,'omitnan');
            B.q1(ii)   = prctile(v,25);
            B.q3(ii)   = prctile(v,75);
            B.iqr(ii)  = B.q3(ii) - B.q1(ii);
            B.sync(ii) = 100 * mean(v <= syncThresh,'omitnan');
        end
    end
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
            per  = (dprev + dnext) / 2;
            f(k) = 1 / per;
        end
    end
end

function p = subsamp(x, idx)
    x   = x(:);
    idx = double(idx(:));
    p   = idx;

    in = idx > 1 & idx < numel(x);
    i  = idx(in);

    ym = x(i-1);
    y0 = x(i);
    yp = x(i+1);

    den   = ym - 2*y0 + yp;
    delta = 0.5*(ym - yp) ./ den;
    delta(~isfinite(delta) | abs(delta) > 0.5) = 0;

    p(in) = i + delta;
end

function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac * prctile(env,95);

    d = diff([0; a; 0]);
    s = find(d ==  1);
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

    bn    = arrayfun(@(x) bestnum(x.name), d);
    [~,ix]= max(bn);

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
