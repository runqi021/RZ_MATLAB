% interv_whisk_detect_QC_RZ.m
% Detect and classify inspiratory vs intervening whisk protraction onsets.
% Produces one 5-panel QC figure — no exports.
%
% Panel 1: example raw traces with event markers
% Panel 2: breath-cycle raster (key panel)
% Panel 3: histogram of onset time relative to inspiration
% Panel 4: # intervening whisks per cycle vs breath frequency
% Panel 5: counts summary bar

% ========================= USER-EDITABLE =========================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

WHISK_EVENT_MODE = "meanLR";  % "meanLR", "left", or "right"

BP         = [5 30];      % whisk bandpass (Hz)
BR_BP      = [2 15];      % breath bandpass (Hz)
PROT_PROM  = 1;           % protraction peak prominence (× std of BP whisk)
MINDIST    = 0.025;       % s, min distance between protraction peaks
BR_PROM    = 0.5;         % inspiration trough prominence (× std of BP breath)

THR_FRAC   = 0.2;         % epoch threshold = THR_FRAC × 95th-pct(envelope)
MIN_DUR    = 0.5;         % s, min epoch duration
MERGE_GAP  = 0.20;        % s, merge epochs separated by less than this

INSP_WHISK_WIN   = [0 0.15];  % s after inspiration onset → inspiratory whisk window
MIN_INTERV_DELAY = 0.05;      % s, intervening whisk must be at least this far after insp onset
MAX_ONSET_TO_PEAK = 0.08;     % s, reject if trough-to-peak too long
MIN_ONSET_TO_PEAK = 0.005;    % s, reject if trough-to-peak too short

EXAMPLE_WIN  = 8;         % s, duration of example trace panel 1
EXAMPLE_SESS = 1;         % which session (1-indexed) to show in panel 1

SAVE_FIGS  = false;
EXCLUDE    = "5840027";
fpsW       = 400;
% =================================================================

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

% ---- session discovery ----
ad = dir(char(dataRoot)); sessAll = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    if any(strcmp(ad(a).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*'));
    [~,o] = sort({rr.name}); rr = rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), ...
                sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sessAll{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sessAll), 'no sessions found in %s', dataRoot);

% ---- pooled event containers ----
events.tInspWhisk       = [];
events.tIntervWhisk     = [];
events.tUnclassifiedWhisk = [];
events.tInsp            = [];
events.breathCycleID    = [];   % global cycle index
events.relTimeToInsp    = [];   % for each onset (insp + interv + unclassified)
events.relType          = [];   % 1=insp, 2=interv, 3=unclassified
events.breathFreq       = [];   % per-cycle breath freq
events.whiskFreq        = [];   % per-onset local whisk freq
events.sessionID        = [];   % which session each cycle belongs to
events.cycleDuration    = [];   % per-cycle duration

% per-cycle aggregates (one row per breath cycle)
cycBreathFreq  = [];
cycNInterv     = [];
cycHasInsp     = [];
cycNextInsp    = [];   % time of next inspiration (for raster end marker)
cycRelTimes    = {};   % {cycIdx} = vector of rel onset times with type labels

globalCycleID = 0;
exampleData   = [];   % store one session for panel 1

% ---- main loop ----
for e = 1:numel(sessAll)
    animal = sessAll{e}{1}; kk = sessAll{e}{2};

    % --- whisk DLC ---
    try
        M = dlc_gate_interp(pick_csv(whiskDir, ...
            sprintf('%s_whisk_n%d', animal, kk)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kk, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)' / fpsW;

    xL = filtfilt(bw,aw, fillmissing(La - mean(La,'omitnan'), 'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra - mean(Ra,'omitnan'), 'linear'));

    switch lower(char(WHISK_EVENT_MODE))
        case 'left',   xW = xL;
        case 'right',  xW = xR;
        otherwise,     xW = mean([xL xR], 2, 'omitnan');
    end

    % --- epoch mask ---
    env = (abs(hilbert(xL)) + abs(hilbert(xR))) / 2;
    ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
    epochMask = false(size(t));
    for q = 1:size(ep,1)
        epochMask = epochMask | (t >= ep(q,1) & t <= ep(q,2));
    end

    % --- breath ---
    try
        Pn = thermal_resolve_paths( ...
            pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);
        if ~isfile(Pn.breath), continue; end
        Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
        tBr = (0:numel(br)-1)' / fb;
        [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
        brf = filtfilt(b2,a2, fillmissing(br,'linear'));
        [~,iloc] = findpeaks(-brf, ...
            'MinPeakProminence', BR_PROM*std(brf), ...
            'MinPeakDistance',   round(fb/BR_BP(2)));
        if numel(iloc) < 3, continue; end
        tInspFull = (subsamp(-brf, iloc) - 1) / fb;
        tInspFull = tInspFull(:);
        brwFull   = interp1(tBr, brf, t, 'linear', NaN);
    catch ME
        warning('breath failed %s n%d: %s', animal, kk, ME.message);
        continue;
    end

    % --- protraction peaks + onsets from xW ---
    stdW = std(xW(epochMask | true));  % use full signal std for prominence
    [~,iPeak] = findpeaks(xW, ...
        'MinPeakProminence', PROT_PROM * stdW, ...
        'MinPeakDistance',   round(MINDIST * fpsW));

    % troughs = peaks of -xW
    [~,iTrough] = findpeaks(-xW, ...
        'MinPeakProminence', PROT_PROM * stdW * 0.5, ...
        'MinPeakDistance',   round(MINDIST * fpsW));

    % for each protraction peak, find nearest preceding trough
    tPeak   = subsamp(xW,  iPeak)   / fpsW;   % sub-frame peak time (seconds)
    tTrough = subsamp(-xW, iTrough) / fpsW;   % sub-frame trough time
    tPeak   = tPeak(:); tTrough = tTrough(:);

    tOnset = nan(numel(tPeak), 1);
    for pk = 1:numel(tPeak)
        % preceding troughs
        prec = tTrough(tTrough < tPeak(pk));
        if isempty(prec), continue; end
        dt = tPeak(pk) - prec(end);
        if dt < MIN_ONSET_TO_PEAK || dt > MAX_ONSET_TO_PEAK, continue; end
        onsetFrame = round(prec(end) * fpsW) + 1;
        if onsetFrame < 1 || onsetFrame > numel(t), continue; end
        if ~epochMask(onsetFrame), continue; end
        tOnset(pk) = prec(end);
    end

    validOnset = ~isnan(tOnset);
    tOnsetV    = tOnset(validOnset);
    tPeakV     = tPeak(validOnset);

    if isempty(tOnsetV), continue; end

    % local whisk frequency at each onset (centered 2-onset interval)
    fwOnset = nan(numel(tOnsetV), 1);
    for k = 2:numel(tOnsetV)-1
        dt = tOnsetV(k+1) - tOnsetV(k-1);
        if dt > 0 && dt <= 0.5
            fwOnset(k) = 2 / dt;
        end
    end

    % --- classify per breath cycle ---
    for cy = 1:numel(tInspFull)-1
        t0 = tInspFull(cy);
        t1 = tInspFull(cy+1);
        dur = t1 - t0;
        if dur <= 0 || dur > 1.5, continue; end   % skip pathological cycles
        bfreq = 1 / dur;

        % onsets in this cycle
        inCycle = tOnsetV >= t0 & tOnsetV < t1;
        if ~any(inCycle), continue; end

        globalCycleID = globalCycleID + 1;
        onsTimes = tOnsetV(inCycle);
        onsFreq  = fwOnset(inCycle);

        % sort by time (should already be, but be safe)
        [onsTimes, so] = sort(onsTimes);
        onsFreq = onsFreq(so);

        relT = onsTimes - t0;

        % find inspiratory whisk: first onset inside INSP_WHISK_WIN
        inspWin = relT >= INSP_WHISK_WIN(1) & relT <= INSP_WHISK_WIN(2);
        hasInsp = false;
        inspIdx = 0;
        if any(inspWin)
            inspIdx = find(inspWin, 1, 'first');
            hasInsp = true;
        end

        for j = 1:numel(onsTimes)
            absT = onsTimes(j);
            rel  = relT(j);
            if hasInsp && j == inspIdx
                typ = 1;  % inspiratory
                events.tInspWhisk(end+1,1)    = absT;
            elseif hasInsp && j > inspIdx && rel >= MIN_INTERV_DELAY
                typ = 2;  % intervening
                events.tIntervWhisk(end+1,1)  = absT;
            else
                typ = 3;  % unclassified
                events.tUnclassifiedWhisk(end+1,1) = absT;
            end
            events.relTimeToInsp(end+1,1)  = rel;
            events.relType(end+1,1)        = typ;
            events.breathCycleID(end+1,1) = globalCycleID;
            events.whiskFreq(end+1,1)     = onsFreq(j);
        end

        events.tInsp(end+1,1)      = t0;
        events.breathFreq(end+1,1) = bfreq;
        events.sessionID(end+1,1)  = e;
        events.cycleDuration(end+1,1) = dur;

        cycBreathFreq(end+1,1) = bfreq; %#ok<AGROW>
        nInterv = sum(events.relType(events.breathCycleID == globalCycleID) == 2);
        cycNInterv(end+1,1)   = nInterv; %#ok<AGROW>
        cycHasInsp(end+1,1)   = hasInsp; %#ok<AGROW>
        cycNextInsp(end+1,1)  = dur;     %#ok<AGROW>   relative duration
    end

    % store first valid session for panel 1
    if isempty(exampleData) && e == EXAMPLE_SESS && ~isempty(tInspFull)
        exampleData.t        = t;
        exampleData.xW       = xW;
        exampleData.brw      = brwFull;
        exampleData.tInsp    = tInspFull;
        exampleData.tOnsetV  = tOnsetV;
        exampleData.sessIdx  = e;
        exampleData.animal   = animal;
        exampleData.kk       = kk;
    end
end

assert(~isempty(events.tInsp), 'no breath cycles processed — check data paths');

% ---- summary counts ----
nCycles  = numel(unique(events.breathCycleID));
nInsp    = numel(events.tInspWhisk);
nInterv  = numel(events.tIntervWhisk);
nUnclass = numel(events.tUnclassifiedWhisk);

% cycles that have ≥1 intervening whisk
cycWithInterv = sum(cycNInterv >= 1);
fracInterv    = 100 * cycWithInterv / numel(cycNInterv);

% median breath freq for cycles WITH intervening whisks
medFreqInterv = median(cycBreathFreq(cycNInterv >= 1), 'omitnan');

fprintf('\n=== Intervening whisk detection ===\n');
fprintf('Sessions processed: %d\n', numel(sessAll));
fprintf('Detected %d inspiratory whisks, %d intervening whisks, %d unclassified whisks across %d breath cycles\n', ...
    nInsp, nInterv, nUnclass, nCycles);
fprintf('Cycles with >=1 intervening whisk: %.1f%%\n', fracInterv);
fprintf('Median breath freq for cycles with intervening whisks: %.2f Hz\n', medFreqInterv);

% ================================================================
%  5-panel QC figure
% ================================================================
cGreen  = [0.13 0.60 0.20];
cMag    = [0.85 0.10 0.55];
cGray   = [0.60 0.60 0.60];
cBlue   = [0.15 0.35 0.75];

fig = figure('Name','Intervening whisk detection QC', ...
    'Color','w', 'Position',[30 30 1300 1000]);
tlo = tiledlayout(3, 2, 'TileSpacing','compact', 'Padding','compact');
suptitle_str = sprintf('Intervening whisk QC  |  %d cycles  |  insp=%d  interv=%d  unclass=%d', ...
    nCycles, nInsp, nInterv, nUnclass);
title(tlo, suptitle_str, 'FontSize',9, 'FontWeight','normal');

% ----------------------------------------------------------------
% Panel 1: example traces  (spans top row, left)
% ----------------------------------------------------------------
ax1 = nexttile(tlo, [1 1]);
if ~isempty(exampleData)
    t0ex = exampleData.tInsp(1) - 0.3;
    t0ex = max(t0ex, 0);
    t1ex = t0ex + EXAMPLE_WIN;
    tW   = exampleData.t;
    xWex = exampleData.xW;
    brex = exampleData.brw;

    tmask = tW >= t0ex & tW <= t1ex;
    tSeg  = tW(tmask) - t0ex;

    % normalize each signal independently to [−1 +1]
    xWseg = xWex(tmask);
    brseg = brex(tmask);
    xWn = normSig(xWseg);
    brn  = normSig(brseg) - 2.5;  % offset below whisk

    hold(ax1,'on');
    plot(ax1, tSeg, xWn,  '-', 'Color',[0.3 0.3 0.3], 'LineWidth',0.8);
    plot(ax1, tSeg, brn,  '-', 'Color',[0.2 0.5 0.9], 'LineWidth',0.8);

    % inspiration onsets
    tInspEx = exampleData.tInsp;
    tInspEx = tInspEx(tInspEx >= t0ex & tInspEx <= t1ex) - t0ex;
    for qi = 1:numel(tInspEx)
        xline(ax1, tInspEx(qi), 'b-', 'Alpha',0.4, 'LineWidth',1.0);
    end

    % classify onsets in this window
    tOnEx = exampleData.tOnsetV;
    tOnEx = tOnEx(tOnEx >= t0ex & tOnEx <= t1ex) - t0ex;
    tInspFull_ex = exampleData.tInsp;

    for qi = 1:numel(tOnEx)
        absT = tOnEx(qi) + t0ex;
        typ = classifyOnset(absT, tInspFull_ex, INSP_WHISK_WIN, MIN_INTERV_DELAY);
        switch typ
            case 1, col = cGreen;
            case 2, col = cMag;
            otherwise, col = cGray;
        end
        plot(ax1, tOnEx(qi), 0.5 + rand()*0.3, 'v', ...
            'Color',col, 'MarkerFaceColor',col, 'MarkerSize',6);
    end

    ylim(ax1, [-4 2]);
    xlabel(ax1,'time (s)'); ylabel(ax1,'normalized');
    title(ax1, sprintf('example traces — %s n%d  (whisk mode: %s)', ...
        exampleData.animal, exampleData.kk, WHISK_EVENT_MODE));
    legend(ax1, {'whisk xW','breath (offset)'}, 'Location','northeast','Box','off');
    box(ax1,'off'); grid(ax1,'on');
else
    text(0.5,0.5,'no example data','HorizontalAlignment','center','Parent',ax1);
end

% ----------------------------------------------------------------
% Panel 2: breath-cycle raster  (spans 2 rows, right column)
% ----------------------------------------------------------------
ax2 = nexttile(tlo, 2, [2 1]);   % tile 2, spanning 2 rows in col 2

% sort cycles by breath duration (ascending = fast breaths at top)
[~,sortOrd] = sort(cycNextInsp, 'descend');   % slow at top → descending dur
nCyc = numel(sortOrd);
cycYmap = zeros(1, nCyc);
for yi = 1:nCyc
    cycYmap(sortOrd(yi)) = yi;
end

hold(ax2,'on');
% gray line at cycle end (next inspiration)
for ci = 1:numel(sortOrd)
    origCy = sortOrd(ci);
    yy = ci;
    dur = cycNextInsp(origCy);
    plot(ax2, dur, yy, '|', 'Color',[0.85 0.85 0.85], 'MarkerSize',5, 'LineWidth',1);
end

% plot events
allCycIDs = unique(events.breathCycleID);
for ci = 1:numel(allCycIDs)
    cid = allCycIDs(ci);
    mask = events.breathCycleID == cid;
    yy   = cycYmap(ci);  % ci == cycle ordering after sort
    rels  = events.relTimeToInsp(mask);
    typs  = events.relType(mask);
    for qi = 1:numel(rels)
        switch typs(qi)
            case 1, col = cGreen; mk = 'o'; sz = 5;
            case 2, col = cMag;   mk = 'o'; sz = 5;
            otherwise, col = cGray; mk = '.'; sz = 8;
        end
        plot(ax2, rels(qi), yy, mk, ...
            'Color',col, 'MarkerFaceColor',col, 'MarkerSize',sz);
    end
end

% legend proxy
plot(ax2, nan,nan,'o','Color',cGreen,'MarkerFaceColor',cGreen,'MarkerSize',5,'DisplayName','inspiratory');
plot(ax2, nan,nan,'o','Color',cMag,  'MarkerFaceColor',cMag,  'MarkerSize',5,'DisplayName','intervening');
plot(ax2, nan,nan,'.','Color',cGray, 'MarkerSize',8,'DisplayName','unclassified');
legend(ax2,'Location','northeast','Box','off');

xline(ax2, INSP_WHISK_WIN(2), 'g--', 'Alpha',0.5, 'LineWidth',1.0, ...
    'Label','insp win end','LabelVerticalAlignment','bottom');

xlabel(ax2,'time from inspiration onset (s)');
ylabel(ax2,'breath cycle (sorted by duration, slow→fast)');
title(ax2, sprintf('breath-cycle raster  (%d cycles)', nCyc));
xlim(ax2, [-0.05, max(cycNextInsp) + 0.05]);
ylim(ax2, [0 nCyc+1]);
box(ax2,'off'); grid(ax2,'on');

% ----------------------------------------------------------------
% Panel 3: histogram of relative onset time
% ----------------------------------------------------------------
ax3 = nexttile(tlo);

relInsp    = events.relTimeToInsp(events.relType == 1);
relInterv  = events.relTimeToInsp(events.relType == 2);
relUnclass = events.relTimeToInsp(events.relType == 3);
maxRelT    = max(cycNextInsp) + 0.02;
edges3     = 0 : 0.02 : min(maxRelT, 1.5);

hold(ax3,'on'); grid(ax3,'on');
if ~isempty(relUnclass)
    histogram(ax3, relUnclass, edges3, 'Normalization','count', ...
        'FaceColor',cGray,'EdgeColor','none','FaceAlpha',0.5);
end
if ~isempty(relInterv)
    histogram(ax3, relInterv, edges3, 'Normalization','count', ...
        'FaceColor',cMag,'EdgeColor','none','FaceAlpha',0.7);
end
if ~isempty(relInsp)
    histogram(ax3, relInsp, edges3, 'Normalization','count', ...
        'FaceColor',cGreen,'EdgeColor','none','FaceAlpha',0.9);
end
xline(ax3, INSP_WHISK_WIN(1), 'g--','LineWidth',1.0);
xline(ax3, INSP_WHISK_WIN(2), 'g--','LineWidth',1.0, ...
    'Label','insp win','LabelVerticalAlignment','bottom');
xlabel(ax3,'time from inspiration onset (s)');
ylabel(ax3,'count');
legend(ax3, {sprintf('unclassified n=%d',numel(relUnclass)), ...
             sprintf('intervening n=%d',numel(relInterv)), ...
             sprintf('inspiratory n=%d',numel(relInsp))}, ...
    'Location','northeast','Box','off');
title(ax3,'protraction onset time relative to inspiration');
box(ax3,'off');

% ----------------------------------------------------------------
% Panel 4: # intervening whisks per cycle vs breath frequency
% ----------------------------------------------------------------
ax4 = nexttile(tlo);

hold(ax4,'on'); grid(ax4,'on');
% scatter with jitter for visibility
jit = 0.04 * (rand(numel(cycBreathFreq),1) - 0.5);
scatter(ax4, cycBreathFreq, cycNInterv + jit, 10, ...
    [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha',0.4);

% mean per freq bin
bfEdges = 0:1:20;
bfCtr   = bfEdges(1:end-1) + 0.5;
mnInterv = nan(numel(bfCtr),1);
for bi = 1:numel(bfCtr)
    ok = cycBreathFreq >= bfEdges(bi) & cycBreathFreq < bfEdges(bi+1);
    if sum(ok) >= 3
        mnInterv(bi) = mean(cycNInterv(ok));
    end
end
plot(ax4, bfCtr, mnInterv, '-o', 'Color',[0.1 0.1 0.1], ...
    'LineWidth',1.6, 'MarkerSize',5, 'MarkerFaceColor',[0.1 0.1 0.1]);

xlabel(ax4,'breathing frequency (Hz)');
ylabel(ax4,'# intervening whisks in cycle');
title(ax4,'intervening whisks per cycle vs breath frequency');
xlim(ax4,[0 18]); ylim(ax4,[-0.3, max(cycNInterv)+0.5]);
box(ax4,'off');

% ----------------------------------------------------------------
% Panel 5: counts summary bar
% ----------------------------------------------------------------
ax5 = nexttile(tlo);

% turn off axis, use text table
axis(ax5,'off');
labels = { ...
    'breath cycles total', ...
    'cycles with insp whisk', ...
    'inspiratory whisk events', ...
    'intervening whisk events', ...
    'unclassified whisk events', ...
    'cycles with ≥1 interv whisk', ...
    'med. breath freq (interv cycles)'};

vals = { ...
    sprintf('%d', nCycles), ...
    sprintf('%d  (%.0f%%)', sum(cycHasInsp), 100*mean(cycHasInsp)), ...
    sprintf('%d', nInsp), ...
    sprintf('%d', nInterv), ...
    sprintf('%d', nUnclass), ...
    sprintf('%d  (%.1f%%)', cycWithInterv, fracInterv), ...
    sprintf('%.2f Hz', medFreqInterv)};

cols5 = {'k','k', cGreen, cMag, cGray, cMag, [0.3 0.3 0.6]};

for ri = 1:numel(labels)
    ypos = 1 - ri/( numel(labels)+1);
    text(ax5, 0.02, ypos, labels{ri}, 'Units','normalized', ...
        'FontSize',10, 'Color',[0.3 0.3 0.3], 'HorizontalAlignment','left');
    text(ax5, 0.98, ypos, vals{ri}, 'Units','normalized', ...
        'FontSize',10, 'Color',cols5{ri}, 'HorizontalAlignment','right', 'FontWeight','bold');
end
title(ax5,'counts summary');

% ================================================================
%  LOCAL FUNCTIONS
% ================================================================

function typ = classifyOnset(absT, tInspFull, INSP_WHISK_WIN, MIN_INTERV_DELAY)
% Classify a single onset given the global inspiration time vector.
    cy = find(tInspFull <= absT, 1, 'last');
    if isempty(cy) || cy >= numel(tInspFull)
        typ = 3; return;
    end
    t0  = tInspFull(cy);
    rel = absT - t0;
    if rel >= INSP_WHISK_WIN(1) && rel <= INSP_WHISK_WIN(2)
        typ = 1;
    elseif rel > INSP_WHISK_WIN(2) && rel >= MIN_INTERV_DELAY
        typ = 2;
    else
        typ = 3;
    end
end

function yn = normSig(x)
    x = x - mean(x,'omitnan');
    rng = max(abs(x));
    if rng > 0, yn = x / rng; else, yn = x; end
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
