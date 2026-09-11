% whisk_meta_shuffle_controls_RZ.m
%
% Meta shuffle-control script for current draft figure set.
%
% Produces:
%   Figure 1: 10-panel L/R |lag| histograms with circular-shift R-whisk null
%   Figure 2: breath-phase histograms with circular-shift whisk-event null
%   Figure 3: peri-inspiration whisk-onset PSTHs with circular-shift whisk-event null
%
% Definitions fixed by Runqi:
%   basal breathing:  <5 Hz
%   sniffing:         >8 Hz
%   slow whisking:    <10 Hz
%   fast whisking:    >20 Hz
%   L/R sync thresh:  7.7 ms
%
% L/R synchrony uses THR_FRAC_LR = 0.2.
% Breath-phase / PSTH control uses THR_FRAC_PHASE = 0.
%
% Null logic:
%   L/R synchrony null:
%       circularly shift R protraction-peak times within each session,
%       re-pair to L peaks, recompute |L-R| lag.
%
%   Breath-phase/PSTH null:
%       circularly shift whisk protraction-onset times within each session
%       relative to thermal inspiration times, recompute phase/PSTH.
%
% This does NOT change your main analysis scripts. It is a control/meta script.

clear; close all; clc;

% ========================= USER SETTINGS =========================
dataRoot   = "D:\260615_thermalNbasler";
whiskDir   = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

EXCLUDE = "5840027";
fpsW    = 400;

% shared processing
BP        = [5 30];
BR_BP     = [2 15];
MINDIST   = 0.025;
PROT_PROM = 1;
BR_PROM   = 0.5;
MAXLAG    = 0.08;

% definitions
BASAL_HZ = 5;
SNIFF_HZ = 8;
SLOW_HZ  = 10;
FAST_HZ  = 20;

SYNC_THRESH_MS = 7.7;
CYC_MAX = 0.5;

% epoch thresholds
THR_FRAC_LR    = 0.2;   % for 10-panel L/R lag figure
THR_FRAC_PHASE = 0.0;   % for whisk-breath phase/PSTH figure

MIN_DUR    = 0.5;
MERGE_GAP  = 0.20;
ENV_SMOOTH = 0;

% whisk local frequency
ITI_MAX_LR    = 0.3;
WHISK_ITI_MAX = 0.25;

% histogram settings
LAG_BIN_MS  = 2.5;

% per-panel lag display windows (ms)
LAG_WIN.all        = 80;
LAG_WIN.slow       = 80;
LAG_WIN.fast       = 80;
LAG_WIN.basal      = 80;
LAG_WIN.sniff      = 80;
LAG_WIN.basal_slow = 80;
LAG_WIN.basal_fast = 80;
LAG_WIN.sniff_slow = 80;
LAG_WIN.sniff_fast = 80;

LAG_XMAX_MS = max(struct2array(LAG_WIN));
lagEdges    = 0:LAG_BIN_MS:LAG_XMAX_MS;

itiEdges    = 0:0.01:ITI_MAX_LR;

NBINS_PH    = 31;
phaseEdges  = linspace(0,1,NBINS_PH+1);

BIN_W       = 0.005;

% per-panel display windows (seconds)
TIME_WIN.pop   = 0.10;
TIME_WIN.slowB = 0.10;
TIME_WIN.fastB = 0.10;
TIME_WIN.slowW = 0.10;
TIME_WIN.sBsW  = 0.10;
TIME_WIN.fBsW  = 0.10;
TIME_WIN.fastW = 0.10;
TIME_WIN.sBfW  = 0.10;
TIME_WIN.fBfW  = 0.10;

PAIR_WIN    = max(struct2array(TIME_WIN));
timeEdges   = -PAIR_WIN:BIN_W:PAIR_WIN;

% shuffle settings
NSHUFF   = 500;
RNG_SEED = 1;

% plotting
YLIM_LAG   = [0 0.4];
YLIM_PHASE = [];
YLIM_TIME  = [];

SAVE_FIGS = false;
OUTDIR    = fullfile(pwd, 'shuffle_control_exports');
% =================================================================

rng(RNG_SEED);

if SAVE_FIGS && ~exist(OUTDIR,'dir')
    mkdir(OUTDIR);
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

grpLag = {'all','slow','fast','basal','sniff','basal_slow','basal_fast','sniff_slow','sniff_fast'};
grpPhase = {'pop','slowB','fastB','slowW','sBsW','fBsW','fastW','sBfW','fBfW'};

% observed event containers
OBS_LAG = emptyStruct(grpLag);
OBS_PHASE = emptyStruct(grpPhase);
OBS_TIME = emptyStruct(grpPhase);

ITI.protL = [];
ITI.protR = [];

% shuffle count matrices; each row = one shuffle replicate
SH_LAG_COUNT   = initCountStruct(grpLag,   NSHUFF, numel(lagEdges)-1);
SH_PHASE_COUNT = initCountStruct(grpPhase, NSHUFF, numel(phaseEdges)-1);
SH_TIME_COUNT  = initCountStruct(grpPhase, NSHUFF, numel(timeEdges)-1);

SH_TIME_VALS = struct();
for i = 1:numel(grpPhase)
    SH_TIME_VALS.(grpPhase{i}) = cell(NSHUFF,1);
end

% ---------------- categorical state palette ----------------
C.pop   = [0.20 0.20 0.20];  % population / all
C.slowB = [0.04 0.17 0.42];  % slow breathing / basal
C.fastB = [0.94 0.29 0.14];  % fast breathing / sniffing
C.slowW = [0.48 0.25 0.70];  % slow whisking
C.fastW = [0.00 0.61 0.48];  % fast whisking
C.sBsW  = [0.13 0.72 0.90];  % slow breath + slow whisk
C.sBfW  = [0.00 0.32 0.80];  % slow breath + fast whisk
C.fBsW  = [0.69 0.07 0.39];  % fast breath + slow whisk
C.fBfW  = [0.95 0.73 0.23];  % fast breath + fast whisk

% ---------------- Figure 1 aliases: L/R lag panels ----------------
C.all        = C.pop;
C.slow       = C.slowW;
C.fast       = C.fastW;
C.basal      = C.slowB;
C.sniff      = C.fastB;
C.basal_slow = C.sBsW;
C.basal_fast = C.sBfW;
C.sniff_slow = C.fBsW;
C.sniff_fast = C.fBfW;

% session discovery
sess = discoverSessions(dataRoot, whiskDir, EXCLUDE);
fprintf('Found %d candidate whisk sessions\n', numel(sess));

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

nSessLR = 0;
nSessPhase = 0;

for e = 1:numel(sess)
    animal = sess{e}{1};
    kk     = sess{e}{2};

    % ------------------------- load whisk -------------------------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
    catch ME
        warning('Skipping %s n%d: whisk load failed: %s', animal, kk, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)' / fpsW;
    Tdur = t(end) - t(1);

    xL = filtfilt(bw,aw, fillmissing(La - mean(La,'omitnan'), 'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra - mean(Ra,'omitnan'), 'linear'));

    env = movmean((abs(hilbert(xL)) + abs(hilbert(xR))) / 2, ...
        max(1, round(ENV_SMOOTH*fpsW)));

    % ------------------------- load breath -------------------------
    hasBreath = false;
    tInsp = [];
    tPeak = [];
    cycdur = [];

    try
        Pn = thermal_resolve_paths( ...
            pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);

        if isfile(Pn.breath)
            Bs = load(Pn.breath);
            br = Bs.breath(:);
            fb = double(Bs.fps);

            [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
            brf = filtfilt(b2,a2, fillmissing(br,'linear'));

            [~,il] = findpeaks(-brf, ...
                'MinPeakProminence', BR_PROM*std(brf), ...
                'MinPeakDistance', round(fb/BR_BP(2)));

            if numel(il) >= 2
                tInsp = (subsamp(-brf, il) - 1) / fb;
                tInsp = tInsp(:);
                cycdur = diff(tInsp);

                tPeak = nan(numel(il)-1,1);
                for k = 1:numel(il)-1
                    [~,rel] = max(brf(il(k):il(k+1)));
                    tPeak(k) = (subsamp(brf, il(k)+rel-1) - 1) / fb;
                end

                hasBreath = ~isempty(cycdur);
            end
        end
    catch ME
        warning('Breath load failed %s n%d: %s', animal, kk, ME.message);
    end

    % ===============================================================
    % PART A — L/R synchrony, THR_FRAC_LR = 0.2
    % ===============================================================
    epLR = detect(env, t, THR_FRAC_LR, MIN_DUR, MERGE_GAP);
    mLR = epochMask(t, epLR);

    if any(mLR)
        [~,iL] = findpeaks(xL, ...
            'MinPeakProminence', PROT_PROM*std(xL), ...
            'MinPeakDistance', round(MINDIST*fpsW));
        [~,iR] = findpeaks(xR, ...
            'MinPeakProminence', PROT_PROM*std(xR), ...
            'MinPeakDistance', round(MINDIST*fpsW));

        iL = iL(mLR(iL));
        iR = iR(mLR(iR));

        if numel(iL) >= 3 && numel(iR) >= 3
            tL = (subsamp(xL, iL) - 1) / fpsW;
            tR = (subsamp(xR, iR) - 1) / fpsW;
            tL = tL(:); tR = tR(:);

            fL = centeredFreqAtPeaks(tL, ITI_MAX_LR);

            % observed L/R pair table
            E = pairLR(tL, tR, fL, MAXLAG, tInsp, cycdur, CYC_MAX);
            OBS_LAG = appendLagGroups(OBS_LAG, E, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ);

            % ITI from protraction onsets, same as 10-panel script
            ITI.protL = [ITI.protL; iti(-xL, mLR, t, PROT_PROM, fpsW, MINDIST, ITI_MAX_LR)];
            ITI.protR = [ITI.protR; iti(-xR, mLR, t, PROT_PROM, fpsW, MINDIST, ITI_MAX_LR)];

            % shuffled R-peak controls
            for s = 1:NSHUFF
                off = randomCircularOffset(Tdur);
                tRsh = circShiftTimes(tR, off, t(1), t(end));

                Esh = pairLR(tL, tRsh, fL, MAXLAG, tInsp, cycdur, CYC_MAX);
                SH_LAG_COUNT = addLagShuffleCounts(SH_LAG_COUNT, Esh, lagEdges, s, ...
                    SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ);
            end

            nSessLR = nSessLR + 1;
        end
    end

    % ===============================================================
    % PART B — whisk-breath phase/PSTH, THR_FRAC_PHASE = 0
    % ===============================================================
    if hasBreath
        epPH = detect(env, t, THR_FRAC_PHASE, MIN_DUR, MERGE_GAP);
        mPH = epochMask(t, epPH);

        if any(mPH)
            [~,oL] = findpeaks(-xL, ...
                'MinPeakProminence', PROT_PROM*std(xL), ...
                'MinPeakDistance', round(MINDIST*fpsW));
            [~,oR] = findpeaks(-xR, ...
                'MinPeakProminence', PROT_PROM*std(xR), ...
                'MinPeakDistance', round(MINDIST*fpsW));

            oL = oL(mPH(oL));
            oR = oR(mPH(oR));

            [twL,fwL] = onsetFreq((subsamp(-xL,oL)-1)/fpsW, WHISK_ITI_MAX);
            [twR,fwR] = onsetFreq((subsamp(-xR,oR)-1)/fpsW, WHISK_ITI_MAX);

            tw = [twL; twR];
            fw = [fwL; fwR];

            good = isfinite(tw) & isfinite(fw);
            tw = tw(good);
            fw = fw(good);

            if ~isempty(tw)
                % observed
                [ph, fwk, brk] = coordPhaseFB_events(tw, fw, tInsp, tPeak, CYC_MAX);
                OBS_PHASE = appendPhaseGroups(OBS_PHASE, ph, fwk, brk, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ);

                Dlag = psthGroups(tw, fw, tInsp, cycdur, CYC_MAX, PAIR_WIN, ...
                    SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ);
                OBS_TIME = appendTimeGroups(OBS_TIME, Dlag);

                % shuffled whisk-onset controls
                for s = 1:NSHUFF
                    off = randomCircularOffset(Tdur);
                    twSh = circShiftTimes(tw, off, t(1), t(end));

                    [phSh, fwSh, brSh] = coordPhaseFB_events(twSh, fw, tInsp, tPeak, CYC_MAX);
                    SH_PHASE_COUNT = addPhaseShuffleCounts(SH_PHASE_COUNT, phSh, fwSh, brSh, ...
                        phaseEdges, s, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ);

                    DlagSh = psthGroups(twSh, fw, tInsp, cycdur, CYC_MAX, PAIR_WIN, ...
                        SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ);
                    SH_TIME_COUNT = addTimeShuffleCounts(SH_TIME_COUNT, DlagSh, timeEdges, s);

                    fn = fieldnames(DlagSh);
                    for ii = 1:numel(fn)
                        g = fn{ii};
                        SH_TIME_VALS.(g){s} = [SH_TIME_VALS.(g){s}; DlagSh.(g)];
                    end
                end

                nSessPhase = nSessPhase + 1;
            end
        end
    end
end

fprintf('\nMeta shuffle control complete.\n');
fprintf('L/R sessions used: %d\n', nSessLR);
fprintf('Phase/PSTH sessions used: %d\n', nSessPhase);
fprintf('NSHUFF = %d\n\n', NSHUFF);

% ============================ FIGURE 1 ============================
fig1 = figure('Color','w','Position',[30 30 1380 1020], ...
    'Name','L/R lag histograms with R-shift null', ...
    'Renderer','painters');
tlo = tiledlayout(fig1,5,2,'TileSpacing','compact','Padding','compact');
title(tlo, sprintf('L/R synchrony with circular-shift R-whisk null | sync threshold %.1f ms | %d shuffles', ...
    SYNC_THRESH_MS, NSHUFF), 'FontSize', 9, 'FontWeight','normal');

ax = nexttile(tlo);
plotITIHist(ax, ITI.protL, ITI.protR, itiEdges, [0.20 0.40 0.80], [0.85 0.10 0.10], [0 ITI_MAX_LR], [0 0.4]);

plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'all',        lagEdges, LAG_WIN.all,        C.all,        'all pairs',        SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'slow',       lagEdges, LAG_WIN.slow,       C.slow,       sprintf('slow whisk <%.0f Hz',SLOW_HZ), SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'fast',       lagEdges, LAG_WIN.fast,       C.fast,       sprintf('fast whisk >%.0f Hz',FAST_HZ), SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'basal',      lagEdges, LAG_WIN.basal,      C.basal,      sprintf('basal breath <%.0f Hz',BASAL_HZ), SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'sniff',      lagEdges, LAG_WIN.sniff,      C.sniff,      sprintf('sniff breath >%.0f Hz',SNIFF_HZ), SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'basal_slow', lagEdges, LAG_WIN.basal_slow, C.basal_slow, 'basal + slow',     SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'basal_fast', lagEdges, LAG_WIN.basal_fast, C.basal_fast, 'basal + fast',     SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'sniff_slow', lagEdges, LAG_WIN.sniff_slow, C.sniff_slow, 'sniff + slow',     SYNC_THRESH_MS, YLIM_LAG);
plotLagPanel(tlo, OBS_LAG, SH_LAG_COUNT, 'sniff_fast', lagEdges, LAG_WIN.sniff_fast, C.sniff_fast, 'sniff + fast',     SYNC_THRESH_MS, YLIM_LAG);

% ============================ FIGURE 2 ============================
fig2 = figure('Color','w','Position',[80 80 1250 900], ...
    'Name','Breath phase histograms with circular-shift whisk null', ...
    'Renderer','painters');
tlo2 = tiledlayout(fig2,3,3,'TileSpacing','compact','Padding','compact');
title(tlo2, sprintf('Whisk protraction-onset breath phase with circular-shift null | %d shuffles', NSHUFF), ...
    'FontSize',9,'FontWeight','normal');

phaseOrder = {'pop','slowB','fastB','slowW','sBsW','fBsW','fastW','sBfW','fBfW'};
phaseLabels = {'population','slow breathing','fast breathing', ...
               'slow whisk','slowB + slowW','fastB + slowW', ...
               'fast whisk','slowB + fastW','fastB + fastW'};

for i = 1:numel(phaseOrder)
    g = phaseOrder{i};
    ax = nexttile(tlo2);
    plotObsNullProb(ax, OBS_PHASE.(g), SH_PHASE_COUNT.(g), phaseEdges, C.(g), phaseLabels{i}, ...
        'breath phase', 'P(phase)', YLIM_PHASE);
    xlim(ax,[0 1]);
    xline(ax,0,'k-','LineWidth',0.8);
    xline(ax,0.5,'k:','LineWidth',0.8);
    xline(ax,1,'k-','LineWidth',0.8);
end

% ============================ FIGURE 3 ============================
fig3 = figure('Color','w','Position',[120 120 1250 900], ...
    'Name','Peri-inspiration PSTH with circular-shift whisk null', ...
    'Renderer','painters');
tlo3 = tiledlayout(fig3,3,3,'TileSpacing','compact','Padding','compact');
title(tlo3, sprintf('Whisk protraction-onset timing relative to inspiration with circular-shift null | %d shuffles', NSHUFF), ...
    'FontSize',9,'FontWeight','normal');

for i = 1:numel(phaseOrder)
    g = phaseOrder{i};
    ax = nexttile(tlo3);

    win      = TIME_WIN.(g);
    edgesNow = -win:BIN_W:win;

    obsNow = OBS_TIME.(g);
    obsNow = obsNow(abs(obsNow) <= win);

    shCountsNow = zeros(NSHUFF, numel(edgesNow)-1);
    for s = 1:NSHUFF
        shVals = SH_TIME_VALS.(g){s};
        shVals = shVals(abs(shVals) <= win);
        shCountsNow(s,:) = histcounts(shVals, edgesNow);
    end

    plotObsNullProb(ax, obsNow, shCountsNow, edgesNow, C.(g), phaseLabels{i}, ...
        'time from inspiration onset (s)', 'P(event lag)', YLIM_TIME);

    xlim(ax,[-win win]);
    xline(ax,0,'k--','LineWidth',1.0);
end

if SAVE_FIGS
    exportVectorFig(fig1, fullfile(OUTDIR, 'shuffle_LR_lag_10panel'));
    exportVectorFig(fig2, fullfile(OUTDIR, 'shuffle_breath_phase'));
    exportVectorFig(fig3, fullfile(OUTDIR, 'shuffle_peri_insp_PSTH'));
end

% ======================================================================
% LOCAL FUNCTIONS
% ======================================================================

function sess = discoverSessions(dataRoot, whiskDir, EXCLUDE)
    ad = dir(char(dataRoot));
    sess = {};
    for a = 1:numel(ad)
        if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
        if any(strcmp(ad(a).name, EXCLUDE)), continue; end
        rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*'));
        [~,o] = sort({rr.name});
        rr = rr(o);
        for kk = 1:numel(rr)
            if ~isempty(dir(fullfile(char(whiskDir), ...
                    sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
                sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
            end
        end
    end
end

function S = emptyStruct(names)
    S = struct();
    for i = 1:numel(names)
        S.(names{i}) = [];
    end
end

function S = initCountStruct(names, nshuff, nbins)
    S = struct();
    for i = 1:numel(names)
        S.(names{i}) = zeros(nshuff, nbins);
    end
end

function m = epochMask(t, ep)
    m = false(size(t));
    for q = 1:size(ep,1)
        m = m | (t >= ep(q,1) & t <= ep(q,2));
    end
end

function E = pairLR(tL, tR, fL, MAXLAG, tInsp, cycdur, CYC_MAX)
    E.absLag = [];
    E.signedLag = [];
    E.wfreq = [];
    E.brate = [];
    E.hasBreath = [];

    hasBreath = ~isempty(tInsp) && ~isempty(cycdur);

    for k = 2:numel(tL)-1
        wfreq = fL(k);
        if isnan(wfreq), continue; end

        [d,j] = min(abs(tR - tL(k)));
        if d > MAXLAG, continue; end

        tPair = 0.5*(tL(k)+tR(j));
        absLag = abs(tL(k)-tR(j))*1000;
        signedLag = (tL(k)-tR(j))*1000;

        brate = NaN;
        hB = false;
        if hasBreath
            bn = discretize(tPair, tInsp);
            if ~isnan(bn) && bn >= 1 && bn <= numel(cycdur)
                bdur = cycdur(bn);
                if bdur > 0 && bdur <= CYC_MAX
                    brate = 1/bdur;
                    hB = true;
                end
            end
        end

        E.absLag(end+1,1) = absLag; %#ok<AGROW>
        E.signedLag(end+1,1) = signedLag; %#ok<AGROW>
        E.wfreq(end+1,1) = wfreq; %#ok<AGROW>
        E.brate(end+1,1) = brate; %#ok<AGROW>
        E.hasBreath(end+1,1) = hB; %#ok<AGROW>
    end
end

function OBS = appendLagGroups(OBS, E, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ)
    if isempty(E.absLag), return; end

    isSlow = E.wfreq < SLOW_HZ;
    isFast = E.wfreq > FAST_HZ;
    isBasal = E.hasBreath & E.brate < BASAL_HZ;
    isSniff = E.hasBreath & E.brate > SNIFF_HZ;

    OBS.all        = [OBS.all;        E.absLag];
    OBS.slow       = [OBS.slow;       E.absLag(isSlow)];
    OBS.fast       = [OBS.fast;       E.absLag(isFast)];
    OBS.basal      = [OBS.basal;      E.absLag(isBasal)];
    OBS.sniff      = [OBS.sniff;      E.absLag(isSniff)];
    OBS.basal_slow = [OBS.basal_slow; E.absLag(isBasal & isSlow)];
    OBS.basal_fast = [OBS.basal_fast; E.absLag(isBasal & isFast)];
    OBS.sniff_slow = [OBS.sniff_slow; E.absLag(isSniff & isSlow)];
    OBS.sniff_fast = [OBS.sniff_fast; E.absLag(isSniff & isFast)];
end

function SH = addLagShuffleCounts(SH, E, edges, s, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ)
    if isempty(E.absLag), return; end

    isSlow = E.wfreq < SLOW_HZ;
    isFast = E.wfreq > FAST_HZ;
    isBasal = E.hasBreath & E.brate < BASAL_HZ;
    isSniff = E.hasBreath & E.brate > SNIFF_HZ;

    SH.all(s,:)        = SH.all(s,:)        + histcounts(E.absLag, edges);
    SH.slow(s,:)       = SH.slow(s,:)       + histcounts(E.absLag(isSlow), edges);
    SH.fast(s,:)       = SH.fast(s,:)       + histcounts(E.absLag(isFast), edges);
    SH.basal(s,:)      = SH.basal(s,:)      + histcounts(E.absLag(isBasal), edges);
    SH.sniff(s,:)      = SH.sniff(s,:)      + histcounts(E.absLag(isSniff), edges);
    SH.basal_slow(s,:) = SH.basal_slow(s,:) + histcounts(E.absLag(isBasal & isSlow), edges);
    SH.basal_fast(s,:) = SH.basal_fast(s,:) + histcounts(E.absLag(isBasal & isFast), edges);
    SH.sniff_slow(s,:) = SH.sniff_slow(s,:) + histcounts(E.absLag(isSniff & isSlow), edges);
    SH.sniff_fast(s,:) = SH.sniff_fast(s,:) + histcounts(E.absLag(isSniff & isFast), edges);
end

function [ph, fwOut, brOut] = coordPhaseFB_events(tEvent, fwEvent, tInsp, tPeak, CYC_MAX)
    ph = nan(size(tEvent));
    fwOut = fwEvent(:);
    brOut = nan(size(tEvent));

    bn = discretize(tEvent, tInsp);

    for i = 1:numel(tEvent)
        b = bn(i);
        if isnan(b) || b < 1 || b > numel(tPeak) || b > numel(tInsp)-1
            continue;
        end

        dur = tInsp(b+1) - tInsp(b);
        if dur <= 0 || dur > CYC_MAX
            continue;
        end

        brOut(i) = 1/dur;

        if tEvent(i) <= tPeak(b)
            den = tPeak(b) - tInsp(b);
            if den > 0
                ph(i) = 0.5 * (tEvent(i) - tInsp(b)) / den;
            end
        else
            den = tInsp(b+1) - tPeak(b);
            if den > 0
                ph(i) = 0.5 + 0.5 * (tEvent(i) - tPeak(b)) / den;
            end
        end
    end

    ok = isfinite(ph) & ph >= 0 & ph <= 1 & isfinite(fwOut) & isfinite(brOut);
    ph = ph(ok);
    fwOut = fwOut(ok);
    brOut = brOut(ok);
end

function OBS = appendPhaseGroups(OBS, ph, fw, br, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ)
    OBS.pop   = [OBS.pop; ph];
    OBS.slowB = [OBS.slowB; ph(br < BASAL_HZ)];
    OBS.fastB = [OBS.fastB; ph(br > SNIFF_HZ)];
    OBS.slowW = [OBS.slowW; ph(fw < SLOW_HZ)];
    OBS.fastW = [OBS.fastW; ph(fw > FAST_HZ)];

    OBS.sBsW = [OBS.sBsW; ph(br < BASAL_HZ & fw < SLOW_HZ)];
    OBS.sBfW = [OBS.sBfW; ph(br < BASAL_HZ & fw > FAST_HZ)];
    OBS.fBsW = [OBS.fBsW; ph(br > SNIFF_HZ & fw < SLOW_HZ)];
    OBS.fBfW = [OBS.fBfW; ph(br > SNIFF_HZ & fw > FAST_HZ)];
end

function SH = addPhaseShuffleCounts(SH, ph, fw, br, edges, s, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ)
    SH.pop(s,:)   = SH.pop(s,:)   + histcounts(ph, edges);
    SH.slowB(s,:) = SH.slowB(s,:) + histcounts(ph(br < BASAL_HZ), edges);
    SH.fastB(s,:) = SH.fastB(s,:) + histcounts(ph(br > SNIFF_HZ), edges);
    SH.slowW(s,:) = SH.slowW(s,:) + histcounts(ph(fw < SLOW_HZ), edges);
    SH.fastW(s,:) = SH.fastW(s,:) + histcounts(ph(fw > FAST_HZ), edges);

    SH.sBsW(s,:) = SH.sBsW(s,:) + histcounts(ph(br < BASAL_HZ & fw < SLOW_HZ), edges);
    SH.sBfW(s,:) = SH.sBfW(s,:) + histcounts(ph(br < BASAL_HZ & fw > FAST_HZ), edges);
    SH.fBsW(s,:) = SH.fBsW(s,:) + histcounts(ph(br > SNIFF_HZ & fw < SLOW_HZ), edges);
    SH.fBfW(s,:) = SH.fBfW(s,:) + histcounts(ph(br > SNIFF_HZ & fw > FAST_HZ), edges);
end

function D = psthGroups(tw, fw, tInsp, cycdur, CYC_MAX, PAIR_WIN, SLOW_HZ, FAST_HZ, BASAL_HZ, SNIFF_HZ)
    D = emptyStruct({'pop','slowB','fastB','slowW','sBsW','fBsW','fastW','sBfW','fBfW'});

    brate = nan(size(tw));
    bn = discretize(tw, tInsp);
    for i = 1:numel(tw)
        b = bn(i);
        if ~isnan(b) && b >= 1 && b <= numel(cycdur)
            dur = cycdur(b);
            if dur > 0 && dur <= CYC_MAX
                brate(i) = 1/dur;
            end
        end
    end

    D.pop = psth(tw, tInsp, PAIR_WIN);

    slowB = brate < BASAL_HZ;
    fastB = brate > SNIFF_HZ;
    slowW = fw < SLOW_HZ;
    fastW = fw > FAST_HZ;

    D.slowB = psth(tw(slowB), tInsp, PAIR_WIN);
    D.fastB = psth(tw(fastB), tInsp, PAIR_WIN);
    D.slowW = psth(tw(slowW), tInsp, PAIR_WIN);
    D.fastW = psth(tw(fastW), tInsp, PAIR_WIN);

    D.sBsW = psth(tw(slowB & slowW), tInsp, PAIR_WIN);
    D.sBfW = psth(tw(slowB & fastW), tInsp, PAIR_WIN);
    D.fBsW = psth(tw(fastB & slowW), tInsp, PAIR_WIN);
    D.fBfW = psth(tw(fastB & fastW), tInsp, PAIR_WIN);
end

function OBS = appendTimeGroups(OBS, D)
    fn = fieldnames(D);
    for i = 1:numel(fn)
        OBS.(fn{i}) = [OBS.(fn{i}); D.(fn{i})];
    end
end

function SH = addTimeShuffleCounts(SH, D, edges, s)
    fn = fieldnames(D);
    for i = 1:numel(fn)
        SH.(fn{i})(s,:) = SH.(fn{i})(s,:) + histcounts(D.(fn{i}), edges);
    end
end

function d = psth(events, refs, win)
    d = [];
    events = events(:);
    refs = refs(:);

    if isempty(events) || isempty(refs)
        return;
    end

    for i = 1:numel(events)
        nearby = refs(abs(refs - events(i)) <= win);
        if isempty(nearby), continue; end
        d = [d; events(i) - nearby(:)]; %#ok<AGROW>
    end

    d = d(abs(d) <= win);
end

function off = randomCircularOffset(Tdur)
    % avoid tiny shifts that act like observed data
    minOff = 0.10 * Tdur;
    off = minOff + rand * (Tdur - 2*minOff);
end

function tsh = circShiftTimes(t, off, t0, tEnd)
    dur = tEnd - t0;
    tsh = mod((t(:) - t0) + off, dur) + t0;
    tsh = sort(tsh);
end

function plotLagPanel(tlo, OBS, SH, g, edges, win, col, label, thresh, yl)
    ax = nexttile(tlo);
    hold(ax,'on'); grid(ax,'on');

    obs      = OBS.(g);
    shCounts = SH.(g);

    % clip edges, data, and shuffle counts to per-panel window
    edgesNow  = edges(edges <= win);
    if isempty(edgesNow) || edgesNow(end) < win
        edgesNow(end+1) = win;
    end
    nBinsNow    = numel(edgesNow) - 1;
    obs         = obs(obs <= win);
    shCountsNow = shCounts(:, 1:nBinsNow);

    if isempty(obs)
        title(ax, sprintf('%s\nobserved n=0', label));
        xlim(ax,[edgesNow(1) edgesNow(end)]);
        axis(ax,'square'); box(ax,'off');
        return;
    end

    plotObsNullProb(ax, obs, shCountsNow, edgesNow, col, label, ...
        '|L-R lag| (ms)', 'probability', yl);

    xline(ax, thresh, 'k--', 'LineWidth', 1.2);

    medObs  = median(obs,'omitnan');
    syncObs = 100*mean(obs <= thresh,'omitnan');

    title(ax, sprintf('%s\nn=%d  med=%.1f ms  sync=%.0f%%', ...
        label, numel(obs), medObs, syncObs));

    xlim(ax,[edgesNow(1) edgesNow(end)]);
    axis(ax,'square'); box(ax,'off');
end

function plotObsNullProb(ax, obs, shCounts, edges, col, label, xlab, ylab, yl)
    hold(ax,'on'); grid(ax,'on');

    obsCounts = histcounts(obs, edges);
    obsProb = obsCounts / max(1,sum(obsCounts));

    shProb = countsToProbRows(shCounts);
    shMed = median(shProb, 1, 'omitnan');
    shLo  = prctile(shProb, 5, 1);
    shHi  = prctile(shProb, 95, 1);

    % observed filled bars
    histogram(ax, obs, edges, ...
        'Normalization','probability', ...
        'FaceColor', col, ...
        'EdgeColor','none');

    % observed stair outline
    stairsX = edges;
    stairsY = [obsProb obsProb(end)];
    stairs(ax, stairsX, stairsY, ...
        'Color', col, ...
        'LineWidth', 2.0);

    % shuffle median and 5-95% envelope as stairs
    stairsMed = [shMed shMed(end)];
    stairsLo  = [shLo  shLo(end)];
    stairsHi  = [shHi  shHi(end)];
    stairs(ax, stairsX, stairsMed, 'k-', 'LineWidth', 1.4);
    stairs(ax, stairsX, stairsLo,  'k:', 'LineWidth', 0.9);
    stairs(ax, stairsX, stairsHi,  'k:', 'LineWidth', 0.9);

    xlabel(ax,xlab);
    ylabel(ax,ylab);

    if ~isempty(yl)
        ylim(ax,yl);
    end

    title(ax, sprintf('%s\nobs n=%d', label, numel(obs)));

    axis(ax,'square'); box(ax,'off');
end

function P = countsToProbRows(C)
    P = nan(size(C));
    for i = 1:size(C,1)
        s = sum(C(i,:));
        if s > 0
            P(i,:) = C(i,:) / s;
        end
    end
end

function plotITIHist(ax, itiL, itiR, edges, colL, colR, xl, yl)
    hold(ax,'on'); grid(ax,'on');

    if ~isempty(itiL)
        histogram(ax, itiL, edges, ...
            'Normalization','probability', ...
            'DisplayStyle','stairs', ...
            'EdgeColor',colL, ...
            'LineWidth',1.6);
    end

    if ~isempty(itiR)
        histogram(ax, itiR, edges, ...
            'Normalization','probability', ...
            'DisplayStyle','stairs', ...
            'EdgeColor',colR, ...
            'LineWidth',1.6);
    end

    if ~isempty(itiL)
        xline(ax, median(itiL,'omitnan'), '-', 'Color',colL,'LineWidth',1.6);
        xline(ax, mean(itiL,'omitnan'), '--', 'Color',colL,'LineWidth',1.0);
    end

    if ~isempty(itiR)
        xline(ax, median(itiR,'omitnan'), '-', 'Color',colR,'LineWidth',1.6);
        xline(ax, mean(itiR,'omitnan'), '--', 'Color',colR,'LineWidth',1.0);
    end

    xlabel(ax,'inter-onset interval (s)');
    ylabel(ax,'probability');
    xlim(ax,xl);
    if ~isempty(yl), ylim(ax,yl); end

    title(ax, sprintf('protraction onset ITI\nL n=%d  R n=%d', numel(itiL), numel(itiR)));
    axis(ax,'square'); box(ax,'off');
end

function [tw, fw] = onsetFreq(tOn, itiMax)
    tOn = tOn(:);
    tOn = sort(tOn);
    n = numel(tOn);
    fw = nan(n,1);

    for k = 2:n-1
        dt = tOn(k+1) - tOn(k-1);
        if dt > 0 && dt <= 2*itiMax
            fw(k) = 2/dt;
        end
    end

    keep = isfinite(fw);
    tw = tOn(keep);
    fw = fw(keep);
end

function f = centeredFreqAtPeaks(t, iti_max)
    n = numel(t);
    f = nan(n,1);

    for k = 2:n-1
        dt = t(k+1) - t(k-1);
        if dt > 0 && dt <= 2*iti_max
            f(k) = 2/dt;
        end
    end
end

function intervals = iti(sig, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)
    if ~any(m)
        intervals = [];
        return;
    end

    sig = sig(:);
    [~,idx] = findpeaks(sig, ...
        'MinPeakProminence', PROT_PROM * std(sig(m)), ...
        'MinPeakDistance', round(MINDIST * fpsW));

    idx = idx(m(idx));

    if numel(idx) < 2
        intervals = [];
        return;
    end

    pk_t = (subsamp(sig, idx) - 1) / fpsW;
    pk_t = pk_t(:);

    ivals = diff(pk_t);
    intervals = ivals(ivals > 0 & ivals < ITI_MAX);
end

function p = subsamp(x, idx)
    x = x(:);
    idx = double(idx(:));
    p = idx;

    in = idx > 1 & idx < numel(x);
    i = idx(in);

    ym = x(i-1);
    y0 = x(i);
    yp = x(i+1);

    den = ym - 2*y0 + yp;
    delta = 0.5*(ym-yp) ./ den;
    delta(~isfinite(delta) | abs(delta) > 0.5) = 0;

    p(in) = i + delta;
end

function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac * prctile(env, 95);
    d = diff([0; a; 0]);
    s = find(d==1);
    e = find(d==-1)-1;

    ep = [t(s) t(e)];

    if isempty(ep)
        return;
    end

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

    bn = arrayfun(@(x) bestnum(x.name), d);
    [~,ix] = max(bn);

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

function exportVectorFig(fig, baseName)
    drawnow;

    set(fig,'Renderer','painters');
    set(fig,'InvertHardcopy','off');
    set(fig,'GraphicsSmoothing','off');

    oldUnits = fig.Units;
    fig.Units = 'inches';
    pos = fig.Position;

    fig.PaperUnits = 'inches';
    fig.PaperSize = [pos(3) pos(4)];
    fig.PaperPosition = [0 0 pos(3) pos(4)];
    fig.PaperPositionMode = 'manual';

    print(fig, baseName + ".pdf", '-dpdf', '-painters');
    print(fig, baseName + ".eps", '-depsc', '-painters');

    fig.Units = oldUnits;

    fprintf('Exported %s.pdf and %s.eps\n', baseName, baseName);
end