% whisk_imbalance_vs_envelope_breath_xcorr_RZ.m
%
% XCORR ONLY:
%   Whisk envelope/amplitude-asymmetry metrics vs breathing, whisk-epoch-gated,
%   with within-epoch shuffle null + inspiration-triggered averages.
%   Three metrics, one 2-panel figure each (vs the same global whisking):
%     imbWave     = env(Lbp-Rbp)              L-R waveform diff env
%     ampAsymRaw  = |Lenv-Renv|               raw L/R amplitude difference (deg)
%     ampAsymNorm = |Lenv-Renv|/(Lenv+Renv)   normalized L/R amplitude asymmetry
%
% Continuous multitaper COHERENCE was split into whisk_breath_coherence_RZ.m.

clear; clc;

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW    = 400;
LIK_THR = 0.6;

% --- filtering ---
WHISK_BP  = [5 60];       % Hz, whisker band
BREATH_BP = [2 20];       % Hz, breathing band
INVERT_BREATH = false;    % loaded breath is ALREADY inhale-up; do NOT re-invert

% --- camera timing ---
EXPOSURE_OFFSET_S = 0.36324e-3;

% --- whisking epoch gating for XCORR ---
USE_EPOCHS = true;
ENV_THR    = 5;           % deg, threshold on smoothed global whisking envelope
ENV_SMOOTH = 0.05;        % s
MIN_DUR    = 2;         % s
MERGE_GAP  = 0.1;        % s

% --- cross-correlation / shuffle stats ---
MAXLAG_S = 0.2;
NSHUF    = 200;
ALPHA    = 0.01;          % p = 0.01

% --- inspiration-triggered average ---
INSP_PROM       = 0.5;
MIN_INSP_DIST_S = 0.05;
TRIG_PRE_S      = 0.2;
TRIG_POST_S     = 0.2;

% --- exclusions ---
EXCLUDE = "5840027";      % whole animal exclusion
EXCLUDE_SESS = {{"5916297",3}, {"5916297",4}, {"5916300",4}};  % specific sessions

SAVE_FIG = false;
% ======================================================================
% Coherence/spectral analysis lives in whisk_breath_coherence_RZ.m (split out).

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

ord = 4;
[b_wb, a_wb] = butter(3,   WHISK_BP/(fpsW/2),  'bandpass');
[b_br, a_br] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

maxLag = round(MAXLAG_S * fpsW);
lagSec = (-maxLag:maxLag)' / fpsW;
nLag   = numel(lagSec);

% --- population accumulators for xcorr ---
allWave  = zeros(nLag,0);
allRaw   = zeros(nLag,0);
allNorm  = zeros(nLag,0);
allWhisk = zeros(nLag,0);

sSurWave  = zeros(nLag,NSHUF);
sSurRaw   = zeros(nLag,NSHUF);
sSurNorm  = zeros(nLag,NSHUF);
sSurWhisk = zeros(nLag,NSHUF);

nSess = 0;
includedSessions = {};

% --- inspiration-triggered accumulators ---
nPre  = round(TRIG_PRE_S*fpsW);
nPost = round(TRIG_POST_S*fpsW);
tRel  = (-nPre:nPost)' / fpsW;

EVwave  = [];
EVraw   = [];
EVnorm  = [];
EVwhisk = [];
EVbr    = [];

% per-session peak summaries
S = struct('animal',{},'run',{}, ...
    'waveR',{},'waveLag',{},'waveSig',{}, ...
    'rawR',{},'rawLag',{},'rawSig',{}, ...
    'normR',{},'normLag',{},'normSig',{}, ...
    'logR',{},'logLag',{},'logSig',{}, ...
    'whiskR',{},'whiskLag',{},'whiskSig',{});

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

fprintf('\n=== SESSION LOOP ===\n');

% ============================ SESSION LOOP ============================
for e = 1:numel(sess)

    animal = sess{e}{1};
    kRun   = sess{e}{2};

    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE))
        fprintf('EXCLUDING animal %s n%d\n', animal, kRun);
        continue;
    end

    if is_excluded_sess(animal, kRun, EXCLUDE_SESS)
        fprintf('EXCLUDING session %s n%d\n', animal, kRun);
        continue;
    end

    fprintf('INCLUDING candidate %s n%d\n', animal, kRun);

    % -------------------- whisker DLC --------------------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, ...
            sprintf('%s_whisk_n%d', animal, kRun)), LIK_THR);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));

    tW = (0:numel(La)-1)' / fpsW;

    La0 = fillmissing(La(:),'linear','EndValues','nearest');
    Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest');

    La0 = La0 - mean(La0,'omitnan');
    Ra0 = Ra0 - mean(Ra0,'omitnan');

    % raw band-passed whisker traces
    LbpFull = filtfilt(b_wb,a_wb, La0);
    RbpFull = filtfilt(b_wb,a_wb, Ra0);

    % Hilbert envelopes for xcorr metrics only
    LenvFull = abs(hilbert(LbpFull));
    RenvFull = abs(hilbert(RbpFull));

    whiskAmpFull    = 0.5 * (LenvFull + RenvFull);
    imbWaveFull     = abs(hilbert(LbpFull - RbpFull));
    ampAsymRawFull  = abs(LenvFull - RenvFull);
    ampAsymNormFull = abs(LenvFull - RenvFull) ./ (LenvFull + RenvFull + eps);
    ampAsymLogFull  = abs(log((LenvFull + eps) ./ (RenvFull + eps)));

    % -------------------- breath --------------------
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, ...
            sprintf('%s_nose_n%d', animal, kRun)), dataRoot);

        if ~isfile(Pn.breath)
            warning('no breath file %s n%d', animal, kRun);
            continue;
        end

        Bs = load(Pn.breath);
        br = Bs.breath(:);
        fb = double(Bs.fps);

        tB = (0:numel(br)-1)' / fb + EXPOSURE_OFFSET_S;
        brwFull = interp1(tB, br, tW, 'linear', NaN);

    catch ME
        warning('breath load failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    % -------------------- crop to shared overlap --------------------
    keep = isfinite(brwFull);

    if nnz(keep) <= 2*maxLag + 1
        warning('breath/whisk overlap too short %s n%d', animal, kRun);
        continue;
    end

    imbWave     = imbWaveFull(keep);
    ampAsymRaw  = ampAsymRawFull(keep);
    ampAsymNorm = ampAsymNormFull(keep);
    ampAsymLog  = ampAsymLogFull(keep);
    whiskAmp    = whiskAmpFull(keep);

    brw = brwFull(keep);

    if INVERT_BREATH
        brw = -brw;
    end

    % raw band-passed breath trace
    Bbp = filtfilt(b_br,a_br, brw);

    N = numel(Bbp);

    % -------------------- whisking mask --------------------
    if USE_EPOCHS
        env = movmean(whiskAmp, max(1, round(ENV_SMOOTH*fpsW)));
        ep  = bool_to_epochs(env > ENV_THR, (0:N-1)'/fpsW, MIN_DUR, MERGE_GAP);

        m = false(N,1);
        for q = 1:size(ep,1)
            s_i = max(1, round(ep(q,1)*fpsW)+1);
            e_i = min(N, round(ep(q,2)*fpsW)+1);
            m(s_i:e_i) = true;
        end
    else
        m = true(N,1);
    end

    nWhiskSec = nnz(m)/fpsW;

    if nnz(m) < 2*maxLag
        warning('too little whisking %s n%d (%.1fs)', animal, kRun, nWhiskSec);
        continue;
    end

    % ============================ XCORR ============================
    [r_wave, sW ] = masked_xcorr(imbWave,     Bbp, m, maxLag, NSHUF);
    [r_raw,  sR ] = masked_xcorr(ampAsymRaw,  Bbp, m, maxLag, NSHUF);
    [r_norm, sN ] = masked_xcorr(ampAsymNorm, Bbp, m, maxLag, NSHUF);
    [r_log,  sL ] = masked_xcorr(ampAsymLog,  Bbp, m, maxLag, NSHUF);
    [r_whisk,sWh] = masked_xcorr(whiskAmp,    Bbp, m, maxLag, NSHUF);

    [wvAbs,wvLag,wvR] = peakAbs(r_wave, lagSec);
    wvSig = wvAbs > prctile(max(abs(sW),[],1),100*(1-ALPHA));

    [rwAbs,rwLag,rwR] = peakAbs(r_raw, lagSec);
    rwSig = rwAbs > prctile(max(abs(sR),[],1),100*(1-ALPHA));

    [nmAbs,nmLag,nmR] = peakAbs(r_norm, lagSec);
    nmSig = nmAbs > prctile(max(abs(sN),[],1),100*(1-ALPHA));

    [lgAbs,lgLag,lgR] = peakAbs(r_log, lagSec);
    lgSig = lgAbs > prctile(max(abs(sL),[],1),100*(1-ALPHA));

    [pwAbs,pwLag,pwR] = peakAbs(r_whisk, lagSec);
    pwSig = pwAbs > prctile(max(abs(sWh),[],1),100*(1-ALPHA));

    allWave(:,end+1)  = r_wave;
    allRaw(:,end+1)   = r_raw;
    allNorm(:,end+1)  = r_norm;
    allWhisk(:,end+1) = r_whisk;

    sSurWave  = sSurWave  + sW;
    sSurRaw   = sSurRaw   + sR;
    sSurNorm  = sSurNorm  + sN;
    sSurWhisk = sSurWhisk + sWh;

    nSess = nSess + 1;
    includedSessions{end+1,1} = sprintf('%s n%d', animal, kRun);

    S(end+1) = struct('animal',animal,'run',kRun, ...
        'waveR',wvR,'waveLag',wvLag,'waveSig',wvSig, ...
        'rawR', rwR,'rawLag', rwLag,'rawSig', rwSig, ...
        'normR',nmR,'normLag',nmLag,'normSig',nmSig, ...
        'logR', lgR,'logLag', lgLag,'logSig', lgSig, ...
        'whiskR',pwR,'whiskLag',pwLag,'whiskSig',pwSig);

    % -------------------- triggered averages --------------------
    [~, insp] = findpeaks(-Bbp, 'MinPeakProminence', INSP_PROM*std(Bbp), ...
                          'MinPeakDistance', round(MIN_INSP_DIST_S*fpsW));

    for ii = insp(:)'
        if ii-nPre < 1 || ii+nPost > N || ~m(ii)
            continue;
        end

        EVwave(end+1,:)  = imbWave(ii-nPre:ii+nPost)';
        EVraw(end+1,:)   = ampAsymRaw(ii-nPre:ii+nPost)';
        EVnorm(end+1,:)  = ampAsymNorm(ii-nPre:ii+nPost)';
        EVwhisk(end+1,:) = whiskAmp(ii-nPre:ii+nPost)';
        EVbr(end+1,:)    = Bbp(ii-nPre:ii+nPost)';
    end

    fprintf(['%s n%d [%4.1fs whisk]:  wave|r|=%.2f@%+.3f(%s)  ' ...
             'raw|r|=%.2f@%+.3f(%s)  norm|r|=%.2f@%+.3f(%s)  ' ...
             'whisk|r|=%.2f@%+.3f(%s)\n'], animal, kRun, nWhiskSec, ...
        wvAbs,wvLag,tern(wvSig,'s','n'), ...
        rwAbs,rwLag,tern(rwSig,'s','n'), ...
        nmAbs,nmLag,tern(nmSig,'s','n'), ...
        pwAbs,pwLag,tern(pwSig,'s','n'));
end

assert(nSess > 0, 'no usable sessions');

fprintf('\n=== INCLUDED USABLE SESSIONS FOR XCORR / MAIN ANALYSIS ===\n');
for i = 1:numel(includedSessions)
    fprintf('  %2d. %s\n', i, includedSessions{i});
end
fprintf('Total usable sessions = %d\n', nSess);

% ============================ XCORR SUMMARIES ============================

muWave = mean(allWave,2);
semWave = std(allWave,0,2) ./ sqrt(nSess);
ciWave = prctile(sSurWave / nSess,[2.5 97.5],2);

muRaw = mean(allRaw,2);
semRaw = std(allRaw,0,2) ./ sqrt(nSess);
ciRaw = prctile(sSurRaw / nSess,[2.5 97.5],2);

muNorm = mean(allNorm,2);
semNorm = std(allNorm,0,2) ./ sqrt(nSess);
ciNorm = prctile(sSurNorm / nSess,[2.5 97.5],2);

muWhisk = mean(allWhisk,2);
semWhisk = std(allWhisk,0,2) ./ sqrt(nSess);
ciWhisk = prctile(sSurWhisk / nSess,[2.5 97.5],2);

waveR  = [S.waveR].';
waveLag = [S.waveLag].';
waveSig = logical([S.waveSig].');

rawR  = [S.rawR].';
rawLag = [S.rawLag].';
rawSig = logical([S.rawSig].');

normR  = [S.normR].';
normLag = [S.normLag].';
normSig = logical([S.normSig].');

logR  = [S.logR].';
logLag = [S.logLag].';
logSig = logical([S.logSig].');

whiskR  = [S.whiskR].';
whiskLag = [S.whiskLag].';
whiskSig = logical([S.whiskSig].');

fprintf('\n=== XCORR SUMMARY: %d sessions ===\n', nSess);
print_summary('global whisking amp',     whiskR, whiskLag, whiskSig, nSess);
print_summary('L-R waveform diff env',    waveR,  waveLag,  waveSig,  nSess);
print_summary('raw L/R amp difference',   rawR,   rawLag,   rawSig,   nSess);
print_summary('norm L/R amp asymmetry',   normR,  normLag,  normSig,  nSess);
print_summary('[diag] log L/R amp asym',  logR,   logLag,   logSig,   nSess);

% ============================ XCORR FIGURES ============================

colImb   = [0.55 0.10 0.65];
colWhisk = [0.00 0.55 0.45];
colBr    = [0.55 0.55 0.55];

C = struct('lagSec',lagSec,'tRel',tRel,'colA',colImb,'colW',colWhisk,'colBr',colBr, ...
    'muW',muWhisk,'semW',semWhisk,'ciW',ciWhisk,'EVW',EVwhisk,'EVbr',EVbr, ...
    'doSave',SAVE_FIG);

% Figure 1: waveform difference envelope
P = C;
P.figName='coupling: L-R waveform diff env';
P.ttl=sprintf('Respiratory coupling: L-R waveform difference envelope vs global whisking   (%d sessions)',nSess);
P.asymName='L-R waveform diff envelope';
P.asymUnitDeg=true;
P.asymYlab='';
P.muA=muWave;
P.semA=semWave;
P.ciA=ciWave;
P.EVA=EVwave;
P.savePath=fullfile(repoRoot,'coupling_waveform_diff_env.png');
coupling_figure(P);

% Figure 2: raw amplitude difference
P = C;
P.figName='coupling: raw L/R amp diff';
P.ttl=sprintf('Respiratory coupling: raw L/R amplitude difference vs global whisking   (%d sessions)',nSess);
P.asymName='raw L/R amplitude difference';
P.asymUnitDeg=true;
P.asymYlab='';
P.muA=muRaw;
P.semA=semRaw;
P.ciA=ciRaw;
P.EVA=EVraw;
P.savePath=fullfile(repoRoot,'coupling_raw_amp_diff.png');
coupling_figure(P);

% Figure 3: normalized amplitude asymmetry
P = C;
P.figName='coupling: norm L/R amp asym';
P.ttl=sprintf('Respiratory coupling: normalized L/R amplitude asymmetry vs global whisking   (%d sessions)',nSess);
P.asymName='normalized L/R amplitude asymmetry';
P.asymUnitDeg=false;
P.asymYlab='normalized L/R amplitude asymmetry (unitless)';
P.muA=muNorm;
P.semA=semNorm;
P.ciA=ciNorm;
P.EVA=EVnorm;
P.savePath=fullfile(repoRoot,'coupling_norm_amp_asym.png');
coupling_figure(P);

% ============================= HELPERS =============================

function [r, sur] = masked_xcorr(x, b, m, maxLag, NSHUF)
    x = x(:);
    b = b(:);
    m = double(m(:));

    nLag = 2*maxLag + 1;

    xm = x .* m;
    xxm = (x.^2) .* m;

    cnt = xcorr(m,   m, maxLag);
    Sx  = xcorr(xm,  m, maxLag);
    Sxx = xcorr(xxm, m, maxLag);

    r = pear_both(b .* m, m, xm, cnt, Sx, Sxx, maxLag, nLag);

    eidx = mask_runs(logical(m));
    sur = zeros(nLag, NSHUF);

    for s = 1:NSHUF
        bs = shift_within_epochs(b, eidx);
        sur(:,s) = pear_both(bs .* m, m, xm, cnt, Sx, Sxx, maxLag, nLag);
    end
end

function rr = pear_both(bm, m, xm, cnt, Sx, Sxx, maxLag, nLag)
    Sb  = xcorr(m,  bm,    maxLag);
    Sbb = xcorr(m,  bm.^2, maxLag);
    Sxb = xcorr(xm, bm,    maxLag);

    num = Sxb - (Sx .* Sb) ./ cnt;
    dX  = Sxx - (Sx.^2)    ./ cnt;
    dB  = Sbb - (Sb.^2)    ./ cnt;

    rr  = zeros(nLag,1);
    ok  = cnt > 1 & dX > 0 & dB > 0;

    rr(ok) = num(ok) ./ sqrt(dX(ok) .* dB(ok));
end

function bs = shift_within_epochs(b, eidx)
    bs = b;

    for k = 1:size(eidx,1)
        i0 = eidx(k,1);
        i1 = eidx(k,2);
        L = i1 - i0 + 1;

        if L < 3
            continue;
        end

        bs(i0:i1) = circshift(b(i0:i1), randi(L-1));
    end
end

function eidx = mask_runs(m)
    m = logical(m(:));
    d = diff([false; m; false]);
    eidx = [find(d == 1), find(d == -1) - 1];
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:));
    d = diff([false; a; false]);

    s = find(d == 1);
    e = find(d == -1) - 1;

    ep = [t(s) t(e)];

    if isempty(ep)
        return;
    end

    ep2 = ep(1,:);

    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= mergeGap
            ep2(end,2) = ep(i,2);
        else
            ep2(end+1,:) = ep(i,:);
        end
    end

    ep = ep2;
    ep = ep(ep(:,2) - ep(:,1) >= minDur, :);
end

function [pk, lagAt, sgn] = peakAbs(r, lagSec)
    [pk, idx] = max(abs(r));
    lagAt = lagSec(idx);
    sgn   = r(idx);
end

function band(ax, x, lo, hi, col)
    fill(ax, [x; flipud(x)], [lo; flipud(hi)], col, ...
        'FaceAlpha',0.15, 'EdgeColor','none');
end

function h = plot_trig(ax, t, X, col, lw)
    if isempty(X)
        h = plot(ax, NaN, NaN, 'Color', col);
        return;
    end

    mu = mean(X,1,'omitnan');
    n  = sum(isfinite(X),1);
    se = std(X,0,1,'omitnan') ./ sqrt(max(n,1));

    fill(ax, [t; flipud(t)], [(mu+se)'; flipud((mu-se)')], col, ...
        'FaceAlpha',0.18, 'EdgeColor','none');

    h = plot(ax, t, mu, 'Color', col, 'LineWidth', lw);
end

function print_summary(name, sgnR, lag, sig, nSess)
    fprintf('%-28s : signed r med=%+.2f   |r| med=%.2f   lag med=%+.3fs   sig %d/%d\n', ...
        name, median(sgnR,'omitnan'), median(abs(sgnR),'omitnan'), ...
        median(lag,'omitnan'), nnz(sig), nSess);
end

function coupling_figure(P)

    f = figure('Color','w','Position',[60 120 1240 520], 'Name', P.figName);
    tlo = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    title(tlo, P.ttl, 'FontWeight','bold','Interpreter','none');

    % Panel A
    axA = nexttile(tlo);
    hold(axA,'on');
    grid(axA,'on');

    band(axA, P.lagSec, P.muA-P.semA, P.muA+P.semA, P.colA);
    band(axA, P.lagSec, P.muW-P.semW, P.muW+P.semW, P.colW);

    plot(axA, P.lagSec, P.ciA, '--', 'Color',P.colA, 'LineWidth',0.8);
    plot(axA, P.lagSec, P.ciW, '--', 'Color',P.colW, 'LineWidth',0.8);

    hA = plot(axA, P.lagSec, P.muA, 'Color',P.colA, 'LineWidth',2);
    hW = plot(axA, P.lagSec, P.muW, 'Color',P.colW, 'LineWidth',2);

    xline(axA,0,'k:');
    yline(axA,0,'k:');

    xlabel(axA,'lag (s)   (+ = whisk signal follows breath)');
    ylabel(axA,'cross-correlation (coeff)');
    title(axA,'A  population mean xcorr  (shaded = \pmSEM; dashed = shuffle 95% CI)');

    legend([hA hW], {P.asymName,'global whisking amplitude'}, ...
        'Location','northeast','Box','off');

    axis(axA,'square');
    box(axA,'off');

    % Panel B
    axB = nexttile(tlo);
    hold(axB,'on');
    grid(axB,'on');

    if P.asymUnitDeg
        yyaxis(axB,'left');
        ha = plot_trig(axB, P.tRel, P.EVA, P.colA, 2.0);
        hw = plot_trig(axB, P.tRel, P.EVW, P.colW, 2.0);
        ylabel(axB,'whisk Hilbert amplitude (deg)');
        set(axB,'YColor','k');

        yyaxis(axB,'right');
        hb = plot_trig(axB, P.tRel, P.EVbr, P.colBr, 1.2);
        ylabel(axB,'breath, filtered (a.u.)');
        set(axB,'YColor',P.colBr);

        hs = [ha hw hb];
        leg = {P.asymName,'global whisking amplitude','breath'};

    else
        yyaxis(axB,'left');
        ha = plot_trig(axB, P.tRel, P.EVA, P.colA, 2.0);
        ylabel(axB, P.asymYlab);
        set(axB,'YColor','k');

        yyaxis(axB,'right');
        hw = plot_trig(axB, P.tRel, P.EVW, P.colW, 2.0);
        ylabel(axB,'global whisking amplitude (deg)');
        set(axB,'YColor',P.colW);

        hs = [ha hw];
        leg = {P.asymName,'global whisking amplitude'};
    end

    xline(axB,0,'k--');
    xlabel(axB,'time from inspiration onset (s)');
    title(axB, sprintf('B  inspiration-triggered average  (%d events)', size(P.EVA,1)));

    legend(hs, leg, 'Location','northeast','Box','off');

    xlim(axB,[P.tRel(1) P.tRel(end)]);
    axis(axB,'square');
    box(axB,'off');

    if P.doSave
        exportgraphics(f, P.savePath, 'Resolution',150);
        fprintf('saved %s\n', P.savePath);
    end
end

function p = signrank_safe(a, b)
    a = a(:);
    b = b(:);
    ok = isfinite(a) & isfinite(b);

    if nnz(ok) < 2
        p = NaN;
        return;
    end

    try
        p = signrank(a(ok), b(ok));
    catch
        p = NaN;
    end
end

function s = tern(c, a, b)
    if c
        s = a;
    else
        s = b;
    end
end

function tf = is_excluded_sess(animal, kRun, excl)
    tf = false;

    for q = 1:numel(excl)
        if strcmp(char(excl{q}{1}), char(animal)) && double(excl{q}{2}) == kRun
            tf = true;
            return;
        end
    end
end

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv'));

    S = {};
    key = {};

    for i = 1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');

        if isempty(tok)
            continue;
        end

        k = sprintf('%s_%s', tok{1}, tok{2});

        if any(strcmp(key,k))
            continue;
        end

        key{end+1} = k;
        S{end+1} = {tok{1}, str2double(tok{2})};
    end

    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S);
        rn = cellfun(@(c) c{2}, S);
        [~,o] = sortrows([an(:) rn(:)]);
        S = S(o);
    end
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