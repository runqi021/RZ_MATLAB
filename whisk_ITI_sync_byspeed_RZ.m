% whisk_ITI_sync_byspeed_RZ  (script)
% Like whisk_ITI_sync_RZ, but SPLITS the L-R protraction-peak synchrony by whisk
% SPEED: each L protraction peak is classified SLOW (centered local whisk freq <
% SLOW_HZ) or FAST (> FAST_HZ) from its bracketing inter-peak intervals (same
% centered-frequency convention as whisk_breath_coord_panels_RZ / onsetFreq),
% paired to the nearest R peak within +/-MAXLAG, and the signed lag (L - R) is
% pooled per speed class. One figure, 1x2: SLOW lag histogram | FAST lag histogram.
% (Whisk epochs only, pooled over sessions; animal EXCLUDE dropped.)

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [5 30];     % fast-whisking bandpass (Hz)
MINDIST   = 0.025;      % s, min distance between peaks/onsets (findpeaks)
PROT_PROM = 1;          % peak/onset prominence (x std of BP whisk)
MAXLAG    = 0.08;       % s, max |L-R| lag to count as a pair
SYNC_FRAC = 0.1;        % "synchronous" = |lag| < SYNC_FRAC * mean whisk cycle (per class)
ITI_MAX   = 0.3;        % s, max bracketing interval for a valid centered cycle
SLOW_HZ   = 10;         % centered whisk freq < this = SLOW whisk
FAST_HZ   = 20;         % centered whisk freq > this = FAST whisk
PHASE_MAX = 1;          % clip synchrony-phase axis to [0 PHASE_MAX]
NBINS     = 30;         % phase histogram bins over [0 PHASE_MAX]
THR_FRAC  = 0; MIN_DUR = 0.5; MERGE_GAP = 0.20; ENV_SMOOTH = 0;   % whisk-epoch gating
EXCLUDE   = "5840027";  % animal ids to drop from the analysis ("0027")
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

S = struct('slow', struct('lag',[],'cyc',[]), 'fast', struct('lag',[],'cyc',[]));
nsess = 0;
ad = dir(char(dataRoot));
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    if any(strcmp(ad(ai).name, EXCLUDE)), continue; end          % skip excluded animals
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);  % lik<0.6 -> linear interp
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        t  = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end

        % ---- protraction PEAKS (local maxima), in whisk epochs ----
        [~,iL] = findpeaks(xL, 'MinPeakProminence',PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
        [~,iR] = findpeaks(xR, 'MinPeakProminence',PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
        iL = iL(m(iL)); iR = iR(m(iR));
        nsess = nsess + 1;
        if numel(iL) < 3 || isempty(iR), continue; end
        tL = (subsamp(xL,iL)-1)/fpsW; tR = (subsamp(xR,iR)-1)/fpsW;   % sub-frame peak times

        % ---- per-L-peak: centered cycle freq + nearest-R lag, binned by speed ----
        dL = diff(tL);
        for k = 2:numel(tL)-1
            dprev = dL(k-1); dnext = dL(k);                  % bracketing intervals
            if dprev > ITI_MAX || dnext > ITI_MAX, continue; end
            per = (dprev+dnext)/2; f = 1/per;                % centered cycle period -> Hz
            [d,j] = min(abs(tR - tL(k)));
            if d > MAXLAG, continue; end
            lag = 1000*(tL(k) - tR(j));                      % signed L - R lag (ms)
            if     f < SLOW_HZ, S.slow.lag(end+1)=lag; S.slow.cyc(end+1)=per; %#ok<AGROW>
            elseif f > FAST_HZ, S.fast.lag(end+1)=lag; S.fast.cyc(end+1)=per; %#ok<AGROW>
            end
        end
    end
end
fprintf('%d sessions  |  paired peaks: slow %d, fast %d\n', ...
    nsess, numel(S.slow.lag), numel(S.fast.lag));
assert(~isempty(S.slow.lag) && ~isempty(S.fast.lag), 'need paired peaks in both speed classes');

% ============================ FIGURE (overlay, PHASE) ============================
% phase = |L-R lag| / (mean whisk cycle of that speed class); absolute, clipped [0 PHASE_MAX].
classes = {'slow','fast'};
labs    = {sprintf('SLOW whisk <%g Hz',SLOW_HZ), sprintf('FAST whisk >%g Hz',FAST_HZ)};
cols    = {[0.2 0.4 0.8],[0.85 0.1 0.1]};
pedges  = linspace(0, PHASE_MAX, NBINS+1);           % shared phase bins

figure('Color','w','Position',[200 180 620 560]); ax = axes; hold(ax,'on'); grid(ax,'on');
hLeg = gobjects(1,2); txt = cell(1,2);
for c = 1:2
    lag = S.(classes{c}).lag(:); cyc = S.(classes{c}).cyc(:);
    meanCycMs = 1000*mean(cyc);
    ph = abs(lag/1000) / mean(cyc);                  % |lag| / average duty-cycle time
    medPh = median(ph); iqrPh = iqr(ph); fracSync = mean(ph < SYNC_FRAC);

    hLeg(c) = histogram(ax, ph, pedges, 'Normalization','probability', ...
        'FaceColor',cols{c}, 'FaceAlpha',0.5, 'EdgeColor',cols{c}, 'LineWidth',0.5);
    txt{c} = sprintf('%s  (n=%d, %.1f Hz): median %.3f, IQR %.3f, sync %.0f%%', ...
        labs{c}, numel(ph), 1000/meanCycMs, medPh, iqrPh, 100*fracSync);
    fprintf('%s: n=%d | cycle %.0f ms (%.1f Hz) | median phase %.3f | IQR %.3f | sync %.0f%%\n', ...
        labs{c}, numel(ph), meanCycMs, 1000/meanCycMs, medPh, iqrPh, 100*fracSync);
end
xline(ax, SYNC_FRAC,'k--','LineWidth',1);
xlim(ax, [0 PHASE_MAX]);
xlabel(ax,'|L-R lag| / whisk cycle  (phase)'); ylabel(ax,'fraction of pairs'); axis(ax,'square');
legend(ax, hLeg, txt, 'Location','northoutside','Box','off');
title(ax, sprintf('L-R protraction peak synchrony phase by whisk speed  (%d sessions)', nsess));

% ================= helpers =================
function p = subsamp(x, idx)
    % parabolic sub-sample refinement of peak locations (removes frame quantization)
    x = x(:); idx = double(idx(:)); p = idx;
    in = idx>1 & idx<numel(x); i = idx(in);
    ym = x(i-1); y0 = x(i); yp = x(i+1);
    den = ym - 2*y0 + yp;                 % <0 at a true maximum
    delta = 0.5*(ym - yp) ./ den;
    delta(~isfinite(delta) | abs(delta)>0.5) = 0;   % guard flat / ill-conditioned
    p(in) = i + delta;
end
function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac*prctile(env,95);
    d = diff([0; a; 0]); s = find(d==1); e = find(d==-1)-1;
    ep = [t(s) t(e)];
    if ~isempty(ep)
        m = ep(1,:);
        for i=2:size(ep,1)
            if ep(i,1)-m(end,2) <= mergeGap, m(end,2)=ep(i,2); else, m(end+1,:)=ep(i,:); end %#ok<AGROW>
        end
        ep = m(m(:,2)-m(:,1) >= minDur, :);
    end
end
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
