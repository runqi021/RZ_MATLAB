% whisk_ITI_sync_setpoint07_RZ  (script)
% Same as whisk_ITI_sync_RZ, but only KEEPS sessions whose per-session L/R
% slow set-point correlation (low-pass <LP_HZ, demeaned) exceeds R_THR.
%   ITI   : inter-onset intervals for protraction (trough) & retraction (peak)
%           onsets, L vs R, plus the pooled population.
%   SYNC  : L-R protraction-PEAK synchrony; signed lag (L - R), median/IQR/% sync.
% 1x2 figure: protraction onset ITI (L+R) | L-R protraction peak lag.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [5 30];     % fast-whisking bandpass (Hz)
LP_HZ     = 1;          % low-pass cutoff for the slow set-point (Hz)  <-- gate signal
R_THR     = 0.8;        % L/R set-point corr threshold
KEEP_ABOVE= false;      % true: keep r>R_THR ; false: keep the SKIPPED (r<=R_THR) sessions
MINDIST   = 0.025;      % s, min distance between peaks/onsets (findpeaks)
PROT_PROM = 1;          % peak/onset prominence (x std of BP whisk)
MAXLAG    = 0.08;       % s, max |L-R| lag to count as a pair
SYNC_FRAC = 0.1;        % "synchronous" = |lag| < SYNC_FRAC * mean whisk cycle (ITI)
ITI_MAX   = 0.3;        % s, drop ITIs longer than this (cross-bout / gap straddles)
NBINS     = 40;         % histogram bins (ITI over [0 ITI_MAX]; lag over +/-MAXLAG)
THR_FRAC  = 0; MIN_DUR = 0.5; MERGE_GAP = 0.20; ENV_SMOOTH = 0;   % whisk-epoch gating
EXCLUDE   = "5840027";  % animal ids to drop from the analysis ("0027")
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
[bl,al] = butter(3, LP_HZ/(fpsW/2),'low');
maxLag = round(MAXLAG*fpsW);

lag_ms = []; iti_s = []; nsess = 0; nkept = 0;
ITI = struct('protL',[],'protR',[]);
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
        nsess = nsess + 1;

        % ---- GATE: per-session slow set-point L/R correlation ----
        sL = filtfilt(bl,al, fillmissing(La,'linear')); sL = sL - mean(sL,'omitnan');
        sR = filtfilt(bl,al, fillmissing(Ra,'linear')); sR = sR - mean(sR,'omitnan');
        gset = isfinite(sL)&isfinite(sR);
        rset = corr(sL(gset), sR(gset));
        pass = (rset > R_THR) == logical(KEEP_ABOVE);   % keep above OR below per toggle
        fprintf('%s n%d: set-point r = %.2f  ->  %s\n', ad(ai).name, kk, rset, ...
            ternary(pass,'KEEP','skip'));
        if ~pass, continue; end
        nkept = nkept + 1;

        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end

        % ---- ITI: protraction (trough) onsets, L & R ----
        ITI.protL = [ITI.protL; iti(-xL, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)]; %#ok<AGROW>
        ITI.protR = [ITI.protR; iti(-xR, m, t, PROT_PROM, fpsW, MINDIST, ITI_MAX)]; %#ok<AGROW>

        % ---- SYNC: protraction PEAKS (local maxima), pair L->nearest R ----
        [~,iL] = findpeaks(xL, 'MinPeakProminence',PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
        [~,iR] = findpeaks(xR, 'MinPeakProminence',PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
        iL = iL(m(iL)); iR = iR(m(iR));
        if isempty(iL) || isempty(iR), continue; end
        dtL=diff(iL)/fpsW; dtR=diff(iR)/fpsW; iti_s=[iti_s; dtL(dtL<=ITI_MAX); dtR(dtR<=ITI_MAX)]; %#ok<AGROW> whisk cycle (same cutoff as ITI plot)
        for k = 1:numel(iL)
            [d,j] = min(abs(iR - iL(k)));
            if d <= maxLag, lag_ms(end+1) = 1000*(iL(k)-iR(j))/fpsW; end %#ok<AGROW>
        end
    end
end
gateStr = ternary(KEEP_ABOVE, sprintf('r>%.2f',R_THR), sprintf('r<=%.2f',R_THR));
fprintf('%d/%d sessions kept (set-point %s)  |  ITI counts: protL %d protR %d\n', ...
    nkept, nsess, gateStr, numel(ITI.protL), numel(ITI.protR));
assert(~isempty(lag_ms),'no paired peaks among kept sessions');
medianLag = median(lag_ms); iqrLag = iqr(lag_ms);
meanITIms = 1000*mean(iti_s); sync_ms = SYNC_FRAC*meanITIms; fracSync = mean(abs(lag_ms) < sync_ms);
fprintf('%d paired peaks  |  mean cycle %.1f ms -> sync window +/-%.1f ms\n', numel(lag_ms), meanITIms, sync_ms);
fprintf('median L-R lag = %.1f ms | IQR = %.1f ms | frac sync = %.2f\n', medianLag, iqrLag, fracSync);

% ============================ FIGURE (1x2) ============================
edges = linspace(0, ITI_MAX, NBINS+1);
colL=[0.85 0.1 0.1]; colR=[0 0.4 0.85];
figure('Color','w','Position',[120 160 1100 460]);

% ---- PROTRACTION onset ITI: L & R, with COMBINED (L+R) mean & median ----
prot = [ITI.protL; ITI.protR]; mP = mean(prot); medP = median(prot);
ax1=subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
hh(ax1, ITI.protL, edges, colL); hh(ax1, ITI.protR, edges, colR);
xline(ax1, mP,'b--','LineWidth',1.4); xline(ax1, medP,'k--','LineWidth',1.4);
xlim(ax1,[0 ITI_MAX]); xlabel(ax1,'inter-onset interval (s)'); ylabel(ax1,'prob');
legend(ax1, {'L','R', sprintf('mean %.0f ms (%.1f Hz)',1000*mP,1/mP), ...
    sprintf('median %.0f ms (%.1f Hz)',1000*medP,1/medP)}, 'Location','northeast');
title(ax1, sprintf('PROTRACTION onset ITI (L+R combined, n=%d)   mean %.0f ms / median %.0f ms', ...
    numel(prot), 1000*mP, 1000*medP));
axis(ax1,'square');

ax4=subplot(1,2,2); hold(ax4,'on'); grid(ax4,'on');
bwid = 2000*MAXLAG/NBINS; nEdge = ceil(max(abs(lag_ms))/bwid) + 1;
ledges = ((-nEdge:nEdge) + 0.5) * bwid;            % a bin CENTERED on 0
histogram(ax4, lag_ms, ledges, 'Normalization','probability','FaceColor',[0.4 0.4 0.4],'EdgeColor','none');
xline(ax4,0,'k--','LineWidth',1);
xline(ax4, sync_ms,'r--','LineWidth',1.2); xline(ax4,-sync_ms,'r--','LineWidth',1.2);
xlabel(ax4,'L peak - R peak lag (ms)'); ylabel(ax4,'% of pairs'); ytickformat(ax4,'percentage');
axis(ax4,'square');
title(ax4, sprintf('L-R protraction peak synchrony  (n=%d pairs, %d/%d sess, set-pt %s)  median %.1f ms, sync %.0f%%', ...
    numel(lag_ms), nkept, nsess, gateStr, medianLag, 100*fracSync));

% ================= helpers =================
function d = iti(sig, m, t, prom, fpsW, mindist, itimax)
    [~,loc] = findpeaks(sig, 'MinPeakProminence', prom*std(sig), 'MinPeakDistance', round(mindist*fpsW));
    loc = loc(m(loc));                       % keep onsets inside whisk epochs
    d = diff((loc-1)/fpsW);                  % inter-onset intervals (s)
    d = d(d>0 & d<=itimax); d = d(:);        % drop cross-bout / gap straddles
end
function hh(ax, x, edges, col)
    histogram(ax, x, edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor',col,'LineWidth',1.4);
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
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
