% whisk_LR_sync_RZ  (script)
% Left-right whisker protraction-PEAK synchrony, pooled over all sessions.
% Core: detect protraction peaks (local maxima of BP angle) on L and R, pair
% each L peak to the nearest R peak within +/-maxLag, collect the signed lag
% (L - R), and report median / IQR / fraction synchronous.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [6 30];     % fast-whisking bandpass (Hz)
MINDIST   = 0.02;       % s, min distance between peaks/onsets (findpeaks)
PROT_PROM = 1;        % peak/onset prominence (x std of BP whisk)
MAXLAG    = 0.07;       % s, max |L-R| lag to count as a pair
SYNC_FRAC = 0.1;        % "synchronous" = |lag| < SYNC_FRAC * mean whisk cycle (ITI)
ITI_MAX   = 0.25;       % s, drop ITIs longer than this (cross-bout / gap straddles)
NBINS     = 40;         % histogram bins (ITI over [0 ITI_MAX]; lag over +/-MAXLAG)
THR_FRAC  = 0; MIN_DUR = 0.5; MERGE_GAP = 0.20; ENV_SMOOTH = 0;   % whisk-epoch gating
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[b,a] = butter(3, BP/(fpsW/2), 'bandpass');
maxLag = round(MAXLAG*fpsW);

lag_ms = []; iti_s = []; nsess = 0;
ad = dir(char(dataRoot));
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);  % lik<0.6 -> linear interp
        % protraction-positive angles (larger = more protracted)
        swL = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        swR = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        xL = filtfilt(b,a, fillmissing(swL-mean(swL,'omitnan'),'linear'));
        xR = filtfilt(b,a, fillmissing(swR-mean(swR,'omitnan'),'linear'));
        % whisk epochs (envelope > THR_FRAC*95pct)
        t   = (0:numel(xL)-1)'/fpsW;
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end
        % protraction peaks = local maxima (larger angle = more protracted), in epochs only
        [~,iL] = findpeaks(xL, 'MinPeakProminence',PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
        [~,iR] = findpeaks(xR, 'MinPeakProminence',PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
        iL = iL(m(iL)); iR = iR(m(iR));
        if isempty(iL) || isempty(iR), continue; end
        dtL=diff(iL)/fpsW; dtR=diff(iR)/fpsW; iti_s=[iti_s; dtL(dtL<0.3); dtR(dtR<0.3)]; %#ok<AGROW> whisk cycle
        % pair each L peak to nearest R peak within maxLag
        for k = 1:numel(iL)
            [d,j] = min(abs(iR - iL(k)));
            if d <= maxLag
                lag_ms(end+1) = 1000*(iL(k)-iR(j))/fpsW; %#ok<AGROW>
            end
        end
        nsess = nsess + 1;
    end
end
assert(~isempty(lag_ms),'no paired onsets');

medianLag = median(lag_ms);
iqrLag    = iqr(lag_ms);
meanITIms = 1000*mean(iti_s);
sync_ms   = SYNC_FRAC*meanITIms;           % synchrony window = +/- 10% of mean cycle
fracSync  = mean(abs(lag_ms) < sync_ms);
fprintf('%d paired peaks, %d sessions\n', numel(lag_ms), nsess);
fprintf('mean whisk cycle = %.1f ms  -> sync window = +/-%.1f ms (%.0f%% of cycle)\n', meanITIms, sync_ms, 100*SYNC_FRAC);
fprintf('median L-R lag = %.1f ms\n', medianLag);
fprintf('IQR = %.1f ms\n', iqrLag);
fprintf('fraction synchronous (|lag|<%.1f ms) = %.2f\n', sync_ms, fracSync);

figure('Color','w','Position',[200 200 640 460]); hold on; grid on;
bwid  = 2000*MAXLAG/NBINS;                          % bin width (ms) across +/-MAXLAG
nEdge = ceil(max(abs(lag_ms))/bwid) + 1;
edges = ((-nEdge:nEdge) + 0.5) * bwid;             % edges at +/-bwid/2, +/-3bwid/2 ... -> a bin CENTERED on 0
histogram(lag_ms, edges, 'Normalization','probability', 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
xline(0,'k--','LineWidth',1);
xline( sync_ms,'r--','LineWidth',1.2); xline(-sync_ms,'r--','LineWidth',1.2);   % +/- sync window
xlabel('L peak - R peak lag (ms)'); ylabel('% of pairs'); ytickformat('percentage');
axis square;
title(sprintf('Left-right protraction peak synchrony  (n=%d, %d sess)  median %.1f ms, IQR %.1f, sync %.0f%% (|lag|<%.1f ms)', ...
    numel(lag_ms), nsess, medianLag, iqrLag, 100*fracSync, sync_ms));

% ================= helpers =================
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
