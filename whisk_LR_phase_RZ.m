% whisk_LR_phase_RZ  (script)
% Left-right whisker synchrony in CYCLE-NORMALIZED phase (not absolute ms).
% Use one whisker's protraction onsets as the cycle clock: consecutive onsets
% bracket a cycle, phase 0 -> 1 by LINEAR interpolation between them. Place each
% onset of the OTHER whisker at its fractional phase within the cycle it falls in.
% Phase is recentred to [-0.5 0.5]: 0 = synchronous, +/-0.5 = anti-phase.
% Both directions pooled (L-in-R cycles + R-in-L cycles). Whisk epochs only.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [6 30];     % fast-whisking bandpass (Hz)
MINDIST   = 0.04;       % s, min distance between onsets (findpeaks)
PROT_PROM = 1;          % protraction trough prominence (x std of BP whisk)
CYC_MAX   = 0.20;       % s, drop reference cycles longer than this (cross-bout gaps)
SYNC_FRAC = 0.10;       % "synchronous" = |phase| < SYNC_FRAC of a cycle
NBINS     = 40;         % phase histogram bins over [-0.5 0.5]
THR_FRAC  = 0.30; MIN_DUR = 1.0; MERGE_GAP = 0.20; ENV_SMOOTH = 0.2;   % whisk-epoch gating
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[b,a] = butter(3, BP/(fpsW/2), 'bandpass');

phL = []; phR = []; nsess = 0;   % phL: L onset phase in R cycle; phR: R onset phase in L cycle
ad = dir(char(dataRoot));
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);  % lik<0.6 -> linear interp
        swL = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        swR = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        xL = filtfilt(b,a, fillmissing(swL-mean(swL,'omitnan'),'linear'));
        xR = filtfilt(b,a, fillmissing(swR-mean(swR,'omitnan'),'linear'));
        t   = (0:numel(xL)-1)'/fpsW;
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end
        [~,iL] = findpeaks(-xL, 'MinPeakProminence',PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
        [~,iR] = findpeaks(-xR, 'MinPeakProminence',PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
        iL = iL(m(iL)); iR = iR(m(iR));
        tL = (iL-1)/fpsW; tR = (iR-1)/fpsW;
        phL = [phL; cyclephase(tL, tR, CYC_MAX)]; %#ok<AGROW>  L onsets within R cycles
        phR = [phR; cyclephase(tR, tL, CYC_MAX)]; %#ok<AGROW>  R onsets within L cycles
        nsess = nsess + 1;
    end
end
ph = [phL; phR];
assert(~isempty(ph),'no phased onsets');
fracSync = mean(abs(ph) < SYNC_FRAC);
fprintf('%d phased onsets (%d L-in-R, %d R-in-L), %d sessions\n', numel(ph), numel(phL), numel(phR), nsess);
fprintf('median |phase| = %.3f cycles | fraction synchronous (|phase|<%.2f) = %.2f\n', median(abs(ph)), SYNC_FRAC, fracSync);

edges = linspace(-0.5, 0.5, NBINS+1);
figure('Color','w','Position',[200 200 660 460]); hold on; grid on;
histogram(ph, edges, 'Normalization','probability','FaceColor',[0.4 0.4 0.4],'EdgeColor','none');
xline(0,'k--','LineWidth',1);
xline( SYNC_FRAC,'r--','LineWidth',1.2); xline(-SYNC_FRAC,'r--','LineWidth',1.2);   % +/- sync window
xlim([-0.5 0.5]); xlabel('other-whisker onset phase within cycle  (0 = sync, \pm0.5 = anti-phase)'); ylabel('prob');
title(sprintf('L-R cycle-normalized onset phase  (n=%d, %d sess)  sync %.0f%% (|phase|<%.2f)', ...
    numel(ph), nsess, 100*fracSync, SYNC_FRAC));

% ================= helpers =================
function p = cyclephase(tEv, tRef, cycMax)
% phase in [-0.5 0.5] of each event time tEv within the reference cycle (tRef(k),tRef(k+1)]
    p = [];
    if numel(tRef) < 2 || isempty(tEv), return; end
    bin = discretize(tEv, tRef);            % which ref cycle each event falls in
    ok  = ~isnan(bin);
    bin = bin(ok); te = tEv(ok);
    D   = tRef(bin+1) - tRef(bin);          % cycle duration
    f   = (te - tRef(bin)) ./ D;            % fractional phase [0,1)
    f   = f(D <= cycMax);                    % drop cross-bout cycles
    p   = mod(f + 0.5, 1) - 0.5;            % recentre: 0 = synchronous
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
