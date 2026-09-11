% raster_L_trig_R_RZ  (script)
% For every RIGHT whisk cycle (consecutive R protraction onsets), show where the
% LEFT protraction onsets fall, NORMALIZED by that cycle's duration:
%   phase 0 = this R onset, phase 1 = next R onset.
% Top: raster (one row per R cycle). Bottom: phase histogram of L onsets.
% Pooled over all sessions, restricted to whisking epochs.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
PROT_PROM = 1;          % protraction trough prominence (x std of BP whisk)
CYC_MAX   = 0.20;       % s, drop R cycles longer than this (gap straddles)
NBINS     = 20;         % phase histogram bins over [0 1]
SORT      = 'phase';    % row order: 'phase' = by earliest L phase | 'none' = occurrence
THR_FRAC  = 0.50; MIN_DUR = 0.5; MERGE_GAP = 0.20; ENV_SMOOTH = 0.2;   % whisk-epoch gating
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

ad = dir(char(dataRoot)); PH={}; KEY=[]; DUR=[]; nsess=0;
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk)))), continue; end
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(a).name,kk)), 0.6);    % lik<0.6 -> linear interp
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        t  = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end
        [~,iL]=findpeaks(-xL,'MinPeakProminence',PROT_PROM*std(xL),'MinPeakDistance',round(0.02*fpsW));
        [~,iR]=findpeaks(-xR,'MinPeakProminence',PROT_PROM*std(xR),'MinPeakDistance',round(0.02*fpsW));
        iL = iL(m(iL)); iR = iR(m(iR));                 % keep onsets inside whisk epochs
        tL = (iL-1)/fpsW; tR = (iR-1)/fpsW;
        for r = 1:numel(tR)-1
            D = tR(r+1) - tR(r);                         % R whisk-cycle duration
            if D<=0 || D>CYC_MAX, continue; end
            ph = (tL - tR(r))/D;                         % L onset phase within this R cycle
            ph = ph(ph>=0 & ph<1);
            if isempty(ph), continue; end
            PH{end+1} = ph(:)'; KEY(end+1) = min(ph); DUR(end+1) = D;  %#ok<AGROW>  earliest L phase = sort key
        end
        nsess = nsess+1;
    end
end
assert(~isempty(PH),'no R cycles');
N = numel(PH);
switch lower(SORT)
    case 'phase', [~,ord] = sort(KEY,'ascend');  % by earliest L phase
    otherwise,    ord = 1:N;                      % occurrence order
end
X=[]; Y=[];
for r = 1:N
    p = PH{ord(r)}; X=[X p]; Y=[Y r*ones(1,numel(p))]; %#ok<AGROW>
end
durMean = 1000*mean(DUR); durMed = 1000*median(DUR); durSD = 1000*std(DUR);
fprintf('%d R cycles, %d sessions  |  cycle dur: mean %.1f ms (%.1f Hz), median %.1f ms, SD %.1f ms\n', ...
    N, nsess, durMean, 1000/durMean, durMed, durSD);

figure('Color','w','Position',[80 50 600 840]);
% ---- raster ----
ax1=subplot(3,1,[1 2]); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, X, Y, '.', 'Color',[0.85 0.1 0.1], 'MarkerSize',3);   % L onsets (red)
xline(ax1,0,'b-','LineWidth',1); xline(ax1,1,'b-','LineWidth',1);  % R onsets (cycle bounds)
xlim(ax1,[0 1]); ylim(ax1,[1 N]); set(ax1,'YDir','reverse');
ylabel(ax1,'R cycle #');
title(ax1,sprintf('L onsets within R whisk cycle (0,1 = R onsets)  —  %d cycles, %d sess  |  mean cycle %.0f ms (%.1f Hz)', ...
    N, nsess, durMean, 1000/durMean));
% ---- histogram ----
ax2=subplot(3,1,3); hold(ax2,'on'); grid(ax2,'on');
edges = linspace(0, 1, NBINS+1);
histogram(ax2, X, edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0.85 0.1 0.1],'LineWidth',1.2);
xline(ax2,0,'b:'); xline(ax2,1,'b:'); xlim(ax2,[0 1]);
xlabel(ax2,'L onset phase within R cycle (0,1 = R onset)'); ylabel(ax2,'prob');

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
