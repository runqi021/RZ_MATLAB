% whisk_ITI_dist_RZ  (script)
% Distribution of inter-onset intervals (ITI) for whisker protraction and
% retraction onsets, Left vs Right, pooled over all sessions (whisk epochs only).
%   protraction onset = trough of BP angle  (findpeaks(-x), start of forward sweep)
%   retraction  onset = peak  of BP angle  (findpeaks(+x), start of backward sweep)
% Two panels: (1) protraction ITI L/R, (2) retraction ITI L/R.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
PROT_PROM = 1;          % trough/peak prominence (x std of BP whisk)
ITI_MAX   = 0.20;       % s, drop ITIs longer than this (cross-bout / gap straddles)
NBINS     = 30;         % histogram bins over [0 ITI_MAX]
THR_FRAC  = 0.50; MIN_DUR = 0.5; MERGE_GAP = 0.20; ENV_SMOOTH = 0.2;   % whisk-epoch gating
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

ad = dir(char(dataRoot)); nsess=0;
ITI = struct('protL',[],'protR',[],'retL',[],'retR',[]);
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
        ITI.protL = [ITI.protL; iti(-xL, m, t, PROT_PROM, fpsW, ITI_MAX)];
        ITI.protR = [ITI.protR; iti(-xR, m, t, PROT_PROM, fpsW, ITI_MAX)];
        ITI.retL  = [ITI.retL;  iti( xL, m, t, PROT_PROM, fpsW, ITI_MAX)];
        ITI.retR  = [ITI.retR;  iti( xR, m, t, PROT_PROM, fpsW, ITI_MAX)];
        nsess = nsess+1;
    end
end
fprintf('%d sessions  |  ITI counts: protL %d protR %d retL %d retR %d\n', ...
    nsess, numel(ITI.protL), numel(ITI.protR), numel(ITI.retL), numel(ITI.retR));

edges = linspace(0, ITI_MAX, NBINS+1);
colL=[0.85 0.1 0.1]; colR=[0 0.4 0.85];
figure('Color','w','Position',[60 120 1500 440]);

ax1=subplot(1,3,1); hold(ax1,'on'); grid(ax1,'on');
hh(ax1, ITI.protL, edges, colL); hh(ax1, ITI.protR, edges, colR);
xlim(ax1,[0 ITI_MAX]); xlabel(ax1,'inter-onset interval (s)'); ylabel(ax1,'prob');
title(ax1, sprintf('PROTRACTION onset ITI   L %.0f ms (%.1f Hz) / R %.0f ms (%.1f Hz)', ...
    1000*median(ITI.protL), 1/median(ITI.protL), 1000*median(ITI.protR), 1/median(ITI.protR)));
legend(ax1,{'L','R'},'Location','northeast');

ax2=subplot(1,3,2); hold(ax2,'on'); grid(ax2,'on');
hh(ax2, ITI.retL, edges, colL); hh(ax2, ITI.retR, edges, colR);
xlim(ax2,[0 ITI_MAX]); xlabel(ax2,'inter-onset interval (s)'); ylabel(ax2,'prob');
title(ax2, sprintf('RETRACTION onset ITI   L %.0f ms (%.1f Hz) / R %.0f ms (%.1f Hz)', ...
    1000*median(ITI.retL), 1/median(ITI.retL), 1000*median(ITI.retR), 1/median(ITI.retR)));
legend(ax2,{'L','R'},'Location','northeast');

% ---- POPULATION: all onsets pooled (prot + ret, L + R) ----
pop = [ITI.protL; ITI.protR; ITI.retL; ITI.retR];
mPop = mean(pop); medPop = median(pop);
ax3=subplot(1,3,3); hold(ax3,'on'); grid(ax3,'on');
histogram(ax3, pop, edges, 'Normalization','probability','FaceColor',[0.4 0.4 0.4],'EdgeColor','none');
xline(ax3, mPop,  'b--','LineWidth',1.4);
xline(ax3, medPop,'r--','LineWidth',1.4);
xlim(ax3,[0 ITI_MAX]); xlabel(ax3,'inter-onset interval (s)'); ylabel(ax3,'prob');
legend(ax3, {'population', sprintf('mean %.0f ms (%.1f Hz)',1000*mPop,1/mPop), ...
    sprintf('median %.0f ms (%.1f Hz)',1000*medPop,1/medPop)}, 'Location','northeast');
title(ax3, sprintf('POPULATION ITI (all onsets, n=%d)   mean %.0f ms / median %.0f ms', ...
    numel(pop), 1000*mPop, 1000*medPop));
linkaxes([ax1 ax2 ax3],'xy');

% ================= helpers =================
function d = iti(sig, m, t, prom, fpsW, itimax)
    [~,loc] = findpeaks(sig, 'MinPeakProminence', prom*std(sig), 'MinPeakDistance', round(0.02*fpsW));
    loc = loc(m(loc));                       % keep onsets inside whisk epochs
    d = diff((loc-1)/fpsW);                  % inter-onset intervals (s)
    d = d(d>0 & d<=itimax);                  % drop cross-bout / gap straddles
    d = d(:);
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
