% whisk_ITI_epoch_RZ  (script)
% Distribution of whisk-cycle ITI (protraction onset -> next protraction onset),
% L and R pooled, restricted to WHISKING EPOCHS, pooled over all sessions.
% Onset = trough of BP angle (findpeaks(-x), prominence + min-distance gated).

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
PROT_PROM = 1;          % protraction trough prominence (x std of BP whisk)
ITI_MAX   = 0.20;       % s, drop ITIs longer than this (cross-bout gaps)
EDGES     = 0:0.005:0.20;  % ITI histogram bins (s)
THR_FRAC  = 0.30; MIN_DUR = 1.0; MERGE_GAP = 0.20; ENV_SMOOTH = 0.2;   % whisk-epoch gating
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

iti = []; nsess = 0;
ad = dir(char(dataRoot));
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);    % lik<0.6 -> linear interp
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        t  = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end
        for s = {xL, xR}
            x = s{1};
            [~,loc] = findpeaks(-x,'MinPeakProminence',PROT_PROM*std(x),'MinPeakDistance',round(0.02*fpsW));
            loc = loc(m(loc));                    % keep onsets inside whisk epochs
            d = diff((loc-1)/fpsW); d = d(d>0 & d<=ITI_MAX);
            iti = [iti; d(:)]; %#ok<AGROW>
        end
        nsess = nsess + 1;
    end
end
assert(~isempty(iti),'no whisk cycles in epochs');
fprintf('%d whisk cycles, %d sessions  |  ITI mean %.0f ms (%.1f Hz), median %.0f ms\n', ...
    numel(iti), nsess, 1000*mean(iti), 1/mean(iti), 1000*median(iti));

figure('Color','w','Position',[200 200 660 440]); hold on; grid on;
histogram(iti, EDGES, 'Normalization','probability','FaceColor',[0.4 0.4 0.4],'EdgeColor','none');
xline(mean(iti),  'b--','LineWidth',1.2);
xline(median(iti),'r--','LineWidth',1.2);
legend({'ITI', sprintf('mean %.0f ms (%.1f Hz)',1000*mean(iti),1/mean(iti)), ...
    sprintf('median %.0f ms (%.1f Hz)',1000*median(iti),1/median(iti))},'Location','northeast');
xlim([EDGES(1) EDGES(end)]); xlabel('whisk-cycle ITI (s)'); ylabel('prob');
title(sprintf('whisk-cycle ITI in whisking epochs (L+R)  —  %d cycles, %d sess  |  mean %.0f ms / median %.0f ms', ...
    numel(iti), nsess, 1000*mean(iti), 1000*median(iti)));

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
