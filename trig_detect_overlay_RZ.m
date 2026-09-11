% trig_detect_overlay_RZ  (script)
% Single raster of the DETECTED events for a RANDOM example session — one plot.
% Rows (ticks at event times): inspiration onset / inspiration peak / R peak /
% L protraction onset / R protraction onset.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
BR_BP     = [2 15];     % breath bandpass (Hz)
PROT_PROM = 1;          % protraction trough prominence (x std of BP whisk)
THR_FRAC  = 0.30; MIN_DUR = 1.0; MERGE_GAP = 0.20; ENV_SMOOTH = 0.20;   % whisk-epoch gating
GATE_EPOCH= true;       % only show events inside whisk epochs (matches the analysis)
ZOOM      = [];         % s, window ([] = full)
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

ad = dir(char(dataRoot)); sess = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess),'no sessions');
pick = sess{randi(numel(sess))}; animal = pick{1}; kRun = pick{2};

% whisk events
M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);    % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
t  = (0:numel(La)-1)'/fpsW;
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);                       % whisk epochs
m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end
phiR = angle(hilbert(xR));
tRpk = t(find(phiR(1:end-1)<0 & phiR(2:end)>=0));                              % R peak
[~,iLp]=findpeaks(-xL,'MinPeakProminence',PROT_PROM*std(xL),'MinPeakDistance',round(0.02*fpsW)); tLon=(iLp-1)/fpsW;
[~,iRp]=findpeaks(-xR,'MinPeakProminence',PROT_PROM*std(xR),'MinPeakDistance',round(0.02*fpsW)); tRon=(iRp-1)/fpsW;

% breath events
tIon=[]; tIpk=[];
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kRun)), dataRoot);
if isfile(Pn.breath)
    Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps);
    [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
    phib=angle(hilbert(brf)); tIpk=(find(phib(1:end-1)<0 & phib(2:end)>=0)-1)/fb;     % inspiration peak
    [~,il]=findpeaks(-brf,'MinPeakProminence',0.5*std(brf),'MinPeakDistance',round(fb/BR_BP(2))); tIon=(il-1)/fb;
end

if GATE_EPOCH
    g = @(e) e(m(min(max(round(e*fpsW)+1,1),numel(m))));   % keep events inside whisk epochs
    tIon=g(tIon); tIpk=g(tIpk); tRpk=g(tRpk); tLon=g(tLon); tRon=g(tRon);
end
zoomwin = ZOOM; if isempty(zoomwin), zoomwin=[t(1) t(end)]; end
ev   = {tIon, tIpk, tRpk, tLon, tRon};
names= {'insp onset','insp peak','R peak','L protract onset','R protract onset'};
cols = {[0.85 0.1 0.1],[0.85 0.1 0.1],[0 0.4 0.85],[0 0.5 0],[0 0.4 0.85]};
fprintf('random session: %s n%d\n', animal, kRun);

figure('Color','w','Position',[60 250 1280 360]); hold on; grid on;
ylr = [0.5 numel(ev)+0.5];
for q=1:size(ep,1)                                           % shade whisk epochs
    if ep(q,2)<zoomwin(1) || ep(q,1)>zoomwin(2), continue; end
    patch(ep(q,[1 2 2 1]), ylr([1 1 2 2]), [0.3 0.75 0.3], 'FaceAlpha',0.12, 'EdgeColor','none');
end
for r = 1:numel(ev)
    e = ev{r}; e = e(e>=zoomwin(1) & e<=zoomwin(2));
    plot(e, r*ones(size(e)), '|', 'Color', cols{r}, 'MarkerSize', 12, 'LineWidth', 1);
end
xlim(zoomwin); ylim(ylr);
yticks(1:numel(ev)); yticklabels(names); set(gca,'YDir','reverse');
xlabel('time (s)');
gs=''; if GATE_EPOCH, gs=' [gated to epochs]'; end
title(sprintf('%s n%d  —  detected events%s  (green = whisk epoch)', animal, kRun, gs));

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
