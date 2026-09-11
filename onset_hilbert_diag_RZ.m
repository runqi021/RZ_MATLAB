% onset_hilbert_diag_RZ  (script)
% Random example session, two panels:
%   TOP : BP whisk L (green) with detected L protraction onsets (troughs) marked.
%   BOT : BP whisk R (blue)  with detected R protraction onsets (troughs) marked.
% Just the detected onsets on the signal — nothing else.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
PROT_PROM = 1;          % protraction trough prominence (x std of BP whisk)
ZOOM      = [20 28];    % s, window ([] = full)
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
assert(~isempty(sess),'no sessions'); pick=sess{randi(numel(sess))}; animal=pick{1}; kRun=pick{2};

M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);    % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
t  = (0:numel(La)-1)'/fpsW;
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
[~,iL]=findpeaks(-xL,'MinPeakProminence',PROT_PROM*std(xL),'MinPeakDistance',round(0.02*fpsW));   % L onsets (trough)
[~,iR]=findpeaks(-xR,'MinPeakProminence',PROT_PROM*std(xR),'MinPeakDistance',round(0.02*fpsW));   % R onsets (trough)

zoomwin = ZOOM; if isempty(zoomwin), zoomwin=[t(1) t(end)]; end
fprintf('random session: %s n%d  |  L onsets %d  R onsets %d\n', animal, kRun, numel(iL), numel(iR));
colL=[0 0.5 0]; colR=[0 0.4 0.85];

figure('Color','w','Position',[60 200 1280 520]);
ax1=subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t, xL, '-', 'Color',colL);
plot(ax1, t(iL), xL(iL), 'v', 'Color',colL, 'MarkerFaceColor',colL, 'MarkerSize',6);
xlim(ax1,zoomwin); ylabel(ax1,'L whisk angle (deg, BP)');
title(ax1, sprintf('%s n%d  —  L protraction onsets (v)', animal,kRun));

ax2=subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t, xR, '-', 'Color',colR);
plot(ax2, t(iR), xR(iR), 'v', 'Color',colR, 'MarkerFaceColor',colR, 'MarkerSize',6);
xlim(ax2,zoomwin); ylabel(ax2,'R whisk angle (deg, BP)'); xlabel(ax2,'time (s)');
title(ax2, sprintf('%s n%d  —  R protraction onsets (v)', animal,kRun));
linkaxes([ax1 ax2],'x');

% ================= helpers =================
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
