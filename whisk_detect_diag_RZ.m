function whisk_detect_diag_RZ()
% whisk_detect_diag_RZ  Inspect the whisking detection for ONE session:
%   (1) WHISK-EPOCH detection: envelope (|Hilbert| of BP L/R, averaged) vs an
%       ABSOLUTE ENV_THR (deg) threshold, with detected epochs shaded (same
%       logic as whisk_imbalance_vs_envelope_breath_xcorr_RZ.m: BP [5 40],
%       ENV_SMOOTH 0.05 s, ENV_THR 5 deg, MIN_DUR 1 s, MERGE_GAP 0).
%   (2) PROTRACTION-ONSET detection: troughs of the BP whisker angle (start of
%       the forward sweep, "large amplitude" = prominence gate), as in the raster.
% Single panel (zoom): BP whisk L/R + (smoothed) envelope +/- + threshold + epoch
% shading; protraction-onset markers optional (SHOW_PEAKS).

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
animal    = "5916297";%"5840027"; %"5916296";
kRun      = 4;
BP        = [5 40];     % whisk bandpass (Hz)  -- matches xcorr scripts
ENV_THR   = 5;          % deg, ABSOLUTE hard envelope threshold -- matches xcorr scripts
ENV_SMOOTH= 0.05;       % s, envelope smoothing (movmean) before thresholding; 0 = none
MIN_DUR   = 0.5;          % s, min epoch duration
MERGE_GAP = 0.00;       % s, merge epochs closer than this
PROT_PROM = 1;        % protraction trough prominence (x std of BP whisk)
SHOW_PEAKS= false;       % overlay protraction-onset markers (true/false)
ZOOM      = [];    % s, zoom window [t0 t1]; [] = full trace
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);    % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
t  = (0:numel(La)-1)'/fpsW;
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
env_raw = (abs(hilbert(xL)) + abs(hilbert(xR)))/2;
env = movmean(env_raw, max(1,round(ENV_SMOOTH*fpsW)));   % smoothed envelope (used for detection)
thr = ENV_THR;                                  % absolute deg threshold
ep  = detect(env, t, ENV_THR, MIN_DUR, MERGE_GAP);
% protraction onsets (troughs of BP angle), L and R
[~,iL] = findpeaks(-xL,'MinPeakProminence',PROT_PROM*std(xL),'MinPeakDistance',round(0.02*fpsW));
[~,iR] = findpeaks(-xR,'MinPeakProminence',PROT_PROM*std(xR),'MinPeakDistance',round(0.02*fpsW));

fprintf('%s n%d: dur=%.0fs | %d whisk epochs (%.0f%% of time) | protraction onsets L %d R %d\n', ...
    animal, kRun, t(end), size(ep,1), 100*sum(ep(:,2)-ep(:,1))/t(end), numel(iL), numel(iR));

figure('Color','w','Position',[60 200 1240 430]);
ax2 = axes; hold(ax2,'on'); grid(ax2,'on');
yl2 = [min([xL;xR;-env]) max([xL;xR;env])];
for q=1:size(ep,1)
    patch(ax2, ep(q,[1 2 2 1]), yl2([1 1 2 2]), [0.3 0.75 0.3], 'FaceAlpha',0.12, 'EdgeColor','none');
end
hL=plot(ax2, t, xL, '-', 'Color',[0 0.5 0]);
hR=plot(ax2, t, xR, '-', 'Color',[0 0.4 0.85]);
he=plot(ax2, t, env, 'k-', 'LineWidth',1.4); plot(ax2, t, -env, 'k-', 'LineWidth',1.4);   % +/- smoothed envelope
yline(ax2, thr, 'r--'); yline(ax2, -thr, 'r--');                                          % detection threshold
if SHOW_PEAKS
    plot(ax2, t(iL), xL(iL), 'v', 'Color',[0 0.5 0],   'MarkerFaceColor',[0 0.5 0],   'MarkerSize',6);
    plot(ax2, t(iR), xR(iR), 'v', 'Color',[0 0.4 0.85],'MarkerFaceColor',[0 0.4 0.85],'MarkerSize',6);
end
zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end   % [] = full trace
xlim(ax2, zoomwin); ylim(ax2, yl2);
xlabel(ax2,'time (s)'); ylabel(ax2,'whisk angle (deg, BP)');
legend(ax2, [hL hR he], {'L (BP)','R (BP)','envelope (\pm)'}, 'Location','northeastoutside');
ttl = sprintf('%s n%d  %g-%g s  (green = whisk epoch; red-- = thr)', animal,kRun,zoomwin(1),zoomwin(2));
if SHOW_PEAKS, ttl = [ttl '  v=protraction onset']; end
title(ax2, ttl);
end

% ================= helpers =================
function ep = detect(env, t, thr, minDur, mergeGap)
    a = env(:) > thr;
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
