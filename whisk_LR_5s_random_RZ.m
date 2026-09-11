% whisk_LR_5s_random_RZ  (script)
% Plot a random 5 s window of LEFT/RIGHT whisking drawn from inside a detected
% whisking EPOCH (so the window is guaranteed to contain active whisking).
% Whisk detection matches whisk_LR_sync_RZ.m / whisk_breath_raster_RZ.m:
%   angle = base->tip (deg), LEFT x mirrored so protraction is +up on both sides
%   BP fast-whisking signal -> epochs from |Hilbert| envelope -> protraction
%   peaks = maxima of BP angle (most-protracted point of each sweep).

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
animal    = "5916296";   % session animal id
kRun      = 1;           % run index (n<k>)
WIN       = 3;           % window length (s)
BP        = [6 30];      % fast-whisking bandpass (Hz)
PROT_PROM = 0.8;         % protraction trough prominence (x std of BP whisk)
MINDIST   = 0.02;        % s, min distance between whisk onsets (findpeaks)
THR_FRAC  = 0.5; MIN_DUR = 1.0; MERGE_GAP = 0.10; ENV_SMOOTH = 0.10;  % whisk-epoch gating
SHOW_PEAKS= true;        % overlay protraction-peak markers
RNG_SEED  = [];          % [] = different window each run; set an integer to reproduce
fpsW      = 400;         % whisk camera fps
% ======================================================================

if ~isempty(RNG_SEED), rng(RNG_SEED); end
addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % (pick_csv lives below)
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

% ---- load whisking, protraction-positive BP angles ----
M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', char(animal), kRun)), 0.6);  % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));   % LEFT  (x mirrored)
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));    % RIGHT
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
t  = (0:numel(xL)-1)'/fpsW;

% ---- whisk epochs (envelope > THR_FRAC*95pct of |Hilbert| of BP L/R, averaged) ----
env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
ep  = ep(ep(:,2)-ep(:,1) >= WIN, :);                 % keep epochs long enough for the window
assert(~isempty(ep), 'no whisking epoch >= %g s in %s n%d (loosen THR_FRAC/MIN_DUR)', WIN, char(animal), kRun);

% ---- pick a random epoch (weighted by length) then a random 5 s window inside it ----
w  = (ep(:,2)-ep(:,1)) - WIN;                          % slack available in each epoch
ei = find(cumsum(w) >= rand*sum(w), 1, 'first');      % longer epochs more likely
t0 = ep(ei,1) + rand*w(ei);                            % window start
t1 = t0 + WIN;

% ---- protraction peaks in window (maxima of BP angle = most protracted) ----
[~,iL] = findpeaks(xL, 'MinPeakProminence', PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
[~,iR] = findpeaks(xR, 'MinPeakProminence', PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
iL = iL(t(iL)>=t0 & t(iL)<=t1);  iR = iR(t(iR)>=t0 & t(iR)<=t1);

fprintf('%s n%d: %d whisk epochs >= %gs | window %.2f-%.2f s (epoch %d/%d) | prot peaks L %d R %d\n', ...
    char(animal), kRun, size(ep,1), WIN, t0, t1, ei, size(ep,1), numel(iL), numel(iR));

% ---- plot ----
cL = [0 0.5 0]; cR = [0 0.4 0.85];
figure('Color','w','Position',[120 200 1100 420]); hold on; grid on;
hL = plot(t, xL, '-', 'Color',cL, 'LineWidth',1.0);
hR = plot(t, xR, '-', 'Color',cR, 'LineWidth',1.0);
if SHOW_PEAKS
    plot(t(iL), xL(iL), '^', 'Color',cL, 'MarkerFaceColor',cL, 'MarkerSize',6);
    plot(t(iR), xR(iR), '^', 'Color',cR, 'MarkerFaceColor',cR, 'MarkerSize',6);
end
xlim([t0 t1]); box on;
xlabel('time (s)'); ylabel('whisk angle (deg, BP)');
legend([hL hR], {'L','R'}, 'Location','northeast');
title(sprintf('%s n%d  —  random %g s from a whisking epoch  (%.2f-%.2f s)', char(animal), kRun, WIN, t0, t1), ...
      'Interpreter','none');

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
