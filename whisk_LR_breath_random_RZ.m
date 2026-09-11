% whisk_LR_breath_random_RZ  (script)
% Random WIN-second window of LEFT/RIGHT whisking drawn from inside a detected
% whisking EPOCH, overlaid with the simultaneous BREATHING trace.
%   left  y-axis : L/R fast-whisk BP angle (deg) + protraction peaks (^)
%   right y-axis : band-passed breathing + inspiration peaks (o)
% Events use the peak convention (whisk_LR_sync_RZ / whisk_breath_raster_RZ):
%   protraction peak = maximum of BP whisker angle (most-protracted point)
%   inspiration peak = maximum of the (inhale-up) breath signal (peak inhalation)
% Set FIX_T0 to pin the window start (e.g. 57) instead of drawing randomly.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
animal    = "5916297";   % session animal id
kRun      = 1;           % run index (n<k>)
WIN       = 3;           % window length (s)
FIX_T0    = [];          % [] = random window inside an epoch; set a number (s) to pin the start
BP        = [6 30];      % fast-whisking bandpass (Hz)
BR_BP     = [2 15];      % breath bandpass before peak detection (Hz)
PROT_PROM = 0.8;         % protraction peak prominence (x std of BP whisk)
BR_PROM   = 0.5;         % inspiration peak prominence (x std of BP breath)
MINDIST   = 0.02;        % s, min distance between whisk peaks (findpeaks)
THR_FRAC  = 0.5; MIN_DUR = 1.0; MERGE_GAP = 0.10; ENV_SMOOTH = 0.10;  % whisk-epoch gating
RNG_SEED  = [];          % [] = different window each run; set an integer to reproduce
fpsW      = 400;         % whisk camera fps
% ======================================================================

if ~isempty(RNG_SEED), rng(RNG_SEED); end
addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % thermal_resolve_paths
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
ep  = ep(ep(:,2)-ep(:,1) >= WIN, :);                  % keep epochs long enough for the window
assert(~isempty(ep), 'no whisking epoch >= %g s in %s n%d (loosen THR_FRAC/MIN_DUR)', WIN, char(animal), kRun);

% ---- choose window ----
if isempty(FIX_T0)
    w  = (ep(:,2)-ep(:,1)) - WIN;                      % slack available in each epoch
    ei = find(cumsum(w) >= rand*sum(w), 1, 'first');   % longer epochs more likely
    t0 = ep(ei,1) + rand*w(ei);
else
    t0 = FIX_T0; ei = NaN;                             % user-pinned start
end
t1 = t0 + WIN;

% ---- protraction peaks in window (maxima of BP angle) ----
[~,iL] = findpeaks(xL, 'MinPeakProminence', PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
[~,iR] = findpeaks(xR, 'MinPeakProminence', PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
iL = iL(t(iL)>=t0 & t(iL)<=t1);  iR = iR(t(iR)>=t0 & t(iR)<=t1);

% ---- load breathing + inspiration peaks (maxima of inhale-up breath) ----
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', char(animal), kRun)), dataRoot);
assert(isfile(Pn.breath), 'no _breath.mat for %s n%d (run drawROI_N_lpsub first)', char(animal), kRun);
Bs  = load(Pn.breath); br = Bs.breath(:); fpsB = double(Bs.fps);
tB  = (0:numel(br)-1)'/fpsB;
[bb,ab] = butter(2, BR_BP/(fpsB/2), 'bandpass');
brf = filtfilt(bb,ab, fillmissing(br,'linear'));
[~,iB] = findpeaks(brf, 'MinPeakProminence', BR_PROM*std(brf), 'MinPeakDistance', round(fpsB/BR_BP(2)));
iB = iB(tB(iB)>=t0 & tB(iB)<=t1);

fprintf('%s n%d: window %.2f-%.2f s | prot peaks L %d R %d | insp peaks %d\n', ...
    char(animal), kRun, t0, t1, numel(iL), numel(iR), numel(iB));

% ---- plot (whisk left axis, breath right axis) ----
cL = [0 0.5 0]; cR = [0 0.4 0.85]; cB = [0.10 0.10 0.10];
figure('Color','w','Position',[120 200 1100 440]);
yyaxis left;  hold on;
hL = plot(t, xL, '-', 'Color',cL, 'LineWidth',1.0);
hR = plot(t, xR, '-', 'Color',cR, 'LineWidth',1.0);
ylabel('whisk angle (deg, BP)'); set(gca,'YColor','k');
yyaxis right; hold on;
hB = plot(tB, brf, '-', 'Color',cB, 'LineWidth',1.3);
ylabel('breathing (BP, inhale up)'); set(gca,'YColor',cB);
xlim([t0 t1]); grid on; box on;
xlabel('time (s)');
legend([hL hR hB], {'L whisk','R whisk','breath'}, 'Orientation','horizontal','Location','northoutside');
title(sprintf('%s n%d  —  %g s whisking + breathing  (%.2f-%.2f s)', char(animal), kRun, WIN, t0, t1), ...
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
