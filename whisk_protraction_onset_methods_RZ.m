% whisk_protraction_onset_methods_RZ  (script)
% THREE random WIN-second whisking windows (drawn from detected whisking epochs),
% each showing protraction-ONSET detection by TWO methods on the same BP trace:
%   (1) findpeaks  -> troughs of BP angle      (v, green)   = start of forward sweep
%   (2) Hilbert phase reset -> wrap +pi->-pi   (o, magenta) = trough of analytic signal
% For a BP whisk cycle the analytic phase is 0 at the protraction PEAK and wraps
% +pi->-pi at the trough, so the phase reset marks the same protraction onset.
% The wrapped Hilbert phase is drawn faintly on the right axis so the resets show.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
animal    = "5916296";   % session animal id
kRun      = 1;           % run index (n<k>)
SIDE      = 'L';         % which whisker to show: 'L' or 'R'
WIN       = 3;           % window length (s)
NPAN      = 3;           % number of random windows (panels)
BP        = [6 30];      % fast-whisking bandpass (Hz)
PROT_PROM = 0.8;         % findpeaks trough prominence (x std of BP whisk)
MINDIST   = 0.03;        % s, min distance between onsets (findpeaks)
THR_FRAC  = 0.5; MIN_DUR = 1.0; MERGE_GAP = 0.10; ENV_SMOOTH = 0.10;  % whisk-epoch gating
RNG_SEED  = [];          % [] = new windows each run; integer to reproduce
fpsW      = 400;         % whisk camera fps
% ======================================================================

if ~isempty(RNG_SEED), rng(RNG_SEED); end
addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(4, BP/(fpsW/2), 'bandpass');

% ---- load whisking, protraction-positive BP angles ----
M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', char(animal), kRun)), 0.6);  % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));   % LEFT  (x mirrored)
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));    % RIGHT
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
if upper(SIDE)=='L', x = xL; else, x = xR; end                       % chosen side
t  = (0:numel(x)-1)'/fpsW;

% ---- whisk epochs (gate detections to active whisking) ----
env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
ep  = ep(ep(:,2)-ep(:,1) >= WIN, :);
assert(~isempty(ep), 'no whisking epoch >= %g s in %s n%d', WIN, char(animal), kRun);
m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end

% ---- METHOD 1: findpeaks troughs (protraction onsets) ----
[~,iTr] = findpeaks(-x, 'MinPeakProminence', PROT_PROM*std(x), 'MinPeakDistance', round(MINDIST*fpsW));
iTr = iTr(m(iTr));

% ---- METHOD 2: Hilbert phase reset (+pi -> -pi wrap = trough) ----
phi  = angle(hilbert(x));                       % wrapped phase, -pi..pi
iWr  = find(diff(phi) < -pi) + 1;               % index just after the +pi->-pi wrap
iWr  = iWr(m(iWr));

% ---- pick NPAN random windows (length-weighted) ----
w = (ep(:,2)-ep(:,1)) - WIN;
t0 = zeros(NPAN,1);
for k = 1:NPAN
    ei = find(cumsum(w) >= rand*sum(w), 1, 'first');
    t0(k) = ep(ei,1) + rand*w(ei);
end

fprintf('%s n%d (%c): findpeaks onsets %d, phase-reset onsets %d (whole trace, in-epoch)\n', ...
    char(animal), kRun, upper(SIDE), numel(iTr), numel(iWr));

% ---- plot ----
cX = [0 0.45 0.75]; cTr = [0 0.55 0]; cWr = [0.85 0.1 0.7]; cP = [0.6 0.6 0.6];
figure('Color','w','Position',[100 80 1050 760]);
tl = tiledlayout(NPAN,1,'TileSpacing','compact','Padding','compact');
for k = 1:NPAN
    a0 = t0(k); a1 = a0 + WIN;
    inT = iTr(t(iTr)>=a0 & t(iTr)<=a1);
    inW = iWr(t(iWr)>=a0 & t(iWr)<=a1);
    nexttile;
    yyaxis right;                                   % faint Hilbert phase
    pp = plot(t, phi, '-', 'Color',cP, 'LineWidth',0.4); pp.Color(4)=0.5;
    ylim([-pi pi]); set(gca,'YColor',cP,'YTick',[-pi 0 pi],'YTickLabel',{'-\pi','0','\pi'}); ylabel('Hilbert phase');
    yyaxis left;  hold on; grid on;
    hX  = plot(t, x, '-', 'Color',cX, 'LineWidth',1.0);
    hTr = plot(t(inT), x(inT), 'v', 'Color',cTr, 'MarkerFaceColor',cTr, 'MarkerSize',7);
    hWr = plot(t(inW), x(inW), 'o', 'Color',cWr, 'MarkerSize',9, 'LineWidth',1.4);
    set(gca,'YColor','k'); xlim([a0 a1]); ylabel('whisk angle (deg, BP)');
    title(sprintf('window %d:  %.2f-%.2f s', k, a0, a1));
    if k==NPAN, xlabel('time (s)'); end
    if k==1, legend([hX hTr hWr], {sprintf('%c BP whisk',upper(SIDE)),'findpeaks onset','phase-reset onset'}, ...
            'Orientation','horizontal','Location','northoutside'); end
end
title(tl, sprintf('%s n%d (%c) — protraction onset: findpeaks vs Hilbert phase reset', char(animal), kRun, upper(SIDE)), ...
      'Interpreter','none','FontWeight','bold');

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
