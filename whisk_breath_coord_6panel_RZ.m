% whisk_breath_coord_6panel_RZ  (script)
% Breathing<->whisking coordination + left-right whisk synchronization for ONE
% session, as an event-timing analysis (NOT continuous correlation/linear fit).
%
% One 2x3 figure:
%   BILATERAL WHISKING
%     (1) L/R fast-whisk traces + detected protraction onsets   (sanity plot)
%     (2) L-R protraction onset-lag histogram                   (main sync result)
%     (3) cycle-by-cycle L-R lag over time                      (stable vs episodic)
%   BREATH <-> WHISK
%     (4) inspiration-triggered protraction raster              (main coordination plot)
%     (5) protraction latency histogram (first prot. after inspiration onset)
%     (6) whisks per breath histogram
% Convention: positive lag / right of t=0 = AFTER the trigger event.
%
% Events:
%   protraction onset = trough of band-passed whisker angle (start of forward sweep)
%   inspiration onset = trough of the (inhale-up) breath signal (start of inhalation)
% Whisk angle is base->tip, LEFT x mirrored so protraction is +up on both sides
% (same convention as breath_whisk_overlay_RZ / whisk_breath_raster_RZ).
% Secondary circular-phase stats (PLV, mean phase) are printed to the console
% only -- the recommended figure stays at 6 event-timing panels.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
animal    = "5916296";   % session animal id
kRun      = 1;           % run index (n<k>)
BP        = [5 30];      % fast-whisking bandpass (Hz)
BR_BP     = [2 15];      % breath bandpass before inspiration-onset detection (Hz)
PROT_PROM = 1.0;         % protraction trough prominence (x std of BP whisk) = "large amplitude"
BR_PROM   = 0.5;         % inspiration trough prominence (x std of BP breath)
MINDIST   = 0.02;        % s, min distance between whisk onsets (findpeaks)
MAXLAG    = 0.07;        % s, max |L-R| lag to count an L/R onset as a pair
SYNC_FRAC = 0.10;        % "synchronous" = |lag| < SYNC_FRAC * mean whisk cycle (ITI)
LAT_WIN   = 0.100;       % s, latency window for "protraction within ___" report
NBINS_LAG = 50;          % L-R lag histogram bins
NBINS_LAT = 40;          % latency histogram bins
fpsW      = 400;         % whisk camera fps
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % thermal_resolve_paths
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

% ---------- load whisking, build protraction-positive BP angles ----------
M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', char(animal), kRun)), 0.6);  % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));   % LEFT  (x mirrored)
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));    % RIGHT
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
tW = (0:numel(xL)-1)'/fpsW;

% protraction onsets = troughs of BP angle (start of forward sweep)
[~,iL] = findpeaks(-xL, 'MinPeakProminence', PROT_PROM*std(xL), 'MinPeakDistance', round(MINDIST*fpsW));
[~,iR] = findpeaks(-xR, 'MinPeakProminence', PROT_PROM*std(xR), 'MinPeakDistance', round(MINDIST*fpsW));
tProtL = (iL-1)/fpsW;  tProtR = (iR-1)/fpsW;
tProt  = sort([tProtL; tProtR]);                       % L+R pooled (for breath-whisk panels)

% ---------- load breathing, detect inspiration onsets ----------
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', char(animal), kRun)), dataRoot);
assert(isfile(Pn.breath), 'no _breath.mat for %s n%d (run drawROI_N_lpsub first)', char(animal), kRun);
Bs  = load(Pn.breath);  br = Bs.breath(:);  fpsB = double(Bs.fps);
[bb,ab] = butter(2, BR_BP/(fpsB/2), 'bandpass');
brf  = filtfilt(bb,ab, fillmissing(br,'linear'));
[~,iInsp] = findpeaks(-brf, 'MinPeakProminence', BR_PROM*std(brf), 'MinPeakDistance', round(fpsB/BR_BP(2)));
tInsp = (iInsp-1)/fpsB;
tBr   = (0:numel(br)-1)'/fpsB;

% ======================================================================
%  BILATERAL: pair each L peak to nearest R peak within MAXLAG
% ======================================================================
maxLag = MAXLAG;
lag_ms = []; lag_t = [];                                % signed lag (L-R) + L peak time
for k = 1:numel(tProtL)
    [d,j] = min(abs(tProtR - tProtL(k)));
    if d <= maxLag
        lag_ms(end+1,1) = 1000*(tProtL(k)-tProtR(j)); %#ok<AGROW>
        lag_t(end+1,1)  = tProtL(k);                  %#ok<AGROW>
    end
end
dtL = diff(tProtL); dtR = diff(tProtR);
iti = [dtL(dtL<0.5); dtR(dtR<0.5)];                    % whisk cycle (drop long gaps)
meanITIms = 1000*mean(iti,'omitnan');
sync_ms   = SYNC_FRAC*meanITIms;
medLag = median(lag_ms); iqrLag = iqr(lag_ms); fracSync = mean(abs(lag_ms)<sync_ms);

% ======================================================================
%  BREATH <-> WHISK: per inspiration cycle [tInsp(i), tInsp(i+1))
% ======================================================================
relProt = {};       % protraction onsets relative to inspiration onset (per breath)
cycDur  = [];       % breath cycle duration (= next insp onset)
firstLat = [];      % latency to first protraction onset after inspiration onset (s)
nPerBr   = [];      % protraction onsets per breath cycle
for i = 1:numel(tInsp)-1
    t0 = tInsp(i); t1 = tInsp(i+1);
    rel = tProt(tProt>=t0 & tProt<t1) - t0;
    relProt{end+1,1} = rel(:)';                        %#ok<AGROW>
    cycDur(end+1,1)  = t1 - t0;                        %#ok<AGROW>
    nPerBr(end+1,1)  = numel(rel);                     %#ok<AGROW>
    if ~isempty(rel), firstLat(end+1,1) = min(rel); else, firstLat(end+1,1) = NaN; end %#ok<AGROW>
end
nBr = numel(relProt);
% raster scatter (y = breath number, temporal order)
Xr = []; Yr = [];
for r = 1:nBr, Xr=[Xr relProt{r}]; Yr=[Yr r*ones(1,numel(relProt{r}))]; end %#ok<AGROW>
latMs   = 1000*firstLat;
medLat  = median(latMs,'omitnan'); iqrLat = iqr(latMs(~isnan(latMs)));
fracLat = mean(firstLat <= LAT_WIN);                   % fraction of breaths with prot. within LAT_WIN

% ---------- secondary circular-phase stats (console only) ----------
% (A) respiratory phase of protraction onsets: insp onset = 0, next = 2*pi
phi = nan(numel(tProt),1);
for q = 1:numel(tProt)
    j = find(tInsp<=tProt(q),1,'last');
    if ~isempty(j) && j<numel(tInsp)
        phi(q) = 2*pi*(tProt(q)-tInsp(j))/(tInsp(j+1)-tInsp(j));
    end
end
phi = phi(~isnan(phi));
PLV_br = abs(mean(exp(1i*phi)));  muPhi = mod(angle(mean(exp(1i*phi))),2*pi);
% (B) L-R Hilbert phase difference (over whisk epochs implicitly via amplitude weight)
phL = angle(hilbert(xL)); phR = angle(hilbert(xR));
dphi = phL - phR;  PLV_lr = abs(mean(exp(1i*dphi)));  muLR = angle(mean(exp(1i*dphi)));

fprintf('\n%s n%d  | whisk %.0fs @%gHz, breath %.0fs @%.3gHz\n', ...
    char(animal), kRun, tW(end), fpsW, tBr(end), fpsB);
fprintf('protraction onsets: L %d, R %d (pooled %d) | inspirations %d (%d cycles)\n', ...
    numel(tProtL), numel(tProtR), numel(tProt), numel(tInsp), nBr);
fprintf('--- bilateral L-R synchrony ---\n');
fprintf('  mean whisk cycle %.1f ms -> sync window +/-%.1f ms (%.0f%% of cycle)\n', meanITIms, sync_ms, 100*SYNC_FRAC);
fprintf('  %d paired onsets | median lag %.1f ms, IQR %.1f ms, sync %.0f%%\n', numel(lag_ms), medLag, iqrLag, 100*fracSync);
fprintf('  L-R Hilbert PLV %.3f, mean offset %.1f deg\n', PLV_lr, rad2deg(muLR));
fprintf('--- breath -> whisk ---\n');
fprintf('  first-protraction latency: median %.1f ms, IQR %.1f ms | within %d ms: %.0f%%\n', ...
    medLat, iqrLat, round(1000*LAT_WIN), 100*fracLat);
fprintf('  whisks/breath: median %.1f (mean %.2f)\n', median(nPerBr), mean(nPerBr));
fprintf('  protraction respiratory-phase PLV %.3f, mean phase %.0f deg\n\n', PLV_br, rad2deg(muPhi));

% ======================================================================
%  FIGURE  (2 x 3 event-timing panels)
% ======================================================================
figure('Color','w','Position',[60 80 1320 760]);
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s n%d  —  breathing<->whisking coordination & L-R whisk synchrony', char(animal), kRun), ...
      'Interpreter','none','FontWeight','bold');
cL = [0 0.5 0]; cR = [0 0.4 0.85];

% (1) L/R fast traces + detected protraction onsets
nexttile; hold on; grid on;
hl=plot(tW, xL, '-', 'Color',cL, 'LineWidth',0.6);
hr=plot(tW, xR, '-', 'Color',cR, 'LineWidth',0.6);
plot(tProtL, xL(iL), 'v', 'Color',cL, 'MarkerFaceColor',cL, 'MarkerSize',4);
plot(tProtR, xR(iR), 'v', 'Color',cR, 'MarkerFaceColor',cR, 'MarkerSize',4);
xlabel('time (s)'); ylabel('whisk angle (deg, BP)');
legend([hl hr], {'L','R'}, 'Location','northeast'); box on;
if tW(end)>8, xlim([0 8]); end                         % zoom for legibility
title('(1) L/R fast whisk + protraction onsets');

% (2) L-R onset-lag histogram (main sync result)  [+lag = L after R]
nexttile; hold on; grid on;
histogram(lag_ms, NBINS_LAG, 'Normalization','probability', 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
xline(0,'k--','LineWidth',1); xline(sync_ms,'r--','LineWidth',1.1); xline(-sync_ms,'r--','LineWidth',1.1);
xlabel('L onset - R onset lag (ms)   [+ = L after R]'); ylabel('fraction of pairs');
title(sprintf('(2) L-R lag  med %.1f ms, IQR %.1f, sync %.0f%%', medLag, iqrLag, 100*fracSync));

% (3) cycle-by-cycle lag over time
nexttile; hold on; grid on;
plot(lag_t, lag_ms, '.-', 'Color',[0.25 0.25 0.25], 'MarkerSize',8);
yline(0,'k--'); yline(sync_ms,'r--'); yline(-sync_ms,'r--');
xlabel('time (s)'); ylabel('L-R lag (ms)');
title('(3) L-R lag over time');

% (4) inspiration-triggered protraction raster (main coordination plot)
nexttile; hold on;
plot(Xr, Yr, '.', 'Color',[0.10 0.55 0.55], 'MarkerSize',4);   % protraction onsets (teal)
plot(cycDur, 1:nBr, '.', 'Color',[0.85 0.1 0.1], 'MarkerSize',4); % next inspiration onset (red)
xline(0,'r-','LineWidth',1);
xmax = prctile(cycDur,99); if isempty(xmax)||xmax<=0, xmax=1; end
xlim([0 xmax]); ylim([1 max(nBr,2)]); set(gca,'YDir','reverse'); box on;
xlabel('time from inspiration onset (s)   [+ = after]'); ylabel('breath number');
title(sprintf('(4) insp-triggered protraction raster  (%d breaths)', nBr));

% (5) protraction latency histogram (first prot. onset after inspiration onset)
nexttile; hold on; grid on;
histogram(latMs, NBINS_LAT, 'Normalization','probability', 'FaceColor',[0.2 0.5 0.7], 'EdgeColor','none');
xline(1000*LAT_WIN,'r--','LineWidth',1.1); xline(medLat,'k-','LineWidth',1.1);
xlim([0 prctile(latMs(~isnan(latMs)),99)]);
xlabel('latency to first protraction onset (ms)'); ylabel('fraction of breaths');
title(sprintf('(5) latency  med %.0f ms, IQR %.0f | <%dms: %.0f%%', medLat, iqrLat, round(1000*LAT_WIN), 100*fracLat));

% (6) whisks per breath histogram
nexttile; hold on; grid on;
edges = -0.5:1:(max(nPerBr)+0.5);
histogram(nPerBr, edges, 'Normalization','probability', 'FaceColor',[0.5 0.4 0.2], 'EdgeColor','none');
xline(median(nPerBr),'k-','LineWidth',1.1);
xlabel('protraction onsets per breath'); ylabel('fraction of breaths');
set(gca,'XTick',0:max(nPerBr));
title(sprintf('(6) whisks per breath  median %.1f', median(nPerBr)));

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
