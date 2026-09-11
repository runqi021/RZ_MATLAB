function whisk_breath_onset_check_RZ()
% whisk_breath_onset_check_RZ
%
% QC / sanity check for the onset detection used by whisk_breath_raster_RZ.
% For ONE session it overlays, on a shared time axis over a whisking window:
%
%   row 1 : RAW nostril temperature (L+R avg, deg C) + slow LP baseline.
%           PHYSICS: inhalation draws cool air -> raw temp FALLS during the
%           inhale. The onset is the START of that fall, so the onset line
%           sits at a raw-temp LOCAL MAX (warmest, end of exhale); raw then
%           dips toward the next trough. (saved breath = -(raw - baseline),
%           verified corr = +1.000, so its trough = raw peak = insp onset.)
%
%   row 2 : inhale-up breath (the saved B.breath, BP 1-20 Hz) with detected
%           INSPIRATION ONSETS (red). Onset = trough of the inhale-up signal
%           = start of the rise. These should sit at the bottom of each cycle,
%           lined up with the START of each dip in row 1.
%
%   row 3 : LEFT whisker BP angle with PROTRACTION ONSETS (markers).
%   row 4 : RIGHT whisker BP angle with PROTRACTION ONSETS (markers).
%           Protraction = increasing angle, so onset = trough (most retracted).
%
% Detection here is IDENTICAL to whisk_breath_raster_RZ. Nothing is written;
% this only draws figures and prints the per-session inverted flag.

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

SESS = {"5916296", 1};   % {animal, run} to inspect
WIN_S = 6;               % seconds to display (centered in the longest whisk epoch)

fpsW = 400;
WHISK_BP  = [5 40];
BREATH_BP = [1 20];

ENV_THR = 5; ENV_SMOOTH = 0.05; MIN_DUR = 0.5; MERGE_GAP = 0.1;
WHISK_PROM = 1; MIN_WHISK_DIST_S = 0.04;
MIN_INSP_DIST_S = 0.06; BREATH_PROM_MAD = 1;
PROTRACTION_POLARITY = -1;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
ord = 4;
[b_wbp,a_wbp] = butter(ord, WHISK_BP/(fpsW/2),'bandpass');
[b_br, a_br ] = butter(ord, BREATH_BP/(fpsW/2),'bandpass');

animal = char(SESS{1}); kk = SESS{2};

% ----- whisker -----
M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kk)), 0.6);
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
tW = (0:numel(La)-1)'/fpsW; Nw = numel(tW);
La0 = fillmissing(La(:),'linear','EndValues','nearest'); La0 = La0 - mean(La0,'omitnan');
Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest'); Ra0 = Ra0 - mean(Ra0,'omitnan');
xL = filtfilt(b_wbp,a_wbp, La0);
xR = filtfilt(b_wbp,a_wbp, Ra0);

env = movmean(abs(hilbert(xL)), max(1,round(ENV_SMOOTH*fpsW)));
ep  = bool_to_epochs(env>ENV_THR, tW, MIN_DUR, MERGE_GAP);
m   = false(Nw,1); for q=1:size(ep,1), m = m | (tW>=ep(q,1)&tW<=ep(q,2)); end

tProtL = detect_large_protractions(xL,tW,m,MIN_WHISK_DIST_S,WHISK_PROM,PROTRACTION_POLARITY);
tProtR = detect_large_protractions(xR,tW,m,MIN_WHISK_DIST_S,WHISK_PROM,PROTRACTION_POLARITY);

% ----- breath (saved inhale-up) + RAW nostril temp for the physics check -----
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kk)), dataRoot);
B  = load(Pn.breath); br = B.breath(:); fb = double(B.fps);
tB = (0:numel(br)-1)'/fb;
invflag = 'NA'; if isfield(B,'inverted'), invflag = mat2str(B.inverted); end
fprintf('%s n%d: inverted=%s, method=%s, fps=%.1f\n', animal, kk, invflag, ...
    char(string(getfield_safe(B,'method','NA'))), fb);

brw = fillmissing(interp1(tB, br, tW, 'linear', NaN),'linear','EndValues','nearest');
brw_bp = filtfilt(b_br,a_br, brw);
prom = BREATH_PROM_MAD * mad(brw_bp,1); if ~isfinite(prom)||prom<=0, prom = 0.5*std(brw_bp,'omitnan'); end
[~,locs] = findpeaks(-brw_bp, 'MinPeakDistance',round(MIN_INSP_DIST_S*fpsW), 'MinPeakProminence',prom);
tInsp = tW(locs);

% RAW nostril temp (L+R avg) + baseline, for "inhale = dip" check
rawA = []; baseA = []; tRaw = [];
try
    D = load(Pn.nostrilROI);
    fbr = double(D.L.fps);
    [lb,la] = butter(2, 1/(fbr/2),'low');     % 1 Hz baseline as in lpsub
    rL = double(D.L.trace(:)); rR = double(D.R.trace(:));
    rawA = (rL + rR)/2; baseA = filtfilt(lb,la, rawA);
    tRaw = (0:numel(rawA)-1)'/fbr;
catch ME
    warning('raw nostrilROI not available: %s', ME.message);
end

% ----- choose display window: centered in the longest whisk epoch -----
if isempty(ep)
    t0 = 0;
else
    [~,iq] = max(ep(:,2)-ep(:,1)); c = mean(ep(iq,:));
    t0 = max(0, c - WIN_S/2);
end
t1 = t0 + WIN_S;
tInspWin = tInsp(tInsp>=t0 & tInsp<=t1);   % only draw onset ref-lines in the window

% ============================== PLOT ==================================
figure('Color','w','Position',[60 60 1180 900], 'Name', sprintf('%s n%d onset check',animal,kk));
tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
ax = gobjects(4,1);

ax(1) = nexttile; hold on; grid on;
if ~isempty(rawA)
    plot(tRaw, rawA, '-', 'Color',[0.4 0.4 0.4]);
    plot(tRaw, baseA, '-', 'Color',[0.85 0.1 0.1], 'LineWidth',1.4);
    for tt = tInspWin', xline(tt,'b:'); end   % insp onsets as ref lines (window only)
    ylabel('raw nostril \circC');
    title('RAW nostril temp (gray) + 1 Hz baseline (red).  inhale cools -> raw FALLS during inhale; blue insp-onset sits at a raw-temp PEAK, then dips');
else
    title('raw nostril temp unavailable'); ylabel('\circC');
end

ax(2) = nexttile; hold on; grid on;
plot(tW, brw_bp, '-', 'Color',[0.1 0.1 0.1]);
plot(tInsp, brw_bp(locs), '.', 'Color',[0.85 0.1 0.1], 'MarkerSize',14);
yline(0,'k:'); ylabel('inhale-up breath (a.u.)');
title(sprintf('inhale-up breath BP %g-%g Hz; red = inspiration onsets (troughs).  inverted=%s', ...
    BREATH_BP(1),BREATH_BP(2), invflag));

ax(3) = nexttile; hold on; grid on;
plot(tW, xL, '-', 'Color',[0.1 0.1 0.1]);
iL = ismember(tW, tProtL); plot(tW(iL), xL(iL), 'v', 'Color',[0.1 0.55 0.55], 'MarkerFaceColor',[0.1 0.55 0.55], 'MarkerSize',6);
for tt = tInspWin', xline(tt,'r:'); end
ylabel('L whisk angle (deg)'); title('LEFT whisker BP angle; teal v = protraction onsets (troughs); red = insp onsets');

ax(4) = nexttile; hold on; grid on;
plot(tW, xR, '-', 'Color',[0.1 0.1 0.1]);
iR = ismember(tW, tProtR); plot(tW(iR), xR(iR), 'v', 'Color',[0.1 0.55 0.55], 'MarkerFaceColor',[0.1 0.55 0.55], 'MarkerSize',6);
for tt = tInspWin', xline(tt,'r:'); end
ylabel('R whisk angle (deg)'); xlabel('time (s)');
title('RIGHT whisker BP angle; teal v = protraction onsets (troughs); red = insp onsets');

linkaxes(ax,'x'); xlim(ax(4),[t0 t1]);
title(tl, sprintf('%s n%d  onset QC  (window %.1f-%.1f s)', animal, kk, t0, t1), 'Interpreter','none');
end

% ============================= HELPERS =============================
function v = getfield_safe(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end

function tProt = detect_large_protractions(x, t, mask, minDistS, whiskProm, polarity)
    x=x(:); t=t(:); mask=logical(mask(:));
    dt=median(diff(t)); minDistSamples=max(1,round(minDistS/dt));
    good=mask&isfinite(x); prom=whiskProm*std(x(good),'omitnan');
    if ~isfinite(prom)||prom<=0, prom=0.5*std(x(isfinite(x)),'omitnan'); end
    if ~isfinite(prom)||prom<=0, tProt=zeros(0,1); return; end
    [~,idx]=findpeaks(polarity*x,'MinPeakProminence',prom,'MinPeakDistance',minDistSamples);
    idx=idx(mask(idx)); tProt=t(idx);
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a=logical(a(:)); d=diff([false;a;false]); s=find(d==1); e=find(d==-1)-1;
    ep=[t(s),t(e)]; if isempty(ep), return; end
    ep2=ep(1,:);
    for i=2:size(ep,1)
        if ep(i,1)-ep2(end,2)<=mergeGap, ep2(end,2)=ep(i,2); else, ep2(end+1,:)=ep(i,:); end %#ok<AGROW>
    end
    ep=ep2; ep=ep(ep(:,2)-ep(:,1)>=minDur,:);
end

function csv = pick_csv(dirPath, prefix)
    d=dir(fullfile(char(dirPath),[char(prefix) '*DLC*.csv']));
    assert(~isempty(d),'no DLC csv matching %s* in %s',prefix,dirPath);
    [~,ix]=max(arrayfun(@(x) bestnum(x.name),d)); csv=fullfile(d(ix).folder,d(ix).name);
end
function n = bestnum(name)
    tok=regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
