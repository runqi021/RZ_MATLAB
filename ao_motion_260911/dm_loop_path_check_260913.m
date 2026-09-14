function R = dm_loop_path_check_260913()
%DM_LOOP_PATH_CHECK_260913  Why does the loop see only 65% of the motion?
%
% THE CONTRADICTION THIS RESOLVES (rig, 2026-09-13):
%   dm_estimator_check : commanded 5.8 um -> estimated 5.90 um.  Gain 0.99.
%   dm_loop_validate   : injected  5.0 um -> measured   3.3 um.  Gain 0.65.
% Same mirror, same Z2C row, same reference, same conversion constant. So it is
% not calibration. The difference is HOW THE FRAME IS READ:
%
%   estimator_check : senddata -> pause 0.3 s -> startGrab(N) -> read AFTER it ends
%   the loop        : senddata -> read DURING a continuous grab, no settle
%
% Two things differ at once there, and they need different fixes:
%   A. NO SETTLE. The DM may not reach the commanded defocus within one 16.7 ms
%      frame. Fix: slower loop, or accept a known actuator lag.
%   B. STALE OR SMEARED FRAMES. lastFrameNumber comes from lastStripeData;
%      the image comes from rollingStripeDataBuffer. Those are different
%      structures, and if the number advances while the image repeats or
%      blends, dz is smeared and its amplitude drops. Fix: read the image
%      differently, or verify each frame is genuinely new.
%
% A STAIRCASE SEPARATES THEM. Hold each level long enough that the DM has
% certainly settled (A is excluded), but read through the LOOP'S path (B is
% still live). Then:
%   plateaus reach the commanded value -> read path is fine, the 0.65 is the
%                                         DM not settling within a frame (A)
%   plateaus come out ~0.65x           -> the read path is the problem (B)
%
% It also reports how many reads returned a genuinely NEW image, which tests B
% directly and independently of the estimator.
%
% ~25 s. Nothing logged, nothing corrected, DM returned to baseline.
%
% Runqi Zhang / 2026-09-13.

%% ---------------- settings ----------------
SI_DIR       = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE     = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_DEFOC    = 3;
CHANNEL      = 3;

REF_ZMIN_UM  = -10;
REF_ZMAX_UM  =  10;
REF_ZSTEP_UM =  1;
REF_NFRAMES  = 30;

STAIR_UM     = [0 5 0 -5 0 3 0 -3 0];
HOLD_S       = 2.0;      % far longer than any plausible DM settling time
GAIN_UM_PER_COEFF = 6.656;
CMD_LIMIT    = 0.95;
SETTLE_S     = 0.3;

%% ---------------- connect ----------------
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;
S = load(Z2C_FILE);  fn = fieldnames(S);  Z2C = S.(fn{1});
zRow = Z2C(ROW_DEFOC,:);
if ~evalin('base', sprintf('exist(''%s'',''var'')', BASELINE_VAR))
    error('dm_path:noBaseline','base workspace has no %s', BASELINE_VAR);
end
baseline = double(evalin('base', BASELINE_VAR));  baseline = baseline(:).';

st = struct('rollingAvg', hSI.hDisplay.displayRollingAverageFactor, ...
            'framesPerSlice', hSI.hStackManager.framesPerSlice, ...
            'numSlices', hSI.hStackManager.numSlices, ...
            'logging', hSI.hChannels.loggingEnable, ...
            'keepResOn', hSI.hScan2D.keepResonantScannerOn);
restore = onCleanup(@() cleanup(hDM, baseline, hSI, st)); %#ok<NASGU>
hSI.hChannels.loggingEnable = false;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = 1;

fps = hSI.hRoiManager.scanFrameRate;
fprintf('\n============ LOOP READ-PATH CHECK ============\n');
fprintf('%.2f Hz | staircase %s um, %.1f s per level\n', ...
        fps, mat2str(STAIR_UM), HOLD_S);

hSI.hScan2D.keepResonantScannerOn = true;
wait_s(1.0);

%% ---------------- reference ----------------
zRef = REF_ZMIN_UM : REF_ZSTEP_UM : REF_ZMAX_UM;
refVol = [];  got = false(1,numel(zRef));
for k = 1:numel(zRef)
    cmd = baseline + (zRef(k)/GAIN_UM_PER_COEFF)*zRow;
    if max(abs(cmd)) > CMD_LIMIT, continue; end
    hDM.senddata(cmd);  pause(SETTLE_S);
    I = grab_after(hSI, REF_NFRAMES, CHANNEL);
    if isempty(refVol), refVol = zeros([size(I) numel(zRef)],'single'); end
    refVol(:,:,k) = I;  got(k) = true; %#ok<AGROW>
end
hDM.senddata(baseline);  pause(SETTLE_S);
refVol = refVol(:,:,got);  zRef = zRef(got);
E = dm_axial_estimator_260913('prep', refVol, zRef, struct('nStrips',1,'pixFrac',0.04));
fprintf('reference: %d planes, %.1f to %.1f um\n', numel(zRef), zRef(1), zRef(end));

%% ---------------- staircase through the LOOP path ----------------
nStep = numel(STAIR_UM);
nFr = round(nStep*HOLD_S*fps);
hSI.hStackManager.framesPerSlice = nFr;

T = nan(1,nFr); DZ = nan(1,nFr); CMDU = nan(1,nFr); NEWIM = false(1,nFr);
lastNo = -1; nSvc = 0; prevIm = [];

hDM.senddata(baseline + (STAIR_UM(1)/GAIN_UM_PER_COEFF)*zRow);
pause(SETTLE_S);
t0 = tic;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_path:grabDidNotStart','startGrab left acqState idle');
end
curUm = STAIR_UM(1);  iStep = 1;

while ~strcmp(hSI.acqState,'idle')
    try, fno = hSI.hDisplay.lastFrameNumber; catch, pause(0.002); continue; end
    if fno == lastNo, pause(0.001); drawnow limitrate; continue; end
    I = read_loop_path(hSI, CHANNEL);
    if isempty(I), pause(0.0005); continue; end
    lastNo = fno;  nSvc = nSvc + 1;
    if nSvc > nFr, break; end

    tNow = toc(t0);
    want = STAIR_UM(min(nStep, floor(tNow/HOLD_S)+1));
    if want ~= curUm
        hDM.senddata(baseline + (want/GAIN_UM_PER_COEFF)*zRow);
        curUm = want;
    end

    % Is this image genuinely NEW, or did the buffer hand back the same one
    % while the frame NUMBER advanced? That is test B, and it needs no
    % estimator at all.
    if ~isempty(prevIm)
        NEWIM(nSvc) = ~isequal(I, prevIm);
    else
        NEWIM(nSvc) = true;
    end
    prevIm = I;

    o = dm_axial_estimator_260913('run', E, single(I), NaN);
    T(nSvc)=tNow; DZ(nSvc)=o.dz(1); CMDU(nSvc)=curUm;
    drawnow limitrate;
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
hDM.senddata(baseline);  pause(SETTLE_S);

k = isfinite(T);
T=T(k); DZ=DZ(k); CMDU=CMDU(k); NEWIM=NEWIM(k);

%% ---------------- report ----------------
fprintf('\nserviced %d frames; %.0f%% of reads returned a NEW image\n', ...
        numel(T), 100*mean(NEWIM));
if mean(NEWIM) < 0.9
    fprintf('  >> THE READ PATH IS REPEATING FRAMES. %.0f%% of reads handed back\n', ...
            100*(1-mean(NEWIM)));
    fprintf('     the SAME image while lastFrameNumber advanced. That alone\n');
    fprintf('     smears dz and shrinks its amplitude.\n');
end

fprintf('\n%10s %12s %10s %10s\n','commanded','plateau est','err','n');
lv = unique(STAIR_UM);
res = struct('cmd',[],'est',[],'n',[]);
for v = lv
    % Use only the second half of each hold, so DM settling cannot contribute.
    sel = false(size(T));
    for s = 1:nStep
        if STAIR_UM(s) ~= v, continue; end
        t1 = (s-1)*HOLD_S + HOLD_S/2;  t2 = s*HOLD_S;
        sel = sel | (T >= t1 & T < t2);
    end
    sel = sel & isfinite(DZ);
    if nnz(sel) < 5, continue; end
    m = median(DZ(sel));
    fprintf('%10.1f %12.2f %10.2f %10d\n', v, m, m-v, nnz(sel));
    res.cmd(end+1)=v; res.est(end+1)=m; res.n(end+1)=nnz(sel);
end

nz = res.cmd ~= 0;
if nnz(nz) >= 2
    g = mean(res.est(nz) ./ res.cmd(nz));
    fprintf('\nplateau gain (settled, read through the LOOP path): %.3f\n', g);
    if g > 0.85
        fprintf('=> THE READ PATH IS FINE. Settled plateaus come out right, so the\n');
        fprintf('   0.65 in the sweep is the DM NOT REACHING the commanded defocus\n');
        fprintf('   within one 16.7 ms frame. That is an actuator/settling limit:\n');
        fprintf('   slow the loop, or calibrate the dynamic gain and correct for it.\n');
    else
        fprintf('=> THE READ PATH IS THE PROBLEM (gain %.3f even fully settled).\n', g);
        fprintf('   The frames the loop reads do not reflect the commanded state.\n');
        fprintf('   Fix the frame source before touching the controller.\n');
    end
end

R = struct('t',T,'dz',DZ,'cmd',CMDU,'newImage',NEWIM,'plateau',res,'zRef',zRef);
fprintf('==============================================\n');
end

%% ---------------- helpers ----------------
function I = read_loop_path(hSI, ch)
% EXACTLY what the closed loop does -- same source, same order, no settle.
I = [];
try
    b = hSI.hDisplay.rollingStripeDataBuffer;
    if ~isempty(b) && ~isempty(b{1}) && ~isempty(b{1}{1}.roiData)
        I = pick_channel(b{1}{1}.roiData{1}, ch);
        if ~isempty(I), return; end
    end
catch
end
try, I = pick_channel(hSI.hDisplay.lastStripeData.roiData{1}, ch); catch, end
end

function I = grab_after(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_path:grabDidNotStart','startGrab left acqState idle');
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
I = read_loop_path(hSI, ch);
assert(~isempty(I), 'could not read channel %d', ch);
hSI.hDisplay.displayRollingAverageFactor = 1;
end

function I = pick_channel(rd, ch)
I = [];
j = find(rd.channels(:).' == ch, 1);
if isempty(j) || j > numel(rd.imageData), return; end
im = rd.imageData{j};
if iscell(im), if isempty(im), return; end, im = im{1}; end
if isempty(im), return; end
if isprop(rd,'transposed') && rd.transposed, im = im.'; end
I = single(im);
end

function wait_s(s)
t0 = tic; while toc(t0) < s, pause(0.05); drawnow limitrate; end
end

function cleanup(hDM, baseline, hSI, st)
try, hDM.senddata(baseline); catch ME, warning('dm_path:dm','%s',ME.message); end
try
    if ~strcmp(hSI.acqState,'idle')
        hSI.abort();
        while ~strcmp(hSI.acqState,'idle'), pause(0.02); end
    end
    hSI.hDisplay.displayRollingAverageFactor = st.rollingAvg;
    hSI.hStackManager.framesPerSlice = st.framesPerSlice;
    hSI.hStackManager.numSlices = st.numSlices;
    hSI.hChannels.loggingEnable = st.logging;
    hSI.hScan2D.keepResonantScannerOn = st.keepResOn;
    fprintf('[cleanup] restored\n');
catch ME, warning('dm_path:si','%s',ME.message); end
end
