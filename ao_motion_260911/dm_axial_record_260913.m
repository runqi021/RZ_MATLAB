function R = dm_axial_record_260913()
%DM_AXIAL_RECORD_260913  Acquire a recording with DM axial correction running.
%
% THE PRODUCTION SCRIPT. Everything else in this folder measures the loop;
% this one uses it. Correction is ON for the whole recording, ScanImage logs
% the frames to disk as usual, and the loop's own log is saved beside them so
% the correction is reconstructable afterwards.
%
%   reference stack -> start logging -> for each frame: register, send -dz
%
% ============================== WHAT IT WRITES ==============================
%   <stem>_00001.tif        the science data, written by ScanImage
%   <stem>_aoloop.mat       dz, command, frame numbers, timestamps, reference
%
% The .mat is not optional bookkeeping. Corrected data has an unknown amount of
% motion removed, and dz(t) is the only record of what was taken out. Analysis
% that assumes a fixed focal plane needs it.
%
% ========================= BEFORE YOU RUN THIS =========================
% 1. GEOMETRY. Frame rate sets the loop delay and delay is the only lever.
%    dm_loop_geometry_260913 sets it. Predicted median removal on the ventral
%    archive: 63% at 512x256/58 Hz, 76% with the gx fix, 86% at 128 lines.
% 2. FOV. The estimator needs structure that CHANGES with z. A vessel running
%    down the z axis looks identical at every plane and carries zero axial
%    information -- measured axial SNR 0.000. The check below refuses a field
%    that fails that test rather than recording 100 s of nothing.
% 3. CHANNEL. Register on a STRUCTURAL channel. On a functional one,
%    breath-locked activity is indistinguishable from breath-locked motion:
%    15% of cells firing put 0.22 um into dz, 40% put in 0.38, and the robust
%    refit makes it WORSE at that point.
% 4. The reference must be taken at the SAME geometry as the recording. A
%    256-line reference cannot register 128-line frames.
%
% ============================== SAFETY ==============================
% DRY_RUN is true by default: the reference is taken and the loop runs, but no
% correction reaches the mirror and nothing is logged, so you can read the
% numbers first. Commands are clamped and a command that would breach the ALPAO
% limit is SKIPPED, not clipped -- a clipped mirror is not producing the
% defocus you asked for. onCleanup restores the DM, the stage z and every
% ScanImage setting, on error or Ctrl-C.
%
% *** NOT YET RUN ON THE RIG. *** Written against the ScanImage 2018b source
% and the estimator and loop proven on 2026-09-13. Read it before you trust it.
%
% Runqi Zhang / 2026-09-13

%% ========================= USER SETTINGS =========================
DRY_RUN      = true;      % <<< false to actually correct and log

DURATION_S   = 100;       % length of the recording
OUT_ROOT     = 'D:\Data\RZ\ao_axial_260913';
FILE_STEM    = '';        % '' = auto from the timestamp

SI_DIR       = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE     = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';   % load your system correction FIRST
ROW_DEFOC    = 3;         % verified from actuator geometry: row 3 is defocus
% *** TWO DIFFERENT CHANNEL SETTINGS. THEY ARE NOT THE SAME THING. ***
% SAVE_CHANNELS is your science data: whatever you want recorded, any
% combination. REG_CHANNEL is the single channel the loop registers on, and it
% must be STRUCTURAL. On a functional channel breath-locked ACTIVITY is
% indistinguishable from breath-locked MOTION -- 15% of cells firing puts
% 0.22 um into dz, 40% puts in 0.38, and the robust refit is WORSE at that
% point because a transient covers a whole cell rather than a few outlier
% pixels. Registering on GCaMP would feed the loop the very signal you are
% imaging.
SAVE_CHANNELS = [1 3];    % logged to the TIFF. Any combination of 1 2 3.
REG_CHANNEL   = 3;        % the loop registers on this one. Structural.

% ---- reference ----------------------------------------------------------
REF_REUSE_FROM = '';      % a previous session dir, or '' to take a fresh one
REUSE_Z_TOL_UM = 1.0;     % refuse a reused reference if the stage z has moved
                          % more than this: its zero would not be your zero
% *** SPAN THE MOTION, NOT THE STROKE. *** Measured 2026-09-13, same mirror,
% rate, read path and estimator, 5 um at 0.5 Hz:
%     reference -10..+10 (21 planes) -> estimator gain 0.938
%     reference -14..+24 (39 planes) -> estimator gain 0.653
% The full-DM-authority stack COSTS 35% of the signal -- its far planes carry no
% usable structure and the coarse tier occasionally locks onto one. Too narrow
% fails the other way: 10 um of motion in a +-10 um reference measured 0.810,
% the estimate saturating at the edges.
REF_MIN_UM   = -12;
REF_MAX_UM   =  12;
REF_STEP_UM  =  1;
% *** AVERAGE WHOLE BREATH CYCLES PER PLANE, NOT A FIXED FRAME COUNT. ***
% A reference plane is only the plane you think it is if the motion averaged out
% while it was taken. A short grab is a SNAPSHOT AT A RANDOM BREATH PHASE, with
% the brain displaced by up to the full motion amplitude -- and a reference built
% from 25 differently-displaced planes is wrong in a way nothing downstream can
% detect. Median breath period on the ventral archive is 1067 ms, so 3.3 s is
% about 3 cycles. dm_loop_validate uses 2.0 s and says so explicitly: that is a
% STATIC-sample number, purely about SNR, and it does not transfer to an animal.
REF_SEC_PER_PLANE = 3.3;  % >= 3 breath cycles. 2.0 is enough on a static prep.

% ---- control law --------------------------------------------------------
LOOP_GAIN    = 0.7;       % fraction of the estimate applied per update
EXTRAP_GAIN  = 1.0;       % FRAMES TO LEAD BY. 0 = react only (always late);
                          % 1.0 = theory; ~1.5 = the delay measured 2026-09-13.
                          % Sweep it with dm_loop_validate before trusting a value.
MAX_STEP_UM  = 2.0;       % rate limit: bounds one bad estimate to one step
DEFOCUS_MIN_UM = -12;     % absolute clamp, inside the reference span
DEFOCUS_MAX_UM =  12;
CMD_LIMIT    = 0.95;      % ALPAO command ceiling; a step past it is skipped

RES_ARM_SEC  = 1.0;       % resonant scanner spin-up before the first frame
%% =================================================================

%% ========================= CONNECT =========================
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;
if ~strcmpi(hSI.acqState,'idle')
    error('dm_rec:notIdle','ScanImage is not idle -- stop the current acquisition.');
end

Z = load(Z2C_FILE); fz = fieldnames(Z); Z2C = Z.(fz{1});
assert(size(Z2C,2)==97, 'expected 97 actuators, got %d', size(Z2C,2));
baseline = get_baseline(BASELINE_VAR);

STAMP = datestr(now,'yymmdd_HHMM');
if isempty(FILE_STEM), FILE_STEM = sprintf('aorec_%s', STAMP); end
sessDir = fullfile(OUT_ROOT, FILE_STEM);
if ~DRY_RUN && ~isfolder(sessDir), mkdir(sessDir); end

siState = snapshot_si(hSI);
p0      = hSI.hMotors.motorPosition;
restore = onCleanup(@() restore_all(hDM, baseline, hSI, siState, p0)); %#ok<NASGU>

% APPLY THE CHANNEL SELECTION, AND MAKE SURE THE LOOP CAN SEE ITS CHANNEL.
% channelSave and channelDisplay are independent: the loop reads frames from the
% DISPLAY, so a registration channel that is saved but not displayed is invisible
% to it and every estimate comes back empty.
avail = 1:hSI.hChannels.channelsAvailable;
bad = setdiff(SAVE_CHANNELS, avail);
assert(isempty(bad), 'dm_rec:chan', 'channel(s) %s do not exist; this rig has %d', ...
       mat2str(bad), hSI.hChannels.channelsAvailable);
assert(ismember(REG_CHANNEL, avail), 'dm_rec:regChan', ...
       'REG_CHANNEL %d does not exist', REG_CHANNEL);
hSI.hChannels.channelSave = SAVE_CHANNELS;
if ~ismember(REG_CHANNEL, hSI.hChannels.channelDisplay)
    hSI.hChannels.channelDisplay = union(hSI.hChannels.channelDisplay, REG_CHANNEL);
    fprintf('   [note] added ch%d to channelDisplay so the loop can read it\n', REG_CHANNEL);
end

rate  = hSI.hRoiManager.scanFrameRate;
nFr   = round(DURATION_S * rate);
fprintf('\n=============== AO-CORRECTED RECORDING ===============\n');
fprintf('geometry   %d x %d px, %.1f Hz -> loop delay %.1f ms at 1.46 frames\n', ...
        hSI.hRoiManager.pixelsPerLine, hSI.hRoiManager.linesPerFrame, rate, 1460/rate);
fprintf('duration   %.0f s = %d frames\n', DURATION_S, nFr);
fprintf('register   ch%d | loop gain %.2f | lead %.2f frames | step limit %.1f um\n', ...
        REG_CHANNEL, LOOP_GAIN, EXTRAP_GAIN, MAX_STEP_UM);
fprintf('baseline   %s, max|cmd| %.3f (leaves %.3f of the +-1 stroke)\n', ...
        BASELINE_VAR, max(abs(baseline)), CMD_LIMIT-max(abs(baseline)));
bytes = nFr * hSI.hRoiManager.pixelsPerLine * hSI.hRoiManager.linesPerFrame * 2 * ...
        numel(SAVE_CHANNELS);
fprintf('save       ch%s -> ~%.1f GB\n', mat2str(SAVE_CHANNELS), bytes/1e9);
if DRY_RUN, fprintf('\n*** DRY RUN: no correction sent, nothing logged ***\n'); end

%% ================= 1. FOV CHECK: does z carry information? =================
% A field whose structure does not change with z cannot be registered in z, and
% no amount of analysis recovers it. Two planes 5 um apart must actually differ.
fprintf('\n--- 1. can this field be registered in z? ---\n');
arm_scanner(hSI, RES_ARM_SEC);
hSI.hChannels.loggingEnable = false;
Ia = grab_at_defocus(hSI, hDM, baseline, Z2C, ROW_DEFOC,  0, 4, REG_CHANNEL, CMD_LIMIT);
Ib = grab_at_defocus(hSI, hDM, baseline, Z2C, ROW_DEFOC,  5, 4, REG_CHANNEL, CMD_LIMIT);
hDM.senddata(baseline); pause(0.2);
r = corr_lin(Ia, Ib);
fprintf('   correlation between planes 5 um apart: %.4f\n', r);
if r > 0.97
    warning('dm_rec:flatField', ...
      ['this field looks nearly IDENTICAL 5 um away (r = %.4f).\n' ...
       'Vessels running along z carry ZERO axial information however bright they\n' ...
       'are -- measured axial SNR 0.000. Move to a field with vessels crossing\n' ...
       'the plane, or branch points. 10 s to Ctrl-C.'], r);
    pause(10);
else
    fprintf('   ok -- the field changes with depth\n');
end

%% ========================= 2. REFERENCE STACK =========================
% The reference defines z = 0, and its zero is wherever the stage was parked
% when it was taken. Record that, and refuse a reused reference from a
% different park -- otherwise the zero shifts silently and dz is meaningless.
if isempty(REF_REUSE_FROM)
    fprintf('\n--- 2. reference stack, stage parked, DM stepping ---\n');
    zRef   = REF_MIN_UM:REF_STEP_UM:REF_MAX_UM;
    nRefFr = max(1, round(REF_SEC_PER_PLANE * rate));
    fprintf('   %d planes x %.1f s (%d frames) = %.0f s of acquisition\n', ...
            numel(zRef), REF_SEC_PER_PLANE, nRefFr, numel(zRef)*REF_SEC_PER_PLANE);
    fprintf('   each plane averages %.1f breath cycles at a 1.07 s period\n', ...
            REF_SEC_PER_PLANE/1.067);
    if REF_SEC_PER_PLANE < 2.0
        warning('dm_rec:refTooShort', ...
          ['%.2f s per plane is under two breath cycles. Each reference plane ' ...
           'will be a snapshot at a random breath phase, displaced by up to the ' ...
           'full motion amplitude. Raise REF_SEC_PER_PLANE unless the prep is static.'], ...
           REF_SEC_PER_PLANE);
    end
    refVol = []; okRef = true(size(zRef));
    for k = 1:numel(zRef)
        [I, sent] = grab_at_defocus(hSI, hDM, baseline, Z2C, ROW_DEFOC, zRef(k), ...
                                    nRefFr, REG_CHANNEL, CMD_LIMIT);
        if ~sent
            fprintf('   z %+5.1f um SKIPPED -- would breach the ALPAO limit\n', zRef(k));
            okRef(k) = false; continue
        end
        if isempty(refVol), refVol = zeros([size(I) numel(zRef)],'single'); end
        refVol(:,:,k) = I; %#ok<AGROW>
        if mod(k,5)==0 || k==numel(zRef), fprintf('   plane %d/%d\n', k, numel(zRef)); end
    end
    hDM.senddata(baseline); pause(0.2);
    % NEVER keep a plane that was not acquired: it is zeros, and the estimator
    % would happily match a frame to it.
    refVol = refVol(:,:,okRef); zRef = zRef(okRef);
    if any(~okRef)
        fprintf('   dropped %d unreachable plane(s); reference spans %+.1f to %+.1f um\n', ...
                nnz(~okRef), zRef(1), zRef(end));
    end
    refZPark = hSI.hMotors.motorPosition(3);
    refGeom  = [hSI.hRoiManager.pixelsPerLine hSI.hRoiManager.linesPerFrame];
else
    fprintf('\n--- 2. reusing reference from %s ---\n', REF_REUSE_FROM);
    Rr = load(fullfile(REF_REUSE_FROM,'reference.mat'));
    refVol = Rr.refVol; zRef = Rr.zRef; refZPark = Rr.refZPark; refGeom = Rr.refGeom;
    dz0 = hSI.hMotors.motorPosition(3) - refZPark;
    if abs(dz0) > REUSE_Z_TOL_UM
        error('dm_rec:refZ', ['reference was taken at stage z = %.2f, stage is now %.2f ' ...
          '(%.2f um away). Its zero is not your zero -- take a fresh one.'], ...
          refZPark, hSI.hMotors.motorPosition(3), dz0);
    end
    if ~isequal(refGeom, [hSI.hRoiManager.pixelsPerLine hSI.hRoiManager.linesPerFrame])
        error('dm_rec:refGeom', ['reference is %dx%d, scan is now %dx%d. A reference ' ...
          'cannot register frames of a different size.'], refGeom(1), refGeom(2), ...
          hSI.hRoiManager.pixelsPerLine, hSI.hRoiManager.linesPerFrame);
    end
end

E = dm_axial_estimator_260913('prep', refVol, zRef, struct('nStrips',1,'pixFrac',0.04));
mid = refVol(:,:,max(1,round(numel(zRef)/2)));
for q = 1:20, dm_axial_estimator_260913('run', E, single(mid), 0); end   % warm up JIT
fprintf('   reference ready: %d planes, %+.1f to %+.1f um, parked at z = %.2f\n', ...
        numel(zRef), zRef(1), zRef(end), refZPark);
if ~DRY_RUN
    save(fullfile(sessDir,'reference.mat'),'refVol','zRef','refZPark','refGeom','-v7.3');
end

%% ========================= 3. THE RECORDING =========================
fprintf('\n--- 3. recording %.0f s with correction %s ---\n', DURATION_S, ...
        ternary(DRY_RUN,'SIMULATED (dry run)','ON'));
hSI.hStackManager.numSlices      = 1;
hSI.hStackManager.framesPerSlice = nFr;
hSI.hChannels.loggingEnable      = ~DRY_RUN;
if ~DRY_RUN
    hSI.hScan2D.logFilePath   = sessDir;
    hSI.hScan2D.logFileStem   = FILE_STEM;
    hSI.hScan2D.logFileCounter = 1;
end
arm_scanner(hSI, RES_ARM_SEC);
hDM.senddata(baseline); pause(0.2);

L = struct('t',nan(1,nFr),'frameNo',nan(1,nFr),'dz',nan(1,nFr),'cmdUm',nan(1,nFr), ...
           'nSvc',0,'nMissed',0,'nSkipCmd',0,'nEmpty',0);
corrUm = 0; uPrev = NaN; fnoPrev = -1; lastNo = -1; nEmpty = 0;
t0 = tic;
hSI.startGrab();
if strcmpi(hSI.acqState,'idle')
    error('dm_rec:grabDidNotStart', ...
      ['startGrab left acqState idle -- the grab never started. Most likely no ' ...
       'channel is selected for Display or Save in ScanImage.']);
end

while ~strcmpi(hSI.acqState,'idle')
    try, fno = hSI.hDisplay.lastFrameNumber; catch, pause(0.002); continue; end
    if fno == lastNo, pause(0.001); drawnow limitrate; continue; end
    [I, ~] = get_disp_frame(hSI, REG_CHANNEL);
    if isempty(I)
        nEmpty = nEmpty + 1;
        if nEmpty > 2000, hSI.abort(); error('dm_rec:noFrames','display gave no frames'); end
        pause(0.002); continue
    end
    nEmpty = 0;
    if lastNo >= 0, L.nMissed = L.nMissed + (fno - lastNo - 1); end
    lastNo = fno;
    L.nSvc = L.nSvc + 1;
    if L.nSvc > nFr, break; end

    dzNow = single_dz(E, I);

    if isfinite(dzNow)
        % EXTRAPOLATE THE DISTURBANCE, NOT THE RESIDUAL, and normalise the
        % velocity by frames actually elapsed. See dm_loop_validate for why:
        % the residual's velocity cancels when tracking works, and a dropped
        % frame otherwise doubles the lead exactly when data is sparse.
        uStar = corrUm - dzNow;
        dn    = max(1, fno - fnoPrev);
        if isfinite(uPrev), pred = uStar + EXTRAP_GAIN*(uStar - uPrev)/dn;
        else,               pred = uStar; end
        uPrev = uStar; fnoPrev = fno;

        want = corrUm + LOOP_GAIN*(pred - corrUm);
        want = max(DEFOCUS_MIN_UM, min(DEFOCUS_MAX_UM, want));
        want = corrUm + max(-MAX_STEP_UM, min(MAX_STEP_UM, want - corrUm));

        cmd = baseline(:).' + (want/gain_um_per_coeff()) * Z2C(ROW_DEFOC,:);
        if max(abs(cmd)) > CMD_LIMIT
            L.nSkipCmd = L.nSkipCmd + 1;        % hold the last good command
        else
            if ~DRY_RUN, hDM.senddata(cmd); end
            corrUm = want;
        end
    end

    L.t(L.nSvc) = toc(t0);  L.frameNo(L.nSvc) = fno;
    L.dz(L.nSvc) = dzNow;   L.cmdUm(L.nSvc) = corrUm;
    if mod(L.nSvc, round(rate*10)) == 0
        fprintf('   %5.0f s  dz rms %5.2f um  cmd %+6.2f um  missed %d\n', ...
                L.t(L.nSvc), rms_nan(L.dz(1:L.nSvc)), corrUm, L.nMissed);
    end
end
hDM.senddata(baseline);

%% ========================= 4. REPORT =========================
n = L.nSvc;
fprintf('\n--- 4. done ---\n');
fprintf('   frames serviced %d of %d (%.1f%%), missed %d, empty %d\n', ...
        n, nFr, 100*n/nFr, L.nMissed, L.nEmpty);
fprintf('   residual dz     %.3f um rms, %.3f um p-p\n', ...
        rms_nan(L.dz(1:n)), max(L.dz(1:n))-min(L.dz(1:n)));
fprintf('   command applied %.3f um rms, %+.2f to %+.2f um\n', ...
        rms_nan(L.cmdUm(1:n)), min(L.cmdUm(1:n)), max(L.cmdUm(1:n)));
if L.nSkipCmd > 0
    fprintf('   !! %d command(s) skipped at the ALPAO limit -- the correction was\n', L.nSkipCmd);
    fprintf('      held, not applied. Narrow DEFOCUS_MIN/MAX or re-centre the baseline.\n');
end
if L.nMissed > 0.02*n
    fprintf('   !! %.1f%% of frames missed -- the loop is not keeping up. The lead\n', 100*L.nMissed/n);
    fprintf('      term divides by elapsed frames so this is handled, but it costs correction.\n');
end
R = struct('log',L,'zRef',zRef,'refZPark',refZPark,'refGeom',refGeom, ...
           'baselineVar',BASELINE_VAR,'loopGain',LOOP_GAIN,'extrapGain',EXTRAP_GAIN, ...
           'saveChannels',SAVE_CHANNELS,'regChannel',REG_CHANNEL,'rate',rate,'nFrames',nFr,'stem',FILE_STEM,'dryRun',DRY_RUN);
if ~DRY_RUN
    out = fullfile(sessDir, [FILE_STEM '_aoloop.mat']);
    save(out,'-struct','R');
    fprintf('\n   science data  %s\n   loop log      %s\n', ...
            fullfile(sessDir,[FILE_STEM '_00001.tif']), out);
else
    fprintf('\n   dry run -- nothing written. Set DRY_RUN = false to record.\n');
end
end

%% ============================ HELPERS ============================
function g = gain_um_per_coeff()
% um of focus per um of Zernike defocus coefficient. MEASURED on the rig
% 2026-09-11: calibration gave 6.656, four stack demos implied 6.76-7.66
% (mean 7.04). Quote ~7 um/um +-5%. Re-measure with dm_motion_calib if the
% objective, the immersion or the system correction changes.
g = 6.9;
end

function [I, sent] = grab_at_defocus(hSI, hDM, baseline, Z2C, row, um, nAvg, ch, lim)
% Park the DM at a commanded defocus and average nAvg frames there.
cmd = baseline(:).' + (um/gain_um_per_coeff()) * Z2C(row,:);
sent = max(abs(cmd)) <= lim;
if ~sent, I = []; return; end
hDM.senddata(cmd); pause(0.25);
I = quick_avg(hSI, nAvg, ch);
end

function I = quick_avg(hSI, n, ch)
hSI.hStackManager.numSlices = 1;
hSI.hStackManager.framesPerSlice = n;
hSI.startGrab();
while ~strcmpi(hSI.acqState,'idle'), pause(0.005); drawnow limitrate; end
I = get_disp_frame(hSI, ch);
end

function [I, why] = get_disp_frame(hSI, ch)
I = []; why = '';
try
    f = hSI.hDisplay.lastFrame;
    if isempty(f), why = 'empty buffer'; return; end
    idx = find(hSI.hDisplay.lastFrameChannels == ch, 1);
    if isempty(idx), why = sprintf('channel %d not displayed', ch); return; end
    I = single(f{idx});
catch ME
    why = ME.message;
end
end

function dz = single_dz(E, I)
o = dm_axial_estimator_260913('run', E, single(I), NaN);   % NaN = coarse then fine
dz = o.dz(1);
end

function b = get_baseline(varName)
% No way to read the current command back off the ALPAO driver, so it comes
% from the base workspace by name -- load your system correction FIRST.
try
    b = evalin('base', varName);
catch
    error('dm_rec:noBaseline', ...
      ['%s is not in the base workspace. Load your system aberration correction ' ...
       'first, e.g.\n   load(''...SystemCorrection_20260117_930_..._water.mat'')'], varName);
end
b = double(b(:)).';
assert(numel(b)==97, '%s is %d long, expected 97', varName, numel(b));
end

function arm_scanner(hSI, sec)
% Hold the resonant scanner on. arm() only ESTIMATES the scan frequency -- the
% measure-after-settle branch is commented out in SI -- so there is no settle
% wait in the grab path and a cold scanner distorts the first frames.
try, hSI.hScan2D.keepResonantScannerOn = true; catch, end
t = tic; while toc(t) < sec, pause(0.05); drawnow limitrate; end
end

function s = snapshot_si(hSI)
s.numSlices      = hSI.hStackManager.numSlices;
s.framesPerSlice = hSI.hStackManager.framesPerSlice;
s.logging        = hSI.hChannels.loggingEnable;
s.chanSave       = hSI.hChannels.channelSave;
s.chanDisplay    = hSI.hChannels.channelDisplay;
s.logFilePath    = hSI.hScan2D.logFilePath;
s.logFileStem    = hSI.hScan2D.logFileStem;
s.logFileCounter = hSI.hScan2D.logFileCounter;
try, s.keepOn    = hSI.hScan2D.keepResonantScannerOn; catch, s.keepOn = false; end
end

function restore_all(hDM, baseline, hSI, s, p0)
try, hDM.senddata(baseline); fprintf('[cleanup] DM back to baseline\n');
catch ME, warning('dm_rec:dmRestore','COULD NOT RESTORE THE DM: %s', ME.message); end
try
    hSI.hStackManager.numSlices      = s.numSlices;
    hSI.hStackManager.framesPerSlice = s.framesPerSlice;
    hSI.hChannels.loggingEnable      = s.logging;
    hSI.hChannels.channelSave        = s.chanSave;
    hSI.hChannels.channelDisplay     = s.chanDisplay;
    hSI.hScan2D.logFilePath          = s.logFilePath;
    hSI.hScan2D.logFileStem          = s.logFileStem;
    hSI.hScan2D.logFileCounter       = s.logFileCounter;
    hSI.hScan2D.keepResonantScannerOn = s.keepOn;
catch ME, warning('dm_rec:siRestore','could not restore ScanImage settings: %s', ME.message); end
try
    p = hSI.hMotors.motorPosition;
    if abs(p(3)-p0(3)) > 0.5, hSI.hMotors.motorPosition = p0; end
catch, end
end

function r = corr_lin(a, b)
a = double(a(:)); b = double(b(:));
a = a - mean(a); b = b - mean(b);
r = (a.'*b) / (norm(a)*norm(b) + eps);
end

function v = rms_nan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = sqrt(mean(x.^2)); end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
