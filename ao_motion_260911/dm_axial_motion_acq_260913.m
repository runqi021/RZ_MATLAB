function dm_axial_motion_acq_260913()
%DM_AXIAL_MOTION_ACQ_260913  Acquire everything needed to measure axial motion.
%
% The go/no-go measurement for DM-based axial motion correction: HOW FAR DOES
% THE BRAIN MOVE IN Z PER BREATH? Nothing is corrected here and no loop is
% closed. This script only ACQUIRES, in four blocks, and every block logs raw
% individual frames to disk so all of the analysis is offline and repeatable.
%
%   A  STAGE z-stack   -- the metric ruler. Encoder-backed, unidirectional.
%   B  DM z-stack      -- the REGISTRATION REFERENCE. Stage parked, defocus only.
%   C  TIME SERIES     -- stage parked, DM at baseline, the animal breathing.
%                         This is the measurement.
%   D  INJECTION       -- a known defocus staircase during a second time series.
%                         Positive control: the estimator must recover a signal
%                         it was told about, on top of breathing.
%
% WHY THE REFERENCE IS TAKEN WITH THE DM AND NOT THE STAGE
% dm_zstack_demo_260911 found peak correlation rising monotonically with
% defocus (0.345 -> 0.526): DM defocus is NOT optically identical to stage
% translation, it carries some spherical aberration, so the PSF differs
% slightly with coefficient. That is harmless for holding a plane but it is
% exactly the mismatch a registration reference must not have. Block B shares
% its optical path with block C, so the estimator never pays for it. Block A
% exists only to put a micron label on the coefficient axis.
%
% WHY ALL FRAMES ARE KEPT
% The reference stack is acquired while the animal breathes, so each plane is
% the true structure convolved with the breathing displacement. The obvious
% worry is that this biases every later dz estimate toward zero.
%
% MEASURED 2026-09-13 (test_axial_estimator_260913, tests 6 and 6b): IT DOES
% NOT. A reference blurred over 4 um still recovers slope 1.004, and a
% reference averaged from only 2 frames per plane still recovers 0.995.
% Smoothing is a deterministic linear operation and least squares with a
% smoothed predictor stays consistent -- regression dilution needs NOISE in the
% predictor, and the fit averages over thousands of pixels, so reference noise
% washes out. Blur costs a constant offset and some sensitivity, not a slope.
%
% So gating is DEFENSIVE, not essential, and nothing here depends on it
% working. Frames are kept individually anyway (logAverageFactor = 1) because
% the gate can then be turned on, off or retuned offline without re-imaging,
% and because the per-frame lateral shifts are where the breath phase comes
% from. Keeping them costs disk and nothing else.
%
% ================================ SAFETY ================================
% Every DM command is baseline + defocus and is checked against the ALPAO +-1
% limit BEFORE being sent; a step that would clip is SKIPPED, never clamped --
% a clamped mirror is not producing the Zernike you asked for. The DM returns
% to your baseline and the stage to its starting position via onCleanup, on
% error and on Ctrl-C, along with every ScanImage setting this touches.
%
% PREREQUISITE: the system aberration correction must already be loaded and
% active, i.e. DMcommand_sys_aber_cor_930 exists in the base workspace. The
% baseline is READ FROM THAT VARIABLE, never from the DM -- DMctrl caches
% nothing and there is genuinely no read-back path on this hardware.
%
% Runqi Zhang / 2026-09-13.  NOT YET RUN -- rig was off when this was written.

%% ========================= USER SETTINGS =========================
SI_DIR       = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE     = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_DEFOC    = 3;

OUT_ROOT     = 'D:\Data\RZ\ao_axial_260913';   % session folders are made here

% WHICH CHANNEL IS THE GUIDE STAR
% Set this to the VESSEL / structural channel, not the GCaMP one.
%
% Measured 2026-09-13 (test_axial_estimator_260913, test 8): with 15% of cells
% firing a 3x transient, activity leaks 0.22 um into dz; at 40% it is 0.38 um.
% The global gain term does not absorb it (a transient is LOCAL) and the robust
% refit barely helps (0.19 um, and it gets WORSE at high activity) because a
% transient covers a whole cell and so is not a trimmable pixel minority.
% On a functional channel, breath-locked ACTIVITY is indistinguishable from
% breath-locked MOTION -- which is the exact thing being measured. A structural
% channel removes the mechanism instead of fighting it.
%
% SI logs every channel in channelSave into the same TIFF, so select BOTH the
% vessel and the GCaMP channel for Save. The analysis picks one by number, and
% running it on each is itself the experiment: if vessel-guided and
% GCaMP-guided dz agree, activity is not contaminating. If they diverge, the
% vessel answer is the trustworthy one.
CHANNEL      = 3;      % vessel / structural channel -- CHECK THIS FOR YOUR PREP

% ---- which blocks to run ------------------------------------------------
DO_STAGE_STACK = true;    % A  ruler
DO_DM_STACK    = true;    % B  reference   <- the one the estimator needs
DO_TIMESERIES  = true;    % C  measurement
DO_INJECTION   = true;    % D  positive control

% ---- block A: stage stack (the ruler) -----------------------------------
% z is 0.3125 um/count, so 1 um is 3.2 counts (-6.25% off) while 1.25 um is
% EXACTLY 4 counts. Use exact steps; the axial PSF is several um so 1.25 is
% already far denser than it needs to be.
STAGE_ZMIN_UM  = -10;
STAGE_ZMAX_UM  =  10;
STAGE_ZSTEP_UM = 1.25;
STAGE_SEC_PER_PLANE = 3.3;   % same breath-cycle argument as block B

% ---- block B: DM reference stack ----------------------------------------
% Same nominal range, in um of focus. Converted to coefficient by GAIN below.
% Keep NFRAMES generous: phase gating throws most of them away, and ~10
% surviving frames per plane is the floor for a usable reference.
% Time per plane is set by BREATH CYCLES, not SNR. Each plane must average over
% the whole cycle so its effective z is the MEAN position; too short and each
% plane lands on whatever phase it caught, offsetting its z by a fraction of
% the excursion -- differently for every plane, which is jitter in the ruler
% itself. (Distinct from reference NOISE, measured 2026-09-13 and found not to
% matter: 2 frames/plane still recovered slope 0.995.)
%
% IN SECONDS, NOT FRAMES: at 30 Hz, 100 frames is 3.3 s (~10 cycles); at 60 Hz
% the same 100 frames is 1.67 s and only ~5. Frame counts are not portable
% across frame rates. The count is derived from the measured fps below.
DM_ZMIN_UM   = -10;
DM_ZMAX_UM   =  10;
DM_ZSTEP_UM  = 1.25;
DM_SEC_PER_PLANE = 3.3;

% ---- block C: time series -----------------------------------------------
TS_SECONDS   = 120;       % at ~30 Hz this is ~3600 frames, ~350 breath cycles

% ---- block D: injected defocus staircase --------------------------------
% Known steps, in um of focus, held for INJ_HOLD_S each, commanded while the
% acquisition runs. Small: the point is to sit inside the linear zone and on
% top of the breathing, not to sweep the range.
INJ_STEPS_UM = [0 -2 0 +2 0 -4 0 +4 0];
INJ_HOLD_S   = 6;

% ---- calibration ---------------------------------------------------------
% um of focus per um of Zernike coefficient (Z2C row 3). Measured 2026-09-11:
% calibration said 6.656, four independent stack demos implied 6.76/6.87/7.66/
% 6.95. Quote ~7 +-5%. Block A is what re-measures this in TISSUE, so this
% number is only used to CHOOSE the steps, never to interpret them -- the
% analysis reads the true um axis back out of the A-vs-B registration.
GAIN_UM_PER_COEFF = 6.656;

AXIAL_PROBE_UM = 5;    % preflight: how far to step to test for axial blindness

SETTLE_S    = 0.3;
RES_ARM_SEC = 1.0;
CMD_LIMIT   = 0.95;
ENC_UM_PER_COUNT = -0.09987;   % signed: encoder counts DOWN as SI z goes up

%% ========================= CONNECT =========================
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;

hLSC = [];
try, hLSC = hSI.hMotors.hMotor(1).hLSC; catch, end

S = load(Z2C_FILE);
fn = fieldnames(S); Z2C = S.(fn{1});
assert(size(Z2C,2) == 97, 'expected 97 actuators, got %d', size(Z2C,2));

baseline = get_baseline(BASELINE_VAR);
p0 = hSI.hMotors.motorPosition;
z0 = p0(3);

STAMP   = datestr(now,'yymmdd_HHMM');
sessDir = fullfile(OUT_ROOT, sprintf('axial_%s', STAMP));
if ~exist(sessDir,'dir'), mkdir(sessDir); end

siState = struct( ...
    'loggingEnable',  hSI.hChannels.loggingEnable, ...
    'framesPerSlice', hSI.hStackManager.framesPerSlice, ...
    'numSlices',      hSI.hStackManager.numSlices, ...
    'rollingAvg',     hSI.hDisplay.displayRollingAverageFactor, ...
    'keepResOn',      hSI.hScan2D.keepResonantScannerOn, ...
    'logFilePath',    hSI.hScan2D.logFilePath, ...
    'logFileStem',    hSI.hScan2D.logFileStem, ...
    'logFileCounter', hSI.hScan2D.logFileCounter, ...
    'logAvgFactor',   hSI.hScan2D.logAverageFactor);
restore = onCleanup(@() restore_all(hDM, baseline, hSI, siState, p0)); %#ok<NASGU>

% RAW FRAMES, not averages -- offline phase gating needs every frame.
hSI.hScan2D.logAverageFactor = 1;
hSI.hScan2D.logFilePath      = sessDir;
hSI.hChannels.loggingEnable  = true;
hSI.hStackManager.numSlices  = 1;
hSI.hDisplay.displayRollingAverageFactor = 1;

%% ===================== PREFLIGHT: fail here, not in 6 minutes =====================
chSave = hSI.hChannels.channelSave;
fprintf('\nchannels saved: %s | registering on channel %d\n', mat2str(chSave(:).'), CHANNEL);
if ~ismember(CHANNEL, chSave)
    error('dm_axial:channelNotSaved', ...
        ['channel %d is not in channelSave (%s).\n' ...
         'The whole session would be logged without the channel the analysis ' ...
         'needs.\nSelect it for Save in ScanImage, or change CHANNEL.'], ...
        CHANNEL, mat2str(chSave(:).'));
end
if isempty(hSI.hChannels.channelDisplay)
    warning('dm_axial:noDisplay', ...
        ['no channel selected for Display. startGrab may refuse to start ' ...
         '(auto_acq_260910 hit exactly this on the rig).']);
end

% AXIAL BLINDNESS CHECK. A vessel running straight down z looks identical at
% every plane, so it carries NO axial information -- measured axial SNR 0.000,
% recovered slope 0.000, for a field of pure penetrating vessels
% (test_axial_estimator_260913, test 9). In-plane vessels gave SNR 4.1.
% There is no way to recover from this in analysis, so it is checked HERE:
% two frames a few um apart must actually look different.
fprintf('checking the field carries axial information...\n');
hSI.hChannels.loggingEnable = false;
arm_scanner(hSI, RES_ARM_SEC, 'the preflight probe');
I1 = quick_avg(hSI, 10, CHANNEL);
goto_z(hSI, p0, z0 + AXIAL_PROBE_UM);  pause(0.3);
I2 = quick_avg(hSI, 10, CHANNEL);
goto_z(hSI, p0, z0);  pause(0.3);
hSI.hChannels.loggingEnable = true;
a = double(I1(:)) - mean(double(I1(:)));
b = double(I2(:)) - mean(double(I2(:)));
rAx = (a.'*b) / (norm(a)*norm(b) + eps);
fprintf('  correlation between planes %g um apart: r = %.3f\n', AXIAL_PROBE_UM, rAx);
if rAx > 0.97
    warning('dm_axial:axiallyBlind', ...
        ['the field barely changes over %g um (r = %.3f). This looks AXIALLY ' ...
         'BLIND -- typically a field of penetrating vessels, which are ' ...
         'invariant under z translation.\nMove to a field with vessels running ' ...
         'ACROSS the plane, or branch points, before spending the animal.'], ...
        AXIAL_PROBE_UM, rAx);
    fprintf('  >> continuing anyway in 10 s; Ctrl-C now to move the FOV\n');
    wait_s(10);
end

% Held on for the whole run, AND re-armed with its settling time at the start of
% every block -- see arm_scanner().
arm_scanner(hSI, RES_ARM_SEC, 'the run');

try
    fps = hSI.hRoiManager.scanFrameRate;
catch
    fps = NaN;
end
fprintf('\nsession %s\n', sessDir);
fprintf('baseline %s, max|cmd| %.3f | frame rate %.2f Hz | channel %d\n', ...
        BASELINE_VAR, max(abs(baseline)), fps, CHANNEL);
fprintf('headroom on top of the baseline: coefficient +-%.2f  (%.1f um at %.2f um/coeff)\n', ...
        (CMD_LIMIT - max(abs(baseline))) / max(abs(Z2C(ROW_DEFOC,:))), ...
        (CMD_LIMIT - max(abs(baseline))) / max(abs(Z2C(ROW_DEFOC,:))) * GAIN_UM_PER_COEFF, ...
        GAIN_UM_PER_COEFF);

meta = struct();
meta.stamp = STAMP;  meta.sessDir = sessDir;  meta.channel = CHANNEL;
meta.baselineVar = BASELINE_VAR;  meta.baseline = baseline;
meta.rowDefoc = ROW_DEFOC;  meta.Z2CFile = Z2C_FILE;
meta.gainUmPerCoeff = GAIN_UM_PER_COEFF;
meta.encUmPerCount = ENC_UM_PER_COUNT;
meta.fpsHeader = fps;  meta.z0 = z0;  meta.p0 = p0;

%% ========================= A. STAGE STACK (ruler) =========================
if DO_STAGE_STACK
    fprintf('\n===== A. STAGE stack: motor moves, DM held at baseline =====\n');
    arm_scanner(hSI, RES_ARM_SEC, 'block A');
    hDM.senddata(baseline); pause(SETTLE_S);

    zList = STAGE_ZMIN_UM : STAGE_ZSTEP_UM : STAGE_ZMAX_UM;
    nZ    = numel(zList);
    STAGE_NFRAMES = max(4, round(STAGE_SEC_PER_PLANE * fps));
    fprintf('  %d planes, %.2f s each = %d frames\n', nZ, STAGE_SEC_PER_PLANE, STAGE_NFRAMES);

    % Run-up below the first plane so the sweep is monotonic and no plane
    % absorbs the 0.205 um reversal backlash.
    goto_z(hSI, p0, z0 + zList(1) - 5*STAGE_ZSTEP_UM);
    wait_s(0.3);

    zEnc  = nan(1,nZ);
    files = cell(1,nZ);
    for k = 1:nZ
        goto_z(hSI, p0, z0 + zList(k));
        pause(0.2);
        zEnc(k) = read_z_counts(hLSC);
        files{k} = log_grab(hSI, STAGE_NFRAMES, sprintf('A_stage_z%+07.2f', zList(k)));
        if mod(k,5) == 0 || k == nZ, fprintf('  plane %d/%d\n', k, nZ); end
    end
    goto_z(hSI, p0, z0); wait_s(0.3);

    if all(isfinite(zEnc)) && range(zEnc) > 0
        zMeas = (zEnc - zEnc(1)) * ENC_UM_PER_COUNT + zList(1);
        fprintf('  encoder: %.4f um/step measured vs %g commanded\n', ...
                (zMeas(end)-zMeas(1))/(nZ-1), STAGE_ZSTEP_UM);
    else
        zMeas = zList;
        warning('dm_axial:noEnc','z encoder unreadable; recording COMMANDED z');
    end

    meta.A = struct('zNom',zList,'zEnc',zEnc,'zMeas',zMeas, ...
                    'nFrames',STAGE_NFRAMES,'files',{files});
    save(fullfile(sessDir,'meta.mat'),'-struct','meta');
end

%% ========================= B. DM STACK (reference) =========================
if DO_DM_STACK
    fprintf('\n===== B. DM stack: STAGE PARKED, defocus only. This is the reference. =====\n');
    arm_scanner(hSI, RES_ARM_SEC, 'block B (the reference stack)');
    zParked = hSI.hMotors.motorPosition(3);
    fprintf('  stage parked at z = %.3f um and left there\n', zParked);

    zList = DM_ZMIN_UM : DM_ZSTEP_UM : DM_ZMAX_UM;
    nZ    = numel(zList);
    DM_NFRAMES = max(4, round(DM_SEC_PER_PLANE * fps));
    fprintf('  %d planes, %.2f s each = %d frames (~%.1f breath cycles at 3 Hz)\n', ...
            nZ, DM_SEC_PER_PLANE, DM_NFRAMES, DM_SEC_PER_PLANE*3);
    cUsed = nan(1,nZ);  got = false(1,nZ);  files = cell(1,nZ);

    for k = 1:nZ
        c   = zList(k) / GAIN_UM_PER_COEFF;
        cmd = baseline(:).' + c * Z2C(ROW_DEFOC,:);
        if max(abs(cmd)) > CMD_LIMIT
            fprintf('  z %+6.2f um (c %+5.2f): SKIPPED, would reach |cmd| %.3f\n', ...
                    zList(k), c, max(abs(cmd)));
            continue
        end
        hDM.senddata(cmd);
        pause(SETTLE_S);
        files{k} = log_grab(hSI, DM_NFRAMES, sprintf('B_dm_z%+07.2f', zList(k)));
        cUsed(k) = c;  got(k) = true;
        if mod(k,5) == 0 || k == nZ, fprintf('  plane %d/%d (c %+5.2f)\n', k, nZ, c); end
    end
    hDM.senddata(baseline); pause(SETTLE_S);

    zAfter = hSI.hMotors.motorPosition(3);
    fprintf('  stage z before %.3f, after %.3f (drift %+.3f um)\n', ...
            zParked, zAfter, zAfter - zParked);
    fprintf('  planes the DM could not reach: %d of %d\n', nnz(~got), nZ);

    meta.B = struct('zNom',zList,'coeff',cUsed,'got',got, ...
                    'nFrames',DM_NFRAMES,'files',{files},'zParked',zParked);
    save(fullfile(sessDir,'meta.mat'),'-struct','meta');
end

%% ========================= C. TIME SERIES =========================
if DO_TIMESERIES
    fprintf('\n===== C. TIME SERIES: stage parked, DM at baseline, animal breathing =====\n');
    arm_scanner(hSI, RES_ARM_SEC, 'block C (the measurement)');
    hDM.senddata(baseline); pause(SETTLE_S);

    if isfinite(fps) && fps > 0
        nFr = max(1, round(TS_SECONDS * fps));
    else
        nFr = round(TS_SECONDS * 30);
        warning('dm_axial:noFps','scanFrameRate unreadable; assuming 30 Hz for the COUNT only');
    end
    fprintf('  %d frames (~%.0f s). Do not touch the rig.\n', nFr, TS_SECONDS);
    t0 = now;
    f = log_grab(hSI, nFr, 'C_timeseries');
    meta.C = struct('nFrames',nFr,'file',f,'tStart',t0,'requestedSec',TS_SECONDS);
    save(fullfile(sessDir,'meta.mat'),'-struct','meta');
end

%% ========================= D. INJECTED STAIRCASE =========================
if DO_INJECTION
    fprintf('\n===== D. INJECTION: known defocus staircase during acquisition =====\n');
    fprintf('  This is the positive control. If the estimator cannot recover a\n');
    fprintf('  step it was TOLD about, it cannot be trusted on breathing.\n');

    % Pre-compute and pre-check every command BEFORE anything is sent, so a
    % clipping step is refused up front rather than mid-acquisition.
    nStep = numel(INJ_STEPS_UM);
    cmds  = zeros(nStep, 97);  cList = zeros(1,nStep);  okStep = false(1,nStep);
    for k = 1:nStep
        cList(k) = INJ_STEPS_UM(k) / GAIN_UM_PER_COEFF;
        cmds(k,:) = baseline(:).' + cList(k) * Z2C(ROW_DEFOC,:);
        okStep(k) = max(abs(cmds(k,:))) <= CMD_LIMIT;
    end
    if ~all(okStep)
        warning('dm_axial:injClip','%d injection step(s) would clip and are REPLACED BY 0 um', nnz(~okStep));
        for k = find(~okStep)
            cmds(k,:) = baseline(:).';  cList(k) = 0;  INJ_STEPS_UM(k) = 0;
        end
    end

    if isfinite(fps) && fps > 0
        nFr = max(1, round(nStep * INJ_HOLD_S * fps));
    else
        nFr = round(nStep * INJ_HOLD_S * 30);
    end

    arm_scanner(hSI, RES_ARM_SEC, 'block D');
    hDM.senddata(cmds(1,:)); pause(SETTLE_S);
    fprintf('  %d steps x %.1f s = %d frames\n', nStep, INJ_HOLD_S, nFr);

    % Start the acquisition, then walk the staircase while it runs. MATLAB
    % pause() jitter is milliseconds against a 6 s hold, so the step TIMES are
    % recorded rather than assumed -- the analysis aligns on the recorded
    % times, it does not reconstruct them from INJ_HOLD_S.
    hSI.hStackManager.framesPerSlice = nFr;
    hSI.hScan2D.logFileStem = sprintf('D_inject_%s', STAMP);
    hSI.hScan2D.logFileCounter = 1;
    tAcq = tic;  tWall0 = now;
    hSI.startGrab();

    tSend = nan(1,nStep);
    for k = 1:nStep
        hDM.senddata(cmds(k,:));
        tSend(k) = toc(tAcq);
        fprintf('    t=%6.2f s  step %+5.1f um (c %+5.2f)\n', tSend(k), INJ_STEPS_UM(k), cList(k));
        tEnd = k * INJ_HOLD_S;
        while toc(tAcq) < tEnd && ~strcmp(hSI.acqState,'idle')
            pause(0.02); drawnow limitrate;
        end
    end
    while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
    hDM.senddata(baseline); pause(SETTLE_S);

    meta.D = struct('stepsUm',INJ_STEPS_UM,'coeff',cList,'holdS',INJ_HOLD_S, ...
                    'tSendSec',tSend,'tWall0',tWall0,'nFrames',nFr, ...
                    'file',fullfile(sessDir, sprintf('D_inject_%s_00001.tif', STAMP)));
    save(fullfile(sessDir,'meta.mat'),'-struct','meta');
end

%% ========================= DONE =========================
save(fullfile(sessDir,'meta.mat'),'-struct','meta');
fprintf('\n=========================================================\n');
fprintf('done. session written to\n  %s\n', sessDir);
fprintf('next:  dm_axial_analyze_260913(''%s'')\n', sessDir);
fprintf('=========================================================\n');
end

%% ========================== HELPERS ==========================
function f = log_grab(hSI, nFrames, stem)
% One logged grab of nFrames RAW frames. Returns the path SI will have written.
%
% THE IDLE GUARD IS NOT OPTIONAL. startGrab() can return with acqState still
% 'idle' -- most often because no channel is selected for Display or Save in
% ScanImage. The obvious wait loop `while ~strcmp(acqState,'idle')` then exits
% IMMEDIATELY, this function returns a filename that was never written, and the
% whole block "completes" in seconds having acquired nothing. auto_acq_260910
% hit exactly this on the rig, which is why it errors here instead.
hSI.hStackManager.framesPerSlice = nFrames;
hSI.hStackManager.numSlices      = 1;
hSI.hScan2D.logFileStem          = stem;
hSI.hScan2D.logFileCounter       = 1;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_axial:grabDidNotStart', ...
        ['startGrab returned with acqState still ''idle'' for "%s".\n' ...
         'The grab never started -- most likely no channel is selected for ' ...
         'Display or Save in ScanImage.\nFix that and re-run.'], stem);
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
f = fullfile(hSI.hScan2D.logFilePath, sprintf('%s_00001.tif', stem));
if ~exist(f,'file')
    warning('dm_axial:noFile','expected "%s" but it is not on disk', f);
end
end

function I = quick_avg(hSI, n, ch)
% Unlogged averaged grab, for the preflight only. Reads the display buffer
% rather than a file, so it leaves nothing on disk.
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_axial:grabDidNotStart', ...
        ['startGrab returned with acqState still ''idle'' during preflight.\n' ...
         'Most likely no channel is selected for Display in ScanImage.']);
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
f = hSI.hDisplay.lastFrame;
chans = hSI.hDisplay.lastFrameChannels;
i = find(chans == ch, 1);
if isempty(i)
    error('dm_axial:noChannel','channel %d not displayed (available %s)', ...
          ch, mat2str(chans));
end
I = single(f{i});
end

function b = get_baseline(varName)
% The baseline is the VARIABLE YOU LOADED, not a read-back. DMctrl caches
% nothing (its properties are only hDM/Gain/serialName) and AOdata's getters
% are broken in this install, so there is no other honest source. Falling back
% to zeros(1,97) would silently measure everything about a FLAT mirror and
% then flatten the DM on cleanup -- so this errors instead.
if ~evalin('base', sprintf('exist(''%s'',''var'')', varName))
    error('dm_axial:noBaseline', ...
        ['base workspace has no ''%s''. Load and activate the system ' ...
         'aberration correction before running this.'], varName);
end
b = double(evalin('base', varName));  b = b(:).';
assert(numel(b) == 97, '%s is %d long, expected 97', varName, numel(b));
end

function goto_z(hSI, p0, zTarget)
p = p0;  p(3) = zTarget;
hSI.hMotors.moveStartRelative(p);
hSI.hMotors.moveWaitForFinish([1 2 3]);
end

function c = read_z_counts(hLSC)
c = NaN;
try
    v = hLSC.getEncoderPosition();
    if numel(v) >= 3, c = v(3); end
catch
end
end

function wait_s(sec)
t0 = tic;  while toc(t0) < sec, pause(0.05); drawnow limitrate; end
end

function arm_scanner(hSI, sec, what)
%ARM_SCANNER  Resonant scanner on and settled before ANY acquisition.
%
% STANDARD PRACTICE ON THIS RIG (RZ, 2026-09-13): the resonant scanner is
% turned on and given ~1 s to reach stable amplitude before every acquisition
% -- the reference stack included, not just the "real" one. A scanner still
% ramping gives a field whose scale and linearity are drifting, which in this
% pipeline corrupts the reference volume itself and therefore every dz estimate
% measured against it. There is no way to detect that after the fact.
%
% Called at the start of every block. keepResonantScannerOn is also held true
% across the whole run, so in practice this is belt and braces -- but a block
% that is skipped, retried or run on its own still gets its arming time.
hSI.hScan2D.keepResonantScannerOn = true;
fprintf('  arming resonant scanner %.1f s before %s...\n', sec, what);
wait_s(sec);
end

function restore_all(hDM, baseline, hSI, s, p0)
try, hDM.senddata(baseline); fprintf('[cleanup] DM restored to baseline\n');
catch ME, warning('dm_axial:noDM','COULD NOT RESTORE THE DM: %s', ME.message); end
try
    hSI.hMotors.moveStartRelative(p0);
    hSI.hMotors.moveWaitForFinish([1 2 3]);
    fprintf('[cleanup] stage returned to z %.3f\n', p0(3));
catch ME, warning('dm_axial:noStage','could not return the stage: %s', ME.message); end
try
    hSI.hChannels.loggingEnable              = s.loggingEnable;
    hSI.hStackManager.framesPerSlice         = s.framesPerSlice;
    hSI.hStackManager.numSlices              = s.numSlices;
    hSI.hDisplay.displayRollingAverageFactor = s.rollingAvg;
    hSI.hScan2D.keepResonantScannerOn        = s.keepResOn;
    hSI.hScan2D.logFilePath                  = s.logFilePath;
    hSI.hScan2D.logFileStem                  = s.logFileStem;
    hSI.hScan2D.logFileCounter               = s.logFileCounter;
    hSI.hScan2D.logAverageFactor             = s.logAvgFactor;
    fprintf('[cleanup] ScanImage settings restored\n');
catch ME, warning('dm_axial:noSI','could not restore SI settings: %s', ME.message); end
end
