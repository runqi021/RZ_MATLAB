function dm_axial_closeloop_260913()
%DM_AXIAL_CLOSELOOP_260913  Closed-loop axial motion correction with the DM.
%
% Image in, DM out. Nothing else enters the loop:
%
%     frame -> compare to reference stack -> dz -> send -(dz) -> repeat
%
% No breath sensor, no camera, no phase model, no periodicity assumed. The
% reference stack is acquired once at the start and is the only thing the loop
% needs.
%
% THE EXPERIMENT, in one acquisition so nothing drifts between conditions:
%   OPEN (correction off)  -> CLOSED (correction on) -> OPEN again
% Every frame is logged, so the residual motion is measured OFFLINE against the
% same reference. That is the honest test: the loop's own dz goes to zero by
% construction when it works, so it cannot be its own evidence. What counts is
% whether the LOGGED FRAMES are more stable in the CLOSED window, which
% dm_axial_analyze_260913 measures without knowing which window was which.
%
% VELOCITY EXTRAPOLATION
% The error is pure delay: an estimate is the mean over its frame and cannot be
% applied until that frame has finished, so it lands one frame period late. At
% 30 Hz that is ~50 ms, which against 3 Hz breathing leaves most of the motion
% and above ~6 Hz makes it worse. So extrapolate with the velocity already in
% hand:   send -( e[k] + gx*(e[k]-e[k-1]) ).
% Two estimates and a subtraction. NOT a breath model -- no periodicity, keeps
% working through sniffs and apnea.
%
% ============================== SAFETY ==============================
% DRY_RUN IS TRUE BY DEFAULT. The first run computes and logs everything and
% sends NOTHING, so you can watch the numbers before the mirror moves.
%
% Three independent clamps, all applied before any senddata:
%   DEFOCUS_MIN/MAX_UM  absolute excursion from the baseline (asymmetric,
%                   because the DM's authority is -14 to +24, not +-N)
%   MAX_STEP_UM     change per update (rate limit -- stops a bad estimate
%                   throwing the focus across the range in one frame)
%   CMD_LIMIT       the ALPAO +-1 actuator limit, checked on the full command
% A command that would breach the actuator limit is SKIPPED, not clamped.
% onCleanup restores the DM baseline and every SI setting on error and Ctrl-C.
%
% SIGN: dz is the DISPLACEMENT of the sample; the correction is -dz. Verified
% both directions in test_axial_estimator_260913. Backwards, the loop drives
% the motion instead of cancelling it and runs away -- which is what MAX_STEP_UM
% and the DEFOCUS_MIN/MAX_UM clamps exist to survive.
%
% Runqi Zhang / 2026-09-13.  NOT YET RUN.

%% ========================= USER SETTINGS =========================
DRY_RUN      = true;    % <<< set false only after watching a dry run

SI_DIR       = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE     = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_DEFOC    = 3;
OUT_ROOT     = 'D:\Data\RZ\ao_axial_260913';

% VESSEL / structural channel, not the GCaMP one: on a functional channel a
% calcium transient is a LOCAL brightness change that the global gain term does
% not absorb, so breath-locked ACTIVITY becomes indistinguishable from
% breath-locked MOTION -- which is the thing being measured.
%
% For THIS test, ch3 alone is all that needs to be saved. Saving the GCaMP
% channel too is only worth it when the question is "is activity contaminating
% the estimate?", answered by running the analysis on each and comparing. That
% is not this run, and dropping it halves the file (~1.8 GB instead of ~3.7 GB
% for 2 min at 60 Hz) and the data rate with it.
CHANNEL      = 3;

% ---- reference stack (acquired inline, stage parked, DM stepping) -------
%
% GEOMETRY IS INHERITED, NOT SET HERE. The reference is acquired at whatever
% zoom and linesPerFrame ScanImage is currently on, which is the same geometry
% the loop then runs at -- consistent by construction. That matters: on a
% resonant rig the fast axis amplitude IS the zoom and the slow axis sets the
% height, so 256 lines is a genuinely SMALLER FIELD, not a coarser one. A
% 512-line reference cannot register 256-line frames. Change lines for a faster
% loop and just re-run; the reference costs under a minute.
%
% NO EXTRA XY MARGIN IS NEEDED. Only z is corrected. Lateral motion is a few px
% out of 512 and is removed by an integer pre-shift before the axial fit; the
% circshift wrap corrupts a band a few px wide, under 1% of the frame, and the
% robust trim absorbs it.
%
% Z EXTENT: take the WHOLE correctable range, -14 to +24 um. The DM cannot
% correct outside it, so a reference that spans it can never be the thing that
% runs out -- and the excursion is not yet known, so a stack sized to a guess
% is a guess you have to re-take. It costs ~45 s once.
%
% NOTE the range is the DM's, and it is ASYMMETRIC because the system
% aberration correction already uses max|cmd| 0.673 of the +-1 stroke. Centre
% your FOV at the focus you want, then this spans everything reachable from
% there. Any plane that would breach CMD_LIMIT is SKIPPED and reported, so if a
% different system correction is loaded (there is a 10x variant) the stack
% simply self-trims to whatever that baseline allows -- no edit needed.
%
% Step 1 um: the DM is continuous so there is no quantisation to respect here
% (unlike the STAGE, where 1 um is 3.2 counts and 1.25 um is exactly 4), and
% 1 um is well inside the axial PSF so the gradient is not aliased.
%
% Once the excursion IS known, narrowing to a few um either side is faster and
% bleaches less -- and the far planes carry a caveat: DM defocus picks up
% spherical aberration that grows with coefficient (peak correlation rose
% 0.345 -> 0.526 across the range on 2026-09-11), so planes far from zero are
% not optically identical to the same plane reached by moving the stage.
%
% FRAMES PER PLANE: 60, and the reason is BREATH CYCLES, not SNR.
% At ~30 Hz, 60 frames is 2 s -- about 6 breath cycles -- so each plane averages
% over the whole cycle and its effective z is the MEAN position. At 20 frames
% (0.67 s, ~2 cycles) each plane instead lands on whatever phase it happened to
% catch, so its effective z is offset by a fraction of the excursion and the
% offset DIFFERS PLANE TO PLANE. That is jitter in the ruler itself.
%
% Note this is a different failure from reference NOISE, which was measured on
% 2026-09-13 and found not to matter (2 frames/plane still recovered slope
% 0.995). More pixels fix noise; only more CYCLES fix phase sampling.
%
% SPECIFIED IN SECONDS, NOT FRAMES, because what buys cycles is TIME. At 30 Hz
% 100 frames is 3.3 s (~10 cycles); at 60 Hz the same 100 frames is 1.67 s and
% only ~5 cycles. Frame count is not portable across frame rates -- seconds
% are. The count is computed below from the ACTUAL scanFrameRate.
%
% 3.3 s is ~10 cycles at 3 Hz. Residual phase bias falls as 1/sqrt(cycles), so
% going below ~2 s (6 cycles) starts to cost; above ~5 s the return is small.
% Worth it on a VESSEL channel, which does not bleach the way GCaMP would;
% drop to 2 s if the dye is dim or you are short of time.
% *** MEASURED 2026-09-13: THE FULL-AUTHORITY STACK IS THE WRONG CHOICE. ***
% Same mirror, rate, read path and estimator, 5 um of motion at 0.5 Hz:
%     reference -10..+10 -> estimator gain 0.938
%     reference -14..+24 -> estimator gain 0.653
% Spanning the DM's whole range costs 35% of the signal, because the far
% positive planes hold no usable structure, the coarse tier sometimes locks
% onto one, and those outliers eat the coherent amplitude.
%
% Too narrow is also wrong: 10 um of motion in a +-10 um reference measured
% 0.810, the estimate saturating at the edges.
%
% RULE: wider than the motion, no wider than the usable signal. For +-5 um of
% breathing, +-12 satisfies both. Widen only if the measured excursion grows,
% and re-check the gain with dm_dynamic_gain_260913 if you do.
REF_ZMIN_UM  = -12;
REF_ZMAX_UM  =  12;
REF_ZSTEP_UM =  1;
REF_SEC_PER_PLANE = 3.3;   % 25 planes -> ~1.5 min at ANY frame rate

% Reuse the stack from an earlier session instead of re-shooting it. '' = take
% a fresh one. Useful while iterating, but ONLY valid if nothing has moved:
% same FOV, same focus, same geometry, same animal, no drift. The script checks
% the frame SIZE matches and refuses otherwise -- it cannot check that the
% sample is still where it was, so that part is on you. Take a fresh stack for
% any run whose numbers you intend to believe.
REF_REUSE_FROM = '';       % e.g. 'D:\Data\RZ\ao_axial_260913\loop_260913_1846'
REUSE_Z_TOL_UM = 1.0;      % refuse a reused reference if the stage has moved
                           % more than this since it was taken -- beyond that
                           % its zero is no longer your zero

% ---- frame rate ----------------------------------------------------------
% [] = inherit whatever ScanImage is set to. A number = set linesPerFrame to
% reach that rate, keeping pixelsPerLine (so the RESOLUTION is unchanged and
% only the FOV HEIGHT shrinks), then restore on exit.
%
% Why bother: the loop's delay is one frame period. At 30 Hz that is ~50 ms,
% which leaves ~91% of 3 Hz motion and AMPLIFIES above ~6 Hz. At 60 Hz the
% break-even moves to ~13 Hz and the residual drops sharply. Lines are the
% currency you buy loop rate with.
TARGET_FPS   = 60;              % [] to inherit

% ---- bidirectional scan phase -------------------------------------------
% -0.7 for this condition (RZ, 2026-09-13). [] = leave whatever is dialled in.
%
% This matters more here than for ordinary imaging. A wrong scan phase shifts
% alternate lines against each other, and that comb is a STRUCTURED, fixed
% pattern sitting on every frame. The estimator would happily fit part of it.
% Worse, if the phase differs between the reference stack and the live frames,
% the difference looks like a real image change and lands in dz.
%
% Reference and live frames come from the same run, so they are consistent
% either way -- but pinning it makes the session reproducible and puts the
% value in the metadata instead of in someone's memory.
%
% *** UNITS: MICROSECONDS. hScan2D.linePhase is in SECONDS. ***
% Measured on the rig 2026-09-13: the working value read back as -7e-07, i.e.
% -0.7 us. Writing -0.7 directly sets -0.7 SECONDS, which ScanImage accepts
% and then floods the console with "Phase is too negative. Adjust the physical
% scan phase on the resonant driver board" -- once per grab -- while imaging
% with a completely wrong phase. A read-back that merely confirms the write
% landed does NOT catch this; the value has to be sanity-checked for
% plausibility as well. (Same trap as positionDeviceUnits being in metres.)
SCAN_PHASE_US = -0.7;           % MICROSECONDS. [] to inherit.

% ---- the run -------------------------------------------------------------
SEG_SECONDS  = [30 60 30];      % OPEN, CLOSED, OPEN
LOOP_GAIN    = 0.7;             % fraction of the correction applied per update
EXTRAP_GAIN  = 1.0;             % gx above; 0 disables extrapolation

% ---- clamps --------------------------------------------------------------
% Sized for the motion actually being corrected, not for a guess. RZ reports
% typical excursions around +-5 um, so:
%   the defocus clamp must be comfortably MORE than the excursion or the loop
%   saturates exactly when the motion is largest. +-5 um was equal to the
%   motion and therefore useless.
%   MAX_STEP_UM is a velocity limit: um per frame x frame rate. 5 um p2p at
%   3 Hz peaks at pi*3*5 = 47 um/s; at 60 Hz, 1 um/frame allows only 60 um/s,
%   so the rate limit would clip the fastest part of every breath. 2 um/frame
%   gives 120 um/s, ~2.5x headroom, while still bounding a bad estimate to one
%   plane spacing per frame.
%
% -13 / +23 (RZ, 2026-09-13): the DM's measured authority is -14 to +24, so this
% is the full range with a micron of margin at each end. Deliberately
% ASYMMETRIC -- a symmetric +-13 would have discarded 10 um of perfectly good
% positive authority, and which way the brain actually moves is not yet known,
% so there is no reason to give any of it away.
% (Re-centring the operating point -- a static defocus bias on the baseline plus
% a compensating stage park -- would buy roughly +-19 either way instead.)
DEFOCUS_MIN_UM = -13.0;
DEFOCUS_MAX_UM =  23.0;
MAX_STEP_UM    = 2.0;
CMD_LIMIT      = 0.95;

GAIN_UM_PER_COEFF = 6.656;      % bench value; only sets the command scale
SETTLE_S    = 0.3;
RES_ARM_SEC = 1.0;

%% ========================= CONNECT =========================
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;

S = load(Z2C_FILE);
fn = fieldnames(S); Z2C = S.(fn{1});
assert(size(Z2C,2) == 97, 'expected 97 actuators, got %d', size(Z2C,2));
zRow = Z2C(ROW_DEFOC,:);

baseline = get_baseline(BASELINE_VAR);
p0 = hSI.hMotors.motorPosition;

STAMP   = datestr(now,'yymmdd_HHMM');
sessDir = fullfile(OUT_ROOT, sprintf('loop_%s', STAMP));
if ~exist(sessDir,'dir'), mkdir(sessDir); end

% EVERY restorable value must be captured BEFORE the onCleanup below, because
% onCleanup takes a snapshot of siState BY VALUE. Fields added to siState after
% the closure exists are invisible to it. That is not hypothetical: on
% 2026-09-13 linePhase was stored after this line, the restore silently never
% ran (no error, just a missing line in the cleanup log), and the rig was left
% at -0.7 SECONDS of line phase.
siState = struct( ...
    'linePhase',            hSI.hScan2D.linePhase, ...
    'linesPerFrame',        hSI.hRoiManager.linesPerFrame, ...
    'forceSquarePixelation',hSI.hRoiManager.forceSquarePixelation, ...
    'forceSquarePixels',    hSI.hRoiManager.forceSquarePixels, ...
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

hSI.hScan2D.logAverageFactor = 1;
hSI.hScan2D.logFilePath      = sessDir;
hSI.hStackManager.numSlices  = 1;
hSI.hDisplay.displayRollingAverageFactor = 1;

% BOTH lists must contain the channel, for two different reasons, and checking
% only one is how the 2026-09-13 dry run "serviced" 7202 frames while doing
% nothing at all:
%   channelSave    -> it reaches the TIFF, which is what the offline evaluation
%                     re-measures dz on;
%   channelDisplay -> it reaches hDisplay.lastFrame, which is the ONLY place the
%                     live loop can read a frame from.
% Saved-but-not-displayed meant get_disp_frame() returned empty every single
% frame, the loop hit `continue` before the estimate, the log and the progress
% print, and the run looked like it had worked.
chSave = hSI.hChannels.channelSave;
chDisp = hSI.hChannels.channelDisplay;
assert(ismember(CHANNEL, chSave), ...
    'channel %d is not in channelSave (%s) -- it would never reach the TIFF', ...
    CHANNEL, mat2str(chSave(:).'));
assert(ismember(CHANNEL, chDisp), ...
    ['channel %d is not in channelDisplay (%s).\nThe live loop reads frames ' ...
     'from the DISPLAY buffer, so it would see nothing at all.\nTick channel ' ...
     '%d for Display in ScanImage.'], ...
    CHANNEL, mat2str(chDisp(:).'), CHANNEL);

%% ===================== FRAME RATE / LINES =====================
% The loop's delay is one frame period, so lines are how you buy loop rate.
%
% *** THE TRAP: forceSquarePixelation SILENTLY UNDOES THIS. ***
% It defaults TRUE, and RoiManager.m:934 resets linesPerFrame back to
% pixelsPerLine whenever they differ. Set 256 lines with it on and you get 512
% back with no error and no warning -- you would simply still be at 30 Hz and
% never know. So it must be turned OFF first, and the result must be READ BACK
% and verified, not assumed.
%
% forceSquarePixels is a DIFFERENT property and must stay ON: it keeps pixels
% square by shrinking scanAngleMultiplierSlow with the line count, so half the
% lines is half the FOV HEIGHT at unchanged um/px. That is the trade we want --
% field of view, not resolution.
if ~isempty(TARGET_FPS)
    lines0 = hSI.hRoiManager.linesPerFrame;
    fps0   = hSI.hRoiManager.scanFrameRate;
    lp     = hSI.hRoiManager.linePeriod;
    assert(isfinite(fps0) && fps0 > 0 && isfinite(lp) && lp > 0, ...
           'cannot read scanFrameRate / linePeriod');

    % Frame period = flyback + lines*linePeriod. Solve for the lines that hit
    % the target, measuring the flyback rather than assuming it is zero.
    flyback = 1/fps0 - lines0*lp;
    want    = round((1/TARGET_FPS - flyback) / lp);
    want    = 2*round(want/2);                 % even: resonant is bidirectional
    want    = max(32, min(lines0, want));      % never ask for MORE lines than set

    hSI.hRoiManager.forceSquarePixelation = false;   % or the next line is undone
    hSI.hRoiManager.forceSquarePixels     = true;    % shrink the FOV, not the pixel
    hSI.hRoiManager.linesPerFrame         = want;

    got = hSI.hRoiManager.linesPerFrame;             % READ BACK. Never assume.
    if got ~= want
        error('dm_loop:linesRefused', ...
            ['asked for %d lines, ScanImage kept %d. forceSquarePixelation is ' ...
             '%d.\nSomething else is constraining the ROI (mroiEnable?) -- set ' ...
             'the frame size by hand and use TARGET_FPS = [].'], ...
            want, got, hSI.hRoiManager.forceSquarePixelation);
    end
    fprintf('\nlines %d -> %d (px/line %d unchanged) : %.2f -> %.2f Hz\n', ...
            lines0, got, hSI.hRoiManager.pixelsPerLine, fps0, hSI.hRoiManager.scanFrameRate);
    fprintf('FOV HEIGHT IS NOW %.0f%% of what it was; pixel size unchanged.\n', ...
            100*got/lines0);
end

%% ===================== SCAN PHASE =====================
% Set whichever property this install actually exposes, then READ BACK. The
% name and units differ between SI versions (scanPhase in GUI units vs
% linePhase in seconds), so this detects rather than assumes -- and if the
% read-back does not match, it says so instead of quietly imaging with a comb.
if ~isempty(SCAN_PHASE_US)
    assert(isprop(hSI.hScan2D,'linePhase'), 'hScan2D has no linePhase');
    phasePrev  = siState.linePhase;              % SECONDS, captured pre-cleanup
    wantPhaseS = SCAN_PHASE_US * 1e-6;

    % PLAUSIBILITY GATE, not just a read-back. A resonant line is ~63 us, so a
    % line phase outside +-20 us is not a phase, it is a units error -- and
    % ScanImage will happily accept it and image with it.
    assert(abs(wantPhaseS) < 20e-6, ...
        ['SCAN_PHASE_US = %g would set linePhase to %g s. linePhase is in ' ...
         'SECONDS and a resonant line is ~63 us, so this is a units error. ' ...
         'Give microseconds.'], SCAN_PHASE_US, wantPhaseS);

    hSI.hScan2D.linePhase = wantPhaseS;
    got = hSI.hScan2D.linePhase;
    fprintf('scan phase: %.3f us -> %.3f us (linePhase %g s)\n', ...
            phasePrev*1e6, got*1e6, got);
    if abs(got - wantPhaseS) > 1e-9
        warning('dm_loop:phaseRefused', ...
            'asked for %.3f us, ScanImage kept %.3f us', SCAN_PHASE_US, got*1e6);
    end
else
    fprintf('scan phase: inheriting %.3f us\n', hSI.hScan2D.linePhase*1e6);
end

fps = hSI.hRoiManager.scanFrameRate;
assert(isfinite(fps) && fps > 0, 'could not read scanFrameRate');
if ~isempty(TARGET_FPS) && abs(fps - TARGET_FPS)/TARGET_FPS > 0.1
    warning('dm_loop:fpsMiss', ...
        'asked for %.1f Hz, got %.2f Hz. Line period may not divide evenly.', ...
        TARGET_FPS, fps);
end

fprintf('\n=================== CLOSED-LOOP AXIAL CORRECTION ===================\n');
if DRY_RUN
    fprintf('*** DRY RUN -- everything is computed and logged, NOTHING is sent ***\n');
end
fprintf('session %s\n', sessDir);
fprintf('frame rate %.2f Hz -> one frame of delay is %.0f ms\n', fps, 1000/fps);
fprintf('channel %d | loop gain %.2f | extrapolation gain %.2f\n', CHANNEL, LOOP_GAIN, EXTRAP_GAIN);
fprintf('clamps: defocus %.1f to %+.1f um, |step| <= %.2f um/frame\n', ...
        DEFOCUS_MIN_UM, DEFOCUS_MAX_UM, MAX_STEP_UM);

% What this rate can and cannot do, stated before the run rather than after.
tau = 1/fps;
fBreak = 1/(6*tau);
fprintf('break-even frequency %.1f Hz: above this a PLAIN loop amplifies.\n', fBreak);
fprintf('  plain      leaves %.0f%% of 3 Hz motion\n', 100*2*sin(pi*3*tau));
fprintf('  extrapolated leaves ~%.0f%% (second order in f*tau)\n', ...
        100*min(1, 0.5*(2*pi*3*tau)^2));

%% ========================= 1. REFERENCE STACK =========================
fprintf('\n--- 1. reference stack (stage parked, DM stepping) ---\n');
hSI.hChannels.loggingEnable = false;

% *** THE REFERENCE DEFINES z = 0, SO IT MUST BE ANCHORED TO A REAL z. ***
% dz is measured against the reference stack, whose zero is wherever the stage
% was parked when that stack was taken. Without recording that park position,
% every new reference silently re-zeros the z coordinate and dz cannot be
% compared across sessions -- and a REUSED reference taken at a different park
% would shift the zero with no complaint at all.
% So: record the stage z with the stack, report absolute depth alongside dz,
% and refuse a reused reference if the stage has moved.
stageZnow = hSI.hMotors.motorPosition(3);

if ~isempty(REF_REUSE_FROM)
    p = fullfile(REF_REUSE_FROM, 'loop_meta.mat');
    assert(exist(p,'file') == 2, 'no loop_meta.mat in %s', REF_REUSE_FROM);
    Mold = load(p, 'refVol', 'zRef', 'linesPerFrame', 'pixelsPerLine', ...
                   'linePhaseS', 'refStageZ');
    refVol = Mold.refVol;  zRef = Mold.zRef;
    if isfield(Mold,'refStageZ') && ~isempty(Mold.refStageZ)
        refStageZ = Mold.refStageZ;
        dStage = stageZnow - refStageZ;
        fprintf('  reference was taken at stage z = %.2f um; stage is now %.2f um\n', ...
                refStageZ, stageZnow);
        assert(abs(dStage) <= REUSE_Z_TOL_UM, ...
            ['the stage has moved %.2f um since that reference was taken ' ...
             '(tolerance %.2f um).\nReusing it would redefine z = 0 and every ' ...
             'dz would be offset by that amount.\nEither return the stage to ' ...
             '%.2f um or take a fresh stack.'], ...
            dStage, REUSE_Z_TOL_UM, refStageZ);
    else
        refStageZ = stageZnow;
        warning('dm_loop:noAnchor', ...
            ['that reference has no recorded stage z, so its zero cannot be ' ...
             'verified. Treating the CURRENT stage z (%.2f um) as its zero.'], stageZnow);
    end
    assert(size(refVol,1) == hSI.hRoiManager.linesPerFrame && ...
           size(refVol,2) == hSI.hRoiManager.pixelsPerLine, ...
        ['reused reference is %dx%d but the current geometry is %dx%d. ' ...
         'A reference cannot register frames of a different size.'], ...
        size(refVol,1), size(refVol,2), ...
        hSI.hRoiManager.linesPerFrame, hSI.hRoiManager.pixelsPerLine);
    fprintf('  REUSING the stack from %s\n', REF_REUSE_FROM);
    fprintf('  %d planes, %.1f to %.1f um, taken at linePhase %.3f us\n', ...
            numel(zRef), zRef(1), zRef(end), Mold.linePhaseS*1e6);
    fprintf('  >> valid ONLY if the sample has not moved since. Take a fresh\n');
    fprintf('     stack for any run whose numbers you intend to believe.\n');
else
refStageZ = stageZnow;
fprintf('  anchoring this reference to stage z = %.2f um\n', refStageZ);
arm_scanner(hSI, RES_ARM_SEC, 'the reference stack');

zRef = REF_ZMIN_UM : REF_ZSTEP_UM : REF_ZMAX_UM;
K = numel(zRef);

% Frames per plane from the MEASURED frame rate, so the number of breath cycles
% averaged is the same whatever geometry is loaded. At 256 lines / 60 Hz this
% is twice the count it would be at 512 lines / 30 Hz -- same seconds, same
% cycles, same wall-clock cost.
REF_NFRAMES = max(4, round(REF_SEC_PER_PLANE * fps));
fprintf('  %d planes, %.2f s per plane = %d frames each (%.2f Hz)\n', ...
        K, REF_SEC_PER_PLANE, REF_NFRAMES, fps);
fprintf('  ~%.1f breath cycles per plane at 3 Hz | stack takes ~%.1f min\n', ...
        REF_SEC_PER_PLANE*3, (K*REF_NFRAMES/fps + K*SETTLE_S)/60);

refVol = [];  got = false(1,K);
for k = 1:K
    c = zRef(k) / GAIN_UM_PER_COEFF;
    cmd = baseline + c*zRow;
    if max(abs(cmd)) > CMD_LIMIT
        fprintf('  z %+5.1f um: SKIPPED (|cmd| %.3f)\n', zRef(k), max(abs(cmd)));
        continue
    end
    hDM.senddata(cmd);  pause(SETTLE_S);
    I = quick_avg(hSI, REF_NFRAMES, CHANNEL);
    if isempty(refVol), refVol = zeros([size(I) K],'single'); end
    refVol(:,:,k) = I; %#ok<AGROW>
    got(k) = true;
end
hDM.senddata(baseline); pause(SETTLE_S);
assert(nnz(got) >= 5, 'only %d reference planes acquired', nnz(got));
refVol = refVol(:,:,got);  zRef = zRef(got);
fprintf('  %d planes, %.1f to %.1f um\n', numel(zRef), zRef(1), zRef(end));
end   % REF_REUSE_FROM

% Axial blindness check -- a field of penetrating vessels carries NO axial
% information (measured axial SNR 0.000) and no loop can work on it.
rMid = corr_lin(refVol(:,:,1), refVol(:,:,end));
fprintf('  end-to-end plane correlation r = %.3f\n', rMid);
if rMid > 0.97
    error('dm_loop:axiallyBlind', ...
        ['the field barely changes over %.0f um (r = %.3f) -- it is axially ' ...
         'BLIND, typically a field of penetrating vessels. Move the FOV.'], ...
        zRef(end)-zRef(1), rMid);
end

% pixFrac 0.04 rather than the 0.08 default: measured 2026-09-13 at 512x512,
% 21 planes -- 11.2 ms vs 7.5 ms per frame with IDENTICAL accuracy (+3.00 um
% recovered either way). In a loop, milliseconds are the product.
E = dm_axial_estimator_260913('prep', refVol, zRef, struct('nStrips',1,'pixFrac',0.04));
fprintf('  estimator ready (whole-frame, 1 strip -- the loop sees whole frames)\n');

% ---- can MATLAB keep up? Answer it BEFORE committing the animal. ----
% Every frame the loop misses doubles its delay, which is the one thing the
% correction is fighting. So time the estimate now, against the real frame
% size, and say plainly whether there is headroom.
tMid = refVol(:,:,max(1,round(numel(zRef)/2)));
dm_axial_estimator_260913('run', E, single(tMid), 0);        % warm up
t0 = tic;
for q = 1:20, dm_axial_estimator_260913('run', E, single(tMid), 0); end
msPerFrame = 1000*toc(t0)/20;
budget = 1000/fps;
fprintf('  estimate takes %.1f ms; frame period is %.1f ms (%.0f%% of budget)\n', ...
        msPerFrame, budget, 100*msPerFrame/budget);
if msPerFrame > 0.5*budget
    warning('dm_loop:tooSlow', ...
        ['the estimate uses %.0f%% of the frame period. Frames WILL be missed, ' ...
         'and each miss doubles the delay.\nReduce opts.pixFrac, or the frame ' ...
         'size, before running closed.'], 100*msPerFrame/budget);
end

%% ========================= 2. THE RUN =========================
fprintf('\n--- 2. OPEN %g s -> CLOSED %g s -> OPEN %g s ---\n', SEG_SECONDS);
arm_scanner(hSI, RES_ARM_SEC, 'the run');
hSI.hChannels.loggingEnable = true;

totalSec = sum(SEG_SECONDS);
nFr = round(totalSec * fps);
tEdges = cumsum(SEG_SECONDS);

hSI.hStackManager.framesPerSlice = nFr;
hSI.hScan2D.logFileStem    = sprintf('loop_%s', STAMP);
hSI.hScan2D.logFileCounter = 1;

L = struct();
L.t = nan(1,nFr); L.frameNo = nan(1,nFr); L.dz = nan(1,nFr);
L.cmdUm = nan(1,nFr); L.closed = false(1,nFr); L.serviced = false(1,nFr);

curUm = 0;  e1 = NaN;  e2 = NaN;  lastNo = -1;  nSvc = 0;  nSkip = 0;
nEmpty = 0;

hDM.senddata(baseline); pause(SETTLE_S);
tAcq = tic;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_loop:grabDidNotStart', ...
        'startGrab returned with acqState idle -- is a channel selected for Display?');
end

while ~strcmp(hSI.acqState,'idle')
    tNow = toc(tAcq);
    try
        fno = hSI.hDisplay.lastFrameNumber;
    catch
        pause(0.002); continue
    end
    if fno == lastNo
        pause(0.001);  drawnow limitrate;  continue
    end

    % ---- get the frame ----
    % Never fail silently here. An empty frame used to `continue` straight past
    % the estimate, the log and the progress print, so a run that did nothing
    % at all still reported thousands of "serviced" frames.
    % RETRY THE SAME FRAME, do not discard it. The stripe's frame number
    % advances before its roiData is filled, so the first read after the number
    % changes often lands on a half-built stripe. Committing lastNo before the
    % read turned that race into a lost frame every single time; leaving lastNo
    % alone means we simply look again a few hundred microseconds later.
    [I, why] = get_disp_frame(hSI, CHANNEL);
    if isempty(I)
        nEmpty = nEmpty + 1;
        if nEmpty > 2000                       % ~1 s of retries, not 30 frames
            hSI.abort();
            error('dm_loop:noFrames', ...
                ['no frame for channel %d after %d attempts.\nReason: %s\n' ...
                 'channelDisplay = %s, channelSave = %s'], ...
                CHANNEL, nEmpty, why, ...
                mat2str(hSI.hChannels.channelDisplay(:).'), ...
                mat2str(hSI.hChannels.channelSave(:).'));
        end
        pause(0.0005);
        continue
    end
    nEmpty = 0;

    % Only now commit to this frame.
    if lastNo >= 0, nSkip = nSkip + (fno - lastNo - 1); end
    lastNo = fno;
    nSvc = nSvc + 1;
    if nSvc > nFr, break; end

    % ---- estimate ----
    % INDEPENDENT EVERY FRAME (zInit = NaN), not carried forward. Carrying it
    % runs the fine gradient tier alone, and on the rig 2026-09-13 that froze:
    % once the estimate drifted outside the linear zone the fit returned ~0,
    % it could never climb back, and it reported a confident CONSTANT while a
    % 5 um sinusoid was being injected.
    % NaN forces coarse-then-fine, which is what dm_estimator_check_260913
    % does -- measured 0.993 slope, 0.05 um residual on real tissue, single
    % frames included. One extra decimated correlation per frame, against a
    % budget only 17% used.
    o = dm_axial_estimator_260913('run', E, single(I), NaN);
    dzNow = o.dz(1);

    % ---- decide whether we are in the closed window ----
    inClosed = tNow > tEdges(1) && tNow <= tEdges(2);

    % ---- correction, with velocity extrapolation ----
    if inClosed && isfinite(dzNow)
        % dz is the DISPLACEMENT measured relative to the CURRENT command, so
        % the absolute target is the current command minus the new error.
        e2 = e1;  e1 = dzNow;
        if isfinite(e2), pred = e1 + EXTRAP_GAIN*(e1 - e2); else, pred = e1; end
        want = curUm - LOOP_GAIN * pred;
        want = max(DEFOCUS_MIN_UM, min(DEFOCUS_MAX_UM, want));          % absolute clamp
        want = curUm + max(-MAX_STEP_UM, min(MAX_STEP_UM, want-curUm)); % rate limit
        cmd  = baseline + (want/GAIN_UM_PER_COEFF)*zRow;
        if max(abs(cmd)) <= CMD_LIMIT
            if ~DRY_RUN, hDM.senddata(cmd); end
            curUm = want;
        end
    elseif ~inClosed && curUm ~= 0
        curUm = 0;  e1 = NaN;  e2 = NaN;
        if ~DRY_RUN, hDM.senddata(baseline); end
    end

    L.t(nSvc)=tNow; L.frameNo(nSvc)=fno; L.dz(nSvc)=dzNow;
    L.cmdUm(nSvc)=curUm; L.closed(nSvc)=inClosed; L.serviced(nSvc)=true;

    if mod(nSvc,60)==0
        % Absolute z, not just dz. refStageZ anchors the reference's zero to a
        % real stage coordinate, so zAbs is comparable across sessions instead
        % of being relative to whichever stack happened to be taken last.
        fprintf('  t=%5.1f s  %s  dz %+6.2f  z %8.2f  cmd %+6.2f um  (missed %d)\n', ...
                tNow, ternary(inClosed,'CLOSED','open  '), dzNow, ...
                refStageZ + dzNow, curUm, nSkip);
    end
    drawnow limitrate;
end
if ~DRY_RUN, hDM.senddata(baseline); end

% THE ACQUISITION MUST BE IDLE BEFORE ANYTHING IS RESTORED. The loop above can
% exit on frame count while the grab is still running, and ScanImage then
% refuses every restore with "Cannot set property ... during an active
% acquisition" -- which onCleanup swallows as warnings, leaving the rig on the
% modified geometry and the caller none the wiser.
if ~strcmp(hSI.acqState,'idle')
    fprintf('  waiting for the acquisition to finish...\n');
    t0 = tic;
    while ~strcmp(hSI.acqState,'idle') && toc(t0) < 10
        pause(0.02); drawnow limitrate;
    end
    if ~strcmp(hSI.acqState,'idle')
        fprintf('  still running after 10 s -- aborting it\n');
        hSI.abort();
        while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
    end
end
pause(SETTLE_S);

%% ========================= 3. REPORT =========================
k = L.serviced;
fprintf('\n--- 3. result ---\n');
fprintf('serviced %d frames, MISSED %d (%.1f%%)\n', nSvc, nSkip, 100*nSkip/max(nSvc+nSkip,1));
if nSkip > 0.1*nSvc
    fprintf('  >> MATLAB is not keeping up. Every missed frame doubles the delay,\n');
    fprintf('     which is exactly what the loop is fighting. Reduce the frame\n');
    fprintf('     size or the pixel fraction before trusting the CLOSED numbers.\n');
end

op = k & ~L.closed;  cl = k & L.closed;
if nnz(op) > 10 && nnz(cl) > 10
    fprintf('\n  dz measured BY THE LOOP (indicative only -- it is its own referee):\n');
    fprintf('    open   RMS %.3f um   (n=%d)\n', rms_nan(L.dz(op)), nnz(op));
    fprintf('    closed RMS %.3f um   (n=%d)\n', rms_nan(L.dz(cl)), nnz(cl));
    fprintf('    correction actually applied, closed: %.3f um RMS, %.2f um p2p\n', ...
            rms_nan(L.cmdUm(cl)), max(L.cmdUm(cl))-min(L.cmdUm(cl)));
    if DRY_RUN && any(L.cmdUm(cl) <= DEFOCUS_MIN_UM + 1e-6 | ...
                      L.cmdUm(cl) >= DEFOCUS_MAX_UM - 1e-6)
        fprintf(['\n  NOTE: the command saturated at the clamp. In a DRY RUN that is\n' ...
                 '  EXPECTED, not a fault: nothing is sent, so the mirror never moves,\n' ...
                 '  so the measured error never responds, so the integrator keeps\n' ...
                 '  winding until it hits the limit. Whether the loop is stable can\n' ...
                 '  only be judged with DRY_RUN = false.\n']);
    end
end
fprintf(['\n  THE REAL TEST IS OFFLINE, on the logged frames:\n' ...
         '    dm_axial_loop_eval_260913(''%s'')\n' ...
         '  It re-measures dz on every logged frame against the same reference,\n' ...
         '  blind to which window was which, and reports CLOSED vs OPEN.\n'], sessDir);

meta = struct('stamp',STAMP,'sessDir',sessDir,'channel',CHANNEL,'fps',fps, ...
    'linesPerFrame',hSI.hRoiManager.linesPerFrame, ...
    'pixelsPerLine',hSI.hRoiManager.pixelsPerLine, ...
    'scanPhaseUs',SCAN_PHASE_US,'linePhaseS',hSI.hScan2D.linePhase, ...
    'linePhaseEntryS',siState.linePhase, ...
    'targetFps',TARGET_FPS, ...
    'segSeconds',SEG_SECONDS,'tEdges',tEdges,'loopGain',LOOP_GAIN, ...
    'extrapGain',EXTRAP_GAIN,'defocusMinUm',DEFOCUS_MIN_UM, ...
    'defocusMaxUm',DEFOCUS_MAX_UM,'maxStepUm',MAX_STEP_UM, ...
    'dryRun',DRY_RUN,'gainUmPerCoeff',GAIN_UM_PER_COEFF,'baseline',baseline, ...
    'baselineVar',BASELINE_VAR,'rowDefoc',ROW_DEFOC,'zRef',zRef,'refVol',refVol, ...
    'refStageZ',refStageZ, ...   % the stage z this reference's zero means
    'stageZatRun',hSI.hMotors.motorPosition(3), ...
    'log',L,'nServiced',nSvc,'nMissed',nSkip, ...
    'file',fullfile(sessDir,sprintf('loop_%s_00001.tif',STAMP)));
save(fullfile(sessDir,'loop_meta.mat'),'-struct','meta','-v7.3');
fprintf('\nsaved %s\n', fullfile(sessDir,'loop_meta.mat'));
fprintf('====================================================================\n');
end

%% ============================ HELPERS ============================
function [I, why] = get_disp_frame(hSI, ch)
%GET_DISP_FRAME  One channel out of the live display buffer.
%
% *** lastFrame IS INDEXED BY channelDisplay, NOT BY lastFrameChannels. ***
% From Display.m's getter:
%     [~,idxs] = ismember(hSI.hChannels.channelDisplay, stripe.roiData{1}.channels);
%     val = stripe.roiData{1}.imageData(idxs);
% so lastFrame{i} is channelDisplay(i). But lastFrameChannels returns
% stripe.roiData{1}.channels, which is a DIFFERENT list. Finding an index in
% one and subscripting the other happens to work when the two lists coincide
% and silently returns the wrong channel -- or nothing -- when they do not.
% ScanImage's own comment on that getter says why they diverge: "fix for
% inconsistent behavior between resscan and linscan when some channels are
% enabled for logging only". Which is exactly the case here.
%
% Two paths are tried, and the reason for failure is RETURNED rather than
% swallowed -- an empty frame used to be indistinguishable from a missing one.
% *** THE RACE THIS EXISTS TO DODGE ***
% lastStripeData advances its frame NUMBER before its roiData is populated, so
% polling the instant the number changes lands on a half-built stripe and
% `stripe.roiData{1}` throws "Cell contents reference from a non-cell array
% object". Every time, because the poll is tight. quick_avg never sees it
% because it reads after the grab has finished.
%
% rollingStripeDataBuffer{zIdx}{1} is the COMPLETED averaged stripe
% (Display.m:891 writes it, :876 checks it), so it is read first. The caller
% also retries the same frame rather than discarding it -- see the loop.
I = [];  why = '';

% 1. the completed frame
try
    b = hSI.hDisplay.rollingStripeDataBuffer;
    if ~isempty(b) && ~isempty(b{1}) && ~isempty(b{1}{1}.roiData)
        rd = b{1}{1}.roiData{1};
        I  = pick_channel(rd, ch);
        if ~isempty(I), return; end
        why = sprintf('rolling buffer has channels %s, wanted %d', ...
                      mat2str(rd.channels(:).'), ch);
    else
        why = 'rolling buffer empty';
    end
catch ME
    why = sprintf('rolling buffer threw: %s', ME.message);
end

% 2. the live stripe (may be mid-build; that is what the retry is for)
try
    rd = hSI.hDisplay.lastStripeData.roiData{1};
    I  = pick_channel(rd, ch);
    if ~isempty(I), why = ''; return; end
catch ME
    why = sprintf('%s | live stripe: %s', why, ME.message);
end

% 3. the convenience wrapper, indexed by channelDisplay (NOT lastFrameChannels)
try
    dispCh = hSI.hChannels.channelDisplay;
    i = find(dispCh(:).' == ch, 1);
    f = hSI.hDisplay.lastFrame;
    if ~isempty(i) && i <= numel(f) && ~isempty(f{i})
        I = single(f{i});  why = '';  return
    end
catch ME
    why = sprintf('%s | lastFrame: %s', why, ME.message);
end
end

function I = pick_channel(rd, ch)
%PICK_CHANNEL  One channel's image out of a roiData, or [] if not ready.
I = [];
j = find(rd.channels(:).' == ch, 1);
if isempty(j) || j > numel(rd.imageData), return; end
im = rd.imageData{j};
if iscell(im), if isempty(im), return; end, im = im{1}; end
if isempty(im), return; end
if isprop(rd,'transposed') && rd.transposed, im = im.'; end
I = single(im);
end

function I = quick_avg(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_loop:grabDidNotStart', ...
        'startGrab returned with acqState idle -- is a channel selected for Display?');
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
[I, why] = get_disp_frame(hSI, ch);
assert(~isempty(I), 'channel %d not readable from the display buffer: %s', ch, why);
I = single(I);
hSI.hDisplay.displayRollingAverageFactor = 1;
end

function b = get_baseline(varName)
if ~evalin('base', sprintf('exist(''%s'',''var'')', varName))
    error('dm_loop:noBaseline', ...
        ['base workspace has no ''%s''. Load and activate the system ' ...
         'aberration correction first -- there is no read-back path on this ' ...
         'hardware, so the loaded variable is the only honest baseline.'], varName);
end
b = double(evalin('base', varName));  b = b(:).';
assert(numel(b) == 97, '%s is %d long, expected 97', varName, numel(b));
end

function arm_scanner(hSI, sec, what)
hSI.hScan2D.keepResonantScannerOn = true;
fprintf('  arming resonant scanner %.1f s before %s...\n', sec, what);
t0 = tic; while toc(t0) < sec, pause(0.05); drawnow limitrate; end
end

function r = corr_lin(a, b)
a = double(a(:)) - mean(double(a(:)));
b = double(b(:)) - mean(double(b(:)));
r = (a.'*b)/(norm(a)*norm(b)+eps);
end

function v = rms_nan(x)
x = x(isfinite(x));  x = x - mean(x);
v = sqrt(mean(x.^2));
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function restore_all(hDM, baseline, hSI, s, p0) %#ok<INUSD>
try, hDM.senddata(baseline); fprintf('[cleanup] DM restored to baseline\n');
catch ME, warning('dm_loop:noDM','COULD NOT RESTORE THE DM: %s', ME.message); end
try
    % ROI geometry first, and forceSquarePixelation LAST of the three -- it
    % rewrites linesPerFrame when it is turned back on, so restoring it before
    % the line count would clobber the value we are trying to put back.
    hSI.hRoiManager.forceSquarePixels     = s.forceSquarePixels;
    hSI.hRoiManager.linesPerFrame         = s.linesPerFrame;
    hSI.hRoiManager.forceSquarePixelation = s.forceSquarePixelation;
    fprintf('[cleanup] frame geometry restored (%d lines)\n', s.linesPerFrame);
catch ME, warning('dm_loop:noRoi','could not restore ROI geometry: %s', ME.message); end
try
    % Unconditional: it was captured before the closure, so it is always here.
    hSI.hScan2D.linePhase = s.linePhase;
    fprintf('[cleanup] linePhase restored to %.3f us\n', s.linePhase*1e6);
catch ME, warning('dm_loop:noPhase','could not restore scan phase: %s', ME.message); end
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
catch ME, warning('dm_loop:noSI','could not restore SI settings: %s', ME.message); end
end
