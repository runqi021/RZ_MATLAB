function dm_motion_calib_260911()
%DM_MOTION_CALIB_260911  Measure how far the sample moves per micron of DM Zernike.
%
% STAGE 0 of the AO motion-correction pipeline.  Produces the only numbers a
% corrector actually needs: the 3x3 matrix G mapping
%
%       [tiltY tipX defocus] coefficients (um)  ->  [dx dy dz] at the sample (um)
%
% so a corrector inverts it:  coeffs = (dr_measured) / G.
%
% It also answers, empirically and against YOUR objective, the two things that
% cannot be read off the Z2C file:
%   - whether the coefficients are um RMS or um PtV (a factor of ~4 in range)
%   - where the response stops being linear and the PSF starts to go
%
% WHAT IT DOES
%   for each of tilt / tip / defocus:
%     save the current DM command, step a coefficient through COEFFS_UM,
%     grab and average NFRAMES at each step, restore the baseline.
%   then: lateral shift by phase correlation against the c = 0 frame;
%         axial shift by matching each frame into a baseline z-stack.
%
% ================================ SAFETY ================================
% Every send is baseline + motion, and the baseline (your system + sample
% aberration correction) is restored in the cleanup handler even on error or
% Ctrl-C.  Commands are checked against the ALPAO +-1 limit BEFORE being sent;
% a step that would clip is skipped, not clamped, because a clamped mirror is
% no longer producing the Zernike you think it is and would silently corrupt
% the fit.  Start with small COEFFS_UM and widen once you have seen it behave.
%
% REQUIRES: ScanImage running (hSI, hAOROIctrl in base), an FOV with structure
% that cross-correlates well (beads, or vasculature), and NO breathing motion
% during the run -- do it on a dead/anaesthetised-still prep or a slide, or the
% motion you are trying to measure gets swamped by the motion you are trying to
% correct.
%
% *** NOT YET RUN ON THE RIG.  Written against the ScanImage 2018b source; the
% *** hardware calls are unexecuted.  Read it before you trust it.
%
% Runqi Zhang / 2026-09-11

%% ========================= USER SETTINGS =========================
SI_DIR     = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE   = fullfile(SI_DIR,'BAX331-Z2C.mat');

% THE BASELINE COMES FROM THE BASE WORKSPACE, BY NAME.
% There is no way to read the current command back off the DM -- DMctrl caches
% nothing (its properties are hDM/Gain/serialName) and hAOdata has no such
% getter. The lab's own scripts already work this way, e.g.
%   Layer4_aberration_correction.m:20  senddata(DMcommand_spherical + DMcommand_sys_aber_cor_930)
%   system_aberration_correction1.m:33 senddata(DMcommand_defocus  + DMcommand_sys_aber_cor_930)
% so load your correction first, then run this. The variable is [1 97] and
% lives in e.g. water_lens_System_aberration_correction_8mm\
%   SystemCorrection_20260117_930_Iter30_FITC_Z60_P3_20K_water.mat
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_TILTY  = 1;        % verified from actuator geometry: row 1 projects 0.988 onto a y-ramp
ROW_TIPX   = 2;        %                                  row 2 projects 0.986 onto an x-ramp
ROW_DEFOC  = 3;        %                                  row 3 projects 0.959 onto r^2

COEFFS_UM  = [-2 -1 -0.5 0 0.5 1 2];   % start small.  Widen after the first run.
NFRAMES    = 15;       % frames averaged per step (sets framesPerSlice AND
                       % displayRollingAverageFactor -- see grab_avg)
SETTLE_S   = 0.3;      % after a send, before grabbing
% WHICH CHANNEL TO REGISTER ON. Name it explicitly -- [] means "take the first
% frame in the buffer", which with 3 channels open is ch1, and ch1 on this rig is
% near-featureless (see project_autostitch_ncc_260910: "no visible cells at all,
% diffuse haze with vessel silhouettes, contrast index 0.02"). The 1813 run did
% exactly that and every correlation came back 0.02. Use the VESSEL channel.
CHANNEL    = 3;

DO_LATERAL = true;
DO_AXIAL   = true;     % needs a z-stack; set false for a quick lateral-only pass
% Total z range for the axial reference stack, centred on the current focus.
% 50 = +-25 um, matching the depth over which the user confirmed good sample on
% 2026-09-11. THE STACK MUST SPAN WHEREVER THE DEFOCUS SWEEP PUSHES THE FOCUS:
% if a coefficient drives the focus past the end of the stack, match_z finds its
% best plane at the boundary and interp1 EXTRAPOLATES, which yields a confident
% but wrong dz. The script now warns when that happens rather than reporting it
% silently -- if you see that warning, narrow COEFFS_UM or widen this.
ZSTACK_UM  = 50;

% *** ZSTEP MUST BE A MULTIPLE OF THE 0.3125 um Z STEPPER QUANTUM ***
% Measured 2026-09-11: z moves in 0.3125 um counts, so a 1 um command is 3.2
% counts and the stage actually steps 0.9375 um -- 6.25% SHORT. This stack is
% the AXIAL RULER for the defocus gain, so that error would land directly in G.
% 1.25 um is exactly 4 counts. Other exact values: 0.9375 (3), 2.5 (8), 5 (16).
% See project_stage_backlash_260911.
ZSTEP_UM   = 1.25;

CMD_LIMIT  = 0.95;     % refuse to send if any actuator would exceed this

% Keep the SWEPT FRAMES in the .mat (7 frames, a few MB). The z-stack is NOT
% saved -- user's call, it is the 40 MB bulk of the file. Set false to store
% nothing but the fit.
SAVE_IMAGES = true;

% HOLD THE RESONANT SCANNER ON ACROSS THE WHOLE RUN.
% This script issues ~48 separate grabs (41 stack planes + one per coefficient).
% With keepResonantScannerOn false the scanner spins down after each and re-arms
% for the next, and arm() only ESTIMATES the scan frequency -- the measure-after-
% settle branch is commented out in SI -- so there is no settle wait in the grab
% path and the first frames of every grab come out distorted. Holding it on
% removes that entirely. Armed once, RES_ARM_SEC before the first grab.
% 1.0 s per the user (2026-09-11); SI's own number is
% resonantScannerWaitSettle(max(2, mdfData.resonantScannerSettleTime)).
RES_ARM_SEC = 1.0;
OUT_DIR    = fileparts(mfilename('fullpath'));

%% ========================= CONNECT =========================
hSI   = evalin('base','hSI');
hAO   = evalin('base','hAOROIctrl');
hDM   = hAO.hDMctrl;

% Raw Galil object: the ONLY route to the z encoder (TP). ScanImage's
% motorPosition is TD, a step-pulse count that echoes the command.
hLSC = [];
try, hLSC = hSI.hMotors.hMotor(1).hLSC; catch, end
if isempty(hLSC)
    warning('dm_calib:noLSC', ...
        'no raw Galil object -- the z ruler will fall back to commanded position');
end

S = load(Z2C_FILE);
fn = fieldnames(S); Z2C = S.(fn{1});
assert(size(Z2C,2) == 97, 'expected 97 actuators, got %d', size(Z2C,2));

rows  = [ROW_TILTY ROW_TIPX ROW_DEFOC];
names = {'tiltY','tipX','defocus'};

% MODES: which of the three to sweep. 3 = defocus only.
% Defocus is the one that matters -- axial motion is the only component that
% post-hoc registration cannot undo, so tip/tilt are diagnostic rather than
% essential. Start with 3 alone: fastest run, fewest DM commands, and it answers
% the actual question (does commanding defocus move the focus, and how much per
% unit coefficient). Widen to [1 2 3] once that behaves, to get the cross-talk.
MODES = 3;

rows  = rows(MODES);
names = names(MODES);
nM    = numel(rows);

% The baseline is whatever the DM is holding now -- your aberration correction.
% There is no "read current command" call on the ALPAO driver, so it is read
% from the BASE WORKSPACE by name and must be re-sent by us at the end.
baseline = get_baseline(BASELINE_VAR);

% grab_avg calls startGrab, which (a) WRITES A TIFF PER GRAB if logging is on --
% ~40 files for one run of this script -- and (b) overwrites framesPerSlice and
% numSlices. This is a calibration, not data. Turn logging off, remember all
% three, and put them back in the SAME cleanup that restores the DM, so Ctrl-C
% leaves the rig as it was found.
siState = struct( ...
    'loggingEnable',  hSI.hChannels.loggingEnable, ...
    'framesPerSlice', hSI.hStackManager.framesPerSlice, ...
    'numSlices',      hSI.hStackManager.numSlices, ...
    'keepResOn',      hSI.hScan2D.keepResonantScannerOn, ...
    'rollingAvg',     hSI.hDisplay.displayRollingAverageFactor);
if siState.loggingEnable
    fprintf('logging was ON -- disabled for this run, restored afterwards\n');
end
hSI.hChannels.loggingEnable = false;

restore = onCleanup(@() restore_all(hDM, baseline, hSI, siState));

% Arm the resonant scanner ONCE and hold it on for the whole run. Note this flag
% is persisted to SI's class data file (it is the cbKeepScannerOn checkbox), so
% the cleanup below matters -- otherwise the rig is left with it changed.
hSI.hScan2D.keepResonantScannerOn = true;
fprintf('arming resonant scanner (%.1f s), held on for the run...\n', RES_ARM_SEC);
t0 = tic; while toc(t0) < RES_ARM_SEC, pause(0.05); drawnow limitrate; end
fprintf('baseline: %s, max|cmd| = %.3f of the +-1 stroke\n', ...
        BASELINE_VAR, max(abs(baseline)));

% NO um/px CONVERSION. Lateral results are reported in PIXELS on purpose.
% The corrector needs "image shift -> Zernike coefficient", so if the gain is
% measured in px/coefficient and the motion estimator also reports px, the unit
% cancels and no ruler is needed. That sidesteps two separate problems:
% objectiveResolution = 38.9 is ONE number for both axes while the stitching
% work measured x and y disagreeing ~6%, and (measured 2026-09-11) the x/y stage
% is open-loop with an unverified 0.78125 um/pulse constant and no encoder, so
% it is not a length standard either. Axial stays in um, where the z stack is
% the ruler and z DOES have a 0.1 um encoder.

%% ========================= AXIAL REFERENCE =========================
zs = []; stack = [];
if DO_AXIAL
    % %g not %d -- ZSTEP_UM is 1.25, and %d on a non-integer prints in exponential
    fprintf('\n[z-stack] baseline reference, %g um at %g um steps (%g stepper counts each)\n', ...
            ZSTACK_UM, ZSTEP_UM, ZSTEP_UM/0.3125);
    [stack, zs, zcnt] = grab_zstack(hSI, hLSC, hDM, baseline, ZSTACK_UM, ZSTEP_UM, ...
                                    NFRAMES, CHANNEL, SETTLE_S);
end

%% ========================= SWEEP =========================
nC = numel(COEFFS_UM);
IM = cell(nM, nC);  applied = nan(nM, nC);

for m = 1:nM
    fprintf('\n[%s] row %d\n', names{m}, rows(m));
    for k = 1:nC
        c   = COEFFS_UM(k);
        cmd = baseline(:).' + c * Z2C(rows(m),:);
        if max(abs(cmd)) > CMD_LIMIT
            fprintf('   c = %+5.2f um  SKIPPED, would reach |cmd| = %.3f\n', c, max(abs(cmd)));
            continue
        end
        hDM.senddata(cmd);
        pause(SETTLE_S);
        IM{m,k} = grab_avg(hSI, NFRAMES, CHANNEL);
        applied(m,k) = c;
        fprintf('   c = %+5.2f um  max|cmd| %.3f  mean %.0f\n', c, max(abs(cmd)), mean(IM{m,k}(:)));
    end
    hDM.senddata(baseline);  pause(SETTLE_S);
end

%% ========================= ANALYSE =========================
% G(m,:) = [dx dy dz] per micron of coefficient m.  Off-diagonals ARE the
% cross-talk, and the corrector must invert the whole matrix -- three separate
% scalars would leave tip-induced defocus uncorrected.
G = nan(nM,3);  fits = struct();
for m = 1:nM
    k0 = find(applied(m,:) == 0, 1);
    assert(~isempty(k0) && ~isempty(IM{m,k0}), 'no c = 0 frame for %s', names{m});
    ref = IM{m,k0};
    dx = nan(1,nC); dy = nan(1,nC); dz = nan(1,nC);
    for k = 1:nC
        if isempty(IM{m,k}), continue; end
        if DO_LATERAL
            d = phase_shift(ref, IM{m,k});
            dx(k) = d(2);  dy(k) = d(1);     % PIXELS. [row col] -> [y x]
        end
        if DO_AXIAL && ~isempty(stack)
            [dz(k), edgeHit, rpk] = match_z(stack, zs, IM{m,k});
            if edgeHit
                warning('dm_calib:zEdge', ...
                    ['%s c = %+g: best plane is at the EDGE of the z-stack. The ' ...
                     'focus moved OUTSIDE it, so dz = %.2f is only a LOWER BOUND. ' ...
                     'Narrow COEFFS_UM or widen ZSTACK_UM.'], names{m}, applied(m,k), dz(k));
            end
            if rpk < 0.3
                warning('dm_calib:weakMatch', ...
                    ['%s c = %+g: best correlation into the stack is only %.2f -- ' ...
                     'weak match, dz is unreliable. Needs more texture in the FOV.'], ...
                    names{m}, applied(m,k), rpk);
            end
        end
    end
    ok = ~isnan(applied(m,:));
    G(m,1) = robust_slope(applied(m,ok), dx(ok));
    G(m,2) = robust_slope(applied(m,ok), dy(ok));
    G(m,3) = robust_slope(applied(m,ok), dz(ok));
    fits.(names{m}) = struct('c',applied(m,:),'dx',dx,'dy',dy,'dz',dz);
    fprintf('\n%-8s -> dx %+7.3f  dy %+7.3f  dz %+7.3f   um per um-coefficient\n', ...
            names{m}, G(m,1), G(m,2), G(m,3));
end

%% ========================= REPORT =========================
fprintf('\n================ GAIN MATRIX G ================\n');
fprintf('            dx        dy        dz     (dx,dy in PX; dz in um; per um coeff)\n');
for m = 1:nM
    fprintf('%-8s %8.3f  %8.3f  %8.3f\n', names{m}, G(m,1), G(m,2), G(m,3));
end
fprintf('\ncondition number %.2f  (large = axes are not separable)\n', cond(G));

% Range: how far can each axis be pushed before an actuator clips, using the
% measured gain.  This is the number that decides whether the DM can cover your
% breathing amplitude at all.
fprintf('\n================ REACHABLE RANGE ================\n');
for m = 1:nM
    cmax = (CMD_LIMIT - max(abs(baseline))) / max(abs(Z2C(rows(m),:)));
    % USE THE AXIAL COLUMN ONLY. G(:,1:2) are PIXELS and G(:,3) is MICRONS, so
    % max(abs(G(m,:))) silently picked whichever was numerically larger and
    % reported a pixel count as microns -- the 1821 run printed "+15.9 um" from
    % a 10.19 px lateral term when the real axial range was 1.56 * 6.35 = 9.9 um.
    fprintf('%-8s coeff +-%.2f um  ->  %+.1f um AXIAL  (lateral %+.1f, %+.1f px)\n', ...
            names{m}, cmax, cmax*G(m,3), cmax*G(m,1), cmax*G(m,2));
end
fprintf('(baseline aberration correction already uses max|cmd| = %.3f of the +-1 stroke)\n', ...
        max(abs(baseline)));

out = fullfile(OUT_DIR, sprintf('dm_motion_calib_%s.mat', datestr(now,'yymmdd_HHMM')));
% SAVE THE IMAGES, NOT JUST THE FIT.
% The 1813 run saved only G and fits. When it turned out to have registered the
% wrong channel on unaveraged frames, there was nothing left to re-analyse -- the
% z-stack and every swept frame had been discarded, and logging was off so no
% TIFFs existed either. The whole run had to be repeated. Same lesson as the
% stage work: keep the raw data, derive at analysis time.
% See feedback_measure_raw_not_converted.
save(out, 'G','fits','COEFFS_UM','rows','names','baseline','BASELINE_VAR', ...
     'Z2C_FILE','ZSTEP_UM','ZSTACK_UM','CHANNEL','NFRAMES','MODES');
if SAVE_IMAGES
    % THE Z-STACK IS DELIBERATELY NOT SAVED (user, 2026-09-11) -- it is ~40 MB
    % per run and by far the bulk of it. Its z AXIS is kept (zs/zcnt, the
    % encoder-measured ruler in um and raw counts), as are the swept frames,
    % which are small and are what you actually inspect when a run looks wrong.
    % Trade-off to know: without the stack, a bad run cannot be re-matched
    % offline against a different z reference -- it has to be repeated.
    fprintf('saving %d swept frames (z-stack intentionally not saved)...\n', ...
            nnz(~cellfun(@isempty,IM)));
    save(out, 'zs','zcnt','IM','applied', '-append');
end
fprintf('\nsaved %s\n', out);
plot_fits(fits, names, G, OUT_DIR);
end

%% ========================== HELPERS ==========================
function b = get_baseline(varName)
% The command the DM is currently holding, taken from the BASE WORKSPACE by name.
%
% There is deliberately NO fallback. An earlier version called
% hAO.hAOdata.getDMdata(), which does not exist anywhere in the install, so the
% try/catch fired every time and silently substituted zeros(1,97). That is the
% worst possible failure: the gain would be measured about a FLAT mirror instead
% of about your correction, the CMD_LIMIT headroom would be computed against
% max|baseline| = 0, and onCleanup would then send zeros and WIPE the aberration
% correction off the DM at the end of the run.
%
% Reading it back from the hardware is not an option: DMctrl caches nothing and
% the ALPAO driver has no read command. The variable you loaded IS the record.
if ~evalin('base', sprintf('exist(''%s'',''var'')', varName))
    error('dm_calib:noBaseline', ...
        ['base workspace has no variable ''%s''.\nLoad and activate your system ' ...
         'aberration correction FIRST -- every command this script sends is\n' ...
         'baseline + motion, and the baseline cannot be read back from the DM.'], ...
        varName);
end
b = double(evalin('base', varName));
b = b(:).';
if numel(b) ~= 97
    error('dm_calib:baselineSize', ...
        '%s is %d long, expected 97 actuators', varName, numel(b));
end
if max(abs(b)) > 1
    error('dm_calib:baselineRange', ...
        '%s already exceeds the ALPAO +-1 limit (max|cmd| = %.3f)', ...
        varName, max(abs(b)));
end
end

function restore_all(hDM, baseline, hSI, s)
% DM FIRST -- that is the safety-critical one. SI settings after, and a failure
% of either must not stop the other being attempted.
safe_restore(hDM, baseline);
try
    hSI.hChannels.loggingEnable                = s.loggingEnable;
    hSI.hStackManager.framesPerSlice           = s.framesPerSlice;
    hSI.hStackManager.numSlices                = s.numSlices;
    hSI.hScan2D.keepResonantScannerOn          = s.keepResOn;
    hSI.hDisplay.displayRollingAverageFactor   = s.rollingAvg;
    fprintf('[cleanup] ScanImage logging/stack/scanner/averaging settings restored\n');
catch ME
    warning('dm_calib:noSIRestore', ...
        'could not restore ScanImage settings: %s', ME.message);
end
end

function safe_restore(hDM, baseline)
try, hDM.senddata(baseline); fprintf('[cleanup] DM restored to baseline\n');
catch ME, warning('dm_calib:noRestore','COULD NOT RESTORE THE DM: %s', ME.message);
end
end

function u = get_um_per_px_mdf(hSI)
% ScanImage's own model: one objectiveResolution (38.9 um/deg) scaling both axes.
% Used only as a sanity check against the measured value -- see stage_ruler.
sf = hSI.hRoiManager.currentRoiGroup.rois(1).get(0);
u  = [sf.sizeXY(1) * hSI.objectiveResolution / sf.pixelResolutionXY(1), ...
      sf.sizeXY(2) * hSI.objectiveResolution / sf.pixelResolutionXY(2)];
end

function u = stage_ruler(hSI, n, ch, stepUm)
% MEASURE um/px with the stage, per axis, instead of trusting the MDF.
%
% Why not the MDF: objectiveResolution is a single number scaling both axes, but
% the stitching work measured x and y disagreeing by ~6% over 144 tile pairs --
% the resonant and galvo axes are not calibrated together.  One scalar is
% therefore ~6% wrong on one axis, and that error would land directly in the DM
% gain.  The stage is the ruler you already trust: 0.7815 um/device-unit, and
% commanded positions land where they are asked.
%
% No calibration grid needed.  Move a known distance, see how far the image went.
p0 = hSI.hMotors.motorPosition;
u  = nan(1,2);
for ax = 1:2
    I = cell(1,2);
    for s = [-1 1]
        p = p0;  p(ax) = p0(ax) + s*stepUm/2;
        hSI.hMotors.motorPosition = p;  pause(0.4);
        I{(s+3)/2} = grab_avg(hSI, n, ch);
    end
    hSI.hMotors.motorPosition = p0;  pause(0.4);
    d = phase_shift(I{1}, I{2});          % [dy dx] in px for a stepUm stage move
    if ax == 1, px = abs(d(2)); else, px = abs(d(1)); end
    u(ax) = stepUm / px;
    fprintf('   axis %d: %.1f um of stage -> %.2f px  =>  %.4f um/px\n', ax, stepUm, px, u(ax));
end
end

function I = grab_avg(hSI, n, ch)
% AVERAGING COMES FROM YOUR DISPLAY SETTING, NOT FROM THIS SCRIPT.
% hDisplay.lastFrame is "the last frame in the image buffer" (Display.m:39) -- a
% SINGLE frame. framesPerSlice alone grabs n frames and hands back the last one,
% unaveraged; displayRollingAverageFactor is what averages. The 1813 run
% registered single noisy frames for exactly this reason.
%
% So the script sets displayRollingAverageFactor = NFRAMES itself, and restores
% whatever you had in the cleanup handler (including on Ctrl-C). Both numbers
% must match: the effective average is min(framesPerSlice, rollingAverage), so
% setting only one of them silently averages over fewer frames than intended.
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;

hSI.hDisplay.displayRollingAverageFactor = n;

hSI.startGrab();
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end

f = hSI.hDisplay.lastFrame;
chans = hSI.hDisplay.lastFrameChannels;
if isempty(ch)
    idx = 1;
else
    idx = find(chans == ch, 1);
    if isempty(idx)
        error('dm_calib:noChannel', ...
            ['channel %d is not being displayed (available: %s). Select it in ' ...
             'ScanImage, or change CHANNEL.'], ch, mat2str(chans));
    end
end
I = single(f{idx});
end

function [stack, zs, zcnt] = grab_zstack(hSI, hLSC, hDM, baseline, range_um, step_um, n, ch, settle)
% Baseline z-stack with the STAGE. This stack is the AXIAL RULER, so its z axis
% must be right.
%
% z is READ FROM THE ENCODER, not assumed from the command. The z axis carries a
% real 0.1 um encoder on Galil TP (measured 2026-09-11) which ScanImage never
% reads -- motorPosition is TD, a step-pulse count that merely echoes the
% command. Reading TP removes BOTH error sources from the ruler at once:
%   - the 0.3125 um stepper quantum (a 1 um command really moves 0.9375)
%   - the 0.205 um backlash on any direction reversal
% and it means step_um no longer has to be an exact multiple of anything.
% Raw counts are returned alongside so the conversion can be revisited without
% re-running -- see feedback_measure_raw_not_converted.
% SIGN MATTERS. The fit on 2026-09-11 gave -10.0127 counts per um of SI z, i.e.
% the encoder counts DOWN as SI z goes UP. Using +0.09987 made the 1821 run's z
% axis run -25 to -74 um instead of -25 to +25 -- backwards, so every absolute
% dz was wrong (the SLOPE, and therefore the gain, survived it).
ENC_UM_PER_COUNT = -0.09987;     % measured; x/y have no encoder, z does

z0    = hSI.hMotors.motorPosition(3);
zsCmd = (-range_um/2 : step_um : range_um/2);
hDM.senddata(baseline); pause(settle);

% Run-up below the first plane, so the sweep is monotonic from there on and the
% first plane is not the one that eats the reversal backlash.
p = hSI.hMotors.motorPosition;  p(3) = z0 + zsCmd(1) - 5*step_um;
hSI.hMotors.motorPosition = p;  pause(0.3);

stack = [];  zcnt = nan(1,numel(zsCmd));
for k = 1:numel(zsCmd)
    p = hSI.hMotors.motorPosition;  p(3) = z0 + zsCmd(k);
    hSI.hMotors.motorPosition = p;  pause(0.2);
    zcnt(k) = read_z_counts(hLSC);
    I = grab_avg(hSI, n, ch);
    if isempty(stack), stack = zeros([size(I) numel(zsCmd)],'single'); end
    stack(:,:,k) = I; %#ok<AGROW>
end
p = hSI.hMotors.motorPosition;  p(3) = z0;  hSI.hMotors.motorPosition = p;
pause(0.3);

if all(isfinite(zcnt)) && range(zcnt) > 0
    zs = (zcnt - zcnt(1)) * ENC_UM_PER_COUNT + zsCmd(1);
    resid = zs - zsCmd;
    fprintf(['   z ruler from ENCODER: commanded %g um/step, measured %.4f um/step\n' ...
             '   max |measured - commanded| = %.3f um over the stack\n'], ...
            step_um, (zs(end)-zs(1))/(numel(zs)-1), max(abs(resid - mean(resid))));
else
    zs = zsCmd;
    warning('dm_calib:noZEncoder', ...
        ['could not read the z encoder; falling back to COMMANDED z. The ruler ' ...
         'is then quantised to 0.3125 um steps and the defocus gain inherits ' ...
         'that error.']);
end
end

function c = read_z_counts(hLSC)
% Galil TP on axis C, raw counts. NaN if unreadable -- never silently fall back
% to TD, which is the command echo and would defeat the point.
c = NaN;
try
    v = hLSC.getEncoderPosition();
    if numel(v) >= 3, c = v(3); end
catch
end
end

function d = phase_shift(ref, im)
% Subpixel translation by phase correlation.  Hann-windowed: the frame edges
% move content in and out of view, and an unwindowed FFT reads that as signal.
%
% SIGN, verified against synthetic translations: if `im` is `ref` displaced by
% +5 px, this returns -5.  So d is the shift that maps `im` BACK onto `ref` --
% the correction, not the displacement.  Get this backwards in a corrector and
% it drives the motion instead of cancelling it, doubling the error rather than
% removing it, and the closed loop runs away.  Measured accuracy: exact at zero
% shift, ~0.17 px at +-11 px (the window loses content at large shifts).
w = hann(size(ref,1)) * hann(size(ref,2)).';
A = fft2(double(ref - mean(ref(:))) .* w);
B = fft2(double(im  - mean(im(:)))  .* w);
R = A .* conj(B);  R = R ./ max(abs(R), eps);
c = real(ifft2(R));
[~, i] = max(c(:));
[r, cc] = ind2sub(size(c), i);
% -1 because MATLAB indexes from 1: the zero-shift peak sits at (1,1) and must
% come out as a shift of 0.  Without it every measurement is off by one pixel,
% which is a constant offset in the fit -- invisible in the slope but wrong in
% every reported displacement.
d = [cent(c, r, 1), cent(c, cc, 2)] - 1;
d = d - (d > size(c)/2) .* size(c);        % wrap to signed shift
end

function v = cent(c, i, dim)
% 3-point parabolic peak interpolation along one axis.
n = size(c, dim);
im1 = mod(i-2, n) + 1;  ip1 = mod(i, n) + 1;
if dim == 1, y = [c(im1,1) c(i,1) c(ip1,1)]; else, y = [c(1,im1) c(1,i) c(1,ip1)]; end
den = y(1) - 2*y(2) + y(3);
if den == 0, v = i; else, v = i + 0.5*(y(1) - y(3))/den; end
end

function [dz, edgeHit, rpk] = match_z(stack, zs, im)
% Which z-plane does this frame look like?  Correlation against every slice,
% parabolic interpolation on the peak for sub-step resolution.
%
% edgeHit flags the case that would otherwise pass silently: the best-matching
% plane is the FIRST or LAST slice, meaning the true focus lies outside the
% stack. interp1 then extrapolates and returns a plausible-looking number that
% is really just a lower bound on the shift.
%
% rpk is the peak correlation. A low value means the match itself is weak --
% worth knowing before trusting the dz, e.g. if the structure is a symmetric
% blur (a bead or a smooth slide), which is degenerate about focus.
n = numel(zs);  r = zeros(1,n);
a = double(im(:)); a = a - mean(a);
for k = 1:n
    b = double(reshape(stack(:,:,k), [], 1)); b = b - mean(b);
    r(k) = (a.'*b) / (norm(a)*norm(b) + eps);
end
[rpk, kmax] = max(r);
edgeHit = (kmax == 1) || (kmax == n);
k = kmax;
if k > 1 && k < n
    den = r(k-1) - 2*r(k) + r(k+1);
    if den ~= 0, k = k + 0.5*(r(k-1) - r(k+1))/den; end
end
dz = interp1(1:n, zs, k, 'linear', 'extrap');
end

function s = robust_slope(x, y)
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 3, s = NaN; return; end
p = polyfit(x(ok), y(ok), 1);
s = p(1);
end

function plot_fits(fits, names, G, outDir)
% nM is NOT in scope here -- plot_fits is its own function. Derive the count
% from names, which is already trimmed to MODES by the caller.
nM = numel(names);
hf = figure('Color','w','Position',[100 100 380*nM+60 340]);
for m = 1:nM
    ax = subplot(1,3,m); hold(ax,'on');
    f = fits.(names{m});
    plot(ax, f.c, f.dx, 'o-', 'DisplayName','dx');
    plot(ax, f.c, f.dy, 's-', 'DisplayName','dy');
    plot(ax, f.c, f.dz, '^-', 'DisplayName','dz');
    xlabel(ax,'coefficient (um)'); ylabel(ax,'sample motion (um)');
    title(ax, sprintf('%s   [%.2f %.2f %.2f]', names{m}, G(m,1), G(m,2), G(m,3)));
    legend(ax,'Location','best'); box(ax,'off'); grid(ax,'on');
end
% exportgraphics only exists from R2020a and the ScanImage MATLAB is older --
% the 1821 run died here with "Undefined function or variable 'exportgraphics'",
% after the save, so nothing was lost but the figure.
png = fullfile(outDir,'dm_motion_calib.png');
if exist('exportgraphics','file') == 2 || exist('exportgraphics','builtin') == 5
    exportgraphics(hf, png, 'Resolution', 150);
else
    set(hf,'PaperPositionMode','auto','InvertHardcopy','off');
    print(hf, png, '-dpng', '-r150');
end
fprintf('figure -> %s\n', png);
end
