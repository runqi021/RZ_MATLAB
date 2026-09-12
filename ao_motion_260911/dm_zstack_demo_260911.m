function dm_zstack_demo_260911()
%DM_ZSTACK_DEMO_260911  Same z-stack twice: once with the STAGE, once with the DM.
%
% The demonstration: take a z-stack by moving the stage, then take the SAME
% stack without moving the stage at all -- only by commanding defocus on the
% deformable mirror. If the DM is really moving the focus, the two stacks show
% the same planes in the same order, and plane k of the DM stack matches plane k
% of the stage stack.
%
% That match is checked quantitatively at the end: every DM plane is correlated
% against the whole stage stack, and the best-matching stage z is reported. A
% working DM gives nominal-z vs matched-z on the identity line.
%
% THE STAGE DOES NOT MOVE DURING THE DM STACK. That is the whole point.
%
% ================================ SAFETY ================================
% Every DM command is baseline + defocus, checked against the ALPAO +-1 limit
% before being sent; a plane that would clip is SKIPPED, not clamped. The DM is
% returned to your baseline and the stage to its starting z by onCleanup, even
% on error or Ctrl-C, along with every ScanImage setting this touches.
%
% STROKE IS ASYMMETRIC. Measured 2026-09-11: usable coefficient runs about
% -2.3 to +4.3, i.e. roughly -15 to +28 um of focus at a gain of 6.656 um per
% um of coefficient. A +-25 um stack therefore CANNOT be reached on the negative
% side -- those planes are skipped and reported. That is a real property of the
% mirror, not a failure of this script.
%
% Runqi Zhang / 2026-09-11.  NOT YET RUN.

%% ========================= USER SETTINGS =========================
SI_DIR     = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE   = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_DEFOC  = 3;

% Range about the CURRENT focus, asymmetric on purpose. The DM's usable
% coefficient runs about -2.3 to +4.3, which at 6.656 um/coeff is roughly -15 to
% +28 um -- so -15 to +25 is the largest window BOTH methods can cover, and
% making them cover the same planes is what makes the comparison meaningful.
ZMIN_UM    = -15;
ZMAX_UM    =  25;
ZSTEP_UM   = 1;        % as requested. NB 1 um is 3.2 stepper counts, so the
                       % STAGE really steps 0.9375 um -- the encoder is read at
                       % every plane so the true z is recorded regardless.
NFRAMES    = 15;
CHANNEL    = 3;

% um of focus per um of Zernike coefficient. From dm_motion_calib_260911_1824.
GAIN_UM_PER_COEFF = 6.656;

SETTLE_S   = 0.3;
RES_ARM_SEC = 1.0;
CMD_LIMIT  = 0.95;
ENC_UM_PER_COUNT = -0.09987;   % signed: encoder counts DOWN as SI z goes up

OUT_DIR    = fileparts(mfilename('fullpath'));
STAMP      = datestr(now,'yymmdd_HHMM');

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

siState = struct( ...
    'loggingEnable',  hSI.hChannels.loggingEnable, ...
    'framesPerSlice', hSI.hStackManager.framesPerSlice, ...
    'numSlices',      hSI.hStackManager.numSlices, ...
    'rollingAvg',     hSI.hDisplay.displayRollingAverageFactor, ...
    'keepResOn',      hSI.hScan2D.keepResonantScannerOn);
hSI.hChannels.loggingEnable = false;
restore = onCleanup(@() restore_all(hDM, baseline, hSI, siState, p0)); %#ok<NASGU>

hSI.hScan2D.keepResonantScannerOn = true;
fprintf('arming resonant scanner (%.1f s), held on for the run...\n', RES_ARM_SEC);
wait_s(RES_ARM_SEC);

zList = (ZMIN_UM : ZSTEP_UM : ZMAX_UM);

% TRIM THE REQUEST TO WHAT THE DM CAN ACTUALLY REACH.
% ZMIN/ZMAX above are the stroke computed against a command limit of 1.0, but the
% guard below is CMD_LIMIT = 0.95, so the two end planes sit in the gap between
% them: on 2026-09-11 they needed 0.9602 and 0.9589, over by 0.010 and 0.009.
% They were skipped -- correctly -- and the pre-allocated zeros went into the
% stack and the TIFF as two black planes. Decide the range here instead, so
% nothing is skipped and there is nothing to leave behind.
cAll   = zList / GAIN_UM_PER_COEFF;
mxAll  = arrayfun(@(cc) max(abs(baseline(:).' + cc*Z2C(ROW_DEFOC,:))), cAll);
reach  = mxAll <= CMD_LIMIT;
if ~all(reach)
    fprintf(['  trimming %d plane(s) the DM cannot reach at CMD_LIMIT %.2f: ' ...
             'z %s um (needed max|cmd| up to %.4f)\n'], nnz(~reach), CMD_LIMIT, ...
             mat2str(zList(~reach)), max(mxAll(~reach)));
    zList = zList(reach);
    fprintf('  z range is now %+.1f to %+.1f um, %d planes\n', ...
            zList(1), zList(end), numel(zList));
end
nZ    = numel(zList);
fprintf('baseline %s, max|cmd| %.3f | %d planes, %g um steps, gain %.3f um/coeff\n', ...
        BASELINE_VAR, max(abs(baseline)), nZ, ZSTEP_UM, GAIN_UM_PER_COEFF);

%% ========================= 1. STAGE STACK =========================
fprintf('\n===== STAGE stack: moving the motor, DM held at baseline =====\n');
hDM.senddata(baseline); pause(SETTLE_S);

% Run-up below the first plane so the sweep is monotonic and the first plane is
% not the one that absorbs the 0.205 um reversal backlash.
goto_z(hSI, p0, z0 + zList(1) - 5*ZSTEP_UM);
wait_s(0.3);

stackStage = [];  zEnc = nan(1,nZ);
for k = 1:nZ
    goto_z(hSI, p0, z0 + zList(k));
    pause(0.2);
    zEnc(k) = read_z_counts(hLSC);
    I = grab_avg(hSI, NFRAMES, CHANNEL);
    if isempty(stackStage), stackStage = zeros([size(I) nZ],'single'); end
    stackStage(:,:,k) = I; %#ok<AGROW>
    if mod(k,10) == 0 || k == nZ, fprintf('  plane %d/%d\n', k, nZ); end
end
goto_z(hSI, p0, z0); wait_s(0.3);

if all(isfinite(zEnc)) && range(zEnc) > 0
    zStageMeas = (zEnc - zEnc(1)) * ENC_UM_PER_COUNT + zList(1);
    fprintf('  encoder: %.3f um/step measured vs %g commanded\n', ...
            (zStageMeas(end)-zStageMeas(1))/(nZ-1), ZSTEP_UM);
else
    zStageMeas = zList;
    warning('dm_demo:noEnc','z encoder unreadable; recording COMMANDED z');
end

%% ========================= 2. DM STACK =========================
fprintf('\n===== DM stack: STAGE DOES NOT MOVE, defocus only =====\n');
zStageNow = hSI.hMotors.motorPosition(3);
fprintf('  stage parked at z = %.3f um and left there\n', zStageNow);

stackDM = zeros(size(stackStage),'single');
cUsed   = nan(1,nZ);  gotDM = false(1,nZ);
for k = 1:nZ
    c   = zList(k) / GAIN_UM_PER_COEFF;
    cmd = baseline(:).' + c * Z2C(ROW_DEFOC,:);
    if max(abs(cmd)) > CMD_LIMIT
        fprintf('  z %+6.1f um (c %+5.2f): SKIPPED, would reach |cmd| %.3f\n', ...
                zList(k), c, max(abs(cmd)));
        continue
    end
    hDM.senddata(cmd);
    pause(SETTLE_S);
    stackDM(:,:,k) = grab_avg(hSI, NFRAMES, CHANNEL);
    cUsed(k) = c;  gotDM(k) = true;
    if mod(k,10) == 0 || k == nZ, fprintf('  plane %d/%d (c %+5.2f)\n', k, nZ, c); end
end
hDM.senddata(baseline); pause(SETTLE_S);

zAfter = hSI.hMotors.motorPosition(3);
fprintf('  stage z before %.3f, after %.3f (drift %+.3f um)\n', ...
        zStageNow, zAfter, zAfter - zStageNow);

%% ========================= 3. DOES IT MATCH? =========================
fprintf('\n===== DM plane -> best-matching STAGE plane =====\n');
fprintf('%10s %10s %12s %10s\n','nominal z','coeff','matched z','peak r');
matched = nan(1,nZ);  rpk = nan(1,nZ);
for k = 1:nZ
    if ~gotDM(k), continue; end
    [matched(k), rpk(k)] = best_plane(stackStage, zStageMeas, stackDM(:,:,k));
    fprintf('%10.1f %10.2f %12.2f %10.3f\n', zList(k), cUsed(k), matched(k), rpk(k));
end

ok = gotDM & isfinite(matched);
if nnz(ok) >= 3
    p = polyfit(zList(ok), matched(ok), 1);
    fprintf(['\nmatched_z = %.3f * nominal_z %+.2f   (slope 1 and offset 0 = the ' ...
             'DM\nreproduces the stage stack exactly)\n'], p(1), p(2));
    fprintf('residual RMS %.2f um over %d planes, median peak r %.3f\n', ...
            sqrt(mean((matched(ok) - polyval(p,zList(ok))).^2)), nnz(ok), median(rpk(ok)));
    fprintf('planes the DM could NOT reach: %d of %d\n', nnz(~gotDM), nZ);
end

%% ========================= 4. SAVE =========================
matOut = fullfile(OUT_DIR, sprintf('dm_zstack_demo_%s.mat', STAMP));
save(matOut, 'stackStage','stackDM','zList','zStageMeas','zEnc','cUsed','gotDM', ...
     'matched','rpk','baseline','BASELINE_VAR','GAIN_UM_PER_COEFF','ZSTEP_UM', ...
     'ZMIN_UM','ZMAX_UM','NFRAMES','CHANNEL','z0','-v7.3');
fprintf('\nsaved %s\n', matOut);

% A plane that was never acquired is still the pre-allocated ZEROS. Writing it
% would put a black slice in the middle of a z-stack, which is worse than a
% missing one: anything that later registers into this stack can match it.
% The trim above should mean this never fires, but a send can still fail.
if ~all(gotDM)
    warning('dm_zstack:skipped', ...
        '%d plane(s) never acquired -- dropping them from the saved stack and TIFF: z %s', ...
        nnz(~gotDM), mat2str(zList(~gotDM)));
end
write_tif(fullfile(OUT_DIR, sprintf('demo_STAGE_zstack_%s.tif', STAMP)), stackStage(:,:,gotDM));
write_tif(fullfile(OUT_DIR, sprintf('demo_DM_zstack_%s.tif',    STAMP)), stackDM(:,:,gotDM));
fprintf(['wrote both stacks as TIFFs -- open them side by side and step through;\n' ...
         'they should show the same structure entering and leaving focus.\n']);
end

%% ========================== HELPERS ==========================
function [zBest, r] = best_plane(stack, zAxis, im)
% Correlate one frame against every plane of the reference stack; return the z
% of the best match, parabolically interpolated, and the peak correlation.
n = size(stack,3);  rr = zeros(1,n);
a = double(im(:));  a = a - mean(a);
na = norm(a);
for k = 1:n
    b = double(reshape(stack(:,:,k),[],1));  b = b - mean(b);
    rr(k) = (a.'*b) / (na*norm(b) + eps);
end
[r, k] = max(rr);
if k > 1 && k < n
    den = rr(k-1) - 2*rr(k) + rr(k+1);
    if den ~= 0, k = k + 0.5*(rr(k-1) - rr(k+1))/den; end
end
zBest = interp1(1:n, zAxis, k, 'linear', 'extrap');
end

function b = get_baseline(varName)
if ~evalin('base', sprintf('exist(''%s'',''var'')', varName))
    error('dm_demo:noBaseline', ...
        ['base workspace has no ''%s''. Load and activate your system ' ...
         'aberration correction first.'], varName);
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

function I = grab_avg(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
f = hSI.hDisplay.lastFrame;
chans = hSI.hDisplay.lastFrameChannels;
idx = find(chans == ch, 1);
if isempty(idx)
    error('dm_demo:noChannel','channel %d not displayed (available %s)', ...
          ch, mat2str(chans));
end
I = single(f{idx});
end

function wait_s(sec)
t0 = tic;  while toc(t0) < sec, pause(0.05); drawnow limitrate; end
end

function write_tif(path, stack)
% uint16 multipage TIFF for viewing. The .mat holds the single-precision
% original; this is scaled to the stack's own full range for display, so DO NOT
% take intensities off the TIFF.
s = stack;  s(~isfinite(s)) = 0;
lo = min(s(:));  hi = max(s(:));
if hi <= lo, hi = lo + 1; end
u = uint16( (s - lo) / (hi - lo) * 65535 );
for k = 1:size(u,3)
    if k == 1, imwrite(u(:,:,k), path);
    else,      imwrite(u(:,:,k), path, 'WriteMode','append');
    end
end
end

function restore_all(hDM, baseline, hSI, s, p0)
try, hDM.senddata(baseline); fprintf('[cleanup] DM restored to baseline\n');
catch ME, warning('dm_demo:noDM','COULD NOT RESTORE THE DM: %s', ME.message); end
try
    hSI.hMotors.moveStartRelative(p0);
    hSI.hMotors.moveWaitForFinish([1 2 3]);
    fprintf('[cleanup] stage returned to z %.3f\n', p0(3));
catch ME, warning('dm_demo:noStage','could not return the stage: %s', ME.message); end
try
    hSI.hChannels.loggingEnable              = s.loggingEnable;
    hSI.hStackManager.framesPerSlice         = s.framesPerSlice;
    hSI.hStackManager.numSlices              = s.numSlices;
    hSI.hDisplay.displayRollingAverageFactor = s.rollingAvg;
    hSI.hScan2D.keepResonantScannerOn        = s.keepResOn;
    fprintf('[cleanup] ScanImage settings restored\n');
catch ME, warning('dm_demo:noSI','could not restore SI settings: %s', ME.message); end
end
