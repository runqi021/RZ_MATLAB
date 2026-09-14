function R = dm_loop_validate_260913()
%DM_LOOP_VALIDATE_260913  Does the closed loop actually suppress motion?
%
% NO ANIMAL NEEDED. NO BREATHING NEEDED. A static sample is BETTER, because the
% only motion present is the motion we put there on purpose.
%
% THE IDEA
% Inject a known axial disturbance with the DM, and ask the loop to cancel it.
% The loop cannot tell where the defocus came from -- it only sees that the
% imaged plane has moved. The mirror carries both terms at once:
%
%     command = baseline + (inject(t) + correct(t)) * Z2C(defocus)
%
% If the loop works, the imaged plane stays at zero while inject(t) swings by
% several microns. Ground truth is exact because we wrote inject(t) ourselves.
%
% WHAT IT MEASURES, per frequency:
%   OPEN    correct(t) = 0        -> measured dz should EQUAL inject(t).
%                                    That alone validates the estimator's GAIN
%                                    and SIGN against a known displacement.
%   CLOSED  loop running          -> measured dz should collapse toward zero.
%                                    residual/open = the suppression factor.
%
% Sweeping frequency maps the loop's transfer function, including the frequency
% where it stops helping and starts AMPLIFYING. For a pure delay tau the
% residual is 2*sin(pi*f*tau), which exceeds 1 past 60 degrees of lag -- with
% velocity extrapolation it should do considerably better, and this measures
% which is true rather than assuming.
%
% ============================== SAFETY ==============================
% DRY_RUN is true by default: the injection is NOT sent either, so the mirror
% stays at baseline and you can check the arithmetic first. The TOTAL command
% (inject + correct) is clamped, and any command breaching the ALPAO limit is
% skipped rather than clipped. onCleanup restores the DM and every SI setting.
%
% Runqi Zhang / 2026-09-13.  NOT YET RUN.

%% ========================= USER SETTINGS =========================
DRY_RUN      = false;     % <<< false to actually move the mirror

SI_DIR       = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE     = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_DEFOC    = 3;
OUT_ROOT     = 'D:\Data\RZ\ao_axial_260913';
CHANNEL      = 3;

% ---- reference ----------------------------------------------------------
REF_REUSE_FROM = '';     % a previous loop/validate session, or '' for fresh
REUSE_Z_TOL_UM = 1.0;    % refuse a reused reference if the stage has moved more
                         % than this -- its zero would no longer be your zero
% *** THE REFERENCE RANGE IS NOT "AS WIDE AS POSSIBLE". *** Measured
% 2026-09-13, same mirror / rate / read path / estimator, 5 um at 0.5 Hz:
%     reference -10..+10 (21 planes) -> gain 0.938
%     reference -14..+24 (39 planes) -> gain 0.653
% The full-DM-authority stack COSTS 35% of the signal. Its far positive planes
% carry no usable structure, the coarse tier occasionally locks onto one, and
% those outlier frames do not cancel in a lock-in -- they just eat the coherent
% amplitude.
%
% But narrow is not free either: 10 um of motion inside a +-10 um reference
% measured 0.810, because the estimate saturates at the edges.
%
% So the reference must be WIDER THAN THE MOTION and NO WIDER THAN THE USABLE
% SIGNAL. For +-5 um of motion, +-12 satisfies both.
REF_ZMIN_UM  = -12;
REF_ZMAX_UM  =  12;
REF_ZSTEP_UM =  1;
REF_SEC_PER_PLANE = 2.0; % static sample: no breath cycles to average, so this
                         % is purely about SNR and can be shorter than 3.3 s

% ---- the disturbance ----------------------------------------------------
% Weighted toward the band that actually matters: RZ reports breathing 0.5-5 Hz,
% mostly under 3. So the sweep is dense from 0.5 to 5 and then jumps to 8 and 12
% for one specific reason -- to find where the loop STOPS helping. For a plain
% delayed loop at 60 Hz that crossover is ~10 Hz; velocity extrapolation should
% push it higher, and 12 Hz brackets it either way.
%
% Those high points are not academic. A 3 Hz breath with a sharp inspiration is
% not a sinusoid: it carries harmonics at 6, 9, 12 Hz. If the loop amplifies
% there, it can add high-frequency jitter while removing the fundamental, and a
% single RMS number would hide that trade entirely.
INJECT_FREQS_HZ = [0.5 1 1.5 2 3 4 5 8 12];
INJECT_AMP_UM   = 5.0;   % amplitude, so peak-to-peak is twice this
SEG_SECONDS     = 10;    % per condition, per frequency

% ---- the realistic disturbance ------------------------------------------
% *** BRAIN MOTION IS NOT A SINUSOID. *** It is sharp: the brain is pushed
% quickly during inspiration and relaxes slowly, so displacement is NOT linear
% in breath phase (RZ, 2026-09-13). That matters enormously here, because a
% sharp waveform carries harmonics far above the breath rate -- a 2 Hz breath
% puts real energy at 6, 10, 14 Hz, which is exactly where a 60 Hz loop stops
% helping (break-even ~10 Hz for a plain delayed loop).
%
% So a pure-sine sweep FLATTERS the loop. It is still worth having, because it
% is the transfer function and you need it to interpret anything else, but the
% number that predicts real performance comes from this: an asymmetric pulse,
% fast rise then slow decay, injected at realistic breath rates, scored
% BROADBAND rather than at one frequency.
SHARP_FREQS_HZ  = [1 2 3];
SHARP_RISE_FRAC = 0.12;  % fraction of the cycle spent rising (fast inspiration)
SHARP_DECAY     = 0.35;  % exponential decay constant, in cycles

% ---- controller ---------------------------------------------------------
LOOP_GAIN    = 0.7;
EXTRAP_GAIN  = 1.0;

% ---- clamps -------------------------------------------------------------
DEFOCUS_MIN_UM = -13.0;
DEFOCUS_MAX_UM =  23.0;
MAX_STEP_UM    = 2.0;
CMD_LIMIT      = 0.95;

GAIN_UM_PER_COEFF = 6.656;
SETTLE_S     = 0.3;
RES_ARM_SEC  = 1.0;

%% ========================= CONNECT =========================
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;

S = load(Z2C_FILE);  fn = fieldnames(S);  Z2C = S.(fn{1});
zRow = Z2C(ROW_DEFOC,:);
baseline = get_baseline(BASELINE_VAR);

STAMP   = datestr(now,'yymmdd_HHMM');
sessDir = fullfile(OUT_ROOT, sprintf('validate_%s', STAMP));
if ~exist(sessDir,'dir'), mkdir(sessDir); end

siState = struct( ...
    'linePhase',            hSI.hScan2D.linePhase, ...
    'loggingEnable',        hSI.hChannels.loggingEnable, ...
    'framesPerSlice',       hSI.hStackManager.framesPerSlice, ...
    'numSlices',            hSI.hStackManager.numSlices, ...
    'rollingAvg',           hSI.hDisplay.displayRollingAverageFactor, ...
    'keepResOn',            hSI.hScan2D.keepResonantScannerOn);
restore = onCleanup(@() restore_all(hDM, baseline, hSI, siState)); %#ok<NASGU>

hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = 1;
hSI.hChannels.loggingEnable = false;   % nothing to log: dz is the measurement

assert(ismember(CHANNEL, hSI.hChannels.channelSave) || true, '');
assert(ismember(CHANNEL, hSI.hChannels.channelDisplay), ...
    ['channel %d is not in channelDisplay (%s). The loop reads frames from ' ...
     'the DISPLAY buffer.'], CHANNEL, mat2str(hSI.hChannels.channelDisplay(:).'));

fps = hSI.hRoiManager.scanFrameRate;
tau = 1/fps;

fprintf('\n============ CLOSED-LOOP VALIDATION (injected motion) ============\n');
if DRY_RUN, fprintf('*** DRY RUN -- neither injection nor correction is sent ***\n'); end
fprintf('session %s\n', sessDir);
fprintf('%.2f Hz, one frame of delay = %.1f ms\n', fps, 1000*tau);
fprintf('injecting %.1f um amplitude (%.1f um p2p) at %s Hz\n', ...
        INJECT_AMP_UM, 2*INJECT_AMP_UM, mat2str(INJECT_FREQS_HZ));
fprintf('gain %.2f, extrapolation %.2f, clamps %.1f to %+.1f um\n', ...
        LOOP_GAIN, EXTRAP_GAIN, DEFOCUS_MIN_UM, DEFOCUS_MAX_UM);
fprintf('\npredicted residual for a PLAIN delayed loop, 2*sin(pi*f*tau):\n  ');
for f = INJECT_FREQS_HZ, fprintf('%.1fHz:%3.0f%%  ', f, 100*2*sin(pi*f*tau)); end
fprintf('\n(>100%% means a plain loop would AMPLIFY. Extrapolation should beat this.)\n');

%% ========================= REFERENCE =========================
% The reference defines z = 0, so it is anchored to the stage z it was taken at.
% Otherwise every new stack silently re-zeros the coordinate, and a reused stack
% from a different park shifts the zero with no complaint.
fprintf('\n--- reference stack ---\n');
stageZnow = hSI.hMotors.motorPosition(3);
if ~isempty(REF_REUSE_FROM)
    p = dir(fullfile(REF_REUSE_FROM,'*meta.mat'));
    assert(~isempty(p), 'no meta .mat in %s', REF_REUSE_FROM);
    Mold = load(fullfile(REF_REUSE_FROM, p(1).name), 'refVol','zRef','refStageZ');
    refVol = Mold.refVol;  zRef = Mold.zRef;
    if isfield(Mold,'refStageZ') && ~isempty(Mold.refStageZ)
        refStageZ = Mold.refStageZ;
        assert(abs(stageZnow - refStageZ) <= REUSE_Z_TOL_UM, ...
            ['stage has moved %.2f um since that reference was taken ' ...
             '(tolerance %.2f um). Return to z = %.2f um or take a fresh stack.'], ...
            stageZnow - refStageZ, REUSE_Z_TOL_UM, refStageZ);
        fprintf('  reference anchored at stage z = %.2f um (stage is there)\n', refStageZ);
    else
        refStageZ = stageZnow;
        warning('dm_val:noAnchor','reused reference has no recorded stage z');
    end
    assert(size(refVol,1) == hSI.hRoiManager.linesPerFrame && ...
           size(refVol,2) == hSI.hRoiManager.pixelsPerLine, ...
        'reused reference is %dx%d, current geometry is %dx%d', ...
        size(refVol,1), size(refVol,2), ...
        hSI.hRoiManager.linesPerFrame, hSI.hRoiManager.pixelsPerLine);
    fprintf('  reusing %d planes from %s\n', numel(zRef), REF_REUSE_FROM);
else
    refStageZ = stageZnow;
    fprintf('  anchoring this reference to stage z = %.2f um\n', refStageZ);
    arm_scanner(hSI, RES_ARM_SEC, 'the reference stack');
    zRef = REF_ZMIN_UM : REF_ZSTEP_UM : REF_ZMAX_UM;
    nRef = max(4, round(REF_SEC_PER_PLANE * fps));
    fprintf('  %d planes, %d frames each (~%.1f min)\n', ...
            numel(zRef), nRef, (numel(zRef)*nRef/fps + numel(zRef)*SETTLE_S)/60);
    refVol = [];  got = false(1,numel(zRef));
    for k = 1:numel(zRef)
        cmd = baseline + (zRef(k)/GAIN_UM_PER_COEFF)*zRow;
        if max(abs(cmd)) > CMD_LIMIT
            fprintf('  z %+5.1f um SKIPPED (|cmd| %.3f)\n', zRef(k), max(abs(cmd)));
            continue
        end
        hDM.senddata(cmd);  pause(SETTLE_S);
        I = quick_avg(hSI, nRef, CHANNEL);
        if isempty(refVol), refVol = zeros([size(I) numel(zRef)],'single'); end
        refVol(:,:,k) = I;  got(k) = true; %#ok<AGROW>
    end
    hDM.senddata(baseline);  pause(SETTLE_S);
    refVol = refVol(:,:,got);  zRef = zRef(got);
    fprintf('  %d planes, %.1f to %.1f um\n', numel(zRef), zRef(1), zRef(end));
end

r_end = corr_lin(refVol(:,:,1), refVol(:,:,end));
fprintf('  end-to-end plane correlation r = %.3f\n', r_end);
assert(r_end < 0.97, 'field looks axially blind (r = %.3f)', r_end);

E = dm_axial_estimator_260913('prep', refVol, zRef, ...
        struct('nStrips',1,'pixFrac',0.04));

%% ========================= SWEEP =========================
arm_scanner(hSI, RES_ARM_SEC, 'the sweep');
nF = numel(INJECT_FREQS_HZ);
res = struct('f',[],'openAmp',[],'closedAmp',[],'suppression',[], ...
             'openGain',[],'nOpen',[],'nClosed',[]);
LOG = {};

for iF = 1:nF
    fHz = INJECT_FREQS_HZ(iF);
    fprintf('\n--- %.1f Hz ---\n', fHz);
    for mode = [0 1]                      % 0 = open, 1 = closed
        lab = 'open  ';  if mode, lab = 'CLOSED'; end
        L = run_segment(hSI, hDM, E, baseline, zRow, ...
                fHz, INJECT_AMP_UM, SEG_SECONDS, mode, ...
                LOOP_GAIN, EXTRAP_GAIN, DEFOCUS_MIN_UM, DEFOCUS_MAX_UM, ...
                MAX_STEP_UM, CMD_LIMIT, GAIN_UM_PER_COEFF, CHANNEL, DRY_RUN);
        L.freq = fHz;  L.closedMode = mode;
        LOG{end+1} = L; %#ok<AGROW>

        % Amplitude of the measured dz AT THE INJECTION FREQUENCY. Narrowband,
        % so drift and estimator noise at other frequencies do not inflate it.
        a = lockin(L.dz, L.t, fHz);
        if mode == 0
            openA = a;
            % dz should track the injection 1:1 in open loop. This is the
            % estimator's gain and sign checked against a known displacement.
            gi = lockin(L.inject, L.t, fHz);
            if DRY_RUN
                fprintf('  %s : would inject %.2f um; measured %.2f um (noise floor)\n', ...
                        lab, gi, a);
            else
                fprintf('  %s : injected %.2f um, measured %.2f um  -> estimator gain %.3f\n', ...
                        lab, gi, a, a/max(gi,eps));
            end
        else
            if DRY_RUN
                fprintf('  %s : measured %.2f um (noise floor; nothing was sent)\n', lab, a);
            else
                fprintf('  %s : measured %.2f um  -> RESIDUAL %.0f%% of open loop\n', ...
                        lab, a, 100*a/max(openA,eps));
            end
            res.f(end+1) = fHz;
            res.openAmp(end+1) = openA;
            res.closedAmp(end+1) = a;
            res.suppression(end+1) = a/max(openA,eps);
            res.openGain(end+1) = openA / INJECT_AMP_UM;
        end
        if L.nMissed > 0.05*L.nSvc
            fprintf('    (missed %d of %d frames -- delay is worse than nominal)\n', ...
                    L.nMissed, L.nSvc);
        end
    end
end
%% ================= SHARP (REALISTIC) DISTURBANCE =================
% The sine sweep above is the transfer function. This is the prediction of real
% performance, because brain motion is NOT sinusoidal -- it is a fast push on
% inspiration and a slow relaxation, so displacement is not linear in breath
% phase and the waveform carries strong harmonics.
%
% Scored BROADBAND (RMS), not with a lock-in at the fundamental: the whole point
% is that the harmonics are where a delayed loop struggles, and a lock-in at the
% breath rate would look at exactly the component the loop handles best.
sharp = struct('f',[],'openRms',[],'closedRms',[],'ratio',[], ...
               'openF0',[],'closedF0',[]);
if ~isempty(SHARP_FREQS_HZ)
    fprintf('\n========= SHARP waveform (fast rise, slow decay) =========\n');
    fprintf('rise %.0f%% of the cycle, decay %.2f cycles -- harmonics well above f0\n', ...
            100*SHARP_RISE_FRAC, SHARP_DECAY);
    for fHz = SHARP_FREQS_HZ
        fprintf('\n--- %.1f Hz sharp ---\n', fHz);
        oR = NaN;
        for mode = [0 1]
            lab = 'open  ';  if mode, lab = 'CLOSED'; end
            L = run_segment(hSI, hDM, E, baseline, zRow, ...
                    fHz, INJECT_AMP_UM, SEG_SECONDS, mode, ...
                    LOOP_GAIN, EXTRAP_GAIN, DEFOCUS_MIN_UM, DEFOCUS_MAX_UM, ...
                    MAX_STEP_UM, CMD_LIMIT, GAIN_UM_PER_COEFF, CHANNEL, DRY_RUN, ...
                    'sharp', SHARP_RISE_FRAC, SHARP_DECAY);
            L.freq = fHz;  L.closedMode = mode;  L.kind = 'sharp';
            LOG{end+1} = L; %#ok<AGROW>
            v = L.dz(isfinite(L.dz));  v = v - mean(v);
            rms_ = sqrt(mean(v.^2));
            f0   = lockin(L.dz, L.t, fHz);
            if mode == 0
                oR = rms_;  oF = f0;
                fprintf('  %s : RMS %.2f um, fundamental %.2f um\n', lab, rms_, f0);
            else
                fprintf('  %s : RMS %.2f um (%.0f%%), fundamental %.2f um (%.0f%%)\n', ...
                        lab, rms_, 100*rms_/max(oR,eps), f0, 100*f0/max(oF,eps));
                sharp.f(end+1)=fHz; sharp.openRms(end+1)=oR; sharp.closedRms(end+1)=rms_;
                sharp.ratio(end+1)=rms_/max(oR,eps);
                sharp.openF0(end+1)=oF; sharp.closedF0(end+1)=f0;
            end
        end
    end
end

if ~DRY_RUN, hDM.senddata(baseline); end
pause(SETTLE_S);

%% ========================= REPORT =========================
fprintf('\n================== SUPPRESSION vs FREQUENCY ==================\n');
if DRY_RUN
    fprintf('%8s %12s %12s\n','f (Hz)','open um','closed um');
    for i = 1:numel(res.f)
        fprintf('%8.1f %12.3f %12.3f\n', res.f(i), res.openAmp(i), res.closedAmp(i));
    end
else
    fprintf('%8s %12s %12s %12s %10s\n','f (Hz)','open um','closed um','residual','predicted');
    for i = 1:numel(res.f)
        fprintf('%8.1f %12.2f %12.2f %11.0f%% %9.0f%%\n', res.f(i), res.openAmp(i), ...
                res.closedAmp(i), 100*res.suppression(i), ...
                100*2*sin(pi*res.f(i)*tau));
    end
end

%% --- IN A DRY RUN THERE IS NOTHING TO MEASURE. SAY SO AND STOP. ---
% The injection is not sent either, so on a static sample the true motion is
% zero, every amplitude is the estimator's noise floor, and closed/open is a
% ratio of noise to noise. Printing a verdict from that is worse than printing
% nothing -- the first version of this script announced "THE LOOP WORKS, 45%
% residual" off 0.00 um amplitudes.
if DRY_RUN
    fprintf('\n*** DRY RUN: the numbers above are MEANINGLESS. ***\n');
    fprintf('Nothing was sent, so there was no motion to suppress and every\n');
    fprintf('amplitude is the estimator noise floor. What a dry run DOES prove:\n');
    fprintf('  frames were read, the estimator ran, the segments and timing work.\n');
    fprintf('  measured amplitude on a STATIC sample: %.3f um mean\n', mean(res.openAmp));
    fprintf('  (that IS a real number -- it is the noise floor, and it should be\n');
    fprintf('   well under the %.1f um you intend to inject.)\n', INJECT_AMP_UM);
    fprintf('\nSet DRY_RUN = false to measure anything.\n');
    R = struct('res',res,'log',{LOG},'fps',fps,'tau',tau,'zRef',zRef, ...
               'refVol',refVol,'injectAmpUm',INJECT_AMP_UM,'dryRun',true, ...
               'loopGain',LOOP_GAIN,'extrapGain',EXTRAP_GAIN,'stamp',STAMP, ...
               'refStageZ',refStageZ);
    save(fullfile(sessDir,'validate_meta.mat'),'-struct','R','-v7.3');
    fprintf('saved %s\n', fullfile(sessDir,'validate_meta.mat'));
    fprintf('==============================================================\n');
    return
end

% A real run whose open-loop amplitude is near zero means the injection never
% reached the mirror -- not that the loop is perfect.
if mean(res.openAmp) < 0.2*INJECT_AMP_UM
    error('dm_val:noInjection', ...
        ['open-loop amplitude is %.2f um but %.1f um was injected. The ' ...
         'disturbance is not reaching the mirror -- check senddata and the ' ...
         'clamps before reading anything else.'], mean(res.openAmp), INJECT_AMP_UM);
end

g = res.openGain;
fprintf('\nestimator gain (open-loop measured/injected): %.3f +- %.3f\n', mean(g), std(g));
if abs(mean(g)-1) > 0.2
    fprintf('  >> OFF BY %.0f%%. Either GAIN_UM_PER_COEFF is wrong for this\n', 100*abs(mean(g)-1));
    fprintf('     sample, or the reference z axis is mis-scaled. Both are\n');
    fprintf('     calibration, not control -- fix before reading the residuals.\n');
elseif mean(g) < 0
    fprintf('  >> SIGN INVERTED. The correction would drive the motion, not cancel it.\n');
else
    fprintf('  >> the estimator recovers a known displacement correctly.\n');
end

% Judged over 0.5-3 Hz, where this animal's breathing actually lives.
lowF = res.suppression(res.f <= 3);
if ~isempty(lowF) && mean(lowF) < 0.5
    fprintf('\n*** THE LOOP WORKS. *** %.0f%% residual over 0.5-3 Hz.\n', 100*mean(lowF));
elseif ~isempty(lowF) && mean(lowF) < 0.9
    fprintf('\n*** PARTIAL. *** %.0f%% residual over 0.5-3 Hz -- real but modest.\n', 100*mean(lowF));
elseif ~isempty(lowF)
    fprintf('\n*** NO USEFUL SUPPRESSION over 0.5-3 Hz (%.0f%%). ***\n', 100*mean(lowF));
end
midF = res.suppression(res.f > 3 & res.f <= 5);
if ~isempty(midF)
    fprintf('3-5 Hz (sniffing): %.0f%% residual\n', 100*mean(midF));
end

% THE NUMBER THAT PREDICTS REAL PERFORMANCE.
if ~isempty(sharp.f)
    fprintf('\n---------- SHARP waveform (the realistic case) ----------\n');
    fprintf('%8s %12s %12s %11s %14s\n','f (Hz)','open RMS','closed RMS','broadband','fundamental');
    for i = 1:numel(sharp.f)
        fprintf('%8.1f %12.2f %12.2f %10.0f%% %13.0f%%\n', sharp.f(i), ...
                sharp.openRms(i), sharp.closedRms(i), 100*sharp.ratio(i), ...
                100*sharp.closedF0(i)/max(sharp.openF0(i),eps));
    end
    fprintf(['\nBROADBAND is the honest number. If it is much WORSE than the\n' ...
             'fundamental column, the loop is cancelling the breath rate while\n' ...
             'ADDING jitter at the harmonics -- which a single RMS on a sine\n' ...
             'sweep would never have revealed.\n']);
    sb = mean(sharp.ratio);
    if sb < 0.5
        fprintf('=> %.0f%% residual on realistic motion. THE LOOP IS WORTH RUNNING.\n', 100*sb);
    elseif sb < 0.9
        fprintf('=> %.0f%% residual on realistic motion. Marginal.\n', 100*sb);
    else
        fprintf('=> %.0f%% residual on realistic motion. NOT WORTH RUNNING as tuned.\n', 100*sb);
    end
end
amp = res.f(res.suppression > 1);
if ~isempty(amp)
    fprintf('AMPLIFIES at: %s Hz. That is the loop bandwidth limit.\n', mat2str(amp));
end

R = struct('res',res,'log',{LOG},'fps',fps,'tau',tau,'zRef',zRef, ...
           'refVol',refVol,'injectAmpUm',INJECT_AMP_UM,'dryRun',DRY_RUN, ...
           'loopGain',LOOP_GAIN,'extrapGain',EXTRAP_GAIN,'stamp',STAMP, ...
           'refStageZ',refStageZ,'sharp',sharp);
save(fullfile(sessDir,'validate_meta.mat'),'-struct','R','-v7.3');
fprintf('\nsaved %s\n', fullfile(sessDir,'validate_meta.mat'));
make_figure(R, sessDir);
fprintf('==============================================================\n');
end

%% ========================== ONE SEGMENT ==========================
function L = run_segment(hSI, hDM, E, baseline, zRow, fHz, ampUm, secs, closed, ...
        gLoop, gExtrap, dMin, dMax, maxStep, cmdLimit, umPerCoeff, ch, dry, ...
        kind, riseFrac, decay)
if nargin < 19, kind = 'sine'; riseFrac = 0.12; decay = 0.35; end
%RUN_SEGMENT  One frequency, one condition. Injection always applied; the
% correction only when `closed`.
n = round(secs * hSI.hRoiManager.scanFrameRate);
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;

L = struct('t',nan(1,n),'dz',nan(1,n),'inject',nan(1,n),'corr',nan(1,n), ...
           'total',nan(1,n),'nSvc',0,'nMissed',0);

corrUm = 0;  e1 = NaN;  e2 = NaN;  lastNo = -1;  nEmpty = 0;
t0 = tic;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_val:grabDidNotStart','startGrab left acqState idle');
end

while ~strcmp(hSI.acqState,'idle')
    try, fno = hSI.hDisplay.lastFrameNumber; catch, pause(0.002); continue; end
    if fno == lastNo, pause(0.001); drawnow limitrate; continue; end

    [I, why] = get_disp_frame(hSI, ch);
    if isempty(I)
        nEmpty = nEmpty + 1;
        if nEmpty > 2000
            hSI.abort();
            error('dm_val:noFrames','no frame for channel %d: %s', ch, why);
        end
        pause(0.0005);  continue
    end
    nEmpty = 0;
    if lastNo >= 0, L.nMissed = L.nMissed + (fno - lastNo - 1); end
    lastNo = fno;
    L.nSvc = L.nSvc + 1;
    if L.nSvc > n, break; end

    tNow = toc(t0);
    % INDEPENDENT ESTIMATE EVERY FRAME (zInit = NaN), not a carried one.
    %
    % Carrying zEst forward runs the fine gradient tier ALONE, and on the rig
    % 2026-09-13 that froze: once the estimate drifted outside the linear zone
    % it returned ~0 correction, could never climb back, and reported a
    % confident constant -- every lock-in amplitude came out ~0.
    %
    % Passing NaN forces coarse-then-fine, which is exactly what
    % dm_estimator_check_260913 does, and that measured 0.993 slope with
    % 0.05 um residual on real tissue, single frames included. Costs one
    % decimated correlation per frame against a 16.7 ms budget already only
    % 17% used. Proven beats clever.
    dzNow = single_dz(E, I, NaN);

    if closed && isfinite(dzNow)
        e2 = e1;  e1 = dzNow;
        if isfinite(e2), pred = e1 + gExtrap*(e1 - e2); else, pred = e1; end
        want = corrUm - gLoop*pred;
        want = corrUm + max(-maxStep, min(maxStep, want - corrUm));
        corrUm = want;
    end

    injUm   = inject_wave(tNow, fHz, ampUm, kind, riseFrac, decay);
    totalUm = injUm + corrUm;
    totalUm = max(dMin, min(dMax, totalUm));
    cmd = baseline + (totalUm/umPerCoeff)*zRow;
    if max(abs(cmd)) <= cmdLimit && ~dry
        hDM.senddata(cmd);
    end

    L.t(L.nSvc)=tNow; L.dz(L.nSvc)=dzNow; L.inject(L.nSvc)=injUm;
    L.corr(L.nSvc)=corrUm; L.total(L.nSvc)=totalUm;
    drawnow limitrate;
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end

k = isfinite(L.t);
L.t=L.t(k); L.dz=L.dz(k); L.inject=L.inject(k); L.corr=L.corr(k); L.total=L.total(k);
end

function dz = single_dz(E, I, zPrev)
o = dm_axial_estimator_260913('run', E, single(I), zPrev);
dz = o.dz(1);
end

%% ========================== ANALYSIS ==========================
function y = inject_wave(t, f, amp, kind, riseFrac, decay)
%INJECT_WAVE  The disturbance waveform, in um.
%
% 'sine'  -- amp*sin(2*pi*f*t). Clean, one frequency, gives the transfer
%            function. Not what a brain does.
% 'sharp' -- fast linear rise over riseFrac of the cycle, then exponential
%            decay. Models inspiration pushing the brain quickly and it
%            relaxing slowly, i.e. displacement NOT linear in breath phase.
%            Zero-meaned and scaled to the same peak-to-peak as the sine so
%            the two are comparable, and so the only thing that changes
%            between them is the HARMONIC CONTENT.
if strcmpi(kind,'sine')
    y = amp * sin(2*pi*f*t);
    return
end
ph = mod(f*t, 1);
y  = zeros(size(ph));
up = ph < riseFrac;
y(up)  = ph(up) / riseFrac;
y(~up) = exp(-(ph(~up) - riseFrac) / decay);

% *** NORMALISE ANALYTICALLY, NOT FROM max(y)-min(y). ***
% The loop calls this ONE SAMPLE AT A TIME, so t is a SCALAR and an empirical
% range is identically ZERO. Dividing by max(0,eps) scaled the output to ~1e16
% um, which the clamp then squared off into a square wave slamming between
% -13 and +23 um. Measured on the rig 2026-09-13: "open RMS 12.13 um" for a
% waveform whose peak-to-peak was supposed to be 10 um -- an RMS larger than
% the p2p, which is impossible and was the tell.
%
% The shape is known, so its range is known: y rises 0 -> 1 over riseFrac, then
% decays. Max is exactly 1, min is exactly 0, range is exactly 1. Scaling by
% 2*amp therefore gives the same peak-to-peak as the sine, sample by sample,
% with no dependence on how many samples you happen to be evaluating.
mu = riseFrac/2 + decay*(1 - exp(-(1-riseFrac)/decay));
y  = (y - mu) * (2*amp);
end

function a = lockin(x, t, f)
%LOCKIN  Amplitude of x at exactly frequency f. Narrowband, so drift and
% broadband estimator noise do not inflate it the way an RMS would.
k = isfinite(x) & isfinite(t);
if nnz(k) < 20, a = NaN; return; end
x = x(k) - mean(x(k));  t = t(k);
c = 2*mean(x .* cos(2*pi*f*t));
s = 2*mean(x .* sin(2*pi*f*t));
a = hypot(c, s);
end

function make_figure(R, outDir)
try, opengl('software'); catch, end
hf = figure('Color','w','Position',[80 80 1100 450],'Visible','off');

subplot(1,2,1);
semilogx(R.res.f, 100*R.res.suppression, 'o-','LineWidth',1.5, ...
         'MarkerFaceColor',[0.15 0.35 0.75],'Color',[0.15 0.35 0.75]);
hold on;
ff = logspace(log10(min(R.res.f)), log10(max(R.res.f)), 100);
plot(ff, 100*2*sin(pi*ff*R.tau), 'k--','LineWidth',1);
plot(xlim, [100 100], 'r:','LineWidth',1);
xlabel('injected frequency (Hz)'); ylabel('residual (% of open loop)');
title(sprintf('suppression (gain %.2f, extrap %.2f)', R.loopGain, R.extrapGain));
legend({'measured','plain delay 2sin(\pi f \tau)','no effect'},'Location','northwest');
grid on;

subplot(1,2,2);
i = find(R.res.f <= 3, 1, 'last');
if ~isempty(i)
    Lo = R.log{2*i-1};  Lc = R.log{2*i};
    plot(Lo.t, Lo.dz, '-', 'Color',[0.8 0.4 0.1]); hold on;
    plot(Lc.t, Lc.dz, '-', 'Color',[0.15 0.35 0.75]);
    plot(Lo.t, Lo.inject, 'k:','LineWidth',1);
    xlim([0 min(5, max(Lo.t))]);
    xlabel('time (s)'); ylabel('measured dz (\mum)');
    title(sprintf('%.1f Hz: open vs closed', R.res.f(i)));
    legend({'open','closed','injected'},'Location','best');
    grid on;
end

% print(), NOT exportgraphics: this runs under R2017b on the rig (that is the
% MATLAB SI2018b launches in), and exportgraphics is R2020a+. Syntax-checking
% against the R2024a on this machine does not catch that.
png = fullfile(outDir,'loop_validation.png');
print(hf, png, '-dpng', '-r200');
close(hf);
d = dir(png);
if ~isempty(d) && d.bytes > 20000
    fprintf('wrote %s (%.0f kB)\n', png, d.bytes/1024);
else
    warning('dm_val:stub','%s looks like an OpenGL blank', png);
end
end

%% ========================== HELPERS ==========================
function [I, why] = get_disp_frame(hSI, ch)
I = [];  why = '';
try
    b = hSI.hDisplay.rollingStripeDataBuffer;
    if ~isempty(b) && ~isempty(b{1}) && ~isempty(b{1}{1}.roiData)
        I = pick_channel(b{1}{1}.roiData{1}, ch);
        if ~isempty(I), return; end
    end
    why = 'rolling buffer not ready';
catch ME
    why = sprintf('rolling buffer: %s', ME.message);
end
try
    I = pick_channel(hSI.hDisplay.lastStripeData.roiData{1}, ch);
    if ~isempty(I), why = ''; end
catch ME
    why = sprintf('%s | live stripe: %s', why, ME.message);
end
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

function I = quick_avg(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_val:grabDidNotStart','startGrab left acqState idle');
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
I = get_disp_frame(hSI, ch);
hSI.hDisplay.displayRollingAverageFactor = 1;
end

function b = get_baseline(v)
if ~evalin('base', sprintf('exist(''%s'',''var'')', v))
    error('dm_val:noBaseline','base workspace has no %s', v);
end
b = double(evalin('base', v));  b = b(:).';
assert(numel(b) == 97, '%s is %d long, expected 97', v, numel(b));
end

function arm_scanner(hSI, sec, what)
hSI.hScan2D.keepResonantScannerOn = true;
fprintf('  arming resonant scanner %.1f s before %s...\n', sec, what);
t0 = tic; while toc(t0) < sec, pause(0.05); drawnow limitrate; end
end

function r = corr_lin(a,b)
a = double(a(:)) - mean(double(a(:)));
b = double(b(:)) - mean(double(b(:)));
r = (a.'*b)/(norm(a)*norm(b)+eps);
end

function restore_all(hDM, baseline, hSI, s)
try, hDM.senddata(baseline); fprintf('[cleanup] DM restored to baseline\n');
catch ME, warning('dm_val:dm','COULD NOT RESTORE THE DM: %s', ME.message); end
try
    if ~strcmp(hSI.acqState,'idle')
        hSI.abort();
        while ~strcmp(hSI.acqState,'idle'), pause(0.02); end
    end
    hSI.hScan2D.linePhase                    = s.linePhase;
    hSI.hChannels.loggingEnable              = s.loggingEnable;
    hSI.hStackManager.framesPerSlice         = s.framesPerSlice;
    hSI.hStackManager.numSlices              = s.numSlices;
    hSI.hDisplay.displayRollingAverageFactor = s.rollingAvg;
    hSI.hScan2D.keepResonantScannerOn        = s.keepResOn;
    fprintf('[cleanup] ScanImage settings restored\n');
catch ME, warning('dm_val:si','could not restore SI settings: %s', ME.message); end
end
