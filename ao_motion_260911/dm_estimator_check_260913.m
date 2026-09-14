function R = dm_estimator_check_260913()
%DM_ESTIMATOR_CHECK_260913  Does the estimator read a KNOWN defocus on REAL data?
%
% The one question that has to be answered before any loop result means
% anything, and it needs no loop, no timing, no injection waveform and no
% animal. ~40 s.
%
% Steps the DM to known offsets, grabs a frame at each, and asks the estimator
% what it sees. Commanded vs estimated is the whole output.
%
% IT SEPARATES THREE FAILURES THAT LOOK IDENTICAL IN A CLOSED LOOP:
%
%   1. THE MIRROR IS NOT MOVING.
%      Reported as the frame-to-frame correlation at each step. If frames at
%      +-8 um are still ~1.0 correlated, senddata is not reaching the DM and
%      nothing downstream can work.
%
%   2. THE MIRROR MOVES BUT THE ESTIMATOR IS BLIND.
%      Frames change, estimate does not. Slope ~0. This is what the rig showed
%      on 2026-09-13: dz sat at a constant while a 5 um sinusoid was injected,
%      giving ~0 lock-in amplitude at every frequency.
%
%   3. IT WORKS ON AVERAGED FRAMES BUT NOT SINGLE ONES.
%      Then the problem is photon noise, not the algorithm, and the fix is
%      averaging / more pixels / a brighter plane -- not the controller. Both
%      are measured here, side by side, which is the point.
%
% Nothing is logged and no correction is applied. The DM returns to baseline.
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

% Deliberately OFF the reference grid, so an estimate that merely snaps to the
% nearest plane is visible as quantisation rather than passing as a good fit.
TEST_Z_UM    = [-6.4 -3.7 -1.3 0 1.6 3.2 5.8];
TEST_NAVG    = 30;       % frames for the averaged test
CMD_LIMIT    = 0.95;
GAIN_UM_PER_COEFF = 6.656;
SETTLE_S     = 0.3;

%% ---------------- connect ----------------
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;

S = load(Z2C_FILE);  fn = fieldnames(S);  Z2C = S.(fn{1});
zRow = Z2C(ROW_DEFOC,:);
if ~evalin('base', sprintf('exist(''%s'',''var'')', BASELINE_VAR))
    error('dm_est:noBaseline','base workspace has no %s', BASELINE_VAR);
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

assert(ismember(CHANNEL, hSI.hChannels.channelDisplay), ...
    'channel %d is not in channelDisplay (%s)', CHANNEL, ...
    mat2str(hSI.hChannels.channelDisplay(:).'));

lp = hSI.hScan2D.linePhase;
assert(abs(lp) < 20e-6, 'linePhase is %.1f us -- units error, frames are scrambled', lp*1e6);

fprintf('\n============== ESTIMATOR CHECK (known defocus) ==============\n');
fprintf('%.2f Hz, %d x %d, linePhase %.3f us, channel %d\n', ...
        hSI.hRoiManager.scanFrameRate, hSI.hRoiManager.pixelsPerLine, ...
        hSI.hRoiManager.linesPerFrame, lp*1e6, CHANNEL);

hSI.hScan2D.keepResonantScannerOn = true;
fprintf('arming resonant scanner 1.0 s...\n');
wait_s(1.0);

%% ---------------- reference ----------------
zRef = REF_ZMIN_UM : REF_ZSTEP_UM : REF_ZMAX_UM;
fprintf('\nreference: %d planes, %.1f to %.1f um, %d frames each\n', ...
        numel(zRef), zRef(1), zRef(end), REF_NFRAMES);
refVol = [];  got = false(1,numel(zRef));
for k = 1:numel(zRef)
    cmd = baseline + (zRef(k)/GAIN_UM_PER_COEFF)*zRow;
    if max(abs(cmd)) > CMD_LIMIT, continue; end
    hDM.senddata(cmd);  pause(SETTLE_S);
    I = grab(hSI, REF_NFRAMES, CHANNEL);
    if isempty(refVol), refVol = zeros([size(I) numel(zRef)],'single'); end
    refVol(:,:,k) = I;  got(k) = true; %#ok<AGROW>
end
hDM.senddata(baseline);  pause(SETTLE_S);
refVol = refVol(:,:,got);  zRef = zRef(got);

E = dm_axial_estimator_260913('prep', refVol, zRef, ...
        struct('nStrips',1,'pixFrac',0.04));

% How much does the image actually change per micron? If this is ~1.0 the
% field carries no axial information and nothing below can work.
rAdj = nan(1,numel(zRef)-1);
for k = 1:numel(zRef)-1
    rAdj(k) = corr_lin(refVol(:,:,k), refVol(:,:,k+1));
end
fprintf('  adjacent-plane correlation (1 um apart): %.4f median\n', median(rAdj));
fprintf('  end-to-end correlation (%.0f um apart) : %.4f\n', ...
        zRef(end)-zRef(1), corr_lin(refVol(:,:,1), refVol(:,:,end)));
if median(rAdj) > 0.999
    fprintf('  >> planes 1 um apart are essentially IDENTICAL. There is no\n');
    fprintf('     axial gradient to fit and no estimator can work here.\n');
end

%% ---------------- the test ----------------
fprintf('\n%8s %10s %12s %12s %10s\n', ...
        'cmd um','corr vs 0','est AVG um','est 1FRAME um','err avg');
n = numel(TEST_Z_UM);
estA = nan(1,n);  estS = nan(1,n);  rr = nan(1,n);
I0 = [];
for k = 1:n
    cmd = baseline + (TEST_Z_UM(k)/GAIN_UM_PER_COEFF)*zRow;
    if max(abs(cmd)) > CMD_LIMIT
        fprintf('%8.1f   SKIPPED (|cmd| %.3f)\n', TEST_Z_UM(k), max(abs(cmd)));
        continue
    end
    hDM.senddata(cmd);  pause(SETTLE_S);

    Iavg = grab(hSI, TEST_NAVG, CHANNEL);   % averaged: noise removed
    Isng = grab(hSI, 1, CHANNEL);           % single frame: what the loop sees

    if TEST_Z_UM(k) == 0, I0 = Iavg; end
    if ~isempty(I0), rr(k) = corr_lin(I0, Iavg); end

    % zInit = NaN so each estimate is independent -- no carried-over state to
    % mask a stuck fine tier.
    oA = dm_axial_estimator_260913('run', E, Iavg, NaN);
    oS = dm_axial_estimator_260913('run', E, Isng, NaN);
    estA(k) = oA.dz(1);  estS(k) = oS.dz(1);

    fprintf('%8.1f %10.3f %12.2f %12.2f %10.2f\n', ...
            TEST_Z_UM(k), rr(k), estA(k), estS(k), estA(k)-TEST_Z_UM(k));
end
hDM.senddata(baseline);  pause(SETTLE_S);

%% ---------------- verdict ----------------
ok = isfinite(estA);
fprintf('\n');
if nnz(ok) < 3
    fprintf('too few points to judge.\n');
else
    pA = polyfit(TEST_Z_UM(ok), estA(ok), 1);
    pS = polyfit(TEST_Z_UM(ok), estS(ok), 1);
    rA = estA(ok) - polyval(pA, TEST_Z_UM(ok));
    fprintf('AVERAGED  : est = %.3f * cmd %+.2f, residual RMS %.2f um\n', ...
            pA(1), pA(2), sqrt(mean(rA.^2)));
    fprintf('SINGLE    : est = %.3f * cmd %+.2f\n', pS(1), pS(2));

    rMin = min(rr(isfinite(rr) & TEST_Z_UM ~= 0));
    fprintf('\nimage correlation at the extremes: %.3f\n', rMin);
    if rMin > 0.98
        fprintf('=> THE MIRROR IS NOT MOVING. The frames barely change over %.1f um.\n', ...
                max(TEST_Z_UM)-min(TEST_Z_UM));
        fprintf('   Fix senddata / the DM before looking at the estimator at all.\n');
    elseif abs(pA(1)) < 0.2
        fprintf('=> THE MIRROR MOVES BUT THE ESTIMATOR IS BLIND (slope %.3f).\n', pA(1));
        fprintf('   The images change, the estimate does not. The fault is in the\n');
        fprintf('   fit, not the optics and not the controller.\n');
    elseif pA(1) < 0
        fprintf('=> SIGN INVERTED (slope %.3f). A loop would drive the motion.\n', pA(1));
    elseif abs(pA(1)-1) > 0.25
        fprintf('=> SCALE IS OFF by %.0f%% (slope %.3f). Calibration, not control:\n', ...
                100*(pA(1)-1), pA(1));
        fprintf('   GAIN_UM_PER_COEFF or the reference z axis is wrong.\n');
    elseif abs(pS(1)-1) > 0.3
        fprintf('=> WORKS ON AVERAGES, FAILS ON SINGLE FRAMES (%.3f vs %.3f).\n', ...
                pA(1), pS(1));
        fprintf('   That is photon noise, not the algorithm. The loop needs a\n');
        fprintf('   brighter plane, more pixels, or frame averaging.\n');
    else
        fprintf('=> THE ESTIMATOR IS FINE on real data, single frames included.\n');
        fprintf('   Any loop failure is in the controller or the timing.\n');
    end
end

R = struct('testZ',TEST_Z_UM,'estAvg',estA,'estSingle',estS,'corr',rr, ...
           'zRef',zRef,'adjCorr',rAdj);
fprintf('============================================================\n');
end

%% ---------------- helpers ----------------
function I = grab(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_est:grabDidNotStart','startGrab left acqState idle');
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
I = [];
try
    b = hSI.hDisplay.rollingStripeDataBuffer;
    if ~isempty(b) && ~isempty(b{1}) && ~isempty(b{1}{1}.roiData)
        I = pick_channel(b{1}{1}.roiData{1}, ch);
    end
catch
end
if isempty(I)
    try, I = pick_channel(hSI.hDisplay.lastStripeData.roiData{1}, ch); catch, end
end
assert(~isempty(I), 'could not read channel %d from the display buffer', ch);
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

function r = corr_lin(a,b)
a = double(a(:)) - mean(double(a(:)));
b = double(b(:)) - mean(double(b(:)));
r = (a.'*b)/(norm(a)*norm(b)+eps);
end

function wait_s(s)
t0 = tic; while toc(t0) < s, pause(0.05); drawnow limitrate; end
end

function cleanup(hDM, baseline, hSI, st)
try, hDM.senddata(baseline); catch ME, warning('dm_est:dm','%s',ME.message); end
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
catch ME, warning('dm_est:si','%s',ME.message); end
end
