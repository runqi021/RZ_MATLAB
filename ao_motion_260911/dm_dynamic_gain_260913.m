function R = dm_dynamic_gain_260913()
%DM_DYNAMIC_GAIN_260913  Why does a continuously-driven DM deliver only 65%?
%
% WHAT IS ALREADY KNOWN (rig, 2026-09-13), so this does not re-test it:
%   estimator, static, averaged frames ............ gain 0.990
%   estimator, static, single frames ............... gain 0.966
%   read path, 2 s plateaus, loop's own frame source gain 1.000, 100% new frames
%   loop, 5 um sinusoid, senddata every frame ...... gain 0.653
%
% So the estimator is right, the frames are right, and a SETTLED mirror is
% right. Only the continuously-driven mirror is short. A first-order settling
% lag does not explain it: 0.65 at 0.5 Hz needs tau ~ 370 ms, which would have
% left the 2 s plateaus at ~97% instead of 1.000 -- and a low-pass falls
% monotonically with frequency, while the measured gain RISES from 0.656 at
% 0.5 Hz to 0.710 at 3 Hz before falling.
%
% THE REAL DIFFERENCE IS THE COMMAND RATE: 8 senddata calls in the staircase
% versus 60 per second in the sweep.
%
% THIS SEPARATES THE THREE CANDIDATES BY THEIR SIGNATURE:
%   slew / rate limit  -> gain falls with AMPLITUDE x FREQUENCY (peak velocity)
%   linear filter      -> gain falls with FREQUENCY alone, amplitude-independent
%   static nonlinearity-> gain falls with AMPLITUDE alone, frequency-independent
%
% It sweeps amplitude at fixed frequency and frequency at fixed amplitude, and
% reports gain against peak velocity (2*pi*f*A) so a rate limit shows up as a
% single curve collapsing both sweeps onto one line.
%
% Open loop throughout -- no correction, nothing logged, DM returned to
% baseline. ~3 min.
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

% (amplitude um, frequency Hz) pairs. Chosen so peak velocity 2*pi*f*A spans
% two decades and the same velocity is reached by different (A,f) combinations
% -- that overlap is what tests the slew hypothesis.
COND = [ 5.0 0.05;     % 1.6 um/s   very slow, large
         5.0 0.20;     % 6.3
         5.0 0.50;     % 15.7       the condition that measured 0.65
         5.0 2.00;     % 62.8
         1.0 0.50;     % 3.1        same f as above, 5x smaller
         2.0 0.50;     % 6.3        same velocity as (5.0, 0.20)
        10.0 0.25;     % 15.7       same velocity as (5.0, 0.50), 2x amplitude
         1.0 2.50];    % 15.7       same velocity again, 5x smaller amplitude
SEG_SECONDS  = 12;

% Command the mirror only every Nth frame, to test the rate directly. 1 = every
% frame (what the loop does).
CMD_EVERY    = 1;

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
    error('dm_dyn:noBaseline','base workspace has no %s', BASELINE_VAR);
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
fprintf('\n============ DM DYNAMIC GAIN ============\n');
fprintf('%.2f Hz, commanding every %d frame(s) = %.1f Hz command rate\n', ...
        fps, CMD_EVERY, fps/CMD_EVERY);

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
fprintf('reference: %d planes, %.1f to %.1f um\n\n', numel(zRef), zRef(1), zRef(end));

%% ---------------- sweep ----------------
n = size(COND,1);
amp = COND(:,1).';  frq = COND(:,2).';
vel = 2*pi*frq.*amp;
g   = nan(1,n);
fprintf('%8s %8s %12s %10s\n','amp um','f Hz','peak um/s','gain');
for i = 1:n
    L = seg(hSI, hDM, E, baseline, zRow, amp(i), frq(i), SEG_SECONDS, ...
            GAIN_UM_PER_COEFF, CMD_LIMIT, CHANNEL, CMD_EVERY);
    a = lockin(L.dz, L.t, frq(i));
    c = lockin(L.cmd, L.t, frq(i));
    g(i) = a / max(c, eps);
    fprintf('%8.1f %8.2f %12.1f %10.3f\n', amp(i), frq(i), vel(i), g(i));
end
hDM.senddata(baseline);  pause(SETTLE_S);

%% ---------------- diagnosis ----------------
fprintf('\n---------------- which signature? ----------------\n');

% Same frequency (0.5 Hz), different amplitudes -> amplitude dependence
sel = frq == 0.5;
if nnz(sel) >= 3
    fprintf('at 0.5 Hz, amplitude %s um -> gain %s\n', ...
            mat2str(amp(sel)), mat2str(round(g(sel),3)));
end
% Same amplitude (5 um), different frequencies -> frequency dependence
sel = amp == 5;
if nnz(sel) >= 3
    fprintf('at 5 um, frequency  %s Hz -> gain %s\n', ...
            mat2str(frq(sel)), mat2str(round(g(sel),3)));
end
% Equal peak velocity, different (A,f) -> the slew test
[~, ord] = sort(vel);
fprintf('\nsorted by peak velocity:\n%10s %8s %8s %8s\n','um/s','amp','f','gain');
for i = ord
    fprintf('%10.1f %8.1f %8.2f %8.3f\n', vel(i), amp(i), frq(i), g(i));
end

sameVel = abs(vel - 15.7) < 1;
if nnz(sameVel) >= 2
    sv = g(sameVel);
    fprintf('\nSAME peak velocity (~15.7 um/s) at amplitudes %s: gain %s\n', ...
            mat2str(amp(sameVel)), mat2str(round(sv,3)));
    if max(sv) - min(sv) < 0.1
        fprintf('  >> gain is set by PEAK VELOCITY, not by amplitude or frequency\n');
        fprintf('     separately. That is a SLEW / RATE LIMIT in the DM or its\n');
        fprintf('     driver. Fix: command smaller excursions, or accept a known\n');
        fprintf('     velocity ceiling and design the loop around it.\n');
    else
        fprintf('  >> gain differs at equal velocity, so it is NOT a slew limit.\n');
    end
end

slow = g(frq <= 0.05);
if ~isempty(slow) && mean(slow) > 0.9
    fprintf('\nAt 0.05 Hz the gain is %.3f -- essentially perfect. So nothing is\n', mean(slow));
    fprintf('wrong with commanding every frame per se; it is the RATE OF CHANGE\n');
    fprintf('that costs amplitude.\n');
elseif ~isempty(slow)
    fprintf('\nEven at 0.05 Hz (peak %.1f um/s) the gain is only %.3f. The loss is\n', ...
            2*pi*0.05*5, mean(slow));
    fprintf('NOT about speed -- calling senddata every frame is itself the problem.\n');
    fprintf('Try CMD_EVERY = 2 or 4 and re-run.\n');
end

R = struct('amp',amp,'freq',frq,'peakVel',vel,'gain',g,'cmdEvery',CMD_EVERY, ...
           'fps',fps,'zRef',zRef);
fprintf('=========================================\n');
end

%% ---------------- one segment ----------------
function L = seg(hSI, hDM, E, baseline, zRow, ampUm, fHz, secs, umPerCoeff, ...
                 cmdLimit, ch, cmdEvery)
n = round(secs * hSI.hRoiManager.scanFrameRate);
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
L = struct('t',nan(1,n),'dz',nan(1,n),'cmd',nan(1,n));
lastNo = -1; nSvc = 0; cur = 0;
hDM.senddata(baseline);  pause(0.2);
t0 = tic;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_dyn:grabDidNotStart','startGrab left acqState idle');
end
while ~strcmp(hSI.acqState,'idle')
    try, fno = hSI.hDisplay.lastFrameNumber; catch, pause(0.002); continue; end
    if fno == lastNo, pause(0.001); drawnow limitrate; continue; end
    I = read_loop_path(hSI, ch);
    if isempty(I), pause(0.0005); continue; end
    lastNo = fno; nSvc = nSvc + 1;
    if nSvc > n, break; end
    tNow = toc(t0);
    o = dm_axial_estimator_260913('run', E, single(I), NaN);
    L.t(nSvc) = tNow;  L.dz(nSvc) = o.dz(1);  L.cmd(nSvc) = cur;
    if mod(nSvc, cmdEvery) == 0
        cur = ampUm * sin(2*pi*fHz*tNow);
        c = baseline + (cur/umPerCoeff)*zRow;
        if max(abs(c)) <= cmdLimit, hDM.senddata(c); end
    end
    drawnow limitrate;
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
k = isfinite(L.t);
L.t=L.t(k); L.dz=L.dz(k); L.cmd=L.cmd(k);
end

%% ---------------- helpers ----------------
function a = lockin(x, t, f)
k = isfinite(x) & isfinite(t);
if nnz(k) < 20, a = NaN; return; end
x = x(k) - mean(x(k));  t = t(k);
a = hypot(2*mean(x.*cos(2*pi*f*t)), 2*mean(x.*sin(2*pi*f*t)));
end

function I = read_loop_path(hSI, ch)
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
    error('dm_dyn:grabDidNotStart','startGrab left acqState idle');
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
try, hDM.senddata(baseline); catch ME, warning('dm_dyn:dm','%s',ME.message); end
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
catch ME, warning('dm_dyn:si','%s',ME.message); end
end
