function dm_defocus_check_260913()
%DM_DEFOCUS_CHECK_260913  30-second check: does the DM actually move the focus?
%
% Grabs one averaged frame at each of a few defocus commands and shows them
% side by side, with the correlation of each against the c=0 frame. No stack,
% no loop, no logging, nothing saved.
%
% Run this whenever the answer to "is the focus moving?" is not obvious by eye.
% A full reference stack takes 2.4 min to tell you the same thing, and if the
% images are bad for an unrelated reason (see the line-phase check below) it
% will tell you nothing at all.
%
% WHAT GOOD LOOKS LIKE
%   correlation falls off smoothly and symmetrically either side of c = 0,
%   reaching well below 1 by +-10 um. Structure visibly enters and leaves focus.
% WHAT BAD LOOKS LIKE
%   correlation stays near 1.0        -> the focus is NOT moving
%   correlation is near 0 everywhere  -> the images are noise; fix the imaging
%                                        before blaming the mirror
%
% Runqi Zhang / 2026-09-13.

%% ---------------- settings ----------------
SI_DIR       = 'C:\Users\dklab\Desktop\SI2018bR1_2018-12-19_4a9264c4fc';
Z2C_FILE     = fullfile(SI_DIR,'BAX331-Z2C.mat');
BASELINE_VAR = 'DMcommand_sys_aber_cor_930';
ROW_DEFOC    = 3;
CHANNEL      = 3;
Z_LIST_UM    = [-10 -5 0 5 10];
NFRAMES      = 30;
GAIN_UM_PER_COEFF = 6.656;
CMD_LIMIT    = 0.95;
SETTLE_S     = 0.3;

%% ---------------- connect ----------------
hSI = evalin('base','hSI');
hAO = evalin('base','hAOROIctrl');
hDM = hAO.hDMctrl;

S = load(Z2C_FILE);  fn = fieldnames(S);  Z2C = S.(fn{1});
zRow = Z2C(ROW_DEFOC,:);
assert(numel(zRow) == 97, 'Z2C row is %d long, expected 97', numel(zRow));

if ~evalin('base', sprintf('exist(''%s'',''var'')', BASELINE_VAR))
    error('dm_check:noBaseline','base workspace has no %s', BASELINE_VAR);
end
baseline = double(evalin('base', BASELINE_VAR));  baseline = baseline(:).';

st = struct('rollingAvg', hSI.hDisplay.displayRollingAverageFactor, ...
            'framesPerSlice', hSI.hStackManager.framesPerSlice, ...
            'numSlices', hSI.hStackManager.numSlices, ...
            'logging', hSI.hChannels.loggingEnable, ...
            'keepResOn', hSI.hScan2D.keepResonantScannerOn);
restore = onCleanup(@() cleanup(hDM, baseline, hSI, st)); %#ok<NASGU>
hSI.hChannels.loggingEnable = false;

%% ---------------- sanity checks first ----------------
fprintf('\n=============== DM DEFOCUS CHECK ===============\n');

% Line phase, because a bad one makes every frame garbage and every
% conclusion below meaningless. A resonant line is ~63 us.
lp = hSI.hScan2D.linePhase;
fprintf('line phase       : %.3f us\n', lp*1e6);
if abs(lp) > 20e-6
    error('dm_check:badLinePhase', ...
        ['linePhase is %.1f us. That is not a phase, it is a units error ' ...
         '(linePhase is in SECONDS, the GUI is in MICROSECONDS).\n' ...
         'Every frame will be scrambled. Fix it first:  ' ...
         'hSI.hScan2D.linePhase = -0.7e-6;'], lp*1e6);
end

chDisp = hSI.hChannels.channelDisplay;
assert(ismember(CHANNEL, chDisp), ...
    'channel %d is not in channelDisplay (%s) -- frames are read from the display', ...
    CHANNEL, mat2str(chDisp(:).'));

fprintf('baseline         : %s, max|cmd| %.3f\n', BASELINE_VAR, max(abs(baseline)));
fprintf('frame rate       : %.2f Hz, %d x %d\n', hSI.hRoiManager.scanFrameRate, ...
        hSI.hRoiManager.pixelsPerLine, hSI.hRoiManager.linesPerFrame);

hSI.hScan2D.keepResonantScannerOn = true;
fprintf('arming resonant scanner 1.0 s...\n');
wait_s(1.0);

%% ---------------- sweep ----------------
n = numel(Z_LIST_UM);
ims = cell(1,n);  ok = false(1,n);
for k = 1:n
    c   = Z_LIST_UM(k) / GAIN_UM_PER_COEFF;
    cmd = baseline + c*zRow;
    if max(abs(cmd)) > CMD_LIMIT
        fprintf('  z %+5.1f um : SKIPPED, |cmd| %.3f\n', Z_LIST_UM(k), max(abs(cmd)));
        continue
    end
    hDM.senddata(cmd);  pause(SETTLE_S);
    ims{k} = grab(hSI, NFRAMES, CHANNEL);
    ok(k)  = true;
    fprintf('  z %+5.1f um (c %+5.2f, max|cmd| %.3f) grabbed\n', ...
            Z_LIST_UM(k), c, max(abs(cmd)));
end
hDM.senddata(baseline);  pause(SETTLE_S);

i0 = find(Z_LIST_UM == 0 & ok, 1);
assert(~isempty(i0), 'the c = 0 frame was not acquired');

%% ---------------- report ----------------
fprintf('\n%10s %12s %12s\n','z (um)','corr vs 0','mean counts');
r = nan(1,n);
for k = 1:n
    if ~ok(k), continue; end
    r(k) = corr_lin(ims{i0}, ims{k});
    fprintf('%10.1f %12.3f %12.0f\n', Z_LIST_UM(k), r(k), mean(ims{k}(:)));
end

rr = r(ok & Z_LIST_UM ~= 0);
fprintf('\n');
if all(rr > 0.95)
    fprintf('*** THE FOCUS IS NOT MOVING. ***\n');
    fprintf('Frames %g um apart are still %.3f correlated. Check, in order:\n', ...
            max(abs(Z_LIST_UM)), min(rr));
    fprintf('  1. is the DM actually in the beam path / powered?\n');
    fprintf('  2. is Z2C row %d really defocus for this mirror?\n', ROW_DEFOC);
    fprintf('  3. does senddata reach the mirror (any error from the driver)?\n');
elseif median(abs(diff(r(ok)))) < 0.01
    fprintf('*** SUSPICIOUS: correlation barely changes between planes. ***\n');
elseif mean(rr) < 0.1
    fprintf('*** THE IMAGES MAY BE NOISE, NOT DEFOCUS. ***\n');
    fprintf('Everything decorrelates to ~0, which is what two noise frames do.\n');
    fprintf('Check the displayed image looks like tissue before trusting this.\n');
else
    fprintf('*** LOOKS RIGHT. *** Correlation falls off with defocus as it should.\n');
    fprintf('Structure should visibly enter and leave focus in the figure.\n');
end

%% ---------------- figure ----------------
f = figure('Color','w','Name','DM defocus check','NumberTitle','off');
for k = 1:n
    if ~ok(k), continue; end
    subplot(1,n,k);
    im = ims{k};
    lo = qtile(im(:),0.02);  hi = qtile(im(:),0.998);
    imagesc(im, [lo max(hi,lo+1)]);  axis image off;  colormap(f, gray);
    title(sprintf('%+g um\nr=%.2f', Z_LIST_UM(k), r(k)));
end
fprintf('================================================\n');
end

%% ---------------- helpers ----------------
function I = grab(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.hDisplay.displayRollingAverageFactor = n;
hSI.startGrab();
if strcmp(hSI.acqState,'idle')
    error('dm_check:grabDidNotStart','startGrab left acqState idle -- channel selected for Display?');
end
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
% lastFrame is indexed by channelDisplay, NOT by lastFrameChannels -- see the
% getter in Display.m. Using the wrong list silently returns the wrong channel
% or nothing at all, and the two lists diverge as soon as any channel is
% enabled for logging only.
dispCh = hSI.hChannels.channelDisplay;
i = find(dispCh(:).' == ch, 1);
assert(~isempty(i), 'channel %d is not in channelDisplay (%s)', ch, mat2str(dispCh(:).'));
fr = hSI.hDisplay.lastFrame;
assert(i <= numel(fr) && ~isempty(fr{i}), ...
    'lastFrame has %d entries, wanted index %d for channel %d', numel(fr), i, ch);
I = single(fr{i});
end

function r = corr_lin(a, b)
a = double(a(:)) - mean(double(a(:)));
b = double(b(:)) - mean(double(b(:)));
r = (a.'*b)/(norm(a)*norm(b)+eps);
end

function q = qtile(x, p)
x = sort(double(x(isfinite(x))));
if isempty(x), q = 0; return; end
q = x(max(1,min(numel(x), round(p*(numel(x)-1))+1)));
end

function wait_s(s)
t0 = tic; while toc(t0) < s, pause(0.05); drawnow limitrate; end
end

function cleanup(hDM, baseline, hSI, st)
try, hDM.senddata(baseline); catch ME, warning('dm_check:dm','%s',ME.message); end
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
catch ME, warning('dm_check:si','%s',ME.message); end
end
