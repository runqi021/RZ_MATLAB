function stage_backlash_260911()
%STAGE_BACKLASH_260911  Direction-dependent step deviation of the Galil stage.
%
% Measures, per axis and per step size, how far the stage actually lands when a
% target is approached going UP versus going DOWN -- i.e. backlash / lost motion
% on reversal, plus per-step scale error.
%
% NO ACQUISITION. Moves the motor and reads position registers. The laser,
% shutter, scanner and DM are never touched.
%
% ----------------------------------------------------------------------------
% READ THIS FIRST -- WHY THERE IS A PROBE MODE
% ----------------------------------------------------------------------------
% `hSI.hMotors.motorPosition` CANNOT see backlash on this rig.  Its read path is
%     motorPosition -> positionAbsoluteRaw -> getStepperCount() -> Galil "TD"
% (+dabs\+galil\DMC4040.m:132-134,164-165).  For a stepper axis TD is the number
% of step pulses the controller GENERATED.  It is an open-loop echo of the
% command: it increments whether or not the stage moved, so it reports the
% target essentially exactly, in both directions, always.  Characterising
% backlash from TD alone would produce a confident table of zeros.
%
% The driver also exposes two registers ScanImage does NOT use for position:
%     TP  getEncoderPosition()    - encoder register.  REAL position, if wired.
%     RP  getReferencePosition()  - commanded reference.
% If this stage has encoders, TP is the ruler and (TD - TP) is precisely the
% deviation we are after.  If it does not, TP will sit at zero or drift with no
% relation to the move, and NO controller register on this machine can measure
% backlash -- it would need the image or an external dial indicator.
%
% MODE = 'probe' answers that question with one small move.  RUN IT FIRST.
%
% ----------------------------------------------------------------------------
% TRAVEL ENVELOPE (user-stated 2026-09-11, and enforced below)
%   start at x,y = 0,0 -> move to CENTER_XY = [500 500] and sweep about there
%   MAX DISTANCE FROM THAT CENTRE:  x +-3000   y +-3000   z +-300
% These are radii from the centre, NOT total spans.  Sweeps are symmetric about
% the centre, so the worst case is half the span plus one run-up step, and that
% is what gets checked.  Every target is validated BEFORE anything moves; a plan
% that does not fit is refused, never clipped -- clipping would quietly change
% the step pattern being characterised.
% ----------------------------------------------------------------------------
% Runqi Zhang / 2026-09-11.  Written against the SI2018b source on this rig.
% NOT YET RUN.  Read it before you trust it.

%% ========================= USER SETTINGS =========================
MODE       = 'full';      % 'probe' first. Then 'full'.

DRY_RUN    = false;         % 'full' only: print the move plan, move nothing.

CENTER_XY  = [500 500];    % move here before sweeping x/y. [] = stay put.

AXES       = [1 2 3];      % 1 = x, 2 = y, 3 = z
% Step sizes to test, per axis. NOTE the device quantum: x/y move in 0.7815 um
% counts, z in 0.3125 um counts, so not every value here is exactly commandable
% -- 0.5 um is BELOW one count on x/y and is 1.6 counts on z. They are left in
% deliberately; the plan printout reports the quantisation for each so the
% result is interpreted with it in view rather than blamed on backlash.
STEPS_UM   = { [0.1 0.5 1 2 5 10 50 100 200 500], ...        % x
               [0.1 0.5 1 2 5 10 50 100 200 500], ...        % y
               [0.1 0.5 1 2 5 10 20 50 100] };      % z  (user: 5 um expected best)

NSTEPS_MAX = 5;            % steps per direction, reduced per step size to fit
NREPS      = 3;            % sweeps per (axis, step size)
SETTLE_S   = 0.5;          % after moveWaitForFinish, before reading registers
REREAD_S   = 0;            % >0: read again after this delay, to catch creep

PROBE_UM   = [20 20 10];   % probe move per axis
PROBE_AXIS = 1;

% HARD limit on DISTANCE FROM THE SWEEP CENTRE, per axis (user-stated).
% Enforced as: |offset from centre| + one run-up step <= SPAN_FRAC*MAX_EXC_UM.
% Checked for every target before the first move; a plan that does not fit is
% refused, never clipped.
% x/y widened to +-3000 on 2026-09-11 (user). NOTE this does NOT make x/y
% measurable: TP returns a constant on those axes, so the result is zeros at any
% range. The wider window only improves statistics on large steps -- at +-1500,
% a 500 um step only fit 2 steps per direction (n=18 instead of 30).
MAX_EXC_UM = [3000 3000 300];
SPAN_FRAC  = 0.9;            % use at most this fraction of it

OUT_DIR    = fileparts(mfilename('fullpath'));

%% ========================= CONNECT =========================
hSI = evalin('base','hSI');
hMotors = hSI.hMotors;

% hMotors.motorPosition is the SI-frame position (relative to the soft origin).
% hLSC is the raw dabs.galil.DMC4040, the only route to TP and RP.
try
    hLSC = hMotors.hMotor(1).hLSC;
catch
    error('stage_backlash:noLSC', ...
        ['Could not reach hSI.hMotors.hMotor(1).hLSC -- without the raw Galil ' ...
         'object there is no way to read TP, and motorPosition alone cannot ' ...
         'measure backlash. See the header.']);
end

% um per device unit. MDF has [0.7815 -0.7815 0.3125] -- note y is NEGATIVE and
% z is finer than x/y. Read it from the object so this never silently drifts
% from the MDF.
% *** positionDeviceUnits IS IN METRES. *** The MDF literal reads
% [0.7815 -0.7815 0.3125] but SI stores it as [7.815e-07 ...] m. Taking it as
% um/count scales every register reading by 1e-6, which makes a real 25-count
% move print as +0.000 and fakes a "the encoder is dead" verdict. Cost one
% bogus probe on 2026-09-11. Convert, then sanity-check the magnitude.
devUm = [0.7815 -0.7815 0.3125];
try
    v = hLSC.positionDeviceUnits;
    if numel(v) >= 3 && all(isfinite(v(1:3))) && all(v(1:3) ~= 0)
        cand = v(1:3) * 1e6;                       % m -> um
        if all(abs(cand) > 0.01) && all(abs(cand) < 100)
            devUm = cand;
        else
            warning('stage_backlash:devUnitsOdd', ...
                ['positionDeviceUnits*1e6 = [%g %g %g] um/count is outside the ' ...
                 'plausible 0.01-100 range; keeping the MDF defaults instead.'], cand);
        end
    end
catch
    warning('stage_backlash:devUnits', ...
        'could not read positionDeviceUnits; assuming [%g %g %g]', devUm);
end
fprintf('device units (um/count): [%.4f %.4f %.4f]\n', devUm);

pStart = hMotors.motorPosition;
fprintf('start position: x %.3f  y %.3f  z %.3f\n', pStart(1), pStart(2), pStart(3));

% Whatever happens -- error, Ctrl+C -- go back to where the user left it. This
% is the lesson from auto_acq, where a Ctrl+C left the rig in a changed state.
restore = onCleanup(@() go_home(hMotors, pStart)); %#ok<NASGU>

% Sweep centre: x/y offset to CENTER_XY, z left where it is.
p0 = pStart;
if ~isempty(CENTER_XY)
    p0(1) = pStart(1) + CENTER_XY(1);
    p0(2) = pStart(2) + CENTER_XY(2);
end

%% ========================= PROBE =========================
if strcmpi(MODE,'probe')
    probe(hMotors, hLSC, pStart, PROBE_AXIS, PROBE_UM(PROBE_AXIS), SETTLE_S, devUm);
    return
end
if ~strcmpi(MODE,'full')
    error('stage_backlash:mode','MODE must be ''probe'' or ''full''');
end

%% ========================= PLAN + BUDGET CHECK =========================
% Sweeps are CENTRED on p0 and symmetric, so the worst-case excursion is half
% the span plus one run-up step. Build every target first and check the budget
% before anything moves, rather than discovering it halfway through a sweep.
plan = {};
for ai = 1:numel(AXES)
    ax = AXES(ai);
    budget = SPAN_FRAC * MAX_EXC_UM(ax);
    for s = STEPS_UM{ax}(:).'
        % largest even nst whose half-span plus a run-up step fits the budget
        nst = min(NSTEPS_MAX, 2*floor((budget/s - 1)));
        nst = nst - mod(nst,2);
        if nst < 2
            warning('stage_backlash:stepTooBig', ...
                ['axis %d step %g um: even 2 steps would need %g um, past the ' ...
                 '%g um budget. SKIPPED.'], ax, s, 2*s, budget);
            continue
        end
        offs   = ((0:nst) - nst/2) * s;
        maxOff = max(abs(offs)) + s;           % + run-up
        if maxOff > budget
            error('stage_backlash:budget', ...
                'axis %d step %g: excursion %g um exceeds budget %g um', ...
                ax, s, maxOff, budget);
        end

        % Device quantisation. The controller can only command whole counts, so
        % a step that is not an integer number of them is rounded -- and that
        % rounding is a fixed, direction-INDEPENDENT error. Reporting it here
        % keeps it from being read as backlash later.
        q      = abs(devUm(ax));
        cnt    = s / q;
        cntR   = round(cnt);
        actual = cntR * q;
        if cntR < 1
            warning('stage_backlash:subCount', ...
                ['axis %d step %g um is BELOW one device count (%g um). The ' ...
                 'stage cannot resolve it -- expect either no motion or an ' ...
                 'erratic 0/1-count stutter, NOT a backlash measurement.'], ...
                ax, s, q);
        elseif abs(cnt-cntR)/cnt > 0.02
            warning('stage_backlash:quantised', ...
                ['axis %d step %g um = %.2f counts of %g um; the controller ' ...
                 'will actually step %g um (%+.1f%%).'], ...
                ax, s, cnt, q, actual, 100*(actual-s)/s);
        end

        plan{end+1} = struct('ax',ax,'step',s,'nst',nst,'offs',offs, ...
                             'maxOff',maxOff,'cnt',cnt,'actual',actual); %#ok<AGROW>
    end
end
if isempty(plan), error('stage_backlash:emptyPlan','nothing fits the budget'); end

nMoves = 0;
for k = 1:numel(plan), nMoves = nMoves + NREPS*(2*(plan{k}.nst+1) + 1); end
fprintf('\nsweep centre: x %.1f  y %.1f  z %.1f\n', p0(1), p0(2), p0(3));
fprintf('plan: %d conditions, %d reps, %d moves\n', numel(plan), NREPS, nMoves);
fprintf('\n%-5s %8s %6s %9s %9s %11s %9s\n', ...
        'axis','step','steps','counts','actual','excursion','budget');
for k = 1:numel(plan)
    p = plan{k};
    flag = '';
    if p.cnt < 1, flag = '  <1 COUNT'; elseif abs(p.actual-p.step) > 0.02*p.step, flag = '  quantised'; end
    fprintf('%-5d %8.2f %6d %9.2f %9.4f %11.1f %9.0f%s\n', ...
            p.ax, p.step, p.nst, p.cnt, p.actual, p.maxOff, ...
            SPAN_FRAC*MAX_EXC_UM(p.ax), flag);
end
fprintf(['counts = step / device quantum (x,y %.4f um; z %.4f um). "actual" is ' ...
         'what\nthe controller can really step. Any gap there is quantisation, ' ...
         'not backlash.\n'], abs(devUm(1)), abs(devUm(3)));

if DRY_RUN
    fprintf(['\nDRY_RUN is true -- nothing moved. Set DRY_RUN = false to run ' ...
             'it for real.\n']);
    return
end

if ~isempty(CENTER_XY)
    fprintf('\nmoving to sweep centre...\n');
    hMotors.moveStartRelative(p0);
    hMotors.moveWaitForFinish([1 2 3]);
    gui_refresh();
end

%% ========================= SWEEP =========================
% Field ORDER must match the struct() built in the loop, or the assignment errors.
R = struct('ax',{},'step',{},'rep',{},'dir',{},'k',{},'target',{}, ...
           'tdCnt',{},'tpCnt',{},'rpCnt',{},'td',{},'tp',{},'rp',{},'si',{});

for k = 1:numel(plan)
    ax = plan{k}.ax;  s = plan{k}.step;  offs = plan{k}.offs;
    fprintf('\n[axis %d, step %g um, %d steps]\n', ax, s, plan{k}.nst);
    for rep = 1:NREPS
        % Run-up below the first target, so the first recorded point is already
        % direction-defined instead of carrying the previous sweep's reversal
        % state.
        goto(hMotors, p0, ax, offs(1)-s, SETTLE_S);

        for dirSign = [+1 -1]
            if dirSign > 0, seq = 1:numel(offs); else, seq = numel(offs):-1:1; end
            for j = seq
                goto(hMotors, p0, ax, offs(j), SETTLE_S);
                [td, tp, rp] = read_raw(hLSC);     % counts
                si = hMotors.motorPosition;
                % Store counts AND um. Counts are the primary record; if the
                % conversion ever needs revisiting the raw data survives.
                R(end+1) = struct('ax',ax,'step',s,'rep',rep, ...
                    'dir',dirSign,'k',j,'target',p0(ax)+offs(j), ...
                    'tdCnt',td(ax),'tpCnt',tp(ax),'rpCnt',rp(ax), ...
                    'td',td(ax)*devUm(ax),'tp',tp(ax)*devUm(ax), ...
                    'rp',rp(ax)*devUm(ax),'si',si(ax)); %#ok<AGROW>
                if REREAD_S > 0
                    pause(REREAD_S);
                    [td2, tp2] = read_raw(hLSC);
                    fprintf('     creep after %.1fs: TD %+0.2f  TP %+0.2f counts\n', ...
                            REREAD_S, td2(ax)-td(ax), tp2(ax)-tp(ax));
                end
            end
        end
        fprintf('  rep %d done\n', rep);
    end
end

%% ========================= REPORT =========================
fprintf('\n================ DIRECTION-DEPENDENT DEVIATION ================\n');
fprintf('error = measured (TP) - commanded target, in um\n\n');
fprintf('%-4s %8s %10s %10s %12s %8s\n','ax','step','err UP','err DOWN','backlash','sd');
for k = 1:numel(plan)
    ax = plan{k}.ax;  s = plan{k}.step;
    m  = [R.ax] == ax & [R.step] == s;
    eu = [R(m & [R.dir] > 0).tp] - [R(m & [R.dir] > 0).target];
    ed = [R(m & [R.dir] < 0).tp] - [R(m & [R.dir] < 0).target];
    fprintf('%-4d %8.2f %10.3f %10.3f %12.3f %8.3f\n', ax, s, ...
            mean(eu), mean(ed), mean(ed)-mean(eu), ...
            sqrt((std(eu)^2 + std(ed)^2)/2));
end
fprintf(['\nbacklash = mean(err DOWN) - mean(err UP). A systematic non-zero ' ...
         'value\nis lost motion on reversal; sd is repeatability.\n']);

out = fullfile(OUT_DIR, sprintf('stage_backlash_%s.mat', datestr(now,'yymmdd_HHMM')));
save(out, 'R','plan','p0','pStart','devUm','STEPS_UM','NREPS','SETTLE_S','MAX_EXC_UM');
fprintf('\nsaved %s\n', out);
plot_result(R, plan, OUT_DIR);
end

%% ========================== HELPERS ==========================
function probe(hMotors, hLSC, p0, ax, dUm, settle, devUm)
% One small move. Decides whether TP is a live encoder or dead weight.
names = {'x','y','z'};
fprintf('\n=== PROBE: axis %s, %+g um ===\n', names{ax}, dUm);

q       = abs(devUm(ax));
expCnt  = dUm / q;                       % counts the command should produce
fprintf('quantum %.4f um/count -> a %+g um command should move %+.2f counts\n', ...
        q, dUm, expCnt);

% Two reads before moving, to show the registers are stable and that a zero
% delta later means "did not move" rather than "read noise swamped it".
[a0,b0,c0] = read_raw(hLSC);
pause(0.2);
[td0,tp0,rp0] = read_raw(hLSC);
fprintf('\n%-8s %14s %14s %14s\n','', 'TD (counts)','TP (counts)','RP (counts)');
fprintf('%-8s %14.1f %14.1f %14.1f\n','rest 1', a0(ax), b0(ax), c0(ax));
fprintf('%-8s %14.1f %14.1f %14.1f\n','rest 2', td0(ax), tp0(ax), rp0(ax));

siBefore = hMotors.motorPosition;
goto(hMotors, p0, ax, dUm, settle);
siAfter  = hMotors.motorPosition;
[td1,tp1,rp1] = read_raw(hLSC);
fprintf('%-8s %14.1f %14.1f %14.1f\n','moved',  td1(ax), tp1(ax), rp1(ax));

dTD = td1(ax)-td0(ax);  dTP = tp1(ax)-tp0(ax);  dRP = rp1(ax)-rp0(ax);
fprintf('%-8s %14.1f %14.1f %14.1f   <- expected %+.2f\n\n', ...
        'DELTA', dTD, dTP, dRP, expCnt);
fprintf('in um:   TD %+8.3f   TP %+8.3f   RP %+8.3f   (commanded %+g)\n', ...
        dTD*q, dTP*q, dRP*q, dUm);
fprintf('SI motorPosition: %.3f -> %.3f  (delta %+.3f um)\n\n', ...
        siBefore(ax), siAfter(ax), siAfter(ax)-siBefore(ax));

% Did the stage move at all? TD is the command echo, so it is the witness for
% "a move was issued", NOT for "the stage arrived".
if abs(dTD) < 0.5*abs(expCnt)
    fprintf(['VERDICT: NO MOVE WAS ISSUED. TD changed by %.1f counts, expected ' ...
        '%.1f.\n  This is a command/plumbing problem, not a stage problem -- ' ...
        'stop and\n  fix it before reading anything into TP.\n'], dTD, expCnt);
    return
end
fprintf('TD moved %.1f of the expected %.1f counts: the command WAS issued.\n', ...
        dTD, expCnt);

if ~isfinite(dTP)
    fprintf(['VERDICT: TP is unreadable (getEncoderPosition threw). No encoder\n' ...
        '  register is reachable -- treat as the no-encoder case below.\n']);
elseif abs(dTP) < 0.1*abs(expCnt)
    fprintf(['VERDICT: TP did NOT follow (%.1f counts vs %.1f expected).\n' ...
        '  No usable encoder on this axis, so NO register on this controller\n' ...
        '  can measure backlash -- TD only echoes the command. Characterising\n' ...
        '  it needs an external ruler: the 2P image itself (phase-correlate a\n' ...
        '  frame before and after a move) or a dial indicator.\n' ...
        '  Do not run MODE = ''full''; it would report zeros.\n'], dTP, expCnt);
else
    fprintf(['VERDICT: TP tracked the move (%.1f counts = %.3f um for a %g um\n' ...
        '  command, error %+.3f um). The encoder is live -- MODE = ''full'' will\n' ...
        '  measure real deviation.\n'], dTP, dTP*q, dUm, dTP*q-dUm);
end
end

function [td, tp, rp] = read_raw(hLSC)
% All three registers as RAW DEVICE COUNTS -- deliberately NOT converted.
% The unit conversion is the one thing that already went wrong here (metres vs
% microns silently scaled a real move to +0.000), so counts are what gets
% reported and compared; um is derived at the point of display only.
% A register that cannot be read comes back NaN rather than silently falling
% back to another register.
td = nan(1,3); tp = nan(1,3); rp = nan(1,3);
try, v = hLSC.getStepperCount();      n = min(3,numel(v)); td(1:n) = v(1:n); catch, end
try, v = hLSC.getEncoderPosition();   n = min(3,numel(v)); tp(1:n) = v(1:n); catch, end
try, v = hLSC.getReferencePosition(); n = min(3,numel(v)); rp(1:n) = v(1:n); catch, end
end

function goto(hMotors, p0, ax, offUm, settle)
% Absolute move in the SI frame. NOTE: despite its name, moveStartRelative on
% the Motors component takes coordinates in the same frame as motorPosition --
% this is the call auto_acq_260910 uses and it is proven on this rig.
t = p0;  t(ax) = p0(ax) + offUm;
hMotors.moveStartRelative(t);
hMotors.moveWaitForFinish([1 2 3]);
if settle > 0
    t0 = tic;
    while toc(t0) < settle, pause(0.02); drawnow limitrate; end
end
gui_refresh();
end

function go_home(hMotors, p0)
try
    hMotors.moveStartRelative(p0);
    hMotors.moveWaitForFinish([1 2 3]);
    gui_refresh();
    fprintf('[cleanup] stage returned to x %.3f y %.3f z %.3f\n', p0(1), p0(2), p0(3));
catch ME
    warning('stage_backlash:noHome','COULD NOT RETURN THE STAGE: %s', ME.message);
end
end

function gui_refresh()
% Push the stage position into the ScanImage Motor Controls boxes.
%
% NOTHING in ScanImage refreshes them in response to a motor move. The display
% is driven by SIController.changedMotorPosition, and that is reached only from
%   SI.m:1034            zzzEndOfAcquisitionMode  - END of every FOCUS/GRAB/LOOP
%   SIController.m:901   a listener on hUserFunctions 'motorPositionUpdate',
%                        which is notified ONLY from StackManager (414/432/526/569)
% So auto_acq appears to track the stage purely because it grabs after every
% tile. A no-acquisition script like this one leaves the display frozen even
% while the stage is moving -- which is what happened on the 921-move run.
%
% We cannot fire the event ourselves. UserFunctions.m:280 declares it
%     events (NotifyAccess = {?scanimage.interfaces.Class})
% so only SI's own classes may notify it; an external script would error.
% Calling the controller method directly is the route the GUI itself uses
% (motorControlsV5.m:36), and it is reachable -- Hidden, but not private.
%
% Cached, revalidated if SI is restarted, and never allowed to break a sweep.
persistent hCtrl tried
if isempty(tried) || (~isempty(hCtrl) && isobject(hCtrl) && ~isvalid(hCtrl))
    tried = true;  hCtrl = [];
    try, hCtrl = evalin('base','hSICtl'); catch, end
end
if isempty(hCtrl), return; end
try, hCtrl.changedMotorPosition(); catch, end
end

function plot_result(R, plan, outDir)
n  = numel(plan);
hf = figure('Color','w','Position',[100 100 min(1500,380*n) 340]);
for k = 1:n
    ax = plan{k}.ax;  s = plan{k}.step;
    a  = subplot(1,n,k); hold(a,'on');
    m  = [R.ax] == ax & [R.step] == s;
    up = m & [R.dir] > 0;  dn = m & [R.dir] < 0;
    plot(a, [R(up).target], [R(up).tp]-[R(up).target], 'o-', 'DisplayName','up');
    plot(a, [R(dn).target], [R(dn).tp]-[R(dn).target], 's-', 'DisplayName','down');
    xlabel(a,'target (um)'); ylabel(a,'TP - target (um)');
    title(a, sprintf('axis %d, step %g um', ax, s));
    legend(a,'Location','best'); box(a,'off'); grid(a,'on');
end
% opengl software: the NVIDIA path on these machines drops ~1% of exportgraphics
% PNGs as blank stubs. See reference_matlab_opengl_export_stubs.
try, opengl('software'); catch, end

% exportgraphics only exists from R2020a. This script runs inside the ScanImage
% MATLAB, which is older than the analysis one -- fall back to print() there
% rather than losing the figure to an undefined-function error at the very end
% of a long sweep.
png = fullfile(outDir,'stage_backlash.png');
if exist('exportgraphics','file') == 2 || exist('exportgraphics','builtin') == 5
    exportgraphics(hf, png, 'Resolution', 150);
else
    set(hf,'PaperPositionMode','auto','InvertHardcopy','off');
    print(hf, png, '-dpng', '-r150');
end
d  = dir(png);
nb = 0;  if ~isempty(d), nb = d.bytes; end
if nb < 20000
    warning('stage_backlash:exportStub', ...
        ['figure PNG is %d bytes -- likely an OpenGL blank stub. The .mat has ' ...
         'the data; re-plot from it.'], nb);
end
end
