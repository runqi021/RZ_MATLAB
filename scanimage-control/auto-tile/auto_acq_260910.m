%% auto_acq_260910.m  --  snake-tile map acquisition
%  Successor to auto_acq_n_shutterOFF.m.  Same grid, same snake, same filename
%  convention (which the stitcher parses), same laser bring-up.  The laser is
%  deliberately left ON and the shutter OPEN at the end, as in the original.
%
%  Run `acq_pause_queue` in another window for live control.
%
% ============================ WHAT CHANGED ============================
% 1. IDLE TIME PER TILE: 10 s -> settleSec (0.5 s), ~13 min back on an 81-tile
%    map.  The old loop did `pause(1)` INSIDE the 10-sub-step move loop, AFTER
%    moveWaitForFinish had already returned -- ten 1 s stops interleaved with the
%    stage stepping in 40 um hops, all of it dead time, with the objective in
%    motion for most of it.  Now the move runs uninterrupted and there is a
%    single settleSec wait at the target with the stage STILL.
%
%    The 10 sub-moves are KEPT, unchanged from the original.  NOT because of
%    motorFastMotionThreshold -- that does nothing on this rig: the Galil's
%    registry entry sets TwoStep.Enable = false
%    (+scanimage/+components/+motors/+MotorRegistry/galil_DMC4040.m:9), so the
%    fast/slow velocity switch is disabled, and the MDF leaves velocitySlow and
%    velocityFast empty in any case.  What the sub-moves actually do is make the
%    traverse ten short accel-decel cycles instead of one continuous 400 um run,
%    which is gentler on the prep.  Kept because the original did it and worked.
%
% 2. ONE MOTOR POSITION QUERY PER TILE INSTEAD OF ELEVEN.  get.motorPosition
%    issues a live `TP` query to the Galil over RS232 and moveStartRelative calls
%    it internally every time -- ScanImage's own source comments that this
%    "might be problematic because we increase numbers of commands sent to motor
%    controller".  Sub-step targets are now interpolated from the previous
%    target, tracked in a variable.  The port is one ScanImage also polls.
%
% 3. EVERY WAIT YIELDS (drawnow), so the ScanImage GUI stays live and the control
%    window stays responsive instead of freezing for the duration.
%
% 4. RESONANT SCANNER held on across the run, spun down while paused, re-armed
%    ResArmSec before acquisition resumes -- including before tile 1.
%
% 5. LIVE CONTROL: scheduled pauses, pause-now, brake, abort.  See below.
%
% 6. ARRIVAL IS CHECKED, and position is re-checked after every pause.
%
% ====================== THE THREE WAYS TO STOP ======================
%   PAUSE (scheduled or PAUSE NOW)  stops at a tile BOUNDARY: the tile being
%       acquired finishes first.  "pause before tile 34" means the stage parks at
%       tile 34, you do what you need, CONTINUE, and 34 is then acquired.
%   BRAKE  aborts the tile IN PROGRESS via hSI.abort().  That tile's file is
%       incomplete, so CONTINUE RE-RUNS IT.
%   ABORT RUN  also stops the tile in progress, and then ends the script.  Same
%       stop as BRAKE, different thing after it -- so it too leaves one short
%       file behind, which is the tile to re-shoot or delete.
%
% DURING A PAUSE the ScanImage GUI is fully usable -- the wait loop yields to the
% event queue, so buttons, focus and the motor controls all work.  The MATLAB
% COMMAND PROMPT is not, because this script still owns the thread.  If you move
% the stage while paused, the resume checks position and drives back to the
% planned target before acquiring.
%
% Runqi Zhang / 2026-09-10

%% ========================= CONFIGURATION =========================
numCols     = 9;        % tiles along X
numRows     = 9;        % tiles along Y
xStep       = 400;      % um
yStep       = 400;      % um
filePrefix  = '260909_ChAT_g8m_shiverer';

pauseEveryN = 15;       % pre-schedule a pause every N tiles.  [] or 0 = none.
                        % Edit the queue freely in acq_pause_queue afterwards.

settleSec   = 0.5;      % stationary wait at the target before the grab
nSubSteps   = 10;       % sub-moves per tile; ten short accel-decel cycles rather
                        % than one 400 um traverse.  See the note above: this is
                        % NOT about motorFastMotionThreshold, which is inert here.
subPause    = 0;        % extra wait between sub-moves; the controller already
                        % reports arrival before it would run
ResArmSec   = 2.0;      % resonant scanner spin-up allowed before each acquisition
posTolUm    = 2.0;      % warn if the stage lands further than this from target
resumeTolUm = 1.0;      % after a pause, drive back if it moved more than this

laserPort   = 'COM4';
laserNm     = 930;
laserGDD    = 12500;
doLaser     = true;     % false if the laser is already up

%% ========================= LASER =========================
if doLaser
    addpath('C:\Users\dklab\Desktop\LaserControl-master\code');
    rehash toolboxcache;
    laserControl.settings.readSettings;
    C = laserControl.chameleon(laserPort);
    C.turnOn(); C.openShutter();
    C.setWavelength(laserNm); waitForLaser(C);
    C.setGDD(laserGDD);       waitForLaser(C);
end

%% ========================= INIT =========================
addpath(genpath('C:\Users\dklab\Desktop\RZ_ScanImage_script'));
hCtrl   = evalin('base','hSICtl');
hSI     = hCtrl.hModel;
hMotors = hSI.hMotors;
hScan2D = hSI.hScan2D;

origin = hMotors.motorPosition;          % the one query needed up front
x0 = origin(1); y0 = origin(2); z0 = origin(3);
fprintf('origin  x %.2f  y %.2f  z %.2f\n', x0, y0, z0);

% NOT onCleanup: in a SCRIPT the cleanup object lands in the base workspace and
% does not run until cleared, so it would restore nothing.  Restored explicitly
% after the loop.  An error mid-run leaves the scanner spinning, which is benign.
prevKeepOn = hScan2D.keepResonantScannerOn;
hScan2D.keepResonantScannerOn = true;

nTiles = numCols*numRows;
ctrl_init(nTiles, pauseEveryN);

% The control window is opened HERE, not left to the user. The two files were
% written as independent scripts and that split was a trap: the run would reach
% a scheduled pause with no window on screen and nothing to press CONTINUE with.
ctrl_open_gui();

% Print the ACTUAL queue, not the requested schedule. ctrl_init MERGES the new
% every-N schedule into whatever pauseAt survived the last run, so the two can
% differ -- printing the request was misleading.
qAt = ctrl_queue();
if isempty(qAt)
    fprintf('pause queue: none\n');
else
    fprintf('pause queue: %s\n', mat2str(qAt));
end

fprintf('arming resonant scanner (%.1f s)...\n', ResArmSec);
waitAndServe(ResArmSec);

%% ========================= TILING LOOP =========================
tile  = 0;  tRun = tic;
xPrev = x0; yPrev = y0;          % where the stage is: tracked, not queried

for row = 0:(numRows-1)
    if mod(row,2) == 0
        colRange = 0:(numCols-1);
    else
        colRange = (numCols-1):-1:0;
    end

    for col = colRange
        tile = tile + 1;
        ctrl_at(tile);                  % window shows it from the moment it starts
        xTarget = x0 + col*xStep;
        yTarget = y0 + row*yStep;

        % ---- move, sub-stepped to stay slow, uninterrupted ----
        for s = 1:nSubSteps
            fr = s/nSubSteps;
            hMotors.moveStartRelative([xPrev + (xTarget-xPrev)*fr, ...
                                       yPrev + (yTarget-yPrev)*fr, z0]);
            hMotors.moveWaitForFinish([1 2 3]);
            if subPause > 0, waitAndServe(subPause); end
        end
        xPrev = xTarget;  yPrev = yTarget;

        pos = hMotors.motorPosition;         % the one query per tile
        err = hypot(pos(1)-xTarget, pos(2)-yTarget);
        if err > posTolUm
            warning('tile %d (col %d row %d): stage landed %.2f um from target', ...
                    tile, col, row, err);
        end

        fprintf('[%3d/%3d] col%02d row%02d  x %.0f y %.0f\n', ...
                tile, nTiles, col, row, pos(1), pos(2));
        waitAndServe(settleSec);

        % ---- acquire, with pause/brake handling and retry ----
        % A braked tile is re-run rather than skipped: hSI.abort() leaves a
        % truncated file on disk that would otherwise look like a complete tile
        % to the stitcher, with fewer z pages than the header claims.
        acquired = false;
        while ~acquired
            if ctrl_should_pause(tile)
                hScan2D.keepResonantScannerOn = false;   % spin down while stopped
                if ctrl_is_retry()
                    fprintf('   PAUSED (braked) at tile %d -- CONTINUE re-runs it\n', tile);
                else
                    fprintf('   PAUSED before tile %d -- press CONTINUE\n', tile);
                end
                ctrl_wait_for_continue(tile);
                if ctrl_aborted(), break; end

                % The stage may have been driven by hand during the pause, and
                % the filename is written from wherever it actually is -- so a
                % tile acquired after an unnoticed nudge would be placed wrong
                % and nothing downstream could tell.  Check, and drive back.
                pos = hMotors.motorPosition;
                d = hypot(pos(1)-xTarget, pos(2)-yTarget) + abs(pos(3)-z0);
                if d > resumeTolUm
                    fprintf('   stage moved %.2f um during the pause -- returning\n', d);
                    hMotors.moveStartRelative([xTarget, yTarget, z0]);
                    hMotors.moveWaitForFinish([1 2 3]);
                    pos = hMotors.motorPosition;
                    fprintf('   back at x %.2f y %.2f z %.2f\n', pos(1), pos(2), pos(3));
                end
                xPrev = pos(1); yPrev = pos(2);

                hScan2D.keepResonantScannerOn = true;
                fprintf('   re-arming scanner (%.1f s)...\n', ResArmSec);
                waitAndServe(ResArmSec);
                ctrl_clear_retry();
            end
            if ctrl_aborted(), break; end

            hScan2D.logFileStem = sprintf('%s_col%02d_row%02d_x%.0f_y%.0f', ...
                                          filePrefix, col, row, pos(1), pos(2));
            hSI.startGrab();
            % SI.start() sets acqState synchronously before it does any hardware
            % work, so by here it is 'grab' -- EXCEPT on the one path that
            % returns early without setting it: no channel selected for Display
            % or Save.  Then the wait loop below would fall straight through and
            % the script would march the stage over all 81 tiles acquiring
            % nothing, silently.  (The original script has the same hole.)
            if strcmp(hSI.acqState,'idle')
                hScan2D.keepResonantScannerOn = prevKeepOn;
                ctrl_finish();
                error('auto_acq:grabDidNotStart', ...
                    ['startGrab returned with acqState still ''idle'' at tile %d.\n' ...
                     'The grab never started -- most likely no channel is selected ' ...
                     'for Display or Save in ScanImage.\nFix that and re-run.'], tile);
            end
            % BRAKE and ABORT both stop the grab where it stands and differ only
            % in what follows.  Either way that tile's file is short, so it is
            % named in the log -- it is the one to re-shoot or delete.
            stopped = '';
            while ~strcmp(hSI.acqState,'idle')
                if ctrl_aborted(), stopped = 'ABORT'; end
                if ctrl_braked(),  stopped = 'BRAKE'; end
                if ~isempty(stopped)
                    fprintf('   %s -- stopping tile %d mid-acquisition (%s is incomplete)\n', ...
                            stopped, tile, hScan2D.logFileStem);
                    hSI.abort();
                    break
                end
                pause(0.05); drawnow limitrate;
            end
            if ~isempty(stopped)
                while ~strcmp(hSI.acqState,'idle'), pause(0.05); drawnow limitrate; end
                if strcmp(stopped,'ABORT')
                    hScan2D.keepResonantScannerOn = false;   % idle, right here
                    break
                end
                ctrl_set_retry();           % CONTINUE comes back to this tile
                continue                    % round again, into the pause branch
            end
            acquired = true;
            ctrl_progress(tile);
        end
        if ctrl_aborted(), break; end
    end
    if ctrl_aborted(), break; end
end

if ctrl_aborted()
    hScan2D.keepResonantScannerOn = false;
    fprintf('\n[ABORTED] at tile %d of %d after %.1f min -- idle\n', tile, nTiles, toc(tRun)/60);
else
    hScan2D.keepResonantScannerOn = prevKeepOn;
    fprintf('\n[DONE] %d of %d tiles in %.1f min\n', tile, nTiles, toc(tRun)/60);
end
ctrl_finish();

%% ========================= HELPERS =========================
function waitAndServe(sec)
% Wait while yielding, so the ScanImage GUI redraws and the control window stays
% responsive.  A bare pause(sec) blocks both for the whole interval.
t0 = tic;
while toc(t0) < sec
    pause(0.05);
    drawnow limitrate;
end
end

%% ---- shared control state, written by acq_pause_queue ----
function ctrl_init(nTiles, everyN)
S = getappdata(0,'rz_acq_ctrl');
if isempty(S) || ~isstruct(S) || ~isfield(S,'pauseAt') || ~isfield(S,'stopTile')
    S = struct('pauseAt',[],'pauseNow',false,'continue_',false,'abort',false, ...
               'brake',false,'retry',false,'running',false,'paused',false, ...
               'tile',0,'atTile',0,'stopTile',0,'nTiles',0);
end
S.running = true; S.paused = false; S.abort = false; S.brake = false;
S.retry = false; S.continue_ = false; S.pauseNow = false;
S.tile = 0; S.atTile = 0; S.stopTile = 0; S.nTiles = nTiles;
if ~isempty(everyN) && everyN > 0
    S.pauseAt = unique([S.pauseAt(:).' everyN:everyN:nTiles]);
end
setappdata(0,'rz_acq_ctrl',S);
end

function ctrl_open_gui()
% Open the live-control window from the run itself.
%
% Resolution order matters. acq_pause_queue.m sits NEXT TO this file, and the
% rig's current folder is not guaranteed to be that folder, so put this file's
% own directory on the path first. mfilename('fullpath') is valid inside a local
% function of a script (verified on R2024a), and returns '' only if the script
% was pasted rather than run by name -- hence the isempty guard.
here = fileparts(mfilename('fullpath'));
if ~isempty(here) && exist(fullfile(here,'acq_pause_queue.m'),'file') == 2
    addpath(here);
end
if exist('acq_pause_queue','file') ~= 2
    warning('auto_acq:noCtrlWindow', ...
        'acq_pause_queue.m not found -- running with NO live control');
    return
end
% A figure failure must never kill an acquisition that is otherwise fine: the
% run degrades to the un-pausable behaviour instead of erroring out mid-map.
try
    acq_pause_queue();
    drawnow;
catch ME
    warning('auto_acq:ctrlWindowFailed', ...
        'control window did not open (%s) -- running with NO live control', ME.message);
end
end

function q = ctrl_queue()
S = getappdata(0,'rz_acq_ctrl');
q = S.pauseAt;
end

function tf = ctrl_should_pause(tile)
S = getappdata(0,'rz_acq_ctrl');
tf = S.pauseNow || S.retry || ismember(tile, S.pauseAt);
end

function tf = ctrl_braked()
S = getappdata(0,'rz_acq_ctrl');
tf = isfield(S,'brake') && S.brake;
end

function ctrl_set_retry()
S = getappdata(0,'rz_acq_ctrl');
S.retry = true; S.brake = false; S.pauseNow = true; S.continue_ = false;
setappdata(0,'rz_acq_ctrl',S);
end

function tf = ctrl_is_retry()
S = getappdata(0,'rz_acq_ctrl');
tf = isfield(S,'retry') && S.retry;
end

function ctrl_clear_retry()
S = getappdata(0,'rz_acq_ctrl'); S.retry = false; setappdata(0,'rz_acq_ctrl',S);
end

function ctrl_wait_for_continue(tile)
S = getappdata(0,'rz_acq_ctrl');
S.paused = true; S.atTile = tile;          % the tile it is parked BEFORE
S.continue_ = false; S.pauseNow = false; S.brake = false;
S.pauseAt(S.pauseAt == tile) = [];      % a queued pause fires once
setappdata(0,'rz_acq_ctrl',S);
while true
    S = getappdata(0,'rz_acq_ctrl');
    if S.continue_ || S.abort, break; end
    pause(0.1); drawnow;                % generous yield: ScanImage stays usable
end
S.paused = false; S.continue_ = false;
setappdata(0,'rz_acq_ctrl',S);
end

function tf = ctrl_aborted()
S = getappdata(0,'rz_acq_ctrl');
tf = isfield(S,'abort') && S.abort;
end

function ctrl_at(tile)
S = getappdata(0,'rz_acq_ctrl'); S.atTile = tile; setappdata(0,'rz_acq_ctrl',S);
end

function ctrl_progress(tile)
S = getappdata(0,'rz_acq_ctrl'); S.tile = tile; setappdata(0,'rz_acq_ctrl',S);
end

function ctrl_finish()
% Back to 0 when nothing is running, so the box only ever shows a tile that is
% actually being worked on.  Where it stopped is kept for the status line.
S = getappdata(0,'rz_acq_ctrl');
S.running = false; S.paused = false;
S.stopTile = S.atTile; S.atTile = 0;
setappdata(0,'rz_acq_ctrl',S);
end

function waitForLaser(C)
while true
    [ready, msg] = C.isReady();
    tuning = C.isTuning();
    lambda = C.readWavelength();
    if ismethod(C,'readGDD'), gdd = C.readGDD(); else, gdd = NaN; end
    if ready && ~tuning
        fprintf('Laser ready, %.1f nm', lambda);
        if ~isnan(gdd), fprintf(', GDD %.0f fs^2', gdd); end
        fprintf('\n'); return
    end
    if tuning
        fprintf('  waiting: tuning (%.1f nm)\n', lambda);
    else
        fprintf('  waiting: %s (%.1f nm)\n', msg, lambda);
    end
    pause(1); drawnow limitrate;
end
end
