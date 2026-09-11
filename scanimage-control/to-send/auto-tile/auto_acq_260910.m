filePrefix  = '260909_ChAT_g8m_shiverer';
numCols     = 9;
numRows     = 9;
xStep       = 400;
yStep       = 400;

pauseEveryN = 15; settleSec   = 0.5;
nSubSteps   = 10; subPause    = 0;

%%
ResArmSec   = 2.0;  posTolUm  = 2.0;  resumeTolUm = 1.0;
laserPort   = 'COM4'; laserNm = 930;  laserGDD    = 12500;  

doLaser     = false;
if doLaser
    addpath('C:\Users\dklab\Desktop\LaserControl-master\code');
    rehash toolboxcache;
    laserControl.settings.readSettings;
    C = laserControl.chameleon(laserPort);
    C.turnOn(); C.openShutter();
    C.setWavelength(laserNm); waitForLaser(C);
    C.setGDD(laserGDD);       waitForLaser(C);
end

addpath(genpath('C:\Users\dklab\Desktop\RZ_ScanImage_script'));
hCtrl   = evalin('base','hSICtl');  hSI = hCtrl.hModel;  hMotors = hSI.hMotors;  hScan2D = hSI.hScan2D;

origin = hMotors.motorPosition;
x0 = origin(1); y0 = origin(2); z0 = origin(3);
fprintf('origin  x %.2f  y %.2f  z %.2f\n', x0, y0, z0);

prevKeepOn = hScan2D.keepResonantScannerOn;
hScan2D.keepResonantScannerOn = true;

nTiles = numCols*numRows;
ctrl_init(nTiles, pauseEveryN);
if ~isempty(pauseEveryN) && pauseEveryN > 0
    fprintf('scheduled pauses every %d tiles: %s\n', pauseEveryN, ...
            mat2str(pauseEveryN:pauseEveryN:nTiles));
end
fprintf('run acq_pause_queue in another window for live control\n');

fprintf('arming resonant scanner (%.1f s)...\n', ResArmSec);
waitAndServe(ResArmSec);

tile = 0;  tRun = tic;  xPrev = x0;  yPrev = y0;

for row = 0:(numRows-1)
    if mod(row,2) == 0
        colRange = 0:(numCols-1);
    else
        colRange = (numCols-1):-1:0;
    end

    for col = colRange
        tile = tile + 1;
        ctrl_at(tile);
        xTarget = x0 + col*xStep;
        yTarget = y0 + row*yStep;

        for s = 1:nSubSteps
            fr = s/nSubSteps;
            hMotors.moveStartRelative([xPrev + (xTarget-xPrev)*fr, ...
                                       yPrev + (yTarget-yPrev)*fr, z0]);
            hMotors.moveWaitForFinish([1 2 3]);
            if subPause > 0, waitAndServe(subPause); end
        end
        xPrev = xTarget;  yPrev = yTarget;

        pos = hMotors.motorPosition;
        err = hypot(pos(1)-xTarget, pos(2)-yTarget);
        if err > posTolUm
            warning('tile %d (col %d row %d): stage landed %.2f um from target', ...
                    tile, col, row, err);
        end

        fprintf('[%3d/%3d] col%02d row%02d  x %.0f y %.0f\n', ...
                tile, nTiles, col, row, pos(1), pos(2));
        waitAndServe(settleSec);

        acquired = false;
        while ~acquired
            if ctrl_should_pause(tile)
                hScan2D.keepResonantScannerOn = false;
                if ctrl_is_retry()
                    fprintf('   PAUSED (braked) at tile %d -- CONTINUE re-runs it\n', tile);
                else
                    fprintf('   PAUSED before tile %d -- press CONTINUE\n', tile);
                end
                ctrl_wait_for_continue(tile);
                if ctrl_aborted(), break; end

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
            if strcmp(hSI.acqState,'idle')
                hScan2D.keepResonantScannerOn = prevKeepOn;
                ctrl_finish();
                error('auto_acq:grabDidNotStart', ...
                    ['startGrab returned with acqState still ''idle'' at tile %d.\n' ...
                     'The grab never started -- most likely no channel is selected ' ...
                     'for Display or Save in ScanImage.\nFix that and re-run.'], tile);
            end
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
                    hScan2D.keepResonantScannerOn = false;
                    break
                end
                ctrl_set_retry();
                continue
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

function waitAndServe(sec)
t0 = tic;
while toc(t0) < sec
    pause(0.05);
    drawnow limitrate;
end
end

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
S.paused = true; S.atTile = tile;
S.continue_ = false; S.pauseNow = false; S.brake = false;
S.pauseAt(S.pauseAt == tile) = [];
setappdata(0,'rz_acq_ctrl',S);
while true
    S = getappdata(0,'rz_acq_ctrl');
    if S.continue_ || S.abort, break; end
    pause(0.1); drawnow;
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
