function acq_pause_queue()
%ACQ_PAUSE_QUEUE  Live control window for auto_acq_260910.
%
% A memo pad for the run: see what pauses are scheduled, delete them, add more,
% stop now, brake out of a tile in progress, continue, abort.
%
% Deliberately tiny.  The acquisition loop runs on MATLAB's main thread, the same
% one that draws ScanImage, so every widget here competes with both.
%
% HOW IT TALKS TO THE LOOP.  Through root appdata ('rz_acq_ctrl') -- not a global
% and not a handle either side has to hold.  The loop polls it; this window only
% writes it.  So they can start in either order, this window can be closed and
% reopened mid-run, and a crash in one cannot take the other down.
%
% THE THREE WAYS TO STOP, which differ in what happens to the current tile:
%   PAUSE NOW   finish the tile being acquired, then stop before the next one.
%               Nothing is lost.
%   BRAKE       abort the tile IN PROGRESS immediately (hSI.abort).  That tile's
%               file is incomplete, so CONTINUE RE-RUNS IT.
%   ABORT RUN   stop the tile in progress too, then end the script.  Like BRAKE,
%               it leaves that tile's file incomplete.  It ENDS the run, and this
%               window closes with it -- there is nothing left to control.
%               Starting again means running auto_acq_260910 again, which starts
%               over at tile 1 with the stage's current position as the origin.
%
% THE LOCK guards those four buttons and nothing else -- the queue stays editable
% while locked, because editing it costs nothing.  It starts LOCKED, since this
% window sits open for the length of a map and a stray click on BRAKE throws away
% a tile.  After that it stays wherever you put it: unlocking is not undone by
% using a button.
%
% Runqi Zhang / 2026-09-10

ctrl_get();                                  % create the shared state if absent

f = findall(0,'Type','figure','Tag','rz_acq_pause_queue');
if ~isempty(f), figure(f(1)); return; end

BG = [0.13 0.13 0.13]; FG = [0.92 0.92 0.92];
YELLOW = [0.95 0.78 0.15]; GREEN = [0.29 0.70 0.36];
msg = '';                                    % last thing clicked, line 2
f = figure('Name','acq control','Tag','rz_acq_pause_queue','NumberTitle','off', ...
    'Units','pixels','Position',[80 300 310 460],'Color',BG, ...
    'MenuBar','none','ToolBar','none','Resize','off');

    function h = lab(y,txt,w)
        h = uicontrol(f,'Style','text','Units','pixels','Position',[10 y w 18], ...
            'String',txt,'HorizontalAlignment','left', ...
            'ForegroundColor',FG,'BackgroundColor',BG);
    end

% Reads 0 until the first tile starts, so 0 means "has not run yet" rather than
% "finished tile 0".  It is the tile being worked on -- moved to, acquired, or
% parked before -- not the last one completed.
lab(430,'current tile',100);
hCur = uicontrol(f,'Style','edit','Units','pixels','Position',[110 430 60 22], ...
    'String','0','Enable','inactive','HorizontalAlignment','center', ...
    'FontWeight','bold');
hOf = lab(430,'of --',115); set(hOf,'Position',[180 432 115 18]);

lab(400,'pause at tile #',100);
hEd = uicontrol(f,'Style','edit','Units','pixels','Position',[110 400 60 22],'String','');
uicontrol(f,'Style','pushbutton','Units','pixels','Position',[180 400 115 22], ...
    'String','Add pause','Callback',@(~,~) addOne());

lab(372,'every',45);
hEvery = uicontrol(f,'Style','edit','Units','pixels','Position',[52 372 38 22],'String','15');
lab2 = lab(372,'up to',40); set(lab2,'Position',[96 372 40 18]);
% How far to schedule.  Auto-fills from the run once one is going; typed in by
% hand otherwise, so the queue can be laid out BEFORE the script starts.
hUpto = uicontrol(f,'Style','edit','Units','pixels','Position',[138 372 40 22],'String','81');
uicontrol(f,'Style','pushbutton','Units','pixels','Position',[186 372 109 22], ...
    'String','Schedule every N','Callback',@(~,~) addEvery());

hList = uicontrol(f,'Style','listbox','Units','pixels','Position',[10 186 285 176], ...
    'String',{},'Max',2,'Min',0);
uicontrol(f,'Style','pushbutton','Units','pixels','Position',[10 154 140 24], ...
    'String','Delete selected','Callback',@(~,~) delSel());
uicontrol(f,'Style','pushbutton','Units','pixels','Position',[155 154 140 24], ...
    'String','Clear all','Callback',@(~,~) clearAll());

hLock = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[10 126 140 22], ...
    'String','LOCK','Callback',@(~,~) setLock(true));
hUnlock = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[155 126 140 22], ...
    'String','UNLOCK','Callback',@(~,~) setLock(false));
defBG = get(hLock,'BackgroundColor');        % the theme's button face

hPause = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[10 86 140 36], ...
    'String','PAUSE NOW','FontWeight','bold', ...
    'BackgroundColor',YELLOW,'ForegroundColor',[0 0 0], ...
    'Callback',@(~,~) setFlag('pauseNow'));
hCont = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[155 86 140 36], ...
    'String','CONTINUE','FontWeight','bold', ...
    'BackgroundColor',GREEN,'ForegroundColor',[0 0 0], ...
    'Callback',@(~,~) doContinue());

hStat = uicontrol(f,'Style','text','Units','pixels','Position',[10 34 285 44], ...
    'String','','HorizontalAlignment','left','FontSize',8, ...
    'ForegroundColor',[.6 .8 1],'BackgroundColor',BG);

% Destructive pair, half height, at the very bottom -- out of the way of the two
% buttons actually used tile to tile.
hBrake = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[10 8 140 18], ...
    'String','BRAKE','FontSize',8,'FontWeight','bold','ForegroundColor',[.7 0 0], ...
    'Callback',@(~,~) doBrake());
hAbort = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[155 8 140 18], ...
    'String','ABORT RUN','FontSize',8,'Callback',@(~,~) doAbort());

t = timer('ExecutionMode','fixedSpacing','Period',0.4,'TimerFcn',@(~,~) refresh());
set(f,'CloseRequestFcn',@(s,~) onClose(s,t));
start(t); setLock(true); refresh();

    function addOne()
        v = str2double(get(hEd,'String'));
        if ~isfinite(v) || v < 1 || v ~= round(v)
            say('tile # must be a positive integer'); return
        end
        addTiles(v); set(hEd,'String','');
    end

    function addEvery()
        N = str2double(get(hEvery,'String'));
        S = ctrl_get();
        if ~isfinite(N) || N < 1 || N ~= round(N), say('N must be a positive integer'); return; end
        % Prefer the live grid size; fall back to the typed one so the queue can
        % be built before the acquisition starts.
        if S.nTiles >= 1
            M = S.nTiles;
        else
            M = str2double(get(hUpto,'String'));
            if ~isfinite(M) || M < 1, say('set "up to" (total tiles) first'); return; end
            M = round(M);
        end
        if N > M, say(sprintf('N=%d is past the last tile (%d)', N, M)); return; end
        addTiles(N:N:M);
    end

    function addTiles(v)
        S = ctrl_get();
        new = setdiff(v(:).', S.pauseAt);
        S.pauseAt = sort([S.pauseAt(:).' new]);
        ctrl_set(S);
        if isempty(new), say('already queued');
        else,            say(sprintf('queued %s', mat2str(new))); end
        refresh();
    end

    function delSel()
        S = ctrl_get();
        k = get(hList,'Value'); k = k(k >= 1 & k <= numel(S.pauseAt));
        if isempty(k), say('nothing selected'); return; end
        gone = S.pauseAt(k); S.pauseAt(k) = []; ctrl_set(S);
        say(sprintf('removed %s', mat2str(gone)));
        set(hList,'Value',[]); refresh();
    end

    function clearAll()
        S = ctrl_get(); S.pauseAt = []; ctrl_set(S);
        say('queue cleared'); set(hList,'Value',[]); refresh();
    end

    function setLock(tf)
        % Enable is the lock -- no state variable to fall out of step with it.
        onoff = 'on'; if tf, onoff = 'off'; end
        set([hPause hCont hBrake hAbort],'Enable',onoff);
        if tf
            set(hLock,'BackgroundColor',[.55 .18 .18],'ForegroundColor',[1 1 1]);
            set(hUnlock,'BackgroundColor',defBG,'ForegroundColor',[0 0 0]);
        else
            set(hLock,'BackgroundColor',defBG,'ForegroundColor',[0 0 0]);
            set(hUnlock,'BackgroundColor',[.20 .45 .25],'ForegroundColor',[1 1 1]);
        end
    end

    function setFlag(name)
        S = ctrl_get();
        if noRun(S), return; end
        S.(name) = true; ctrl_set(S);
        say(sprintf('%s -- takes effect after the current tile', name));
    end

    function doBrake()
        S = ctrl_get();
        if noRun(S), return; end
        S.brake = true; S.pauseNow = true; S.continue_ = false;
        ctrl_set(S);
        say('BRAKE: stopping this tile now.');
    end

    function doContinue()
        S = ctrl_get();
        if noRun(S), return; end
        S.continue_ = true; S.pauseNow = false; S.brake = false;
        ctrl_set(S);
        say('continuing');
    end

    function doAbort()
        S = ctrl_get();
        if noRun(S), return; end
        S.abort = true; S.pauseNow = false; S.continue_ = true; S.brake = false;
        ctrl_set(S);
        % The flag lives in root appdata, not in this figure, so the loop still
        % sees it after the window is gone.  Nothing left here to control.
        fprintf('ABORT sent -- closing the control window.\n');
        close(f);
    end

    function tf = noRun(S)
        % These four buttons only set flags for a loop to read.  With no loop
        % running the click would vanish silently, which reads as a broken
        % button -- so say what is actually the case.
        tf = ~S.running;
        if tf
            say('no run in progress -- these do nothing. Start auto_acq_260910.');
        end
    end

    function say(t2), msg = t2; refresh(); end

    function s = stateLine(S)
        % What the run is doing and what to press -- the line that matters when
        % you walk back to the rig and the window has been sitting there.
        if S.running && S.paused && S.retry
            s = sprintf('BRAKED at %d -- UNLOCK, then CONTINUE to redo it', S.atTile);
        elseif S.running && S.paused
            s = sprintf('PAUSED before %d -- UNLOCK to continue', S.atTile);
        elseif ~S.running && S.abort
            s = sprintf('ABORTED at %d -- idle', S.stopTile);
        elseif ~S.running && S.stopTile > 0
            s = sprintf('finished %d tiles -- idle', S.nTiles);
        elseif S.running
            s = sprintf('running -- tile %d / %d', S.atTile, S.nTiles);
        else
            s = 'not running';
        end
    end

    function refresh()
        % Stop the timer if the window is gone.  delete(f) bypasses
        % CloseRequestFcn, so a timer that only returns here would keep firing
        % against a dead handle forever -- and once the .m file is edited, error
        % on every tick.
        if ~isvalid(f)
            try, stop(t); delete(t); catch, end
            return
        end
        S = ctrl_get();
        if S.nTiles >= 1, set(hUpto,'String',num2str(S.nTiles)); end
        set(hCur,'String',num2str(S.atTile));
        if S.nTiles >= 1
            set(hOf,'String',sprintf('of %d', S.nTiles));
        else
            set(hOf,'String','of --');
        end
        set(hStat,'String',{stateLine(S); msg});
        items = arrayfun(@(v) sprintf('  pause before tile %d', v), S.pauseAt, 'uni', 0);
        if isempty(items), items = {'  (nothing queued)'}; end
        set(hList,'String',items);
        set(f,'Name',sprintf('acq control  |  %s', stateLine(S)));
    end

    function onClose(s, tm)
        try, stop(tm); delete(tm); catch, end
        delete(s);
    end
end

%% ------------------------------------------------------------------
function S = ctrl_get()
S = getappdata(0, 'rz_acq_ctrl');
if isempty(S) || ~isstruct(S) || ~isfield(S,'pauseAt') || ~isfield(S,'stopTile')
    S = struct('pauseAt',[],'pauseNow',false,'continue_',false,'abort',false, ...
               'brake',false,'retry',false,'running',false,'paused',false, ...
               'tile',0,'atTile',0,'stopTile',0,'nTiles',0);
    setappdata(0,'rz_acq_ctrl',S);
end
end

function ctrl_set(S)
setappdata(0, 'rz_acq_ctrl', S);
end
