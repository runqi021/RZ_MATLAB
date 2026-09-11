function acq_pause_queue()

ctrl_get();

f = findall(0,'Type','figure','Tag','rz_acq_pause_queue');
if ~isempty(f), figure(f(1)); return; end

BG = [0.13 0.13 0.13]; FG = [0.92 0.92 0.92];
YELLOW = [0.95 0.78 0.15]; GREEN = [0.29 0.70 0.36];
msg = '';
f = figure('Name','acq control','Tag','rz_acq_pause_queue','NumberTitle','off', ...
    'Units','pixels','Position',[80 300 310 460],'Color',BG, ...
    'MenuBar','none','ToolBar','none','Resize','off');

    function h = lab(y,txt,w)
        h = uicontrol(f,'Style','text','Units','pixels','Position',[10 y w 18], ...
            'String',txt,'HorizontalAlignment','left', ...
            'ForegroundColor',FG,'BackgroundColor',BG);
    end

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
defBG = get(hLock,'BackgroundColor');

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
        fprintf('ABORT sent -- closing the control window.\n');
        close(f);
    end

    function tf = noRun(S)
        tf = ~S.running;
        if tf
            say('no run in progress -- these do nothing. Start auto_acq_260910.');
        end
    end

    function say(t2), msg = t2; refresh(); end

    function s = stateLine(S)
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
