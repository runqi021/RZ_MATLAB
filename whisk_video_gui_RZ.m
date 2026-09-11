function whisk_video_gui_RZ()
% whisk_video_gui_RZ  Scrub the whisker video alongside L/R whisk-angle traces,
% frame-locked, and export a clip in the same layout.
%   Top    : video frame (DLC *_labeled.mp4; raw .avi is FFV1 = unreadable here)
%   Bottom : L (green) / R (blue) BP whisk angle with a playhead at the frame
% Navigate sessions with the dropdown / Prev / Next. Type a start and duration
% and press "Export clip" to render that chunk (video + trace) to an mp4.
% Frame k in the video == row k in the DLC csv == xL(k)/xR(k).

% ============================ USER-EDITABLE ============================
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BP       = [6 30];      % fast-whisking bandpass (Hz)
fpsW     = 400;         % whisk camera fps
WINSEC   = 1;         % trace window shown around the playhead (s)
T0       = 39;          % initial playhead time (s)
EXPORT_FPS = 100;       % exported-clip playback fps (the "speed" box sets slow-mo: 0.25 = 4x slower)
EXPORT_DIR = "";        % output folder ("" = this script's folder)
EXPORT_WIN = 1.0;       % exported trace: rolling window width (s), like the GUI
TRACE_TALL = 1.4;       % exported trace panel height factor (vertical; 1 = default)
% --- export layout knobs (tweak to reposition; all in the exported frame) ---
EXP_BASE_H  = 240;      % base trace-panel height in px (final height = EXP_BASE_H * TRACE_TALL)
EXP_PPS     = 600;      % rolling strip resolution (px per second of trace)
EXP_YAXIS_FRAC = 0.10;  % y-axis (label+ticks) width as a fraction of the VISIBLE trace panel; raise for more room
EXP_RMARGIN = 0.01;     % RIGHT margin (fraction)
EXP_AXBOT   = 0.18;     % axes BOTTOM (fraction) = room for x tick numbers; raise for more room
EXP_AXTOP   = 0.05;     % axes TOP margin (fraction)
EXP_CAP_H   = 30;       % 'Time (s)' caption band height px (0 = no caption); font scales with it
EXP_SIDEGAP = 0.03;     % white gap on LEFT & RIGHT of the trace panel (fraction of panel width)
EXP_CAP_Y   = 0.93;     % 'Time (s)' vertical center as fraction of trace height (SMALLER = higher / toward top)
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');
cL = [0 0.55 0]; cR = [0 0.4 0.85];

% ---- enumerate sessions (one per readable *_labeled.mp4) ----
sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_labeled.mp4 in %s', whiskDir);
labels = cellfun(@(c) sprintf('%s  n%d', c{1}, c{2}), sess, 'UniformOutput',false);

% ---- state (filled by loadSession) ----
v=[]; xL=[]; xR=[]; t=[]; n=0; cf=max(1,round(T0*fpsW)+1); yl=[0 1];
animal=""; kRun=0; tmr=[]; si=1; loading=false;
him=[]; pL=[]; pR=[]; hPlay=[];
PERIOD = round(1000/60)/1000;   % timer period (s), ms-rounded to avoid precision warning

% ---- UI ----
f = figure('Color','w','Position',[80 60 1180 800],'Name','whisk video + angle', ...
    'NumberTitle','off','CloseRequestFcn',@onClose);
axV = axes(f,'Units','normalized','Position',[0.05 0.46 0.90 0.50]);
axT = axes(f,'Units','normalized','Position',[0.07 0.18 0.88 0.21]);

% nav row
uicontrol(f,'Style','text','Units','normalized','Position',[0.07 0.108 0.05 0.026], ...
    'String','session','BackgroundColor','w','HorizontalAlignment','left');
pop = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.12 0.10 0.18 0.04], ...
    'String',labels,'Value',1,'Callback',@onPick);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.31 0.10 0.05 0.04], ...
    'String','< Prev','Callback',@(s,e)onStep(-1));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.36 0.10 0.05 0.04], ...
    'String','Next >','Callback',@(s,e)onStep(+1));

% transport row
sld = uicontrol(f,'Style','slider','Units','normalized','Position',[0.07 0.065 0.52 0.028], ...
    'Min',1,'Max',2,'Value',1,'Callback',@onSlide);
btn = uicontrol(f,'Style','togglebutton','Units','normalized','Position',[0.61 0.06 0.07 0.038], ...
    'String','Play','FontWeight','bold','Callback',@onPlay);
uicontrol(f,'Style','text','Units','normalized','Position',[0.69 0.063 0.05 0.026], ...
    'String','speed','BackgroundColor','w','HorizontalAlignment','left');
spd = uicontrol(f,'Style','edit','Units','normalized','Position',[0.74 0.06 0.05 0.038],'String','1');

% export row
uicontrol(f,'Style','text','Units','normalized','Position',[0.07 0.022 0.055 0.026], ...
    'String','start (s)','BackgroundColor','w','HorizontalAlignment','left');
eStart = uicontrol(f,'Style','edit','Units','normalized','Position',[0.125 0.018 0.05 0.04],'String',num2str(T0));
uicontrol(f,'Style','text','Units','normalized','Position',[0.185 0.022 0.05 0.026], ...
    'String','dur (s)','BackgroundColor','w','HorizontalAlignment','left');
eDur = uicontrol(f,'Style','edit','Units','normalized','Position',[0.235 0.018 0.05 0.04],'String','2');
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.30 0.016 0.11 0.044], ...
    'String','Export clip','FontWeight','bold','Callback',@onExport);
txt = uicontrol(f,'Style','text','Units','normalized','Position',[0.43 0.02 0.55 0.03], ...
    'String','','BackgroundColor','w','HorizontalAlignment','left','FontName','Consolas');

loadSession(1);

% ================= session loading =================
    function loadSession(idx)
        if loading, return; end          % ignore reentrant calls (drawnow can fire callbacks)
        loading = true; cleanup = onCleanup(@() setLoading(false));
        stopplay(); set(btn,'Value',0,'String','Play');
        si = idx; animal = sess{idx}{1}; kRun = sess{idx}{2};
        M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', char(animal), kRun)), 0.6);
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));   % LEFT  (x mirrored)
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));    % RIGHT
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        t  = (0:numel(xL)-1)'/fpsW;
        v  = VideoReader(pick_labeled_mp4(whiskDir, sprintf('%s_whisk_n%d', char(animal), kRun))); %#ok<TNMLP>
        n  = min(v.NumFrames, numel(xL));
        cf = max(1, min(n, cf));
        % video
        cla(axV); him = imshow(read(v,cf),'Parent',axV);
        % traces
        cla(axT); hold(axT,'on'); grid(axT,'on');
        pL = plot(axT, t(1:n), xL(1:n), '-', 'Color',cL, 'LineWidth',0.9);
        pR = plot(axT, t(1:n), xR(1:n), '-', 'Color',cR, 'LineWidth',0.9);
        ylabel(axT,'Whisker angle (degree)'); xlabel(axT,'time (s)');
        legend([pL pR],{'L whisk','R whisk'},'Orientation','horizontal','Location','northoutside');
        yl = [min([xL(1:n);xR(1:n)]) max([xL(1:n);xR(1:n)])]; yl = yl + 0.05*diff(yl)*[-1 1];
        ylim(axT, yl);
        hPlay = line(axT,[t(cf) t(cf)], yl, 'Color',[0.85 0.1 0.1],'LineWidth',1.2,'HandleVisibility','off');
        set(sld,'Min',1,'Max',n,'Value',cf,'SliderStep',[1/(n-1) 30/(n-1)]);
        set(pop,'Value',si);
        redraw();
    end
    function setLoading(tf), loading = tf; end

% ================= navigation =================
    function onPick(src,~), loadSession(get(src,'Value')); end
    function onStep(d)
        idx = min(numel(sess), max(1, si + d));
        if idx~=si, loadSession(idx); end
    end

% ================= transport =================
    function onSlide(~,~), stopplay(); cf = round(get(sld,'Value')); redraw(); end
    function onPlay(src,~)
        if get(src,'Value')
            set(src,'String','Pause');
            v.CurrentTime = max(0,(cf-1)/fpsW);
            tmr = timer('ExecutionMode','fixedRate','Period',PERIOD, ...
                'BusyMode','drop','TimerFcn',@tick);
            start(tmr);
        else
            stopplay();
        end
    end
    function tick(~,~)
        if loading || isempty(him) || ~isgraphics(him), return; end
        s = str2double(get(spd,'String')); if isnan(s)||s<=0, s=1; end
        step = max(1, round(s*fpsW*PERIOD));
        img = [];
        for i = 1:step
            if cf < n && hasFrame(v), img = readFrame(v); cf = cf + 1; else, break; end
        end
        if isempty(img) || cf >= n
            if ~isempty(img), showFrame(img); end
            stopplay(); set(btn,'Value',0,'String','Play'); return;
        end
        showFrame(img);
    end
    function stopplay()
        if ~isempty(tmr) && isvalid(tmr), stop(tmr); delete(tmr); end
        tmr = [];
    end
    function showFrame(img)            % playback: decoded frame, no seek, no limitrate
        if isempty(him) || ~isgraphics(him), return; end
        cf = max(1, min(n, cf)); set(him,'CData', img); overlay(); drawnow;
    end
    function redraw()                  % scrub: random-access read + throttled draw
        if isempty(him) || ~isgraphics(him), return; end
        cf = max(1, min(n, cf)); set(him,'CData', read(v,cf)); overlay(); drawnow limitrate;
    end
    function overlay()
        if ~isgraphics(hPlay), return; end
        tc = t(cf); set(hPlay,'XData',[tc tc]); xlim(axT, tc + WINSEC*[-0.5 0.5]);
        set(sld,'Value',cf);
        set(txt,'String', sprintf('frame %d / %d   t=%.3f s   xL=%+6.2f   xR=%+6.2f', cf, n, tc, xL(cf), xR(cf)));
    end

% ================= export =================
    function onExport(~,~)
        t0 = str2double(get(eStart,'String')); dur = str2double(get(eDur,'String'));
        if isnan(t0)||isnan(dur)||dur<=0, errordlg('Enter numeric start and positive duration.','export'); return; end
        stopplay(); set(btn,'Value',0,'String','Play');
        f0 = max(1, round(t0*fpsW)+1); f1 = min(n, f0 + round(dur*fpsW) - 1);
        if f1 <= f0, errordlg('Requested range is outside this recording.','export'); return; end
        s = str2double(get(spd,'String')); if isnan(s)||s<=0, s=1; end   % speed box: 0.1 = 10x slower
        dstep  = max(1, round(fpsW*s/EXPORT_FPS));   % drop frames only if needed to stay near EXPORT_FPS
        outfps = max(1, (fpsW/dstep)*s);             % => clip duration = dur/s exactly (any speed)
        if strlength(EXPORT_DIR)>0, od = char(EXPORT_DIR); else, od = fileparts(mfilename('fullpath')); end
        outfile = fullfile(od, sprintf('%s_n%d_t%05.1f-%05.1fs_%gx_traceclip.mp4', char(animal), kRun, t0, t0+dur, s));
        vw = VideoWriter(outfile,'MPEG-4'); vw.FrameRate = outfps; vw.Quality = 95; open(vw);
        % ROLLING trace like the GUI: pre-render the WHOLE clip as one wide strip
        % (with the y-axis on the left), then per frame crop the EXPORT_WIN window
        % and paste it next to the static y-axis. No per-frame getframe (fast).
        Htr = round(EXP_BASE_H*TRACE_TALL); Htr = Htr - mod(Htr,2); PPS = EXP_PPS;  % trace panel height
        ws_t0 = t0 - EXPORT_WIN/2; ws_t1 = t0 + dur + EXPORT_WIN/2;       % render with half-window margins
        i0 = max(1, round(ws_t0*fpsW)+1); i1 = min(n, round(ws_t1*fpsW)+1);
        Rmn = EXP_RMARGIN;                                               % right margin (fraction of strip)
        winFrac = EXPORT_WIN/(ws_t1-ws_t0);                              % visible window / full rendered span
        % choose strip left margin so the y-axis is EXP_YAXIS_FRAC of the VISIBLE panel after cropping
        Lmn = EXP_YAXIS_FRAC*winFrac*(1-Rmn) / (1 - EXP_YAXIS_FRAC + EXP_YAXIS_FRAC*winFrac);
        Lmn = min(max(Lmn,0.005),0.30);
        Wfig = max(400, round((ws_t1-ws_t0)*PPS));
        tf = figure('Color','w','Position',[60 60 Wfig EXP_BASE_H],'Name','rendering trace...','NumberTitle','off');
        tax = axes(tf,'Units','normalized','Position',[Lmn EXP_AXBOT 1-Lmn-Rmn 1-EXP_AXBOT-EXP_AXTOP]); hold(tax,'on'); grid(tax,'on');
        plot(tax, t(i0:i1), xL(i0:i1), '-', 'Color',cL, 'LineWidth',1);
        plot(tax, t(i0:i1), xR(i0:i1), '-', 'Color',cR, 'LineWidth',1);
        xlim(tax,[ws_t0 ws_t1]); ylim(tax, yl);
        xlabel(tax,''); ylabel(tax,'Whisker angle (degree)');   % x label added later as a STATIC caption
        drawnow; frb = getframe(tf); strip = frb.cdata; axFont = get(tax,'FontName'); axFsz = get(tax,'FontSize'); delete(tf);
        Hs = size(strip,1); Ws = size(strip,2); Lpx = round(Lmn*Ws); plotW = Ws-Lpx-round(Rmn*Ws);
        colOf = @(tt) Lpx + round((tt-ws_t0)/(ws_t1-ws_t0)*plotW);       % time -> column in strip
        wcols = colOf(t0+EXPORT_WIN/2) - colOf(t0-EXPORT_WIN/2);         % fixed window width (px)
        gapL = max(2, round(EXP_SIDEGAP*(Lpx+wcols))); gapR = round(1.5*gapL);   % right gap 1.5x left
        Wc = gapL + gapR + Lpx + wcols + 1;                             % gapL | y-axis | window | gapR
        Wout = round(Htr * Wc / Hs); Wout = Wout - mod(Wout,2);          % width from crop aspect -> UNIFORM resize (no stretch)
        rTop = max(1,round(EXP_AXTOP*Hs)); rBot = min(Hs,round((1-EXP_AXBOT)*Hs));  % playhead spans the axes box only
        vH = round(v.Height*Wout/v.Width); vH = vH - mod(vH,2);          % video same width (aligned), aspect kept
        Hcap = EXP_CAP_H - mod(EXP_CAP_H,2);                             % caption text-band height (even)
        cap = []; capMask = []; capR0 = [];
        if Hcap > 0                                                     % render 'Time (s)' in MATLAB's font (matches axes)
            cf = figure('Color','w','Position',[60 60 Wout Hcap],'Name','caption','NumberTitle','off');
            ca = axes(cf,'Position',[0 0 1 1]); axis(ca,'off');
            text(ca,0.5,0.5,'Time (s)','Units','normalized','HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName',axFont,'FontSize',max(8,round(axFsz*Htr/EXP_BASE_H)));
            drawnow; cg = getframe(cf); cap = imresize(cg.cdata,[Hcap Wout]); delete(cf);
            capMask = repmat(any(cap<160,3),1,1,3);                      % text pixels only (don't paint over ticks)
            capR0 = min(Htr-Hcap+1, max(1, round(EXP_CAP_Y*Htr - Hcap/2)));   % top row inside the trace panel
        end
        wb = waitbar(0,sprintf('exporting %.1f-%.1f s ...',t0,t0+dur));
        v.CurrentTime = (f0-1)/fpsW; kk = f0; i = 0;
        while kk <= f1 && hasFrame(v)
            img = readFrame(v);
            if mod(i,dstep)==0
                if size(img,3)==1, img = repmat(img,1,1,3); end
                vb = imresize(img,[vH Wout]);                            % video aligned with trace, aspect preserved
                c0 = max(Lpx+1, colOf(t(kk)-EXPORT_WIN/2)); c1 = min(Ws, c0+wcols);
                ws = [uint8(255*ones(Hs,gapL,3)), strip(:,1:Lpx,:), strip(:,c0:c1,:), uint8(255*ones(Hs,gapR,3))];
                xc = gapL + Lpx + round((c1-c0)/2); cc = max(1,xc-1):min(size(ws,2),xc+1);
                ws(rTop:rBot,cc,1) = 230; ws(rTop:rBot,cc,2) = 30; ws(rTop:rBot,cc,3) = 30;   % red playhead, axes-box height
                tp = imresize(ws,[Htr Wout]);                           % trace panel
                if ~isempty(cap)                                        % overlay 'Time (s)' text inside the bottom margin
                    reg = tp(capR0:capR0+Hcap-1,:,:); reg(capMask) = cap(capMask); tp(capR0:capR0+Hcap-1,:,:) = reg;
                end
                writeVideo(vw, [vb; tp]);                               % video | rolling trace (caption inside)
                if isvalid(wb), waitbar((kk-f0)/(f1-f0), wb); end
            end
            kk = kk + 1; i = i + 1;
        end
        close(vw); if isvalid(wb), delete(wb); end
        if isgraphics(txt), set(txt,'String', sprintf('exported: %s', outfile)); end
        fprintf('exported %s\n', outfile);
    end

    function onClose(~,~), stopplay(); delete(f); end
end

% ================= helpers =================
function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*labeled*.mp4'));
    S = {}; key = {};
    for i = 1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1} = k; S{end+1} = {tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S); rn = cellfun(@(c) c{2}, S);
        [~,o] = sortrows([an(:) rn(:)]); S = S(o);
    end
end
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function mp4 = pick_labeled_mp4(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*labeled*.mp4']));
    assert(~isempty(d), ['no *_labeled.mp4 matching %s* in %s\n' ...
        '(raw .avi is FFV1 and unreadable by MATLAB; create a labeled/transcoded mp4)'], prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    mp4 = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
