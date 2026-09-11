function whisk_breath_video_gui_RZ()
% whisk_breath_video_gui_RZ  Scrub the whisker video alongside L/R whisk-angle
% traces AND the thermal nostril breathing trace, all frame-locked, and export a
% clip in the same layout.
%   Top    : video frame (DLC *_labeled.mp4; raw .avi is FFV1 = unreadable here)
%   Middle : L (green) / R (blue) BP whisk angle with a playhead at the frame
%   Bottom : breathing (inhale UP), detrended + band-passed thermal nostril trace
% Navigate sessions with the dropdown / Prev / Next. Type a start and duration
% and press "Export clip" to render that chunk (video + both traces) to an mp4.
% Frame k in the video == row k in the DLC csv == xL(k)/xR(k).
%
% Breathing is loaded from <stem>_breath.mat (resolved from the matching NOSE
% DLC csv via thermal_resolve_paths). The thermal and whisk cameras share the
% WFG trigger and both start at t=0, so breath is aligned by REAL SECONDS (no
% resampling) -- it is simply x-linked to the whisk time axis.
%
% Sign: drawROI_N_lpsub saves _breath.mat already inhale-up (INVERT=true,
% recorded in B.inverted), so the panel is oriented inhale-up directly. Set
% BR_SIGN=-1 to flip if a session reads the other way.
% This is a NEW file -- whisk_video_gui_RZ.m is left untouched.

% ============================ USER-EDITABLE ============================
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
dataRoot = "D:\260615_thermalNbasler";
BP       = [5 60];      % fast-whisking bandpass (Hz)
BR_BP    = [1 20];      % breathing bandpass (Hz); the 1 Hz HP detrends
BR_SIGN  = +1;          % breath display sign: +1 = inhale up (canonical); -1 flips
fpsW     = 400;         % whisk camera fps
WINSEC   = 1;         % trace window shown around the playhead (s)
T0       = 39;          % initial playhead time (s)
EXPORT_FPS = 100;       % exported-clip playback fps (the "speed" box sets slow-mo: 0.25 = 4x slower)
EXPORT_DIR = "";        % output folder ("" = this script's folder)
EXPORT_WIN = 1.0;       % exported trace: rolling window width (s), like the GUI
TRACE_TALL = 1.4;       % exported WHISK trace panel height factor (vertical; 1 = default)
BR_TALL    = 1.0;       % exported BREATH trace panel height factor (vertical; 1 = default)
% --- export layout knobs (tweak to reposition; all in the exported frame) ---
EXP_BASE_H  = 240;      % base trace-panel height in px (final height = EXP_BASE_H * *_TALL)
EXP_PPS     = 600;      % rolling strip resolution (px per second of trace)
EXP_YAXIS_FRAC = 0.10;  % y-axis (label+ticks) width as a fraction of the VISIBLE trace panel; raise for more room
EXP_RMARGIN = 0.01;     % RIGHT margin (fraction)
EXP_AXBOT   = 0.18;     % axes BOTTOM (fraction) = room for x tick numbers; raise for more room
EXP_AXTOP   = 0.05;     % axes TOP margin (fraction)
EXP_CAP_H   = 30;       % 'Time (s)' caption band height px (0 = no caption); font scales with it
EXP_SIDEGAP = 0.03;     % white gap on LEFT & RIGHT of the trace panel (fraction of panel width)
EXP_CAP_Y   = 0.93;     % 'Time (s)' vertical center as fraction of trace height (SMALLER = higher / toward top)
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % thermal_resolve_paths
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');
cL = [0 0.55 0]; cR = [0 0.4 0.85]; cBr = [0.5 0.2 0.7];

% ---- enumerate sessions (one per readable *_labeled.mp4) ----
sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_labeled.mp4 in %s', whiskDir);
labels = cellfun(@(c) sprintf('%s  n%d', c{1}, c{2}), sess, 'UniformOutput',false);

% ---- state (filled by loadSession) ----
v=[]; xL=[]; xR=[]; t=[]; n=0; cf=max(1,round(T0*fpsW)+1); yl=[0 1];
animal=""; kRun=0; tmr=[]; si=1; loading=false;
him=[]; pL=[]; pR=[]; hPlay=[];
brf=[]; tBr=[]; fpsB=0; hasBreath=false; pBr=[]; hPlayB=[]; ylB=[-1 1];
PERIOD = round(1000/60)/1000;   % timer period (s), ms-rounded to avoid precision warning

% ---- UI ----
f = figure('Color','w','Position',[80 50 1180 860],'Name','whisk video + angle + breath', ...
    'NumberTitle','off','CloseRequestFcn',@onClose);
axV = axes(f,'Units','normalized','Position',[0.05 0.55 0.90 0.42]);
axT = axes(f,'Units','normalized','Position',[0.07 0.355 0.88 0.155]);
axB = axes(f,'Units','normalized','Position',[0.07 0.175 0.88 0.155]);

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
        loading = true; cleanup = onCleanup(@() setLoading(false)); %#ok<NASGU>
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

        % ---- breath: thermal nostril trace (aligned by shared trigger, t=0) ----
        loadBreath();

        % video
        cla(axV); him = imshow(read(v,cf),'Parent',axV);
        % whisk traces (middle) -- x ticks hidden; breath panel carries the time axis
        cla(axT); hold(axT,'on'); grid(axT,'on');
        pL = plot(axT, t(1:n), xL(1:n), '-', 'Color',cL, 'LineWidth',0.9);
        pR = plot(axT, t(1:n), xR(1:n), '-', 'Color',cR, 'LineWidth',0.9);
        ylabel(axT,'Whisker angle (degree)');
        legend([pL pR],{'L whisk','R whisk'},'Orientation','horizontal','Location','northeast');
        yl = [min([xL(1:n);xR(1:n)]) max([xL(1:n);xR(1:n)])]; yl = yl + 0.05*diff(yl)*[-1 1];
        ylim(axT, yl);
        hPlay = line(axT,[t(cf) t(cf)], yl, 'Color',[0.85 0.1 0.1],'LineWidth',1.2,'HandleVisibility','off');

        % breath trace (bottom)
        cla(axB); hold(axB,'on'); grid(axB,'on');
        if hasBreath
            pBr = plot(axB, tBr, brf, '-', 'Color',cBr, 'LineWidth',0.9);
            ylB = [min(brf) max(brf)]; if ~(diff(ylB)>0), ylB=[-1 1]; end
            ylB = ylB + 0.08*diff(ylB)*[-1 1]; ylim(axB, ylB);
            hPlayB = line(axB,[t(cf) t(cf)], ylB, 'Color',[0.85 0.1 0.1],'LineWidth',1.2,'HandleVisibility','off');
        else
            pBr = []; hPlayB = [];
            text(axB,0.5,0.5,sprintf('no breath signal for %s n%d',char(animal),kRun), ...
                'Units','normalized','HorizontalAlignment','center','Color',[0.6 0.1 0.1]);
        end
        ylabel(axB,'Breath (inhale \uparrow)'); xlabel(axB,'time (s)');
        set(axT,'XTickLabel',[]);                 % middle panel: no x numbers
        linkaxes([axT axB],'x');                  % breath follows whisk x-window

        set(sld,'Min',1,'Max',n,'Value',cf,'SliderStep',[1/(n-1) 30/(n-1)]);
        set(pop,'Value',si);
        redraw();
    end
    function setLoading(tf), loading = tf; end

    function loadBreath()
        hasBreath = false; brf = []; tBr = []; fpsB = 0; ylB = [-1 1];
        try
            noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', char(animal), kRun));
            Pn = thermal_resolve_paths(noseCsv, dataRoot);
            if ~isfile(Pn.breath)
                fprintf('breath: no _breath.mat for %s n%d (run drawROI_N_lpsub)\n', char(animal), kRun); return;
            end
            Bs = load(Pn.breath); br = Bs.breath(:); fpsB = double(Bs.fps);
            if isfield(Bs,'inverted') && ~Bs.inverted, br = -br; end   % normalize to inhale-up
            bpHi = min(BR_BP(2), 0.45*fpsB);                            % keep below Nyquist
            [bb,ab] = butter(2, [BR_BP(1) bpHi]/(fpsB/2), 'bandpass');
            brf = BR_SIGN * filtfilt(bb,ab, fillmissing(br,'linear'));  % inhale-up, detrended
            tBr = (0:numel(brf)-1)'/fpsB;
            hasBreath = true;
        catch ME
            fprintf('breath load skipped (%s n%d): %s\n', char(animal), kRun, ME.message);
        end
    end

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
        tc = t(cf); set(hPlay,'XData',[tc tc]);
        if ~isempty(hPlayB) && isgraphics(hPlayB), set(hPlayB,'XData',[tc tc]); end
        xlim(axT, tc + WINSEC*[-0.5 0.5]);   % axB follows via linkaxes
        set(sld,'Value',cf);
        bstr = '';
        if hasBreath
            jb = max(1, min(numel(brf), round(tc*fpsB)+1)); bstr = sprintf('   br=%+6.2f', brf(jb));
        end
        set(txt,'String', sprintf('frame %d / %d   t=%.3f s   xL=%+6.2f   xR=%+6.2f%s', cf, n, tc, xL(cf), xR(cf), bstr));
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
        outfile = fullfile(od, sprintf('%s_n%d_t%05.1f-%05.1fs_%gx_breathwhisk.mp4', char(animal), kRun, t0, t0+dur, s));
        vw = VideoWriter(outfile,'MPEG-4'); vw.FrameRate = outfps; vw.Quality = 95; open(vw);

        % ROLLING traces like the GUI: pre-render the WHOLE clip as wide strips
        % (with the y-axis on the left), then per frame crop the EXPORT_WIN window
        % and paste it. WHISK (middle) + BREATH (bottom) share horizontal geometry
        % so their rolling windows stay column-aligned. No per-frame getframe.
        HtrW  = round(EXP_BASE_H*TRACE_TALL); HtrW  = HtrW  - mod(HtrW,2);   % whisk panel height
        HtrBr = round(EXP_BASE_H*BR_TALL);    HtrBr = HtrBr - mod(HtrBr,2);  % breath panel height
        ws_t0 = t0 - EXPORT_WIN/2; ws_t1 = t0 + dur + EXPORT_WIN/2;          % render with half-window margins
        i0 = max(1, round(ws_t0*fpsW)+1); i1 = min(n, round(ws_t1*fpsW)+1);
        Rmn = EXP_RMARGIN;                                                  % right margin (fraction of strip)
        winFrac = EXPORT_WIN/(ws_t1-ws_t0);                                 % visible window / full rendered span
        % choose strip left margin so the y-axis is EXP_YAXIS_FRAC of the VISIBLE panel after cropping
        Lmn = EXP_YAXIS_FRAC*winFrac*(1-Rmn) / (1 - EXP_YAXIS_FRAC + EXP_YAXIS_FRAC*winFrac);
        Lmn = min(max(Lmn,0.005),0.30);
        Wfig = max(400, round((ws_t1-ws_t0)*EXP_PPS));

        % --- whisk strip (x numbers hidden when breath present: breath carries them) ---
        W = render_trace_strip(t(i0:i1), {xL(i0:i1), xR(i0:i1)}, {cL,cR}, ...
            'Whisker angle (degree)', yl, ws_t0, ws_t1, Wfig, EXP_BASE_H, Lmn, Rmn, EXP_AXBOT, EXP_AXTOP, ~hasBreath);
        Hs = W.Hs; Ws = W.Ws; Lpx = W.Lpx; plotW = W.plotW;
        colOf = @(tt) Lpx + round((tt-ws_t0)/(ws_t1-ws_t0)*plotW);          % time -> column in strip
        wcols = colOf(t0+EXPORT_WIN/2) - colOf(t0-EXPORT_WIN/2);            % fixed window width (px)
        gapL = max(2, round(EXP_SIDEGAP*(Lpx+wcols))); gapR = round(1.5*gapL);   % right gap 1.5x left
        Wc = gapL + gapR + Lpx + wcols + 1;                                % gapL | y-axis | window | gapR
        Wout = round(HtrW * Wc / Hs); Wout = Wout - mod(Wout,2);           % width from whisk crop aspect (master)
        rTop = W.rTop; rBot = W.rBot;                                      % playhead spans the axes box only
        vH = round(v.Height*Wout/v.Width); vH = vH - mod(vH,2);            % video same width (aligned), aspect kept

        % --- breath strip (same horizontal geometry; bottom panel shows x numbers) ---
        Br = struct(); rTopBr = 0; rBotBr = 0;
        if hasBreath
            jb0 = max(1, floor(ws_t0*fpsB)+1); jb1 = min(numel(brf), ceil(ws_t1*fpsB)+1);
            Br = render_trace_strip(tBr(jb0:jb1), {brf(jb0:jb1)}, {cBr}, ...
                'Breath (inhale \uparrow)', ylB, ws_t0, ws_t1, Wfig, EXP_BASE_H, Lmn, Rmn, EXP_AXBOT, EXP_AXTOP, true);
            if Br.Ws ~= Ws, Br.strip = imresize(Br.strip,[size(Br.strip,1) Ws]); Br.Hs = size(Br.strip,1); end
            rTopBr = max(1,round(EXP_AXTOP*Br.Hs)); rBotBr = min(Br.Hs, round((1-EXP_AXBOT)*Br.Hs));
        end

        % --- 'Time (s)' caption on the BOTTOM panel (breath if present, else whisk) ---
        Hcap = EXP_CAP_H - mod(EXP_CAP_H,2);
        capPanelH = HtrW; if hasBreath, capPanelH = HtrBr; end
        cap = []; capMask = []; capR0 = [];
        if Hcap > 0                                                        % render 'Time (s)' in MATLAB's font (matches axes)
            capF = figure('Color','w','Position',[60 60 Wout Hcap],'Name','caption','NumberTitle','off');
            ca = axes(capF,'Position',[0 0 1 1]); axis(ca,'off');
            text(ca,0.5,0.5,'Time (s)','Units','normalized','HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName',W.axFont,'FontSize',max(8,round(W.axFsz*capPanelH/EXP_BASE_H)));
            drawnow; cg = getframe(capF); cap = imresize(cg.cdata,[Hcap Wout]); delete(capF);
            capMask = repmat(any(cap<160,3),1,1,3);                        % text pixels only (don't paint over ticks)
            capR0 = min(capPanelH-Hcap+1, max(1, round(EXP_CAP_Y*capPanelH - Hcap/2)));   % top row inside the bottom panel
        end

        wb = waitbar(0,sprintf('exporting %.1f-%.1f s ...',t0,t0+dur));
        v.CurrentTime = (f0-1)/fpsW; kk = f0; i = 0;
        while kk <= f1 && hasFrame(v)
            img = readFrame(v);
            if mod(i,dstep)==0
                if size(img,3)==1, img = repmat(img,1,1,3); end
                vb = imresize(img,[vH Wout]);                              % video aligned with traces, aspect preserved
                c0 = max(Lpx+1, colOf(t(kk)-EXPORT_WIN/2)); c1 = min(Ws, c0+wcols);
                xc = gapL + Lpx + round((c1-c0)/2); cc = max(1,xc-1):xc+1;  % red playhead columns (in the cropped panel)

                % whisk panel
                wsW = [uint8(255*ones(Hs,gapL,3)), W.strip(:,1:Lpx,:), W.strip(:,c0:c1,:), uint8(255*ones(Hs,gapR,3))];
                ccW = cc(cc<=size(wsW,2));
                wsW(rTop:rBot,ccW,1) = 230; wsW(rTop:rBot,ccW,2) = 30; wsW(rTop:rBot,ccW,3) = 30;
                tpW = imresize(wsW,[HtrW Wout]);
                if ~hasBreath && ~isempty(cap)                             % caption rides the whisk panel only if no breath
                    reg = tpW(capR0:capR0+Hcap-1,:,:); reg(capMask) = cap(capMask); tpW(capR0:capR0+Hcap-1,:,:) = reg;
                end

                if hasBreath
                    % breath panel (same crop columns -> aligned with whisk)
                    wsB = [uint8(255*ones(Br.Hs,gapL,3)), Br.strip(:,1:Lpx,:), Br.strip(:,c0:c1,:), uint8(255*ones(Br.Hs,gapR,3))];
                    ccB = cc(cc<=size(wsB,2));
                    wsB(rTopBr:rBotBr,ccB,1) = 230; wsB(rTopBr:rBotBr,ccB,2) = 30; wsB(rTopBr:rBotBr,ccB,3) = 30;
                    tpBr = imresize(wsB,[HtrBr Wout]);
                    if ~isempty(cap)
                        reg = tpBr(capR0:capR0+Hcap-1,:,:); reg(capMask) = cap(capMask); tpBr(capR0:capR0+Hcap-1,:,:) = reg;
                    end
                    writeVideo(vw, [vb; tpW; tpBr]);                       % video | whisk | breath
                else
                    writeVideo(vw, [vb; tpW]);                             % video | whisk
                end
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
function R = render_trace_strip(tvec, sigs, cols, ylab, yl, ws_t0, ws_t1, Wfig, baseH, Lmn, Rmn, axbot, axtop, showX)
% Pre-render one trace panel as a wide strip; return the image + geometry so the
% per-frame rolling crop can map time -> column. WHISK and BREATH call this with
% identical horizontal args (Wfig, Lmn, Rmn) so their crops stay column-aligned.
    tf = figure('Color','w','Position',[60 60 Wfig baseH],'Name','rendering trace...','NumberTitle','off');
    tax = axes(tf,'Units','normalized','Position',[Lmn axbot 1-Lmn-Rmn 1-axbot-axtop]); hold(tax,'on'); grid(tax,'on');
    for s = 1:numel(sigs), plot(tax, tvec, sigs{s}, '-', 'Color',cols{s}, 'LineWidth',1); end
    xlim(tax,[ws_t0 ws_t1]); if ~isempty(yl) && diff(yl)>0, ylim(tax, yl); end
    xlabel(tax,''); ylabel(tax, ylab);
    if ~showX, set(tax,'XTickLabel',[]); end
    drawnow; frb = getframe(tf);
    R.strip = frb.cdata; R.axFont = get(tax,'FontName'); R.axFsz = get(tax,'FontSize'); delete(tf);
    R.Hs = size(R.strip,1); R.Ws = size(R.strip,2);
    R.Lpx = round(Lmn*R.Ws); R.plotW = R.Ws - R.Lpx - round(Rmn*R.Ws);
    R.rTop = max(1,round(axtop*R.Hs)); R.rBot = min(R.Hs, round((1-axbot)*R.Hs));
end
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
