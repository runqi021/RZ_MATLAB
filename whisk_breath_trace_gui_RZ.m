function whisk_breath_trace_gui_RZ()
% whisk_breath_trace_gui_RZ  Browse sessions and export a breathing + whisking
% trace figure over a chosen window.
%   Top    : L (green) / R (blue) BP whisk angle (deg)
%   Bottom : BP breathing (nose thermal), inspiration up
% Pick a session (dropdown / Prev / Next), type START (s) and DURATION (s),
% preview updates live; "Export" saves the window as a vector PDF (+ PNG).

% ============================ USER-EDITABLE ============================
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
dataRoot = "D:\260615_thermalNbasler";
BP       = [6 30];      % whisk bandpass (Hz)
BR_BP    = [2 15];      % breath bandpass (Hz)
fpsW     = 400;         % whisk camera fps
T0       = 39;          % initial start time (s)
DUR0     = 3;           % initial duration (s)
EXPORT_DIR = "";        % output folder ("" = this script's folder)
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');
cL=[0 0.55 0]; cR=[0 0.4 0.85]; cB=[0.1 0.1 0.1];

sess = list_sessions(whiskDir);
assert(~isempty(sess),'no *_whisk_n*.csv in %s', whiskDir);
labels = cellfun(@(c) sprintf('%s  n%d', c{1}, c{2}), sess, 'UniformOutput',false);

% state
xL=[]; xR=[]; t=[]; brf=[]; tB=[]; animal=""; kRun=0; si=1; loading=false;

f = figure('Color','w','Position',[80 90 1180 720],'Name','breath + whisk trace export', ...
    'NumberTitle','off','CloseRequestFcn',@(s,e)delete(s));
axT = axes(f,'Units','normalized','Position',[0.08 0.26 0.86 0.58]);

% nav row
uicontrol(f,'Style','text','Units','normalized','Position',[0.08 0.11 0.05 0.03],'String','session', ...
    'BackgroundColor','w','HorizontalAlignment','left');
pop = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.13 0.105 0.18 0.04], ...
    'String',labels,'Value',1,'Callback',@onPick);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.32 0.105 0.05 0.04],'String','< Prev','Callback',@(s,e)onStep(-1));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.37 0.105 0.05 0.04],'String','Next >','Callback',@(s,e)onStep(+1));

% window + export row
uicontrol(f,'Style','text','Units','normalized','Position',[0.08 0.045 0.06 0.03],'String','start (s)', ...
    'BackgroundColor','w','HorizontalAlignment','left');
eStart = uicontrol(f,'Style','edit','Units','normalized','Position',[0.14 0.04 0.06 0.045],'String',num2str(T0),'Callback',@(s,e)applyWin());
uicontrol(f,'Style','text','Units','normalized','Position',[0.21 0.045 0.06 0.03],'String','dur (s)', ...
    'BackgroundColor','w','HorizontalAlignment','left');
eDur = uicontrol(f,'Style','edit','Units','normalized','Position',[0.27 0.04 0.06 0.045],'String',num2str(DUR0),'Callback',@(s,e)applyWin());
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.35 0.038 0.11 0.05],'String','Export PDF', ...
    'FontWeight','bold','Callback',@onExport);
txt = uicontrol(f,'Style','text','Units','normalized','Position',[0.48 0.045 0.48 0.035],'String','', ...
    'BackgroundColor','w','HorizontalAlignment','left','FontName','Consolas');

loadSession(1);

% ================= session loading =================
    function loadSession(idx)
        if loading, return; end
        loading = true; cu = onCleanup(@() setLoading(false));
        si = idx; animal = sess{idx}{1}; kRun = sess{idx}{2};
        M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',char(animal),kRun)), 0.6);
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        t  = (0:numel(xL)-1)'/fpsW;
        brf = []; tB = [];
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',char(animal),kRun)), dataRoot);
            if isfile(Pn.breath)
                Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps);
                [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
                tB=(0:numel(brf)-1)'/fb;
            end
        catch
        end
        drawTraces(axT);
        set(pop,'Value',si); applyWin();
    end
    function setLoading(tf), loading = tf; end

    function drawTraces(ax)
        cla(ax,'reset');
        yyaxis(ax,'left'); hold(ax,'on');
        pL = plot(ax, t, xL, '-', 'Color',cL, 'LineWidth',0.9);
        pR = plot(ax, t, xR, '-', 'Color',cR, 'LineWidth',0.9);
        ylabel(ax,'Whisker angle (degree)'); set(ax,'YColor','k');
        yyaxis(ax,'right'); hold(ax,'on'); pB = [];
        if ~isempty(brf), pB = plot(ax, tB, brf, '-', 'Color',cB, 'LineWidth',1.1); end
        ylabel(ax,'Breathing (BP)'); set(ax,'YColor',cB);
        grid(ax,'on'); box(ax,'on'); xlabel(ax,'time (s)');
        hobj=[pL pR]; lab={'L whisk','R whisk'};
        if ~isempty(pB), hobj=[hobj pB]; lab=[lab {'breath'}]; end
        legend(ax, hobj, lab, 'Orientation','horizontal','Location','northoutside');
    end

% ================= navigation / window =================
    function onPick(src,~), loadSession(get(src,'Value')); end
    function onStep(d), idx=min(numel(sess),max(1,si+d)); if idx~=si, loadSession(idx); end, end
    function [t0,t1] = winRange()
        t0 = str2double(get(eStart,'String')); dur = str2double(get(eDur,'String'));
        if isnan(t0), t0=0; end, if isnan(dur)||dur<=0, dur=3; end
        t1 = t0 + dur;
    end
    function applyWin()
        [t0,t1] = winRange();
        xlim(axT,[t0 t1]);
        set(txt,'String', sprintf('%s n%d   window %.2f - %.2f s', char(animal), kRun, t0, t1));
    end

% ================= export =================
    function onExport(~,~)
        [t0,t1] = winRange();
        if strlength(EXPORT_DIR)>0, od=char(EXPORT_DIR); else, od=fileparts(mfilename('fullpath')); end
        stem = fullfile(od, sprintf('%s_n%d_t%05.1f-%05.1fs_breathwhisk', char(animal), kRun, t0, t1));
        ef = figure('Color','w','Position',[100 100 1100 460]);
        eax = axes(ef);
        drawTraces(eax);
        xlim(eax,[t0 t1]);
        title(eax, sprintf('%s n%d   %.2f - %.2f s', char(animal), kRun, t0, t1), 'Interpreter','none');
        drawnow;
        try
            exportgraphics(ef, [stem '.pdf'], 'ContentType','vector', 'BackgroundColor','none');
            exportgraphics(ef, [stem '.png'], 'Resolution',300);
            set(txt,'String', sprintf('exported: %s.pdf', stem));
            fprintf('exported %s.pdf (+ .png)\n', stem);
        catch ME
            set(txt,'String', ['export error: ' ME.message]);
        end
        if isgraphics(ef), delete(ef); end
    end
end

% ================= helpers =================
function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv'));
    S = {}; key = {};
    for i = 1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1}=k; S{end+1}={tok{1}, str2double(tok{2})}; %#ok<AGROW>
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
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
