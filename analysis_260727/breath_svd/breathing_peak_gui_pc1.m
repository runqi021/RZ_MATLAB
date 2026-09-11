function breathing_peak_gui_pc1()
%BREATHING_PEAK_GUI_PC1  Interactive inspiration-peak detector on SVD PC1.
%
%   Identical workflow to breathing_peak_gui.m, but the breathing signal is
%   PC1 from breath_pc1.mat (written by breath_svd_pc1.m) instead of being
%   rebuilt from DLC dots. PC1 is detrended + z-scored so the same peak
%   thresholds apply; SVD sign is arbitrary so an Invert toggle is provided.
%
%   Sessions = folders containing a breath_pc1.mat file.
%
%   Output (per session):
%       breath_peak_pc1.mat
%           insp_onset_idx, insp_onsets_t, insp_onsets_amp, insp_onsets_train,
%           breath, t_breath, findpeak_params   (same fields as the DLC version,
%           so breathing_trough_gui_pc1.m consumes it unchanged)
%       breath_peak_pc1_figure.fig / .png
%
%   Requires: MATLAB R2021a+

%% ── Shared state ─────────────────────────────────────────────────────────
S.pc1raw     = [];      % [F×1] raw PC1 from the mat
S.eigImg     = [];      % PC1 spatial eigenimage (for confirming the ROI)
S.diffImg    = [];      % inspiration-baseline diff map (top100 - bottom100 frames)
S.t          = [];      % [F×1] seconds
S.breath     = [];      % [F×1] processed (detrended, z-scored, optionally inverted)
S.fps        = NaN;
S.peaks      = [];      % [K×1] frame indices of accepted peaks
S.folderPath = '';
S.label      = '';
S.editMode   = 'none';
S.masterFolder = '';
S.sessionList  = {};
S.sessionIdx   = 0;

DEFAULTS.minProm   = 1.00;
DEFAULTS.minDist   = 0.20;
DEFAULTS.minHeight = 0.50;
DEFAULTS.minWidth  = 0.05;
DEFAULTS.snapWin   = 0.10;
DEFAULTS.invert    = false;

%% ── Figure ──────────────────────────────────────────────────────────────
fig = uifigure('Name','Breathing Peak Sorter (SVD PC1)', ...
               'Position',[40 30 1560 930], 'WindowState','maximized');
rootGL = uigridlayout(fig, [1 2]);
rootGL.ColumnWidth = {320, '1x'};
rootGL.Padding = [6 6 6 6]; rootGL.ColumnSpacing = 8;

ctrlPan = uipanel(rootGL, 'Title','Controls','FontSize',12,'FontWeight','bold');
ctrlPan.Layout.Column = 1;
nCtrlRows = 28;
cGL = uigridlayout(ctrlPan, [nCtrlRows 2]);
cGL.RowHeight = repmat({'fit'}, 1, nCtrlRows);
cGL.ColumnWidth = {'fit','1x'}; cGL.Padding = [8 8 8 8]; cGL.RowSpacing = 4;
r = 0;

% ── Session Navigation ──
r=r+1; h=uilabel(cGL,'Text','── Session Navigation ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnLoad=uibutton(cGL,'Text','Browse Master Folder','ButtonPushedFcn',@cb_loadMaster); btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];
r=r+1; lblMasterPath=uilabel(cGL,'Text','No folder selected','WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]); lblMasterPath.Layout.Row=r; lblMasterPath.Layout.Column=[1 2];
r=r+1; lblCurrentFolder=uilabel(cGL,'Text','—','FontWeight','bold','WordWrap','on','FontSize',10); lblCurrentFolder.Layout.Row=r; lblCurrentFolder.Layout.Column=[1 2];
r=r+1; lblNavStatus=uilabel(cGL,'Text','0 / 0  folders  |  0 complete','FontSize',10,'FontColor',[0.25 0.65 0.25]); lblNavStatus.Layout.Row=r; lblNavStatus.Layout.Column=[1 2];
r=r+1; btnPrev=uibutton(cGL,'Text','← Prev','ButtonPushedFcn',@cb_prevSession); btnPrev.Layout.Row=r; btnPrev.Layout.Column=1;
       btnNext=uibutton(cGL,'Text','Next →','ButtonPushedFcn',@cb_nextSession); btnNext.Layout.Row=r; btnNext.Layout.Column=2;

% ── Signal (PC1) ──
r=r+1; h=uilabel(cGL,'Text','── Signal (PC1) ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; h=uilabel(cGL,'Text','fps (from calcium):'); h.Layout.Row=r; h.Layout.Column=1;
       efFps=uieditfield(cGL,'numeric','Value',30,'Limits',[1 2000],'ValueChangedFcn',@cb_updateTrace); efFps.Layout.Row=r; efFps.Layout.Column=2;
r=r+1; cbInvert=uicheckbox(cGL,'Text','Invert signal (SVD sign is arbitrary)','Value',false,'ValueChangedFcn',@cb_updateTrace); cbInvert.Layout.Row=r; cbInvert.Layout.Column=[1 2];
r=r+1; lblInfo=uilabel(cGL,'Text','','FontSize',9,'FontColor',[0.3 0.3 0.6],'WordWrap','on'); lblInfo.Layout.Row=r; lblInfo.Layout.Column=[1 2];
r=r+1; btnTrace=uibutton(cGL,'Text','Update Trace','ButtonPushedFcn',@cb_updateTrace); btnTrace.Layout.Row=r; btnTrace.Layout.Column=[1 2];

% ── Peak Detection ──
r=r+1; h=uilabel(cGL,'Text','── Peak Detection ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
paramLabels = {'Min Prominence:','Min Distance (s):','Min Height:','Min Width (s):'};
paramDefs   = {DEFAULTS.minProm, DEFAULTS.minDist, DEFAULTS.minHeight, DEFAULTS.minWidth};
efParam = gobjects(4,1);
for p = 1:4
    r=r+1; h=uilabel(cGL,'Text',paramLabels{p}); h.Layout.Row=r; h.Layout.Column=1;
    efParam(p)=uieditfield(cGL,'numeric','Value',paramDefs{p},'Limits',[-Inf Inf],'ValueChangedFcn',@(~,~) saveWorkingParams());
    efParam(p).Layout.Row=r; efParam(p).Layout.Column=2;
end
r=r+1; h=uilabel(cGL,'Text','Snap window (s):'); h.Layout.Row=r; h.Layout.Column=1;
       efSnap=uieditfield(cGL,'numeric','Value',DEFAULTS.snapWin,'Limits',[0 10],'ValueChangedFcn',@(~,~) saveWorkingParams()); efSnap.Layout.Row=r; efSnap.Layout.Column=2;
r=r+1; btnDetect=uibutton(cGL,'Text','Detect Peaks','ButtonPushedFcn',@cb_detect,'BackgroundColor',[0.18 0.65 0.30],'FontColor','white'); btnDetect.Layout.Row=r; btnDetect.Layout.Column=[1 2];

% ── Manual Editing ──
r=r+1; h=uilabel(cGL,'Text','── Manual Editing ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; lblMode=uilabel(cGL,'Text','Mode: Navigate','FontWeight','bold','FontSize',11,'FontColor',[0.90 0.80 0.10]); lblMode.Layout.Row=r; lblMode.Layout.Column=[1 2];
r=r+1; btnAdd=uibutton(cGL,'Text','+ Add Mode','ButtonPushedFcn',@cb_setAdd,'BackgroundColor',[0.20 0.42 0.85],'FontColor','white'); btnAdd.Layout.Row=r; btnAdd.Layout.Column=1;
       btnDel=uibutton(cGL,'Text','- Delete Mode','ButtonPushedFcn',@cb_setDelete,'BackgroundColor',[0.82 0.22 0.18],'FontColor','white'); btnDel.Layout.Row=r; btnDel.Layout.Column=2;
r=r+1; btnNav=uibutton(cGL,'Text','Navigate (N)','ButtonPushedFcn',@cb_setNav); btnNav.Layout.Row=r; btnNav.Layout.Column=1;
       btnClear=uibutton(cGL,'Text','Clear All Peaks','ButtonPushedFcn',@cb_clearPeaks,'BackgroundColor',[0.40 0.10 0.10],'FontColor','white'); btnClear.Layout.Row=r; btnClear.Layout.Column=2;
r=r+1; lblCount=uilabel(cGL,'Text','Peaks: 0','FontWeight','bold','FontSize',12); lblCount.Layout.Row=r; lblCount.Layout.Column=[1 2];

% ── Save ──
r=r+1; h=uilabel(cGL,'Text','── Save ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnSave=uibutton(cGL,'Text','Save breath_peak_pc1.mat','ButtonPushedFcn',@cb_save,'BackgroundColor',[0.65 0.45 0.05],'FontColor','white'); btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];
r=r+1; btnReset=uibutton(cGL,'Text','Reset to Defaults','ButtonPushedFcn',@cb_reset,'BackgroundColor',[0.30 0.30 0.30],'FontColor','white'); btnReset.Layout.Row=r; btnReset.Layout.Column=[1 2];

if ispref('BreathingPeakGUI_PC1','workingParams')
    try, applyWorkingParams(getpref('BreathingPeakGUI_PC1','workingParams')); catch, end
end

%% ── Right: plots ────────────────────────────────────────────────────────
rightGL = uigridlayout(rootGL, [2 1]); rightGL.Layout.Column = 2;
rightGL.RowHeight = {'1x', 300}; rightGL.Padding = [0 0 0 0]; rightGL.RowSpacing = 6;

axMain = uiaxes(rightGL); axMain.Layout.Row = 1;
axMain.Title.String = 'PC1 breathing trace  |  Left-click: add/delete peaks  |  N: navigate';
axMain.Title.Interpreter='none'; axMain.XLabel.String='Time (s)'; axMain.XLabel.Interpreter='none';
axMain.YLabel.String='PC1 (z-score)'; axMain.YLabel.Interpreter='none';
hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');

imgGL = uigridlayout(rightGL, [1 2]); imgGL.Layout.Row = 2;
imgGL.Padding = [0 0 0 0]; imgGL.ColumnSpacing = 8;
axEig = uiaxes(imgGL); axEig.Layout.Column = 1;
axEig.Title.String = 'PC1 eigenimage'; axEig.Title.Interpreter='none';
axis(axEig,'image'); axEig.XTick=[]; axEig.YTick=[];
axDiff = uiaxes(imgGL); axDiff.Layout.Column = 2;
axDiff.Title.String = 'top100 - bottom100 (inspiration - baseline)'; axDiff.Title.Interpreter='none';
axis(axDiff,'image'); axDiff.XTick=[]; axDiff.YTick=[];

hPeakSc = gobjects(0);
fig.KeyPressFcn = @cb_key;

%% ═══════════════════════════════ CALLBACKS ═══════════════════════════════
    function cb_loadMaster(~,~)
        master = uigetdir(S.masterFolder,'Select Master Folder');
        if isequal(master,0), return; end
        S.masterFolder = master;
        sessions = scanForSessions(master);
        if isempty(sessions)
            uialert(fig,'No breath_pc1.mat found in any subfolder.','No Sessions'); return;
        end
        S.sessionList = sessions; S.sessionIdx = 1;
        [~, mName] = fileparts(master); lblMasterPath.Text = mName;
        loadCurrentSession();
    end
    function cb_prevSession(~,~)
        if S.sessionIdx <= 1, return; end
        S.sessionIdx = S.sessionIdx - 1; loadCurrentSession();
    end
    function cb_nextSession(~,~)
        if S.sessionIdx >= numel(S.sessionList), return; end
        S.sessionIdx = S.sessionIdx + 1; loadCurrentSession();
    end

    function loadCurrentSession()
        sess = S.sessionList{S.sessionIdx};
        S.folderPath = sess.folder; S.label = sess.label; S.peaks = [];
        lblCurrentFolder.Text = sess.label; updateNavDisplay();

        try
            PC = load(sess.pc1Mat);
        catch ME
            uialert(fig, ME.message,'Read Error'); return;
        end
        if ~isfield(PC,'breathTrace')
            uialert(fig,'breath_pc1.mat has no breathTrace (run breath_svd_pc1.m).','Bad file'); return;
        end
        S.pc1raw = double(PC.breathTrace(:));
        if isfield(PC,'t'), S.t = double(PC.t(:)); else, S.t = (0:numel(S.pc1raw)-1)'/30; end
        if isfield(PC,'fps') && PC.fps > 0, S.fps = double(PC.fps); else, S.fps = 1/median(diff(S.t)); end
        efFps.Value = S.fps;
        if isfield(PC,'eigImg'), S.eigImg = PC.eigImg; else, S.eigImg = []; end
        S.diffImg = compute_diff_map(PC);
        pcVar = NaN; if isfield(PC,'varExp') && isfield(PC,'pc'), pcVar = 100*PC.varExp(PC.pc); end
        lblInfo.Text = sprintf('fps=%.2f  frames=%d%s', S.fps, numel(S.pc1raw), ...
            ternary(isnan(pcVar),'',sprintf('  PC var=%.0f%%',pcVar)));

        % preload existing results
        if isfile(sess.matFile)
            try
                saved = load(sess.matFile,'insp_onset_idx','findpeak_params');
                if isfield(saved,'findpeak_params')
                    fp = saved.findpeak_params;
                    if isfield(fp,'camera_fps'),    efFps.Value      = fp.camera_fps;    end
                    if isfield(fp,'minProminence'), efParam(1).Value = fp.minProminence; end
                    if isfield(fp,'minDist_s'),     efParam(2).Value = fp.minDist_s;     end
                    if isfield(fp,'minHeight'),     efParam(3).Value = fp.minHeight;     end
                    if isfield(fp,'minWidth_s'),    efParam(4).Value = fp.minWidth_s;    end
                    if isfield(fp,'snap_win_s'),    efSnap.Value     = fp.snap_win_s;    end
                    if isfield(fp,'inverted'),      cbInvert.Value   = fp.inverted;      end
                end
                cb_updateTrace();
                if isfield(saved,'insp_onset_idx') && ~isempty(saved.insp_onset_idx)
                    S.peaks = saved.insp_onset_idx(:); redraw_peaks();
                end
                lblCurrentFolder.Text = [sess.label '  ✓']; updateNavDisplay();
                return;
            catch
            end
        end
        cb_updateTrace();
    end

    function updateNavDisplay()
        n = numel(S.sessionList);
        nDone = sum(cellfun(@(s) isfile(s.matFile), S.sessionList));
        lblNavStatus.Text = sprintf('%d / %d  folders  |  %d complete', S.sessionIdx, n, nDone);
    end

    function cb_updateTrace(~,~)
        if isempty(S.pc1raw), return; end
        S.fps = efFps.Value;
        b = detrend(S.pc1raw);                 % remove slow linear drift
        sd = std(b,0,'omitnan'); if sd > 0, b = b / sd; end   % z-score
        if cbInvert.Value, b = -b; end
        S.breath = b;
        saveWorkingParams();
        refresh_imgs(); refresh_mainPlot();
    end

    function cb_detect(~,~)
        if isempty(S.breath), uialert(fig,'Load a session first.','No Data'); return; end
        S.peaks = find_peaks_local(S.breath, efFps.Value, ...
            efParam(1).Value, efParam(2).Value, efParam(3).Value, efParam(4).Value);
        saveWorkingParams(); redraw_peaks(); updatePeakCount();
    end

    function cb_axClick(~,~)
        if isempty(S.breath), return; end
        cp = axMain.CurrentPoint; clickT = cp(1,1);
        if clickT < S.t(1) || clickT > S.t(end), return; end
        [~, nearIdx] = min(abs(S.t - clickT));
        switch S.editMode
            case 'add'
                snapFr = max(1, round(efSnap.Value * efFps.Value));
                lo = max(1, nearIdx - snapFr); hi = min(length(S.breath), nearIdx + snapFr);
                [~, imax] = max(S.breath(lo:hi)); peakIdx = lo + imax - 1;
                if isempty(S.peaks) || min(abs(S.peaks - peakIdx)) > 2
                    S.peaks = sort([S.peaks; peakIdx]);
                end
            case 'delete'
                if isempty(S.peaks), return; end
                [~, ci] = min(abs(S.peaks - nearIdx)); S.peaks(ci) = [];
        end
        redraw_peaks(); updatePeakCount();
    end

    function cb_key(~, event)
        switch lower(event.Key)
            case 'n', cb_setNav(); case 'a', cb_setAdd(); case 'd', cb_setDelete();
        end
    end
    function cb_setAdd(~,~)
        S.editMode='add'; lblMode.Text='Mode: + ADD  (click → snap to local max)';
        lblMode.FontColor=[0.30 0.62 1.00]; disableDefaultInteractivity(axMain); axMain.ButtonDownFcn=@cb_axClick;
    end
    function cb_setDelete(~,~)
        S.editMode='delete'; lblMode.Text='Mode: - DELETE  (click → remove nearest)';
        lblMode.FontColor=[1.00 0.35 0.25]; disableDefaultInteractivity(axMain); axMain.ButtonDownFcn=@cb_axClick;
    end
    function cb_setNav(~,~)
        S.editMode='none'; lblMode.Text='Mode: Navigate'; lblMode.FontColor=[0.90 0.80 0.10];
        axMain.ButtonDownFcn=[]; enableDefaultInteractivity(axMain);
    end
    function cb_reset(~,~)
        applyWorkingParams(DEFAULTS); saveWorkingParams(); cb_updateTrace();
    end
    function cb_clearPeaks(~,~)
        S.peaks = []; redraw_peaks(); updatePeakCount();
    end

    function applyWorkingParams(p)
        if isfield(p,'minProm'),   efParam(1).Value = p.minProm;   end
        if isfield(p,'minDist'),   efParam(2).Value = p.minDist;   end
        if isfield(p,'minHeight'), efParam(3).Value = p.minHeight; end
        if isfield(p,'minWidth'),  efParam(4).Value = p.minWidth;  end
        if isfield(p,'snapWin'),   efSnap.Value     = p.snapWin;   end
        if isfield(p,'invert'),    cbInvert.Value   = p.invert;    end
    end
    function saveWorkingParams()
        wp.minProm=efParam(1).Value; wp.minDist=efParam(2).Value;
        wp.minHeight=efParam(3).Value; wp.minWidth=efParam(4).Value;
        wp.snapWin=efSnap.Value; wp.invert=cbInvert.Value;
        setpref('BreathingPeakGUI_PC1','workingParams',wp);
    end

    function cb_save(~,~)
        if isempty(S.folderPath), uialert(fig,'No session loaded.','Error'); return; end
        insp_onset_idx = S.peaks;
        insp_onsets_t = []; insp_onsets_amp = []; insp_onsets_train = zeros(size(S.breath));
        if ~isempty(S.peaks) && ~isempty(S.t)
            insp_onsets_t = S.t(S.peaks); insp_onsets_amp = S.breath(S.peaks);
            insp_onsets_train(S.peaks) = 1;
        end
        breath = S.breath; t_breath = S.t; %#ok<NASGU>
        findpeak_params.minProminence = efParam(1).Value;
        findpeak_params.minDist_s     = efParam(2).Value;
        findpeak_params.minHeight     = efParam(3).Value;
        findpeak_params.minWidth_s    = efParam(4).Value;
        findpeak_params.snap_win_s    = efSnap.Value;
        findpeak_params.camera_fps    = efFps.Value;
        findpeak_params.inverted      = cbInvert.Value;
        findpeak_params.source        = 'svd_pc1';
        findpeak_params.pc1_file      = S.sessionList{S.sessionIdx}.pc1Mat;
        findpeak_params.saved_at      = datestr(now); %#ok<TNOW1,DATST>
        sess = S.sessionList{S.sessionIdx};
        save(sess.matFile,'insp_onset_idx','insp_onsets_t','insp_onsets_amp', ...
            'insp_onsets_train','breath','t_breath','findpeak_params');
        saveScreenshots(sess.figBase);
        uialert(fig, sprintf('Saved %d peaks →\n%s\n+ .fig and .png', numel(S.peaks), sess.matFile), ...
            'Saved','Icon','success');
        updateNavDisplay();
    end

    function saveScreenshots(outBase)
        hFig = figure('Visible','off','Color','white','Position',[100 100 1400 460]);
        ax = axes(hFig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        copyobj(axMain.Children, ax);
        ax.XLim = axMain.XLim; ax.YLim = axMain.YLim;
        ax.XLabel.String = axMain.XLabel.String; ax.YLabel.String = axMain.YLabel.String;
        ax.Title.String = [strrep(S.label,'_',' ') '  —  PC1 breathing'];
        ax.Title.Interpreter='none'; ax.XLabel.Interpreter='none'; ax.YLabel.Interpreter='none';
        savefig(hFig, [outBase '.fig']);
        exportgraphics(hFig, [outBase '.png'], 'Resolution',150);
        close(hFig);
    end

%% ═══════════════════════════════ PLOT HELPERS ═══════════════════════════
    function refresh_imgs()
        cla(axEig);
        if isempty(S.eigImg)
            title(axEig,'(no eigenimage)','Interpreter','none');
        else
            imagesc(axEig, S.eigImg); axis(axEig,'image'); axEig.XTick=[]; axEig.YTick=[];
            colormap(axEig, parula); title(axEig,'PC1 eigenimage','Interpreter','none');
        end
        cla(axDiff);
        if isempty(S.diffImg)
            title(axDiff,'(no diff map)','Interpreter','none');
        else
            imagesc(axDiff, S.diffImg); axis(axDiff,'image'); axDiff.XTick=[]; axDiff.YTick=[];
            colormap(axDiff, parula);
            title(axDiff,'top100 - bottom100 (inspiration - baseline)','Interpreter','none');
        end
    end
    function refresh_mainPlot()
        cla(axMain); hold(axMain,'on');
        if isempty(S.t) || isempty(S.breath), return; end
        h = plot(axMain, S.t, S.breath, 'Color',[0 0 0], 'LineWidth',0.9); h.HitTest='off';
        h2 = yline(axMain, 0, '--','Color',[0.5 0.5 0.5],'LineWidth',0.6); h2.HitTest='off';
        hPeakSc = scatter(axMain, [], [], 75, 'v','filled', ...
            'MarkerFaceColor',[1.00 0.28 0.28],'MarkerEdgeColor','none','Tag','peak_marker');
        hPeakSc.HitTest='off';
        xlim(axMain,[S.t(1) S.t(end)]);
        redraw_peaks(); updatePeakCount();
    end
    function redraw_peaks()
        valid = [];
        if ~isempty(S.peaks) && ~isempty(S.breath)
            valid = S.peaks(S.peaks >= 1 & S.peaks <= length(S.breath));
        end
        if isgraphics(hPeakSc)
            if isempty(valid), set(hPeakSc,'XData',[],'YData',[]);
            else, set(hPeakSc,'XData',S.t(valid),'YData',S.breath(valid)); end
        end
        updatePeakCount();
    end
    function updatePeakCount()
        lblCount.Text = sprintf('Peaks: %d', numel(S.peaks));
    end
end  % breathing_peak_gui_pc1

%% ═══════════════════════════════ LOCAL HELPERS ══════════════════════════
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function dimg = compute_diff_map(PC)
% Inspiration-baseline difference map (top-100 highest-signal frames minus
% bottom-100). Uses a cached full-res PC.diffImg if present; otherwise
% reconstructs it from the SVD factors: U * (sv .* (mean V_top - mean V_bot)).
dimg = [];
if isfield(PC,'diffImg') && ~isempty(PC.diffImg), dimg = double(PC.diffImg); return; end
if ~all(isfield(PC,{'U','V','sv','breathTrace','eigImg'})), return; end
bt = double(PC.breathTrace(:));
if numel(bt) < 10 || std(bt) == 0, return; end
z = (bt-mean(bt))/std(bt); if mean(z.^3) < 0, bt = -bt; end   % inspiration = high side
nf = min(100, floor(numel(bt)/2));
[~, ord] = sort(bt, 'descend'); top = ord(1:nf); bot = ord(end-nf+1:end);
U = double(PC.U); sv = double(PC.sv(:)); V = double(PC.V);
k = min([size(U,2), numel(sv), size(V,2)]);
dV = (mean(V(top,1:k),1) - mean(V(bot,1:k),1))';
dimg = reshape(U(:,1:k) * (sv(1:k).*dV), size(PC.eigImg));
end

function sessions = scanForSessions(masterFolder)
%SCANFORSESSIONS  Recursively find every breath_pc1.mat under masterFolder.
sessions = {};
hits = dir(fullfile(masterFolder,'**','breath_pc1.mat'));
if isempty(hits), return; end
for i = 1:numel(hits)
    h = hits(i);
    [~, leaf] = fileparts(h.folder);
    sessions{end+1} = struct( ...                                   %#ok<AGROW>
        'folder',  h.folder, ...
        'pc1Mat',  fullfile(h.folder, h.name), ...
        'label',   leaf, ...
        'matFile', fullfile(h.folder, 'breath_peak_pc1.mat'), ...
        'figBase', fullfile(h.folder, 'breath_peak_pc1_figure'));
end
[~, ord] = sort(cellfun(@(s) s.folder, sessions, 'UniformOutput',false));
sessions = sessions(ord);
end

function locs = find_peaks_local(sig, fps, minProm, minDistS, minH, minWS)
%FIND_PEAKS_LOCAL  Local-max finder w/ prominence/distance/height/width filters.
sig = sig(:); n = length(sig);
isMax = false(n,1);
isMax(2:end-1) = sig(2:end-1) > sig(1:end-2) & sig(2:end-1) > sig(3:end);
locs = find(isMax);
if isempty(locs), return; end
if minH ~= 0
    locs = locs(sig(locs) >= minH); if isempty(locs), return; end
end
if minProm > 0
    prom = compute_prominence(sig, locs); locs = locs(prom >= minProm);
    if isempty(locs), return; end
end
minDistFr = max(1, round(minDistS * fps));
if minDistFr > 1 && length(locs) > 1
    [~, ord] = sort(sig(locs),'descend'); sorted = locs(ord);
    keep = true(length(sorted),1);
    for i = 1:length(sorted)
        if ~keep(i), continue; end
        tooClose = abs(sorted - sorted(i)) < minDistFr; tooClose(i) = false;
        keep(tooClose) = false;
    end
    locs = sort(sorted(keep));
end
if isempty(locs), return; end
if minWS > 0
    minWFr = max(1, round(minWS * fps)); prom = compute_prominence(sig, locs);
    keep = true(length(locs),1);
    for i = 1:length(locs)
        pk = locs(i); half_lvl = sig(pk) - prom(i)/2;
        lc = find(sig(1:pk) < half_lvl, 1,'last'); if isempty(lc), lc = 1; end
        rc = find(sig(pk:end) < half_lvl, 1,'first'); if isempty(rc), rc = n - pk + 1; end
        keep(i) = ((pk - lc) + (rc - 1)) >= minWFr;
    end
    locs = locs(keep);
end
locs = locs(:);
end

function prom = compute_prominence(sig, locs)
n = length(sig); np = length(locs); prom = zeros(np,1);
for i = 1:np
    pk = locs(i);
    left_higher = locs(sig(locs) > sig(pk) & locs < pk);
    li = 1; if ~isempty(left_higher), li = left_higher(end); end
    right_higher = locs(sig(locs) > sig(pk) & locs > pk);
    ri = n; if ~isempty(right_higher), ri = right_higher(1); end
    prom(i) = sig(pk) - max(min(sig(li:pk)), min(sig(pk:ri)));
end
end
