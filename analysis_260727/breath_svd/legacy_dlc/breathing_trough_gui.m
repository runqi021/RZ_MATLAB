function breathing_trough_gui()
%BREATHING_TROUGH_GUI  Interactive trough detector on the breath waveform
%                     already produced by breathing_peak_gui.m.
%
%   This GUI does NOT rebuild the breath signal from DLC. It loads each
%   session's existing *_breath_peak_data.mat and uses BP.breath and
%   BP.t_breath verbatim, plus camera_fps and the source CSV stem. You
%   only edit TROUGH thresholds + the troughs themselves; everything
%   waveform-related stays locked to the peak data.
%
%   Sessions = folders containing a *_breath_peak_data.mat file.
%
%   Output (per session, alongside the peak .mat):
%       [csvStem]_breath_trough_data.mat
%           exp_trough_idx, exp_troughs_t, exp_troughs_amp,
%           exp_troughs_train, breath, t_breath, findtrough_params,
%           peak_file (path to the peak .mat that produced this waveform)
%
%   Requires: MATLAB R2021a+

%% ── Shared state ─────────────────────────────────────────────────────────
S.breath       = [];     % [F×1] waveform from peak .mat (verbatim)
S.t            = [];     % [F×1] seconds (from peak .mat)
S.fps          = NaN;
S.peakIdx      = [];     % [P×1] peak frame indices (reference overlay)
S.troughs      = [];     % [K×1] trough frame indices (editable)
S.editMode     = 'none';
S.masterFolder = '';
S.sessionList  = {};
S.sessionIdx   = 0;
S.peakInfo     = '';     % multiline summary of how the waveform was built

DEFAULTS.steep_frac     = 0.30;
DEFAULTS.max_lookback_s = 1.00;
DEFAULTS.smooth_win_s   = 0.05;
DEFAULTS.snapWin        = 0.10;

%% ── Figure ──────────────────────────────────────────────────────────────
fig = uifigure('Name','Breathing Inspiration-Start Sorter (foot of rising flank)', ...
               'Position',[40 30 1560 930], ...
               'WindowState','maximized');

rootGL = uigridlayout(fig, [1 2]);
rootGL.ColumnWidth   = {340, '1x'};
rootGL.Padding       = [6 6 6 6];
rootGL.ColumnSpacing = 8;

%% ── Left: control panel ─────────────────────────────────────────────────
ctrlPan = uipanel(rootGL, 'Title','Controls', 'FontSize',12,'FontWeight','bold');
ctrlPan.Layout.Column = 1;

nCtrlRows = 26;
cGL = uigridlayout(ctrlPan, [nCtrlRows 2]);
cGL.RowHeight   = repmat({'fit'}, 1, nCtrlRows);
cGL.ColumnWidth = {'fit','1x'};
cGL.Padding     = [8 8 8 8];
cGL.RowSpacing  = 4;

r = 0;

% ── Session Navigation ─────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Session Navigation ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnLoad = uibutton(cGL,'Text','Browse Master Folder','ButtonPushedFcn',@cb_loadMaster);
       btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];
r=r+1; lblMasterPath = uilabel(cGL,'Text','No folder selected','WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]);
       lblMasterPath.Layout.Row=r; lblMasterPath.Layout.Column=[1 2];
r=r+1; lblCurrentFolder = uilabel(cGL,'Text','—','FontWeight','bold','WordWrap','on','FontSize',10);
       lblCurrentFolder.Layout.Row=r; lblCurrentFolder.Layout.Column=[1 2];
r=r+1; lblNavStatus = uilabel(cGL,'Text','0 / 0  folders  |  0 complete','FontSize',10,'FontColor',[0.25 0.65 0.25]);
       lblNavStatus.Layout.Row=r; lblNavStatus.Layout.Column=[1 2];
r=r+1; btnPrev = uibutton(cGL,'Text','← Prev','ButtonPushedFcn',@cb_prevSession);
       btnPrev.Layout.Row=r; btnPrev.Layout.Column=1;
       btnNext = uibutton(cGL,'Text','Next →','ButtonPushedFcn',@cb_nextSession);
       btnNext.Layout.Row=r; btnNext.Layout.Column=2;

% ── Peak-data summary (READ-ONLY; locked to the peak .mat) ────────────
r=r+1; h=uilabel(cGL,'Text','── Source waveform (locked) ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; taInfo = uitextarea(cGL,'Value',{'(no session loaded)'}, ...
            'Editable','off','FontSize',9,'FontName','Consolas');
       taInfo.Layout.Row=r; taInfo.Layout.Column=[1 2];

% ── Inspiration-Start Detection ──────────────────────────────────────
% Algorithm: foot of the rising flank.  For each frame i, compute slope
% d(i) = (breath(i+1) - breath(i)) * fps.  Mark i as an "insp start" when
% d(i) crosses from below the slope threshold to above.  Frame i is the
% LAST sample sitting on the baseline before the steep rise (matches the
% red-circled points in user's reference image).
r=r+1; h=uilabel(cGL,'Text','── Inspiration Start Detection ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

paramLabels = {'Steep frac of max slope:','Max lookback (s):','Smooth win (s):'};
paramDefs   = {0.30, 1.00, 0.05};
efParam     = gobjects(3,1);
for p = 1:3
    r=r+1;
    h=uilabel(cGL,'Text',paramLabels{p});
    h.Layout.Row=r; h.Layout.Column=1;
    efParam(p) = uieditfield(cGL,'numeric','Value',paramDefs{p},'Limits',[-Inf Inf], ...
        'ValueChangedFcn',@(~,~) saveWorkingParams());
    efParam(p).Layout.Row=r; efParam(p).Layout.Column=2;
end
r=r+1; h=uilabel(cGL,'Text','Snap window (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efSnap = uieditfield(cGL,'numeric','Value',DEFAULTS.snapWin,'Limits',[0 10], ...
           'ValueChangedFcn',@(~,~) saveWorkingParams());
       efSnap.Layout.Row=r; efSnap.Layout.Column=2;
r=r+1; btnDetect = uibutton(cGL,'Text','Detect Insp Starts','ButtonPushedFcn',@cb_detect, ...
           'BackgroundColor',[0.18 0.40 0.85],'FontColor','white');
       btnDetect.Layout.Row=r; btnDetect.Layout.Column=[1 2];
r=r+1; btnGap = uibutton(cGL,'Text','Find Next Gap →','ButtonPushedFcn',@cb_find_gap, ...
           'BackgroundColor',[0.85 0.55 0.10],'FontColor','white');
       btnGap.Layout.Row=r; btnGap.Layout.Column=[1 2];

% ── Manual Editing ────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Manual Editing ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; lblMode = uilabel(cGL,'Text','Mode: Navigate','FontWeight','bold','FontSize',11, ...
           'FontColor',[0.90 0.80 0.10]);
       lblMode.Layout.Row=r; lblMode.Layout.Column=[1 2];
r=r+1; btnAdd = uibutton(cGL,'Text','+ Add Mode','ButtonPushedFcn',@cb_setAdd, ...
           'BackgroundColor',[0.20 0.42 0.85],'FontColor','white');
       btnAdd.Layout.Row=r; btnAdd.Layout.Column=1;
       btnDel = uibutton(cGL,'Text','- Delete Mode','ButtonPushedFcn',@cb_setDelete, ...
           'BackgroundColor',[0.82 0.22 0.18],'FontColor','white');
       btnDel.Layout.Row=r; btnDel.Layout.Column=2;
r=r+1; btnNav = uibutton(cGL,'Text','Navigate (N)','ButtonPushedFcn',@cb_setNav);
       btnNav.Layout.Row=r; btnNav.Layout.Column=1;
       btnClear = uibutton(cGL,'Text','Clear All Troughs','ButtonPushedFcn',@cb_clearTroughs, ...
           'BackgroundColor',[0.40 0.10 0.10],'FontColor','white');
       btnClear.Layout.Row=r; btnClear.Layout.Column=2;
r=r+1; lblCount = uilabel(cGL,'Text','Insp Starts: 0  |  Peaks: 0','FontWeight','bold','FontSize',12);
       lblCount.Layout.Row=r; lblCount.Layout.Column=[1 2];

% ── Save ──────────────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Save ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnSave = uibutton(cGL,'Text','Save breath_insp_start_data.mat','ButtonPushedFcn',@cb_save, ...
           'BackgroundColor',[0.05 0.45 0.65],'FontColor','white');
       btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];
r=r+1; btnReset = uibutton(cGL,'Text','Reset Trough Params','ButtonPushedFcn',@cb_reset, ...
           'BackgroundColor',[0.30 0.30 0.30],'FontColor','white');
       btnReset.Layout.Row=r; btnReset.Layout.Column=[1 2];

% NOTE: do NOT restore persisted prefs on launch. Code defaults win every
% time, so stale values from prior runs can't override them. Per-session
% changes are still written via saveWorkingParams so they survive within
% a running GUI instance (and can be inspected via getpref if needed).
if ispref('BreathingTroughGUI','workingParams')
    try, rmpref('BreathingTroughGUI','workingParams'); catch, end
end

%% ── Right: plot panel ───────────────────────────────────────────────────
rightGL = uigridlayout(rootGL, [1 1]);
rightGL.Layout.Column = 2;
rightGL.Padding       = [0 0 0 0];

axMain = uiaxes(rightGL);
axMain.Title.String         = 'Breath waveform (from peak .mat)  |  ▲ peak ref  |  ▼ insp start (foot of rising flank)';
axMain.Title.Interpreter    = 'none';
axMain.XLabel.String        = 'Time (s)';
axMain.XLabel.Interpreter   = 'none';
axMain.YLabel.String        = 'Signal (z-score)';
axMain.YLabel.Interpreter   = 'none';
hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');

hTroughSc = gobjects(0);
hPeakSc   = gobjects(0);

fig.KeyPressFcn = @cb_key;

%% ═══════════════════════════════════════════════════════════════════════
%% CALLBACKS
%% ═══════════════════════════════════════════════════════════════════════

    function cb_loadMaster(~,~)
        master = uigetdir(S.masterFolder,'Select Master Folder');
        if isequal(master,0), return; end
        S.masterFolder = master;
        sessions = scanForSessions(master);
        if isempty(sessions)
            uialert(fig,'No *_breath_peak_data.mat found in any subfolder.','No Sessions'); return;
        end
        S.sessionList = sessions; S.sessionIdx = 1;
        [~, mName] = fileparts(master);
        lblMasterPath.Text = mName;
        loadCurrentSession();
    end

    function cb_prevSession(~,~)
        if S.sessionIdx <= 1, return; end
        S.sessionIdx = S.sessionIdx - 1;
        loadCurrentSession();
    end

    function cb_nextSession(~,~)
        if S.sessionIdx >= numel(S.sessionList), return; end
        S.sessionIdx = S.sessionIdx + 1;
        loadCurrentSession();
    end

    function loadCurrentSession()
        sess     = S.sessionList{S.sessionIdx};
        S.troughs = [];

        lblCurrentFolder.Text = sess.csvStem;
        updateNavDisplay();

        % --- load the peak .mat (waveform + peaks + how-it-was-built) ---
        try
            PK = load(sess.peakMat);
        catch ME
            uialert(fig, ME.message,'Load Error'); return;
        end
        assert(isfield(PK,'breath') && isfield(PK,'t_breath'), ...
               'Peak .mat missing "breath" or "t_breath".');

        S.breath = double(PK.breath(:));
        S.t      = double(PK.t_breath(:));
        if isfield(PK,'findpeak_params') && isfield(PK.findpeak_params,'camera_fps')
            S.fps = PK.findpeak_params.camera_fps;
        elseif numel(S.t) > 1
            S.fps = 1 / median(diff(S.t));
        else
            S.fps = NaN;
        end
        if isfield(PK,'insp_onset_idx')
            S.peakIdx = sort(double(PK.insp_onset_idx(:)));
        else
            S.peakIdx = [];
        end

        % --- summary string of how the waveform was built ---
        info = {sprintf('csv  : %s', sess.csvStem)};
        if isfield(PK,'findpeak_params')
            fp = PK.findpeak_params;
            if isfield(fp,'camera_fps'), info{end+1} = sprintf('fps  : %.3f', fp.camera_fps); end
            if isfield(fp,'dot_selection')
                info{end+1} = sprintf('dots : %s', strjoin(fp.dot_selection, ', '));
            end
            if isfield(fp,'coord_types')
                info{end+1} = sprintf('coord: %s', strjoin(fp.coord_types, ', '));
            end
            if isfield(fp,'combine_method'), info{end+1} = sprintf('combine: %s', fp.combine_method); end
            if isfield(fp,'inverted'),        info{end+1} = sprintf('inverted: %d', fp.inverted); end
            if isfield(fp,'likelihood_thr'),  info{end+1} = sprintf('lk_thr: %.2f', fp.likelihood_thr); end
        end
        info{end+1} = sprintf('peaks: %d', numel(S.peakIdx));
        taInfo.Value = info;
        S.peakInfo   = strjoin(info, newline);

        % --- preload existing trough data if present ---
        if isfile(sess.matFile)
            try
                saved = load(sess.matFile,'insp_start_idx','findstart_params');
                if isfield(saved,'findstart_params')
                    fsp = saved.findstart_params;
                    if isfield(fsp,'slope_frac'),       efParam(1).Value = fsp.slope_frac;       end
                    if isfield(fsp,'max_lookback_s'),   efParam(2).Value = fsp.max_lookback_s;   end
                    if isfield(fsp,'smooth_win_s'),     efParam(3).Value = fsp.smooth_win_s;     end
                    if isfield(fsp,'snap_win_s'),       efSnap.Value     = fsp.snap_win_s;       end
                end
                if isfield(saved,'insp_start_idx')
                    S.troughs = double(saved.insp_start_idx(:));
                end
                lblCurrentFolder.Text = [sess.csvStem '  ✓'];
                updateNavDisplay();
            catch
            end
        end

        refresh_mainPlot();
    end

    function updateNavDisplay()
        n = numel(S.sessionList);
        nDone = sum(cellfun(@(s) isfile(s.matFile), S.sessionList));
        lblNavStatus.Text = sprintf('%d / %d  folders  |  %d complete', S.sessionIdx, n, nDone);
    end

    function cb_detect(~,~)
        if isempty(S.breath) || isempty(S.peakIdx)
            uialert(fig,'Need breath waveform AND validated peaks.','No Data'); return;
        end
        % LINE-FIT method: fit a line through the steepest portion of each
        % rising flank, extrapolate down to the baseline, take the
        % intersection as the foot. Geometry of the rise itself sets the
        % answer -- no threshold tuning of "how shallow is shallow".
        S.troughs = find_feet_lineFit(S.breath, S.fps, S.peakIdx, ...
            efParam(1).Value, efParam(2).Value, efParam(3).Value);
        saveWorkingParams();
        redraw_marks(); updateCount();
    end

    function cb_find_gap(~,~)
        % Find next inter-peak interval that has no insp-start, then pan
        % the x-axis to show it. "Next" = first one whose midpoint is to
        % the right of the current view's center. Wraps to start if none.
        if numel(S.peakIdx) < 2 || isempty(S.t)
            lblMode.Text = 'no peaks loaded'; return;
        end
        gaps = [];
        for k = 1:numel(S.peakIdx)-1
            n_in = sum(S.troughs > S.peakIdx(k) & S.troughs < S.peakIdx(k+1));
            if n_in == 0, gaps(end+1) = k; end %#ok<AGROW>
        end
        if isempty(gaps)
            lblMode.Text = 'No gaps -- 1/cycle ✓';
            lblMode.FontColor = [0.10 0.55 0.10];
            return;
        end
        % which gap to jump to: first gap whose midpoint lies right of current center
        cur_center = mean(axMain.XLim);
        gap_mid    = arrayfun(@(g) (S.t(S.peakIdx(g)) + S.t(S.peakIdx(g+1)))/2, gaps);
        nextGap    = find(gap_mid > cur_center, 1, 'first');
        if isempty(nextGap), nextGap = 1; end   % wrap to start
        g = gaps(nextGap);
        p1 = S.peakIdx(g); p2 = S.peakIdx(g+1);
        t_lo = S.t(p1); t_hi = S.t(p2);
        % show 4x the gap width so neighboring cycles are visible too
        pad = max(0.5, (t_hi - t_lo) * 1.5);
        axMain.XLim = [max(S.t(1), t_lo - pad), min(S.t(end), t_hi + pad)];
        lblMode.Text = sprintf('Gap %d / %d -- between peaks %d & %d', nextGap, numel(gaps), g, g+1);
        lblMode.FontColor = [0.85 0.45 0.10];
    end

    function cb_axClick(~,~)
        if isempty(S.breath), return; end
        cp     = axMain.CurrentPoint;
        clickT = cp(1,1);
        if clickT < S.t(1) || clickT > S.t(end), return; end
        [~, nearIdx] = min(abs(S.t - clickT));

        switch S.editMode
            case 'add'
                % Snap: scan within the window for the FOOT of the rising flank
                % (= first frame within window where d(breath)/dt > slope_thresh,
                %  then step back one frame so we sit on the baseline sample).
                snapFr   = max(1, round(efSnap.Value * S.fps));
                lo       = max(1, nearIdx - snapFr);
                hi       = min(length(S.breath)-1, nearIdx + snapFr);
                slopeTh  = efParam(1).Value;
                trIdx    = [];
                if hi > lo
                    d_win = (S.breath(lo+1:hi+1) - S.breath(lo:hi)) * S.fps;
                    rise  = find(d_win > slopeTh, 1, 'first');
                    if ~isempty(rise), trIdx = lo + rise - 1; end
                end
                if isempty(trIdx)
                    % fallback: just use the local min in the window
                    [~, imin] = min(S.breath(lo:hi));
                    trIdx = lo + imin - 1;
                end
                if isempty(S.troughs) || min(abs(S.troughs - trIdx)) > 2
                    S.troughs = sort([S.troughs; trIdx]);
                end
            case 'delete'
                if isempty(S.troughs), return; end
                [~, ci] = min(abs(S.troughs - nearIdx));
                S.troughs(ci) = [];
        end
        redraw_marks(); updateCount();
    end

    function cb_key(~, event)
        switch lower(event.Key)
            case 'n', cb_setNav();
            case 'a', cb_setAdd();
            case 'd', cb_setDelete();
        end
    end

    function cb_setAdd(~,~)
        S.editMode = 'add';
        lblMode.Text = 'Mode: + ADD  (click → snap to local MIN)';
        lblMode.FontColor = [0.30 0.62 1.00];
        disableDefaultInteractivity(axMain);
        axMain.ButtonDownFcn = @cb_axClick;
    end
    function cb_setDelete(~,~)
        S.editMode = 'delete';
        lblMode.Text = 'Mode: - DELETE  (click → remove nearest)';
        lblMode.FontColor = [1.00 0.35 0.25];
        disableDefaultInteractivity(axMain);
        axMain.ButtonDownFcn = @cb_axClick;
    end
    function cb_setNav(~,~)
        S.editMode = 'none';
        lblMode.Text = 'Mode: Navigate';
        lblMode.FontColor = [0.90 0.80 0.10];
        axMain.ButtonDownFcn = [];
        enableDefaultInteractivity(axMain);
    end

    function cb_reset(~,~)
        applyWorkingParams(DEFAULTS);
        saveWorkingParams();
    end

    function applyWorkingParams(p)
        if isfield(p,'steep_frac'),     efParam(1).Value = p.steep_frac;     end
        if isfield(p,'max_lookback_s'), efParam(2).Value = p.max_lookback_s; end
        if isfield(p,'smooth_win_s'),   efParam(3).Value = p.smooth_win_s;   end
        if isfield(p,'snapWin'),        efSnap.Value     = p.snapWin;        end
    end

    function saveWorkingParams()
        wp.steep_frac     = efParam(1).Value;
        wp.max_lookback_s = efParam(2).Value;
        wp.smooth_win_s   = efParam(3).Value;
        wp.snapWin        = efSnap.Value;
        setpref('BreathingTroughGUI','workingParams',wp);
    end

    function cb_clearTroughs(~,~)
        S.troughs = []; redraw_marks(); updateCount();
    end

    function cb_save(~,~)
        if isempty(S.sessionList), uialert(fig,'No session loaded.','Error'); return; end
        sess = S.sessionList{S.sessionIdx};

        insp_start_idx    = S.troughs;
        insp_starts_t     = []; insp_starts_amp = [];
        insp_starts_train = zeros(size(S.breath));
        if ~isempty(S.troughs) && ~isempty(S.t)
            insp_starts_t   = S.t(S.troughs);
            insp_starts_amp = S.breath(S.troughs);
            insp_starts_train(S.troughs) = 1;
        end
        breath   = S.breath;       %#ok<NASGU>  copied verbatim from peak .mat
        t_breath = S.t;            %#ok<NASGU>

        findstart_params.slope_frac     = efParam(1).Value;
        findstart_params.max_lookback_s = efParam(2).Value;
        findstart_params.smooth_win_s   = efParam(3).Value;
        findstart_params.snap_win_s     = efSnap.Value;
        findstart_params.camera_fps     = S.fps;
        findstart_params.peak_file      = sess.peakMat;
        findstart_params.saved_at       = datestr(now);     %#ok<TNOW1,DATST>

        peak_file = sess.peakMat;    %#ok<NASGU>

        save(sess.matFile, ...
             'insp_start_idx','insp_starts_t','insp_starts_amp','insp_starts_train', ...
             'breath','t_breath','findstart_params','peak_file');

        saveScreenshots(sess.figBase);

        uialert(fig, sprintf('Saved %d troughs →\n%s\n+ .fig and .png', ...
                numel(S.troughs), sess.matFile), 'Saved','Icon','success');
        updateNavDisplay();
    end

    function saveScreenshots(outBase)
        hFig = figure('Visible','off','Color','white','Position',[100 100 1400 460]);
        ax = axes(hFig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        copyobj(axMain.Children, ax);
        ax.XLim = axMain.XLim; ax.YLim = axMain.YLim;
        ax.XLabel.String = axMain.XLabel.String;
        ax.YLabel.String = axMain.YLabel.String;
        ax.Title.String  = axMain.Title.String;
        ax.Title.Interpreter='none'; ax.XLabel.Interpreter='none'; ax.YLabel.Interpreter='none';
        savefig(hFig,  [outBase '.fig']);
        exportgraphics(hFig, [outBase '.png'], 'Resolution',150);
        close(hFig);
    end

%% ═══════════════════════════════════════════════════════════════════════
%% PLOT HELPERS
%% ═══════════════════════════════════════════════════════════════════════

    function refresh_mainPlot()
        cla(axMain); hold(axMain,'on');
        if isempty(S.t) || isempty(S.breath), return; end
        h = plot(axMain, S.t, S.breath, 'Color',[0 0 0], 'LineWidth',0.9);
        h.HitTest = 'off';
        h2 = yline(axMain, 0, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.6);
        h2.HitTest = 'off';
        % Peak reference (orange ▲, not editable here)
        hPeakSc = scatter(axMain, [], [], 25, '^', 'filled', ...
            'MarkerFaceColor',[0.95 0.50 0.05], 'MarkerEdgeColor','none', 'Tag','peak_marker');
        hPeakSc.HitTest = 'off';
        % Trough markers (blue ▼, editable)
        hTroughSc = scatter(axMain, [], [], 75, 'v', 'filled', ...
            'MarkerFaceColor',[0 0.55 1], 'MarkerEdgeColor','none', 'Tag','trough_marker');
        hTroughSc.HitTest = 'off';
        xlim(axMain,[S.t(1) S.t(end)]);
        redraw_marks(); updateCount();
    end

    function redraw_marks()
        % Peaks
        if isgraphics(hPeakSc)
            v = S.peakIdx;
            v(v < 1 | v > length(S.breath)) = [];
            if isempty(v)
                set(hPeakSc, 'XData',[], 'YData',[]);
            else
                set(hPeakSc, 'XData',S.t(v), 'YData',S.breath(v));
            end
        end
        % Troughs
        if isgraphics(hTroughSc)
            v = S.troughs;
            v(v < 1 | v > length(S.breath)) = [];
            if isempty(v)
                set(hTroughSc, 'XData',[], 'YData',[]);
            else
                set(hTroughSc, 'XData',S.t(v), 'YData',S.breath(v));
            end
        end
        updateCount();
    end

    function updateCount()
        nGaps = 0; nExtra = 0;
        if numel(S.peakIdx) >= 2 && ~isempty(S.troughs)
            for k = 1:numel(S.peakIdx)-1
                n_in = sum(S.troughs > S.peakIdx(k) & S.troughs < S.peakIdx(k+1));
                if n_in == 0, nGaps  = nGaps  + 1; end
                if n_in > 1,  nExtra = nExtra + 1; end
            end
        end
        if nGaps + nExtra == 0
            lblCount.Text = sprintf('Insp Starts: %d  |  Peaks: %d  |  ✓ 1/cycle', ...
                numel(S.troughs), numel(S.peakIdx));
            lblCount.FontColor = [0.10 0.55 0.10];
        else
            lblCount.Text = sprintf('Insp Starts: %d  |  Peaks: %d  |  gaps: %d  extra: %d', ...
                numel(S.troughs), numel(S.peakIdx), nGaps, nExtra);
            lblCount.FontColor = [0.85 0.20 0.10];
        end
    end

end  % breathing_trough_gui


%% ═══════════════════════════════════════════════════════════════════════
%% LOCAL HELPERS
%% ═══════════════════════════════════════════════════════════════════════
function feet = find_feet_lineFit(breath, fps, peak_idx, steep_frac, max_lookback_s, smooth_win_s)
% FIND_FEET_LINEFIT  Geometric foot detector.
%
%   For each validated peak, fit a straight line to the steep portion of
%   the rising flank (frames where the smoothed slope is above
%   `steep_frac` of the max slope on that flank), then extrapolate the
%   line backward until it hits the local baseline.  The intersection
%   point is the foot.
%
%   Rationale: the rise is approximately linear once it starts moving, and
%   the foot is by definition "where that linear rise meets the resting
%   baseline".  This is just geometry -- no slope threshold tuning.
%
%   steep_frac     : fraction of max-slope above which frames count as
%                    part of the "rising line" used for the fit (e.g. 0.30)
%   max_lookback_s : how far back to look from each peak
%   smooth_win_s   : moving-mean window on slope before picking steep frames
%
%   Returns one foot frame index per peak (NaN if not solvable).
    breath          = breath(:);
    feet            = nan(size(peak_idx));
    max_lookback_fr = max(3, round(max_lookback_s * fps));
    smooth_fr       = max(1, round(smooth_win_s * fps));
    for k = 1:numel(peak_idx)
        p  = peak_idx(k);
        lo = max(1, p - max_lookback_fr);
        if p <= lo+2, continue; end

        idx_w = (lo:p)';                                  % window indices
        seg   = breath(idx_w);                            % signal there
        d     = (seg(2:end) - seg(1:end-1)) * fps;        % slope, length p-lo
        d     = movmean(d, smooth_fr);
        max_d = max(d);
        if ~isfinite(max_d) || max_d <= 0
            feet(k) = lo; continue;
        end

        % steep portion = slopes above steep_frac * max
        steep = d > steep_frac * max_d;

        % take only the LATEST contiguous steep run (the one leading INTO
        % the peak), in case noise produced earlier steep blips.
        s_idx = find(steep);
        if isempty(s_idx), feet(k) = lo; continue; end
        gaps  = find(diff(s_idx) > 1);
        if isempty(gaps)
            run = s_idx;
        else
            run = s_idx(gaps(end)+1:end);
        end
        if numel(run) < 2
            % too short to fit a line: fall back to the run's first frame
            feet(k) = idx_w(run(1)); continue;
        end

        % fit line on the steep run: signal(frame) ~ m*frame + c
        % frames here are absolute (in 'breath' indexing). Use midpoint
        % between consecutive samples to align with the slope grid.
        frames_fit = idx_w(run) + 0.5;
        signal_fit = (seg(run) + seg(run+1)) / 2;
        pf      = polyfit(frames_fit, signal_fit, 1);
        m_slope = pf(1);
        c_int   = pf(2);
        if abs(m_slope) < eps, feet(k) = lo; continue; end

        % baseline = median of the pre-steep portion of the window
        pre_end = run(1) - 1;
        if pre_end >= 1
            baseline = median(seg(1:pre_end));
        else
            baseline = seg(1);
        end

        % foot = x where the rising line crosses the baseline
        x_foot  = (baseline - c_int) / m_slope;
        x_foot  = max(lo, min(p-1, x_foot));
        feet(k) = round(x_foot);
    end
    feet = sort(feet(~isnan(feet)));
end

function feet = find_feet_walkback(breath, fps, peak_idx, frac, max_lookback_s, smooth_win_s)
% FIND_FEET_WALKBACK  For each validated peak, walk backwards along the
% rising flank until the (smoothed) slope drops below an adaptive
% threshold = frac * max(slope on that flank).  That last frame is the
% foot (= start of inspiration).  Guarantees one foot per peak.
%
%   peak_idx       : validated peak frame indices (from breath_peak_data)
%   frac           : threshold as fraction of this cycle's max slope
%                     (e.g. 0.10 -> stop walking back when slope falls to
%                      10% of the steepest part of this rise)
%   max_lookback_s : how far back to look from each peak (seconds)
%   smooth_win_s   : moving-mean window for the slope before walking
%
%   Returns a sorted column of foot frame indices (one per peak).
    breath          = breath(:);
    feet            = nan(size(peak_idx));
    max_lookback_fr = max(2, round(max_lookback_s * fps));
    smooth_fr       = max(1, round(smooth_win_s * fps));
    for k = 1:numel(peak_idx)
        p  = peak_idx(k);
        lo = max(1, p - max_lookback_fr);
        if p <= lo+1, continue; end
        % smoothed slope on [lo, p]  (length p-lo)
        seg    = breath(lo:p);
        d      = (seg(2:end) - seg(1:end-1)) * fps;
        d      = movmean(d, smooth_fr);
        max_d  = max(d);
        if ~isfinite(max_d) || max_d <= 0
            feet(k) = lo; continue;
        end
        thresh = frac * max_d;
        % walk back from the peak end of d until slope drops below thresh
        cur = numel(d);
        while cur >= 1 && d(cur) > thresh
            cur = cur - 1;
        end
        if cur >= 1
            feet(k) = lo + cur - 1;   % LAST baseline frame before steep rise
        else
            feet(k) = lo;             % rising flank reached lookback limit
        end
    end
    feet = sort(feet(~isnan(feet)));
end

function trs = enforce_one_per_cycle_latest(starts, peaks)
% Keep at most one insp-start between each pair of consecutive peaks.
% When multiple candidates fall in the same inter-peak interval, keep the
% LATEST one (the foot closest to the next peak's steep rise).
starts = sort(starts(:));
peaks  = sort(peaks(:));
if isempty(peaks) || isempty(starts), trs = starts; return; end
keep = false(size(starts));
pre = starts < peaks(1);
if any(pre), idx = find(pre); keep(idx(end)) = true; end
for k = 1:numel(peaks)-1
    in_int = starts > peaks(k) & starts < peaks(k+1);
    if any(in_int), idx = find(in_int); keep(idx(end)) = true; end
end
post = starts > peaks(end);
if any(post), idx = find(post); keep(idx(end)) = true; end
trs = starts(keep);
end

function trs = enforce_one_per_cycle(troughs, peaks, breath)
% Keep at most one trough between each pair of consecutive peaks. When
% multiple troughs fall in the same inter-peak interval, keep the one with
% the lowest breath value (deepest minimum). Also handles troughs before
% the first peak and after the last peak (one of each kept, if any).
    troughs = sort(troughs(:));
    peaks   = sort(peaks(:));
    if isempty(peaks) || isempty(troughs), trs = troughs; return; end
    keep = false(size(troughs));
    pre = troughs < peaks(1);
    if any(pre)
        idx = find(pre); [~, mi] = min(breath(troughs(idx))); keep(idx(mi)) = true;
    end
    for k = 1:numel(peaks)-1
        in_int = troughs > peaks(k) & troughs < peaks(k+1);
        if any(in_int)
            idx = find(in_int); [~, mi] = min(breath(troughs(idx))); keep(idx(mi)) = true;
        end
    end
    post = troughs > peaks(end);
    if any(post)
        idx = find(post); [~, mi] = min(breath(troughs(idx))); keep(idx(mi)) = true;
    end
    trs = troughs(keep);
end

function locs = find_peaks_local(sig, fps, minProm, minDistS, minH, minWS)
sig  = sig(:); n = length(sig);
isMax = false(n,1);
isMax(2:end-1) = sig(2:end-1) > sig(1:end-2) & sig(2:end-1) > sig(3:end);
locs = find(isMax);
if isempty(locs), return; end
if minH ~= 0
    locs = locs(sig(locs) >= minH);
    if isempty(locs), return; end
end
if minProm > 0
    prom = compute_prominence(sig, locs);
    locs = locs(prom >= minProm);
    if isempty(locs), return; end
end
minDistFr = max(1, round(minDistS * fps));
if minDistFr > 1 && length(locs) > 1
    [~, ord]  = sort(sig(locs),'descend');
    sorted    = locs(ord);
    keep      = true(length(sorted),1);
    for i = 1:length(sorted)
        if ~keep(i), continue; end
        tooClose = abs(sorted - sorted(i)) < minDistFr;
        tooClose(i) = false;
        keep(tooClose) = false;
    end
    locs = sort(sorted(keep));
end
if isempty(locs), return; end
if minWS > 0
    minWFr = max(1, round(minWS * fps));
    prom   = compute_prominence(sig, locs);
    keep   = true(length(locs),1);
    for i = 1:length(locs)
        pk       = locs(i);
        half_lvl = sig(pk) - prom(i)/2;
        lc = find(sig(1:pk) < half_lvl, 1,'last');
        if isempty(lc), lc = 1; end
        rc = find(sig(pk:end) < half_lvl, 1,'first');
        if isempty(rc), rc = n - pk + 1; end
        width   = (pk - lc) + (rc - 1);
        keep(i) = width >= minWFr;
    end
    locs = locs(keep);
end
locs = locs(:);
end

function prom = compute_prominence(sig, locs)
n = length(sig); np = length(locs);
prom = zeros(np,1);
for i = 1:np
    pk = locs(i);
    left_higher = locs(sig(locs) > sig(pk) & locs < pk);
    li = 1; if ~isempty(left_higher), li = left_higher(end); end
    right_higher = locs(sig(locs) > sig(pk) & locs > pk);
    ri = n; if ~isempty(right_higher), ri = right_higher(1); end
    base = max(min(sig(li:pk)), min(sig(pk:ri)));
    prom(i) = sig(pk) - base;
end
end

function sessions = scanForSessions(masterFolder)
%SCANFORSESSIONS  Find one session per recording folder.
%   Recursively finds folders containing a DLC CSV.  For each, looks for
%   the matching [csvStem]_breath_peak_data.mat first, then falls back to
%   the bare breath_peak_data.mat (legacy naming).  A folder without any
%   matching peak .mat is skipped.
sessions = {};
% accept all DLC CSV variants (best-10, best-20, snapshot_200, etc.)
hits = dir(fullfile(masterFolder,'**','*DLC*.csv'));
if isempty(hits), return; end
seen_folders = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(hits)
    h         = hits(i);
    if seen_folders.isKey(h.folder), continue; end   % one session per folder
    [~, stem] = fileparts(h.name);
    % try stem-prefixed first, then bare filename
    candidates = { ...
        fullfile(h.folder, [stem '_breath_peak_data.mat']), ...
        fullfile(h.folder, 'breath_peak_data.mat')};
    peakPath = '';
    for c = 1:numel(candidates)
        if isfile(candidates{c}), peakPath = candidates{c}; break; end
    end
    if isempty(peakPath), continue; end
    seen_folders(h.folder) = true;
    % name outputs to mirror whichever peak filename was found:
    %  - stem-prefixed peak -> stem-prefixed insp_start
    %  - bare peak           -> bare insp_start (breath_insp_start_data.mat)
    [~, pkBase] = fileparts(peakPath);
    matStem = erase(pkBase, '_breath_peak_data');
    if isempty(matStem)     % i.e. peakPath was just 'breath_peak_data.mat'
        matFile = fullfile(h.folder, 'breath_insp_start_data.mat');
        figBase = fullfile(h.folder, 'breath_insp_start_figure');
    else
        matFile = fullfile(h.folder, [matStem '_breath_insp_start_data.mat']);
        figBase = fullfile(h.folder, [matStem '_breath_insp_start_figure']);
    end
    sessions{end+1} = struct( ...                                    %#ok<AGROW>
        'folder',  h.folder, ...
        'csvStem', stem, ...
        'csvFile', fullfile(h.folder, h.name), ...
        'peakMat', peakPath, ...
        'matFile', matFile, ...
        'figBase', figBase);
end
if isempty(sessions), return; end
[~, ord] = sort(cellfun(@(s) s.peakMat, sessions, 'UniformOutput',false));
sessions = sessions(ord);
end
