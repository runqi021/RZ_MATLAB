function breathing_is_motion_gui()
%BREATHING_IS_MOTION_GUI  Interactive breathing peak detection from MC |shift|
%
%   "Breathing is motion": instead of DeepLabCut tracking of a behavior video,
%   this GUI derives the breathing trace directly from the 2-photon
%   motion-correction (NoRMCorre) shifts. Brain motion from breathing shows up
%   in the per-frame rigid shift, so |shift| (microns) is a built-in breathing
%   surrogate — no DLC pipeline needed for these recordings.
%
%   Mirrors breathing_peak_gui.m (same navigation, peak detection, manual
%   add/delete editor, and save format) but the front-end reads
%   *_preproc_MC_output.mat instead of a DLC CSV.
%
%   Usage:
%       breathing_is_motion_gui()
%
%   Output (saved next to each MC-output file, one per session):
%       [stem]_breath_peak_data.mat  - insp_onset_idx, insp_onsets_t,
%               insp_onsets_amp, insp_onsets_train, breath, t_breath,
%               findpeak_params   (same fields as breathing_peak_gui)
%
%   Requires: MATLAB R2021a+

%% -- Shared state (accessible by all nested callbacks) -------------------
S.shiftYX1   = [];      % [F x 2] first-pass MC shifts [y x] (pixels)
S.shiftYX2   = [];      % [F x 2] second-pass (residual) MC shifts, or []
S.pixelSize  = NaN;     % microns / pixel (from meta)
S.t          = [];      % [F x 1] seconds
S.breath     = [];      % [F x 1] processed (detrended, z-scored) signal
S.compTraces = [];      % [F x 3] component signals (shift_x, shift_y, |shift|) in um
S.peaks        = [];    % [K x 1] frame indices of accepted peaks
S.folderPath   = '';
S.mcFile       = '';
S.editMode     = 'none';  % 'none' | 'add' | 'delete'
S.masterFolder = '';
S.sessionList  = {};    % cell array of structs: .folder .mcFile .stem .matFile .figBase
S.sessionIdx   = 0;

% -- Factory defaults (restored by Reset button) ---------------------------
DEFAULTS.source       = 'MC1';
DEFAULTS.axis         = 'magnitude';
DEFAULTS.invert       = false;
DEFAULTS.minProm      = 1.00;
DEFAULTS.minDist      = 0.20;
DEFAULTS.minHeight     = 0.50;
DEFAULTS.minWidth      = 0.05;
DEFAULTS.snapWin       = 0.10;

%% -- Figure ---------------------------------------------------------------
fig = uifigure('Name','Breathing-is-Motion Peak Sorter (|shift|)', ...
               'Position',[40 30 1560 930], ...
               'WindowState','maximized');

rootGL = uigridlayout(fig, [1 2]);
rootGL.ColumnWidth   = {320, '1x'};
rootGL.Padding       = [6 6 6 6];
rootGL.ColumnSpacing = 8;

%% -- Left: control panel --------------------------------------------------
ctrlPan = uipanel(rootGL, 'Title','Controls', ...
    'FontSize',12,'FontWeight','bold');
ctrlPan.Layout.Column = 1;

nCtrlRows = 32;
cGL = uigridlayout(ctrlPan, [nCtrlRows 2]);
cGL.RowHeight   = repmat({'fit'}, 1, nCtrlRows);
cGL.ColumnWidth = {'fit','1x'};
cGL.Padding     = [8 8 8 8];
cGL.RowSpacing  = 4;

r = 0;   % row counter helper

% -- Session Navigation ----------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Session Navigation --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnLoad = uibutton(cGL,'Text','Browse Master Folder', ...
           'ButtonPushedFcn',@cb_loadMaster);
       btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];

r=r+1; lblMasterPath = uilabel(cGL,'Text','No folder selected', ...
           'WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]);
       lblMasterPath.Layout.Row=r; lblMasterPath.Layout.Column=[1 2];

r=r+1; lblCurrentFolder = uilabel(cGL,'Text','-', ...
           'FontWeight','bold','WordWrap','on','FontSize',10);
       lblCurrentFolder.Layout.Row=r; lblCurrentFolder.Layout.Column=[1 2];

r=r+1; lblNavStatus = uilabel(cGL,'Text','0 / 0  folders  |  0 complete', ...
           'FontSize',10,'FontColor',[0.25 0.65 0.25]);
       lblNavStatus.Layout.Row=r; lblNavStatus.Layout.Column=[1 2];

r=r+1; btnPrev = uibutton(cGL,'Text','< Prev','ButtonPushedFcn',@cb_prevSession);
       btnPrev.Layout.Row=r; btnPrev.Layout.Column=1;
       btnNext = uibutton(cGL,'Text','Next >','ButtonPushedFcn',@cb_nextSession);
       btnNext.Layout.Row=r; btnNext.Layout.Column=2;

% -- Imaging ----------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Imaging --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Imaging FPS:');
       h.Layout.Row=r; h.Layout.Column=1;
       efFps = uieditfield(cGL,'numeric','Value',30, ...
           'Limits',[1 2000],'ValueChangedFcn',@cb_updateTrace);
       efFps.Layout.Row=r; efFps.Layout.Column=2;

r=r+1; lblPx = uilabel(cGL,'Text','px size: -', ...
           'FontSize',9,'FontColor',[0.5 0.5 0.5]);
       lblPx.Layout.Row=r; lblPx.Layout.Column=[1 2];

% -- Shift Source -----------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Shift Source --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Source:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddSource = uidropdown(cGL,'Items',{'MC1','MC1+MC2 (sum)'}, ...
           'ItemsData',{'MC1','MC1+MC2'},'Value','MC1', ...
           'ValueChangedFcn',@cb_updateTrace);
       ddSource.Layout.Row=r; ddSource.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Axis:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddAxis = uidropdown(cGL,'Items',{'magnitude','x','y','-x','-y'}, ...
           'Value','magnitude','ValueChangedFcn',@cb_updateTrace);
       ddAxis.Layout.Row=r; ddAxis.Layout.Column=2;

r=r+1; lblSrcNote = uilabel(cGL, ...
           'Text','magnitude rectifies motion (freq-doubles); x/y keep sign', ...
           'FontSize',9,'FontColor',[0.7 0.5 0.1],'WordWrap','on');
       lblSrcNote.Layout.Row=r; lblSrcNote.Layout.Column=[1 2];

r=r+1; cbInvert = uicheckbox(cGL,'Text','Invert signal','Value',false, ...
           'ValueChangedFcn',@cb_updateTrace);
       cbInvert.Layout.Row=r; cbInvert.Layout.Column=[1 2];

r=r+1; btnTrace = uibutton(cGL,'Text','Update Trace', ...
           'ButtonPushedFcn',@cb_updateTrace);
       btnTrace.Layout.Row=r; btnTrace.Layout.Column=[1 2];

% -- Peak Detection ---------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Peak Detection --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

paramLabels = {'Min Prominence:','Min Distance (s):','Min Height:','Min Width (s):'};
paramDefs   = {1.00, 0.20, 0.50, 0.05};
efParam     = gobjects(4,1);
for p = 1:4
    r=r+1;
    h=uilabel(cGL,'Text',paramLabels{p});
    h.Layout.Row=r; h.Layout.Column=1;
    efParam(p) = uieditfield(cGL,'numeric','Value',paramDefs{p},'Limits',[-Inf Inf], ...
        'ValueChangedFcn',@(~,~) saveWorkingParams());
    efParam(p).Layout.Row=r; efParam(p).Layout.Column=2;
end

r=r+1; h=uilabel(cGL,'Text','Snap window (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efSnap = uieditfield(cGL,'numeric','Value',0.10,'Limits',[0 10], ...
           'ValueChangedFcn',@(~,~) saveWorkingParams());
       efSnap.Layout.Row=r; efSnap.Layout.Column=2;

r=r+1; btnDetect = uibutton(cGL,'Text','Detect Peaks', ...
           'ButtonPushedFcn',@cb_detect, ...
           'BackgroundColor',[0.18 0.65 0.30],'FontColor','white');
       btnDetect.Layout.Row=r; btnDetect.Layout.Column=[1 2];

% -- Manual Editing ---------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Manual Editing --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; lblMode = uilabel(cGL,'Text','Mode: Navigate', ...
           'FontWeight','bold','FontSize',11, ...
           'FontColor',[0.90 0.80 0.10]);
       lblMode.Layout.Row=r; lblMode.Layout.Column=[1 2];

r=r+1; btnAdd = uibutton(cGL,'Text','+ Add Mode','ButtonPushedFcn',@cb_setAdd, ...
           'BackgroundColor',[0.20 0.42 0.85],'FontColor','white');
       btnAdd.Layout.Row=r; btnAdd.Layout.Column=1;
       btnDel = uibutton(cGL,'Text','- Delete Mode','ButtonPushedFcn',@cb_setDelete, ...
           'BackgroundColor',[0.82 0.22 0.18],'FontColor','white');
       btnDel.Layout.Row=r; btnDel.Layout.Column=2;

r=r+1; btnNav = uibutton(cGL,'Text','Navigate (N)', ...
           'ButtonPushedFcn',@cb_setNav);
       btnNav.Layout.Row=r; btnNav.Layout.Column=1;
       btnClear = uibutton(cGL,'Text','Clear All Peaks', ...
           'ButtonPushedFcn',@cb_clearPeaks, ...
           'BackgroundColor',[0.40 0.10 0.10],'FontColor','white');
       btnClear.Layout.Row=r; btnClear.Layout.Column=2;

r=r+1; lblCount = uilabel(cGL,'Text','Peaks: 0', ...
           'FontWeight','bold','FontSize',12);
       lblCount.Layout.Row=r; lblCount.Layout.Column=[1 2];

% -- Save -------------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Save --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnSave = uibutton(cGL,'Text','Save breath_peak_data.mat', ...
           'ButtonPushedFcn',@cb_save, ...
           'BackgroundColor',[0.65 0.45 0.05],'FontColor','white');
       btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];

r=r+1; btnReset = uibutton(cGL,'Text','Reset to Defaults', ...
           'ButtonPushedFcn',@cb_reset, ...
           'BackgroundColor',[0.30 0.30 0.30],'FontColor','white');
       btnReset.Layout.Row=r; btnReset.Layout.Column=[1 2];

% -- Load persisted working params (overrides hardcoded defaults) ----------
if ispref('BreathingIsMotionGUI','workingParams')
    try, applyWorkingParams(getpref('BreathingIsMotionGUI','workingParams')); catch, end
end

%% -- Right: plot panel ----------------------------------------------------
rightGL = uigridlayout(rootGL, [2 1]);
rightGL.Layout.Column = 2;
rightGL.RowHeight     = {'1x','1x'};
rightGL.Padding       = [0 0 0 0];
rightGL.RowSpacing    = 6;

% Main axes: processed trace + peaks
axMain = uiaxes(rightGL);
axMain.Layout.Row           = 1;
axMain.Title.String         = 'Breathing Trace (from |shift|)  |  Left-click: add/delete peaks  |  N: navigate mode';
axMain.Title.Interpreter    = 'none';
axMain.XLabel.String        = 'Time (s)';
axMain.XLabel.Interpreter   = 'none';
axMain.YLabel.String        = 'Signal (z-score)';
axMain.YLabel.Interpreter   = 'none';
hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');

% Component traces axes (shift_x, shift_y, |shift| in um)
axComp = uiaxes(rightGL);
axComp.Layout.Row           = 2;
axComp.Title.String         = 'Shift components (um): shift_x, shift_y, |shift| (offset)';
axComp.Title.Interpreter    = 'none';
axComp.XLabel.String        = 'Time (s)';
axComp.XLabel.Interpreter   = 'none';
axComp.YLabel.String        = '';
hold(axComp,'on'); grid(axComp,'on'); box(axComp,'on');

% Link x-axes
linkaxes([axMain axComp],'x');

% Persistent peak overlay handles (updated in-place - no object churn)
hPeakSc = gobjects(0);   % scatter on axMain
hPeakVl = gobjects(0);   % NaN-separated line on axComp

% Keyboard shortcut on figure
fig.KeyPressFcn = @cb_key;

%% =======================================================================
%% CALLBACKS
%% =======================================================================

    % -- Browse master folder and build session list ----------------------
    function cb_loadMaster(~,~)
        master = uigetdir(S.masterFolder,'Select Master Folder');
        if isequal(master,0), return; end
        S.masterFolder = master;

        sessions = scanForSessions(master);
        if isempty(sessions)
            uialert(fig,'No *_preproc_MC_output.mat files found in any subfolder.','No Sessions');
            return;
        end
        S.sessionList = sessions;
        S.sessionIdx  = 1;

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

    % -- Load the session at S.sessionIdx ----------------------------------
    function loadCurrentSession()
        sess         = S.sessionList{S.sessionIdx};
        S.folderPath = sess.folder;
        S.mcFile     = sess.mcFile;
        S.peaks      = [];
        S.shiftYX1   = [];
        S.shiftYX2   = [];

        lblCurrentFolder.Text = sess.stem;
        updateNavDisplay();

        % Load first-pass MC shifts -> [F x 2] pixels [y x]
        try
            S.shiftYX1 = load_mc_shifts(sess.mcFile);
        catch ME
            uialert(fig, ME.message,'MC Load Error'); return;
        end

        % Optional second-pass (residual) shifts for the sum option
        if ~isempty(sess.mcFile2) && isfile(sess.mcFile2)
            try
                S.shiftYX2 = load_mc_shifts(sess.mcFile2);
                if size(S.shiftYX2,1) ~= size(S.shiftYX1,1)
                    S.shiftYX2 = [];   % length mismatch -> ignore
                end
            catch
                S.shiftYX2 = [];
            end
        end
        if isempty(S.shiftYX2)
            ddSource.Enable = 'off';   ddSource.Value = 'MC1';
        else
            ddSource.Enable = 'on';
        end

        % Sibling meta for fps + pixel size
        S.pixelSize = NaN;
        if ~isempty(sess.metaFile) && isfile(sess.metaFile)
            try
                M = load(sess.metaFile);
                if isfield(M,'fps') && isfinite(M.fps),         efFps.Value = M.fps; end
                if isfield(M,'pixelSize_um'),                   S.pixelSize = M.pixelSize_um; end
            catch
            end
        end
        if isfinite(S.pixelSize)
            lblPx.Text = sprintf('px size: %.4f um', S.pixelSize);
        else
            lblPx.Text = 'px size: unknown (using 1 px units)';
        end

        % Restore previously saved peaks/params if present
        matFile = sess.matFile;
        if isfile(matFile)
            try
                saved = load(matFile,'insp_onset_idx','findpeak_params');
                if isfield(saved,'findpeak_params')
                    fp = saved.findpeak_params;
                    if isfield(fp,'fps') && isfinite(fp.fps), efFps.Value = fp.fps; end
                    if isfield(fp,'minProminence'), efParam(1).Value = fp.minProminence; end
                    if isfield(fp,'minDist_s'),     efParam(2).Value = fp.minDist_s;     end
                    if isfield(fp,'minHeight'),     efParam(3).Value = fp.minHeight;     end
                    if isfield(fp,'minWidth_s'),    efParam(4).Value = fp.minWidth_s;    end
                    if isfield(fp,'snap_win_s'),    efSnap.Value     = fp.snap_win_s;    end
                    if isfield(fp,'axis'),          ddAxis.Value     = fp.axis;          end
                    if isfield(fp,'shift_source') && strcmp(ddSource.Enable,'on')
                        ddSource.Value = fp.shift_source;
                    end
                    if isfield(fp,'inverted'),      cbInvert.Value   = fp.inverted;      end
                end
                cb_updateTrace();
                if isfield(saved,'insp_onset_idx') && ~isempty(saved.insp_onset_idx)
                    S.peaks = saved.insp_onset_idx(:);
                    redraw_peaks();
                end
                lblCurrentFolder.Text = [sess.stem '  (done)'];
                updateNavDisplay();
                return;
            catch
                % Corrupted mat - fall through to fresh trace
            end
        end

        cb_updateTrace();
    end

    function updateNavDisplay()
        n = numel(S.sessionList);
        nDone = sum(cellfun(@(s) isfile(s.matFile), S.sessionList));
        lblNavStatus.Text = sprintf('%d / %d  folders  |  %d complete', ...
                                    S.sessionIdx, n, nDone);
    end

    % -- Recompute breathing trace from shifts -----------------------------
    function cb_updateTrace(~,~)
        if isempty(S.shiftYX1), return; end

        fps  = efFps.Value;
        nF   = size(S.shiftYX1,1);
        S.t  = (0:nF-1)' / fps;

        % Select source: MC1 or MC1+MC2 (cumulative displacement)
        if strcmp(ddSource.Value,'MC1+MC2') && ~isempty(S.shiftYX2)
            yx = S.shiftYX1 + S.shiftYX2;
        else
            yx = S.shiftYX1;
        end

        % Pixels -> microns (fall back to pixel units if px size unknown)
        sc = S.pixelSize;  if ~isfinite(sc), sc = 1; end
        ys = yx(:,1) * sc;     % y in um
        xs = yx(:,2) * sc;     % x in um
        mg = sqrt(xs.^2 + ys.^2);

        % Component display traces: shift_x, shift_y, |shift| (um)
        S.compTraces = [xs, ys, mg];

        % Build the 1-D breathing signal from the selected axis.
        % NOTE: 'magnitude' rectifies the back-and-forth breathing motion
        % (frequency-doubles, smears phase). Signed x/y preserve the rhythm.
        switch ddAxis.Value
            case 'magnitude', sig =  mg;
            case 'x',         sig =  xs;
            case 'y',         sig =  ys;
            case '-x',        sig = -xs;
            case '-y',        sig = -ys;
        end

        % Detrend + z-score (same as breathing_peak_gui)
        sig = detrend(sig);
        mu  = mean(sig,'omitnan');
        sd  = std(sig,0,'omitnan');
        if sd > 0, sig = (sig - mu) / sd; end
        if cbInvert.Value, sig = -sig; end
        S.breath = sig;

        saveWorkingParams();
        refresh_compPlot();
        refresh_mainPlot();
    end

    % -- Detect peaks ------------------------------------------------------
    function cb_detect(~,~)
        if isempty(S.breath)
            uialert(fig,'Load and process data first.','No Data'); return;
        end
        % Custom finder avoids Chronux findpeaks shadowing MATLAB's version
        S.peaks = find_peaks_local(S.breath, efFps.Value, ...
            efParam(1).Value, efParam(2).Value, ...
            efParam(3).Value, efParam(4).Value);
        saveWorkingParams();
        redraw_peaks();
        updatePeakCount();
    end

    % -- Axes click: add or delete peak ------------------------------------
    function cb_axClick(~,~)
        if isempty(S.breath), return; end

        cp     = axMain.CurrentPoint;
        clickT = cp(1,1);
        if clickT < S.t(1) || clickT > S.t(end), return; end

        [~, nearIdx] = min(abs(S.t - clickT));

        switch S.editMode
            case 'add'
                snapFr = max(1, round(efSnap.Value * efFps.Value));
                lo = max(1, nearIdx - snapFr);
                hi = min(length(S.breath), nearIdx + snapFr);
                [~, imax] = max(S.breath(lo:hi));
                peakIdx = lo + imax - 1;
                if isempty(S.peaks) || min(abs(S.peaks - peakIdx)) > 2
                    S.peaks = sort([S.peaks; peakIdx]);
                end
            case 'delete'
                if isempty(S.peaks), return; end
                [~, ci] = min(abs(S.peaks - nearIdx));
                S.peaks(ci) = [];
        end

        redraw_peaks();
        updatePeakCount();
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
        lblMode.Text      = 'Mode: + ADD  (click -> snap to local max)';
        lblMode.FontColor = [0.30 0.62 1.00];
        disableDefaultInteractivity(axMain);
        axMain.ButtonDownFcn = @cb_axClick;
    end

    function cb_setDelete(~,~)
        S.editMode = 'delete';
        lblMode.Text      = 'Mode: - DELETE  (click -> remove nearest)';
        lblMode.FontColor = [1.00 0.35 0.25];
        disableDefaultInteractivity(axMain);
        axMain.ButtonDownFcn = @cb_axClick;
    end

    function cb_setNav(~,~)
        S.editMode = 'none';
        lblMode.Text      = 'Mode: Navigate';
        lblMode.FontColor = [0.90 0.80 0.10];
        axMain.ButtonDownFcn = [];
        enableDefaultInteractivity(axMain);
    end

    function cb_reset(~,~)
        applyWorkingParams(DEFAULTS);
        saveWorkingParams();
        cb_updateTrace();
    end

    function applyWorkingParams(p)
        if isfield(p,'source') && strcmp(ddSource.Enable,'on'), ddSource.Value = p.source; end
        if isfield(p,'axis'),   ddAxis.Value   = p.axis;   end
        if isfield(p,'invert'), cbInvert.Value = p.invert; end
        efParam(1).Value = p.minProm;
        efParam(2).Value = p.minDist;
        efParam(3).Value = p.minHeight;
        efParam(4).Value = p.minWidth;
        efSnap.Value     = p.snapWin;
    end

    function saveWorkingParams()
        wp.source    = ddSource.Value;
        wp.axis      = ddAxis.Value;
        wp.invert    = cbInvert.Value;
        wp.minProm   = efParam(1).Value;
        wp.minDist   = efParam(2).Value;
        wp.minHeight = efParam(3).Value;
        wp.minWidth  = efParam(4).Value;
        wp.snapWin   = efSnap.Value;
        setpref('BreathingIsMotionGUI','workingParams',wp);
    end

    function cb_clearPeaks(~,~)
        S.peaks = [];
        redraw_peaks();
        updatePeakCount();
    end

    % -- Save --------------------------------------------------------------
    function cb_save(~,~)
        if isempty(S.folderPath)
            uialert(fig,'No session loaded.','Error'); return;
        end

        insp_onset_idx    = S.peaks;                 % [K x 1] frame indices
        insp_onsets_t     = [];                      % [K x 1] seconds
        insp_onsets_amp   = [];                      % [K x 1] z-score at peak
        insp_onsets_train = zeros(size(S.breath));   % [F x 1] binary train
        if ~isempty(S.peaks) && ~isempty(S.t)
            insp_onsets_t   = S.t(S.peaks);
            insp_onsets_amp = S.breath(S.peaks);
            insp_onsets_train(S.peaks) = 1;
        end
        breath   = S.breath;       %#ok<NASGU>
        t_breath = S.t;            %#ok<NASGU>

        findpeak_params.minProminence  = efParam(1).Value;
        findpeak_params.minDist_s      = efParam(2).Value;
        findpeak_params.minHeight      = efParam(3).Value;
        findpeak_params.minWidth_s     = efParam(4).Value;
        findpeak_params.snap_win_s     = efSnap.Value;
        findpeak_params.fps            = efFps.Value;
        findpeak_params.shift_source   = ddSource.Value;     % 'MC1' | 'MC1+MC2'
        findpeak_params.axis           = ddAxis.Value;       % 'magnitude'|'x'|'y'|...
        findpeak_params.inverted       = cbInvert.Value;
        findpeak_params.pixelSize_um   = S.pixelSize;
        findpeak_params.mc_file        = S.mcFile;
        findpeak_params.saved_at       = datestr(now);       %#ok<TNOW1,DATST>

        sess    = S.sessionList{S.sessionIdx};
        outFile = sess.matFile;
        save(outFile,'insp_onset_idx','insp_onsets_t','insp_onsets_amp','insp_onsets_train', ...
                     'breath','t_breath','findpeak_params');

        saveScreenshots(sess.figBase);

        uialert(fig, sprintf('Saved %d peaks ->\n%s\n+ .fig and .png', numel(S.peaks), outFile), ...
            'Saved','Icon','success');
        updateNavDisplay();
    end

    function saveScreenshots(outBase)
        sess = S.sessionList{S.sessionIdx};
        sessName = sess.stem;

        hFig = figure('Visible','off','Color','white','Position',[100 100 1400 820]);

        ax1 = subplot(2,1,1,'Parent',hFig);
        hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
        copyobj(axMain.Children, ax1);
        ax1.XLim          = axMain.XLim;
        ax1.YLim          = axMain.YLim;
        ax1.XLabel.String       = axMain.XLabel.String;
        ax1.YLabel.String       = axMain.YLabel.String;
        ax1.Title.String        = [strrep(sessName,'_',' ') '  -  Breathing Trace (|shift|)'];
        ax1.Title.Interpreter   = 'none';
        ax1.XLabel.Interpreter  = 'none';
        ax1.YLabel.Interpreter  = 'none';

        ax2 = subplot(2,1,2,'Parent',hFig);
        hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
        copyobj(axComp.Children, ax2);
        ax2.XLim                = axComp.XLim;
        ax2.YLim                = axComp.YLim;
        ax2.XLabel.String       = axComp.XLabel.String;
        ax2.YLabel.String       = '';
        ax2.YTick               = [];
        ax2.Title.String        = [strrep(sessName,'_',' ') '  -  Shift components (um)'];
        ax2.Title.Interpreter   = 'none';
        ax2.XLabel.Interpreter  = 'none';

        savefig(hFig,  [outBase '.fig']);
        exportgraphics(hFig, [outBase '.png'], 'Resolution',150);
        close(hFig);
    end

%% =======================================================================
%% PLOT HELPERS
%% =======================================================================

    function refresh_compPlot()
        cla(axComp); hold(axComp,'on');
        axComp.YLimMode = 'auto';
        if isempty(S.t) || isempty(S.compTraces), return; end

        names  = {'shift_x','shift_y','|shift|'};
        cols   = {[0.2 0.2 0.8],[0.8 0.2 0.2],[0 0 0]};
        offset = 0;
        for d = size(S.compTraces,2):-1:1
            sig = S.compTraces(:,d);
            if all(isnan(sig)), continue; end
            rng = max(sig) - min(sig);
            if rng < eps, rng = 1; end
            sig_n = (sig - min(sig)) / rng + offset;
            h = plot(axComp, S.t, sig_n, 'Color', cols{d}, 'LineWidth',0.85);
            h.HitTest = 'off';
            text(axComp, S.t(1), offset+0.5, [' ' names{d}], ...
                'Color', cols{d}, 'FontSize',8, 'VerticalAlignment','middle', ...
                'HitTest','off');
            offset = offset + 1.3;
        end
        legend(axComp,'off');
        axComp.YTick = [];
        if ~isempty(S.t), xlim(axComp,[S.t(1) S.t(end)]); end
        if offset > 0, axComp.YLim = [-0.2, offset]; end
        axComp.YLimMode = 'manual';

        hPeakVl = plot(axComp, NaN, NaN, '-r', 'LineWidth',0.6, 'Tag','peak_vline');
        hPeakVl.HitTest = 'off';
    end

    function refresh_mainPlot()
        cla(axMain); hold(axMain,'on');
        if isempty(S.t) || isempty(S.breath), return; end

        h = plot(axMain, S.t, S.breath, 'Color',[0 0 0], 'LineWidth',0.9);
        h.HitTest = 'off';

        h2 = yline(axMain, 0, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.6);
        h2.HitTest = 'off';

        hPeakSc = scatter(axMain, [], [], 75, 'v', 'filled', ...
            'MarkerFaceColor',[1.00 0.28 0.28], 'MarkerEdgeColor','none', ...
            'Tag','peak_marker');
        hPeakSc.HitTest = 'off';

        xlim(axMain,[S.t(1) S.t(end)]);
        redraw_peaks();
        updatePeakCount();
    end

    function redraw_peaks()
        valid = [];
        if ~isempty(S.peaks) && ~isempty(S.breath)
            valid = S.peaks(S.peaks >= 1 & S.peaks <= length(S.breath));
        end

        if isgraphics(hPeakSc)
            if isempty(valid)
                set(hPeakSc, 'XData',[], 'YData',[]);
            else
                set(hPeakSc, 'XData',S.t(valid), 'YData',S.breath(valid));
            end
        end

        if isgraphics(hPeakVl)
            if isempty(valid) || isempty(S.t)
                set(hPeakVl, 'XData',NaN, 'YData',NaN);
            else
                tpk = S.t(valid);   n = length(tpk);
                yl  = axComp.YLim;
                xv  = reshape([tpk'; tpk'; NaN(1,n)], 1, []);
                yv  = repmat([yl(1) yl(2) NaN], 1, n);
                set(hPeakVl, 'XData',xv, 'YData',yv);
            end
        end

        updatePeakCount();
    end

    function updatePeakCount()
        lblCount.Text = sprintf('Peaks: %d', numel(S.peaks));
    end

end  % breathing_is_motion_gui


%% =======================================================================
%% LOCAL HELPER - load NoRMCorre rigid shifts as [F x 2] = [y x] pixels
%% =======================================================================
function yx = load_mc_shifts(mcFile)
L = load(mcFile);
if ~isfield(L,'mcOut') || ~isfield(L.mcOut,'shifts')
    error('No mcOut.shifts in %s', mcFile);
end
sh = L.mcOut.shifts;
nF = numel(sh);
yx = zeros(nF, 2);
for k = 1:nF
    s = squeeze(sh(k).shifts);   % rigid [y x] in pixels
    yx(k,:) = s(:)';
end
end  % load_mc_shifts


%% =======================================================================
%% LOCAL HELPER - custom peak finder
%%   Avoids Chronux findpeaks shadowing the Signal Processing Toolbox
%% =======================================================================
function locs = find_peaks_local(sig, fps, minProm, minDistS, minH, minWS)
%FIND_PEAKS_LOCAL  Local-max finder with prominence / distance / height / width filters.

sig  = sig(:);
n    = length(sig);
locs = [];

% -- 1. Local maxima (strict: larger than both neighbours) ----------------
isMax = false(n,1);
isMax(2:end-1) = sig(2:end-1) > sig(1:end-2) & sig(2:end-1) > sig(3:end);
locs = find(isMax);
if isempty(locs), return; end

% -- 2. Height threshold ----------------------------------------------------
if minH ~= 0
    locs = locs(sig(locs) >= minH);
    if isempty(locs), return; end
end

% -- 3. Prominence filter ---------------------------------------------------
if minProm > 0
    prom = compute_prominence(sig, locs);
    locs = locs(prom >= minProm);
    if isempty(locs), return; end
end

% -- 4. Min distance: greedily keep highest peak, remove neighbours --------
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

% -- 5. Min width at half-prominence ---------------------------------------
if minWS > 0
    minWFr = max(1, round(minWS * fps));
    prom   = compute_prominence(sig, locs);
    keep   = true(length(locs),1);
    for i = 1:length(locs)
        pk       = locs(i);
        half_lvl = sig(pk) - prom(i)/2;
        left_seg  = sig(1:pk);
        lc = find(left_seg < half_lvl, 1,'last');
        if isempty(lc), lc = 1; end
        right_seg = sig(pk:end);
        rc = find(right_seg < half_lvl, 1,'first');
        if isempty(rc), rc = n - pk + 1; end
        width   = (pk - lc) + (rc - 1);
        keep(i) = width >= minWFr;
    end
    locs = locs(keep);
end

locs = locs(:);
end  % find_peaks_local


function prom = compute_prominence(sig, locs)
%COMPUTE_PROMINENCE  Prominence of each peak in locs relative to sig.
n    = length(sig);
np   = length(locs);
prom = zeros(np,1);
for i = 1:np
    pk = locs(i);
    left_higher = locs(sig(locs) > sig(pk) & locs < pk);
    li = 1;
    if ~isempty(left_higher), li = left_higher(end); end
    right_higher = locs(sig(locs) > sig(pk) & locs > pk);
    ri = n;
    if ~isempty(right_higher), ri = right_higher(1); end
    left_min  = min(sig(li:pk));
    right_min = min(sig(pk:ri));
    base      = max(left_min, right_min);
    prom(i)   = sig(pk) - base;
end
end  % compute_prominence


function sessions = scanForSessions(masterFolder)
%SCANFORSESSIONS  Recursively find all first-pass MC outputs under masterFolder.
%   One session per *_preproc_MC_output.mat (EXCLUDING *_preproc_MC_MC_output.mat).
%   Returns a cell array of structs with fields:
%     .folder    - directory containing the MC output
%     .mcFile    - full path to the first-pass MC output
%     .mcFile2   - full path to the second-pass (_MC_MC_output) if present, else ''
%     .metaFile  - full path to sibling *_meta.mat if present, else ''
%     .stem      - MC filename without extension (used for output naming)
%     .matFile   - output [stem]_breath_peak_data.mat
%     .figBase   - output figure base (no extension)

sessions = {};

hits = dir(fullfile(masterFolder,'**','*_preproc_MC_output.mat'));
% Exclude the second-pass residual outputs (*_preproc_MC_MC_output.mat)
hits = hits(~contains({hits.name}, '_MC_MC_output'));
if isempty(hits), return; end

for i = 1:numel(hits)
    h       = hits(i);
    mcPath  = fullfile(h.folder, h.name);
    [~, stem] = fileparts(h.name);

    % second-pass sibling: replace _preproc_MC_output -> _preproc_MC_MC_output
    mc2 = strrep(mcPath, '_preproc_MC_output.mat', '_preproc_MC_MC_output.mat');
    if ~isfile(mc2), mc2 = ''; end

    % meta sibling: prefer the matching *_meta.mat in the same folder
    metaFile = '';
    mhits = dir(fullfile(h.folder, '*_meta.mat'));
    if ~isempty(mhits)
        metaFile = fullfile(mhits(1).folder, mhits(1).name);
    end

    sessions{end+1} = struct( ...                            %#ok<AGROW>
        'folder',  h.folder, ...
        'mcFile',  mcPath, ...
        'mcFile2', mc2, ...
        'metaFile',metaFile, ...
        'stem',    stem, ...
        'matFile', fullfile(h.folder, [stem '_breath_peak_data.mat']), ...
        'figBase', fullfile(h.folder, [stem '_breath_peak_figure']));
end

[~, ord] = sort(cellfun(@(s) s.mcFile, sessions, 'UniformOutput',false));
sessions = sessions(ord);
end  % scanForSessions
