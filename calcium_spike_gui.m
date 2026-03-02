function calcium_spike_gui()
%CALCIUM_SPIKE_GUI  Interactive GUI for calcium spike detection with per-ROI inspection
%
%   Scans a master folder for *_dFF.mat files (output of helper.dFF_RZ),
%   pairs them with *_cpSAM_output.mat (ROI masks) and *_minusDark_MC.tif
%   (motion-corrected movie), and lets you detect / manually edit calcium
%   spikes on each ROI individually.
%
%   Usage:
%       calcium_spike_gui()
%
%   Output (saved to each session folder):
%       ca_spike_data.mat  -- roi_spikes struct array, spike_params
%
%   Master output (saved to master folder):
%       ca_spike_master.mat -- master struct indexing all sessions
%
%   Requires: MATLAB R2021a+ with Image Processing Toolbox

%% -- Shared state (accessible by all nested callbacks) --------------------
S.dFF         = [];      % [T x N] dF/F matrix for current session
S.t_img       = [];      % [T x 1] time vector (seconds)
S.fps         = 30;      % imaging frame rate
S.nROI        = 0;       % number of ROIs in current session
S.curROI      = 1;       % currently displayed ROI index
S.spikes      = {};      % {N x 1} cell array of spike frame indices per ROI
S.maskL       = [];      % [H x W] label mask from cpSAM
S.meanImg     = [];      % [H x W] mean projection for video crop
S.hasMask     = false;
S.hasTif      = false;
S.tifFile     = '';
S.samFile     = '';
S.dffFile     = '';
S.folderPath  = '';
S.editMode    = 'none';  % 'none' | 'add' | 'delete'
S.masterFolder = '';
S.sessionList  = {};     % cell array of structs: .folder .dffFile
S.sessionIdx   = 0;
S.dFFout       = [];     % full dFFout struct (may contain spikes_oasis)

% -- Factory defaults (restored by Reset button) --------------------------
DEFAULTS.fps        = 30;
DEFAULTS.method     = 'Threshold';
DEFAULTS.threshold  = 2.0;
DEFAULTS.minDist    = 0.5;
DEFAULTS.minWidth   = 0.1;
DEFAULTS.snapWin    = 0.20;
DEFAULTS.diameter   = 30;
DEFAULTS.zoomFactor = 4;

%% -- Figure ---------------------------------------------------------------
fig = uifigure('Name','Calcium Spike Sorter', ...
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

r=r+1; lblCurrentFolder = uilabel(cGL,'Text','---', ...
           'FontWeight','bold','WordWrap','on','FontSize',10);
       lblCurrentFolder.Layout.Row=r; lblCurrentFolder.Layout.Column=[1 2];

r=r+1; lblNavStatus = uilabel(cGL,'Text','0 / 0  sessions  |  0 complete', ...
           'FontSize',10,'FontColor',[0.25 0.65 0.25]);
       lblNavStatus.Layout.Row=r; lblNavStatus.Layout.Column=[1 2];

r=r+1; btnPrev = uibutton(cGL,'Text','<- Prev Session','ButtonPushedFcn',@cb_prevSession);
       btnPrev.Layout.Row=r; btnPrev.Layout.Column=1;
       btnNext = uibutton(cGL,'Text','Next Session ->','ButtonPushedFcn',@cb_nextSession);
       btnNext.Layout.Row=r; btnNext.Layout.Column=2;

% -- ROI Navigation --------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- ROI Navigation --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; lblROI = uilabel(cGL,'Text','ROI: 0 / 0','FontWeight','bold','FontSize',12);
       lblROI.Layout.Row=r; lblROI.Layout.Column=[1 2];

r=r+1; btnPrevROI = uibutton(cGL,'Text','<- Prev ROI','ButtonPushedFcn',@cb_prevROI);
       btnPrevROI.Layout.Row=r; btnPrevROI.Layout.Column=1;
       btnNextROI = uibutton(cGL,'Text','Next ROI ->','ButtonPushedFcn',@cb_nextROI);
       btnNextROI.Layout.Row=r; btnNextROI.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Jump to ROI:');
       h.Layout.Row=r; h.Layout.Column=1;
       spnROI = uispinner(cGL,'Value',1,'Limits',[1 1],'Step',1, ...
           'ValueChangedFcn',@cb_jumpROI);
       spnROI.Layout.Row=r; spnROI.Layout.Column=2;

% -- Imaging Parameters ---------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Parameters --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Imaging FPS:');
       h.Layout.Row=r; h.Layout.Column=1;
       efFps = uieditfield(cGL,'numeric','Value',DEFAULTS.fps, ...
           'Limits',[1 2000],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efFps.Layout.Row=r; efFps.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Diameter (px):');
       h.Layout.Row=r; h.Layout.Column=1;
       efDiameter = uieditfield(cGL,'numeric','Value',DEFAULTS.diameter, ...
           'Limits',[1 500],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efDiameter.Layout.Row=r; efDiameter.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Zoom Factor:');
       h.Layout.Row=r; h.Layout.Column=1;
       efZoom = uieditfield(cGL,'numeric','Value',DEFAULTS.zoomFactor, ...
           'Limits',[1 20],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efZoom.Layout.Row=r; efZoom.Layout.Column=2;

% -- Spike Detection -------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Spike Detection --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Method:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddMethod = uidropdown(cGL,'Items',{'Threshold','OASIS'},'Value','Threshold', ...
           'ValueChangedFcn',@(~,~) saveWorkingParams());
       ddMethod.Layout.Row=r; ddMethod.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Threshold (dF/F):');
       h.Layout.Row=r; h.Layout.Column=1;
       efThresh = uieditfield(cGL,'numeric','Value',DEFAULTS.threshold, ...
           'Limits',[-Inf Inf],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efThresh.Layout.Row=r; efThresh.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Min Distance (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efMinDist = uieditfield(cGL,'numeric','Value',DEFAULTS.minDist, ...
           'Limits',[0 Inf],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efMinDist.Layout.Row=r; efMinDist.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Min Width (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efMinWidth = uieditfield(cGL,'numeric','Value',DEFAULTS.minWidth, ...
           'Limits',[0 Inf],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efMinWidth.Layout.Row=r; efMinWidth.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Snap window (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efSnap = uieditfield(cGL,'numeric','Value',DEFAULTS.snapWin, ...
           'Limits',[0 10],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efSnap.Layout.Row=r; efSnap.Layout.Column=2;

r=r+1; btnDetect = uibutton(cGL,'Text','Detect Spikes (Current ROI)', ...
           'ButtonPushedFcn',@cb_detect, ...
           'BackgroundColor',[0.18 0.65 0.30],'FontColor','white');
       btnDetect.Layout.Row=r; btnDetect.Layout.Column=[1 2];

r=r+1; btnDetectAll = uibutton(cGL,'Text','Detect Spikes (All ROIs)', ...
           'ButtonPushedFcn',@cb_detectAll, ...
           'BackgroundColor',[0.10 0.50 0.22],'FontColor','white');
       btnDetectAll.Layout.Row=r; btnDetectAll.Layout.Column=[1 2];

% -- Manual Editing --------------------------------------------------------
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
       btnClearSpikes = uibutton(cGL,'Text','Clear ROI Spikes', ...
           'ButtonPushedFcn',@cb_clearSpikes, ...
           'BackgroundColor',[0.40 0.10 0.10],'FontColor','white');
       btnClearSpikes.Layout.Row=r; btnClearSpikes.Layout.Column=2;

r=r+1; lblCount = uilabel(cGL,'Text','Spikes: 0', ...
           'FontWeight','bold','FontSize',12);
       lblCount.Layout.Row=r; lblCount.Layout.Column=[1 2];

% -- Save ------------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Save --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnSave = uibutton(cGL,'Text','Save ca_spike_data.mat', ...
           'ButtonPushedFcn',@cb_save, ...
           'BackgroundColor',[0.65 0.45 0.05],'FontColor','white');
       btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];

r=r+1; btnReset = uibutton(cGL,'Text','Reset to Defaults', ...
           'ButtonPushedFcn',@cb_reset, ...
           'BackgroundColor',[0.30 0.30 0.30],'FontColor','white');
       btnReset.Layout.Row=r; btnReset.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','-- Clear Data --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnClearAll = uibutton(cGL,'Text','Clear All Sorted Data (All Folders)', ...
           'ButtonPushedFcn',@cb_clearAllData, ...
           'BackgroundColor',[0.55 0.08 0.08],'FontColor','white');
       btnClearAll.Layout.Row=r; btnClearAll.Layout.Column=[1 2];

% Sanity check: r must equal nCtrlRows
assert(r == nCtrlRows, 'nCtrlRows mismatch: r=%d, nCtrlRows=%d', r, nCtrlRows);

% -- Load persisted working params (overrides hardcoded defaults) ----------
if ispref('CalciumSpikeGUI','workingParams')
    try, applyWorkingParams(getpref('CalciumSpikeGUI','workingParams')); catch, end
end

%% -- Right: plot panel ----------------------------------------------------
rightGL = uigridlayout(rootGL, [3 1]);
rightGL.Layout.Column = 2;
rightGL.RowHeight     = {'0.3x','0.4x','0.3x'};
rightGL.Padding       = [0 0 0 0];
rightGL.RowSpacing    = 6;

% Row 1: ROI video crop
axCrop = uiaxes(rightGL);
axCrop.Layout.Row            = 1;
axCrop.Title.String          = 'ROI Crop';
axCrop.Title.Interpreter     = 'none';
axCrop.XLabel.String         = '';
axCrop.YLabel.String         = '';
colormap(axCrop, 'gray');
axis(axCrop,'image');

% Row 2: dF/F trace with spikes
axTrace = uiaxes(rightGL);
axTrace.Layout.Row           = 2;
axTrace.Title.String         = 'dF/F Trace  |  Left-click: add/delete spikes  |  N: navigate mode';
axTrace.Title.Interpreter    = 'none';
axTrace.XLabel.String        = 'Time (s)';
axTrace.XLabel.Interpreter   = 'none';
axTrace.YLabel.String        = 'dF/F';
axTrace.YLabel.Interpreter   = 'none';
hold(axTrace,'on'); grid(axTrace,'on'); box(axTrace,'on');

% Row 3: placeholder for future raster
lblFuture = uilabel(rightGL,'Text','[ Future: population raster ]', ...
    'HorizontalAlignment','center','FontSize',14,'FontColor',[0.5 0.5 0.5]);
lblFuture.Layout.Row = 3;

% Persistent overlay handles (updated in-place)
hTraceLine = gobjects(0);   % dF/F trace line
hSpikeSc   = gobjects(0);   % spike scatter markers
hThreshLn  = gobjects(0);   % threshold horizontal line

% Keyboard shortcut on figure
fig.KeyPressFcn = @cb_key;

%% =========================================================================
%% CALLBACKS
%% =========================================================================

    % -- Browse master folder and build session list -----------------------
    function cb_loadMaster(~,~)
        master = uigetdir(S.masterFolder,'Select Master Folder');
        if isequal(master,0), return; end
        S.masterFolder = master;

        sessions = scanForDFFSessions(master);
        if isempty(sessions)
            uialert(fig,'No *_dFF.mat files found in any subfolder.','No Sessions');
            return;
        end
        S.sessionList = sessions;
        S.sessionIdx  = 1;

        [~, mName] = fileparts(master);
        lblMasterPath.Text = mName;
        loadCurrentSession();
    end

    % -- Navigate to previous session --------------------------------------
    function cb_prevSession(~,~)
        if S.sessionIdx <= 1, return; end
        S.sessionIdx = S.sessionIdx - 1;
        loadCurrentSession();
    end

    % -- Navigate to next session ------------------------------------------
    function cb_nextSession(~,~)
        if S.sessionIdx >= numel(S.sessionList), return; end
        S.sessionIdx = S.sessionIdx + 1;
        loadCurrentSession();
    end

    % -- Load the session at S.sessionIdx ----------------------------------
    function loadCurrentSession()
        sess         = S.sessionList{S.sessionIdx};
        S.folderPath = sess.folder;
        S.dffFile    = sess.dffFile;
        S.samFile    = sess.samFile;
        S.tifFile    = sess.tifFile;
        S.hasMask    = ~isempty(sess.samFile) && isfile(sess.samFile);
        S.hasTif     = ~isempty(sess.tifFile) && isfile(sess.tifFile);
        S.maskL      = [];
        S.meanImg    = [];
        S.dFFout     = [];

        % Show current folder name
        [~, folderName] = fileparts(S.folderPath);
        lblCurrentFolder.Text = folderName;

        % Load dFF data
        try
            loaded = load(S.dffFile);
        catch ME
            uialert(fig, ME.message, 'Load Error'); return;
        end

        if isfield(loaded,'dFF')
            S.dFF = loaded.dFF;
        else
            uialert(fig,'No dFF variable in file.','Load Error'); return;
        end

        [T, N] = size(S.dFF);
        S.nROI = N;

        % Extract fps from params if available
        if isfield(loaded,'params') && isfield(loaded.params,'FPS')
            S.fps = loaded.params.FPS;
            efFps.Value = S.fps;
        else
            S.fps = efFps.Value;
        end

        % Reconstruct time vector
        if isfield(loaded,'dFFout') && isfield(loaded.dFFout,'t_dff')
            S.t_img  = loaded.dFFout.t_dff(:);
            S.dFFout = loaded.dFFout;
        else
            S.t_img = (0:T-1)' / S.fps;
        end

        % Sanity: trim t_img if size mismatch
        if length(S.t_img) ~= T
            S.t_img = (0:T-1)' / S.fps;
        end

        % Load cpSAM mask
        if S.hasMask
            try
                samLoaded = load(S.samFile, 'maskL');
                if isfield(samLoaded,'maskL')
                    S.maskL = samLoaded.maskL;
                else
                    S.hasMask = false;
                end
            catch
                S.hasMask = false;
            end
        end

        % Load mean image from TIF (first frame as fallback for crop display)
        if S.hasTif
            try
                tifInfo = imfinfo(S.tifFile);
                S.meanImg = double(imread(S.tifFile, 1, 'Info', tifInfo));
            catch
                S.hasTif = false;
            end
        end

        % Initialize spike storage
        S.spikes = cell(N, 1);
        S.curROI = 1;

        % Update ROI spinner limits
        spnROI.Limits = [1 max(1, N)];
        spnROI.Value  = 1;

        % Check for previously saved results
        matFile = fullfile(S.folderPath, 'ca_spike_data.mat');
        if isfile(matFile)
            try
                saved = load(matFile, 'roi_spikes', 'spike_params');
                if isfield(saved,'roi_spikes')
                    rs = saved.roi_spikes;
                    for ii = 1:min(numel(rs), N)
                        if isfield(rs(ii),'spike_idx') && ~isempty(rs(ii).spike_idx)
                            S.spikes{ii} = rs(ii).spike_idx(:);
                        end
                    end
                end
                if isfield(saved,'spike_params')
                    sp = saved.spike_params;
                    if isfield(sp,'fps'),       efFps.Value      = sp.fps;       end
                    if isfield(sp,'threshold'),  efThresh.Value   = sp.threshold;  end
                    if isfield(sp,'minDist'),     efMinDist.Value  = sp.minDist;    end
                    if isfield(sp,'minWidth'),    efMinWidth.Value = sp.minWidth;   end
                    if isfield(sp,'snapWin'),     efSnap.Value     = sp.snapWin;    end
                    if isfield(sp,'method'),      ddMethod.Value   = sp.method;     end
                    if isfield(sp,'diameter'),    efDiameter.Value = sp.diameter;   end
                    if isfield(sp,'zoomFactor'),  efZoom.Value     = sp.zoomFactor; end
                end
                lblCurrentFolder.Text = [folderName '  (saved)'];
            catch
                % Corrupted mat -- fall through
            end
        end

        updateNavDisplay();
        refreshAll();
    end

    % -- Update folder counter display -------------------------------------
    function updateNavDisplay()
        n = numel(S.sessionList);
        nDone = sum(cellfun(@(s) isfile(fullfile(s.folder,'ca_spike_data.mat')), ...
                            S.sessionList));
        lblNavStatus.Text = sprintf('%d / %d  sessions  |  %d complete', ...
                                    S.sessionIdx, n, nDone);
    end

    % -- ROI navigation callbacks ------------------------------------------
    function cb_prevROI(~,~)
        if S.curROI <= 1, return; end
        S.curROI = S.curROI - 1;
        spnROI.Value = S.curROI;
        refreshAll();
    end

    function cb_nextROI(~,~)
        if S.curROI >= S.nROI, return; end
        S.curROI = S.curROI + 1;
        spnROI.Value = S.curROI;
        refreshAll();
    end

    function cb_jumpROI(~,~)
        newVal = round(spnROI.Value);
        if newVal < 1, newVal = 1; end
        if newVal > S.nROI, newVal = S.nROI; end
        S.curROI = newVal;
        spnROI.Value = newVal;
        refreshAll();
    end

    % -- Detect spikes (current ROI) ---------------------------------------
    function cb_detect(~,~)
        if isempty(S.dFF)
            uialert(fig,'Load data first.','No Data'); return;
        end
        detectSpikesForROI(S.curROI);
        saveWorkingParams();
        refreshTrace();
    end

    % -- Detect spikes (all ROIs) ------------------------------------------
    function cb_detectAll(~,~)
        if isempty(S.dFF)
            uialert(fig,'Load data first.','No Data'); return;
        end
        for ii = 1:S.nROI
            detectSpikesForROI(ii);
        end
        saveWorkingParams();
        refreshTrace();
    end

    % -- Core spike detection for one ROI ----------------------------------
    function detectSpikesForROI(roiIdx)
        trace = S.dFF(:, roiIdx);
        method = ddMethod.Value;

        switch method
            case 'Threshold'
                thr      = efThresh.Value;
                minDistS = efMinDist.Value;
                minWS    = efMinWidth.Value;
                fps      = efFps.Value;
                spk = find_spikes_threshold(trace, fps, thr, minDistS, minWS);

            case 'OASIS'
                % Try loading pre-computed spikes from dFFout
                spk = [];
                if ~isempty(S.dFFout) && isfield(S.dFFout,'spikes_oasis') ...
                        && ~isempty(S.dFFout.spikes_oasis)
                    oasisSpk = S.dFFout.spikes_oasis(:, roiIdx);
                    thr = efThresh.Value;
                    % Threshold the OASIS spike train
                    candidates = find(oasisSpk > thr);
                    if ~isempty(candidates)
                        % Apply min distance: keep only local maxima
                        minDistFr = max(1, round(efMinDist.Value * efFps.Value));
                        spk = enforce_min_distance(oasisSpk, candidates, minDistFr);
                    end
                end
                if isempty(spk)
                    % Fallback: run threshold on dF/F trace
                    thr      = efThresh.Value;
                    minDistS = efMinDist.Value;
                    minWS    = efMinWidth.Value;
                    fps      = efFps.Value;
                    spk = find_spikes_threshold(trace, fps, thr, minDistS, minWS);
                end
        end

        S.spikes{roiIdx} = spk(:);
    end

    % -- Axes click: add or delete spike -----------------------------------
    function cb_axClick(~,~)
        if isempty(S.dFF) || S.nROI == 0, return; end

        cp     = axTrace.CurrentPoint;
        clickT = cp(1,1);
        if clickT < S.t_img(1) || clickT > S.t_img(end), return; end

        % Map time to nearest frame index
        [~, nearIdx] = min(abs(S.t_img - clickT));
        trace = S.dFF(:, S.curROI);
        curSpk = S.spikes{S.curROI};

        switch S.editMode
            case 'add'
                % Snap to local maximum within snap window
                snapFr = max(1, round(efSnap.Value * efFps.Value));
                lo = max(1, nearIdx - snapFr);
                hi = min(length(trace), nearIdx + snapFr);
                [~, imax] = max(trace(lo:hi));
                peakIdx = lo + imax - 1;

                % Add only if not a duplicate (within 2 frames)
                if isempty(curSpk) || min(abs(curSpk - peakIdx)) > 2
                    curSpk = sort([curSpk; peakIdx]);
                end

            case 'delete'
                if isempty(curSpk), return; end
                [~, ci] = min(abs(curSpk - nearIdx));
                curSpk(ci) = [];
        end

        S.spikes{S.curROI} = curSpk;
        redraw_spikes();
        updateSpikeCount();
    end

    % -- Keyboard shortcuts ------------------------------------------------
    function cb_key(~, event)
        switch lower(event.Key)
            case 'n', cb_setNav();
            case 'a', cb_setAdd();
            case 'd', cb_setDelete();
            case 'rightarrow', cb_nextROI();
            case 'leftarrow',  cb_prevROI();
        end
    end

    % -- Edit mode setters -------------------------------------------------
    function cb_setAdd(~,~)
        S.editMode = 'add';
        lblMode.Text      = 'Mode: + ADD  (click to snap to local max)';
        lblMode.FontColor = [0.30 0.62 1.00];
        disableDefaultInteractivity(axTrace);
        axTrace.ButtonDownFcn = @cb_axClick;
    end

    function cb_setDelete(~,~)
        S.editMode = 'delete';
        lblMode.Text      = 'Mode: - DELETE  (click to remove nearest)';
        lblMode.FontColor = [1.00 0.35 0.25];
        disableDefaultInteractivity(axTrace);
        axTrace.ButtonDownFcn = @cb_axClick;
    end

    function cb_setNav(~,~)
        S.editMode = 'none';
        lblMode.Text      = 'Mode: Navigate';
        lblMode.FontColor = [0.90 0.80 0.10];
        axTrace.ButtonDownFcn = [];
        enableDefaultInteractivity(axTrace);
    end

    function cb_reset(~,~)
        applyWorkingParams(DEFAULTS);
        saveWorkingParams();
        if ~isempty(S.dFF)
            refreshAll();
        end
    end

    function cb_clearSpikes(~,~)
        if S.nROI == 0, return; end
        S.spikes{S.curROI} = [];
        redraw_spikes();
        updateSpikeCount();
    end

    function cb_clearAllData(~,~)
        if isempty(S.sessionList)
            uialert(fig,'No master folder loaded.','Error'); return;
        end
        n = numel(S.sessionList);
        answer = uiconfirm(fig, ...
            sprintf(['Delete ca_spike_data.mat from all %d session folders?\n' ...
                     'Also remove entries from ca_spike_master.mat.\n\n' ...
                     'This cannot be undone.'], n), ...
            'Confirm Clear All Data', ...
            'Options',    {'Delete All','Cancel'}, ...
            'DefaultOption','Cancel', ...
            'CancelOption','Cancel', ...
            'Icon','warning');
        if ~strcmp(answer,'Delete All'), return; end

        nDeleted = 0;
        for i = 1:n
            folder = S.sessionList{i}.folder;
            fp = fullfile(folder, 'ca_spike_data.mat');
            if isfile(fp)
                delete(fp);
                nDeleted = nDeleted + 1;
            end
        end

        % Remove master file
        masterFile = fullfile(S.masterFolder, 'ca_spike_master.mat');
        if isfile(masterFile)
            delete(masterFile);
            nDeleted = nDeleted + 1;
        end

        % Clear current session spikes
        S.spikes = cell(S.nROI, 1);
        redraw_spikes();
        updateSpikeCount();
        updateNavDisplay();

        uialert(fig, sprintf('Deleted %d file(s) across %d folders.', nDeleted, n), ...
            'Done','Icon','success');
    end

    function applyWorkingParams(p)
        efFps.Value      = p.fps;
        ddMethod.Value   = p.method;
        efThresh.Value   = p.threshold;
        efMinDist.Value  = p.minDist;
        efMinWidth.Value = p.minWidth;
        efSnap.Value     = p.snapWin;
        efDiameter.Value = p.diameter;
        efZoom.Value     = p.zoomFactor;
    end

    function saveWorkingParams()
        wp.fps        = efFps.Value;
        wp.method     = ddMethod.Value;
        wp.threshold  = efThresh.Value;
        wp.minDist    = efMinDist.Value;
        wp.minWidth   = efMinWidth.Value;
        wp.snapWin    = efSnap.Value;
        wp.diameter   = efDiameter.Value;
        wp.zoomFactor = efZoom.Value;
        setpref('CalciumSpikeGUI','workingParams',wp);
    end

    % -- Save --------------------------------------------------------------
    function cb_save(~,~)
        if isempty(S.folderPath) || S.nROI == 0
            uialert(fig,'No data loaded.','Error'); return;
        end

        T = size(S.dFF, 1);
        fps = efFps.Value;

        % Build roi_spikes struct array
        roi_spikes = struct('spike_idx',{},'spike_t',{},'spike_amp',{}, ...
                            'spike_train',{},'n_spikes',{});
        for ii = 1:S.nROI
            spk = S.spikes{ii};
            rs = struct();
            if isempty(spk)
                rs.spike_idx   = [];
                rs.spike_t     = [];
                rs.spike_amp   = [];
                rs.spike_train = zeros(T, 1);
                rs.n_spikes    = 0;
            else
                spk = spk(spk >= 1 & spk <= T);
                rs.spike_idx   = spk(:);
                rs.spike_t     = S.t_img(spk);
                rs.spike_amp   = S.dFF(spk, ii);
                train          = zeros(T, 1);
                train(spk)     = 1;
                rs.spike_train = train;
                rs.n_spikes    = numel(spk);
            end
            roi_spikes(ii) = rs; %#ok<AGROW>
        end

        % Build spike_params struct
        spike_params = struct();
        spike_params.method     = ddMethod.Value;
        spike_params.threshold  = efThresh.Value;
        spike_params.minDist    = efMinDist.Value;
        spike_params.minWidth   = efMinWidth.Value;
        spike_params.snapWin    = efSnap.Value;
        spike_params.fps        = fps;
        spike_params.diameter   = efDiameter.Value;
        spike_params.zoomFactor = efZoom.Value;
        spike_params.dff_file   = S.dffFile;
        spike_params.saved_at   = datestr(now); %#ok<TNOW1,DATST>

        % Save per-session file
        outFile = fullfile(S.folderPath, 'ca_spike_data.mat');
        save(outFile, 'roi_spikes', 'spike_params', '-v7.3');

        % Update master file
        updateMasterFile(roi_spikes, spike_params);

        [~, folderName] = fileparts(S.folderPath);
        lblCurrentFolder.Text = [folderName '  (saved)'];
        updateNavDisplay();

        uialert(fig, sprintf('Saved %d ROIs to:\n%s', S.nROI, outFile), ...
            'Saved','Icon','success');
    end

    function updateMasterFile(roi_spikes, spike_params)
        masterFile = fullfile(S.masterFolder, 'ca_spike_master.mat');
        [~, sessName] = fileparts(S.folderPath);

        entry = struct();
        entry.folder       = S.folderPath;
        entry.session_name = sessName;
        entry.n_rois       = S.nROI;
        entry.roi_spikes   = roi_spikes;
        entry.spike_params = spike_params;

        if isfile(masterFile)
            try
                loaded = load(masterFile, 'master');
                master = loaded.master;
            catch
                master = struct('sessions',{});
            end
        else
            master = struct('sessions',{});
        end

        % Find existing entry for this folder, or append
        found = false;
        if isfield(master,'sessions') && ~isempty(master.sessions)
            for ii = 1:numel(master.sessions)
                if strcmp(master.sessions(ii).folder, S.folderPath)
                    master.sessions(ii) = entry;
                    found = true;
                    break;
                end
            end
        end
        if ~found
            if isempty(master.sessions)
                master.sessions = entry;
            else
                master.sessions(end+1) = entry;
            end
        end

        save(masterFile, 'master', '-v7.3');
    end

%% =========================================================================
%% PLOT HELPERS
%% =========================================================================

    function refreshAll()
        refreshCrop();
        refreshTrace();
    end

    function refreshCrop()
        cla(axCrop);
        if S.nROI == 0, return; end

        [~, sessName] = fileparts(S.folderPath);
        roiIdx = S.curROI;

        % Update ROI label
        lblROI.Text = sprintf('ROI: %d / %d', roiIdx, S.nROI);

        if ~S.hasMask || isempty(S.maskL)
            axCrop.Title.String = sprintf('ROI %d  --  %s  (no mask)', ...
                roiIdx, strrep(sessName,'_',' '));
            axCrop.Title.Interpreter = 'none';
            return;
        end

        % Get ROI centroid from mask
        labels = unique(S.maskL(:));
        labels(labels == 0) = [];
        if roiIdx > numel(labels)
            axCrop.Title.String = sprintf('ROI %d  --  %s  (out of range)', ...
                roiIdx, strrep(sessName,'_',' '));
            axCrop.Title.Interpreter = 'none';
            return;
        end

        roiMask = (S.maskL == labels(roiIdx));
        props = regionprops(roiMask, 'Centroid');
        if isempty(props)
            axCrop.Title.String = sprintf('ROI %d  --  %s  (empty mask)', ...
                roiIdx, strrep(sessName,'_',' '));
            axCrop.Title.Interpreter = 'none';
            return;
        end

        cx = round(props(1).Centroid(1));
        cy = round(props(1).Centroid(2));

        % Choose image source: mean projection if available, else mask
        if S.hasTif && ~isempty(S.meanImg)
            baseImg = S.meanImg;
        else
            baseImg = double(S.maskL > 0) * 100;
        end

        [H, W] = size(baseImg);
        halfSide = round(efDiameter.Value * efZoom.Value / 2);
        r1 = max(1, cy - halfSide);
        r2 = min(H, cy + halfSide);
        c1 = max(1, cx - halfSide);
        c2 = min(W, cx + halfSide);

        cropImg  = baseImg(r1:r2, c1:c2);
        cropMask = roiMask(r1:r2, c1:c2);

        imagesc(axCrop, cropImg);
        colormap(axCrop, 'gray');
        axis(axCrop, 'image');
        hold(axCrop, 'on');

        % Overlay ROI contour in red
        if any(cropMask(:))
            bnd = bwboundaries(cropMask);
            for b = 1:numel(bnd)
                plot(axCrop, bnd{b}(:,2), bnd{b}(:,1), 'r-', 'LineWidth', 1.5);
            end
        end

        hold(axCrop, 'off');
        axCrop.Title.String = sprintf('ROI %d  --  %s', ...
            roiIdx, strrep(sessName,'_',' '));
        axCrop.Title.Interpreter = 'none';
    end

    function refreshTrace()
        cla(axTrace); hold(axTrace,'on');
        if isempty(S.dFF) || S.nROI == 0, return; end

        roiIdx = S.curROI;
        trace  = S.dFF(:, roiIdx);

        % Update ROI label
        lblROI.Text = sprintf('ROI: %d / %d', roiIdx, S.nROI);

        % Plot dF/F trace
        hTraceLine = plot(axTrace, S.t_img, trace, 'Color',[0 0 0], 'LineWidth',0.9);
        hTraceLine.HitTest = 'off';

        % Threshold dashed line
        hThreshLn = yline(axTrace, efThresh.Value, '--', ...
            'Color',[0.85 0.15 0.15], 'LineWidth',0.8, 'Label','threshold', ...
            'LabelHorizontalAlignment','left');
        hThreshLn.HitTest  = 'off';

        % Pre-allocate spike scatter (empty data)
        hSpikeSc = scatter(axTrace, [], [], 75, 'v', 'filled', ...
            'MarkerFaceColor',[1.00 0.28 0.28], 'MarkerEdgeColor','none', ...
            'Tag','spike_marker');
        hSpikeSc.HitTest = 'off';

        xlim(axTrace, [S.t_img(1) S.t_img(end)]);
        [~, sessName] = fileparts(S.folderPath);
        axTrace.Title.String = sprintf('ROI %d / %d  --  %s', ...
            roiIdx, S.nROI, strrep(sessName,'_',' '));
        axTrace.Title.Interpreter = 'none';

        redraw_spikes();
        updateSpikeCount();
    end

    function redraw_spikes()
        % Fast in-place update: just set XData/YData
        if S.nROI == 0, return; end

        roiIdx = S.curROI;
        trace  = S.dFF(:, roiIdx);
        curSpk = S.spikes{roiIdx};

        valid = [];
        if ~isempty(curSpk)
            valid = curSpk(curSpk >= 1 & curSpk <= length(trace));
        end

        if isgraphics(hSpikeSc)
            if isempty(valid)
                set(hSpikeSc, 'XData',[], 'YData',[]);
            else
                set(hSpikeSc, 'XData',S.t_img(valid), 'YData',trace(valid));
            end
        end

        updateSpikeCount();
    end

    function updateSpikeCount()
        if S.nROI == 0
            lblCount.Text = 'Spikes: 0';
            return;
        end
        curSpk = S.spikes{S.curROI};
        lblCount.Text = sprintf('Spikes: %d', numel(curSpk));
    end

end  % calcium_spike_gui


%% =========================================================================
%% LOCAL HELPER -- threshold-based spike finder
%% =========================================================================
function locs = find_spikes_threshold(trace, fps, threshold, minDistS, minWS)
%FIND_SPIKES_THRESHOLD  Detect calcium spikes by threshold crossing + local max.
%
%   1. Find contiguous above-threshold regions.
%   2. Within each region, take the local maximum as the spike location.
%   3. Filter by minimum distance between spikes.
%   4. Filter by minimum width of the above-threshold region.

trace = trace(:);
T     = length(trace);
locs  = [];

if T < 3, return; end

% -- 1. Find above-threshold crossings ------------------------------------
above = trace > threshold;
if ~any(above), return; end

% Find contiguous regions
d_above    = diff([0; above; 0]);
region_on  = find(d_above == 1);
region_off = find(d_above == -1) - 1;  % last index in region

nRegions = numel(region_on);
if nRegions == 0, return; end

% -- 2. Min width filter ---------------------------------------------------
minWFr = max(1, round(minWS * fps));
keep   = (region_off - region_on + 1) >= minWFr;
region_on  = region_on(keep);
region_off = region_off(keep);

if isempty(region_on), return; end

% -- 3. Local max within each region --------------------------------------
nReg = numel(region_on);
locs = zeros(nReg, 1);
for i = 1:nReg
    seg = trace(region_on(i):region_off(i));
    [~, imax] = max(seg);
    locs(i) = region_on(i) + imax - 1;
end

% -- 4. Min distance: greedily keep highest, remove neighbours -------------
minDistFr = max(1, round(minDistS * fps));
if minDistFr > 1 && numel(locs) > 1
    [~, ord]  = sort(trace(locs), 'descend');
    sorted    = locs(ord);
    keepFlag  = true(length(sorted), 1);
    for i = 1:length(sorted)
        if ~keepFlag(i), continue; end
        tooClose = abs(sorted - sorted(i)) < minDistFr;
        tooClose(i) = false;
        keepFlag(tooClose) = false;
    end
    locs = sort(sorted(keepFlag));
end

locs = locs(:);
end  % find_spikes_threshold


%% =========================================================================
%% LOCAL HELPER -- enforce minimum distance on OASIS spike candidates
%% =========================================================================
function locs = enforce_min_distance(spkTrain, candidates, minDistFr)
%ENFORCE_MIN_DISTANCE  From candidate indices, keep highest-amplitude with min gap.

candidates = candidates(:);
amps       = spkTrain(candidates);
[~, ord]   = sort(amps, 'descend');
sorted     = candidates(ord);
keepFlag   = true(length(sorted), 1);

for i = 1:length(sorted)
    if ~keepFlag(i), continue; end
    tooClose = abs(sorted - sorted(i)) < minDistFr;
    tooClose(i) = false;
    keepFlag(tooClose) = false;
end

locs = sort(sorted(keepFlag));
locs = locs(:);
end  % enforce_min_distance


%% =========================================================================
%% LOCAL HELPER -- scan for dFF sessions
%% =========================================================================
function sessions = scanForDFFSessions(masterFolder)
%SCANFORDFF SESSIONS  Recursively find all subfolders containing a *_dFF.mat file.
%   Returns a cell array of structs with fields:
%       .folder   -- path to the subfolder
%       .dffFile  -- path to the *_dFF.mat file
%       .samFile  -- path to the *_cpSAM_output.mat file (or '' if not found)
%       .tifFile  -- path to the *_minusDark_MC.tif file (or '' if not found)
%   Sorted alphabetically by folder path.

sessions = {};

% Recursive glob for dFF.mat files
hits = dir(fullfile(masterFolder, '**', '*_dFF.mat'));
if isempty(hits), return; end

% One entry per unique folder (take first dFF if multiple exist)
folders  = {hits.folder};
[uFolders, ia] = unique(folders, 'stable');

for i = 1:numel(uFolders)
    h = hits(ia(i));
    folder  = h.folder;
    dffFile = fullfile(folder, h.name);

    % Look for paired cpSAM_output.mat
    samHits = dir(fullfile(folder, '*_cpSAM_output.mat'));
    if ~isempty(samHits)
        samFile = fullfile(folder, samHits(1).name);
    else
        samFile = '';
    end

    % Look for paired minusDark_MC.tif
    tifHits = dir(fullfile(folder, '*_minusDark_MC.tif'));
    if ~isempty(tifHits)
        tifFile = fullfile(folder, tifHits(1).name);
    else
        tifFile = '';
    end

    sessions{end+1} = struct('folder', folder, ...
                             'dffFile', dffFile, ...
                             'samFile', samFile, ...
                             'tifFile', tifFile); %#ok<AGROW>
end
end  % scanForDFFSessions
