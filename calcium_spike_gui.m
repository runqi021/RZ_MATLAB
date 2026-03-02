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

% Video / sum-projection state
S.tifInfo      = [];    % cached imfinfo result (so it's only called once per session)
S.sumImg       = [];    % [H x W] mean projection over sampled frames
S.vidRawFrames = [];    % [Hc x Wc x nF] raw double crops (for contrast sliders)
S.vidBndMask   = [];    % [Hc x Wc] logical ROI boundary pixels
S.vidVmin      = 0;     % contrast min for video (from slider)
S.vidVmax      = 1;     % contrast max for video (from slider)
S.vidTimer     = [];    % timer object for video playback
S.vidFrIdx     = 1;     % current video frame index
S.hVidImg      = [];    % handle to the image object in axVideo (for fast CData updates)
S.nDropFrames  = 0;     % frames dropped in dFF, used to align TIF start
S.dFFraw       = [];    % copy of dFF as loaded (before optional detrend)
S.detrend      = false; % is detrend currently applied?
S.hTraceCursor = [];    % blue vertical xline on axTrace (current video frame)
S.hOasisLine   = [];    % orange overlay line for OASIS deconvolved trace

% -- Factory defaults (restored by Reset button) --------------------------
DEFAULTS.fps        = 30;
DEFAULTS.method     = 'Raw dFF';
DEFAULTS.threshold  = 2.0;
DEFAULTS.minDist    = 0.5;
DEFAULTS.minWidth   = 0.1;
DEFAULTS.snapWin    = 0.20;
DEFAULTS.diameter   = 30;
DEFAULTS.zoomFactor = 4;
DEFAULTS.prominence = 0;      % MinPeakProminence for findpeaks (0 = disabled)
DEFAULTS.tauDecay   = 0.45;   % OASIS AR(1) decay constant (s)
DEFAULTS.pythonExe  = '';      % filled at runtime by detectPyExe()

% Auto-detect Python executable
DEFAULTS.pythonExe = detectPyExe();

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

nCtrlRows = 37;
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

r=r+1; h=uilabel(cGL,'Text','Python Exe:');
       h.Layout.Row=r; h.Layout.Column=1;
       efPythonExe = uieditfield(cGL,'text','Value',DEFAULTS.pythonExe, ...
           'ValueChangedFcn',@(~,~) saveWorkingParams());
       efPythonExe.Layout.Row=r; efPythonExe.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Tau Decay (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efTauDecay = uieditfield(cGL,'numeric','Value',DEFAULTS.tauDecay, ...
           'Limits',[0.01 10],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efTauDecay.Layout.Row=r; efTauDecay.Layout.Column=2;

% -- Spike Detection -------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Spike Detection --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Detect on:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddMethod = uidropdown(cGL,'Items',{'Raw dFF','OASIS deconv'},'Value','Raw dFF', ...
           'ValueChangedFcn',@(~,~) saveWorkingParams());
       ddMethod.Layout.Row=r; ddMethod.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Threshold (dF/F):');
       h.Layout.Row=r; h.Layout.Column=1;
       efThresh = uieditfield(cGL,'numeric','Value',DEFAULTS.threshold, ...
           'Limits',[-Inf Inf],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efThresh.Layout.Row=r; efThresh.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Min Prominence:');
       h.Layout.Row=r; h.Layout.Column=1;
       efProminence = uieditfield(cGL,'numeric','Value',DEFAULTS.prominence, ...
           'Limits',[0 Inf],'ValueChangedFcn',@(~,~) saveWorkingParams());
       efProminence.Layout.Row=r; efProminence.Layout.Column=2;

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

r=r+1; btnRunOasis = uibutton(cGL,'Text','Run OASIS Deconv (All ROIs)', ...
           'ButtonPushedFcn',@cb_runOasis, ...
           'BackgroundColor',[0.55 0.30 0.05],'FontColor','white');
       btnRunOasis.Layout.Row=r; btnRunOasis.Layout.Column=[1 2];

r=r+1; btnDetrend = uibutton(cGL,'Text','Detrend: OFF', ...
           'ButtonPushedFcn',@cb_toggleDetrend, ...
           'BackgroundColor',[0.25 0.25 0.55],'FontColor','white');
       btnDetrend.Layout.Row=r; btnDetrend.Layout.Column=[1 2];

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
rightGL = uigridlayout(rootGL, [4 1]);
rightGL.Layout.Column = 2;
rightGL.RowHeight     = {'0.27x', 45, '0.40x', '0.27x'};
rightGL.Padding       = [0 0 0 0];
rightGL.RowSpacing    = 6;

% Row 1: Nested 1x4 grid for Sum Projection + contrast + Video + contrast
topGL = uigridlayout(rightGL, [1 4]);
topGL.Layout.Row    = 1;
topGL.ColumnWidth   = {'1x', 65, '1x', 65};
topGL.Padding       = [0 0 0 0];
topGL.ColumnSpacing = 4;

axSum = uiaxes(topGL);
axSum.Layout.Column          = 1;
axSum.Title.String           = 'Avg Projection';
axSum.Title.Interpreter      = 'none';
axSum.XLabel.String          = '';
axSum.YLabel.String          = '';
colormap(axSum, 'gray');
axis(axSum, 'image');

% Sum contrast sliders (col 2)
sumSldGL = uigridlayout(topGL, [4 1]);
sumSldGL.Layout.Column = 2;
sumSldGL.RowHeight = {'fit','1x','fit','1x'};
sumSldGL.Padding = [2 4 2 4];

uilabel(sumSldGL,'Text','Max','FontSize',8,'HorizontalAlignment','center');
sldSumMax = uislider(sumSldGL, 'Value',1, 'Limits',[0 1], ...
    'ValueChangedFcn', @(~,~) cb_sumContrastChanged());
sldSumMax.Layout.Row = 2;

uilabel(sumSldGL,'Text','Min','FontSize',8,'HorizontalAlignment','center');
sldSumMin = uislider(sumSldGL, 'Value',0, 'Limits',[0 1], ...
    'ValueChangedFcn', @(~,~) cb_sumContrastChanged());
sldSumMin.Layout.Row = 4;

axVideo = uiaxes(topGL);
axVideo.Layout.Column        = 3;
axVideo.Title.String         = 'Video';
axVideo.Title.Interpreter    = 'none';
axVideo.XLabel.String        = '';
axVideo.YLabel.String        = '';
colormap(axVideo, 'gray');
axis(axVideo, 'image');

% Video contrast sliders (col 4)
vidSldGL = uigridlayout(topGL, [4 1]);
vidSldGL.Layout.Column = 4;
vidSldGL.RowHeight = {'fit','1x','fit','1x'};
vidSldGL.Padding = [2 4 2 4];

uilabel(vidSldGL,'Text','Max','FontSize',8,'HorizontalAlignment','center');
sldVidMax = uislider(vidSldGL, 'Value',1, 'Limits',[0 1], ...
    'ValueChangedFcn', @(~,~) cb_vidContrastChanged());
sldVidMax.Layout.Row = 2;

uilabel(vidSldGL,'Text','Min','FontSize',8,'HorizontalAlignment','center');
sldVidMin = uislider(vidSldGL, 'Value',0, 'Limits',[0 1], ...
    'ValueChangedFcn', @(~,~) cb_vidContrastChanged());
sldVidMin.Layout.Row = 4;

% Row 2: Time scrub slider
scrubGL = uigridlayout(rightGL, [1 2]);
scrubGL.Layout.Row = 2;
scrubGL.ColumnWidth = {'fit', '1x'};
scrubGL.Padding = [4 2 4 2];

uilabel(scrubGL, 'Text','Time:', 'FontSize',9);
sldTime = uislider(scrubGL, 'Value',1, 'Limits',[1 2], ...
    'ValueChangedFcn', @cb_scrubTime, ...
    'MajorTicksMode','auto', 'MinorTicksMode','auto');
sldTime.Layout.Column = 2;

% Row 3: dF/F trace with spikes
axTrace = uiaxes(rightGL);
axTrace.Layout.Row           = 3;
axTrace.Title.String         = 'dF/F Trace  |  Left-click: add/delete spikes  |  N: navigate mode';
axTrace.Title.Interpreter    = 'none';
axTrace.XLabel.String        = 'Time (s)';
axTrace.XLabel.Interpreter   = 'none';
axTrace.YLabel.String        = 'dF/F';
axTrace.YLabel.Interpreter   = 'none';
hold(axTrace,'on'); grid(axTrace,'on'); box(axTrace,'on');

% Row 4: placeholder for future raster
lblFuture = uilabel(rightGL,'Text','[ Future: population raster ]', ...
    'HorizontalAlignment','center','FontSize',14,'FontColor',[0.5 0.5 0.5]);
lblFuture.Layout.Row = 4;

% Persistent overlay handles (updated in-place)
hTraceLine = gobjects(0);   % dF/F trace line
hSpikeSc   = gobjects(0);   % spike scatter markers
hThreshLn  = gobjects(0);   % threshold horizontal line

% Keyboard shortcut on figure
fig.KeyPressFcn = @cb_key;

% Close callback: stop video timer before destroying figure
fig.CloseRequestFcn = @cb_close;

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
        S.sumImg     = [];
        S.tifInfo    = [];
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

        % Extract nDropFrames from saved params (aligns TIF with dFF)
        S.nDropFrames = 0;
        if isfield(loaded,'dFFout') && isfield(loaded.dFFout,'params') ...
                && isfield(loaded.dFFout.params,'nDropFrames')
            S.nDropFrames = loaded.dFFout.params.nDropFrames;
        end

        % Load TIF: cache imfinfo, compute avg projection from sampled frames
        S.tifInfo = [];
        S.sumImg  = [];
        S.meanImg = [];
        if S.hasTif
            try
                S.tifInfo = imfinfo(S.tifFile);
                nTotal   = numel(S.tifInfo);
                startFr  = S.nDropFrames + 1;
                if startFr > nTotal, startFr = 1; end
                nAvail   = nTotal - startFr + 1;
                nSamp    = min(50, nAvail);
                sampIdxs = unique(round(linspace(startFr, nTotal, nSamp)));
                fr1      = double(imread(S.tifFile, sampIdxs(1), 'Info', S.tifInfo));
                acc      = zeros(size(fr1));
                for ki = sampIdxs
                    acc = acc + double(imread(S.tifFile, ki, 'Info', S.tifInfo));
                end
                S.sumImg  = acc / numel(sampIdxs);   % average projection
                S.meanImg = S.sumImg;
            catch
                S.hasTif  = false;
                S.tifInfo = [];
            end
        end

        % Store raw dFF and reset detrend
        S.dFFraw  = S.dFF;
        S.detrend = false;
        btnDetrend.Text = 'Detrend: OFF';
        btnDetrend.BackgroundColor = [0.25 0.25 0.55];

        % Reset OASIS button appearance
        btnRunOasis.Text            = 'Run OASIS Deconv (All ROIs)';
        btnRunOasis.BackgroundColor = [0.55 0.30 0.05];

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
                    if isfield(sp,'fps'),        efFps.Value       = sp.fps;        end
                    if isfield(sp,'threshold'),   efThresh.Value    = sp.threshold;   end
                    if isfield(sp,'minDist'),      efMinDist.Value   = sp.minDist;     end
                    if isfield(sp,'minWidth'),     efMinWidth.Value  = sp.minWidth;    end
                    if isfield(sp,'snapWin'),      efSnap.Value      = sp.snapWin;     end
                    if isfield(sp,'method')
                        mval = sp.method;
                        if strcmp(mval,'Threshold'), mval = 'Raw dFF'; end
                        if strcmp(mval,'OASIS'),     mval = 'OASIS deconv'; end
                        if ismember(mval, ddMethod.Items), ddMethod.Value = mval; end
                    end
                    if isfield(sp,'diameter'),     efDiameter.Value  = sp.diameter;    end
                    if isfield(sp,'zoomFactor'),   efZoom.Value      = sp.zoomFactor;  end
                    if isfield(sp,'prominence'),   efProminence.Value = sp.prominence; end
                    if isfield(sp,'tauDecay'),     efTauDecay.Value  = sp.tauDecay;    end
                    if isfield(sp,'pythonExe'),    efPythonExe.Value = sp.pythonExe;   end
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
        fps     = efFps.Value;
        thr     = efThresh.Value;
        minDist = efMinDist.Value;
        minW    = efMinWidth.Value;
        prom    = efProminence.Value;
        detOn   = ddMethod.Value;   % 'Raw dFF' or 'OASIS deconv'

        switch detOn
            case 'OASIS deconv'
                if ~isempty(S.dFFout) && isfield(S.dFFout,'dFF_oasis_deconv') ...
                        && ~isempty(S.dFFout.dFF_oasis_deconv) ...
                        && size(S.dFFout.dFF_oasis_deconv,2) >= roiIdx
                    trace_det = S.dFFout.dFF_oasis_deconv(:, roiIdx);
                elseif ~isempty(S.dFFout) && isfield(S.dFFout,'spikes_oasis') ...
                        && ~isempty(S.dFFout.spikes_oasis) ...
                        && size(S.dFFout.spikes_oasis,2) >= roiIdx
                    trace_det = S.dFFout.spikes_oasis(:, roiIdx);
                else
                    uialert(fig, ...
                        'No OASIS data found. Press "Run OASIS Deconv" first.', ...
                        'OASIS Not Run');
                    return;
                end
            otherwise   % 'Raw dFF'
                trace_det = S.dFF(:, roiIdx);
        end

        spk = find_spikes_threshold(trace_det, fps, thr, minDist, minW, prom);
        S.spikes{roiIdx} = spk(:);
    end

    % -- Run OASIS deconvolution on all ROIs -------------------------------
    function cb_runOasis(~,~)
        if isempty(S.dFF) || S.nROI == 0
            uialert(fig,'Load data first.','No Data'); return;
        end
        if isempty(S.dFFout) || ~isfield(S.dFFout,'F_dff') || isempty(S.dFFout.F_dff)
            uialert(fig, ...
                ['No raw fluorescence (F_dff) found in the loaded dFF file.\n' ...
                 'Re-run dFF_RZ with UseOASIS=false to ensure F_dff is saved.'], ...
                'F_dff Missing');
            return;
        end

        pyExe = efPythonExe.Value;
        if isempty(pyExe), pyExe = detectPyExe(); end
        if isempty(pyExe) || ~isfile(pyExe)
            uialert(fig, ...
                'Python executable not found. Set it in the Parameters section.', ...
                'Python Not Found');
            return;
        end

        fps      = efFps.Value;
        tauDecay = efTauDecay.Value;
        g = round(exp(-(1/fps)/tauDecay), 2);

        btnRunOasis.Text   = 'Running OASIS...';
        btnRunOasis.Enable = 'off';
        drawnow;
        try
            [F_oasis, dFF_oasis, spikes_oasis, baseline_oasis] = ...
                helper.oasis_deconv_and_dff_AR1(S.dFFout.F_dff, g, 'PythonExe', pyExe);

            S.dFFout.F_oasis_deconv   = F_oasis;
            S.dFFout.dFF_oasis_deconv = dFF_oasis;
            S.dFFout.spikes_oasis     = spikes_oasis;
            S.dFFout.baseline_oasis   = baseline_oasis;
            S.dFFout.g_AR1            = g;

            btnRunOasis.Text            = 'OASIS: Done';
            btnRunOasis.BackgroundColor = [0.10 0.60 0.45];
            btnRunOasis.Enable          = 'on';
            drawnow;

            refreshTrace();   % show orange overlay
        catch ME
            btnRunOasis.Text            = 'Run OASIS Deconv (All ROIs)';
            btnRunOasis.BackgroundColor = [0.55 0.30 0.05];
            btnRunOasis.Enable          = 'on';
            uialert(fig, sprintf('OASIS failed:\n%s', ME.message), 'OASIS Error');
        end
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

    % -- Detrend toggle callback -------------------------------------------
    function cb_toggleDetrend(~,~)
        if isempty(S.dFF) || isempty(S.dFFraw), return; end
        S.detrend = ~S.detrend;
        if S.detrend
            S.dFF = detrend(S.dFFraw);
            btnDetrend.Text = 'Detrend: ON';
            btnDetrend.BackgroundColor = [0.10 0.60 0.45];
        else
            S.dFF = S.dFFraw;
            btnDetrend.Text = 'Detrend: OFF';
            btnDetrend.BackgroundColor = [0.25 0.25 0.55];
        end
        refreshTrace();
    end

    function applyWorkingParams(p)
        efFps.Value      = p.fps;
        mval = p.method;
        if strcmp(mval,'Threshold'), mval = 'Raw dFF'; end
        if strcmp(mval,'OASIS'),     mval = 'OASIS deconv'; end
        if ismember(mval, ddMethod.Items), ddMethod.Value = mval; end
        efThresh.Value   = p.threshold;
        efMinDist.Value  = p.minDist;
        efMinWidth.Value = p.minWidth;
        efSnap.Value     = p.snapWin;
        efDiameter.Value = p.diameter;
        efZoom.Value     = p.zoomFactor;
        if isfield(p,'prominence'),  efProminence.Value = p.prominence;  end
        if isfield(p,'tauDecay'),    efTauDecay.Value   = p.tauDecay;    end
        if isfield(p,'pythonExe'),   efPythonExe.Value  = p.pythonExe;   end
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
        wp.prominence = efProminence.Value;
        wp.tauDecay   = efTauDecay.Value;
        wp.pythonExe  = efPythonExe.Value;
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
                            'spike_train',{},'n_spikes',{},'ifSpike',{});
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
            rs.ifSpike = rs.n_spikes > 0;   % flag: skip in future analysis if false
            roi_spikes(ii) = rs; %#ok<AGROW>
        end

        % Convenience logical vector: ifSpike(ii) = true if ROI ii has any spikes
        ifSpike = logical([roi_spikes.n_spikes] > 0);

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
        spike_params.prominence = efProminence.Value;
        spike_params.tauDecay   = efTauDecay.Value;
        spike_params.dff_file   = S.dffFile;
        spike_params.saved_at   = datestr(now); %#ok<TNOW1,DATST>
        spike_params.oasis_run     = ~isempty(S.dFFout) && isfield(S.dFFout,'spikes_oasis') ...
                                      && ~isempty(S.dFFout.spikes_oasis);
        spike_params.oasis_tauDecay = efTauDecay.Value;
        spike_params.oasis_g        = '';
        if ~isempty(S.dFFout) && isfield(S.dFFout,'g_AR1')
            spike_params.oasis_g = S.dFFout.g_AR1;
        end
        spike_params.detect_on = ddMethod.Value;

        % Save per-session file
        outFile = fullfile(S.folderPath, 'ca_spike_data.mat');
        save(outFile, 'roi_spikes', 'spike_params', 'ifSpike', '-v7.3');

        % Update master file
        updateMasterFile(roi_spikes, spike_params);

        [~, folderName] = fileparts(S.folderPath);
        lblCurrentFolder.Text = [folderName '  (saved)'];
        updateNavDisplay();

        % Save composite video for the current ROI
        vidFile = saveROIVideo();
        fullVidFile = saveFullFrameVideo();

        msg = sprintf('Saved %d ROIs to:\n%s', S.nROI, outFile);
        if ~isempty(vidFile),     msg = [msg sprintf('\n\nCrop video:\n%s',      vidFile)];     end
        if ~isempty(fullVidFile), msg = [msg sprintf('\n\nFull-frame video:\n%s', fullVidFile)]; end
        uialert(fig, msg, 'Saved','Icon','success');
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

        % Load or initialise master struct
        master = struct();
        master.sessions = [];   % plain [] avoids the 0x0 struct bug

        if isfile(masterFile)
            try
                loaded = load(masterFile, 'master');
                if isstruct(loaded.master) && isfield(loaded.master, 'sessions')
                    master = loaded.master;
                end
            catch
                % keep fresh master on corruption
            end
        end

        % Upsert: find existing entry for this folder or append
        nSess = numel(master.sessions);
        found = false;
        for ii = 1:nSess
            if strcmp(master.sessions(ii).folder, S.folderPath)
                master.sessions(ii) = entry;
                found = true;
                break;
            end
        end
        if ~found
            if nSess == 0
                master.sessions = entry;           % first session
            else
                master.sessions(end+1) = entry;    % append
            end
        end

        save(masterFile, 'master', '-v7.3');
    end

    % -- Figure close: stop timer before deleting figure -------------------
    function cb_close(~,~)
        stopVideoTimer();
        delete(fig);
    end

%% =========================================================================
%% PLOT HELPERS
%% =========================================================================

    function refreshAll()
        refreshCrop();
        refreshTrace();
    end

    function refreshCrop()
        % Stop any running video playback
        stopVideoTimer();

        % Clear video state and cursor
        S.vidRawFrames = [];
        S.vidBndMask   = [];
        S.hVidImg      = [];
        S.hTraceCursor = [];
        S.vidFrIdx     = 1;

        cla(axSum);
        cla(axVideo);

        if S.nROI == 0, return; end

        [~, sessName] = fileparts(S.folderPath);
        roiIdx = S.curROI;

        % Update ROI label
        lblROI.Text = sprintf('ROI: %d / %d', roiIdx, S.nROI);

        if ~S.hasMask || isempty(S.maskL)
            axSum.Title.String = sprintf('ROI %d  --  %s  (no mask)', ...
                roiIdx, strrep(sessName,'_',' '));
            axSum.Title.Interpreter = 'none';
            axVideo.Title.String = 'Video';
            axVideo.Title.Interpreter = 'none';
            return;
        end

        % Get ROI centroid from mask
        labels = unique(S.maskL(:));
        labels(labels == 0) = [];
        if roiIdx > numel(labels)
            axSum.Title.String = sprintf('ROI %d  --  %s  (out of range)', ...
                roiIdx, strrep(sessName,'_',' '));
            axSum.Title.Interpreter = 'none';
            axVideo.Title.String = 'Video';
            axVideo.Title.Interpreter = 'none';
            return;
        end

        roiMask = (S.maskL == labels(roiIdx));
        props = regionprops(roiMask, 'Centroid');
        if isempty(props)
            axSum.Title.String = sprintf('ROI %d  --  %s  (empty mask)', ...
                roiIdx, strrep(sessName,'_',' '));
            axSum.Title.Interpreter = 'none';
            axVideo.Title.String = 'Video';
            axVideo.Title.Interpreter = 'none';
            return;
        end

        cx = round(props(1).Centroid(1));
        cy = round(props(1).Centroid(2));

        % Choose image source for avg projection
        if ~isempty(S.sumImg)
            baseImg = S.sumImg;
        elseif ~isempty(S.meanImg)
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

        % Compute boundary mask (logical) for the crop region
        bndMask = false(size(cropMask));
        if any(cropMask(:))
            bnd = bwboundaries(cropMask);
            for b = 1:numel(bnd)
                for bp = 1:size(bnd{b},1)
                    br = bnd{b}(bp,1);
                    bc = bnd{b}(bp,2);
                    if br >= 1 && br <= size(bndMask,1) && bc >= 1 && bc <= size(bndMask,2)
                        bndMask(br, bc) = true;
                    end
                end
            end
        end

        % === Left panel: Avg Projection ===
        imagesc(axSum, cropImg);
        colormap(axSum, 'gray');
        axis(axSum, 'image');
        hold(axSum, 'on');
        if any(bndMask(:))
            bnd = bwboundaries(cropMask);
            for b = 1:numel(bnd)
                plot(axSum, bnd{b}(:,2), bnd{b}(:,1), 'r-', 'LineWidth', 1.5);
            end
        end
        hold(axSum, 'off');
        axSum.Title.String = sprintf('Avg Proj  |  ROI %d', roiIdx);
        axSum.Title.Interpreter = 'none';

        % === Right panel: Video ===
        if S.hasTif && ~isempty(S.tifInfo)
            axVideo.Title.String = 'Video (loading...)';
            axVideo.Title.Interpreter = 'none';
            drawnow;

            startFr = S.nDropFrames + 1;
            nTotal  = numel(S.tifInfo);
            if startFr > nTotal, startFr = 1; end
            nAvail  = nTotal - startFr + 1;
            nVid    = min(500, nAvail);
            vidIdxs = startFr - 1 + round(linspace(1, nAvail, nVid));
            nVid    = numel(vidIdxs);

            % Load raw grayscale crops
            cropH = r2 - r1 + 1;
            cropW = c2 - c1 + 1;
            rawFrames = zeros(cropH, cropW, nVid, 'double');
            for ki = 1:nVid
                try
                    fr = double(imread(S.tifFile, vidIdxs(ki), 'Info', S.tifInfo));
                    rawFrames(:,:,ki) = fr(r1:r2, c1:c2);
                catch
                    % If frame read fails, leave as zeros
                end
            end

            % Compute normalization range (1st-99th percentile of first 10 frames)
            calibData = rawFrames(:,:,1:min(10,nVid));
            vmin = prctile(calibData(:), 1);
            vmax = prctile(calibData(:), 99);
            if vmax <= vmin, vmax = vmin + 1; end

            % Store raw frames and boundary mask for contrast sliders
            S.vidRawFrames = rawFrames;
            S.vidBndMask   = bndMask;
            S.vidVmin      = vmin;
            S.vidVmax      = vmax;

            % Set slider ranges and values
            globalMin = min(rawFrames(:));
            globalMax = max(rawFrames(:));
            if globalMax == globalMin, globalMax = globalMin + 1; end

            % Sum sliders
            sldSumMin.Limits = [globalMin globalMax];
            sldSumMax.Limits = [globalMin globalMax];
            sumVmin = prctile(cropImg(:), 1);
            sumVmax = prctile(cropImg(:), 99);
            if sumVmax <= sumVmin, sumVmax = sumVmin + 1; end
            sldSumMin.Value = max(globalMin, min(globalMax - eps, sumVmin));
            sldSumMax.Value = min(globalMax, max(globalMin + eps, sumVmax));

            % Video sliders
            sldVidMin.Limits = [globalMin globalMax];
            sldVidMax.Limits = [globalMin globalMax];
            sldVidMin.Value = max(globalMin, min(globalMax - eps, vmin));
            sldVidMax.Value = min(globalMax, max(globalMin + eps, vmax));

            % Apply clim to axSum
            clim(axSum, [sldSumMin.Value, sldSumMax.Value]);

            % Display first frame as RGB with boundary overlay
            S.hVidImg = image(axVideo, applyContrast(rawFrames(:,:,1), vmin, vmax, bndMask));
            axis(axVideo, 'image');

            % Set scrub slider range
            sldTime.Limits = [1 max(2, nVid)];
            sldTime.Value  = 1;
            S.vidFrIdx     = 1;

            axVideo.Title.String = sprintf('Video  |  ROI %d  (%d frames)', roiIdx, nVid);
            axVideo.Title.Interpreter = 'none';

            % Start playback timer
            startVideoTimer();
        else
            axVideo.Title.String = 'Video (no TIF)';
            axVideo.Title.Interpreter = 'none';
        end
    end

    % -- Contrast helper: raw frame -> uint8 RGB with red boundary ---------
    function rgb = applyContrast(rawFrame, vmin, vmax, bndMask)
        fr8 = uint8(min(255, max(0, (rawFrame - vmin) / (vmax - vmin) * 255)));
        R = fr8; G = fr8; B = fr8;
        R(bndMask) = 255; G(bndMask) = 0; B(bndMask) = 0;
        rgb = cat(3, R, G, B);
    end

    % -- Contrast slider callbacks -----------------------------------------
    function cb_sumContrastChanged()
        lo = sldSumMin.Value;  hi = sldSumMax.Value;
        if lo >= hi
            % clamp: keep a small gap
            if lo >= sldSumMax.Limits(2) - eps
                lo = sldSumMax.Limits(2) - (sldSumMax.Limits(2) - sldSumMax.Limits(1))*0.01;
            end
            hi = max(lo + eps, hi);
        end
        if isgraphics(axSum)
            clim(axSum, [lo hi]);
        end
    end

    function cb_vidContrastChanged()
        lo = sldVidMin.Value;  hi = sldVidMax.Value;
        if lo >= hi
            hi = max(lo + eps, hi);
        end
        S.vidVmin = lo;  S.vidVmax = hi;
        % Immediately re-render current frame
        if ~isempty(S.vidRawFrames) && S.vidFrIdx >= 1 && S.vidFrIdx <= size(S.vidRawFrames,3)
            fr = S.vidRawFrames(:,:,S.vidFrIdx);
            if isgraphics(S.hVidImg)
                S.hVidImg.CData = applyContrast(fr, lo, hi, S.vidBndMask);
            end
        end
    end

    % -- Time scrub slider callback ----------------------------------------
    function cb_scrubTime(~,~)
        if isempty(S.vidRawFrames), return; end
        nF = size(S.vidRawFrames, 3);
        newIdx = max(1, min(nF, round(sldTime.Value)));
        S.vidFrIdx = newIdx;
        fr = S.vidRawFrames(:,:,newIdx);
        if isgraphics(S.hVidImg)
            S.hVidImg.CData = applyContrast(fr, S.vidVmin, S.vidVmax, S.vidBndMask);
        end
        updateTraceCursor();
        drawnow limitrate;
    end

    % -- Video playback helpers --------------------------------------------
    function startVideoTimer()
        stopVideoTimer();
        if isempty(S.vidRawFrames) || size(S.vidRawFrames,3) < 2, return; end
        period = 1 / min(efFps.Value, 15);
        S.vidTimer = timer('ExecutionMode','fixedRate', ...
                           'Period', max(0.04, period), ...
                           'TimerFcn', @(~,~) stepVideoFrame());
        start(S.vidTimer);
    end

    function stopVideoTimer()
        if ~isempty(S.vidTimer) && isvalid(S.vidTimer)
            stop(S.vidTimer);
            delete(S.vidTimer);
        end
        S.vidTimer = [];
    end

    function stepVideoFrame()
        if isempty(S.vidRawFrames) || ~isgraphics(S.hVidImg), return; end
        nF = size(S.vidRawFrames, 3);
        S.vidFrIdx = mod(S.vidFrIdx, nF) + 1;

        % Update video display with current contrast
        fr = S.vidRawFrames(:,:,S.vidFrIdx);
        S.hVidImg.CData = applyContrast(fr, S.vidVmin, S.vidVmax, S.vidBndMask);

        % Update scrub slider position
        if isgraphics(sldTime) && sldTime.Limits(2) >= S.vidFrIdx
            sldTime.Value = S.vidFrIdx;
        end

        % Update blue trace cursor
        updateTraceCursor();
        drawnow limitrate;
    end

    % -- Blue cursor on axTrace showing current video frame time -----------
    function updateTraceCursor()
        if isempty(S.vidRawFrames) || isempty(S.t_img), return; end
        if ~isgraphics(axTrace), return; end

        % vidIdxs maps video frame ki -> TIF frame number -> t_img index
        nFrames  = numel(S.tifInfo);
        nAvail   = nFrames - S.nDropFrames;
        nVid     = size(S.vidRawFrames, 3);
        vidIdxs  = round(linspace(1, nAvail, nVid));   % indices into dFF (1-based)

        fi = min(S.vidFrIdx, numel(vidIdxs));
        dffIdx = vidIdxs(fi);
        dffIdx = max(1, min(length(S.t_img), dffIdx));
        t_cur  = S.t_img(dffIdx);

        if isgraphics(S.hTraceCursor)
            S.hTraceCursor.Value = t_cur;
        elseif isgraphics(axTrace)
            hold(axTrace, 'on');
            S.hTraceCursor = xline(axTrace, t_cur, 'b-', 'LineWidth', 1.5, ...
                'Label','frame', 'LabelHorizontalAlignment','right');
            S.hTraceCursor.HitTest = 'off';
        end
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

        % OASIS deconvolved overlay (orange) -- shown if data is available
        S.hOasisLine = gobjects(0);
        if ~isempty(S.dFFout) && isfield(S.dFFout,'dFF_oasis_deconv') ...
                && ~isempty(S.dFFout.dFF_oasis_deconv) ...
                && size(S.dFFout.dFF_oasis_deconv,2) >= roiIdx ...
                && length(S.dFFout.dFF_oasis_deconv(:,roiIdx)) == length(S.t_img)
            S.hOasisLine = plot(axTrace, S.t_img, S.dFFout.dFF_oasis_deconv(:,roiIdx), ...
                'Color',[1.0 0.55 0.10], 'LineWidth', 0.8);
            S.hOasisLine.HitTest = 'off';
        end

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

        % Reset and re-create trace cursor
        S.hTraceCursor = [];
        updateTraceCursor();
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

    % -- Save composite video for current ROI (crop video + trace w/ time cursor) --
    function vidFile = saveROIVideo()
        vidFile = '';
        if isempty(S.vidRawFrames) || isempty(S.dFF) || S.nROI == 0, return; end

        roiIdx  = S.curROI;
        nVid    = size(S.vidRawFrames, 3);
        nFrames = numel(S.tifInfo);
        nAvail  = nFrames - S.nDropFrames;
        vidIdxs = round(linspace(1, nAvail, nVid));

        trace = S.dFF(:, roiIdx);
        t     = S.t_img;
        spk   = S.spikes{roiIdx};

        % Panel dimensions
        panelH = 400;
        traceW = 800;

        % Scale crop to panelH (preserve aspect, force even width)
        cropH_src = size(S.vidRawFrames, 1);
        cropW_src = size(S.vidRawFrames, 2);
        sf    = panelH / cropH_src;
        cropH = panelH;
        cropW = max(2, 2 * round(cropW_src * sf / 2));

        outFile = fullfile(S.folderPath, sprintf('roi%d_ca_video.mp4', roiIdx));
        vw  = [];
        hFig = [];
        try
            vw = VideoWriter(outFile, 'MPEG-4'); %#ok<TNMLP>
            vw.FrameRate = min(efFps.Value, 15);
            open(vw);

            % Build trace figure (hidden) -- rendered once, time cursor updated per frame
            hFig = figure('Visible','off', 'Color','white', ...
                          'Units','pixels', 'Position',[0 0 traceW panelH], ...
                          'MenuBar','none', 'ToolBar','none');
            ax = axes(hFig, 'Units','normalized', 'Position',[0.10 0.15 0.86 0.76]);
            plot(ax, t, trace, 'Color',[0 0 0], 'LineWidth', 0.9);
            hold(ax, 'on');

            % Static spike markers
            if ~isempty(spk)
                valid = spk(spk >= 1 & spk <= length(trace));
                if ~isempty(valid)
                    scatter(ax, t(valid), trace(valid), 40, 'v', 'filled', ...
                        'MarkerFaceColor',[1.00 0.28 0.28], 'MarkerEdgeColor','none');
                end
            end

            yLims = ylim(ax);
            % Time cursor as a two-point vertical line
            hXL = plot(ax, [t(1) t(1)], yLims, 'r-', 'LineWidth', 2);
            xlim(ax, [t(1) t(end)]);  ylim(ax, yLims);
            xlabel(ax, 'Time (s)', 'FontSize', 9);
            ylabel(ax, 'dF/F',     'FontSize', 9);
            grid(ax, 'on');  box(ax, 'on');
            title(ax, sprintf('ROI %d', roiIdx), 'FontSize', 10, 'Interpreter','none');
            drawnow;

            for ki = 1:nVid
                % Left panel: render crop frame with current contrast settings
                cropRGB    = applyContrast(S.vidRawFrames(:,:,ki), ...
                    S.vidVmin, S.vidVmax, S.vidBndMask);
                cropScaled = imresize(cropRGB, [cropH cropW]);

                % Right panel: advance time cursor and capture trace
                dffIdx = min(vidIdxs(ki), length(t));
                t_cur = t(dffIdx);
                set(hXL, 'XData', [t_cur t_cur]);
                drawnow limitrate;
                fr = getframe(hFig);
                traceImg = fr.cdata;
                if size(traceImg,1) ~= panelH || size(traceImg,2) ~= traceW
                    traceImg = imresize(traceImg, [panelH traceW]);
                end

                writeVideo(vw, [cropScaled, traceImg]);
            end

            close(vw);   vw   = [];
            close(hFig); hFig = [];
            vidFile = outFile;
        catch ME
            if ~isempty(vw)  && isopen(vw),  close(vw);  end
            if ~isempty(hFig) && isgraphics(hFig), close(hFig); end
            warning('saveROIVideo: %s', ME.message);
        end
    end

    % -- Save full-frame video for current ROI -----------------------------
    function vidFile = saveFullFrameVideo()
        vidFile = '';
        if ~S.hasTif || isempty(S.tifInfo) || isempty(S.dFF) || S.nROI == 0, return; end
        if ~S.hasMask || isempty(S.maskL), return; end

        roiIdx = S.curROI;

        % ROI full-frame mask + boundary
        labels = unique(S.maskL(:)); labels(labels == 0) = [];
        if roiIdx > numel(labels), return; end
        roiMask  = (S.maskL == labels(roiIdx));
        bndSegs  = bwboundaries(roiMask);
        [H, W]   = size(roiMask);

        % Same frame sampling used by the GUI video
        nFrames  = numel(S.tifInfo);
        nAvail   = max(1, nFrames - S.nDropFrames);
        nVid     = min(500, nAvail);
        startFr  = S.nDropFrames + 1;
        vidIdxs  = startFr - 1 + round(linspace(1, nAvail, nVid));
        nVid     = numel(vidIdxs);

        trace = S.dFF(:, roiIdx);
        t     = S.t_img;
        spk   = S.spikes{roiIdx};

        % Panel sizes (scale full frame to max 512 px height)
        panelH = min(H, 512);
        panelW = max(2, 2*round(W * panelH / H / 2));   % keep aspect, even width
        traceW = 800;

        % Normalization from avg projection
        if ~isempty(S.sumImg)
            vmin = double(prctile(S.sumImg(:), 1));
            vmax = double(prctile(S.sumImg(:), 99));
        elseif ~isempty(S.meanImg)
            vmin = double(prctile(S.meanImg(:), 1));
            vmax = double(prctile(S.meanImg(:), 99));
        else
            vmin = 0; vmax = 1000;
        end
        if vmax <= vmin, vmax = vmin + 1; end

        outFile = fullfile(S.folderPath, sprintf('roi%d_fullframe_video.mp4', roiIdx));
        vw   = [];
        hFig2 = [];
        try
            vw = VideoWriter(outFile, 'MPEG-4'); %#ok<TNMLP>
            vw.FrameRate = min(efFps.Value, 15);
            open(vw);

            % Build hidden trace figure
            hFig2 = figure('Visible','off','Color','white','Units','pixels', ...
                           'Position',[0 0 traceW panelH],'MenuBar','none','ToolBar','none');
            ax2 = axes(hFig2,'Units','normalized','Position',[0.10 0.15 0.86 0.76]);
            plot(ax2, t, trace,'Color',[0 0 0],'LineWidth',0.9);
            hold(ax2,'on');
            if ~isempty(spk)
                valid = spk(spk >= 1 & spk <= length(trace));
                if ~isempty(valid)
                    scatter(ax2, t(valid), trace(valid), 40, 'v', 'filled', ...
                        'MarkerFaceColor',[1 0.28 0.28],'MarkerEdgeColor','none');
                end
            end
            yLims = ylim(ax2);
            hXL   = plot(ax2,[t(1) t(1)], yLims,'r-','LineWidth',2);
            xlim(ax2,[t(1) t(end)]); ylim(ax2, yLims);
            xlabel(ax2,'Time (s)','FontSize',9);
            ylabel(ax2,'dF/F','FontSize',9);
            grid(ax2,'on'); box(ax2,'on');
            title(ax2, sprintf('ROI %d  (full frame)', roiIdx),'FontSize',10,'Interpreter','none');
            drawnow;

            for ki = 1:nVid
                % Load full frame from TIF
                try
                    fr = double(imread(S.tifFile, vidIdxs(ki), 'Info', S.tifInfo));
                catch
                    fr = zeros(H, W);
                end

                % Normalise -> uint8 RGB
                fr8 = uint8(min(255, max(0, (fr - vmin)/(vmax - vmin)*255)));
                R = fr8; G = fr8; B = fr8;

                % Draw ROI boundary red on full frame
                for b = 1:numel(bndSegs)
                    pts = bndSegs{b};
                    for pi = 1:size(pts,1)
                        rr = pts(pi,1); cc = pts(pi,2);
                        if rr>=1&&rr<=H&&cc>=1&&cc<=W
                            R(rr,cc)=255; G(rr,cc)=0; B(rr,cc)=0;
                        end
                    end
                end
                fullRGB    = imresize(cat(3,R,G,B), [panelH panelW]);

                % Update red time cursor
                dffIdx = max(1, min(length(t), vidIdxs(ki) - S.nDropFrames));
                set(hXL,'XData',[t(dffIdx) t(dffIdx)]);
                drawnow limitrate;
                fr2      = getframe(hFig2);
                traceImg = fr2.cdata;
                if size(traceImg,1) ~= panelH || size(traceImg,2) ~= traceW
                    traceImg = imresize(traceImg, [panelH traceW]);
                end

                writeVideo(vw, [fullRGB, traceImg]);
            end

            close(vw);    vw    = [];
            close(hFig2); hFig2 = [];
            vidFile = outFile;
        catch ME
            if ~isempty(vw)    && isopen(vw),       close(vw);    end
            if ~isempty(hFig2) && isgraphics(hFig2), close(hFig2); end
            warning('saveFullFrameVideo: %s', ME.message);
        end
    end

end  % calcium_spike_gui


%% =========================================================================
%% LOCAL HELPER -- findpeaks-based spike finder
%% =========================================================================
function locs = find_spikes_threshold(trace, fps, threshold, minDistS, minWS, minProm)
%FIND_SPIKES_THRESHOLD  Detect calcium spikes using MATLAB findpeaks.
%   Uses MinPeakHeight, MinPeakDistance, MinPeakWidth, MinPeakProminence.

trace = trace(:);
locs  = [];
if length(trace) < 3, return; end

minDistFr = max(1, round(minDistS * fps));
minWFr    = max(0, round(minWS * fps));

pkArgs = {'MinPeakHeight',    threshold, ...
          'MinPeakDistance',  minDistFr};
if minWFr > 0
    pkArgs = [pkArgs, {'MinPeakWidth', minWFr}];
end
if minProm > 0
    pkArgs = [pkArgs, {'MinPeakProminence', minProm}];
end

[~, locs] = findpeaks(trace, pkArgs{:});
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


%% =========================================================================
%% LOCAL HELPER -- detect Python executable for cellpose-gpu environment
%% =========================================================================
function pyExe = detectPyExe()
%DETECTPYEXE  Find cellpose-gpu Python executable under USERPROFILE.
    candidates = { ...
        fullfile(getenv('USERPROFILE'), '.conda',     'envs','cellpose-gpu','python.exe'), ...
        fullfile(getenv('USERPROFILE'), 'anaconda3',  'envs','cellpose-gpu','python.exe'), ...
        fullfile(getenv('USERPROFILE'), 'miniconda3', 'envs','cellpose-gpu','python.exe') };
    pyExe = '';
    for i = 1:numel(candidates)
        if isfile(candidates{i}), pyExe = candidates{i}; return; end
    end
end
