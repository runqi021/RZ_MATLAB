function calcium_spike_gui()
%CALCIUM_SPIKE_GUI  Interactive GUI for calcium spike detection with per-ROI inspection
%
%   Scans a master folder for *_dFF.mat files (output of helper.dFF_RZ),
%   pairs them with *_cpSAM_output.mat (ROI masks), and lets you detect /
%   manually edit calcium spikes on each ROI individually.
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
S.meanImg     = [];      % [H x W] mean projection
S.hasMask     = false;
S.samFile     = '';
S.dffFile     = '';
S.folderPath  = '';
S.editMode    = 'none';  % 'none' | 'add' | 'delete'
S.masterFolder = '';
S.sessionList  = {};     % cell array of structs: .folder .dffFile
S.sessionIdx   = 0;
S.dFFout       = [];     % full dFFout struct (may contain spikes_oasis)

S.sumImg        = [];    % [H x W] mean projection (from TIF)
S.F_raw         = [];    % [T x N] raw F from SAM (already mean-per-pixel, all frames)
S.hasTif        = false;
S.tifFile       = '';
S.dFFraw        = [];    % copy of dFF as loaded (before optional detrend)
S.detrend       = false; % is detrend currently applied?
S.trimInfo      = [];    % struct with trim metadata (frame_idx, trimStart_sec, etc.)
S.hOasisLine    = [];    % orange overlay line for OASIS deconvolved trace
S.hRawFLine     = [];    % line handle for raw F panel
S.showOasis     = false; % toggle: overlay OASIS dFF on trace plot


% -- Factory defaults (restored by Reset button) --------------------------
DEFAULTS.BaselineWinSec  = 20;
DEFAULTS.fps        = 30;
DEFAULTS.method     = 'Raw dFF';
DEFAULTS.threshold  = 2.0;
DEFAULTS.minDist    = 0.5;
DEFAULTS.minWidth   = 0.1;
DEFAULTS.snapWin    = 0.10;
DEFAULTS.diameter   = 40;
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

nCtrlRows = 34;
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

r=r+1; btnDetect = uibutton(cGL,'Text','Detect Spikes (Current ROI)', ...
           'ButtonPushedFcn',@cb_detect, ...
           'BackgroundColor',[0.18 0.65 0.30],'FontColor','white');
       btnDetect.Layout.Row=r; btnDetect.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Trim Start (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efTrimStart = uieditfield(cGL,'numeric','Value',0,'Limits',[0 Inf]);
       efTrimStart.Layout.Row=r; efTrimStart.Layout.Column=2;
       efTrimStart.ValueChangedFcn = @(~,~) refreshRawF();

r=r+1; h=uilabel(cGL,'Text','Trim End (s):');
       h.Layout.Row=r; h.Layout.Column=1;
       efTrimEnd = uieditfield(cGL,'numeric','Value',Inf,'Limits',[0 Inf]);
       efTrimEnd.Layout.Row=r; efTrimEnd.Layout.Column=2;
       efTrimEnd.ValueChangedFcn = @(~,~) refreshRawF();

r=r+1; btnApplyTrim = uibutton(cGL,'Text','Apply Trim (crop data)', ...
           'ButtonPushedFcn',@cb_applyTrim, ...
           'BackgroundColor',[0.85 0.55 0.05],'FontColor','white');
       btnApplyTrim.Layout.Row=r; btnApplyTrim.Layout.Column=[1 2];

r=r+1; btnPlayMovie = uibutton(cGL,'Text','Play ROI Movie', ...
           'ButtonPushedFcn',@cb_playMovie, ...
           'BackgroundColor',[0.15 0.45 0.70],'FontColor','white');
       btnPlayMovie.Layout.Row=r; btnPlayMovie.Layout.Column=[1 2];

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
rightGL = uigridlayout(rootGL, [3 1]);
rightGL.Layout.Column = 2;
rightGL.RowHeight     = {'0.22x', '0.42x', '0.36x'};
rightGL.Padding       = [0 0 0 0];
rightGL.RowSpacing    = 6;

% Row 1: Avg Projection
axSum = uiaxes(rightGL);
axSum.Layout.Row             = 1;
axSum.Title.String           = 'Avg Projection';
axSum.Title.Interpreter      = 'none';
axSum.XLabel.String          = '';
axSum.YLabel.String          = '';
axis(axSum, 'image');

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

% Row 3: Raw F panel
axRawF = uiaxes(rightGL);
axRawF.Layout.Row            = 3;
axRawF.Title.String          = 'Raw F';
axRawF.Title.Interpreter     = 'none';
axRawF.XLabel.String         = 'Time (s)';
axRawF.XLabel.Interpreter    = 'none';
axRawF.YLabel.String         = 'F (a.u.)';
axRawF.YLabel.Interpreter    = 'none';
hold(axRawF,'on'); grid(axRawF,'on'); box(axRawF,'on');

% Persistent overlay handles (updated in-place)
hTraceLine = gobjects(0);   % dF/F trace line
hSpikeSc   = gobjects(0);   % spike scatter markers
hThreshLn  = gobjects(0);   % threshold horizontal line

% Keyboard shortcut on figure
fig.KeyPressFcn = @cb_key;

% Close callback
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
            uialert(fig,'No *_dFF.mat or *_cpSAM_output.mat files found in any subfolder.','No Sessions');
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
        S.hasMask    = ~isempty(sess.samFile) && isfile(sess.samFile);
        S.maskL      = [];
        S.meanImg    = [];
        S.sumImg     = [];
        S.dFFout     = [];
        S.showOasis  = false;

        % Show current folder name
        [~, folderName] = fileparts(S.folderPath);
        lblCurrentFolder.Text = folderName;

        % Load dFF data — or compute from SAM F if *_dFF.mat is missing
        if ~isempty(S.dffFile) && isfile(S.dffFile)
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
        elseif S.hasMask
            % No *_dFF.mat but SAM exists — compute dFF from F
            try
                samData = load(S.samFile);
            catch ME
                uialert(fig, ME.message, 'Load SAM Error'); return;
            end
            if ~isfield(samData, 'F')
                uialert(fig, 'SAM file has no F variable.', 'Load Error'); return;
            end

            fps_now    = efFps.Value;
            tossFrames = 30;               % pipeline default: TossFrames = 30
            F_roi_raw  = double(samData.F);
            F_roi      = F_roi_raw;
            if tossFrames > 0 && size(F_roi,1) > tossFrames
                F_roi(1:tossFrames,:) = [];
            end

            fprintf('Computing dFF from SAM F (%d frames, %d ROIs, toss=%d)...\n', ...
                size(F_roi,1), size(F_roi,2), tossFrames);
            dFFout_comp = dFF_RZ_dispatch(F_roi);
            S.dFF = dFFout_comp.dFF;

            % Build params struct matching Batch_BG_MC_SAM save format
            params = struct();
            params.sam_mat    = S.samFile;
            params.FPS        = fps_now;
            params.tossFrames = tossFrames;

            % Derive output filename from SAM filename stem
            [~, samStem] = fileparts(S.samFile);
            samStem = regexprep(samStem, '_cpSAM_output$', '');
            outDffMat = fullfile(S.folderPath, [samStem '_dFF.mat']);

            dFF     = S.dFF;          %#ok<NASGU>  % for save
            dFFout  = dFFout_comp;     %#ok<NASGU>
            save(outDffMat, 'dFF','dFFout','F_roi_raw','F_roi','params', '-v7.3');
            fprintf('Saved computed dFF: %s\n', outDffMat);

            % Update session entry so subsequent loads find the file
            S.dffFile = outDffMat;
            S.sessionList{S.sessionIdx}.dffFile = outDffMat;

            % Set loaded struct for downstream use
            loaded = struct('dFF', S.dFF, 'dFFout', dFFout_comp, 'params', params);
        else
            uialert(fig, 'No *_dFF.mat and no *_cpSAM_output.mat found.', 'Load Error');
            return;
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

        % Load raw F from SAM file; compute pixel count per ROI from maskL
        S.F_raw         = [];
        S.roiPixelCount = [];
        S.sumImg        = [];
        S.meanImg       = [];
        S.hasTif        = ~isempty(sess.tifFile) && isfile(sess.tifFile);
        S.tifFile       = sess.tifFile;

        if S.hasMask
            try
                samAll = load(S.samFile);
                if isfield(samAll, 'F')
                    % F is already mean-per-pixel (cp_sam_extract_F_cli divides by npix)
                    % Store all frames; frame alignment handled via dFFout.F_dff below
                    S.F_raw = double(samAll.F);
                end
            catch
            end
        end

        % Prefer the cellpose input AVG image if available
        avgHits = dir(fullfile(S.folderPath, '*_AVG_for_CP.tif'));
        if ~isempty(avgHits)
            try
                S.sumImg  = double(imread(fullfile(S.folderPath, avgHits(1).name)));
                S.meanImg = S.sumImg;
            catch
            end
        end

        % Fallback: mean projection from TIF (sample 50 frames, compute average)
        if isempty(S.sumImg) && S.hasTif
            try
                tifInfo = imfinfo(S.tifFile);
                nTotal  = numel(tifInfo);
                nDrop   = max(0, nTotal - size(S.dFF, 1));
                startFr = nDrop + 1;
                nAvail  = nTotal - startFr + 1;
                nSamp   = min(50, nAvail);
                sampIdx = unique(round(linspace(startFr, nTotal, nSamp)));
                acc = zeros(tifInfo(1).Height, tifInfo(1).Width);
                for ki = 1:numel(sampIdx)
                    acc = acc + double(imread(S.tifFile, sampIdx(ki), 'Info', tifInfo));
                end
                S.sumImg  = acc / numel(sampIdx);
                S.meanImg = S.sumImg;
            catch
            end
        end

        % Fallback projection: mean of F_raw per ROI painted onto mask (F already per-pixel)
        if isempty(S.sumImg) && ~isempty(S.F_raw) && S.hasMask && ~isempty(S.maskL)
            labels_fb = unique(S.maskL(:)); labels_fb(labels_fb == 0) = [];
            projImg   = zeros(size(S.maskL), 'double');
            nR_fb = min(numel(labels_fb), size(S.F_raw, 2));
            for k = 1:nR_fb
                projImg(S.maskL == labels_fb(k)) = mean(S.F_raw(:, k));
            end
            S.sumImg  = projImg;
            S.meanImg = projImg;
        end

        % Store raw dFF and reset detrend
        S.dFFraw  = S.dFF;
        S.detrend = false;
        btnDetrend.Text = 'Detrend: OFF';
        btnDetrend.BackgroundColor = [0.25 0.25 0.55];

        % Reset trim per FOV
        efTrimStart.Value = 0;
        efTrimEnd.Value   = Inf;
        S.trimInfo = [];

        % Reset OASIS state

        % Load previously saved OASIS results if available
        oasisFile = fullfile(S.folderPath, 'oasis.mat');
        if isfile(oasisFile)
            try
                oas = load(oasisFile);
                if isfield(oas,'dFF_oasis_deconv')
                    S.dFFout.F_oasis_deconv   = oas.F_oasis_deconv;
                    S.dFFout.dFF_oasis_deconv = oas.dFF_oasis_deconv;
                    S.dFFout.spikes_oasis     = oas.spikes_oasis;
                    S.dFFout.baseline_oasis   = oas.baseline_oasis;
                    if isfield(oas,'g_AR1'), S.dFFout.g_AR1 = oas.g_AR1; end
                    % OASIS data loaded from file
                end
            catch
            end
        end

        % Initialize spike storage
        S.spikes = cell(N, 1);
        S.curROI = 1;

        % Update ROI spinner limits
        % ROI range: 1..N

        % Check for previously saved results
        matFile = fullfile(S.folderPath, 'ca_spike_data.mat');
        if isfile(matFile)
            try
                saved = load(matFile, 'roi_spikes', 'spike_params', 'trim_info');
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
                    % method, snapWin, tauDecay, pythonExe: use hardcoded defaults
                    if isfield(sp,'diameter'),     efDiameter.Value  = sp.diameter;    end
                    if isfield(sp,'zoomFactor'),   efZoom.Value      = sp.zoomFactor;  end
                    if isfield(sp,'prominence'),   efProminence.Value = sp.prominence; end
                    if isfield(sp,'tauDecay'),     DEFAULTS.tauDecay  = sp.tauDecay;    end
                    if isfield(sp,'pythonExe'),    DEFAULTS.pythonExe = sp.pythonExe;   end
                end
                % Re-apply saved trim to freshly loaded full-length data
                if isfield(saved,'trim_info') && isstruct(saved.trim_info) ...
                        && isfield(saved.trim_info,'frame_idx')
                    ti = saved.trim_info;
                    idx = ti.frame_idx(:);
                    idx = idx(idx >= 1 & idx <= size(S.dFF,1));
                    if ~isempty(idx)
                        S.dFF    = S.dFF(idx, :);
                        S.dFFraw = S.dFF;
                        S.t_img  = (0:size(S.dFF,1)-1)' / S.fps;
                        if ~isempty(S.F_raw) && size(S.F_raw,1) >= max(idx)
                            S.F_raw = S.F_raw(idx, :);
                        end
                        if ~isempty(S.dFFout) && isstruct(S.dFFout)
                            oF = {'dFF_oasis_deconv','F_oasis_deconv', ...
                                  'spikes_oasis','baseline_oasis','F_dff'};
                            for ff = 1:numel(oF)
                                fn = oF{ff};
                                if isfield(S.dFFout,fn) && ~isempty(S.dFFout.(fn)) ...
                                        && size(S.dFFout.(fn),1) >= max(idx)
                                    S.dFFout.(fn) = S.dFFout.(fn)(idx, :);
                                end
                            end
                        end
                        S.trimInfo = ti;
                    end
                end
                lblCurrentFolder.Text = [folderName '  (saved)'];
            catch
                % Corrupted mat -- fall through
            end
        end

        % Auto-detect FPS and zoom from TIF metadata / _meta.mat
        % (runs AFTER saved params restore so detected values win)
        try
            [det_fps, det_meta] = detect_session_fps(S.folderPath, S.fps);
            S.fps = det_fps;  efFps.Value = det_fps;
            if isfield(det_meta,'zoomFactor') && ~isempty(det_meta.zoomFactor) ...
                    && isfinite(det_meta.zoomFactor) && det_meta.zoomFactor > 0
                efZoom.Value = det_meta.zoomFactor;
            end
        catch
        end

        % Fallback: parse zoom from folder name (e.g. '12x' in '...z0_12x_00001')
        if efZoom.Value == DEFAULTS.zoomFactor
            [~, fName] = fileparts(S.folderPath);
            tok = regexp(fName, '(\d+)x', 'tokens', 'once');
            if ~isempty(tok)
                zf = str2double(tok{1});
                if isfinite(zf) && zf >= 1 && zf <= 50
                    efZoom.Value = zf;
                end
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
        refreshAll();
    end

    function cb_nextROI(~,~)
        if S.curROI >= S.nROI, return; end
        S.curROI = S.curROI + 1;
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
        trace_det = S.dFF(:, roiIdx);

        spk = find_spikes_threshold(trace_det, fps, thr, minDist, minW, prom);

        % Filter spikes to trim window
        tStart = efTrimStart.Value;
        tEnd   = efTrimEnd.Value;
        if ~isempty(spk) && ~isempty(S.t_img)
            spk_t = S.t_img(min(spk, numel(S.t_img)));
            keep = spk_t >= tStart & spk_t <= tEnd;
            spk = spk(keep);
        end

        S.spikes{roiIdx} = spk(:);
    end

    % -- Run OASIS deconvolution on all ROIs -------------------------------
    function cb_runOasis(~,~)
        if isempty(S.dFF) || S.nROI == 0
            uialert(fig,'Load data first.','No Data'); return;
        end
        % Use dFFout.F_dff (frame-aligned, first 30 frames already tossed);
        % fall back to S.F_raw with manual tossing
        fps = efFps.Value;
        nToss = round(fps);  % toss first 30 frames (1 sec at 30 fps)
        if ~isempty(S.dFFout) && isfield(S.dFFout,'F_dff') && ~isempty(S.dFFout.F_dff)
            F_for_oasis = S.dFFout.F_dff;   % already frame-aligned
        elseif ~isempty(S.F_raw)
            if size(S.F_raw, 1) > nToss
                F_for_oasis = S.F_raw(nToss+1:end, :);
            else
                F_for_oasis = S.F_raw;
            end
        else
            uialert(fig, ...
                'No raw fluorescence available. Ensure *_cpSAM_output.mat is present.', ...
                'F Missing');
            return;
        end

        pyExe = DEFAULTS.pythonExe;
        if isempty(pyExe), pyExe = detectPyExe(); end
        if isempty(pyExe) || ~isfile(pyExe)
            uialert(fig, ...
                'Python executable not found. Set it in the Parameters section.', ...
                'Python Not Found');
            return;
        end

        tauDecay = DEFAULTS.tauDecay;
        %g = round(exp(-(1/fps)/tauDecay), 2);
        g = 0.93;
        drawnow;
        try
            [F_oasis, dFF_oasis, spikes_oasis, baseline_oasis] = ...
                helper.oasis_deconv_and_dff_AR1(F_for_oasis, g, 'PythonExe', pyExe);

            S.dFFout.F_oasis_deconv   = F_oasis;
            S.dFFout.dFF_oasis_deconv = dFF_oasis;
            S.dFFout.spikes_oasis     = spikes_oasis;
            S.dFFout.baseline_oasis   = baseline_oasis;
            S.dFFout.g_AR1            = g;

            % Save OASIS results to .mat
            oasis_save = struct();
            oasis_save.F_oasis_deconv   = F_oasis;
            oasis_save.dFF_oasis_deconv = dFF_oasis;
            oasis_save.spikes_oasis     = spikes_oasis;
            oasis_save.baseline_oasis   = baseline_oasis;
            oasis_save.g_AR1            = g;
            oasis_save.tauDecay         = tauDecay;
            oasis_save.fps              = fps;
            oasisPath = fullfile(S.folderPath, 'oasis.mat');
            save(oasisPath, '-struct', 'oasis_save', '-v7.3');
            fprintf('OASIS results saved: %s\n', oasisPath);

            S.showOasis = true;
            refreshTrace();
        catch ME
            uialert(fig, sprintf('OASIS failed:\n%s', ME.message), 'OASIS Error');
        end
    end

    % -- Toggle OASIS overlay on/off ----------------------------------------
    function cb_toggleOasisOverlay(~,~)
        if ~isfield(S.dFFout,'dFF_oasis_deconv') || isempty(S.dFFout.dFF_oasis_deconv)
            uialert(fig,'Run OASIS first.','No OASIS Data'); return;
        end
        S.showOasis = ~S.showOasis;
        refreshTrace();
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
                snapFr = max(1, round(DEFAULTS.snapWin * efFps.Value));
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
        efThresh.Value   = p.threshold;
        efMinDist.Value  = p.minDist;
        efMinWidth.Value = p.minWidth;
        efDiameter.Value = p.diameter;
        efZoom.Value     = p.zoomFactor;
        if isfield(p,'prominence'),  efProminence.Value = p.prominence;  end
    end

    function saveWorkingParams()
        wp.fps        = efFps.Value;
        wp.method     = 'Raw dFF';
        wp.threshold  = efThresh.Value;
        wp.minDist    = efMinDist.Value;
        wp.minWidth   = efMinWidth.Value;
        wp.snapWin    = DEFAULTS.snapWin;
        wp.diameter   = efDiameter.Value;
        wp.zoomFactor = efZoom.Value;
        wp.prominence = efProminence.Value;
        wp.tauDecay   = DEFAULTS.tauDecay;
        wp.pythonExe  = DEFAULTS.pythonExe;
        setpref('CalciumSpikeGUI','workingParams',wp);
    end

    % -- Apply Trim (hard crop) --------------------------------------------
    function cb_applyTrim(~,~)
        if isempty(S.dFF) || S.nROI == 0
            uialert(fig,'Load data first.','No Data'); return;
        end

        tStart = efTrimStart.Value;
        tEnd   = efTrimEnd.Value;
        T_orig = size(S.dFF, 1);

        % If no meaningful trim, do nothing
        if tStart <= S.t_img(1) && (isinf(tEnd) || tEnd >= S.t_img(end))
            uialert(fig,'Trim range covers the entire recording. Nothing to trim.','No Trim');
            return;
        end

        % Compute frame indices to keep
        keepMask = S.t_img >= tStart;
        if isfinite(tEnd)
            keepMask = keepMask & S.t_img <= tEnd;
        end
        keepIdx = find(keepMask);

        if isempty(keepIdx)
            uialert(fig,'Trim range is empty. No frames to keep.','Error'); return;
        end

        T_new = numel(keepIdx);

        % Track original frame indices (cumulative if trim already applied)
        if ~isempty(S.trimInfo) && isfield(S.trimInfo,'frame_idx')
            origIdx = S.trimInfo.frame_idx(keepIdx);
        else
            origIdx = keepIdx;
        end

        % Build old→new index map for spike re-indexing
        oldToNew = zeros(T_orig, 1);
        oldToNew(keepIdx) = 1:T_new;

        % Crop dFF and dFFraw
        S.dFF    = S.dFF(keepIdx, :);
        S.dFFraw = S.dFFraw(keepIdx, :);

        % Crop F_raw
        if ~isempty(S.F_raw) && size(S.F_raw,1) >= max(keepIdx)
            S.F_raw = S.F_raw(keepIdx, :);
        end

        % Crop OASIS / dFFout subfields
        if ~isempty(S.dFFout) && isstruct(S.dFFout)
            oFields = {'dFF_oasis_deconv','F_oasis_deconv','spikes_oasis', ...
                       'baseline_oasis','F_dff'};
            for ff = 1:numel(oFields)
                fn = oFields{ff};
                if isfield(S.dFFout, fn) && ~isempty(S.dFFout.(fn)) ...
                        && size(S.dFFout.(fn),1) >= max(keepIdx)
                    S.dFFout.(fn) = S.dFFout.(fn)(keepIdx, :);
                end
            end
        end

        % Reset time vector to start from 0
        S.t_img = (0:T_new-1)' / S.fps;

        % Re-index all spike indices
        for ii = 1:S.nROI
            spk = S.spikes{ii};
            if isempty(spk), continue; end
            spk = spk(spk >= 1 & spk <= T_orig);
            newSpk = oldToNew(spk);
            S.spikes{ii} = newSpk(newSpk > 0);
        end

        % Store trim metadata
        S.trimInfo = struct();
        S.trimInfo.frame_idx      = origIdx(:);
        S.trimInfo.trimStart_sec  = tStart;
        S.trimInfo.trimEnd_sec    = tEnd;
        S.trimInfo.orig_T         = T_orig;
        S.trimInfo.applied_at     = datestr(now); %#ok<TNOW1,DATST>

        % Reset trim fields (trim has been applied)
        efTrimStart.Value = 0;
        efTrimEnd.Value   = Inf;

        % Refresh display
        refreshTrace();
        refreshRawF();

        uialert(fig, sprintf('Trim applied: kept frames %d–%d (%d of %d).\nNew T = %d frames (%.1f s).', ...
            keepIdx(1), keepIdx(end), T_new, T_orig, T_new, T_new/S.fps), ...
            'Trim Applied', 'Icon','success');
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
            rs.ifSpike = rs.n_spikes > 0;
            roi_spikes(ii) = rs; %#ok<AGROW>
        end

        ifSpike = logical([roi_spikes.n_spikes] > 0);

        % Build spike_params struct
        spike_params = struct();
        spike_params.method     = 'Raw dFF';
        spike_params.threshold  = efThresh.Value;
        spike_params.minDist    = efMinDist.Value;
        spike_params.minWidth   = efMinWidth.Value;
        spike_params.snapWin    = DEFAULTS.snapWin;
        spike_params.fps        = fps;
        spike_params.diameter   = efDiameter.Value;
        spike_params.zoomFactor = efZoom.Value;
        spike_params.prominence = efProminence.Value;
        spike_params.tauDecay   = DEFAULTS.tauDecay;
        spike_params.dff_file   = S.dffFile;
        spike_params.saved_at   = datestr(now); %#ok<TNOW1,DATST>
        spike_params.oasis_run     = ~isempty(S.dFFout) && isfield(S.dFFout,'spikes_oasis') ...
                                      && ~isempty(S.dFFout.spikes_oasis);
        spike_params.oasis_tauDecay = DEFAULTS.tauDecay;
        spike_params.oasis_g        = '';
        if ~isempty(S.dFFout) && isfield(S.dFFout,'g_AR1')
            spike_params.oasis_g = S.dFFout.g_AR1;
        end
        spike_params.detect_on  = 'Raw dFF';
        spike_params.trimStart  = efTrimStart.Value;
        spike_params.trimEnd    = efTrimEnd.Value;

        % Include trim info (empty if no trim applied)
        trim_info = S.trimInfo; %#ok<NASGU>

        outFile = fullfile(S.folderPath, 'ca_spike_data.mat');
        save(outFile, 'roi_spikes', 'spike_params', 'ifSpike', 'trim_info', '-v7.3');
        updateMasterFile(roi_spikes, spike_params);

        [~, folderName] = fileparts(S.folderPath);
        lblCurrentFolder.Text = [folderName '  (saved)'];
        updateNavDisplay();

        uialert(fig, sprintf('Saved %d ROIs to:\n%s', S.nROI, outFile), 'Saved', 'Icon','success');
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
        entry.trim_info    = S.trimInfo;

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

        % Backfill trim_info on old sessions that lack the field
        if ~isempty(master.sessions) && ~isfield(master.sessions, 'trim_info')
            [master.sessions.trim_info] = deal([]);
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

    % -- Figure close -------------------------------------------------------
    function cb_close(~,~)
        delete(fig);
    end

%% =========================================================================
%% PLOT HELPERS
%% =========================================================================

    function refreshAll()
        refreshCrop();
        refreshTrace();
        refreshRawF();
    end

    function refreshCrop()
        cla(axSum);
        axSum.Color = [0 0 0];
        axSum.XLim = [0 1]; axSum.YLim = [0 1];

        if S.nROI == 0, return; end

        [~, sessName] = fileparts(S.folderPath);
        roiIdx = S.curROI;

        % Update ROI label
        lblROI.Text = sprintf('ROI: %d / %d', roiIdx, S.nROI);

        if ~S.hasMask || isempty(S.maskL)
            axSum.Title.String = sprintf('ROI %d  --  %s  (no mask)', ...
                roiIdx, strrep(sessName,'_',' '));
            axSum.Title.Interpreter = 'none';
            return;
        end

        % Get ROI centroid from mask
        labels = unique(S.maskL(:));
        labels(labels == 0) = [];
        if roiIdx > numel(labels)
            axSum.Title.String = sprintf('ROI %d  --  %s  (out of range)', ...
                roiIdx, strrep(sessName,'_',' '));
            axSum.Title.Interpreter = 'none';
            return;
        end

        roiMask = (S.maskL == labels(roiIdx));
        props = regionprops(roiMask, 'Centroid');
        if isempty(props)
            axSum.Title.String = sprintf('ROI %d  --  %s  (empty mask)', ...
                roiIdx, strrep(sessName,'_',' '));
            axSum.Title.Interpreter = 'none';
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
        um_per_px = 1.7778 / efZoom.Value;          % ScanImage default
        halfSide  = round(30 / um_per_px / 2);       % 30 um crop box
        boxSz     = 2 * halfSide + 1;                % guaranteed odd
        r1 = max(1, cy - halfSide);
        r2 = min(H, cy + halfSide);
        c1 = max(1, cx - halfSide);
        c2 = min(W, cx + halfSide);

        cropImg  = zeros(boxSz, boxSz);
        cropMask = false(boxSz, boxSz);
        pr1 = halfSide - (cy - r1) + 1;  % paste offset in padded box
        pc1 = halfSide - (cx - c1) + 1;
        cropImg(pr1:pr1+(r2-r1), pc1:pc1+(c2-c1))  = baseImg(r1:r2, c1:c2);
        cropMask(pr1:pr1+(r2-r1), pc1:pc1+(c2-c1)) = roiMask(r1:r2, c1:c2);

        % === Summary-style display: percentile clip + gamma + burned-in ROI outline ===
        cropImg = double(cropImg);
        validPx = cropImg(cropImg > 0);
        if isempty(validPx), validPx = cropImg(:); end
        lo = prctile(validPx, 0.5);
        hi = prctile(validPx, 99.5);
        if hi <= lo, hi = lo + 1; end
        cropImg = (cropImg - lo) / (hi - lo);
        cropImg = max(0, min(1, cropImg));
        cropImg = cropImg .^ 0.6;

        % Build RGB and burn yellow ROI perimeter
        patchRGB = repmat(cropImg, [1 1 3]);
        roi_perim = bwperim(cropMask);
        yellow = [1 1 0];
        for ch = 1:3
            plane = patchRGB(:,:,ch);
            plane(roi_perim) = yellow(ch);
            patchRGB(:,:,ch) = plane;
        end

        image(axSum, patchRGB);
        axis(axSum, 'image');
        set(axSum, 'XTick', [], 'YTick', []);
        axSum.Title.String = sprintf('Avg Proj  |  ROI %d', roiIdx);
        axSum.Title.Interpreter = 'none';

    end

    function cb_playMovie(~,~)
        % Open per-ROI MP4 in system video player
        if isempty(S.folderPath) || S.nROI == 0, return; end
        movieDir = fullfile(S.folderPath, 'Fmovie_perROI');
        vidPath = fullfile(movieDir, sprintf('ROI_%02d_movie.mp4', S.curROI));
        if isfile(vidPath)
            if ispc
                winopen(vidPath);
            else
                system(sprintf('open "%s" &', vidPath));
            end
        else
            uialert(fig, sprintf('Video not found:\n%s\nRun Batch_dffQC with MakeMontageVideo=true first.', vidPath), 'No Movie');
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

        % OASIS deconvolved overlay (orange) -- shown only when toggle is ON
        S.hOasisLine = gobjects(0);
        if S.showOasis && ~isempty(S.dFFout) ...
                && isfield(S.dFFout,'dFF_oasis_deconv') ...
                && ~isempty(S.dFFout.dFF_oasis_deconv) ...
                && size(S.dFFout.dFF_oasis_deconv,2) >= roiIdx
            oasisTrace = S.dFFout.dFF_oasis_deconv(:, roiIdx);
            nOasis = numel(oasisTrace);
            nT     = numel(S.t_img);
            % Use matching time vector (truncate to shorter length)
            nPlot = min(nOasis, nT);
            S.hOasisLine = plot(axTrace, S.t_img(1:nPlot), oasisTrace(1:nPlot), ...
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
        ylo = min(trace); yhi = max(trace);
        if yhi <= ylo, yhi = ylo + 1; end
        pad = (yhi - ylo) / 0.6 - (yhi - ylo);  % 20% padding each side
        ylim(axTrace, [ylo - pad/2, yhi + pad/2]);
        [~, sessName] = fileparts(S.folderPath);
        axTrace.Title.String = sprintf('ROI %d / %d  --  %s', ...
            roiIdx, S.nROI, strrep(sessName,'_',' '));
        axTrace.Title.Interpreter = 'none';

        redraw_spikes();
        updateSpikeCount();

        % Also update the Raw F panel when trace refreshes
        refreshRawF();
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

    function refreshRawF()
        cla(axRawF); hold(axRawF,'on');
        S.hRawFLine = [];

        if S.nROI == 0
            hold(axRawF,'off'); return;
        end

        % Prefer dFFout.F_dff — already frame-aligned and per-pixel (from dFF_RZ)
        % Fall back to S.F_raw (all frames, also per-pixel from cp_sam_extract_F_cli)
        if ~isempty(S.dFFout) && isfield(S.dFFout,'F_dff') ...
                && ~isempty(S.dFFout.F_dff) && S.curROI <= size(S.dFFout.F_dff, 2)
            rawF = S.dFFout.F_dff(:, S.curROI);
            t_f  = S.t_img(1:min(numel(rawF), numel(S.t_img)));
        elseif ~isempty(S.F_raw) && S.curROI <= size(S.F_raw, 2)
            rawF = S.F_raw(:, S.curROI);   % already per-pixel — no further division
            t_f  = (0:numel(rawF)-1)' / S.fps;
        else
            axRawF.Title.String = 'Raw F  (not available)';
            axRawF.Title.Interpreter = 'none';
            hold(axRawF,'off'); return;
        end

        S.hRawFLine = plot(axRawF, t_f, rawF(1:numel(t_f)), 'Color', [0 0 0], 'LineWidth', 0.9);
        S.hRawFLine.HitTest = 'off';

        xlim(axRawF, [t_f(1) t_f(end)]);
        fVals = rawF(1:numel(t_f));
        flo = min(fVals); fhi = max(fVals);
        if fhi <= flo, fhi = flo + 1; end
        fpad = (fhi - flo) / 0.6 - (fhi - flo);
        ylim(axRawF, [flo - fpad/2, fhi + fpad/2]);
        % Trim lines (yellow) on Raw F
        tStart = efTrimStart.Value;
        tEnd   = efTrimEnd.Value;
        if tStart > t_f(1)
            xline(axRawF, tStart, '-', 'Color', [0.9 0.8 0.0], 'LineWidth', 1.5, ...
                'Label', 'trim start', 'LabelHorizontalAlignment', 'right');
        end
        if isfinite(tEnd) && tEnd < t_f(end)
            xline(axRawF, tEnd, '-', 'Color', [0.9 0.8 0.0], 'LineWidth', 1.5, ...
                'Label', 'trim end', 'LabelHorizontalAlignment', 'left');
        end

        axRawF.Title.String = sprintf('Raw F  |  ROI %d', S.curROI);
        axRawF.Title.Interpreter = 'none';
        axRawF.YLabel.String = 'F (mean pixel intensity)';
        axRawF.YLabel.Interpreter = 'none';
        hold(axRawF,'off');
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
%SCANFORDFF SESSIONS  Recursively find all subfolders containing a *_dFF.mat
%   or *_cpSAM_output.mat (with F but no dFF yet).
%   Returns a cell array of structs with fields:
%       .folder   -- path to the subfolder
%       .dffFile  -- path to the *_dFF.mat file ('' if needs computing)
%       .samFile  -- path to the *_cpSAM_output.mat file (or '')
%       .tifFile  -- path to the *_minusDark_MC.tif file (or '')
%   Sorted alphabetically by folder path.

sessions = {};
folderSet = containers.Map('KeyType','char','ValueType','logical');

% 1) Recursive glob for dFF.mat files
hits = dir(fullfile(masterFolder, '**', '*_dFF.mat'));
if ~isempty(hits)
    folders  = {hits.folder};
    [uFolders, ia] = unique(folders, 'stable');
    for i = 1:numel(uFolders)
        h = hits(ia(i));
        folder  = h.folder;
        dffFile = fullfile(folder, h.name);

        samHits = dir(fullfile(folder, '*_cpSAM_output.mat'));
        samFile = '';
        if ~isempty(samHits), samFile = fullfile(folder, samHits(1).name); end

        tifFile = find_mc_tif(folder);

        sessions{end+1} = struct('folder',  folder, ...
                                 'dffFile', dffFile, ...
                                 'samFile', samFile, ...
                                 'tifFile', tifFile); %#ok<AGROW>
        folderSet(folder) = true;
    end
end

% 2) Also find folders with *_cpSAM_output.mat but NO *_dFF.mat
samHitsAll = dir(fullfile(masterFolder, '**', '*_cpSAM_output.mat'));
if ~isempty(samHitsAll)
    folders  = {samHitsAll.folder};
    [uFolders, ia] = unique(folders, 'stable');
    for i = 1:numel(uFolders)
        folder = uFolders{i};
        if folderSet.isKey(folder), continue; end   % already have a dFF entry

        % Verify the SAM file actually contains F
        samFile = fullfile(folder, samHitsAll(ia(i)).name);
        try
            w = whos('-file', samFile, 'F');
            if isempty(w), continue; end   % no F variable — skip
        catch
            continue;
        end

        tifFile = find_mc_tif(folder);

        sessions{end+1} = struct('folder',  folder, ...
                                 'dffFile', '', ...
                                 'samFile', samFile, ...
                                 'tifFile', tifFile); %#ok<AGROW>
        folderSet(folder) = true;
    end
end
end  % scanForDFFSessions


%% =========================================================================
%% LOCAL HELPER -- dFF computation (dispatch to helper or local fallback)
%% =========================================================================
function dFFout = dFF_RZ_dispatch(F)
%DFF_RZ_DISPATCH  Compute dF/F using helper.dFF_RZ if available, else local fallback.
%   F: [T x N] single/double ok
try
    dFFout = helper.dFF_RZ(F, 'BaselineWinSec', 20);
catch
    dFFout = dFF_RZ_local(F);
end
end

function out = dFF_RZ_local(F)
%DFF_RZ_LOCAL  Simple robust dF/F fallback.
%   F0 per ROI = 20th percentile over time (robust baseline)
%   dFF = (F - F0) ./ max(F0, eps)
F = double(F);
F0 = prctile(F, 20, 1);
F0(F0 <= 0) = eps;
dFF = (F - F0) ./ F0;

out = struct();
out.F0    = F0;
out.dFF   = dFF;
out.F_dff = F;   % frame-aligned raw F (after toss)
out.t_dff = [];  % caller should set from fps
end


%% =========================================================================
%% LOCAL HELPER -- detect Python executable for OASIS / cellpose-gpu environment
%% =========================================================================
function pyExe = detectPyExe()
%DETECTPYEXE  Find Python executable (cellpose-gpu has both cellpose + oasis).
    up = getenv('USERPROFILE');
    candidates = { ...
        fullfile(up, '.conda',     'envs','cellpose-gpu','python.exe'), ...
        fullfile(up, 'anaconda3',  'envs','cellpose-gpu','python.exe'), ...
        fullfile(up, 'miniconda3', 'envs','cellpose-gpu','python.exe'), ...
        fullfile(up, '.conda',     'envs','oasis','python.exe'), ...
        fullfile(up, 'anaconda3',  'envs','oasis','python.exe'), ...
        fullfile(up, 'miniconda3', 'envs','oasis','python.exe') };
    pyExe = '';
    for i = 1:numel(candidates)
        if isfile(candidates{i}), pyExe = candidates{i}; return; end
    end
end


%% =========================================================================
%% LOCAL HELPER -- find MC'd single-channel TIFF
%% =========================================================================
function tifFile = find_mc_tif(folder)
%FIND_MC_TIF  Locate the best MC TIFF for movie playback.
%   Priority: preproc_MC_MC > preproc_MC > minusDark_MC > any *MC*.tif
tifFile = '';
patterns = {"*_preproc_MC_MC.tif", "*_preproc_MC.tif", "*_minusDark_MC.tif"};
for ip = 1:numel(patterns)
    hits = dir(fullfile(folder, patterns{ip}));
    % Filter out AVG/QC/shifts/label TIFs
    hits = hits(~contains({hits.name}, {'_AVG','_QC','_shifts','_label','_SELROI'}, 'IgnoreCase', true));
    if ~isempty(hits)
        tifFile = fullfile(hits(end).folder, hits(end).name);
        return;
    end
end
end
