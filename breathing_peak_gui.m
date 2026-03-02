function breathing_peak_gui()
%BREATHING_PEAK_GUI  Interactive GUI for breathing inspiration peak detection
%
%   Reads DLC snapshot-best-20 CSV files, lets you select which dots
%   (dot1/dot2/dot3) and coordinate (x / y / magnitude) to build the
%   breathing signal from, detects peaks with findpeaks, and lets you
%   add/delete peaks interactively by clicking on the plot.
%
%   Usage:
%       breathing_peak_gui()
%
%   Output (saved to session folder):
%       breath_peak_data.mat  – insp_onset_idx, insp_onsets_t,
%                               breath, t_breath, findpeak_params
%
%   Requires: MATLAB R2021a+

%% ── Shared state (accessible by all nested callbacks) ───────────────────
S.rawData    = [];      % [F×10]: frameIdx, d1x d1y d1p, d2x d2y d2p, d3x d3y d3p
S.t          = [];      % [F×1] seconds
S.breath     = [];      % [F×1] processed (detrended, z-scored) signal
S.dotTraces  = [];      % [F×3] individual dot signals (normalised)
S.peaks        = [];    % [K×1] frame indices of accepted peaks
S.folderPath   = '';
S.csvFile      = '';
S.editMode     = 'none';  % 'none' | 'add' | 'delete'
S.masterFolder = '';
S.sessionList  = {};    % cell array of structs: .folder .csvFile
S.sessionIdx   = 0;

% ── Factory defaults (restored by Reset button) ───────────────────────────
DEFAULTS.fps          = 30;
DEFAULTS.dotsOn       = [true true true false];
DEFAULTS.coordTypes   = {'y','y','y','y'};
DEFAULTS.likelihood   = 0.5;
DEFAULTS.combine      = 'average';
DEFAULTS.invert       = false;
DEFAULTS.minProm      = 1.00;
DEFAULTS.minDist      = 0.20;
DEFAULTS.minHeight    = 0.50;
DEFAULTS.minWidth     = 0.05;
DEFAULTS.snapWin      = 0.10;

nDots      = 4;                                    % max supported
dotNames   = {'dot1','dot2','dot3','dot4'};        % updated on CSV load
dotColors  = {[0 0 0]; [0 0 0]; [0 0 0]; [0 0 0]};
% rawData column indices for each dot: [xCol yCol pCol]
dotCols    = [2 3 4; 5 6 7; 8 9 10; 11 12 13];

%% ── Figure ──────────────────────────────────────────────────────────────
fig = uifigure('Name','Breathing Peak Sorter', ...
               'Position',[40 30 1560 930], ...
               'WindowState','maximized');

rootGL = uigridlayout(fig, [1 2]);
rootGL.ColumnWidth   = {320, '1x'};
rootGL.Padding       = [6 6 6 6];
rootGL.ColumnSpacing = 8;

%% ── Left: control panel ─────────────────────────────────────────────────
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

% ── Session Navigation ────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Session Navigation ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnLoad = uibutton(cGL,'Text','Browse Master Folder', ...
           'ButtonPushedFcn',@cb_loadMaster);
       btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];

r=r+1; lblMasterPath = uilabel(cGL,'Text','No folder selected', ...
           'WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]);
       lblMasterPath.Layout.Row=r; lblMasterPath.Layout.Column=[1 2];

r=r+1; lblCurrentFolder = uilabel(cGL,'Text','—', ...
           'FontWeight','bold','WordWrap','on','FontSize',10);
       lblCurrentFolder.Layout.Row=r; lblCurrentFolder.Layout.Column=[1 2];

r=r+1; lblNavStatus = uilabel(cGL,'Text','0 / 0  folders  |  0 complete', ...
           'FontSize',10,'FontColor',[0.25 0.65 0.25]);
       lblNavStatus.Layout.Row=r; lblNavStatus.Layout.Column=[1 2];

r=r+1; btnPrev = uibutton(cGL,'Text','← Prev','ButtonPushedFcn',@cb_prevSession);
       btnPrev.Layout.Row=r; btnPrev.Layout.Column=1;
       btnNext = uibutton(cGL,'Text','Next →','ButtonPushedFcn',@cb_nextSession);
       btnNext.Layout.Row=r; btnNext.Layout.Column=2;

% ── Camera ────────────────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Camera ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Camera FPS:');
       h.Layout.Row=r; h.Layout.Column=1;
       efFps = uieditfield(cGL,'numeric','Value',30, ...
           'Limits',[1 2000],'ValueChangedFcn',@cb_updateTrace);
       efFps.Layout.Row=r; efFps.Layout.Column=2;

% ── Dot Selection ─────────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Dot Selection ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

dotCheck = gobjects(nDots,1);
dotCoord = gobjects(nDots,1);
for d = 1:nDots
    r=r+1;
    dotCheck(d) = uicheckbox(cGL,'Text',dotNames{d},'Value',d<=3, ...
        'ValueChangedFcn',@cb_updateTrace);
    dotCheck(d).Layout.Row=r; dotCheck(d).Layout.Column=1;

    dotCoord(d) = uidropdown(cGL,'Items',{'y','x','magnitude','-y','-x'},'Value','y', ...
        'ValueChangedFcn',@cb_updateTrace);
    dotCoord(d).Layout.Row=r; dotCoord(d).Layout.Column=2;
end

r=r+1; lblDotsNote = uilabel(cGL,'Text','', ...
           'FontSize',9,'FontColor',[0.7 0.5 0.1],'WordWrap','on');
       lblDotsNote.Layout.Row=r; lblDotsNote.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Likelihood thresh:');
       h.Layout.Row=r; h.Layout.Column=1;
       efLk = uieditfield(cGL,'numeric','Value',0.5,'Limits',[0 1], ...
           'ValueChangedFcn',@cb_updateTrace);
       efLk.Layout.Row=r; efLk.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','Combine dots:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddComb = uidropdown(cGL,'Items',{'average','sum'},'Value','average', ...
           'ValueChangedFcn',@cb_updateTrace);
       ddComb.Layout.Row=r; ddComb.Layout.Column=2;

r=r+1; cbInvert = uicheckbox(cGL,'Text','Invert signal','Value',false, ...
           'ValueChangedFcn',@cb_updateTrace);
       cbInvert.Layout.Row=r; cbInvert.Layout.Column=[1 2];

r=r+1; btnTrace = uibutton(cGL,'Text','Update Trace', ...
           'ButtonPushedFcn',@cb_updateTrace);
       btnTrace.Layout.Row=r; btnTrace.Layout.Column=[1 2];

% ── Peak Detection ────────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Peak Detection ──','FontWeight','bold');
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

% ── Manual Editing ────────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Manual Editing ──','FontWeight','bold');
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

% ── Save ──────────────────────────────────────────────────────────────────
r=r+1; h=uilabel(cGL,'Text','── Save ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnSave = uibutton(cGL,'Text','Save breath_peak_data.mat', ...
           'ButtonPushedFcn',@cb_save, ...
           'BackgroundColor',[0.65 0.45 0.05],'FontColor','white');
       btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];

r=r+1; btnReset = uibutton(cGL,'Text','Reset to Defaults', ...
           'ButtonPushedFcn',@cb_reset, ...
           'BackgroundColor',[0.30 0.30 0.30],'FontColor','white');
       btnReset.Layout.Row=r; btnReset.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','── Clear Data ──','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnClearAll = uibutton(cGL,'Text','Clear All Sorted Data (All Folders)', ...
           'ButtonPushedFcn',@cb_clearAllData, ...
           'BackgroundColor',[0.55 0.08 0.08],'FontColor','white');
       btnClearAll.Layout.Row=r; btnClearAll.Layout.Column=[1 2];

% ── Load persisted working params (overrides hardcoded defaults) ──────────
if ispref('BreathingPeakGUI','workingParams')
    try, applyWorkingParams(getpref('BreathingPeakGUI','workingParams')); catch, end
end

%% ── Right: plot panel ────────────────────────────────────────────────────
rightGL = uigridlayout(rootGL, [2 1]);
rightGL.Layout.Column = 2;
rightGL.RowHeight     = {'1x','1x'};
rightGL.Padding       = [0 0 0 0];
rightGL.RowSpacing    = 6;

% Main axes: processed trace + peaks
axMain = uiaxes(rightGL);
axMain.Layout.Row           = 1;
axMain.Title.String         = 'Breathing Trace  |  Left-click: add/delete peaks  |  N: navigate mode';
axMain.Title.Interpreter    = 'none';
axMain.XLabel.String        = 'Time (s)';
axMain.XLabel.Interpreter   = 'none';
axMain.YLabel.String        = 'Signal (z-score)';
axMain.YLabel.Interpreter   = 'none';
hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');

% Dot traces axes
axDots = uiaxes(rightGL);
axDots.Layout.Row           = 2;
axDots.Title.String         = 'Individual Dot Traces (normalised, offset)';
axDots.Title.Interpreter    = 'none';
axDots.XLabel.String        = 'Time (s)';
axDots.XLabel.Interpreter   = 'none';
axDots.YLabel.String        = '';
hold(axDots,'on'); grid(axDots,'on'); box(axDots,'on');

% Link x-axes
linkaxes([axMain axDots],'x');

% Persistent peak overlay handles (updated in-place — no object churn)
hPeakSc = gobjects(0);   % scatter on axMain
hPeakVl = gobjects(0);   % NaN-separated line on axDots

% Keyboard shortcut on figure
fig.KeyPressFcn = @cb_key;

%% ═══════════════════════════════════════════════════════════════════════
%% CALLBACKS
%% ═══════════════════════════════════════════════════════════════════════

    % ── Load file ────────────────────────────────────────────────────────
    % ── Browse master folder and build session list ───────────────────────
    function cb_loadMaster(~,~)
        master = uigetdir(S.masterFolder,'Select Master Folder');
        if isequal(master,0), return; end
        S.masterFolder = master;

        sessions = scanForSessions(master);
        if isempty(sessions)
            uialert(fig,'No DLC CSV files found in any subfolder.','No Sessions');
            return;
        end
        S.sessionList = sessions;
        S.sessionIdx  = 1;

        [~, mName] = fileparts(master);
        lblMasterPath.Text = mName;
        loadCurrentSession();
    end

    % ── Navigate to previous session ─────────────────────────────────────
    function cb_prevSession(~,~)
        if S.sessionIdx <= 1, return; end
        S.sessionIdx = S.sessionIdx - 1;
        loadCurrentSession();
    end

    % ── Navigate to next session ──────────────────────────────────────────
    function cb_nextSession(~,~)
        if S.sessionIdx >= numel(S.sessionList), return; end
        S.sessionIdx = S.sessionIdx + 1;
        loadCurrentSession();
    end

    % ── Load the session at S.sessionIdx ─────────────────────────────────
    function loadCurrentSession()
        sess         = S.sessionList{S.sessionIdx};
        S.folderPath = sess.folder;
        S.csvFile    = sess.csvFile;
        S.peaks      = [];
        S.rawData    = [];

        % Show current folder name
        [~, folderName] = fileparts(S.folderPath);
        lblCurrentFolder.Text = folderName;
        updateNavDisplay();

        % Read DLC CSV
        try
            raw = readmatrix(S.csvFile,'NumHeaderLines',3);
        catch ME
            uialert(fig, ME.message,'Read Error'); return;
        end
        S.rawData = raw;
        detectAndApplyDots(S.csvFile);

        % Check for previously saved results
        matFile = fullfile(S.folderPath,'breath_peak_data.mat');
        if isfile(matFile)
            try
                saved = load(matFile,'insp_onset_idx','findpeak_params');
                if isfield(saved,'findpeak_params')
                    fp = saved.findpeak_params;
                    if isfield(fp,'camera_fps'),     efFps.Value       = fp.camera_fps;     end
                    if isfield(fp,'minProminence'),  efParam(1).Value  = fp.minProminence;  end
                    if isfield(fp,'minDist_s'),      efParam(2).Value  = fp.minDist_s;      end
                    if isfield(fp,'minHeight'),      efParam(3).Value  = fp.minHeight;      end
                    if isfield(fp,'minWidth_s'),     efParam(4).Value  = fp.minWidth_s;     end
                    if isfield(fp,'snap_win_s'),     efSnap.Value      = fp.snap_win_s;     end
                    if isfield(fp,'likelihood_thr'), efLk.Value        = fp.likelihood_thr; end
                    if isfield(fp,'combine_method'), ddComb.Value      = fp.combine_method; end
                    if isfield(fp,'inverted'),       cbInvert.Value    = fp.inverted;       end
                    if isfield(fp,'dot_selection') && isfield(fp,'coord_types')
                        for d = 1:nDots
                            dotCheck(d).Value = ismember(dotNames{d}, fp.dot_selection);
                            if d <= numel(fp.coord_types)
                                dotCoord(d).Value = fp.coord_types{d};
                            end
                        end
                    end
                end
                cb_updateTrace();
                if isfield(saved,'insp_onset_idx') && ~isempty(saved.insp_onset_idx)
                    S.peaks = saved.insp_onset_idx(:);
                    redraw_peaks();
                end
                lblCurrentFolder.Text = [folderName '  ✓'];
                updateNavDisplay();
                return;
            catch
                % Corrupted mat — fall through to fresh trace
            end
        end

        cb_updateTrace();
    end

    % ── Update folder counter display ─────────────────────────────────────
    % ── Detect dots from CSV header and configure UI ──────────────────────
    function detectAndApplyDots(csvFile)
        try
            hdr   = readlines(csvFile);
            parts = strsplit(strtrim(hdr(2)), ',');
            % parts: 'bodyparts', dot1,dot1,dot1, dot2,dot2,dot2, ...
            bpNames      = strtrim(parts(2:end));
            detectedNames = bpNames(1:3:end);   % one name per dot
            nDet          = numel(detectedNames);
        catch
            nDet          = 3;                  % fallback if header unreadable
            detectedNames = {'dot1','dot2','dot3'};
        end

        nShow = min(nDet, 4);   % cap at 4

        for d = 1:4
            if d <= nShow
                dotCheck(d).Text   = detectedNames{d};
                dotCheck(d).Enable = 'on';
                dotCheck(d).Value  = true;
                dotNames{d}        = detectedNames{d};
            else
                dotCheck(d).Text   = sprintf('dot%d', d);
                dotCheck(d).Enable = 'off';
                dotCheck(d).Value  = false;
            end
        end

        if nDet > 4
            lblDotsNote.Text = sprintf('⚠ %d dots detected — showing first 4', nDet);
        else
            lblDotsNote.Text = sprintf('%d dot(s) detected', nDet);
        end
    end

    function updateNavDisplay()
        n = numel(S.sessionList);
        nDone = sum(cellfun(@(s) isfile(fullfile(s.folder,'breath_peak_data.mat')), ...
                            S.sessionList));
        lblNavStatus.Text = sprintf('%d / %d  folders  |  %d complete', ...
                                    S.sessionIdx, n, nDone);
    end

    % ── Recompute breathing trace ─────────────────────────────────────────
    function cb_updateTrace(~,~)
        if isempty(S.rawData), return; end

        fps    = efFps.Value;
        nF     = size(S.rawData,1);
        S.t    = (0:nF-1)' / fps;
        lk_thr = efLk.Value;
        comb   = ddComb.Value;

        allTraces = NaN(nF, nDots);   % all available dots  → individual plot
        traces    = NaN(nF, nDots);   % checked dots only   → breath signal

        for d = 1:nDots
            if dotCols(d,1) > size(S.rawData,2), continue; end  % dot not in CSV

            xc = S.rawData(:, dotCols(d,1));
            yc = S.rawData(:, dotCols(d,2));
            pc = S.rawData(:, dotCols(d,3));
            bad = pc < lk_thr;

            switch dotCoord(d).Value
                case 'x',         sig =  xc;
                case 'y',         sig =  yc;
                case '-x',        sig = -xc;
                case '-y',        sig = -yc;
                case 'magnitude', sig = sqrt(xc.^2 + yc.^2);
            end
            sig(bad) = NaN;

            % Linear interpolation over low-likelihood gaps
            if sum(~isnan(sig)) >= 2
                sig = fillmissing(sig,'linear','EndValues','nearest');
            end
            allTraces(:,d) = sig;          % always store for display
            if dotCheck(d).Value
                traces(:,d) = sig;         % only checked → breath signal
            end
        end

        S.dotTraces = allTraces;           % individual plot shows all dots

        % Combine checked dots
        active = arrayfun(@(d) dotCheck(d).Value, 1:nDots);
        if ~any(active)
            S.breath = zeros(nF,1);
        else
            sub = traces(:,active);
            switch comb
                case 'average', S.breath = mean(sub,2,'omitnan');
                case 'sum',     S.breath = sum(sub,2,'omitnan');
            end
        end

        % Detrend + z-score
        S.breath = detrend(S.breath);
        mu = mean(S.breath,'omitnan');
        sd = std(S.breath,0,'omitnan');
        if sd > 0, S.breath = (S.breath - mu) / sd; end

        if cbInvert.Value, S.breath = -S.breath; end

        saveWorkingParams();
        refresh_dotPlot();
        refresh_mainPlot();
    end

    % ── Detect peaks ─────────────────────────────────────────────────────
    function cb_detect(~,~)
        if isempty(S.breath)
            uialert(fig,'Load and process data first.','No Data'); return;
        end
        % Use custom finder to avoid Chronux findpeaks shadowing MATLAB's version
        S.peaks = find_peaks_local(S.breath, efFps.Value, ...
            efParam(1).Value, efParam(2).Value, ...
            efParam(3).Value, efParam(4).Value);
        saveWorkingParams();
        redraw_peaks();
        updatePeakCount();
    end

    % ── Axes click: add or delete peak ───────────────────────────────────
    function cb_axClick(~,~)
        if isempty(S.breath), return; end

        % CurrentPoint in data coordinates — only gate on the time range
        cp     = axMain.CurrentPoint;
        clickT = cp(1,1);
        if clickT < S.t(1) || clickT > S.t(end), return; end

        % Map time → nearest frame index
        [~, nearIdx] = min(abs(S.t - clickT));

        switch S.editMode
            case 'add'
                % Snap to local maximum within snap window
                snapFr = max(1, round(efSnap.Value * efFps.Value));
                lo = max(1, nearIdx - snapFr);
                hi = min(length(S.breath), nearIdx + snapFr);
                [~, imax] = max(S.breath(lo:hi));
                peakIdx = lo + imax - 1;

                % Add only if not a duplicate (within 2 frames)
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

    % ── Keyboard shortcuts ────────────────────────────────────────────────
    function cb_key(~, event)
        switch lower(event.Key)
            case 'n', cb_setNav();
            case 'a', cb_setAdd();
            case 'd', cb_setDelete();
        end
    end

    % ── Edit mode setters ─────────────────────────────────────────────────
    function cb_setAdd(~,~)
        S.editMode = 'add';
        lblMode.Text      = 'Mode: + ADD  (click → snap to local max)';
        lblMode.FontColor = [0.30 0.62 1.00];
        disableDefaultInteractivity(axMain);
        axMain.ButtonDownFcn = @cb_axClick;
    end

    function cb_setDelete(~,~)
        S.editMode = 'delete';
        lblMode.Text      = 'Mode: - DELETE  (click → remove nearest)';
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

    function cb_clearAllData(~,~)
        if isempty(S.sessionList)
            uialert(fig,'No master folder loaded.','Error'); return;
        end
        n = numel(S.sessionList);
        answer = uiconfirm(fig, ...
            sprintf(['Delete breath_peak_data.mat, breath_peak_figure.fig, ' ...
                     'and breath_peak_figure.png\nfrom all %d session folders?\n\n' ...
                     'This cannot be undone.'], n), ...
            'Confirm Clear All Data', ...
            'Options',    {'Delete All','Cancel'}, ...
            'DefaultOption','Cancel', ...
            'CancelOption','Cancel', ...
            'Icon','warning');
        if ~strcmp(answer,'Delete All'), return; end

        targets   = {'breath_peak_data.mat','breath_peak_figure.fig','breath_peak_figure.png'};
        nDeleted  = 0;
        for i = 1:n
            folder = S.sessionList{i}.folder;
            for f = 1:numel(targets)
                fp = fullfile(folder, targets{f});
                if isfile(fp)
                    delete(fp);
                    nDeleted = nDeleted + 1;
                end
            end
        end

        % Clear current session peaks from view
        S.peaks = [];
        redraw_peaks();
        updatePeakCount();
        updateNavDisplay();

        uialert(fig, sprintf('Deleted %d file(s) across %d folders.', nDeleted, n), ...
            'Done','Icon','success');
    end

    function applyWorkingParams(p)
        efFps.Value    = p.fps;
        for d = 1:nDots
            if d <= numel(p.dotsOn),    dotCheck(d).Value = p.dotsOn(d);    end
            if d <= numel(p.coordTypes), dotCoord(d).Value = p.coordTypes{d}; end
        end
        efLk.Value       = p.likelihood;
        ddComb.Value     = p.combine;
        cbInvert.Value   = p.invert;
        efParam(1).Value = p.minProm;
        efParam(2).Value = p.minDist;
        efParam(3).Value = p.minHeight;
        efParam(4).Value = p.minWidth;
        efSnap.Value     = p.snapWin;
    end

    function saveWorkingParams()
        wp.fps        = efFps.Value;
        wp.dotsOn     = arrayfun(@(d) dotCheck(d).Value, 1:nDots);
        wp.coordTypes = arrayfun(@(d) dotCoord(d).Value, 1:nDots, 'UniformOutput',false);
        wp.likelihood = efLk.Value;
        wp.combine    = ddComb.Value;
        wp.invert     = cbInvert.Value;
        wp.minProm    = efParam(1).Value;
        wp.minDist    = efParam(2).Value;
        wp.minHeight  = efParam(3).Value;
        wp.minWidth   = efParam(4).Value;
        wp.snapWin    = efSnap.Value;
        setpref('BreathingPeakGUI','workingParams',wp);
    end

    function cb_clearPeaks(~,~)
        S.peaks = [];
        redraw_peaks();
        updatePeakCount();
    end

    % ── Save ──────────────────────────────────────────────────────────────
    function cb_save(~,~)
        if isempty(S.folderPath)
            uialert(fig,'No folder loaded.','Error'); return;
        end

        insp_onset_idx    = S.peaks;            % [K×1] frame indices
        insp_onsets_t     = [];                 % [K×1] seconds
        insp_onsets_amp   = [];                 % [K×1] z-score amplitude at each peak
        insp_onsets_train = zeros(size(S.breath)); % [F×1] binary train
        if ~isempty(S.peaks) && ~isempty(S.t)
            insp_onsets_t   = S.t(S.peaks);
            insp_onsets_amp = S.breath(S.peaks);
            insp_onsets_train(S.peaks) = 1;
        end
        breath   = S.breath;
        t_breath = S.t;

        findpeak_params.minProminence  = efParam(1).Value;
        findpeak_params.minDist_s      = efParam(2).Value;
        findpeak_params.minHeight      = efParam(3).Value;
        findpeak_params.minWidth_s     = efParam(4).Value;
        findpeak_params.snap_win_s     = efSnap.Value;
        findpeak_params.camera_fps     = efFps.Value;
        findpeak_params.dot_selection  = dotNames(arrayfun(@(d) dotCheck(d).Value, 1:nDots));
        findpeak_params.coord_types    = arrayfun(@(d) dotCoord(d).Value, 1:nDots, ...
                                             'UniformOutput',false);
        findpeak_params.combine_method = ddComb.Value;
        findpeak_params.inverted       = cbInvert.Value;
        findpeak_params.likelihood_thr = efLk.Value;
        findpeak_params.csv_file       = S.csvFile;
        findpeak_params.saved_at       = datestr(now);        %#ok<TNOW1,DATST>

        outFile = fullfile(S.folderPath,'breath_peak_data.mat');
        save(outFile,'insp_onset_idx','insp_onsets_t','insp_onsets_amp','insp_onsets_train', ...
                     'breath','t_breath','findpeak_params');

        % Save figure screenshots
        figBase = fullfile(S.folderPath,'breath_peak_figure');
        saveScreenshots(figBase);

        uialert(fig, sprintf('Saved %d peaks →\n%s\n+ .fig and .png', numel(S.peaks), outFile), ...
            'Saved','Icon','success');
        if ~isempty(S.sessionList), updateNavDisplay(); end
    end

    function saveScreenshots(outBase)
        [~, sessName] = fileparts(S.folderPath);

        hFig = figure('Visible','off','Color','white','Position',[100 100 1400 820]);

        % ── Top: main breathing trace ─────────────────────────────────────
        ax1 = subplot(2,1,1,'Parent',hFig);
        hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
        copyobj(axMain.Children, ax1);
        ax1.XLim          = axMain.XLim;
        ax1.YLim          = axMain.YLim;
        ax1.XLabel.String       = axMain.XLabel.String;
        ax1.YLabel.String       = axMain.YLabel.String;
        ax1.Title.String        = [strrep(sessName,'_',' ') '  —  Breathing Trace'];
        ax1.Title.Interpreter   = 'none';
        ax1.XLabel.Interpreter  = 'none';
        ax1.YLabel.Interpreter  = 'none';

        % ── Bottom: individual dot traces ─────────────────────────────────
        ax2 = subplot(2,1,2,'Parent',hFig);
        hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
        copyobj(axDots.Children, ax2);
        ax2.XLim                = axDots.XLim;
        ax2.YLim                = axDots.YLim;
        ax2.XLabel.String       = axDots.XLabel.String;
        ax2.YLabel.String       = '';
        ax2.YTick               = [];
        ax2.Title.String        = [strrep(sessName,'_',' ') '  —  Individual Dot Traces'];
        ax2.Title.Interpreter   = 'none';
        ax2.XLabel.Interpreter  = 'none';

        savefig(hFig,  [outBase '.fig']);
        exportgraphics(hFig, [outBase '.png'], 'Resolution',150);
        close(hFig);
    end

%% ═══════════════════════════════════════════════════════════════════════
%% PLOT HELPERS
%% ═══════════════════════════════════════════════════════════════════════

    function refresh_dotPlot()
        cla(axDots); hold(axDots,'on');
        axDots.YLimMode = 'auto';
        if isempty(S.t) || isempty(S.dotTraces), return; end

        offset = 0;
        for d = nDots:-1:1
            sig = S.dotTraces(:,d);
            if all(isnan(sig)), continue; end
            rng = max(sig) - min(sig);
            if rng < eps, rng = 1; end
            sig_n = (sig - min(sig)) / rng + offset;
            h = plot(axDots, S.t, sig_n, 'Color', dotColors{d}, 'LineWidth',0.85);
            h.HitTest = 'off';
            offset = offset + 1.3;
        end
        legend(axDots,'off');
        axDots.YTick = [];
        if ~isempty(S.t), xlim(axDots,[S.t(1) S.t(end)]); end

        % Set y-range explicitly from tracked offset — no auto-scale guessing
        if offset > 0
            axDots.YLim = [-0.2, offset];
        end
        axDots.YLimMode = 'manual';

        % Create single persistent vline object
        hPeakVl = plot(axDots, NaN, NaN, '-r', 'LineWidth',0.6, 'Tag','peak_vline');
        hPeakVl.HitTest = 'off';
    end

    function refresh_mainPlot()
        cla(axMain); hold(axMain,'on');
        if isempty(S.t) || isempty(S.breath), return; end

        h = plot(axMain, S.t, S.breath, 'Color',[0 0 0], 'LineWidth',0.9);
        h.HitTest = 'off';

        h2 = yline(axMain, 0, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.6);
        h2.HitTest = 'off';

        % Re-create scatter handle after cla (pre-allocated, empty data)
        hPeakSc = scatter(axMain, [], [], 75, 'v', 'filled', ...
            'MarkerFaceColor',[1.00 0.28 0.28], 'MarkerEdgeColor','none', ...
            'Tag','peak_marker');
        hPeakSc.HitTest = 'off';

        xlim(axMain,[S.t(1) S.t(end)]);
        redraw_peaks();
        updatePeakCount();
    end

    function redraw_peaks()
        % Fast in-place update: just set XData/YData — no object churn
        valid = [];
        if ~isempty(S.peaks) && ~isempty(S.breath)
            valid = S.peaks(S.peaks >= 1 & S.peaks <= length(S.breath));
        end

        % Peak markers on main plot
        if isgraphics(hPeakSc)
            if isempty(valid)
                set(hPeakSc, 'XData',[], 'YData',[]);
            else
                set(hPeakSc, 'XData',S.t(valid), 'YData',S.breath(valid));
            end
        end

        % Vertical lines on dot plot (single NaN-separated line)
        if isgraphics(hPeakVl)
            if isempty(valid) || isempty(S.t)
                set(hPeakVl, 'XData',NaN, 'YData',NaN);
            else
                tpk = S.t(valid);   n = length(tpk);
                yl  = axDots.YLim;
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

end  % breathing_peak_gui


%% ═══════════════════════════════════════════════════════════════════════
%% LOCAL HELPER — custom peak finder
%%   Avoids Chronux findpeaks shadowing the Signal Processing Toolbox
%% ═══════════════════════════════════════════════════════════════════════
function locs = find_peaks_local(sig, fps, minProm, minDistS, minH, minWS)
%FIND_PEAKS_LOCAL  Local-max finder with prominence / distance / height / width filters.

sig  = sig(:);
n    = length(sig);
locs = [];

% ── 1. Local maxima (strict: larger than both neighbours) ────────────────
isMax = false(n,1);
isMax(2:end-1) = sig(2:end-1) > sig(1:end-2) & sig(2:end-1) > sig(3:end);
locs = find(isMax);
if isempty(locs), return; end

% ── 2. Height threshold ───────────────────────────────────────────────────
if minH ~= 0
    locs = locs(sig(locs) >= minH);
    if isempty(locs), return; end
end

% ── 3. Prominence filter ─────────────────────────────────────────────────
if minProm > 0
    prom = compute_prominence(sig, locs);
    locs = locs(prom >= minProm);
    if isempty(locs), return; end
end

% ── 4. Min distance: greedily keep highest peak, remove neighbours ────────
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

% ── 5. Min width at half-prominence ──────────────────────────────────────
if minWS > 0
    minWFr = max(1, round(minWS * fps));
    prom   = compute_prominence(sig, locs);
    keep   = true(length(locs),1);
    for i = 1:length(locs)
        pk       = locs(i);
        half_lvl = sig(pk) - prom(i)/2;
        % left crossing
        left_seg  = sig(1:pk);
        lc = find(left_seg < half_lvl, 1,'last');
        if isempty(lc), lc = 1; end
        % right crossing
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
    % Left boundary: nearest higher peak to the left, or signal start
    left_higher = locs(sig(locs) > sig(pk) & locs < pk);
    li = 1;
    if ~isempty(left_higher), li = left_higher(end); end
    % Right boundary: nearest higher peak to the right, or signal end
    right_higher = locs(sig(locs) > sig(pk) & locs > pk);
    ri = n;
    if ~isempty(right_higher), ri = right_higher(1); end
    % Prominence = peak height above highest of the two valley minima
    left_min  = min(sig(li:pk));
    right_min = min(sig(pk:ri));
    base      = max(left_min, right_min);
    prom(i)   = sig(pk) - base;
end
end  % compute_prominence


function sessions = scanForSessions(masterFolder)
%SCANFORSESSIONS  Recursively find all subfolders containing a DLC CSV file.
%   Returns a cell array of structs with fields .folder and .csvFile,
%   sorted alphabetically by folder path.

sessions = {};

% Recursive glob for DLC snapshot CSV files
hits = dir(fullfile(masterFolder,'**','*snapshot_best-20.csv'));
if isempty(hits)
    hits = dir(fullfile(masterFolder,'**','*DLC*.csv'));
end
if isempty(hits), return; end

% One entry per unique folder (take first CSV if multiple exist)
folders  = {hits.folder};
[uFolders, ia] = unique(folders,'stable');
for i = 1:numel(uFolders)
    h = hits(ia(i));
    sessions{end+1} = struct('folder', h.folder, ...        %#ok<AGROW>
                             'csvFile', fullfile(h.folder, h.name));
end
end  % scanForSessions
