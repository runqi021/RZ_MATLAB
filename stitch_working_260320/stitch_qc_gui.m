function stitch_qc_gui()
%STITCH_QC_GUI  Quality-control GUI for tile stitching pipeline
%
%   Loads a stitch_coords .mat (output of NChan_vol_stitch_zUse), lets the
%   user browse stitched Z slices, toss bad tiles at specific Z planes, and
%   save the cleaned output as multi-page TIFFs.
%
%   MULTI-CHANNEL ARCHITECTURE:
%   ---------------------------
%   The GUI auto-discovers all channel folders (ch1/, ch2/, ch3/, etc.)
%   under the dataset root by scanning for sibling chN/avgz/ directories.
%   Each channel has its own set of tile TIFF files; the stitch coordinates
%   (x, y offsets) are shared because they were computed on the reference
%   channel.
%
%   Toss matrix is 3D: [nTiles x nZ x nChannels].
%     - Toss operations (toss at Z, toss range, toss ALL Z) write across
%       ALL channels simultaneously (:, in the 3rd dim), because a bad
%       tile is usually bad on every channel.
%     - "Add back (this ch)" un-tosses a tile on the CURRENT channel only,
%       allowing per-channel exceptions (e.g. a tile is noisy on ch2 but
%       usable on ch3).
%     - Stitching, projection, and overlay display all read the toss state
%       for S.curChannel, so switching channels shows channel-specific
%       toss decisions.
%
%   "Save all" iterates over EVERY channel:
%     - Switches tile handles to that channel
%     - Re-stitches all Z slices using that channel's toss state
%     - Saves stitched TIFF and projection TIFF with channel tag in filename
%     - Saves toss_info.mat once (contains the full 3D toss matrix + channelNames)
%     - Restores the original channel view when done
%
%   PROJECTION:
%   -----------
%   updateProjection() uses buildValidMask() to determine which pixels are
%   valid at each Z. A pixel is valid if at least one NON-TOSSED tile covers
%   it (using ground-truth tile geometry from S.place, not pixel values).
%   This prevents tossed slices from contributing zeros to the projection:
%     - Max projection: tossed regions don't suppress the max
%     - Mean projection: tossed regions don't dilute the average
%
%   Usage:
%       stitch_qc_gui()
%
%   Requires: MATLAB R2021a+ with Image Processing Toolbox
%
%   Output files (per channel):
%       stitched_qc_ch2.tif           -- multi-page stitched TIFF (cleaned)
%       projection_max_qc_ch2.tif     -- single-page max/mean projection
%       toss_info.mat                 -- 3D toss matrix + log (all channels)

%% -- Shared state --------------------------------------------------------
S = struct();
S.coordFile   = '';
S.coordDir    = '';
S.zUse        = [];
S.nZ          = 0;
S.curZ        = 1;
S.nTiles      = 0;
S.H0          = 0;
S.W0          = 0;
S.outH        = 0;
S.outW        = 0;
S.x           = [];
S.y           = [];
S.place       = [];
S.tileHandles = {};     % {nTiles x 1} Tiff objects for current channel
S.tileFiles   = [];     % string array [nTiles x 1] for current channel
S.tossMatrix  = [];     % logical [nTiles x nZ x nChannels]  (3D: per-channel toss)
S.tossLog     = {};     % cell array of structs: .action, .tile, .kzRange, .timestamp
S.zCache      = [];     % single [outH x outW x nZ]
S.projImg     = [];     % single [outH x outW]
S.projType    = 'Max';
S.inClass     = 'uint16';
S.nRows       = 0;
S.nCols       = 0;
S.c_of_t      = [];
S.rB_of_t     = [];
S.selectedTile = 0;     % 0 = none
S.loaded      = false;

% Channel info
S.channelNames   = {};   % e.g. {'ch2','ch1+ch3'}
S.channelFiles   = {};   % cell of string arrays, one per channel
S.curChannel     = 1;

% Coord file fields to pass through on save
S.tileFileRef   = [];
S.tileFileOther = [];
S.mode          = '';

%% -- Figure --------------------------------------------------------------
fig = uifigure('Name','Stitch QC', ...
               'Position',[40 30 1700 950], ...
               'WindowState','maximized', ...
               'CloseRequestFcn',@cb_closeFig);

rootGL = uigridlayout(fig, [1 2]);
rootGL.ColumnWidth   = {300, '1x'};
rootGL.Padding       = [4 4 4 4];
rootGL.ColumnSpacing = 6;

%% -- Left: control panel -------------------------------------------------
ctrlPan = uipanel(rootGL, 'Title','Controls', ...
    'FontSize',12,'FontWeight','bold');
ctrlPan.Layout.Column = 1;

nCtrlRows = 35;
cGL = uigridlayout(ctrlPan, [nCtrlRows 2]);
cGL.RowHeight   = repmat({'fit'}, 1, nCtrlRows);
cGL.ColumnWidth = {'fit','1x'};
cGL.Padding     = [6 6 6 6];
cGL.RowSpacing  = 4;

r = 0;

% -- Load -----------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Load --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnLoad = uibutton(cGL,'Text','Load stitch_coords .mat', ...
           'ButtonPushedFcn',@cb_load);
       btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];

r=r+1; lblPath = uilabel(cGL,'Text','No file loaded', ...
           'WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]);
       lblPath.Layout.Row=r; lblPath.Layout.Column=[1 2];

% -- Channel --------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Channel --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Channel:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddChannel = uidropdown(cGL,'Items',{'(none)'},'Enable','off', ...
           'ValueChangedFcn',@cb_channelChanged);
       ddChannel.Layout.Row=r; ddChannel.Layout.Column=2;

% -- Z Navigation ---------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Z Navigation --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; lblZ = uilabel(cGL,'Text','Z: 0/0 (slice -)','FontWeight','bold','FontSize',12);
       lblZ.Layout.Row=r; lblZ.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Z slice:');
       h.Layout.Row=r; h.Layout.Column=1;
       spnZ = uispinner(cGL,'Value',1,'Limits',[1 1],'Step',1, ...
           'ValueChangedFcn',@cb_zSpinnerChanged);
       spnZ.Layout.Row=r; spnZ.Layout.Column=2;

r=r+1; sldZ = uislider(cGL,'Limits',[1 2],'Value',1, ...
           'ValueChangedFcn',@cb_zSliderChanged, ...
           'MajorTicks',[],'MinorTicks',[]);
       sldZ.Layout.Row=r; sldZ.Layout.Column=[1 2];

% -- Tile Info ------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Tile Info --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; lblTile = uilabel(cGL,'Text','Selected: None', ...
           'FontSize',10,'FontColor',[0.1 0.4 0.7]);
       lblTile.Layout.Row=r; lblTile.Layout.Column=[1 2];

% -- Toss Controls --------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Toss Controls --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnTossCurZ = uibutton(cGL,'Text','Toss at current Z', ...
           'ButtonPushedFcn',@cb_tossCurZ, ...
           'BackgroundColor',[0.85 0.25 0.20],'FontColor','white');
       btnTossCurZ.Layout.Row=r; btnTossCurZ.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Toss Z range:');
       h.Layout.Row=r; h.Layout.Column=1;
       pnlZRange = uigridlayout(cGL,[1 3]);
       pnlZRange.Layout.Row=r; pnlZRange.Layout.Column=2;
       pnlZRange.ColumnWidth = {'1x','fit','1x'};
       pnlZRange.Padding = [0 0 0 0];
       spnZStart = uispinner(pnlZRange,'Value',1,'Limits',[1 1],'Step',1);
       spnZStart.Layout.Row=1; spnZStart.Layout.Column=1;
       hTo = uilabel(pnlZRange,'Text','to');
       hTo.Layout.Row=1; hTo.Layout.Column=2;
       spnZEnd = uispinner(pnlZRange,'Value',1,'Limits',[1 1],'Step',1);
       spnZEnd.Layout.Row=1; spnZEnd.Layout.Column=3;

r=r+1; btnTossRange = uibutton(cGL,'Text','Toss range', ...
           'ButtonPushedFcn',@cb_tossRange, ...
           'BackgroundColor',[0.75 0.30 0.25],'FontColor','white');
       btnTossRange.Layout.Row=r; btnTossRange.Layout.Column=[1 2];

r=r+1; btnTossAllZ = uibutton(cGL,'Text','Toss ALL Z', ...
           'ButtonPushedFcn',@cb_tossAllZ, ...
           'BackgroundColor',[0.65 0.15 0.10],'FontColor','white');
       btnTossAllZ.Layout.Row=r; btnTossAllZ.Layout.Column=[1 2];

r=r+1; btnAddBack = uibutton(cGL,'Text','Add back (this ch)', ...
           'ButtonPushedFcn',@cb_addBackThisCh, ...
           'BackgroundColor',[0.20 0.65 0.35],'FontColor','white');
       btnAddBack.Layout.Row=r; btnAddBack.Layout.Column=[1 2];

r=r+1; btnUndo = uibutton(cGL,'Text','Undo last', ...
           'ButtonPushedFcn',@cb_undoLast);
       btnUndo.Layout.Row=r; btnUndo.Layout.Column=1;
       btnClearToss = uibutton(cGL,'Text','Clear all tosses', ...
           'ButtonPushedFcn',@cb_clearToss);
       btnClearToss.Layout.Row=r; btnClearToss.Layout.Column=2;

r=r+1; lblToss = uilabel(cGL,'Text','Tossed: 0 tile-Z entries', ...
           'FontSize',10,'FontColor',[0.6 0.1 0.1],'WordWrap','on');
       lblToss.Layout.Row=r; lblToss.Layout.Column=[1 2];

% -- Projection -----------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Projection --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','Proj type:');
       h.Layout.Row=r; h.Layout.Column=1;
       ddProj = uidropdown(cGL,'Items',{'Max','Mean'},'Value','Max', ...
           'ValueChangedFcn',@cb_projChanged);
       ddProj.Layout.Row=r; ddProj.Layout.Column=2;

% -- Save -----------------------------------------------------------------
r=r+1; h=uilabel(cGL,'Text','-- Save --','FontWeight','bold');
       h.Layout.Row=r; h.Layout.Column=[1 2];

r=r+1; btnSaveStitched = uibutton(cGL,'Text','Save stitched TIFF', ...
           'ButtonPushedFcn',@cb_saveStitched);
       btnSaveStitched.Layout.Row=r; btnSaveStitched.Layout.Column=[1 2];

r=r+1; btnSaveProj = uibutton(cGL,'Text','Save projection TIFF', ...
           'ButtonPushedFcn',@cb_saveProj);
       btnSaveProj.Layout.Row=r; btnSaveProj.Layout.Column=[1 2];

r=r+1; btnSaveToss = uibutton(cGL,'Text','Save toss info .mat', ...
           'ButtonPushedFcn',@cb_saveToss);
       btnSaveToss.Layout.Row=r; btnSaveToss.Layout.Column=[1 2];

r=r+1; btnSaveAll = uibutton(cGL,'Text','Save all', ...
           'ButtonPushedFcn',@cb_saveAll, ...
           'BackgroundColor',[0.18 0.55 0.30],'FontColor','white');
       btnSaveAll.Layout.Row=r; btnSaveAll.Layout.Column=[1 2];

r=r+1; lblStatus = uilabel(cGL,'Text','Ready.', ...
           'FontSize',10,'FontColor',[0.3 0.3 0.3],'WordWrap','on');
       lblStatus.Layout.Row=r; lblStatus.Layout.Column=[1 2];

% Pad remaining rows
for ir = (r+1):nCtrlRows
    h2 = uilabel(cGL,'Text','');
    h2.Layout.Row = ir; h2.Layout.Column = [1 2];
end

%% -- Right: image panels -------------------------------------------------
rightGL = uigridlayout(rootGL, [3 2]);
rightGL.Layout.Column = 2;
rightGL.RowHeight   = {'1x', 'fit', 100};
rightGL.ColumnWidth = {'1x', '1x'};
rightGL.Padding     = [2 2 2 2];
rightGL.RowSpacing  = 4;
rightGL.ColumnSpacing = 6;

% --- Top row: image axes ---
axLeft = uiaxes(rightGL);
axLeft.Layout.Row = 1; axLeft.Layout.Column = 1;
title(axLeft, 'Per-Z Stitched', 'Interpreter','none');
colormap(axLeft, gray(256));
axLeft.XTick = []; axLeft.YTick = [];
axLeft.DataAspectRatio = [1 1 1];
axLeft.ButtonDownFcn = @cb_imageClick;

axRight = uiaxes(rightGL);
axRight.Layout.Row = 1; axRight.Layout.Column = 2;
title(axRight, 'Projection', 'Interpreter','none');
colormap(axRight, gray(256));
axRight.XTick = []; axRight.YTick = [];
axRight.DataAspectRatio = [1 1 1];

% --- Middle row: slider labels ---
midGL = uigridlayout(rightGL, [1 2]);
midGL.Layout.Row = 2; midGL.Layout.Column = [1 2];
midGL.ColumnWidth = {'1x','1x'};
midGL.Padding = [0 0 0 0];

% Left B/C controls
bcLeftGL = uigridlayout(midGL, [2 3]);
bcLeftGL.Layout.Column = 1;
bcLeftGL.ColumnWidth = {'fit','1x','fit'};
bcLeftGL.RowHeight = {'fit','fit'};
bcLeftGL.Padding = [4 0 4 0];

lblMinLeft = uilabel(bcLeftGL,'Text','Min:');
lblMinLeft.Layout.Row=1; lblMinLeft.Layout.Column=1;
sldMinLeft = uislider(bcLeftGL,'Limits',[0 65535],'Value',0, ...
    'ValueChangedFcn',@(~,~) cb_updateBC('left'));
sldMinLeft.Layout.Row=1; sldMinLeft.Layout.Column=2;
sldMinLeft.MajorTicks=[]; sldMinLeft.MinorTicks=[];
btnAutoLeft = uibutton(bcLeftGL,'Text','Auto', ...
    'ButtonPushedFcn',@(~,~) cb_autoBC('left'));
btnAutoLeft.Layout.Row=1; btnAutoLeft.Layout.Column=3;

lblMaxLeft = uilabel(bcLeftGL,'Text','Max:');
lblMaxLeft.Layout.Row=2; lblMaxLeft.Layout.Column=1;
sldMaxLeft = uislider(bcLeftGL,'Limits',[0 65535],'Value',65535, ...
    'ValueChangedFcn',@(~,~) cb_updateBC('left'));
sldMaxLeft.Layout.Row=2; sldMaxLeft.Layout.Column=2;
sldMaxLeft.MajorTicks=[]; sldMaxLeft.MinorTicks=[];

% Right B/C controls
bcRightGL = uigridlayout(midGL, [2 3]);
bcRightGL.Layout.Column = 2;
bcRightGL.ColumnWidth = {'fit','1x','fit'};
bcRightGL.RowHeight = {'fit','fit'};
bcRightGL.Padding = [4 0 4 0];

lblMinRight = uilabel(bcRightGL,'Text','Min:');
lblMinRight.Layout.Row=1; lblMinRight.Layout.Column=1;
sldMinRight = uislider(bcRightGL,'Limits',[0 65535],'Value',0, ...
    'ValueChangedFcn',@(~,~) cb_updateBC('right'));
sldMinRight.Layout.Row=1; sldMinRight.Layout.Column=2;
sldMinRight.MajorTicks=[]; sldMinRight.MinorTicks=[];
btnAutoRight = uibutton(bcRightGL,'Text','Auto', ...
    'ButtonPushedFcn',@(~,~) cb_autoBC('right'));
btnAutoRight.Layout.Row=1; btnAutoRight.Layout.Column=3;

lblMaxRight = uilabel(bcRightGL,'Text','Max:');
lblMaxRight.Layout.Row=2; lblMaxRight.Layout.Column=1;
sldMaxRight = uislider(bcRightGL,'Limits',[0 65535],'Value',65535, ...
    'ValueChangedFcn',@(~,~) cb_updateBC('right'));
sldMaxRight.Layout.Row=2; sldMaxRight.Layout.Column=2;
sldMaxRight.MajorTicks=[]; sldMaxRight.MinorTicks=[];

% --- Bottom row: histograms ---
axHistLeft = uiaxes(rightGL);
axHistLeft.Layout.Row = 3; axHistLeft.Layout.Column = 1;
title(axHistLeft, 'Histogram (left)', 'Interpreter','none','FontSize',9);
axHistLeft.YTick = [];

axHistRight = uiaxes(rightGL);
axHistRight.Layout.Row = 3; axHistRight.Layout.Column = 2;
title(axHistRight, 'Histogram (right)', 'Interpreter','none','FontSize',9);
axHistRight.YTick = [];

% Handles for overlay graphics (populated after load)
hImgLeft      = [];   % image handle
hImgRight     = [];   % image handle
hTileBordersL = {};   % cell of line handles (left panel)
hTileBordersR = {};   % cell of line handles (right panel)
hSelectBorder = [];   % line handle for selected tile (cyan)
hTossPatchesL = {};   % cell of patch handles for tossed tiles (left)

%% ========================================================================
%  CALLBACKS
%% ========================================================================

% -- Load -----------------------------------------------------------------
    function cb_load(~,~)
        [fname, fdir] = uigetfile('*.mat','Select stitch_coords .mat');
        if isequal(fname,0), return; end
        coordFile = fullfile(fdir, fname);
        loadCoordFile(coordFile);
    end

    function loadCoordFile(coordFile)
        setStatus('Loading coordinates...');
        drawnow;

        d = load(coordFile);

        % Validate required fields
        reqFields = {'x','y','tileFileRef','nRows','nCols','zUse','H0','W0','outH','outW'};
        for fi = 1:numel(reqFields)
            if ~isfield(d, reqFields{fi})
                uialert(fig, sprintf('Missing field: %s', reqFields{fi}), 'Load Error');
                return;
            end
        end

        S.coordFile   = coordFile;
        S.coordDir    = fileparts(coordFile);
        S.x           = d.x(:);
        S.y           = d.y(:);
        S.nTiles      = numel(d.x);
        S.H0          = d.H0;
        S.W0          = d.W0;
        S.outH        = d.outH;
        S.outW        = d.outW;
        S.zUse        = d.zUse(:).';
        S.nZ          = numel(S.zUse);
        S.nRows       = d.nRows;
        S.nCols       = d.nCols;
        S.mode        = '';
        if isfield(d,'mode'),       S.mode      = d.mode; end
        if isfield(d,'c_of_t'),     S.c_of_t    = d.c_of_t(:); end
        if isfield(d,'rB_of_t'),    S.rB_of_t   = d.rB_of_t(:); end

        % --- Dynamic path resolution -----------------------------------------
        % Derive root directory from .mat location:
        %   <root>/matlab_stitch/stitch_coords_from_ch*.mat
        rootDir = fileparts(S.coordDir);

        % Parse ref channel ID from coord filename (e.g. stitch_coords_from_ch2.mat -> 2)
        [~, coordBN] = fileparts(coordFile);
        tokCoord = regexp(coordBN, 'ch(\d+)', 'tokens', 'once');

        % Rebuild tileFileRef paths relative to current rootDir
        storedRef = d.tileFileRef(:);
        refChID = parseChIDFromPath(storedRef(1));
        if isnan(refChID) && ~isempty(tokCoord)
            refChID = str2double(tokCoord{1});
        end
        S.tileFileRef = rebuildTilePaths(storedRef, rootDir, refChID);

        % Rebuild tileFileOther if present
        S.tileFileOther = [];
        if isfield(d,'tileFileOther') && ~isempty(d.tileFileOther)
            storedOther = d.tileFileOther(:);
            otherChID = parseChIDFromPath(storedOther(1));
            S.tileFileOther = rebuildTilePaths(storedOther, rootDir, otherChID);
        end

        % Validate that rebuilt ref tile files exist
        missingIdx = [];
        for ti = 1:numel(S.tileFileRef)
            if ~isfile(S.tileFileRef(ti))
                missingIdx(end+1) = ti; %#ok<AGROW>
            end
        end
        if ~isempty(missingIdx)
            uialert(fig, sprintf('Tile file not found:\n%s\n(%d of %d tiles missing)', ...
                char(S.tileFileRef(missingIdx(1))), numel(missingIdx), numel(S.tileFileRef)), ...
                'Path Error');
            return;
        end

        % --- Discover all available channels from disk ------------------------
        S.channelNames = {};
        S.channelFiles = {};

        % Ref channel first
        refChName = sprintf('ch%d', refChID);
        S.channelNames{end+1} = refChName;
        S.channelFiles{end+1} = S.tileFileRef;

        % Scan rootDir for other ch*/avgz/ folders
        chDirs = dir(fullfile(rootDir, 'ch*'));
        chDirs = chDirs([chDirs.isdir]);
        for ci = 1:numel(chDirs)
            tok = regexp(chDirs(ci).name, '^ch(\d+)$', 'tokens', 'once');
            if isempty(tok), continue; end
            thisChID = str2double(tok{1});
            if thisChID == refChID, continue; end  % already added
            avgzDir = fullfile(rootDir, chDirs(ci).name, 'avgz');
            if ~isfolder(avgzDir), continue; end
            % Try to build tile list by matching basenames from ref
            otherFiles = rebuildTilePaths(storedRef, rootDir, thisChID);
            % Verify at least the first tile exists
            if isfile(otherFiles(1))
                S.channelNames{end+1} = sprintf('ch%d', thisChID);
                S.channelFiles{end+1} = otherFiles;
            end
        end

        % Channel dropdown
        ddChannel.Items = S.channelNames;
        ddChannel.Value = S.channelNames{1};
        ddChannel.Enable = 'on';
        if numel(S.channelNames) < 2
            ddChannel.Enable = 'off';
        end
        S.curChannel = 1;

        % Close any previously open tile handles
        closeTileHandles();

        % Open Tiff handles for current channel
        S.tileFiles = S.channelFiles{S.curChannel};
        openTileHandles();

        % Determine input class from first tile
        tTest = Tiff(S.tileFiles(1),'r');
        testPage = tTest.read();
        S.inClass = class(testPage);
        tTest.close();

        % Compute place struct
        computePlace();

        % Init toss matrix
        S.nCh        = numel(S.channelNames);
        S.tossMatrix = false(S.nTiles, S.nZ, S.nCh);
        S.tossLog    = {};

        % Update Z controls
        S.curZ = 1;
        spnZ.Limits = [1 S.nZ];
        spnZ.Value  = 1;
        if S.nZ > 1
            sldZ.Limits = [1 S.nZ];
        else
            sldZ.Limits = [1 1.001]; % uislider needs range > 0
        end
        sldZ.Value = 1;
        spnZStart.Limits = [1 S.nZ]; spnZStart.Value = 1;
        spnZEnd.Limits   = [1 S.nZ]; spnZEnd.Value   = S.nZ;

        S.selectedTile = 0;

        % Update path label
        lblPath.Text = S.coordFile;

        % Stitch all Z slices into cache
        stitchAllZWithWaitbar();

        % Compute projection
        updateProjection();

        S.loaded = true;

        % Display
        showLeftPanel();
        showRightPanel();
        updateHistograms();
        cb_autoBC('left');
        cb_autoBC('right');
        updateZLabel();
        updateTileLabel();
        updateTossLabel();

        setStatus(sprintf('Loaded: %d tiles, %d Z slices, %dx%d mosaic', ...
            S.nTiles, S.nZ, S.outH, S.outW));
    end

% -- Channel changed ------------------------------------------------------
    function cb_channelChanged(~,~)
        if ~S.loaded, return; end
        idx = find(strcmp(ddChannel.Value, S.channelNames), 1);
        if isempty(idx), return; end
        S.curChannel = idx;

        closeTileHandles();
        S.tileFiles = S.channelFiles{S.curChannel};
        openTileHandles();

        % Re-stitch
        stitchAllZWithWaitbar();
        updateProjection();
        showLeftPanel();
        showRightPanel();
        updateHistograms();
        cb_autoBC('left');
        cb_autoBC('right');
        setStatus(sprintf('Switched to channel: %s', ddChannel.Value));
    end

% -- Z navigation ---------------------------------------------------------
    function cb_zSpinnerChanged(~,~)
        if ~S.loaded, return; end
        S.curZ = round(spnZ.Value);
        sldZ.Value = S.curZ;
        updateZLabel();
        showLeftPanel();
        updateHistLeft();
    end

    function cb_zSliderChanged(~,~)
        if ~S.loaded, return; end
        S.curZ = round(sldZ.Value);
        spnZ.Value = S.curZ;
        updateZLabel();
        showLeftPanel();
        updateHistLeft();
    end

    function updateZLabel()
        if S.nZ > 0
            lblZ.Text = sprintf('Z: %d/%d (slice %d)', S.curZ, S.nZ, S.zUse(S.curZ));
        else
            lblZ.Text = 'Z: 0/0 (slice -)';
        end
    end

% -- Image click (tile selection) -----------------------------------------
    function cb_imageClick(~, evt)
        if ~S.loaded, return; end
        cp = evt.IntersectionPoint;
        cx = cp(1); cy = cp(2);

        % Find which tile contains this point; pick closest center if overlap
        bestT = 0;
        bestDist = Inf;
        for t = 1:S.nTiles
            rr = S.place(t).rr;
            cc = S.place(t).cc;
            if cy >= rr(1) && cy <= rr(end) && cx >= cc(1) && cx <= cc(end)
                % Compute distance to tile center
                ctrY = (rr(1) + rr(end)) / 2;
                ctrX = (cc(1) + cc(end)) / 2;
                dd = (cx - ctrX)^2 + (cy - ctrY)^2;
                if dd < bestDist
                    bestDist = dd;
                    bestT = t;
                end
            end
        end

        S.selectedTile = bestT;
        updateTileLabel();
        drawTileOverlays(axLeft, 'left');
    end

    function updateTileLabel()
        if S.selectedTile > 0
            t = S.selectedTile;
            gridStr = '';
            if ~isempty(S.c_of_t) && ~isempty(S.rB_of_t)
                gridStr = sprintf(', grid(%d,%d)', S.rB_of_t(t), S.c_of_t(t));
            end
            isTossed = '';
            if S.loaded && S.tossMatrix(t, S.curZ, S.curChannel)
                isTossed = ' [TOSSED]';
            end
            lblTile.Text = sprintf('Selected: Tile %d%s%s', t, gridStr, isTossed);
        else
            lblTile.Text = 'Selected: None';
        end
    end

% -- Toss callbacks -------------------------------------------------------
    % Toss operations write (:) in the channel dimension — a tossed tile is
    % tossed on ALL channels by default.  Use "Add back (this ch)" to
    % selectively un-toss on one channel.
    function cb_tossCurZ(~,~)
        if ~S.loaded || S.selectedTile == 0
            uialert(fig, 'Select a tile first by clicking on the left image.', 'No Tile');
            return;
        end
        t = S.selectedTile;
        kz = S.curZ;
        if all(S.tossMatrix(t, kz, :))
            setStatus(sprintf('Tile %d already tossed at Z=%d (all ch)', t, kz));
            return;
        end
        S.tossMatrix(t, kz, :) = true;
        logToss('toss', t, kz, kz);
        restitchZRange(kz, kz);
        updateProjection();
        refreshDisplay();
        setStatus(sprintf('Tossed tile %d at Z=%d', t, kz));
    end

    function cb_tossRange(~,~)
        if ~S.loaded || S.selectedTile == 0
            uialert(fig, 'Select a tile first.', 'No Tile');
            return;
        end
        t = S.selectedTile;
        kzA = round(spnZStart.Value);
        kzB = round(spnZEnd.Value);
        if kzA > kzB
            tmp = kzA; kzA = kzB; kzB = tmp;
        end
        anyNew = false;
        for kz = kzA:kzB
            if ~all(S.tossMatrix(t, kz, :))
                S.tossMatrix(t, kz, :) = true;
                anyNew = true;
            end
        end
        if anyNew
            logToss('toss_range', t, kzA, kzB);
            restitchZRange(kzA, kzB);
            updateProjection();
            refreshDisplay();
        end
        setStatus(sprintf('Tossed tile %d at Z=%d:%d', t, kzA, kzB));
    end

    function cb_tossAllZ(~,~)
        if ~S.loaded || S.selectedTile == 0
            uialert(fig, 'Select a tile first.', 'No Tile');
            return;
        end
        t = S.selectedTile;
        S.tossMatrix(t, :, :) = true;
        logToss('toss_all_z', t, 1, S.nZ);
        restitchZRange(1, S.nZ);
        updateProjection();
        refreshDisplay();
        setStatus(sprintf('Tossed tile %d at ALL Z', t));
    end

    function cb_addBackThisCh(~,~)
        % Un-toss the selected tile at the current Z for the CURRENT channel
        % only.  This is the per-channel exception mechanism: all toss ops
        % write to every channel (`:` in 3rd dim), but this lets the user
        % selectively recover a tile on one channel while keeping it tossed
        % on others.  Saves prevState for full undo support.
        if ~S.loaded || S.selectedTile == 0
            uialert(fig, 'Select a tile first by clicking on the left image.', 'No Tile');
            return;
        end
        t = S.selectedTile;
        kz = S.curZ;
        ch = S.curChannel;
        if ~S.tossMatrix(t, kz, ch)
            setStatus(sprintf('Tile %d not tossed at Z=%d on %s', t, kz, S.channelNames{ch}));
            return;
        end
        prevState = S.tossMatrix;
        S.tossMatrix(t, kz, ch) = false;
        entry.action    = 'add_back';
        entry.tile      = t;
        entry.kzStart   = kz;
        entry.kzEnd     = kz;
        entry.timestamp = datetime('now');
        entry.prevState = prevState;
        S.tossLog{end+1} = entry;
        restitchZRange(kz, kz);
        updateProjection();
        refreshDisplay();
        setStatus(sprintf('Added back tile %d at Z=%d on %s only', t, kz, S.channelNames{ch}));
    end

    function cb_undoLast(~,~)
        if ~S.loaded || isempty(S.tossLog)
            setStatus('Nothing to undo.');
            return;
        end
        entry = S.tossLog{end};
        S.tossLog(end) = [];

        % Reverse the action
        t = entry.tile;
        kzA = entry.kzStart;
        kzB = entry.kzEnd;
        switch entry.action
            case {'toss','toss_range','toss_all_z'}
                S.tossMatrix(t, kzA:kzB, :) = false;
            case {'add_back','clear_all'}
                % Restore: we saved the previous state
                if isfield(entry, 'prevState')
                    S.tossMatrix = entry.prevState;
                    kzA = 1; kzB = S.nZ;
                end
        end
        restitchZRange(kzA, kzB);
        updateProjection();
        refreshDisplay();
        setStatus(sprintf('Undid: %s tile %d', entry.action, t));
    end

    function cb_clearToss(~,~)
        if ~S.loaded, return; end
        if ~any(S.tossMatrix(:))
            setStatus('No tosses to clear.');
            return;
        end
        prevState = S.tossMatrix;
        S.tossMatrix(:) = false;
        entry.action = 'clear_all';
        entry.tile = 0;
        entry.kzStart = 1;
        entry.kzEnd = S.nZ;
        entry.timestamp = datetime('now');
        entry.prevState = prevState;
        S.tossLog{end+1} = entry;
        restitchZRange(1, S.nZ);
        updateProjection();
        refreshDisplay();
        setStatus('Cleared all tosses.');
    end

    function logToss(action, tile, kzStart, kzEnd)
        entry.action    = action;
        entry.tile      = tile;
        entry.kzStart   = kzStart;
        entry.kzEnd     = kzEnd;
        entry.timestamp = datetime('now');
        S.tossLog{end+1} = entry;
        updateTossLabel();
    end

    function updateTossLabel()
        n = nnz(S.tossMatrix(:,:,S.curChannel));
        lblToss.Text = sprintf('Tossed: %d tile-Z entries (this ch)', n);
    end

% -- Projection type ------------------------------------------------------
    function cb_projChanged(~,~)
        if ~S.loaded, return; end
        S.projType = ddProj.Value;
        updateProjection();
        showRightPanel();
        updateHistRight();
    end

% -- B/C controls ---------------------------------------------------------
    function cb_updateBC(panel)
        if ~S.loaded, return; end
        switch panel
            case 'left'
                cmin = sldMinLeft.Value;
                cmax = sldMaxLeft.Value;
                if cmax <= cmin, cmax = cmin + 1; end
                axLeft.CLim = [cmin cmax];
            case 'right'
                cmin = sldMinRight.Value;
                cmax = sldMaxRight.Value;
                if cmax <= cmin, cmax = cmin + 1; end
                axRight.CLim = [cmin cmax];
        end
    end

    function cb_autoBC(panel)
        if ~S.loaded, return; end
        switch panel
            case 'left'
                img = S.zCache(:,:,S.curZ);
                vals = img(img > 0);
                if isempty(vals), vals = img(:); end
                lo = double(prctile(vals, 1));
                hi = double(prctile(vals, 99.5));
                if hi <= lo, hi = lo + 1; end
                sldMinLeft.Limits = [0 max(hi*2, 1)];
                sldMaxLeft.Limits = [0 max(hi*2, 1)];
                sldMinLeft.Value = lo;
                sldMaxLeft.Value = hi;
                axLeft.CLim = [lo hi];
            case 'right'
                img = S.projImg;
                vals = img(img > 0);
                if isempty(vals), vals = img(:); end
                lo = double(prctile(vals, 1));
                hi = double(prctile(vals, 99.5));
                if hi <= lo, hi = lo + 1; end
                sldMinRight.Limits = [0 max(hi*2, 1)];
                sldMaxRight.Limits = [0 max(hi*2, 1)];
                sldMinRight.Value = lo;
                sldMaxRight.Value = hi;
                axRight.CLim = [lo hi];
        end
    end

% -- Save callbacks -------------------------------------------------------
    function cb_saveStitched(~,~)
        if ~S.loaded, return; end
        [fname, fdir] = uiputfile('*.tif','Save stitched TIFF', ...
            fullfile(S.coordDir, 'stitched_qc.tif'));
        if isequal(fname,0), return; end
        outFile = fullfile(fdir, fname);
        saveStitchedTiff(outFile);
    end

    function cb_saveProj(~,~)
        if ~S.loaded, return; end
        [fname, fdir] = uiputfile('*.tif','Save projection TIFF', ...
            fullfile(S.coordDir, sprintf('projection_%s_qc.tif', lower(S.projType))));
        if isequal(fname,0), return; end
        outFile = fullfile(fdir, fname);
        saveProjectionTiff(outFile);
    end

    function cb_saveToss(~,~)
        if ~S.loaded, return; end
        [fname, fdir] = uiputfile('*.mat','Save toss info', ...
            fullfile(S.coordDir, 'toss_info.mat'));
        if isequal(fname,0), return; end
        outFile = fullfile(fdir, fname);
        saveTossInfo(outFile);
    end

    function cb_saveAll(~,~)
        % Save stitched TIFFs and projections for EVERY channel, not just
        % the one currently displayed.  For each channel we:
        %   1. Switch tile file handles to that channel's TIFFs
        %   2. Re-stitch all Z slices (respecting that channel's toss state
        %      via S.tossMatrix(:,:,chi))
        %   3. Recompute the projection (with toss-aware valid mask)
        %   4. Write stitched_qc_<chTag>.tif and projection_<type>_qc_<chTag>.tif
        % After the loop, toss_info.mat is saved once (it contains the full
        % 3D toss matrix shared across channels).  Finally the GUI is
        % restored to whichever channel the user was viewing before saving.
        if ~S.loaded, return; end
        setStatus('Saving all...');
        drawnow;

        origChannel = S.curChannel;

        for chi = 1:numel(S.channelNames)
            chTag = S.channelNames{chi};  % e.g. 'ch2'

            % Switch to this channel
            S.curChannel = chi;
            closeTileHandles();
            S.tileFiles = S.channelFiles{chi};
            openTileHandles();
            stitchAllZWithWaitbar();
            updateProjection();

            % Save with channel tag in filename
            saveStitchedTiff(fullfile(S.coordDir, sprintf('stitched_qc_%s.tif', chTag)));
            saveProjectionTiff(fullfile(S.coordDir, sprintf('projection_%s_qc_%s.tif', lower(S.projType), chTag)));
        end

        % Save toss info once (shared across channels)
        saveTossInfo(fullfile(S.coordDir, 'toss_info.mat'));

        % Restore original channel
        S.curChannel = origChannel;
        closeTileHandles();
        S.tileFiles = S.channelFiles{S.curChannel};
        openTileHandles();
        stitchAllZWithWaitbar();
        updateProjection();
        refreshDisplay();

        setStatus('All files saved (all channels).');
    end

    function saveStitchedTiff(outFile)
        setStatus(sprintf('Writing stitched TIFF (%d pages)...', S.nZ));
        drawnow;

        [bitsPerSample, sampleFormat] = class_to_tiff_format(S.inClass);

        bytesPerPx = bitsPerSample / 8;
        estBytes = double(S.outH) * double(S.outW) * double(S.nZ) * double(bytesPerPx);
        useBigTiff = estBytes > 3.5e9;

        if isfile(outFile), delete(outFile); end
        if useBigTiff
            tout = Tiff(outFile, 'w8');
        else
            tout = Tiff(outFile, 'w');
        end
        c_cleanup = onCleanup(@() safe_close_tiff(tout));

        tag.ImageLength       = S.outH;
        tag.ImageWidth        = S.outW;
        tag.Photometric       = Tiff.Photometric.MinIsBlack;
        tag.SamplesPerPixel   = 1;
        tag.BitsPerSample     = bitsPerSample;
        tag.SampleFormat      = sampleFormat;
        tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tag.Compression       = Tiff.Compression.None;
        tag.RowsPerStrip      = 64;
        tag.Software          = 'stitch_qc_gui';

        wb = uiprogressdlg(fig,'Title','Saving stitched TIFF','Message','Writing...', ...
            'Cancelable','on');
        for kz = 1:S.nZ
            if wb.CancelRequested, break; end
            wb.Value = kz / S.nZ;
            wb.Message = sprintf('Z %d / %d', kz, S.nZ);

            page = S.zCache(:,:,kz);
            page = cast(round(double(page)), S.inClass);

            tout.setTag(tag);
            tout.write(page);
            if kz < S.nZ
                tout.writeDirectory();
            end
        end
        close(wb);
        setStatus(sprintf('Saved: %s', outFile));
    end

    function saveProjectionTiff(outFile)
        setStatus('Writing projection TIFF...');
        drawnow;

        [bitsPerSample, sampleFormat] = class_to_tiff_format(S.inClass);

        if isfile(outFile), delete(outFile); end
        tout = Tiff(outFile, 'w');
        c_cleanup = onCleanup(@() safe_close_tiff(tout));

        tag.ImageLength       = S.outH;
        tag.ImageWidth        = S.outW;
        tag.Photometric       = Tiff.Photometric.MinIsBlack;
        tag.SamplesPerPixel   = 1;
        tag.BitsPerSample     = bitsPerSample;
        tag.SampleFormat      = sampleFormat;
        tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tag.Compression       = Tiff.Compression.None;
        tag.RowsPerStrip      = 64;
        tag.Software          = sprintf('stitch_qc_gui %s projection', S.projType);

        page = cast(round(double(S.projImg)), S.inClass);
        tout.setTag(tag);
        tout.write(page);

        setStatus(sprintf('Saved: %s', outFile));
    end

    function saveTossInfo(outFile)
        setStatus('Saving toss info...');
        drawnow;

        tossMatrix    = S.tossMatrix;
        tossLog       = S.tossLog;
        zUse          = S.zUse;
        tileFileRef   = S.tileFileRef;
        tileFileOther = S.tileFileOther;
        coordFile     = S.coordFile;
        nRows         = S.nRows;
        nCols         = S.nCols;
        c_of_t        = S.c_of_t;
        rB_of_t       = S.rB_of_t;
        channelNames  = S.channelNames;
        projType      = S.projType;
        timestamp     = datetime('now');

        save(outFile, 'tossMatrix','tossLog','zUse','tileFileRef','tileFileOther', ...
            'coordFile','nRows','nCols','c_of_t','rB_of_t','channelNames','projType','timestamp');

        setStatus(sprintf('Saved: %s', outFile));
    end

% -- Close ----------------------------------------------------------------
    function cb_closeFig(~,~)
        closeTileHandles();
        delete(fig);
    end

%% ========================================================================
%  INTERNAL HELPERS
%% ========================================================================

    function openTileHandles()
        S.tileHandles = cell(S.nTiles, 1);
        for t = 1:S.nTiles
            if ~isfile(S.tileFiles(t))
                uialert(fig, sprintf('Tile file not found:\n%s', S.tileFiles(t)), 'Missing File');
                error('stitch_qc_gui:missingTile', 'Tile file not found: %s', S.tileFiles(t));
            end
            S.tileHandles{t} = Tiff(S.tileFiles(t), 'r');
        end
    end

    function closeTileHandles()
        if ~isempty(S.tileHandles)
            for t = 1:numel(S.tileHandles)
                if ~isempty(S.tileHandles{t})
                    try S.tileHandles{t}.close(); catch, end
                end
            end
            S.tileHandles = {};
        end
    end

    function computePlace()
        wTile = linear_blend_weights(S.H0, S.W0);
        S.place = struct([]);
        for t = 1:S.nTiles
            xi = S.x(t); yi = S.y(t);
            xInt = floor(xi); yInt = floor(yi);
            fx = xi - xInt; fy = yi - yInt;
            r0 = yInt + 1; c0 = xInt + 1;
            rr = r0:(r0 + S.H0 - 1);
            cc = c0:(c0 + S.W0 - 1);
            if (abs(fx) > 1e-6 || abs(fy) > 1e-6)
                w2 = shift_bilinear(wTile, fx, fy);
            else
                w2 = wTile;
            end
            S.place(t).fx = fx;
            S.place(t).fy = fy;
            S.place(t).rr = rr;
            S.place(t).cc = cc;
            S.place(t).w2 = w2;
        end
    end

    function img = stitchOneZ(kz)
        % Stitch one Z slice using linear-blend weights.  Skips tiles that
        % are tossed for the CURRENT channel (3rd dim of tossMatrix), so
        % different channels can have different toss states.
        iz = S.zUse(kz);
        acc  = zeros(S.outH, S.outW, 'double');
        wacc = zeros(S.outH, S.outW, 'double');
        for t = 1:S.nTiles
            if S.tossMatrix(t, kz, S.curChannel), continue; end
            S.tileHandles{t}.setDirectory(iz);
            I = double(S.tileHandles{t}.read());
            fx = S.place(t).fx; fy = S.place(t).fy;
            if abs(fx) > 1e-6 || abs(fy) > 1e-6
                I2 = shift_bilinear(I, fx, fy);
            else
                I2 = I;
            end
            rr = S.place(t).rr;
            cc = S.place(t).cc;
            w2 = S.place(t).w2;
            acc(rr, cc)  = acc(rr, cc)  + I2 .* w2;
            wacc(rr, cc) = wacc(rr, cc) + w2;
        end
        img = single(acc ./ max(wacc, eps));
    end

    function stitchAllZWithWaitbar()
        S.zCache = zeros(S.outH, S.outW, S.nZ, 'single');
        wb = uiprogressdlg(fig,'Title','Stitching Z slices','Message','Please wait...', ...
            'Cancelable','on');
        for kz = 1:S.nZ
            if wb.CancelRequested, break; end
            wb.Value = kz / S.nZ;
            wb.Message = sprintf('Z %d / %d', kz, S.nZ);
            S.zCache(:,:,kz) = stitchOneZ(kz);
        end
        close(wb);
    end

    function restitchZRange(kzA, kzB)
        for kz = kzA:kzB
            S.zCache(:,:,kz) = stitchOneZ(kz);
        end
    end

    function updateProjection()
        % Build max or mean projection, excluding tossed tile regions.
        %
        % Why not just max/mean over zCache directly?
        %   Tossed tiles leave 0 in zCache (no tiles contributed).  For
        %   mean projection these zeros dilute the average.  For max it's
        %   mostly harmless but conceptually wrong -- a pixel that was
        %   "not measured" at a given Z should not participate at all.
        %
        % We use buildValidMask(kz) to get a per-pixel logical mask based
        % on ground-truth tile geometry (S.place footprints) rather than
        % pixel values, because real data CAN contain legitimate zeros.
        % The mask is channel-aware via S.tossMatrix(:,:,S.curChannel).
        switch S.projType
            case 'Max'
                projMax = -inf(S.outH, S.outW, 'single');
                for kz = 1:S.nZ
                    slice = S.zCache(:,:,kz);
                    valid = buildValidMask(kz);
                    update = valid & (slice > projMax);
                    projMax(update) = slice(update);
                end
                projMax(isinf(projMax)) = 0;
                S.projImg = projMax;
            case 'Mean'
                projSum = zeros(S.outH, S.outW, 'single');
                projCnt = zeros(S.outH, S.outW, 'single');
                for kz = 1:S.nZ
                    slice = S.zCache(:,:,kz);
                    valid = buildValidMask(kz);
                    projSum = projSum + slice .* single(valid);
                    projCnt = projCnt + single(valid);
                end
                S.projImg = projSum ./ max(projCnt, 1);
        end
    end

    function valid = buildValidMask(kz)
        % Per-pixel validity mask for Z slice kz on the current channel.
        % A pixel is valid if at least one NON-TOSSED tile covers it,
        % determined from tile footprints (S.place(t).rr/cc), NOT from
        % pixel values.  This is channel-aware: reads S.curChannel.
        valid = false(S.outH, S.outW);
        for t = 1:S.nTiles
            if ~S.tossMatrix(t, kz, S.curChannel)
                valid(S.place(t).rr, S.place(t).cc) = true;
            end
        end
    end

    function showLeftPanel()
        if ~S.loaded, return; end
        img = S.zCache(:,:,S.curZ);

        if isempty(hImgLeft) || ~isvalid(hImgLeft)
            hImgLeft = imagesc(axLeft, img);
            colormap(axLeft, gray(256));
            axLeft.DataAspectRatio = [1 1 1];
            axLeft.XTick = []; axLeft.YTick = [];
            hImgLeft.ButtonDownFcn = @cb_imageClick;
            axLeft.ButtonDownFcn = @cb_imageClick;
        else
            hImgLeft.CData = img;
            axLeft.XLim = [0.5, S.outW + 0.5];
            axLeft.YLim = [0.5, S.outH + 0.5];
        end

        nTossed = nnz(S.tossMatrix(:, S.curZ, S.curChannel));
        title(axLeft, sprintf('Z=%d/%d (slice %d) | %d tossed', ...
            S.curZ, S.nZ, S.zUse(S.curZ), nTossed), 'Interpreter','none');

        drawTileOverlays(axLeft, 'left');
    end

    function showRightPanel()
        if ~S.loaded, return; end
        img = S.projImg;

        if isempty(hImgRight) || ~isvalid(hImgRight)
            hImgRight = imagesc(axRight, img);
            colormap(axRight, gray(256));
            axRight.DataAspectRatio = [1 1 1];
            axRight.XTick = []; axRight.YTick = [];
        else
            hImgRight.CData = img;
            axRight.XLim = [0.5, S.outW + 0.5];
            axRight.YLim = [0.5, S.outH + 0.5];
        end

        title(axRight, sprintf('%s Projection', S.projType), 'Interpreter','none');

        drawTileOverlays(axRight, 'right');
    end

    function drawTileOverlays(ax, panel)
        % Clear existing overlays for this panel
        switch panel
            case 'left'
                cellfun(@(h) delete(h(isvalid(h))), hTileBordersL, 'UniformOutput', false);
                hTileBordersL = {};
                cellfun(@(h) delete(h(isvalid(h))), hTossPatchesL, 'UniformOutput', false);
                hTossPatchesL = {};
                if ~isempty(hSelectBorder) && isvalid(hSelectBorder)
                    delete(hSelectBorder);
                end
                hSelectBorder = [];
            case 'right'
                cellfun(@(h) delete(h(isvalid(h))), hTileBordersR, 'UniformOutput', false);
                hTileBordersR = {};
        end

        hold(ax, 'on');

        for t = 1:S.nTiles
            rr = S.place(t).rr;
            cc = S.place(t).cc;
            r1 = rr(1) - 0.5; r2 = rr(end) + 0.5;
            c1 = cc(1) - 0.5; c2 = cc(end) + 0.5;

            % Tile boundary (thin gray)
            hLine = line(ax, [c1 c2 c2 c1 c1], [r1 r1 r2 r2 r1], ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HitTest','off');

            switch panel
                case 'left'
                    hTileBordersL{end+1} = hLine; %#ok<AGROW>

                    % Tossed overlay (red semi-transparent)
                    if S.tossMatrix(t, S.curZ, S.curChannel)
                        hp = patch(ax, [c1 c2 c2 c1], [r1 r1 r2 r2], ...
                            'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', ...
                            'LineWidth', 1.5, 'HitTest','off');
                        hTossPatchesL{end+1} = hp; %#ok<AGROW>
                    end
                case 'right'
                    hTileBordersR{end+1} = hLine; %#ok<AGROW>
            end
        end

        % Selected tile highlight (left panel only)
        if strcmp(panel, 'left') && S.selectedTile > 0
            t = S.selectedTile;
            rr = S.place(t).rr;
            cc = S.place(t).cc;
            r1 = rr(1) - 0.5; r2 = rr(end) + 0.5;
            c1 = cc(1) - 0.5; c2 = cc(end) + 0.5;
            hSelectBorder = line(ax, [c1 c2 c2 c1 c1], [r1 r1 r2 r2 r1], ...
                'Color', 'c', 'LineWidth', 2.5, 'HitTest','off');
        end

        hold(ax, 'off');
    end

    function updateHistograms()
        updateHistLeft();
        updateHistRight();
    end

    function updateHistLeft()
        if ~S.loaded, return; end
        img = S.zCache(:,:,S.curZ);
        vals = double(img(img > 0));
        if isempty(vals), vals = double(img(:)); end
        cla(axHistLeft);
        if ~isempty(vals)
            histogram(axHistLeft, vals, 200, 'EdgeColor','none','FaceColor',[0.3 0.3 0.8]);
        end
        axHistLeft.YTick = [];
        title(axHistLeft, 'Hist (left)', 'Interpreter','none','FontSize',9);
    end

    function updateHistRight()
        if ~S.loaded, return; end
        img = S.projImg;
        vals = double(img(img > 0));
        if isempty(vals), vals = double(img(:)); end
        cla(axHistRight);
        if ~isempty(vals)
            histogram(axHistRight, vals, 200, 'EdgeColor','none','FaceColor',[0.3 0.3 0.8]);
        end
        axHistRight.YTick = [];
        title(axHistRight, 'Hist (right)', 'Interpreter','none','FontSize',9);
    end

    function refreshDisplay()
        showLeftPanel();
        showRightPanel();
        updateHistograms();
        updateTileLabel();
        updateTossLabel();
    end

    function setStatus(msg)
        lblStatus.Text = msg;
    end

    function name = extractChannelName(filepath, fallback)
        % Try to parse chN from filename or parent folder
        [~, bn] = fileparts(filepath);
        tok = regexp(char(bn), 'ch(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            name = sprintf('ch%s', tok{1});
        else
            % Try parent folder name
            parts = split(string(filepath), filesep);
            for k = numel(parts):-1:1
                tok2 = regexp(char(parts(k)), 'ch(\d+)', 'tokens', 'once');
                if ~isempty(tok2)
                    name = sprintf('ch%s', tok2{1});
                    return;
                end
            end
            name = fallback;
        end
    end

    function chID = parseChIDFromPath(filepath)
        % Parse channel ID from path: .../ch2/avgz/file.tif -> 2
        parts = split(string(filepath), filesep);
        for k = numel(parts):-1:1
            tok = regexp(char(parts(k)), '^ch(\d+)$', 'tokens', 'once');
            if ~isempty(tok)
                chID = str2double(tok{1});
                return;
            end
        end
        % Fallback: parse from filename
        [~, bn] = fileparts(filepath);
        tok = regexp(char(bn), 'ch(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            chID = str2double(tok{1});
        else
            chID = NaN;
        end
    end

    function newPaths = rebuildTilePaths(storedPaths, rootDir, chID)
        % Rebuild tile file paths relative to rootDir using channel ID.
        % For each stored path, extract basename and reconstruct under
        % rootDir/ch<chID>/avgz/<basename>
        %
        % Also handles channel substitution: if the stored path has a
        % different chN in the filename, replace it with the target chID.
        nTiles = numel(storedPaths);
        newPaths = strings(nTiles, 1);
        chFolder = sprintf('ch%d', chID);
        for ti = 1:nTiles
            [~, bn, ext] = fileparts(storedPaths(ti));
            bn = char(bn); ext = char(ext);
            % Replace channel ID in filename if different
            bn = regexprep(bn, 'ch\d+', sprintf('ch%d', chID));
            newPaths(ti) = fullfile(rootDir, chFolder, 'avgz', [bn ext]);
        end
    end

end % main function

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%% ========================================================================

function w = linear_blend_weights(H, W)
%LINEAR_BLEND_WEIGHTS  Ramp weights: distance from nearest edge + 1
[xg, yg] = meshgrid(1:W, 1:H);
dx = min(xg-1, W-xg);
dy = min(yg-1, H-yg);
w = min(dx, dy) + 1;
w = double(w);
end

function I2 = shift_bilinear(I, fx, fy)
%SHIFT_BILINEAR  Sub-pixel shift via bilinear interpolation
[H, W] = size(I);
[xg, yg] = meshgrid(1:W, 1:H);
I2 = interp2(xg, yg, double(I), xg - fx, yg - fy, 'linear', 0);
end

function [bitsPerSample, sampleFormat] = class_to_tiff_format(cls)
%CLASS_TO_TIFF_FORMAT  Map MATLAB class name to TIFF tag values
switch char(cls)
    case 'uint8',  bitsPerSample = 8;  sampleFormat = Tiff.SampleFormat.UInt;
    case 'uint16', bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.UInt;
    case 'int16',  bitsPerSample = 16; sampleFormat = Tiff.SampleFormat.Int;
    case 'single', bitsPerSample = 32; sampleFormat = Tiff.SampleFormat.IEEEFP;
    case 'double', bitsPerSample = 64; sampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise, error('stitch_qc_gui:unsupportedClass', 'Unsupported TIFF class: %s', cls);
end
end

function safe_close_tiff(tt)
%SAFE_CLOSE_TIFF  Close a Tiff handle, ignoring errors
try
    tt.close();
catch
end
end
