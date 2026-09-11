function ca_recheck_gui_260901(bundleDir)
%CA_RECHECK_GUI_260901  Per-CELL, per-RECORDING re-check of calcium event detection.
%
%   ca_recheck_gui_260901()            % uses the default bundle folder below
%   ca_recheck_gui_260901(bundleDir)   % a bundle built by ca_recheck_build_bundle_260901
%
%   THE UNIT IS THE CELL. calcium_spike_gui.m walks an animal's folders and, in
%   each, its ROIs -- so one neuron imaged in three recordings appears as three
%   unrelated ROIs in three unrelated sessions, and there is nowhere to say "this
%   cell is fine but that one recording of it is not". Here the outer loop is the
%   cell and the inner loop is that cell's recordings, which is the level the
%   judgement is actually made at.
%
%   TOSSING A RECORDING is per (cell, recording), not per recording: a folder
%   holds many cells, and a trace can be unusable for one of them while being
%   fine for the others. TOSSING THE WHOLE CELL is a separate, cell-level flag,
%   for a cell that should not enter the population at all; it overrides the
%   per-recording flags without erasing them, so un-tossing the cell gives back
%   the recording decisions that were already made. A tossed recording drops out of this cell's pooled
%   heatmap and pooled mean immediately, so the panel shows what the per-cell
%   summary figure will look like without it, and downstream (heatmap, polar)
%   picks it up through ca_recheck_load_curation_260901.
%
%   NOTHING IN THE ARCHIVE IS TOUCHED. Events and toss flags are written to
%   <bundle>\curated\<stem>_curated.mat plus two CSVs. The archive's
%   ca_spike_data.mat files are left exactly as they are; pushing the curation
%   back into them is a separate, explicit step.
%
%   DETECTION IS THE OLD GUI'S, not a re-implementation of it. It calls
%   find_spikes_threshold_260901.m, which is a verbatim copy of the local
%   function at the bottom of calcium_spike_gui.m, through the same wrapper
%   calcium_spike_gui's detectSpikesForROI uses: the same Imaging FPS /
%   threshold / min distance / min width / prominence fields, the same optional
%   linear Detrend, and the same trim-window filter on the detected frames. The
%   fields are seeded from the same saved CalciumSpikeGUI working parameters,
%   and the same click-to-add / click-to-delete editing sits on top of it.
%
%   KEYS
%     n / a / d      navigate | add event | delete event   (click in the dF/F axes)
%     left / right   previous | next recording of this cell
%     , / .          previous | next cell
%     t              toss / un-toss this recording
%     shift+t        toss / un-toss the WHOLE cell
%     r              re-detect this recording with the current parameters
%     s              save this cell
%
%   Runqi Zhang / 2026-09-01

if nargin < 1 || isempty(bundleDir)
    bundleDir = ['D:\Ventral_surface_summary\per-cell-summary_active_260812\' ...
                 'spike_recheck_260901'];
end
assert(isfolder(bundleDir), 'Bundle folder not found:\n  %s', bundleDir);

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fileparts(thisDir));

%% ---------------- shared state ----------------
S = struct();
S.bundleDir  = bundleDir;
S.cellDir    = fullfile(bundleDir,'cells');
S.curDir     = fullfile(bundleDir,'curated');
S.CL         = [];        % cell list (struct array)
S.shown      = [];        % indices of CL currently in the listbox
S.ci         = 0;         % index into S.CL of the current cell
S.C          = [];        % loaded cell: cellInfo
S.R          = [];        % loaded cell: RECc struct array
S.ri         = 1;         % current recording within the cell
S.ev         = {};        % {nRec x 1} working event frame indices
S.ev0        = {};        % {nRec x 1} events as loaded from the archive
S.toss       = false(0);  % [nRec x 1]
S.note       = strings(0);
S.src        = strings(0);% per-recording provenance of S.ev
S.reviewed   = false;
S.dirty      = false;
S.editMode   = 'none';
S.detrend    = false;     % same linear detrend toggle as the old GUI
S.cellToss   = false;     % the whole cell is thrown away

if ~isfolder(S.curDir), mkdir(S.curDir); end

%% ---------------- detection defaults ----------------
% The same four parameters calcium_spike_gui.m uses, with the same meanings and
% the same units, seeded from the same saved working parameters -- so this GUI
% opens with the detection the old one was last run with.
D.fps        = 30;      % Imaging FPS field, exactly as in the old GUI
D.threshold  = 0;       % MinPeakHeight, absolute dF/F
D.minDist    = 0.50;    % s
D.minWidth   = 0.10;    % s
D.prominence = 0;       % MinPeakProminence, absolute dF/F (0 = off)
D.snapWin    = 0.10;    % s, half-window a click snaps to the local max in
if ispref('CalciumSpikeGUI','workingParams')
    try
        wp = getpref('CalciumSpikeGUI','workingParams');
        if isfield(wp,'fps'),        D.fps        = wp.fps;        end
        if isfield(wp,'threshold'),  D.threshold  = wp.threshold;  end
        if isfield(wp,'minDist'),    D.minDist    = wp.minDist;    end
        if isfield(wp,'minWidth'),   D.minWidth   = wp.minWidth;   end
        if isfield(wp,'prominence'), D.prominence = wp.prominence; end
        if isfield(wp,'snapWin'),    D.snapWin    = wp.snapWin;    end
    catch
    end
end

%% ================= FIGURE =================
fig = uifigure('Name','Calcium re-check  --  cell / recording', ...
               'Position',[40 30 1680 960], 'WindowState','maximized', ...
               'KeyPressFcn',@cb_key);

rootGL = uigridlayout(fig,[1 2]);
rootGL.ColumnWidth   = {410,'1x'};
rootGL.Padding       = [6 6 6 6];
rootGL.ColumnSpacing = 8;

% NOTHING ON THE LEFT SCROLLS. Every panel is a fixed row of the column, so no
% control can end up below a fold -- which is how the toss button and the
% detection parameters each went missing once.
leftGL = uigridlayout(rootGL,[4 1]);
leftGL.Layout.Column = 1;
leftGL.RowHeight     = {'1x', 186, 236, 168};
leftGL.Padding       = [0 0 0 0];
leftGL.RowSpacing    = 6;

cellPan = uipanel(leftGL,'Title','Cell','FontWeight','bold','Tooltip',bundleDir);
cellPan.Layout.Row = 1;
recPan  = uipanel(leftGL,'Title','Recordings of this cell','FontWeight','bold');
recPan.Layout.Row = 2;
detPan  = uipanel(leftGL,'Title','Detection','FontWeight','bold');
detPan.Layout.Row = 3;
savePan = uipanel(leftGL,'Title','Edit / Save','FontWeight','bold');
savePan.Layout.Row = 4;
plotPan = uipanel(rootGL,'Title','','BorderType','none');
plotPan.Layout.Column = 2;

%% ---------------- left, fixed: cell ----------------
cGL = uigridlayout(cellPan,[7 2]);
cGL.RowHeight     = {22, 20, '1x', 24, 30, 26, 16};
cGL.ColumnWidth   = {60,'1x'};
cGL.Padding       = [8 4 8 4];
cGL.RowSpacing    = 3;
cGL.ColumnSpacing = 5;

h=uilabel(cGL,'Text','Group:'); h.Layout.Row=1; h.Layout.Column=1;
ddGroup = uidropdown(cGL,'Items',{'all'},'ValueChangedFcn',@(~,~) refreshCellList());
ddGroup.Layout.Row=1; ddGroup.Layout.Column=2;
cbTodo = uicheckbox(cGL,'Text','only cells not yet reviewed', ...
         'ValueChangedFcn',@(~,~) refreshCellList());
cbTodo.Layout.Row=2; cbTodo.Layout.Column=[1 2];
lbCell = uilistbox(cGL,'Items',{'(loading)'},'ValueChangedFcn',@cb_pickCell, ...
         'FontName','Consolas','FontSize',10);
lbCell.Layout.Row=3; lbCell.Layout.Column=[1 2];
navC = uigridlayout(cGL,[1 2]); navC.Layout.Row=4; navC.Layout.Column=[1 2];
navC.Padding=[0 0 0 0]; navC.ColumnSpacing=5; navC.ColumnWidth={'1x','1x'};
bPrevC = uibutton(navC,'Text','<- Prev cell  (,)','ButtonPushedFcn',@(~,~) stepCell(-1));
bNextC = uibutton(navC,'Text','Next cell (.) ->','ButtonPushedFcn',@(~,~) stepCell(+1));
lblCell = uilabel(cGL,'Text','---','FontWeight','bold','WordWrap','on','FontSize',10);
lblCell.Layout.Row=5; lblCell.Layout.Column=[1 2];
btnTossCell = uibutton(cGL,'state','Text','TOSS ENTIRE CELL  (shift+T)', ...
              'ValueChangedFcn',@cb_tossCell);
btnTossCell.Layout.Row=6; btnTossCell.Layout.Column=[1 2];
lblProgress = uilabel(cGL,'Text','0 / 0 cells reviewed','FontSize',9, ...
              'FontColor',[0.20 0.55 0.20]);
lblProgress.Layout.Row=7; lblProgress.Layout.Column=[1 2];

%% ---------------- left, fixed: recordings + toss ----------------
rGL = uigridlayout(recPan,[4 2]);
rGL.RowHeight     = {'1x', 24, 26, 24};
rGL.ColumnWidth   = {60,'1x'};
rGL.Padding       = [8 4 8 4];
rGL.RowSpacing    = 3;
rGL.ColumnSpacing = 5;

lbRec = uilistbox(rGL,'Items',{'---'},'ValueChangedFcn',@cb_pickRec, ...
        'FontName','Consolas','FontSize',10);
lbRec.Layout.Row=1; lbRec.Layout.Column=[1 2];
navR = uigridlayout(rGL,[1 2]); navR.Layout.Row=2; navR.Layout.Column=[1 2];
navR.Padding=[0 0 0 0]; navR.ColumnSpacing=5; navR.ColumnWidth={'1x','1x'};
bPrevR = uibutton(navR,'Text','<- Prev rec','ButtonPushedFcn',@(~,~) stepRec(-1));
bNextR = uibutton(navR,'Text','Next rec ->','ButtonPushedFcn',@(~,~) stepRec(+1));
btnToss = uibutton(rGL,'state','Text','TOSS this recording  (t)','ValueChangedFcn',@cb_toss);
btnToss.Layout.Row=3; btnToss.Layout.Column=[1 2];
h=uilabel(rGL,'Text','Note:'); h.Layout.Row=4; h.Layout.Column=1;
efNote = uieditfield(rGL,'text','Value','','ValueChangedFcn',@cb_note);
efNote.Layout.Row=4; efNote.Layout.Column=2;

%% ---------------- left, fixed: detection ----------------
% Every parameter that shapes a detection lives here, in a panel that does not
% scroll. In the scrolling list, min width and prominence sat below the fold --
% a parameter you cannot see silently deciding the result.
dGL = uigridlayout(detPan,[8 4]);
dGL.RowHeight     = repmat({22},1,8);
dGL.ColumnWidth   = {88,'1x',92,'1x'};
dGL.Padding       = [8 4 8 4];
dGL.RowSpacing    = 3;
dGL.ColumnSpacing = 5;

h=uilabel(dGL,'Text','Imaging FPS:');   h.Layout.Row=1; h.Layout.Column=1;
efFps = uieditfield(dGL,'numeric','Value',D.fps,'Limits',[1 2000], ...
        'ValueChangedFcn',@(~,~) refreshTrace());
efFps.Layout.Row=1; efFps.Layout.Column=2;
h=uilabel(dGL,'Text','Threshold dF/F:'); h.Layout.Row=1; h.Layout.Column=3;
efThr = uieditfield(dGL,'numeric','Value',D.threshold,'Limits',[-Inf Inf], ...
        'ValueChangedFcn',@(~,~) refreshTrace());
efThr.Layout.Row=1; efThr.Layout.Column=4;

h=uilabel(dGL,'Text','Min dist (s):');  h.Layout.Row=2; h.Layout.Column=1;
efDist = uieditfield(dGL,'numeric','Value',D.minDist,'Limits',[0 Inf]);
efDist.Layout.Row=2; efDist.Layout.Column=2;
h=uilabel(dGL,'Text','Min width (s):'); h.Layout.Row=2; h.Layout.Column=3;
efWid = uieditfield(dGL,'numeric','Value',D.minWidth,'Limits',[0 Inf]);
efWid.Layout.Row=2; efWid.Layout.Column=4;

h=uilabel(dGL,'Text','Prominence:');    h.Layout.Row=3; h.Layout.Column=1;
efProm = uieditfield(dGL,'numeric','Value',D.prominence,'Limits',[0 Inf]);
efProm.Layout.Row=3; efProm.Layout.Column=2;
btnDetrend = uibutton(dGL,'Text','Detrend: OFF','ButtonPushedFcn',@cb_toggleDetrend, ...
             'BackgroundColor',[0.25 0.25 0.55],'FontColor','white');
btnDetrend.Layout.Row=3; btnDetrend.Layout.Column=[3 4];

h=uilabel(dGL,'Text','Trim start (s):'); h.Layout.Row=4; h.Layout.Column=1;
efTrimStart = uieditfield(dGL,'numeric','Value',0,'Limits',[0 Inf]);
efTrimStart.Layout.Row=4; efTrimStart.Layout.Column=2;
h=uilabel(dGL,'Text','Trim end (s):');   h.Layout.Row=4; h.Layout.Column=3;
efTrimEnd = uieditfield(dGL,'numeric','Value',Inf,'Limits',[0 Inf]);
efTrimEnd.Layout.Row=4; efTrimEnd.Layout.Column=4;

lblInfo = uilabel(dGL,'Text','','FontSize',9,'FontColor',[0.35 0.35 0.35]);
lblInfo.Layout.Row=5; lblInfo.Layout.Column=[1 4];

% @(~,~) detectRec(S.ri,...) would freeze S.ri at its value here (1) -- an
% anonymous function captures workspace variables BY VALUE when it is created.
% A nested function reads the live S instead.
bDet1 = uibutton(dGL,'Text','Detect this rec (r)','ButtonPushedFcn',@cb_detectThis);
bDet1.Layout.Row=6; bDet1.Layout.Column=[1 2];
bDetC = uibutton(dGL,'Text','Detect whole cell','ButtonPushedFcn',@cb_detectCell);
bDetC.Layout.Row=6; bDetC.Layout.Column=[3 4];
bRestore = uibutton(dGL,'Text','Restore archive events','ButtonPushedFcn',@cb_restore);
bRestore.Layout.Row=7; bRestore.Layout.Column=[1 2];
bDetAll = uibutton(dGL,'Text','Batch-detect ALL cells...','ButtonPushedFcn',@cb_detectAll, ...
          'FontColor',[0.60 0.25 0.05]);
bDetAll.Layout.Row=7; bDetAll.Layout.Column=[3 4];
lblCount = uilabel(dGL,'Text','events: 0','FontSize',10);
lblCount.Layout.Row=8; lblCount.Layout.Column=[1 4];

%% ---------------- left, fixed: edit + save bar ----------------
% Editing and saving live OUTSIDE the scrolling control list, because in the
% scroll pane they sat below the fold and were effectively invisible.
sGL = uigridlayout(savePan,[6 2]);
sGL.RowHeight   = repmat({'fit'},1,6);
sGL.ColumnWidth = {'1x','1x'};
sGL.Padding     = [8 6 8 6];
sGL.RowSpacing  = 4;

lblMode = uilabel(sGL,'Text','Mode: navigate','FontWeight','bold');
lblMode.Layout.Row = 1; lblMode.Layout.Column = [1 2];
bNav = uibutton(sGL,'Text','Navigate (n)','ButtonPushedFcn',@(~,~) setMode('none'));
bNav.Layout.Row = 2; bNav.Layout.Column = 1;
bAdd = uibutton(sGL,'Text','+ Add (a)','ButtonPushedFcn',@(~,~) setMode('add'));
bAdd.Layout.Row = 2; bAdd.Layout.Column = 2;
bDel = uibutton(sGL,'Text','- Delete (d)','ButtonPushedFcn',@(~,~) setMode('delete'));
bDel.Layout.Row = 3; bDel.Layout.Column = 1;
bClr = uibutton(sGL,'Text','Clear this rec','ButtonPushedFcn',@cb_clearRec);
bClr.Layout.Row = 3; bClr.Layout.Column = 2;

cbReviewed = uicheckbox(sGL,'Text','cell reviewed','ValueChangedFcn',@cb_reviewed);
cbReviewed.Layout.Row = 4; cbReviewed.Layout.Column = [1 2];
bSave = uibutton(sGL,'Text','Save cell (s)','ButtonPushedFcn',@(~,~) saveCell(true), ...
                 'BackgroundColor',[0.20 0.55 0.25],'FontColor','white','FontWeight','bold');
bSave.Layout.Row = 5; bSave.Layout.Column = 1;
bCsv = uibutton(sGL,'Text','Write CSVs','ButtonPushedFcn',@(~,~) writeCsvs(true));
bCsv.Layout.Row = 5; bCsv.Layout.Column = 2;
lblSave = uilabel(sGL,'Text','','FontSize',9,'FontColor',[0.45 0.45 0.45],'WordWrap','on');
lblSave.Layout.Row = 6; lblSave.Layout.Column = [1 2];

%% ---------------- right: plots ----------------
% Same three panels, same proportions and same styling as calcium_spike_gui.m:
% avg projection crop on top, dF/F trace in the middle, raw F at the bottom.
pGL = uigridlayout(plotPan,[3 1]);
pGL.RowHeight  = {'0.22x','0.42x','0.36x'};
pGL.Padding    = [4 4 4 4];
pGL.RowSpacing = 6;

axSum = uiaxes(pGL);
axSum.Layout.Row        = 1;
axSum.Title.String      = 'Avg Projection';
axSum.Title.Interpreter = 'none';
axis(axSum,'image');

axDff = uiaxes(pGL);
axDff.Layout.Row         = 2;
axDff.Title.String       = ['dF/F Trace  |  Left-click: add/delete events  |  ' ...
                            'N: navigate mode'];
axDff.Title.Interpreter  = 'none';
axDff.XLabel.String      = 'Time (s)';
axDff.YLabel.String      = 'dF/F';
hold(axDff,'on'); grid(axDff,'on'); box(axDff,'on');

axRawF = uiaxes(pGL);
axRawF.Layout.Row        = 3;
axRawF.Title.String      = 'Raw F';
axRawF.Title.Interpreter = 'none';
axRawF.XLabel.String     = 'Time (s)';
axRawF.YLabel.String     = 'F (mean pixel intensity)';
hold(axRawF,'on'); grid(axRawF,'on'); box(axRawF,'on');

hDffLine = gobjects(0); hSpk = gobjects(0); hThrLine = gobjects(0);

%% ---------------- go ----------------
scanBundle();
if ~isempty(S.CL), selectCell(1); end

%% =====================================================================
%%  BUNDLE / CELL LIST
%% =====================================================================
    function scanBundle()
        d = dir(fullfile(S.cellDir,'*.mat'));
        assert(~isempty(d), 'No cell files in %s -- run ca_recheck_build_bundle_260901 first', S.cellDir);
        CL = struct('stem',{},'file',{},'group',{},'cell',{},'nRec',{}, ...
                    'reviewed',{},'nToss',{},'nEv',{},'cellToss',{},'curFile',{});
        for k = 1:numel(d)
            stem = erase(d(k).name,'.mat');
            parts = regexp(stem,'_','split');
            g = parts{1};
            cn = str2double(regexprep(stem,'.*cell(\d+)$','$1'));
            CL(k) = struct('stem',stem, 'file',fullfile(d(k).folder,d(k).name), ...
                           'group',g, 'cell',cn, 'nRec',NaN, 'reviewed',false, ...
                           'nToss',0, 'nEv',NaN, 'cellToss',false, ...
                           'curFile',fullfile(S.curDir,[stem '_curated.mat']));
        end
        % fold in whatever curation already exists, so the list shows state
        for k = 1:numel(CL)
            if isfile(CL(k).curFile)
                q = load(CL(k).curFile,'CUR');
                if isfield(q,'CUR')
                    CL(k).reviewed = logical(q.CUR.reviewed);
                    CL(k).nToss    = nnz([q.CUR.rec.toss]);
                    CL(k).nEv      = sum([q.CUR.rec.n_events]);
                    CL(k).nRec     = numel(q.CUR.rec);
                    if isfield(q.CUR,'cell_toss')
                        CL(k).cellToss = logical(q.CUR.cell_toss);
                    end
                end
            end
        end
        [~,ord] = sort([CL.cell]);  S.CL = CL(ord);
        ddGroup.Items = [{'all'}, unique({S.CL.group},'stable')];
        refreshCellList();
    end

    function refreshCellList()
        keep = true(numel(S.CL),1);
        if ~strcmp(ddGroup.Value,'all')
            keep = keep & strcmp({S.CL.group}', ddGroup.Value);
        end
        if cbTodo.Value
            keep = keep & ~[S.CL.reviewed]';
        end
        S.shown = find(keep);
        if isempty(S.shown)
            lbCell.Items = {'(nothing matches)'}; lbCell.ItemsData = 0;
        else
            items = cell(numel(S.shown),1);
            for q = 1:numel(S.shown)
                items{q} = cellLabel(S.CL(S.shown(q)));
            end
            lbCell.Items = items; lbCell.ItemsData = S.shown;
            if S.ci > 0 && ismember(S.ci, S.shown), lbCell.Value = S.ci; end
        end
        nRev = nnz([S.CL.reviewed]);
        lblProgress.Text = sprintf('%d / %d reviewed  |  %d tossed recs  |  %d tossed cells', ...
            nRev, numel(S.CL), sum([S.CL.nToss]), nnz([S.CL.cellToss]));
    end

    function s = cellLabel(c)
        mark = '  ';
        if c.reviewed, mark = 'x '; end
        if c.nToss > 0, mark(2) = 'T'; end
        if c.cellToss, mark = 'XX'; end
        if isnan(c.nEv), s = sprintf('%s%-24s', mark, c.stem);
        else,            s = sprintf('%s%-24s %3d ev', mark, c.stem, c.nEv);
        end
    end

    function cb_pickCell(~,~)
        v = lbCell.Value;
        if isempty(v) || v == 0, return; end
        selectCell(v);
    end

    function stepCell(step)
        if isempty(S.shown), return; end
        p = find(S.shown == S.ci, 1);
        if isempty(p), p = 1; else, p = min(max(p+step,1), numel(S.shown)); end
        selectCell(S.shown(p));
    end

%% =====================================================================
%%  LOAD ONE CELL
%% =====================================================================
    function selectCell(ci)
        if S.dirty, saveCell(false); end
        S.ci = ci;
        L = load(S.CL(ci).file,'cellInfo','RECc');
        S.C = L.cellInfo;  S.R = L.RECc;
        n  = numel(S.R);
        S.ev0 = arrayfun(@(x) x.spike_idx(:), S.R, 'UniformOutput',false);
        S.ev  = S.ev0;
        S.toss = false(n,1);
        S.note = strings(n,1);
        S.src  = repmat("archive", n, 1);
        S.reviewed = false;
        S.cellToss = false;

        % curation, if any, wins over the archive snapshot
        cf = S.CL(ci).curFile;
        if isfile(cf)
            q = load(cf,'CUR');
            if isfield(q,'CUR')
                CUR = q.CUR;
                for i = 1:n
                    m = find(strcmp({CUR.rec.folder}, S.R(i).folder) & ...
                             [CUR.rec.roi] == S.R(i).roi, 1);
                    if isempty(m), continue; end
                    S.ev{i}   = CUR.rec(m).spike_idx(:);
                    S.toss(i) = logical(CUR.rec(m).toss);
                    S.note(i) = string(CUR.rec(m).note);
                    S.src(i)  = string(CUR.rec(m).source);
                end
                S.reviewed = logical(CUR.reviewed);
                if isfield(CUR,'cell_toss'), S.cellToss = logical(CUR.cell_toss); end
            end
        end
        S.CL(ci).nRec = n;
        S.dirty = false;
        S.ri = 1;
        cbReviewed.Value  = S.reviewed;
        btnTossCell.Value = S.cellToss;
        styleTossCell();
        lblCell.Text = sprintf('%s  |  %d rec  |  %.1f ev/min pooled\n%s', ...
            S.C.stem, n, S.C.rate_min, S.C.key);
        if ismember(ci, S.shown), lbCell.Value = ci; end
        refreshRecList();
        selectRec(1);
        lblSave.Text = '';
    end

    function refreshRecList()
        items = cell(numel(S.R),1);
        for i = 1:numel(S.R)
            tag = '  ';
            if S.toss(i), tag = 'T '; end
            items{i} = sprintf('%s%d) roi%-3d %5d ev  %5.1fs  %s', tag, i, ...
                S.R(i).roi, numel(S.ev{i}), S.R(i).dur_s, shortName(S.R(i).recName));
        end
        lbRec.Items = items; lbRec.ItemsData = 1:numel(S.R);
        lbRec.Value = S.ri;
    end

    function s = shortName(nm)
        s = char(nm);
        if numel(s) > 34, s = ['...' s(end-30:end)]; end
    end

    function cb_pickRec(~,~)
        v = lbRec.Value; if isempty(v), return; end
        selectRec(v);
    end

    function stepRec(step)
        selectRec(min(max(S.ri+step,1), numel(S.R)));
    end

    function selectRec(i)
        S.ri = i;
        lbRec.Value = i;
        btnToss.Value = S.toss(i);
        styleToss();
        efNote.Value = char(S.note(i));
        refreshTrace();
    end

%% =====================================================================
%%  TRACE PANELS
%% =====================================================================
    function x = workTrace(i)
        % The trace detection and display both run on -- S.dFF in the old GUI,
        % which is the raw dF/F unless the Detrend toggle is on, in which case it
        % is detrend(dFF). detrend() is column-wise, so detrending this one column
        % gives the same numbers as detrending the whole matrix and taking it.
        x = S.R(i).dff;
        if S.detrend, x = detrend(x); end
    end

    function cb_toggleDetrend(~,~)
        S.detrend = ~S.detrend;
        if S.detrend
            btnDetrend.Text = 'Detrend: ON';
            btnDetrend.BackgroundColor = [0.10 0.60 0.45];
        else
            btnDetrend.Text = 'Detrend: OFF';
            btnDetrend.BackgroundColor = [0.25 0.25 0.55];
        end
        refreshTrace();
    end

    function refreshTrace()
        cla(axDff);
        if isempty(S.R), return; end
        r = S.R(S.ri);
        x = workTrace(S.ri);
        t = (0:r.T-1)'/efFps.Value;

        hDffLine = plot(axDff, t, x, 'Color',[0 0 0], 'LineWidth',0.9);
        hDffLine.HitTest = 'off';
        hold(axDff,'on');
        thr = efThr.Value;
        hThrLine = yline(axDff, thr, '--', 'Color',[0.85 0.15 0.15], 'LineWidth',0.8, ...
                         'Label','threshold', 'LabelHorizontalAlignment','left');
        hThrLine.HitTest = 'off';
        hSpk = scatter(axDff, [], [], 75, 'v', 'filled', ...
            'MarkerFaceColor',[1.00 0.28 0.28], 'MarkerEdgeColor','none', ...
            'Tag','spike_marker');
        hSpk.HitTest = 'off';

        xlim(axDff, [t(1) t(end)]);
        ylo = min(x); yhi = max(x);
        if yhi <= ylo, yhi = ylo + 1; end
        pad = (yhi - ylo) / 0.6 - (yhi - ylo);      % 20% padding each side
        ylim(axDff, [ylo - pad/2, yhi + pad/2]);

        ttl = sprintf('%s   rec %d/%d   roi %d   --   %s', S.C.stem, S.ri, ...
                      numel(S.R), r.roi, r.recName);
        if S.toss(S.ri), ttl = ['[TOSSED]  ' ttl]; end
        if S.cellToss,   ttl = ['[CELL TOSSED]  ' ttl]; end
        axDff.Title.String = ttl;
        axDff.Title.Color  = ternary(S.toss(S.ri) || S.cellToss, [0.75 0.15 0.15], [0 0 0]);
        lblInfo.Text = sprintf('%d frames  |  file fps %.2f%s', r.T, r.fps, ...
            ternary(abs(r.fps - efFps.Value) > 0.01, ...
                    sprintf('  (field says %.2f)', efFps.Value), ''));
        setMode(S.editMode);
        redrawSpikes();
        refreshRawF();
        refreshCrop();
    end

    function refreshRawF()
        cla(axRawF);
        if isempty(S.R), return; end
        r = S.R(S.ri);
        if all(r.F_roi == 0)
            axRawF.Title.String = 'Raw F  (not available)';
            return
        end
        t = (0:r.T-1)'/efFps.Value;
        h = plot(axRawF, t, r.F_roi, 'Color',[0 0 0], 'LineWidth',0.9);
        h.HitTest = 'off';
        xlim(axRawF, [t(1) t(end)]);
        flo = min(r.F_roi); fhi = max(r.F_roi);
        if fhi <= flo, fhi = flo + 1; end
        fpad = (fhi - flo)/0.6 - (fhi - flo);
        ylim(axRawF, [flo - fpad/2, fhi + fpad/2]);
        axRawF.Title.String = sprintf('Raw F  |  ROI %d', r.roi);
    end

    function refreshCrop()
        % 30 um box around the ROI centroid off the *_AVG_for_CP.tif -- the frame
        % Cellpose actually segmented, so it matches maskL by construction. Same
        % percentile clip, gamma 0.6 and burned-in yellow outline as the old GUI.
        cla(axSum); axSum.Color = [0 0 0];
        axSum.XTick = []; axSum.YTick = [];
        if isempty(S.R), return; end
        r = S.R(S.ri);
        [maskL, avgimg] = cropSource(r.folder);
        if isempty(maskL)
            axSum.Title.String = sprintf('Avg Proj  |  ROI %d  (no mask)', r.roi);
            return
        end
        labels = unique(maskL(:)); labels(labels==0) = [];
        if r.roi > numel(labels)
            axSum.Title.String = sprintf('Avg Proj  |  ROI %d  (out of range)', r.roi);
            return
        end
        roiMask = (maskL == labels(r.roi));
        pr = regionprops(roiMask,'Centroid');
        if isempty(pr)
            axSum.Title.String = sprintf('Avg Proj  |  ROI %d  (empty mask)', r.roi);
            return
        end
        if isempty(avgimg), avgimg = double(maskL>0)*100; end
        cx = round(pr(1).Centroid(1)); cy = round(pr(1).Centroid(2));
        [H,W] = size(avgimg);
        px_um = r.px_um; if ~isfinite(px_um) || px_um<=0, px_um = 1.7778/2; end
        half  = max(4, round(30/px_um/2));          % 30 um box
        box   = 2*half + 1;
        r1 = max(1, cy-half); r2 = min(H, cy+half);
        c1 = max(1, cx-half); c2 = min(W, cx+half);
        crop  = zeros(box,box); cmask = false(box,box);
        pr1 = half-(cy-r1)+1;  pc1 = half-(cx-c1)+1;
        crop( pr1:pr1+(r2-r1), pc1:pc1+(c2-c1)) = avgimg(r1:r2, c1:c2);
        cmask(pr1:pr1+(r2-r1), pc1:pc1+(c2-c1)) = roiMask(r1:r2, c1:c2);
        v = crop(crop>0); if isempty(v), v = crop(:); end
        lo = prctile(v,0.5); hi = prctile(v,99.5);
        if hi <= lo, hi = lo + 1; end
        crop = max(0, min(1, (crop-lo)/(hi-lo))) .^ 0.6;
        rgb = repmat(crop,[1 1 3]);
        per = bwperim(cmask);
        for ch = 1:3
            pl = rgb(:,:,ch); pl(per) = double(ch<3); rgb(:,:,ch) = pl;   % yellow
        end
        image(axSum, rgb); axis(axSum,'image');
        axSum.XTick = []; axSum.YTick = [];
        axSum.Title.String = sprintf('Avg Proj  |  ROI %d', r.roi);
    end

    function [maskL, avgimg] = cropSource(folder)
        % One folder's mask and projection, cached: switching between recordings
        % of a cell would otherwise re-read the same TIFF every time.
        maskL = []; avgimg = [];
        if isfield(S,'cropCache') && ~isempty(S.cropCache) && ...
                strcmp(S.cropCache.folder, folder)
            maskL = S.cropCache.maskL; avgimg = S.cropCache.avgimg; return
        end
        sam = dir(fullfile(folder,'*_cpSAM_output.mat'));
        if ~isempty(sam)
            q = load(fullfile(sam(1).folder, sam(1).name),'maskL');
            if isfield(q,'maskL'), maskL = q.maskL; end
        end
        av = dir(fullfile(folder,'*_AVG_for_CP.tif'));
        if ~isempty(av)
            try, avgimg = double(imread(fullfile(av(1).folder, av(1).name))); catch, end
        end
        if isempty(avgimg)
            av = dir(fullfile(folder,'AVG_*_MC_MC.tif'));
            if ~isempty(av)
                try, avgimg = double(imread(fullfile(av(1).folder, av(1).name))); catch, end
            end
        end
        if ~isempty(maskL) && ~isempty(avgimg) && ~isequal(size(avgimg), size(maskL))
            avgimg = [];
        end
        S.cropCache = struct('folder',folder,'maskL',maskL,'avgimg',avgimg);
    end

    function redrawSpikes()
        if isempty(S.R) || ~isgraphics(hSpk), return; end
        r = S.R(S.ri);
        x = workTrace(S.ri);
        e = S.ev{S.ri};
        e = e(e>=1 & e<=r.T);
        if isempty(e)
            set(hSpk,'XData',[],'YData',[]);
        else
            set(hSpk,'XData',(e-1)/efFps.Value,'YData',x(e));
        end
        lblCount.Text = sprintf('events: %d   (%.1f/min, source: %s)', numel(e), ...
            numel(e)/max(r.dur_s,eps)*60, S.src(S.ri));
        items = lbRec.Items;
        if numel(items) >= S.ri
            tag = '  '; if S.toss(S.ri), tag = 'T '; end
            items{S.ri} = sprintf('%s%d) roi%-3d %5d ev  %5.1fs  %s', tag, S.ri, ...
                r.roi, numel(e), r.dur_s, shortName(r.recName));
            lbRec.Items = items; lbRec.Value = S.ri;
        end
    end


%% =====================================================================
%%  DETECTION
%% =====================================================================



    function detectRec(i, doDraw)
        % calcium_spike_gui.m's detectSpikesForROI, step for step: read the five
        % fields, run find_spikes_threshold on the dF/F trace, then drop anything
        % outside the trim window.
        if isempty(S.R), return; end
        fps     = efFps.Value;
        thr     = efThr.Value;
        minDist = efDist.Value;
        minW    = efWid.Value;
        prom    = efProm.Value;
        trace_det = workTrace(i);

        spk = find_spikes_threshold_260901(trace_det, fps, thr, minDist, minW, prom);

        tStart = efTrimStart.Value;
        tEnd   = efTrimEnd.Value;
        t_img  = (0:numel(trace_det)-1)'/fps;
        if ~isempty(spk)
            spk_t = t_img(min(spk, numel(t_img)));
            keep  = spk_t >= tStart & spk_t <= tEnd;
            spk   = spk(keep);
        end

        S.ev{i}  = spk(:);
        S.src(i) = "redetected";
        S.dirty  = true;
        if doDraw, redrawSpikes(); end
    end

    function cb_detectThis(~,~)
        detectRec(S.ri, true);
    end

    function cb_detectCell(~,~)
        for i = 1:numel(S.R), detectRec(i, false); end
        refreshRecList(); redrawSpikes();
    end

    function cb_restore(~,~)
        S.ev = S.ev0;  S.src(:) = "archive";  S.dirty = true;
        refreshRecList(); redrawSpikes();
    end

    function cb_detectAll(~,~)
        ch = uiconfirm(fig, sprintf(['Re-detect events in EVERY cell of the bundle ' ...
            'with the current parameters?\n\nfps %g,  threshold %.4g dF/F,  ' ...
            'min dist %.2f s,  min width %.2f s,  prominence %.4g,  detrend %s' ...
            '\n\nThis writes a curated file for each cell. The archive is not ' ...
            'touched.'], efFps.Value, efThr.Value, efDist.Value, efWid.Value, ...
            efProm.Value, ternary(S.detrend,'ON','OFF')), ...
            'Batch detect', 'Options',{'All cells','Only cells not yet reviewed','Cancel'}, ...
            'DefaultOption',3, 'CancelOption',3);
        if strcmp(ch,'Cancel'), return; end
        onlyTodo = strcmp(ch,'Only cells not yet reviewed');
        if S.dirty, saveCell(false); end

        dlg = uiprogressdlg(fig,'Title','Batch detect','Message','...','Cancelable','on');
        tgt = 1:numel(S.CL);
        if onlyTodo, tgt = tgt(~[S.CL.reviewed]); end
        keepCi = S.ci;
        for q = 1:numel(tgt)
            if dlg.CancelRequested, break; end
            k = tgt(q);
            dlg.Value = q/numel(tgt);
            dlg.Message = sprintf('%d / %d   %s', q, numel(tgt), S.CL(k).stem);
            L = load(S.CL(k).file,'cellInfo','RECc');
            Rk = L.RECc; nk = numel(Rk);
            evk = cell(nk,1); srck = repmat("redetected", nk, 1);
            tossk = false(nk,1); notek = strings(nk,1); ctossk = false;
            if isfile(S.CL(k).curFile)      % keep existing toss flags and notes
                qq = load(S.CL(k).curFile,'CUR');
                if isfield(qq,'CUR')
                    if isfield(qq.CUR,'cell_toss'), ctossk = logical(qq.CUR.cell_toss); end
                    for i = 1:nk
                        m = find(strcmp({qq.CUR.rec.folder}, Rk(i).folder) & ...
                                 [qq.CUR.rec.roi] == Rk(i).roi, 1);
                        if ~isempty(m)
                            tossk(i) = logical(qq.CUR.rec(m).toss);
                            notek(i) = string(qq.CUR.rec(m).note);
                        end
                    end
                end
            end
            for i = 1:nk
                xk = Rk(i).dff;
                if S.detrend, xk = detrend(xk); end
                sp = find_spikes_threshold_260901(xk, efFps.Value, efThr.Value, ...
                                                  efDist.Value, efWid.Value, efProm.Value);
                tk = (0:numel(xk)-1)'/efFps.Value;
                if ~isempty(sp)
                    spt = tk(min(sp, numel(tk)));
                    sp  = sp(spt >= efTrimStart.Value & spt <= efTrimEnd.Value);
                end
                evk{i} = sp(:);
            end
            writeCurated(S.CL(k), L.cellInfo, Rk, evk, tossk, notek, srck, ...
                         S.CL(k).reviewed, ctossk);
            S.CL(k).nEv   = sum(cellfun(@numel, evk));
            S.CL(k).nToss    = nnz(tossk);
            S.CL(k).cellToss = ctossk;
            S.CL(k).nRec  = nk;
        end
        close(dlg);
        refreshCellList();
        writeCsvs(false);
        selectCell(keepCi);
        lblSave.Text = sprintf('batch detect done (%d cells)', numel(tgt));
    end

%% =====================================================================
%%  EDITING
%% =====================================================================
    function setMode(m)
        S.editMode = m;
        switch m
            case 'add'
                lblMode.Text = 'Mode: + ADD   (click snaps to local max)';
                lblMode.FontColor = [0.20 0.45 0.90];
                disableDefaultInteractivity(axDff); axDff.ButtonDownFcn = @cb_click;
            case 'delete'
                lblMode.Text = 'Mode: - DELETE   (click removes nearest)';
                lblMode.FontColor = [0.85 0.20 0.20];
                disableDefaultInteractivity(axDff); axDff.ButtonDownFcn = @cb_click;
            otherwise
                lblMode.Text = 'Mode: navigate   (zoom / pan)';
                lblMode.FontColor = [0.25 0.25 0.25];
                enableDefaultInteractivity(axDff); axDff.ButtonDownFcn = '';
        end
    end

    function cb_click(~,~)
        % Same click behaviour as calcium_spike_gui.m: add snaps to the local max
        % within snapWin and refuses a duplicate within 2 frames, delete removes
        % the nearest event.
        if isempty(S.R), return; end
        r  = S.R(S.ri);
        x  = workTrace(S.ri);
        cp = axDff.CurrentPoint;
        fr = round(cp(1,1)*efFps.Value) + 1;
        if fr < 1 || fr > r.T, return; end
        e = S.ev{S.ri};
        switch S.editMode
            case 'add'
                sn = max(1, round(D.snapWin*efFps.Value));
                lo = max(1, fr-sn); hi = min(r.T, fr+sn);
                [~,im] = max(x(lo:hi));
                q = lo + im - 1;
                if isempty(e) || min(abs(e-q)) > 2, e = sort([e(:); q]); end
            case 'delete'
                if isempty(e), return; end
                [~,ci] = min(abs(e-fr)); e(ci) = [];
            otherwise
                return
        end
        S.ev{S.ri} = e; S.src(S.ri) = "edited"; S.dirty = true;
        redrawSpikes();
    end

    function cb_clearRec(~,~)
        S.ev{S.ri} = zeros(0,1);  S.src(S.ri) = "cleared";  S.dirty = true;
        redrawSpikes();
    end

%% =====================================================================
%%  TOSS / NOTE / KEYS
%% =====================================================================
    function cb_toss(~,~)
        S.toss(S.ri) = btnToss.Value;
        S.dirty = true;
        styleToss(); refreshRecList(); refreshTrace();    end

    function cb_tossCell(~,~)
        % A cell-level throw-away, kept SEPARATE from the per-recording flags:
        % un-tossing the cell must give back exactly the per-recording decisions
        % that were made before, not clear them. Downstream, a recording counts
        % as kept only if neither its own flag nor its cell's is set.
        S.cellToss = btnTossCell.Value;
        S.dirty = true;
        styleTossCell(); refreshCellList(); refreshTrace();
    end

    function styleTossCell()
        if btnTossCell.Value
            btnTossCell.BackgroundColor = [0.75 0.15 0.15];
            btnTossCell.FontColor = [1 1 1];
            btnTossCell.Text = 'CELL TOSSED -- click to keep  (shift+T)';
        else
            btnTossCell.BackgroundColor = [0.94 0.94 0.94];
            btnTossCell.FontColor = [0 0 0];
            btnTossCell.Text = 'TOSS ENTIRE CELL  (shift+T)';
        end
    end

    function styleToss()
        if btnToss.Value
            btnToss.BackgroundColor = [0.90 0.35 0.35];
            btnToss.Text = 'TOSSED -- click to keep  (t)';
        else
            btnToss.BackgroundColor = [0.94 0.94 0.94];
            btnToss.Text = 'TOSS this recording  (t)';
        end
    end

    function cb_note(~,~)
        S.note(S.ri) = string(efNote.Value); S.dirty = true;
    end

    function cb_reviewed(~,~)
        S.reviewed = cbReviewed.Value; S.dirty = true;
    end

    function cb_key(~, ev)
        switch lower(ev.Key)
            case 'n', setMode('none');
            case 'a', setMode('add');
            case 'd', setMode('delete');
            case 't'
                if any(strcmpi(ev.Modifier,'shift'))
                    btnTossCell.Value = ~btnTossCell.Value; cb_tossCell();
                else
                    btnToss.Value = ~btnToss.Value; cb_toss();
                end
            case 'r', detectRec(S.ri, true);
            case 's', saveCell(true);
            case 'rightarrow', stepRec(+1);
            case 'leftarrow',  stepRec(-1);
            case 'comma',  stepCell(-1);
            case 'period', stepCell(+1);
        end
    end

%% =====================================================================
%%  SAVING
%% =====================================================================
    function saveCell(loud)
        if isempty(S.C), return; end
        writeCurated(S.CL(S.ci), S.C, S.R, S.ev, S.toss, S.note, S.src, S.reviewed, ...
                     S.cellToss);
        S.CL(S.ci).reviewed = S.reviewed;
        S.CL(S.ci).cellToss = S.cellToss;
        S.CL(S.ci).nEv      = sum(cellfun(@numel, S.ev));
        S.CL(S.ci).nToss    = nnz(S.toss);
        S.CL(S.ci).nRec     = numel(S.R);
        S.dirty = false;
        refreshCellList();
        if loud
            lblSave.Text = sprintf('saved %s  (%d ev, %d tossed)  %s', S.C.stem, ...
                S.CL(S.ci).nEv, S.CL(S.ci).nToss, datestr(now,'HH:MM:SS')); %#ok<TNOW1,DATST>
        end
    end

    function writeCurated(CLk, cellInfo, R, ev, toss, note, src, reviewed, cellToss)
        rec = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{}, ...
                     'fps',{},'T',{},'dur_s',{},'spike_idx',{},'n_events',{}, ...
                     'rate_per_min',{},'toss',{},'note',{},'source',{});
        for i = 1:numel(R)
            e = ev{i}; e = e(e>=1 & e<=R(i).T);
            rec(i) = struct('folder',R(i).folder, 'roi',R(i).roi, ...
                'recName',R(i).recName, 'group',R(i).group, 'recDate',R(i).recDate, ...
                'fps',R(i).fps, 'T',R(i).T, 'dur_s',R(i).dur_s, ...
                'spike_idx',e(:), 'n_events',numel(e), ...
                'rate_per_min',numel(e)/max(R(i).dur_s,eps)*60, ...
                'toss',logical(toss(i)), 'note',char(note(i)), 'source',char(src(i)));
        end
        % A tossed CELL keeps nothing, whatever its per-recording flags say.
        kept = ~[rec.toss] & ~logical(cellToss);
        CUR = struct('stem',cellInfo.stem, 'cell',cellInfo.cell, 'key',cellInfo.key, ...
            'group',cellInfo.group, 'recDate',cellInfo.recDate, ...
            'rec',rec, 'reviewed',logical(reviewed), 'cell_toss',logical(cellToss), ...
            'n_rec',numel(rec), 'n_rec_kept',nnz(kept), ...
            'pooled_n_events',sum([rec(kept).n_events]), ...
            'pooled_dur_s',sum([rec(kept).dur_s]), ...
            'params',struct('fps',efFps.Value, 'threshold',efThr.Value, ...
                            'minDist',efDist.Value, 'minWidth',efWid.Value, ...
                            'prominence',efProm.Value, 'detrend',S.detrend, ...
                            'trimStart',efTrimStart.Value, 'trimEnd',efTrimEnd.Value), ...
            'bundle_file',CLk.file, 'saved',datestr(now)); %#ok<TNOW1,DATST>
        save(CLk.curFile, 'CUR', '-v7.3');
    end

    function writeCsvs(loud)
        rows = {};
        for k = 1:numel(S.CL)
            if ~isfile(S.CL(k).curFile), continue; end
            q = load(S.CL(k).curFile,'CUR');
            if ~isfield(q,'CUR'), continue; end
            C = q.CUR;
            for i = 1:numel(C.rec)
                ct = false;
                if isfield(C,'cell_toss'), ct = logical(C.cell_toss); end
                rows(end+1,:) = {C.cell, string(C.stem), string(C.group), ...
                    string(C.recDate), i, string(C.rec(i).recName), ...
                    string(C.rec(i).folder), C.rec(i).roi, C.rec(i).fps, ...
                    C.rec(i).dur_s, C.rec(i).n_events, C.rec(i).rate_per_min, ...
                    C.rec(i).toss, ct, ~C.rec(i).toss && ~ct, ...
                    string(C.rec(i).source), string(C.rec(i).note), ...
                    C.reviewed}; %#ok<AGROW>
            end
        end
        if isempty(rows)
            if loud, lblSave.Text = 'nothing curated yet'; end
            return
        end
        T = cell2table(rows, 'VariableNames', {'cell_idx','stem','group','recDate', ...
            'rec_in_cell','recName','folder','roi','fps','dur_s','n_events', ...
            'rate_per_min','toss','cell_toss','kept','source','note','cell_reviewed'});
        T = sortrows(T,{'cell_idx','rec_in_cell'});
        writetable(T, fullfile(S.bundleDir,'recording_curation.csv'));
        writetable(T(~T.kept,:), fullfile(S.bundleDir,'toss_recordings.csv'));
        if loud
            lblSave.Text = sprintf('wrote recording_curation.csv (%d rows, %d dropped: %d rec, %d whole cells)', ...
                height(T), nnz(~T.kept), nnz(T.toss & ~T.cell_toss), ...
                numel(unique(T.cell_idx(T.cell_toss))));
        end
    end
end


function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
