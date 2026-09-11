function stitch_manual_gui_260910()
%STITCH_MANUAL_GUI_260910  Refine an auto-stitched mosaic by hand.
%
% Simplified successor to stitch_manual_gui_fast_260818.m.  Same layout -- control
% column on the left, mosaic on the right -- and the same nudge/zoom/pan feel.
% What changed is the job it does.
%
% THE WORKFLOW IS NOW: run NChan_vol_stitch_zUse_260910 FIRST, then open this to
% fix whatever it got wrong.  That single assumption removes most of the old
% controls, because the things they existed to decide have already been decided:
%
%   REMOVED  Place dropdown (saved/stage/grid).  The starting placement is the
%            auto-stitch solution; stage placement is only the fallback when no
%            coords file exists yet.
%   REMOVED  Flip X / Flip Y.  In stage coordinates the sign convention is
%            measured, not chosen -- a flip here silently mirrors the whole
%            mosaic, which cost 1197 px of Y error once already.
%   REMOVED  stepX_um / stepY_um / umPerPxBase / zoom.  The step is solved from
%            the images by the stitcher and read back from its coords file.
%   REMOVED  "Start fresh", which duplicated Reset ALL.
%
%   ADDED    Channel selector -- the map GUI never had one, so refining on the
%            vessel channel and exporting the cell channel meant editing the file.
%   ADDED    Seam view: for the selected tile, each overlap as a RED/GREEN overlay
%            in X-Y and X-Z.  Yellow = registered.  Judging a seam in the blended
%            mosaic is judging it through the very blend designed to hide it.
%   ADDED    Black background, and Auto B/C that writes the LEVELS it chose back
%            into the manual boxes, so the automatic choice is visible and can be
%            taken over rather than being an opaque one-shot.
%
% Runqi Zhang / 2026-09-10

%% ========================= USER SETTINGS =========================
S.datasetRoot = "E:\260804_shiver_dbh\map";   % holds ch<N>/avgz and matlab_stitch
S.refChID     = 3;            % channel shown first; the popup lists what exists

S.zUse_um     = [-300 0];     % microns of motor z. [] = whole stack.
S.zUseProj    = [];           % page indices, used only when zUse_um is []
S.projMode    = 'avg';        % 'avg' | 'max'
S.workSpace   = "raw";        % "raw" | "rescaled"
                              % which of the stitcher's two coordinate spaces to
                              % nudge in.  "raw" = the acquired pixels, nothing
                              % interpolated.  "rescaled" = 1 px is 1 um.

S.clipLo      = 1;            % Auto B/C percentiles -> absolute levels in the boxes
S.clipHi      = 99.5;
S.gamma       = 0.8;          % display-only, <1 brighter

S.moveStep    = 5;            % px per arrow press (Shift = 10x)
S.feather_um  = 30;           % seam-view blend ramp
% NO projection is written -- only the avgz VOLUME.  Projections are made in
% ImageJ, where the TIFF scale and units have to be set by hand anyway, so one
% written here would just be a second copy to keep in sync.

S.dsFast      = 4;            % preview downsample while nudging
S.dsIdle      = 1;            % 1 = full resolution once idle
S.idleSec     = 0.30;
%% =================================================================

repoRoot = fileparts(fileparts(mfilename('fullpath'))); addpath(repoRoot);
here = fileparts(mfilename('fullpath')); addpath(here);

BG   = [0.13 0.13 0.13];
FG   = [0.92 0.92 0.92];
DIM  = [0.55 0.55 0.55];
COL  = [ 42 120 214; 235 104  52] / 255;   % seam-view channel hues (CVD-checked)

%% ---------------------------- state ----------------------------
S.chList = []; S.files = strings(0,1); S.nTiles = 0;
S.H0 = 0; S.W0 = 0; S.Hn = 0; S.Wn = 0; S.inClass = 'uint16';
S.colIdx = []; S.rowIdx = []; S.nRows = 0; S.nCols = 0;
S.tAvg = {}; S.tMax = {}; S.tXZ = {};        % per-tile projections (full res)
S.tD = {}; S.dsCache = NaN; S.wD = [];       % downsampled + its blend weights
S.x = []; S.y = []; S.manX = []; S.manY = [];
S.baseX = []; S.baseY = [];                  % placement before any nudge
S.zPages = []; S.zHow = '';
S.sel = 0; S.moveMode = 'Tile';
S.clim = [0 1]; S.loaded = false; S.pan = [];
S.outW = 0; S.outH = 0; S.dsPreview = S.dsFast;
S.idleTimer = []; S.src = '';
S.rescaled_um = [];
S.sxPx = NaN; S.syPx = NaN; S.rescaleTarget = [];
S.rawH = 0; S.rawW = 0; S.umH = 0; S.umW = 0;

seamFig = []; seamAx = []; seamIm = []; seamKey = '';

%% ========================= FIGURE =========================
fig = figure('Name','Manual Stitch 260910','NumberTitle','off', ...
    'Units','pixels','Position',[60 60 1500 950], 'Color',BG, ...
    'MenuBar','none','ToolBar','none','Renderer','opengl', ...
    'GraphicsSmoothing','off','DoubleBuffer','on', ...
    'KeyPressFcn',@cb_key,'WindowScrollWheelFcn',@cb_scroll, ...
    'WindowButtonDownFcn',@cb_btndown,'WindowButtonMotionFcn',@cb_btnmotion, ...
    'WindowButtonUpFcn',@cb_btnup,'CloseRequestFcn',@(s,~) closeGui(s));
try, fig.WindowState = 'maximized'; catch, end

ctrl = uipanel(fig,'Title','Controls','FontWeight','bold', ...
    'BackgroundColor',BG,'ForegroundColor',FG, ...
    'Units','normalized','Position',[0.005 0.005 0.185 0.99]);

ax = axes('Parent',fig,'Units','normalized','Position',[0.205 0.04 0.79 0.92], ...
          'Color','k','XColor',DIM,'YColor',DIM);
colormap(ax, gray(256));
ax.XTick=[]; ax.YTick=[]; ax.YDir='reverse'; ax.DataAspectRatio=[1 1 1]; ax.Box='on';
hTitle = title(ax,'Manual stitch','Interpreter','none','FontSize',10,'Color',FG);
hImg = []; hBord = []; hSel = [];

yy = 0.965; H = 0.026; G = 0.009;
    function p = row(n), p = [0.05 yy-n*H 0.90 n*H]; end
    function adv(n), yy = yy - n*H - G; end
    function lab(t)
        uicontrol(ctrl,'Style','text','Units','normalized','Position',row(1), ...
            'String',t,'HorizontalAlignment','left','FontWeight','bold', ...
            'ForegroundColor',[.55 .75 .95],'BackgroundColor',BG); adv(1);
    end
    function h = txt(t, p, fs, col)
        h = uicontrol(ctrl,'Style','text','Units','normalized','Position',p, ...
            'String',t,'HorizontalAlignment','left','FontSize',fs, ...
            'ForegroundColor',col,'BackgroundColor',BG);
    end
    function h = btn(t, cb, p)
        h = uicontrol(ctrl,'Style','pushbutton','Units','normalized','Position',p, ...
            'String',t,'Callback',cb,'KeyPressFcn',@cb_key);
    end

lab('-- Load --');
btn('Load dataset', @(~,~) cb_load(), [0.05 yy-H 0.90 H]); adv(1);
lblPath = txt(char(S.datasetRoot), row(2), 7, DIM); adv(2);

lab('-- Channel / Z --');
txt('Channel', [0.05 yy-H 0.36 H], 8, FG);
ddChan = uicontrol(ctrl,'Style','popupmenu','Units','normalized', ...
    'Position',[0.42 yy-H 0.53 H],'String',{'(load)'},'Callback',@cb_chan); adv(1);
txt('Z proj', [0.05 yy-H 0.36 H], 8, FG);
ddProj = uicontrol(ctrl,'Style','popupmenu','Units','normalized', ...
    'Position',[0.42 yy-H 0.53 H],'String',{'avg','max'}, ...
    'Value',max(1,find(strcmp({'avg','max'},S.projMode),1)),'Callback',@cb_proj); adv(1);
lblSrc = txt('', row(3), 7, [.6 .85 .6]); adv(3);

lab('-- Manual nudge --');
lblSel = txt('Selected: none (click a tile)', row(2), 8, [.55 .75 .95]); adv(2);
txt('Move:', [0.05 yy-H 0.42 H], 8, FG);
ddMove = uicontrol(ctrl,'Style','popupmenu','Units','normalized', ...
    'Position',[0.48 yy-H 0.47 H],'String',{'Tile','Row','Col'},'Value',1, ...
    'Callback',@(s,~) setMove(s)); adv(1);
txt('Move step (px):', [0.05 yy-H 0.60 H], 8, FG);
efMove = uicontrol(ctrl,'Style','edit','Units','normalized','Position',[0.66 yy-H 0.29 H], ...
    'String',num2str(S.moveStep),'Callback',@(s,~) setMoveStep(s)); adv(1);
txt('Arrow = step   Shift+Arrow = 10x', row(2), 7, DIM); adv(2);
btn('Reset this tile', @(~,~) cb_reset(false), [0.05 yy-H 0.44 H]);
btn('Reset ALL',       @(~,~) cb_reset(true),  [0.51 yy-H 0.44 H]); adv(1);

lab('-- Display --');
cbSeams = uicontrol(ctrl,'Style','checkbox','Units','normalized', ...
    'Position',[0.05 yy-H 0.50 H],'String','Borders','Value',1, ...
    'ForegroundColor',FG,'BackgroundColor',BG, ...
    'Callback',@(~,~) setBorders(),'KeyPressFcn',@cb_key);
btn('Auto B/C', @(~,~) autoBC(true), [0.57 yy-H 0.38 H]); adv(1);
% Black / White hold the ABSOLUTE levels, not percentiles.  Auto B/C fills them
% in, so the automatic choice is visible and can be taken over by typing --
% rather than being a one-shot you cannot see or adjust from.
txt('Black', [0.05 yy-H 0.26 H], 8, FG);
efLo = uicontrol(ctrl,'Style','edit','Units','normalized','Position',[0.30 yy-H 0.20 H], ...
    'String','0','Callback',@(s,~) setClim());
txt('White', [0.53 yy-H 0.26 H], 8, FG);
efHi = uicontrol(ctrl,'Style','edit','Units','normalized','Position',[0.76 yy-H 0.19 H], ...
    'String','1','Callback',@(s,~) setClim()); adv(1);
txt('Gamma', [0.05 yy-H 0.26 H], 8, FG);
efGam = uicontrol(ctrl,'Style','edit','Units','normalized','Position',[0.30 yy-H 0.20 H], ...
    'String',num2str(S.gamma),'Callback',@(~,~) setGamma());
txt('(<1 brighter)', [0.53 yy-H 0.42 H], 7, DIM); adv(1);
btn('Reset view', @(~,~) resetView(),      [0.05 yy-H 0.44 H]);
btn('Seam view',  @(~,~) showSeams(true),  [0.51 yy-H 0.44 H]); adv(1);

lab('-- Save --');
btn('Save volumes + coords', @(~,~) cb_save(), [0.05 yy-H 0.90 H]); adv(1);
lblStat = txt('Ready. Click "Load dataset".', row(4), 8, FG); adv(4);
txt(sprintf(['scroll = zoom\nright drag = pan\nr = reset view   b = borders\n' ...
             'v = seam view   1..9 = pick tile']), row(5), 7, DIM);

hFocus = uicontrol(ctrl,'Style','pushbutton','Units','pixels', ...
    'Position',[-20 -20 1 1],'KeyPressFcn',@cb_key);

cb_load();

%% ============================== LOAD ==============================
    function cb_load()
        root = S.datasetRoot;
        d = dir(fullfile(root, "ch*"));
        d = d([d.isdir]);
        chs = [];
        for k = 1:numel(d)
            t = regexp(d(k).name, '^ch(\d+)$', 'tokens', 'once');
            if ~isempty(t) && isfolder(fullfile(root, d(k).name, "avgz"))
                chs(end+1) = str2double(t{1}); %#ok<AGROW>
            end
        end
        assert(~isempty(chs), 'no ch<N>/avgz folders under %s', root);
        S.chList = sort(chs);
        if ~ismember(S.refChID, S.chList), S.refChID = S.chList(1); end
        set(ddChan, 'String', arrayfun(@(c) sprintf('ch %d', c), S.chList, 'uni', 0), ...
                    'Value', find(S.chList == S.refChID, 1));
        set(lblPath, 'String', char(root));

        loadChannel();
        S.loaded = true;
        autoBC(true);
        resetView();
        redraw();
        showSeams(false);
    end

    function loadChannel()
        setStat('listing tiles ...'); drawnow;
        avgzDir = fullfile(S.datasetRoot, sprintf("ch%d", S.refChID), "avgz");
        c = [dir(fullfile(avgzDir,"*.tif")); dir(fullfile(avgzDir,"*.tiff"))];
        assert(~isempty(c), 'no TIFFs in %s', avgzDir);
        S.files = string(fullfile({c.folder}, {c.name})).';
        S.nTiles = numel(S.files);

        [S.colIdx, S.rowIdx] = parse_grid(S.files);
        uC = unique(S.colIdx(~isnan(S.colIdx)));
        uR = unique(S.rowIdx(~isnan(S.rowIdx)));
        S.nCols = numel(uC); S.nRows = numel(uR);

        [S.H0, S.W0, S.inClass] = first_page(S.files(1));
        nZ = count_pages(S.files(1));
        [S.zPages, S.zHow] = resolve_z(S.datasetRoot, S.files(1), nZ, S.zUse_um, S.zUseProj);

        % Projections: avg, max AND an X-Z view, all in one pass over the pages.
        % The X-Z costs one max() per page and ~10 MB for the whole map, and it is
        % the only way the seam view can show a depth offset.
        S.tAvg = cell(S.nTiles,1); S.tMax = cell(S.nTiles,1); S.tXZ = cell(S.nTiles,1);
        for t = 1:S.nTiles
            [S.tAvg{t}, S.tMax{t}, S.tXZ{t}] = read_proj(S.files(t), S.zPages, S.H0, S.W0);
            if mod(t,10)==0 || t==S.nTiles
                setStat(sprintf('reading %d/%d tiles ...', t, S.nTiles)); drawnow;
            end
        end
        S.dsCache = NaN;

        loadPlacement();
    end

    function loadPlacement()
        % Placement comes from the auto-stitcher's coords file when there is one.
        % Refining is the whole point of this GUI, so starting anywhere else would
        % mean throwing away the solve you just ran.
        S.rescaled_um = [];
        cf = "";
        for k = [S.refChID, S.chList]
            f = fullfile(S.datasetRoot, "matlab_stitch", sprintf("stitch_coords_from_ch%d.mat", k));
            if isfile(f), cf = f; break; end
        end
        if cf ~= ""
            C = load(cf);
            names = base_noch(S.files);
            cnames = base_noch(C.tileFileRef(:));
            [tf, loc] = ismember(names, cnames);
            if all(tf) && numel(C.x) == S.nTiles
                % The stitcher saves BOTH spaces: raw pixels (what it solved in)
                % and rescaled microns (what it wrote).  workSpace picks which one
                % this GUI nudges in, so a nudge of 1 px means the same thing here
                % as in the file that gets saved back.
                useUm = S.workSpace == "rescaled" && isfield(C,'x_um') && ...
                        isfield(C,'rescaled_um') && ~isempty(C.rescaled_um);
                if useUm
                    S.baseX = C.x_um(loc);  S.baseY = C.y_um(loc);
                    S.Hn = C.Hn;  S.Wn = C.Wn;  S.rescaled_um = C.rescaled_um;
                    spaceTxt = sprintf('rescaled %g um/px', C.rescaled_um);
                elseif isfield(C,'x_raw')
                    S.baseX = C.x_raw(loc); S.baseY = C.y_raw(loc);
                    S.Hn = C.H0;  S.Wn = C.W0;  S.rescaled_um = [];
                    spaceTxt = 'raw pixels';
                else                      % coords from an older run
                    S.baseX = C.x(loc); S.baseY = C.y(loc);
                    if isfield(C,'Hn'), S.Hn = C.Hn; S.Wn = C.Wn;
                    else,               S.Hn = S.H0; S.Wn = S.W0; end
                    if isfield(C,'rescaled_um'), S.rescaled_um = C.rescaled_um; end
                    spaceTxt = 'legacy coords';
                end
                % Keep BOTH conversion factors, whichever space is being
                % worked in, so Save can emit the other one without re-solving.
                if isfield(C,'sxPx'), S.sxPx = C.sxPx; S.syPx = C.syPx;
                else,                 S.sxPx = NaN;    S.syPx = NaN; end
                if isfield(C,'rescaled_um'), S.rescaleTarget = C.rescaled_um;
                else,                        S.rescaleTarget = []; end
                if isfield(C,'H0'), S.rawH = C.H0; S.rawW = C.W0;
                else,               S.rawH = S.H0; S.rawW = S.W0; end
                if isfield(C,'Hn'), S.umH = C.Hn; S.umW = C.Wn;
                else,               S.umH = S.rawH; S.umW = S.rawW; end
                [~, nm, ex] = fileparts(cf);
                S.src = sprintf('auto-stitch %s%s  [%s]', nm, ex, spaceTxt);
            else
                cf = "";
            end
        end
        if cf == ""
            % Fallback: stage positions, at the header pixel size.
            [sx, sy, ok, how] = stage_positions(S.files, S.datasetRoot);
            acqPx = header_um_per_px(S.datasetRoot, S.files(1));
            assert(ok && isfinite(acqPx), ['no coords file and no usable stage ' ...
                'positions -- run NChan_vol_stitch_zUse_260910 first.']);
            S.baseX = (sx - min(sx)) / acqPx;
            S.baseY = (max(sy) - sy) / acqPx;
            S.Hn = S.H0; S.Wn = S.W0;
            S.src = sprintf('stage (%s), %.4f um/px -- NO auto-stitch coords found', how, acqPx);
        end

        % Resample the tiles into the coords' pixel space, so a nudge of 1 px in
        % this GUI is 1 px of the mosaic the stitcher built.
        if S.Hn ~= S.H0 || S.Wn ~= S.W0
            for t = 1:S.nTiles
                S.tAvg{t} = imresize(S.tAvg{t}, [S.Hn S.Wn], 'bicubic');
                S.tMax{t} = imresize(S.tMax{t}, [S.Hn S.Wn], 'bicubic');
                S.tXZ{t}  = imresize(S.tXZ{t}, [size(S.tXZ{t},1) S.Wn], 'bicubic');
            end
            S.dsCache = NaN;
        end

        S.manX = zeros(S.nTiles,1); S.manY = zeros(S.nTiles,1);
        S.sel = 0;
        set(lblSrc, 'String', S.src);
        setStat(sprintf('%d tiles, %dx%d grid, %s', S.nTiles, S.nRows, S.nCols, S.zHow));
    end

%% ============================ COMPOSE =============================
    function P = proj(t)
        if strcmp(S.projMode,'max'), P = S.tMax{t}; else, P = S.tAvg{t}; end
    end

    function buildDS(ds)
        if isequaln(S.dsCache, ds), return; end
        S.tD = cell(S.nTiles,1);
        for t = 1:S.nTiles
            if ds == 1, S.tD{t} = proj(t); else, S.tD{t} = imresize(proj(t), 1/ds, 'bilinear'); end
        end
        [h, w] = size(S.tD{1});
        S.wD = feather(h, w, max(1, (S.feather_um / max(1,scale_um())) / ds));
        S.dsCache = ds;
    end

    function u = scale_um()
        if ~isempty(S.rescaled_um) && S.rescaled_um > 0, u = S.rescaled_um; else, u = 1; end
    end

    function M = mosaic(ds)
        buildDS(ds);
        X = round((S.baseX + S.manX) / ds);  Y = round((S.baseY + S.manY) / ds);
        X = X - min(X);  Y = Y - min(Y);
        [h, w] = size(S.tD{1});
        Wd = max(X) + w;  Hd = max(Y) + h;
        num = zeros(Hd, Wd, 'single'); den = zeros(Hd, Wd, 'single');
        for t = 1:S.nTiles
            rr = Y(t)+(1:h); cc = X(t)+(1:w);
            num(rr,cc) = num(rr,cc) + single(S.tD{t}) .* S.wD;
            den(rr,cc) = den(rr,cc) + S.wD;
        end
        den(den==0) = 1;  M = num ./ den;
        S.outW = Wd*ds; S.outH = Hd*ds;
    end

    function redraw()
        if ~S.loaded, return; end
        M = mosaic(S.dsPreview);
        if isempty(hImg) || ~isvalid(hImg)
            hImg = imagesc(ax, M); axis(ax,'image');
            ax.XTick=[]; ax.YTick=[]; ax.YDir='reverse';
            set(ax,'Color','k','XColor',DIM,'YColor',DIM);
        else
            set(hImg,'CData',M,'XData',[1 size(M,2)],'YData',[1 size(M,1)]);
        end
        applyGamma();
        caxis(ax, S.clim);
        drawBorders();
        updateSeams();
        % The title is re-created, not cached.  imagesc() calls newplot(), which
        % CLEARS the axes and deletes its title object -- so a handle grabbed
        % before the first draw is already dead by the time it is written to.
        hTitle = title(ax, sprintf('%s  |  ch%d  |  %s  |  ds %d', ...
            char(S.datasetRoot), S.refChID, S.projMode, S.dsPreview), ...
            'Interpreter','none','FontSize',10,'Color',FG);
    end

    function drawBorders()
        delete(hBord(isgraphics(hBord))); hBord = gobjects(0);
        delete(hSel(isgraphics(hSel)));   hSel  = gobjects(0);
        if ~get(cbSeams,'Value'), return; end
        ds = S.dsPreview;
        X = round((S.baseX + S.manX)/ds); Y = round((S.baseY + S.manY)/ds);
        X = X - min(X); Y = Y - min(Y);
        [h, w] = size(S.tD{1});
        hold(ax,'on');
        hBord = gobjects(S.nTiles,1);
        for t = 1:S.nTiles
            c = [0.35 0.75 1]; lw = 0.5;
            if t == S.sel, c = [1 0.85 0.1]; lw = 1.75; end
            hBord(t) = rectangle(ax,'Position',[X(t)+0.5 Y(t)+0.5 w h],'EdgeColor',c,'LineWidth',lw);
        end
        hold(ax,'off');
    end

%% =========================== DISPLAY ==============================
    function autoBC(write)
        v = [];
        for t = 1:S.nTiles
            P = proj(t);
            v = [v; single(reshape(P(1:6:end,1:6:end), [], 1))]; %#ok<AGROW>
        end
        q = double(prctile(v, [S.clipLo S.clipHi]));
        if q(2) <= q(1), q(2) = q(1) + 1; end
        S.clim = q;
        if write
            set(efLo,'String',num2str(q(1),'%.6g'));
            set(efHi,'String',num2str(q(2),'%.6g'));
        end
        if ~isempty(hImg) && isvalid(hImg), caxis(ax, S.clim); end
    end

    function setClim()
        lo = str2double(get(efLo,'String'));  hi = str2double(get(efHi,'String'));
        if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
            setStat('black/white must be finite with white > black');
            set(efLo,'String',num2str(S.clim(1),'%.6g'));
            set(efHi,'String',num2str(S.clim(2),'%.6g'));
        else
            S.clim = [lo hi];  caxis(ax, S.clim);
        end
        uicontrol(hFocus);
    end

    function applyGamma()
        g = S.gamma; if ~isfinite(g) || g <= 0, g = 1; end
        cm = gray(256); colormap(ax, cm .^ g);
    end

    function setGamma()
        g = str2double(get(efGam,'String'));
        if isfinite(g) && g > 0, S.gamma = g; else, set(efGam,'String',num2str(S.gamma)); end
        applyGamma(); uicontrol(hFocus);
    end

    function setBorders(), drawBorders(); end

    function resetView()
        if isempty(hImg) || ~isvalid(hImg), return; end
        xlim(ax, [0.5 size(get(hImg,'CData'),2)+0.5]);
        ylim(ax, [0.5 size(get(hImg,'CData'),1)+0.5]);
    end

%% ========================= SEAM VIEW ==============================
    function P = seamPairs(k)
        X = round(S.baseX + S.manX); Y = round(S.baseY + S.manY);
        P = struct('j',{},'rk',{},'ck',{},'rj',{},'cj',{});
        for j = 1:S.nTiles
            if j == k, continue; end
            c0 = max(X(k),X(j)); c1 = min(X(k)+S.Wn, X(j)+S.Wn) - 1;
            r0 = max(Y(k),Y(j)); r1 = min(Y(k)+S.Hn, Y(j)+S.Hn) - 1;
            if c1 - c0 < 8 || r1 - r0 < 8, continue; end
            P(end+1) = struct('j',j, 'rk',(r0:r1)-Y(k)+1, 'ck',(c0:c1)-X(k)+1, ...
                                     'rj',(r0:r1)-Y(j)+1, 'cj',(c0:c1)-X(j)+1); %#ok<AGROW>
        end
    end

    function showSeams(force)
        if isempty(seamFig) || ~isvalid(seamFig)
            if ~force, return; end
            seamFig = figure('Name','seam view','NumberTitle','off','Color',BG, ...
                'Units','pixels','Position',[80 80 1150 650],'MenuBar','none', ...
                'ToolBar','none','KeyPressFcn',@cb_key);
            seamKey = '';
        end
        updateSeams();
    end

    function updateSeams()
        if isempty(seamFig) || ~isvalid(seamFig) || ~S.loaded, return; end
        set(seamFig,'Name',sprintf('seam view  |  ch %d  (%s)', S.refChID, S.projMode));
        if S.sel < 1
            clf(seamFig);
            annotation(seamFig,'textbox',[0.1 0.45 0.8 0.1],'String', ...
                'click a tile to see its seams','Color',FG,'EdgeColor','none', ...
                'HorizontalAlignment','center');
            seamKey = ''; return
        end
        P = seamPairs(S.sel);
        key = sprintf('%d:%s:%d', S.sel, mat2str([P.j]), S.refChID);
        nS = numel(P);
        if nS == 0
            clf(seamFig);
            annotation(seamFig,'textbox',[0.1 0.45 0.8 0.1],'String', ...
                sprintf('tile %d overlaps nothing', S.sel),'Color',FG, ...
                'EdgeColor','none','HorizontalAlignment','center');
            seamKey = key; return
        end
        if ~strcmp(key, seamKey)
            clf(seamFig);
            tl = tiledlayout(seamFig, nS, 4, 'TileSpacing','compact','Padding','compact');
            seamAx = gobjects(nS,4); seamIm = gobjects(nS,4);
            for s = 1:nS
                for q = 1:4
                    seamAx(s,q) = nexttile(tl, (s-1)*4+q);
                    seamIm(s,q) = imagesc(seamAx(s,q), 0);
                    axis(seamAx(s,q),'image');
                    set(seamAx(s,q),'XTick',[],'YTick',[],'Color','k', ...
                        'XColor',DIM,'YColor',DIM,'Box','on');
                    if q <= 2, colormap(seamAx(s,q), gray); end
                end
            end
            seamKey = key;
        end
        for s = 1:nS
            p = P(s); k = S.sel; j = p.j;
            A = proj(k); B = proj(j);
            a = A(p.rk, p.ck);  b = B(p.rj, p.cj);
            az = S.tXZ{k}(:, p.ck);  bz = S.tXZ{j}(:, p.cj);
            set(seamIm(s,1),'CData',a); caxis(seamAx(s,1), S.clim);
            set(seamIm(s,2),'CData',b); caxis(seamAx(s,2), S.clim);
            set(seamIm(s,3),'CData',rg(a,b));
            set(seamIm(s,4),'CData',rg(az,bz));
            for q = 1:4
                sz = size(get(seamIm(s,q),'CData'));
                set(seamIm(s,q),'XData',[1 sz(2)],'YData',[1 sz(1)]);
                xlim(seamAx(s,q),[0.5 sz(2)+0.5]); ylim(seamAx(s,q),[0.5 sz(1)+0.5]);
            end
            title(seamAx(s,1), sprintf('tile %d (red)', k), 'Color',COL(1,:)*0+[1 .5 .5], 'FontWeight','normal');
            title(seamAx(s,2), sprintf('tile %d (green)', j), 'Color',[.5 1 .5], 'FontWeight','normal');
            title(seamAx(s,3), sprintf('X-Y overlap  %dx%d px', numel(p.ck), numel(p.rk)), ...
                  'Color',FG,'FontWeight','normal');
            title(seamAx(s,4), 'X-Z overlap', 'Color',FG,'FontWeight','normal');
        end
    end

%% ======================== INTERACTION =============================
    function cb_key(~, e)
        if ~S.loaded, return; end
        sh = any(strcmp(e.Modifier,'shift'));
        st = S.moveStep * (1 + 9*sh);
        switch e.Key
            case 'leftarrow',  nudge(-st, 0);
            case 'rightarrow', nudge( st, 0);
            case 'uparrow',    nudge(0, -st);
            case 'downarrow',  nudge(0,  st);
            case 'r', resetView();
            case 'b', set(cbSeams,'Value',~get(cbSeams,'Value')); drawBorders();
            case 'v', showSeams(true);
            otherwise
                k = str2double(e.Key);
                if ~isnan(k) && k >= 1 && k <= S.nTiles, S.sel = k; fastRedraw(); end
        end
    end

    function nudge(dx, dy)
        if S.sel < 1, setStat('click a tile first'); return; end
        switch S.moveMode
            case 'Tile', idx = S.sel;
            case 'Row',  idx = find(S.rowIdx == S.rowIdx(S.sel));
            case 'Col',  idx = find(S.colIdx == S.colIdx(S.sel));
        end
        S.manX(idx) = S.manX(idx) + dx;
        S.manY(idx) = S.manY(idx) + dy;
        fastRedraw();
    end

    function fastRedraw()
        % Coarse while keys are moving, sharp once they stop.  Blending the full
        % mosaic costs ~300 ms, too slow to nudge in; /4 is fast but soft, so the
        % sharp version is what you judge on and the coarse one is what you drive.
        S.dsPreview = S.dsFast;
        redraw(); setSelLabel();
        if ~isempty(S.idleTimer) && isvalid(S.idleTimer), stop(S.idleTimer); delete(S.idleTimer); end
        S.idleTimer = timer('StartDelay',S.idleSec,'TimerFcn',@(~,~) sharpen());
        start(S.idleTimer);
    end

    function sharpen()
        S.dsPreview = S.dsIdle;
        try, redraw(); catch, end
    end

    function cb_btndown(~,~)
        if ~S.loaded, return; end
        switch get(fig,'SelectionType')
            case 'normal'
                cp = get(ax,'CurrentPoint'); px = cp(1,1); py = cp(1,2);
                ds = S.dsPreview;
                X = round((S.baseX+S.manX)/ds); Y = round((S.baseY+S.manY)/ds);
                X = X - min(X); Y = Y - min(Y);
                [h, w] = size(S.tD{1});
                for t = S.nTiles:-1:1
                    if px>=X(t) && px<=X(t)+w && py>=Y(t) && py<=Y(t)+h
                        S.sel = t; fastRedraw(); return
                    end
                end
            case 'open', resetView();
            otherwise
                S.pan = struct('p0',get(ax,'CurrentPoint'),'xl',xlim(ax),'yl',ylim(ax));
        end
    end

    function cb_btnmotion(~,~)
        if isempty(S.pan), return; end
        cp = get(ax,'CurrentPoint');
        xlim(ax, S.pan.xl + (S.pan.p0(1,1)-cp(1,1)));
        ylim(ax, S.pan.yl + (S.pan.p0(1,2)-cp(1,2)));
    end

    function cb_btnup(~,~), S.pan = []; end

    function cb_scroll(~, e)
        if ~S.loaded, return; end
        cp = get(ax,'CurrentPoint'); cx = cp(1,1); cy = cp(1,2);
        f = 1.2 ^ double(e.VerticalScrollCount);
        xlim(ax, cx + (xlim(ax)-cx)*f);  ylim(ax, cy + (ylim(ax)-cy)*f);
    end

    function cb_chan(~,~)
        S.refChID = S.chList(get(ddChan,'Value'));
        setStat(sprintf('loading ch%d ...', S.refChID)); drawnow;
        mx = S.manX; my = S.manY;        % nudges are geometry, not pixels -- keep them
        loadChannel();
        if numel(mx) == S.nTiles, S.manX = mx; S.manY = my; end
        autoBC(true); redraw(); uicontrol(hFocus);
    end

    function cb_proj(~,~)
        v = get(ddProj,'String'); S.projMode = v{get(ddProj,'Value')};
        S.dsCache = NaN; autoBC(true); redraw(); uicontrol(hFocus);
    end

    function cb_reset(all)
        if all, S.manX(:) = 0; S.manY(:) = 0;
        elseif S.sel >= 1, S.manX(S.sel) = 0; S.manY(S.sel) = 0; end
        fastRedraw();
    end

    function setMove(s)
        v = get(s,'String'); S.moveMode = v{get(s,'Value')}; uicontrol(hFocus);
    end

    function setMoveStep(s)
        v = str2double(get(s,'String'));
        if isfinite(v) && v >= 1, S.moveStep = round(v); end
        set(s,'String',num2str(S.moveStep)); uicontrol(hFocus);
    end

    function setSelLabel()
        if S.sel < 1, set(lblSel,'String','Selected: none (click a tile)'); return; end
        [~, nm] = fileparts(S.files(S.sel));
        set(lblSel,'String',sprintf('#%d col%d row%d   dx%+g dy%+g\n%s', ...
            S.sel, S.colIdx(S.sel), S.rowIdx(S.sel), S.manX(S.sel), S.manY(S.sel), nm));
    end

    function setStat(t), set(lblStat,'String',t); end

%% ============================= SAVE ===============================
    function cb_save()
        outDir = fullfile(S.datasetRoot, "matlab_stitch");
        if ~isfolder(outDir), mkdir(outDir); end

        % Curated coords in the space they were nudged in, plus the factors
        % needed to move to the other one -- so nothing here has to be re-derived.
        xw = S.baseX + S.manX;  yw = S.baseY + S.manY;
        xw = xw - min(xw);      yw = yw - min(yw);
        workSpace = S.workSpace;  sxPx = S.sxPx;  syPx = S.syPx; %#ok<NASGU>
        rescaleTarget = S.rescaleTarget; %#ok<NASGU>
        x = xw; y = yw; %#ok<NASGU>
        files = S.files; H0 = S.rawH; W0 = S.rawW; Hn = S.umH; Wn = S.umW; %#ok<NASGU>
        manX = S.manX; manY = S.manY; baseX = S.baseX; baseY = S.baseY; %#ok<NASGU>
        colIdx = S.colIdx; rowIdx = S.rowIdx; %#ok<NASGU>
        refChID = S.refChID; zPages = S.zPages; zUse_um = S.zUse_um; %#ok<NASGU>
        srcPlacement = S.src; %#ok<NASGU>
        savedOn = datestr(now,'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST,NASGU>

        cf = fullfile(outDir, "stitch_manual_coords.mat");
        if isfile(cf)
            copyfile(cf, fullfile(outDir, sprintf("stitch_manual_coords_backup_%s.mat", ...
                     datestr(now,'yymmdd_HHMMSS')))); %#ok<TNOW1,DATST>
        end
        save(cf, 'x','y','xw','yw','workSpace','sxPx','syPx','rescaleTarget', ...
                 'files','H0','W0','Hn','Wn','manX','manY','baseX','baseY', ...
                 'colIdx','rowIdx','refChID','zPages','zUse_um','srcPlacement','savedOn');
        setStat(sprintf('coords -> %s', cf)); drawnow;

        % BOTH SPACES, EVERY ACQUIRED CHANNEL, from the one curated placement.
        % The nudges are geometry, not pixels: converting them to the other space
        % is a multiplication, so there is no reason to make you re-curate or
        % re-run to get the other rendering.
        % THIS IS AN UPDATE STEP.  It overwrites the stitcher's own outputs under
        % the SAME names, so a curated placement never sits beside a stale mosaic
        % for a downstream script to pick the wrong one of.
        %
        % It therefore writes the SAME ARTIFACT the stitcher wrote -- the z-page
        % VOLUME, not a projection.  Overwriting a 31-page volume with a 1-page
        % projection under the same filename would look like a clean save and
        % quietly destroy 30 planes.
        SP = spaceList();
        for sp = 1:numel(SP)
            for ch = S.chList
                fn = fullfile(outDir, sprintf("stitched_ch%d_avgz_%s.tif", ch, SP(sp).tag));
                setStat(sprintf('re-rendering ch%d %s (%d pages) ...', ...
                        ch, SP(sp).tag, numel(S.zPages))); drawnow;
                write_volume(fn, ch, SP(sp));
                fprintf('updated %s\n', fn);
            end
        end
        setStat(sprintf('updated %d volumes (%d ch x %d space) + coords', ...
                numel(SP)*numel(S.chList), numel(S.chList), numel(SP)));
    end

    function SP = spaceList()
        % The curated coords live in workSpace; the other space is reached by the
        % ratio of the measured raw pixel size to the rescale target.  Per axis,
        % because x and y differ by ~6% and a single factor would shear the mosaic.
        xw = S.baseX + S.manX;  yw = S.baseY + S.manY;
        xw = xw - min(xw);      yw = yw - min(yw);
        haveUm = ~isempty(S.rescaleTarget) && S.rescaleTarget > 0 && ...
                 isfinite(S.sxPx) && isfinite(S.syPx) && S.umH > 0;
        % Tags carry the pixel size, matching the stitcher exactly so the GUI
        % overwrites the files it is meant to update rather than creating a
        % parallel set under near-identical names.
        rawTag = sprintf("raw_%gum", round(mean([S.sxPx S.syPx]), 1));
        if ~isfinite(S.sxPx), rawTag = "raw"; end
        if S.workSpace == "rescaled"
            umTag = sprintf("rescaled_%gum", S.rescaleTarget);
            fx = S.rescaleTarget / S.sxPx;   fy = S.rescaleTarget / S.syPx;
            SP = struct('tag',umTag,'x',xw,'y',yw,'H',S.umH,'W',S.umW);
            SP(2) = struct('tag',rawTag,'x',xw*fx,'y',yw*fy,'H',S.rawH,'W',S.rawW);
        else
            SP = struct('tag',rawTag,'x',xw,'y',yw,'H',S.rawH,'W',S.rawW);
            if haveUm
                umTag = sprintf("rescaled_%gum", S.rescaleTarget);
                fx = S.sxPx / S.rescaleTarget;  fy = S.syPx / S.rescaleTarget;
                SP(2) = struct('tag',umTag,'x',xw*fx,'y',yw*fy,'H',S.umH,'W',S.umW);
            end
        end
    end

    function write_volume(fn, ch, SP)
        % Streams the z pages the GUI is showing, blends each into the mosaic and
        % writes it out -- the same thing the stitcher does, so the file it
        % replaces is of the same kind.  Written page by page: the full volume in
        % memory would be several GB on a 9x9 map.
        X = round(SP.x); Y = round(SP.y); X = X - min(X); Y = Y - min(Y);
        Wo = max(X) + SP.W;  Ho = max(Y) + SP.H;
        w2 = feather(SP.H, SP.W, max(1, S.feather_um / max(1, scale_um())));
        rd = cell(S.nTiles,1);
        for t = 1:S.nTiles
            rd{t} = Tiff(swap_ch(S.files(t), S.refChID, ch), 'r');
        end
        cl = onCleanup(@() cellfun(@(o) o.close(), rd));
        if isfile(fn), delete(fn); end
        tw = Tiff(fn, 'w8');  cw = onCleanup(@() tw.close());
        tag.ImageLength = Ho; tag.ImageWidth = Wo;
        tag.Photometric = Tiff.Photometric.MinIsBlack;
        tag.BitsPerSample = 16; tag.SamplesPerPixel = 1;
        tag.SampleFormat = Tiff.SampleFormat.UInt;
        tag.RowsPerStrip = Ho;        % one strip: MATLAB pads a partial strip
        tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tag.Compression = Tiff.Compression.None;
        tag.Software = 'stitch_manual_gui_260910';
        for zi = 1:numel(S.zPages)
            num = zeros(Ho, Wo, 'single'); den = zeros(Ho, Wo, 'single');
            for t = 1:S.nTiles
                rd{t}.setDirectory(S.zPages(zi));
                P = single(rd{t}.read());
                if SP.H ~= S.H0 || SP.W ~= S.W0
                    P = imresize(P, [SP.H SP.W], 'bicubic');
                end
                rr = Y(t)+(1:SP.H); cc = X(t)+(1:SP.W);
                num(rr,cc) = num(rr,cc) + P .* w2;
                den(rr,cc) = den(rr,cc) + w2;
            end
            den(den==0) = 1;
            tw.setTag(tag);
            tw.write(uint16(max(0, min(65535, round(num ./ den)))));
            if zi < numel(S.zPages), tw.writeDirectory(); end
        end
    end

    function M = render_full(ch, x, y, Ht, Wt)
        X = round(x); Y = round(y); X = X - min(X); Y = Y - min(Y);
        Wo = max(X) + Wt; Ho = max(Y) + Ht;
        w2 = feather(Ht, Wt, max(1, S.feather_um / max(1, scale_um())));
        num = zeros(Ho, Wo, 'single'); den = zeros(Ho, Wo, 'single');
        for t = 1:S.nTiles
            f = swap_ch(S.files(t), S.refChID, ch);
            [a, m] = read_proj(f, S.zPages, S.H0, S.W0);
            if strcmp(S.projMode,'max'), P = m; else, P = a; end
            if Ht ~= S.H0 || Wt ~= S.W0, P = imresize(P, [Ht Wt], 'bicubic'); end
            rr = Y(t)+(1:Ht); cc = X(t)+(1:Wt);
            num(rr,cc) = num(rr,cc) + P .* w2;
            den(rr,cc) = den(rr,cc) + w2;
        end
        den(den==0) = 1;  M = num ./ den;
    end

    function closeGui(s)
        if ~isempty(S.idleTimer) && isvalid(S.idleTimer)
            stop(S.idleTimer); delete(S.idleTimer);
        end
        if ~isempty(seamFig) && isvalid(seamFig), delete(seamFig); end
        delete(s);
    end
end

%% ============================ HELPERS =============================
function C = rg(A, B)
% Red = selected tile, green = its neighbour; yellow where they agree.
% Each is normalised to its OWN percentiles: on a shared window the brighter
% tile paints the whole overlap its own colour and says nothing about alignment.
C = cat(3, norm01(A), norm01(B), zeros(size(A),'single'));
end

function y = norm01(X)
v = single(X(:));  p = prctile(v, [2 99.5]);
if p(2) <= p(1), p(2) = p(1) + 1; end
y = max(0, min(1, (single(X) - p(1)) / (p(2) - p(1))));
end

function W = feather(h, w, fpx)
cx = min(1:w, w:-1:1); cy = min(1:h, h:-1:1);
W = max(min(min(1,single(cy(:))/fpx), min(1,single(cx(:).')/fpx)), 1e-3);
end

function [col, row] = parse_grid(files)
n = numel(files); col = nan(n,1); row = nan(n,1);
for k = 1:n
    [~, bn] = fileparts(files(k));
    a = regexp(bn,'col(\d+)','tokens','once');  if ~isempty(a), col(k)=str2double(a{1}); end
    b = regexp(bn,'row(\d+)','tokens','once');  if ~isempty(b), row(k)=str2double(b{1}); end
end
if any(isnan(col))     % tile_### fallback: index only, no row/col grouping
    for k = 1:n
        [~, bn] = fileparts(files(k));
        a = regexp(bn,'(?i)tile[_-]?(\d+)','tokens','once');
        if ~isempty(a), col(k) = str2double(a{1}); row(k) = 1; end
    end
end
end

function b = base_noch(files)
b = strings(numel(files),1);
for k = 1:numel(files)
    [~, nm] = fileparts(files(k));
    b(k) = regexprep(nm, '_ch\d+_avgz$', '');
end
end

function f = swap_ch(f, fromCh, toCh)
if fromCh == toCh, return; end
f = replace(f, sprintf("ch%d", fromCh), sprintf("ch%d", toCh));
end

function [H, W, cls] = first_page(fn)
t = Tiff(fn,'r'); A = t.read(); t.close();
H = size(A,1); W = size(A,2); cls = class(A);
end

function n = count_pages(fn)
t = Tiff(fn,'r'); n = 1;
while ~t.lastDirectory(), t.nextDirectory(); n = n + 1; end
t.close();
end

function [A, M, XZ] = read_proj(fn, pages, H, W)
t = Tiff(fn,'r'); c = onCleanup(@() t.close());
A = zeros(H,W,'single'); M = -inf(H,W,'single'); XZ = zeros(numel(pages), W, 'single');
for i = 1:numel(pages)
    t.setDirectory(pages(i));
    P = single(t.read());
    A = A + P;  M = max(M, P);  XZ(i,:) = max(P, [], 1);
end
A = A / max(1,numel(pages));
M(~isfinite(M)) = 0;
end

function [pages, how] = resolve_z(root, refTile, nZ, zUse_um, zUseProj)
if isempty(zUse_um)
    if isempty(zUseProj)
        pages = 1:nZ; how = sprintf('all %d pages', nZ); return
    end
    pages = zUseProj(:).'; pages = pages(pages>=1 & pages<=nZ);
    assert(~isempty(pages), 'zUseProj empty after clipping to 1..%d', nZ);
    how = sprintf('pages %d-%d (indices)', min(pages), max(pages)); return
end
zs = raw_zs(root, refTile);
assert(~isempty(zs), ['zUse_um is in MICRONS and needs SI.hStackManager.zs from a ' ...
    'raw tile; none found in %s. Use zUseProj instead.'], root);
assert(numel(zs) == nZ, 'raw header lists %d slices but the avgz has %d pages', numel(zs), nZ);
pages = find(zs >= min(zUse_um) & zs <= max(zUse_um)).';
assert(~isempty(pages), 'no slice in [%g %g] um', min(zUse_um), max(zUse_um));
how = sprintf('%g..%g um -> pages %d-%d', min(zUse_um), max(zUse_um), min(pages), max(pages));
end

function zs = raw_zs(root, refTile)
zs = [];
c = dir(fullfile(root,"*.tif"));
if isempty(c), return; end
files = string(fullfile({c.folder},{c.name})); names = string({c.name});
pick = 1;
[~, bn] = fileparts(refTile);
stem = regexprep(bn, '_ch\d+_avgz$', '');
hit = find(startsWith(names, stem + "."), 1);
if ~isempty(hit), pick = hit; end
zs = si_field(files(pick), "SI.hStackManager.zs");
end

function px = header_um_per_px(root, refTile)
px = NaN;
c = dir(fullfile(root,"*.tif"));
if isempty(c), return; end
f = string(fullfile(c(1).folder, c(1).name));
[~, bn] = fileparts(refTile);
stem = regexprep(bn, '_ch\d+_avgz$', '');
hit = find(startsWith(string({c.name}), stem + "."), 1);
if ~isempty(hit), f = string(fullfile(c(hit).folder, c(hit).name)); end
n  = si_field(f, "SI.hRoiManager.pixelsPerLine");
fv = si_field(f, "SI.hRoiManager.imagingFovUm");
if numel(fv) >= 4 && ~isempty(n)
    px = (max(fv(1:2:end)) - min(fv(1:2:end))) / n(1);
end
end

function v = si_field(f, key)
v = [];
try
    t = Tiff(f,'r'); s = string(t.getTag('Software')); t.close();
catch
    return
end
tok = regexp(s, key + "\s*=\s*([^\n\r]*)", 'tokens', 'once');
if isempty(tok), return; end
raw = erase(string(tok{1}), ["[","]","{","}","'"]);
raw = replace(replace(raw, ";", " "), ",", " ");
d = str2double(split(strtrim(raw)));
v = d(~isnan(d));
end

function [sx, sy, ok, how] = stage_positions(files, root)
n = numel(files); sx = nan(n,1); sy = nan(n,1); how = 'filename';
c = dir(fullfile(root,"*.tif"));
if numel(c) >= n
    names = string({c.name}); paths = string(fullfile({c.folder},{c.name}));
    good = true;
    for k = 1:n
        [~, bn] = fileparts(files(k));
        stem = regexprep(bn, '_ch\d+_avgz$', '');
        hit = find(startsWith(names, stem + "."), 1);
        if isempty(hit), good = false; break; end
        mp = si_field(paths(hit), "SI.hMotors.motorPosition");
        if numel(mp) < 2, good = false; break; end
        sx(k) = mp(1); sy(k) = mp(2);
    end
    if good, ok = true; how = 'header'; return; end
    sx(:) = NaN; sy(:) = NaN;
end
for k = 1:n
    [~, bn] = fileparts(files(k));
    a = regexp(bn,'_x(-?\d+(\.\d+)?)','tokens','once'); if ~isempty(a), sx(k)=str2double(a{1}); end
    b = regexp(bn,'_y(-?\d+(\.\d+)?)','tokens','once'); if ~isempty(b), sy(k)=str2double(b{1}); end
end
ok = all(isfinite(sx)) && all(isfinite(sy));
end
