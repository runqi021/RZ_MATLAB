function stack_stitch_gui_260910()
%STACK_STITCH_GUI_260910  Manual stitcher for a few overlapping 3-D stacks.
%
%   stack_stitch_gui_260910
%
% For the case the tile-grid stitcher does not fit: a HANDFUL of big volumes
% that overlap, already resampled to 1 um isotropic by struct_stack_prep_260910.
% Stage coordinates place them; you nudge; it saves the coords; and
% stack_stitch_apply_260910 renders the full-resolution mosaic for any channel.
%
% WHY NOT stitch_manual_gui_fast_260818.  That one is built for a lattice of
% many tiles indexed by col/row, each reduced to one avg-projection image. Here
% there are three volumes, no lattice, and z matters -- so this shows an X-Y
% mosaic AND an X-Z mosaic and lets you nudge in all three axes. The maths and
% the save-then-apply split are the same idea.
%
% ONE MICRON IS ONE PIXEL. The inputs are 1 um isotropic, so every offset in
% this GUI is simultaneously pixels and microns, and a z nudge is a whole plane.
% That is deliberate: chasing a sub-micron z offset by interpolating would blur
% every tile to correct a shift smaller than the axial PSF.
%
% PLACEMENT COMES FROM THE STAGE, and the sign convention is measured, not
% guessed: image COLUMN runs along +stage x, image ROW along -stage y, so
%       x_px = (stageX - min stageX) / umPerPx
%       y_px = (max stageY - stageY) / umPerPx
% There is deliberately NO flip control. A flip that "fixes" the orientation
% here would mirror the whole mosaic; if the strip comes out reversed, the
% convention is what is wrong and it needs measuring again, not a checkbox.
%
% CONTROLS
%   click a tile          select it (its border turns yellow)
%   arrows                nudge selected in x/y by Step (Shift = 10x)
%   [  and  ]             nudge selected in z by Step planes
%   1..9                  select tile by number
%   r                     reset the selected tile's nudge;  R resets all
%   scroll / right-drag   zoom / pan whichever view the cursor is over
%   a                     auto brightness/contrast
%   v                     open the seam (overlap) window
%   channel popup         switch the displayed channel.  It lists the SI channel
%                         numbers the tiles actually saved, intersected across
%                         tiles -- not a fixed 1/2/3.  Display only: the apply
%                         step stitches every shared channel from the one set of
%                         coords, so curation is never per-channel.
%   s                     save coords
%
% THE SEAM WINDOW is where the alignment is actually judged. For the selected
% tile it shows, per overlapping neighbour, the overlap region from each tile in
% grey and then a RED/GREEN overlay of the two -- in X-Y and again in X-Z.
% Aligned structure goes yellow; a misregistration splits every soma into a red
% ghost and a green one, and which side the red sits on tells you the direction
% to nudge. Far more sensitive than looking for a seam in the blended mosaic,
% because the feather blend is designed to HIDE exactly that.
%
% Runqi Zhang / 2026-09-10

%% ========================= USER SETTINGS =========================
S.root      = 'C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer';
S.stems     = {};            % {} = every <stem>_stackmeta.mat under root, sorted
S.refCh     = 1;             % channel shown in the GUI (both get stitched later)
S.projMode  = 'max';         % 'max' | 'avg' -- preview projection
S.zUse_um   = [];            % depth slab for the X-Y preview, um, [] = all
S.feather   = 30;            % blend ramp, um (preview AND apply)
S.moveStep  = 2;             % nudge, px (= um)
S.clipLo    = 1;             % Auto B/C percentiles
S.clipHi    = 99.5;
%% =================================================================

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);

S.n = 0; S.sel = 1; S.pan = [];
S.x = []; S.y = []; S.z = [];
S.manX = []; S.manY = []; S.manZ = [];
S.XY = {}; S.XZ = {}; S.names = {}; S.metaFiles = {};
S.tileH = []; S.tileW = []; S.tileD = [];
S.clim = [0 1];

%% ========================= FIGURE =========================
f = figure('Name','stack_stitch_gui_260910','Color',[.13 .13 .13], ...
    'NumberTitle','off','Units','pixels','Position',[60 60 1450 900], ...
    'KeyPressFcn',@cb_key, 'WindowScrollWheelFcn',@cb_scroll, ...
    'WindowButtonDownFcn',@cb_btndown, 'WindowButtonMotionFcn',@cb_btnmotion, ...
    'WindowButtonUpFcn',@cb_btnup, 'MenuBar','none','ToolBar','none');

% The X-Z panel is as wide as the mosaic but 531 planes deep, and `axis image`
% keeps it square-pixelled, so it is HEIGHT-limited: give it nearly as much of
% the canvas as the X-Y view or it renders as a thumbnail.  The figure cannot
% simply be made taller -- a MATLAB canvas is capped at the screen, and asking
% for more silently returns less (see reference_matlab_opengl_export_stubs).
axXY = axes(f,'Units','normalized','Position',[0.04 0.47 0.80 0.45], ...
            'Color','k','XColor','w','YColor','w');
axXZ = axes(f,'Units','normalized','Position',[0.04 0.05 0.80 0.38], ...
            'Color','k','XColor','w','YColor','w');
hImXY = []; hImXZ = []; hBox = []; hBoxZ = [];
hStep = []; hFocus = []; hChan = [];
projCache = containers.Map();   % 'ch<N>_<mode>' -> {XY, XZ}, so switching back is instant
seamFig = [];      % the live overlap window
seamAx  = [];      % its axes, [nSeam x 4]
seamIm  = [];      % its images, so only CData is touched on a nudge
seamKey = '';      % which (tile, partner-set) the panel grid was built for

pnl = @(yy,txt,cb) uicontrol(f,'Style','pushbutton','Units','normalized', ...
        'Position',[0.86 yy 0.12 0.045],'String',txt,'Callback',cb,'KeyPressFcn',@cb_key);
hStat = uicontrol(f,'Style','text','Units','normalized','Position',[0.04 0.965 0.94 0.03], ...
        'BackgroundColor',[.13 .13 .13],'ForegroundColor','w','HorizontalAlignment','left', ...
        'FontName','Consolas','FontSize',9);
lbl = @(yy,txt) uicontrol(f,'Style','text','Units','normalized', ...
        'Position',[0.86 yy 0.12 0.028],'String',txt,'BackgroundColor',[.13 .13 .13], ...
        'ForegroundColor','w','HorizontalAlignment','left','FontSize',9);

% The right-hand column is laid out ONCE, top to bottom, with no gaps left for
% something to be dropped into later.  The step box and the seam button were
% first placed straight on top of hInfo; uicontrols do not complain about
% overlapping, the later one simply draws over the earlier, so the info panel
% was being half-hidden with nothing to show for it.
pnl(0.930, 'Auto B/C',        @(~,~) autoBC());
pnl(0.875, 'Reset tile (r)',  @(~,~) resetTile(S.sel));
pnl(0.820, 'Reset ALL (R)',   @(~,~) resetTile(0));

lbl(0.780, 'channel');
hChan = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.86 0.738 0.12 0.04], ...
        'String',{'(loading)'},'Callback',@cb_chan,'KeyPressFcn',@cb_key);

% Nudge step, typed.  Floored at 1 because the voxel IS a micron: a step below
% one pixel cannot move anything, it would just accumulate a fractional offset
% that round() throws away at composite time, so the tile would appear stuck.
lbl(0.692, 'nudge step (um)');
hStep = uicontrol(f,'Style','edit','Units','normalized','Position',[0.86 0.650 0.055 0.04], ...
        'String',num2str(S.moveStep),'BackgroundColor','w','Callback',@cb_step);
uicontrol(f,'Style','text','Units','normalized','Position',[0.921 0.647 0.06 0.032], ...
        'String','Shift=10x','BackgroundColor',[.13 .13 .13],'ForegroundColor',[.6 .6 .6], ...
        'HorizontalAlignment','left','FontSize',8);

pnl(0.590, 'Seam view (v)',   @(~,~) showSeams(true));
pnl(0.535, 'proj max/avg',    @(~,~) toggleProj());

hInfo = uicontrol(f,'Style','text','Units','normalized','Position',[0.86 0.255 0.13 0.26], ...
        'BackgroundColor',[.13 .13 .13],'ForegroundColor',[.8 .9 1],'HorizontalAlignment','left', ...
        'FontName','Consolas','FontSize',9);

pnl(0.195, 'Save coords (s)', @(~,~) cb_save());
pnl(0.140, 'Apply -> volume', @(~,~) cb_apply());

% Focus sink: after typing in the edit box the arrow keys belong to the caret,
% not the tiles.  Parking focus on an off-screen control that forwards to cb_key
% hands the keyboard back without the user having to click the image first.
hFocus = uicontrol(f,'Style','pushbutton','Units','pixels','Position',[-30 -30 1 1], ...
        'KeyPressFcn',@cb_key);

loadTiles();
redraw();
showSeams(true);        % the seam window is the point -- open it up front

%% ========================= LOAD =========================
    function loadTiles()
        stems = S.stems;
        if isempty(stems)
            d = dir(fullfile(S.root, '**', '*_stackmeta.mat'));
            assert(~isempty(d), 'no *_stackmeta.mat under %s -- run struct_stack_prep_260910 first', S.root);
            [~, ord] = sort({d.name});  d = d(ord);
            S.metaFiles = arrayfun(@(k) fullfile(d(k).folder, d(k).name), (1:numel(d)).', 'uni', 0);
        else
            S.metaFiles = cellfun(@(s) fullfile(S.root, s, [s '_stackmeta.mat']), stems(:), 'uni', 0);
        end
        S.n = numel(S.metaFiles);
        setStat(sprintf('loading %d tiles ...', S.n)); drawnow;

        stageX = zeros(S.n,1); stageY = zeros(S.n,1); d0 = zeros(S.n,1);
        for i = 1:S.n
            M = load(S.metaFiles{i}, 'meta'); m = M.meta;
            [~, S.names{i}] = fileparts(m.tifPath);
            if i == 1
                S.pxOut = m.pxOut; S.channels = m.channels(:).';
            else
                assert(abs(m.pxOut - S.pxOut) < 1e-9, 'tile %d has a different voxel size', i);
                % INTERSECTION, not the first tile's list: a channel one tile
                % never saved cannot be stitched, and offering it would only
                % fail later, in the apply step, after the curation is done.
                S.channels = intersect(S.channels, m.channels(:).');
            end
            stageX(i) = m.motorPosition(1);
            stageY(i) = m.motorPosition(2);
            d0(i)     = m.depth(1);
            S.tileH(i) = m.size(1); S.tileW(i) = m.size(2); S.tileD(i) = m.size(3);
        end
        assert(~isempty(S.channels), 'the tiles share no saved channel');

        % Offer only the SI channel numbers actually present.  Which numbers
        % those are is a property of the acquisition (channelSave), not a fixed
        % 1/2/3 list -- this session saved [1 3], so 2 must not appear.
        if ~ismember(S.refCh, S.channels), S.refCh = S.channels(1); end
        set(hChan, 'String', arrayfun(@(c) sprintf('ch %d', c), S.channels, 'uni', 0), ...
                   'Value', find(S.channels == S.refCh, 1));
        loadProjections();

        % Measured rig convention: column -> +stage x, row -> -stage y.
        S.x = (stageX - min(stageX)) / S.pxOut;
        S.y = (max(stageY) - stageY) / S.pxOut;
        S.z = (d0 - min(d0)) / S.pxOut;          % depth increases downward
        S.manX = zeros(S.n,1); S.manY = zeros(S.n,1); S.manZ = zeros(S.n,1);
        S.stageX = stageX; S.stageY = stageY;

        prior = fullfile(S.root, 'matlab_stitch', 'stack_stitch_coords.mat');
        if isfile(prior)
            q = questdlg(sprintf('Load curated coords from\n%s ?', prior), ...
                         'stack_stitch', 'Load', 'Start from stage', 'Load');
            if strcmp(q, 'Load')
                C = load(prior, 'stitch');
                if numel(C.stitch.manX) == S.n
                    S.manX = C.stitch.manX(:); S.manY = C.stitch.manY(:); S.manZ = C.stitch.manZ(:);
                end
            end
        end
        S.sel = 1;
    end

%% ========================= MOSAIC =========================
    function [M, org] = composite(P, px, py, hh, ww, fx, fy)
        % Weighted blend of the tile images P at integer offsets, in one canvas.
        % Everything is forced to a COLUMN first: the offsets are columns but the
        % tile sizes are collected as rows, and `ox + ww` on mixed orientations
        % broadcasts to an n-by-n matrix instead of erroring, so the canvas size
        % silently comes out as a vector.
        px = px(:); py = py(:); hh = hh(:); ww = ww(:);
        ox = round(px - min(px));  oy = round(py - min(py));
        W = max(ox + ww);  H = max(oy + hh);
        num = zeros(H, W, 'single');  den = zeros(H, W, 'single');
        for k = 1:numel(P)
            wk = feather2(hh(k), ww(k), fx, fy);
            rr = oy(k)+(1:hh(k));  cc = ox(k)+(1:ww(k));
            num(rr,cc) = num(rr,cc) + single(P{k}) .* wk;
            den(rr,cc) = den(rr,cc) + wk;
        end
        den(den == 0) = 1;
        M = num ./ den;
        org = [min(px) min(py)];
    end

    function redraw()
        fpx = max(1, S.feather / S.pxOut);
        [MXY, ~] = composite(S.XY, S.x+S.manX, S.y+S.manY, S.tileH, S.tileW, fpx, fpx);
        % X-Z mosaic: rows are depth, columns are the same x axis as above.  The
        % blend ramps in x ONLY here (fy = 1 -> flat): a depth ramp would fade
        % the top and bottom of the side view, which is exactly the part being
        % judged when nudging z.
        [MXZ, ~] = composite(S.XZ, S.x+S.manX, S.z+S.manZ, S.tileD, S.tileW, fpx, 1);

        if isempty(hImXY) || ~isvalid(hImXY)
            hImXY = imagesc(axXY, MXY); axis(axXY,'image'); colormap(axXY, gray);
            set(axXY,'XColor','w','YColor','w','Box','on');
            title(axXY, 'X-Y   (click a tile to select)', 'Color','w','FontWeight','normal');
            hold(axXY,'on'); hBox = gobjects(0); hold(axXY,'off');
            hImXZ = imagesc(axXZ, MXZ); axis(axXZ,'image'); colormap(axXZ, gray);
            set(axXZ,'XColor','w','YColor','w','Box','on');
            title(axXZ, 'X-Z   (depth down; [ and ] nudge z)', 'Color','w','FontWeight','normal');
            autoBC();
        else
            set(hImXY,'CData',MXY); set(hImXZ,'CData',MXZ);
            set(axXY,'XLimMode','manual'); % keep the user's zoom
        end
        set(hImXY,'XData',[1 size(MXY,2)],'YData',[1 size(MXY,1)]);
        set(hImXZ,'XData',[1 size(MXZ,2)],'YData',[1 size(MXZ,1)]);
        caxis(axXY, S.clim); caxis(axXZ, S.clim);
        drawBoxes();
        updateSeams();
        setInfo();
    end

    function drawBoxes()
        delete(hBox(isgraphics(hBox)));  delete(hBoxZ(isgraphics(hBoxZ)));
        hBox = gobjects(S.n,1); hBoxZ = gobjects(S.n,1);
        ox = round(S.x + S.manX - min(S.x + S.manX));
        oy = round(S.y + S.manY - min(S.y + S.manY));
        oz = round(S.z + S.manZ - min(S.z + S.manZ));
        hold(axXY,'on'); hold(axXZ,'on');
        for k = 1:S.n
            col = [0.35 0.75 1]; lw = 0.75;
            if k == S.sel, col = [1 0.85 0.1]; lw = 1.75; end
            hBox(k)  = rectangle(axXY,'Position',[ox(k)+0.5 oy(k)+0.5 S.tileW(k) S.tileH(k)], ...
                                 'EdgeColor',col,'LineWidth',lw);
            hBoxZ(k) = rectangle(axXZ,'Position',[ox(k)+0.5 oz(k)+0.5 S.tileW(k) S.tileD(k)], ...
                                 'EdgeColor',col,'LineWidth',lw);
        end
        hold(axXY,'off'); hold(axXZ,'off');
    end

%% ========================= INTERACTION =========================
    function cb_key(~, e)
        sh = any(strcmp(e.Modifier, 'shift'));
        st = S.moveStep * (1 + 9*sh);
        switch e.Key
            case 'leftarrow',  S.manX(S.sel) = S.manX(S.sel) - st; redraw();
            case 'rightarrow', S.manX(S.sel) = S.manX(S.sel) + st; redraw();
            case 'uparrow',    S.manY(S.sel) = S.manY(S.sel) - st; redraw();
            case 'downarrow',  S.manY(S.sel) = S.manY(S.sel) + st; redraw();
            case 'leftbracket',  S.manZ(S.sel) = S.manZ(S.sel) - st; redraw();
            case 'rightbracket', S.manZ(S.sel) = S.manZ(S.sel) + st; redraw();
            case 'a', autoBC(); redraw();
            case 'v', showSeams(true);
            case 's', cb_save();
            case 'r'
                if sh, resetTile(0); else, resetTile(S.sel); end
            otherwise
                k = str2double(e.Key);
                if ~isnan(k) && k >= 1 && k <= S.n, S.sel = k; redraw(); end
        end
    end

    function cb_btndown(~, ~)
        switch get(f, 'SelectionType')
            case 'normal'
                cp = get(axXY, 'CurrentPoint');
                px = cp(1,1); py = cp(1,2);
                ox = round(S.x + S.manX - min(S.x + S.manX));
                oy = round(S.y + S.manY - min(S.y + S.manY));
                % Topmost match wins, so overlapping tiles stay selectable.
                for k = S.n:-1:1
                    if px >= ox(k) && px <= ox(k)+S.tileW(k) && py >= oy(k) && py <= oy(k)+S.tileH(k)
                        S.sel = k; redraw(); return
                    end
                end
            case 'alt'
                ax = axUnderPointer();
                S.pan = struct('ax', ax, 'p0', get(ax,'CurrentPoint'), ...
                               'xl', xlim(ax), 'yl', ylim(ax));
        end
    end

    function cb_btnmotion(~, ~)
        if isempty(S.pan), return; end
        cp = get(S.pan.ax, 'CurrentPoint');
        dx = S.pan.p0(1,1) - cp(1,1);  dy = S.pan.p0(1,2) - cp(1,2);
        xlim(S.pan.ax, S.pan.xl + dx);  ylim(S.pan.ax, S.pan.yl + dy);
    end

    function cb_btnup(~, ~), S.pan = []; end

    function ax = axUnderPointer()
        % Zoom and pan act on whichever view the cursor is over, so the X-Z
        % panel is inspectable at the seams without a second set of controls.
        cp = get(f, 'CurrentPoint');  fp = get(f, 'Position');
        rel = [cp(1)/fp(3), cp(2)/fp(4)];
        pz  = get(axXZ, 'Position');
        if rel(2) >= pz(2) && rel(2) <= pz(2)+pz(4), ax = axXZ; else, ax = axXY; end
    end

    function cb_scroll(~, e)
        ax = axUnderPointer();
        cp = get(ax,'CurrentPoint'); cx = cp(1,1); cy = cp(1,2);
        fac = 1.2 ^ double(e.VerticalScrollCount);
        xl = xlim(ax); yl = ylim(ax);
        xlim(ax, cx + (xl - cx)*fac);  ylim(ax, cy + (yl - cy)*fac);
    end

    function cb_step(~, ~)
        v = str2double(get(hStep, 'String'));
        if ~isfinite(v) || v < 1
            v = max(1, S.moveStep);
            setStat('step must be >= 1 um -- the voxel is 1 um, a smaller step moves nothing');
        end
        S.moveStep = round(v);
        set(hStep, 'String', num2str(S.moveStep));
        uicontrol(hFocus);              % hand the arrow keys back to the tiles
        setInfo();
    end

%% ========================= SEAM (OVERLAP) VIEW =========================
    function P = seamPairs(k)
        % Every tile that actually overlaps tile k, with the intersection
        % expressed in BOTH tiles' local pixel coordinates.  Recomputed on each
        % redraw because a nudge changes the overlap.
        ox = round(S.x + S.manX);  oy = round(S.y + S.manY);  oz = round(S.z + S.manZ);
        P = struct('j',{},'rk',{},'ck',{},'rj',{},'cj',{},'zk',{},'zj',{},'ov',{});
        for j = 1:S.n
            if j == k, continue; end
            % x/y intersection in canvas coordinates
            c0 = max(ox(k), ox(j));  c1 = min(ox(k)+S.tileW(k), ox(j)+S.tileW(j)) - 1;
            r0 = max(oy(k), oy(j));  r1 = min(oy(k)+S.tileH(k), oy(j)+S.tileH(j)) - 1;
            if c1 < c0 || r1 < r0, continue; end
            % z intersection, for the depth overlay
            z0 = max(oz(k), oz(j));  z1 = min(oz(k)+S.tileD(k), oz(j)+S.tileD(j)) - 1;
            if z1 < z0, z0 = oz(k); z1 = oz(k)+S.tileD(k)-1; end
            P(end+1) = struct('j', j, ...
                'rk', (r0:r1) - oy(k) + 1, 'ck', (c0:c1) - ox(k) + 1, ...
                'rj', (r0:r1) - oy(j) + 1, 'cj', (c0:c1) - ox(j) + 1, ...
                'zk', (z0:z1) - oz(k) + 1, 'zj', (z0:z1) - oz(j) + 1, ...
                'ov', (c1-c0+1)*(r1-r0+1)); %#ok<AGROW>
        end
        [~, ord] = sort([P.ov], 'descend');  P = P(ord);
    end

    function showSeams(force)
        if isempty(seamFig) || ~isvalid(seamFig)
            if ~force, return; end        % only opened on demand, or kept alive
            seamFig = figure('Name','stack_stitch seam view','NumberTitle','off', ...
                'Color',[.13 .13 .13],'Units','pixels','Position',[80 80 1150 620], ...
                'MenuBar','none','ToolBar','none','KeyPressFcn',@cb_key);
            seamKey = '';
        end
        updateSeams();
    end

    function updateSeams()
        if isempty(seamFig) || ~isvalid(seamFig), return; end
        % Name carries the channel: the overlay is the one panel where you could
        % not otherwise tell which label you are judging the alignment on.
        set(seamFig, 'Name', sprintf('stack_stitch seam view  |  ch %d  (%s)', S.refCh, S.projMode));
        P = seamPairs(S.sel);
        key = sprintf('%d:%s', S.sel, mat2str([P.j]));
        nS = numel(P);
        if nS == 0
            clf(seamFig);
            annotation(seamFig,'textbox',[0.1 0.45 0.8 0.1],'String', ...
                sprintf('tile %d overlaps nothing at the current placement', S.sel), ...
                'Color','w','EdgeColor','none','HorizontalAlignment','center');
            seamKey = key; seamAx = []; seamIm = [];
            return
        end

        % Rebuild the panel grid only when the seam SET changes; a nudge just
        % pushes new CData, which is what keeps this live rather than flickery.
        if ~strcmp(key, seamKey)
            clf(seamFig);
            tl = tiledlayout(seamFig, nS, 4, 'TileSpacing','compact','Padding','compact');
            seamAx = gobjects(nS,4); seamIm = gobjects(nS,4);
            for s = 1:nS
                for q = 1:4
                    seamAx(s,q) = nexttile(tl, (s-1)*4 + q);
                    seamIm(s,q) = imagesc(seamAx(s,q), 0);
                    axis(seamAx(s,q), 'image');
                    set(seamAx(s,q), 'XTick', [], 'YTick', [], 'Color','k', ...
                        'XColor',[.5 .5 .5], 'YColor',[.5 .5 .5], 'Box','on');
                    if q <= 2, colormap(seamAx(s,q), gray); end
                end
            end
            seamKey = key;
        end

        for s = 1:nS
            p = P(s);  k = S.sel;  j = p.j;
            A  = S.XY{k}(p.rk, p.ck);   B  = S.XY{j}(p.rj, p.cj);
            Az = S.XZ{k}(p.zk, p.ck);   Bz = S.XZ{j}(p.zj, p.cj);

            set(seamIm(s,1), 'CData', A);  caxis(seamAx(s,1), S.clim);
            set(seamIm(s,2), 'CData', B);  caxis(seamAx(s,2), S.clim);
            % RED = selected tile, GREEN = partner.  Aligned structure goes
            % YELLOW; a misregistration splits every soma into a red ghost and a
            % green one, and the split direction tells you which way to nudge.
            set(seamIm(s,3), 'CData', rg(A, B, S.clim));
            set(seamIm(s,4), 'CData', rg(Az, Bz, S.clim));
            for q = 1:4
                sz = size(get(seamIm(s,q),'CData'));
                set(seamIm(s,q), 'XData', [1 sz(2)], 'YData', [1 sz(1)]);
                xlim(seamAx(s,q), [0.5 sz(2)+0.5]); ylim(seamAx(s,q), [0.5 sz(1)+0.5]);
            end
            title(seamAx(s,1), sprintf('tile %d  (red)', k), 'Color',[1 .5 .5], 'FontWeight','normal');
            title(seamAx(s,2), sprintf('tile %d  (green)', j), 'Color',[.5 1 .5], 'FontWeight','normal');
            title(seamAx(s,3), sprintf('X-Y overlap  %d x %d px', numel(p.ck), numel(p.rk)), ...
                  'Color','w', 'FontWeight','normal');
            title(seamAx(s,4), sprintf('X-Z overlap  dz %+g', S.manZ(j)-S.manZ(k)), ...
                  'Color','w', 'FontWeight','normal');
        end
    end

    function resetTile(k)
        if k == 0
            S.manX(:) = 0; S.manY(:) = 0; S.manZ(:) = 0;
        else
            S.manX(k) = 0; S.manY(k) = 0; S.manZ(k) = 0;
        end
        redraw();
    end

    function loadProjections()
        % Projections for the current channel + projection mode, cached.  They
        % are 2-D, so holding every combination costs a few MB, and switching
        % channel to compare the label against the vasculature is then instant
        % instead of another pass over 531 pages per tile.
        key = sprintf('ch%d_%s', S.refCh, S.projMode);
        if isKey(projCache, key)
            c = projCache(key);  S.XY = c{1};  S.XZ = c{2};  return
        end
        for i = 1:S.n
            setStat(sprintf('projecting ch%d (%s)  %d/%d  %s', S.refCh, S.projMode, i, S.n, S.names{i}));
            drawnow;
            M = load(S.metaFiles{i}, 'meta');
            [S.XY{i}, S.XZ{i}] = tile_projections(M.meta, S.refCh, S.zUse_um, S.projMode);
        end
        projCache(key) = {S.XY, S.XZ};
    end

    function cb_chan(~, ~)
        S.refCh = S.channels(get(hChan, 'Value'));
        loadProjections();
        seamKey = '';          % force the seam grid to rebuild on the new channel
        autoBC(); redraw();
        uicontrol(hFocus);
        setStat(sprintf('channel %d', S.refCh));
    end

    function toggleProj()
        if strcmp(S.projMode,'max'), S.projMode = 'avg'; else, S.projMode = 'max'; end
        loadProjections();
        autoBC(); redraw();
    end

    function autoBC()
        v = [];
        for k = 1:S.n, v = [v; single(reshape(S.XY{k}(1:4:end,1:4:end), [], 1))]; end %#ok<AGROW>
        S.clim = double(prctile(v, [S.clipLo S.clipHi]));
        if S.clim(2) <= S.clim(1), S.clim(2) = S.clim(1) + 1; end
        if ~isempty(hImXY) && isvalid(hImXY), caxis(axXY, S.clim); caxis(axXZ, S.clim); end
    end

    function setStat(t), set(hStat, 'String', t); end

    function setInfo()
        L = {sprintf('ch %d  proj %s', S.refCh, S.projMode), ...
             sprintf('voxel %g um', S.pxOut), sprintf('feather %g um', S.feather), ...
             sprintf('step %g um', S.moveStep), ''};
        for k = 1:S.n
            mark = '  '; if k == S.sel, mark = '> '; end
            L{end+1} = sprintf('%s%d %s', mark, k, S.names{k}); %#ok<AGROW>
            L{end+1} = sprintf('   dx%+5g dy%+5g dz%+5g', S.manX(k), S.manY(k), S.manZ(k)); %#ok<AGROW>
        end
        set(hInfo, 'String', L);
        setStat(sprintf(['tile %d/%d  %s   |   nudge  dx %+g  dy %+g  dz %+g px(um)   |   ' ...
                         'stage x %.1f y %.1f um'], S.sel, S.n, S.names{S.sel}, ...
                         S.manX(S.sel), S.manY(S.sel), S.manZ(S.sel), ...
                         S.stageX(S.sel), S.stageY(S.sel)));
    end

%% ========================= SAVE / APPLY =========================
    function cb_save()
        outDir = fullfile(S.root, 'matlab_stitch');
        if ~isfolder(outDir), mkdir(outDir); end
        stitch = struct('root', S.root, 'metaFiles', {S.metaFiles}, 'names', {S.names}, ...
            'x', S.x, 'y', S.y, 'z', S.z, 'manX', S.manX, 'manY', S.manY, 'manZ', S.manZ, ...
            'tileH', S.tileH, 'tileW', S.tileW, 'tileD', S.tileD, 'pxOut', S.pxOut, ...
            'channels', S.channels, 'stageX', S.stageX, 'stageY', S.stageY, ...
            'feather', S.feather, 'projMode', S.projMode, 'refCh', S.refCh, ...
            'savedOn', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
        % manZ is folded into z at apply time, so the saved base placement stays
        % exactly what the stage said and the curation stays inspectable.
        fn = fullfile(outDir, 'stack_stitch_coords.mat');
        save(fn, 'stitch');
        setStat(sprintf('saved %s', fn));
    end

    function cb_apply()
        cb_save();
        fn = fullfile(S.root, 'matlab_stitch', 'stack_stitch_coords.mat');
        setStat('rendering the mosaic -- see the command window'); drawnow;
        stack_stitch_apply_260910(fn, 'Feather', S.feather);
        setStat('mosaic written to matlab_stitch\');
    end
end

% =============================================================================
function [XY, XZ] = tile_projections(meta, ch, zUse_um, mode)
%TILE_PROJECTIONS  One X-Y and one X-Z projection of a tile, read once.
%
% The X-Z view is what makes a z nudge checkable: two tiles can look perfectly
% registered in X-Y and still be a slab apart in depth, and an X-Y preview
% cannot show that at all.

ci = find(meta.channels == ch, 1);
assert(~isempty(ci), 'channel %d is not in this tile', ch);
tif = meta.outTif{ci};
if ~isfile(tif)
    [~, nm, ex] = fileparts(tif);
    tif = fullfile(fileparts(meta.tifPath), [nm ex]);
end
assert(isfile(tif), 'missing tile stack %s', tif);

depth = meta.depth(:);
if isempty(zUse_um)
    keep = true(size(depth));
else
    keep = depth >= min(zUse_um) & depth <= max(zUse_um);
end

t = Tiff(tif, 'r'); c = onCleanup(@() close(t));
A = t.read();  [h, w] = size(A);
nz = numel(depth);
XY = -inf(h, w, 'single');  if strcmpi(mode,'avg'), XY = zeros(h, w, 'single'); end
XZ = zeros(nz, w, 'single');
nk = 0;
t.setDirectory(1);
for k = 1:nz
    P = single(t.read());
    XZ(k,:) = max(P, [], 1);              % X-Z is always a MAX along rows
    if keep(k)
        if strcmpi(mode,'avg'), XY = XY + P; else, XY = max(XY, P); end
        nk = nk + 1;
    end
    if k < nz, t.nextDirectory(); end
end
if strcmpi(mode,'avg') && nk > 0, XY = XY / nk; end
if ~isfinite(XY(1)), XY(~isfinite(XY)) = 0; end
end

function R = feather2(h, w, fx, fy)
% Separable ramp: 1 in the middle, falling to ~0 within fx/fy pixels of the
% edge.  fx = fy = 1 gives a flat weight, i.e. a plain average in the overlap.
cx = min(1:w, w:-1:1);  cy = min(1:h, h:-1:1);
R = max(min(min(1, single(cy(:))/fy), min(1, single(cx(:).')/fx)), 1e-3);
end

function C = rg(A, B, ~)
%RG  Two images as the red and green channels of one RGB, each on its OWN window.
%
% Per-image normalisation, deliberately, and it is the opposite of what the grey
% panels beside it do.  Those keep the shared display window so a real
% brightness difference between the tiles stays visible.  This panel has one job
% -- is the structure in the same PLACE -- and on a shared window the brighter
% tile simply paints the whole overlap its own colour: fov2 outshines fov1 here
% and the overlay came out uniformly green, which says nothing about alignment.
% Equalising them first is what makes "yellow = registered" actually work.
a = norm01(A);  b = norm01(B);
C = cat(3, a, b, zeros(size(a), 'single'));
end

function y = norm01(X)
v = single(X(:));
p = prctile(v, [2 99.5]);
if p(2) <= p(1), p(2) = p(1) + 1; end
y = max(0, min(1, (single(X) - p(1)) / (p(2) - p(1))));
end
