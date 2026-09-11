function stitch_manual_gui()
%STITCH_MANUAL_GUI  Manual 2D tile stitcher for ScanImage overview maps (ch1).
%
%   Built for high-zoom maps where tile overlap is too small (~12%) for the
%   phase-correlation auto-stitcher. Tiles are placed on a REGULAR GRID by
%   their (col,row) index times an adjustable STEP SIZE (in pixels). You set
%   the step to line up the seams, then fine-tune individual tiles by hand.
%
%   * One channel only (ch1). Each tile = the AVG PROJECTION of its avgz
%     stack (mean over Z pages -> one 2D image). No single-Z browsing.
%   * The step that spaces all tiles is NOT a GUI control: it comes from
%     S.stepX_um / S.stepY_um + S.zoom below. It is only auto-derived from
%     the filename micron coords when the names lack colXX/rowYY indices.
%     Flip X / Flip Y fix a mirrored mosaic.
%   * Click a tile, then arrow keys nudge it (1 px; Shift = 10 px).
%     Move mode (Tile/Row/Col) shifts one tile, a row, or a column.
%   * Live preview is DOWNSAMPLED for speed; full resolution is rendered
%     only on Save.
%
%   Output (in <root>/matlab_stitch/):
%     stitched_ch1_avgproj.tif   -- full-res manually-aligned mosaic
%     stitch_manual_coords.mat   -- step, flips, positions, nudges (reloadable)
%
%   Standalone. MATLAB R2021a+. Runqi Zhang / 2026.

%% ========================= USER SETTINGS =========================
S.datasetRoot = "C:\260804_shiver_dbh\map";  % has ch1/avgz, ch3/avgz
S.refChID     = 3;            % channel to stitch

S.umPerPxBase = 1.7778;       % microscope constant @ zoom 1
S.zoom        = 1.2;            % -> 1.4815 um/px
S.stepX_um    = 600;          % stage step (um) -> default step in px
S.stepY_um    = 600;

S.zUseProj    = [];           % Z pages for avg proj ([]=all)
S.flipX       = false;
S.flipY       = true;         % stage-Y usually opposite image rows
S.dsPreview   = 4;            % live-preview downsample factor (speed)

%% ========================= STATE =========================
umPerPx       = S.umPerPxBase / S.zoom;
S.stepXpx     = round(S.stepX_um / umPerPx);   % ~450
S.stepYpx     = round(S.stepY_um / umPerPx);

S.files   = strings(0,1);     % ch1 avgz tile paths
S.nTiles  = 0;
S.H0 = 0; S.W0 = 0; S.inClass = 'uint16';
S.colIdx = []; S.rowIdx = []; S.nRows = 0; S.nCols = 0;
S.tileImg   = {};             % full-res avg proj (double) per tile
S.tileImgDS = {};             % downsampled (single) per tile, for preview
S.Hd = 0; S.Wd = 0;           % downsampled tile size
S.manX = []; S.manY = [];     % manual nudges (px)
S.x = []; S.y = [];           % full-res top-left placement (px)
S.outW = 0; S.outH = 0;
S.dsR0 = []; S.dsC0 = [];     % per-tile DS top-left (row0,col0)
S.mosaic = [];                % downsampled preview mosaic
S.selectedTile = 0;
S.moveMode = 'Tile';
S.moveStep = 5;               % arrow-key nudge increment (px); Shift = x10
S.loaded = false;
S.coordDir = '';

%% ========================= FIGURE =========================
fig = uifigure('Name','Manual Stitch (ch1, avg proj)', ...
    'Position',[40 40 1500 950],'WindowState','maximized', ...
    'CloseRequestFcn',@(src,~) delete(src),'WindowKeyPressFcn',@cb_key);

rootGL = uigridlayout(fig,[1 2]);
rootGL.ColumnWidth = {280,'1x'}; rootGL.Padding=[4 4 4 4]; rootGL.ColumnSpacing=6;

ctrlPan = uipanel(rootGL,'Title','Controls','FontWeight','bold'); ctrlPan.Layout.Column=1;
nR=26; cGL=uigridlayout(ctrlPan,[nR 2]);
cGL.RowHeight=repmat({'fit'},1,nR); cGL.ColumnWidth={'fit','1x'};
cGL.Padding=[6 6 6 6]; cGL.RowSpacing=4; r=0;

r=r+1; h=uilabel(cGL,'Text','-- Load --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnLoad=uibutton(cGL,'Text','Load dataset (ch1/avgz)','ButtonPushedFcn',@cb_load);
       btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];
r=r+1; lblPath=uilabel(cGL,'Text',S.datasetRoot,'WordWrap','on','FontSize',9,'FontColor',[.5 .5 .5]);
       lblPath.Layout.Row=r; lblPath.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','-- Layout --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; cbFlipX=uicheckbox(cGL,'Text','Flip X','Value',S.flipX,'ValueChangedFcn',@cb_flip);
       cbFlipX.Layout.Row=r; cbFlipX.Layout.Column=1;
       cbFlipY=uicheckbox(cGL,'Text','Flip Y','Value',S.flipY,'ValueChangedFcn',@cb_flip);
       cbFlipY.Layout.Row=r; cbFlipY.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','-- Manual nudge --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; lblSel=uilabel(cGL,'Text','Selected: none (click a tile)','FontColor',[.1 .4 .7]);
       lblSel.Layout.Row=r; lblSel.Layout.Column=[1 2];
r=r+1; h=uilabel(cGL,'Text','Move:'); h.Layout.Row=r; h.Layout.Column=1;
       ddMove=uidropdown(cGL,'Items',{'Tile','Row','Col'},'Value','Tile', ...
           'ValueChangedFcn',@(s,~) assignMove(s.Value)); ddMove.Layout.Row=r; ddMove.Layout.Column=2;
r=r+1; h=uilabel(cGL,'Text','Move step (px):'); h.Layout.Row=r; h.Layout.Column=1;
       efMove=uieditfield(cGL,'numeric','Value',S.moveStep,'ValueDisplayFormat','%.0f', ...
           'Limits',[1 Inf],'ValueChangedFcn',@(s,~) assignMoveStep(s.Value));
       efMove.Layout.Row=r; efMove.Layout.Column=2;
r=r+1; h=uilabel(cGL,'Text','Arrow = step   Shift+Arrow = 10x step','FontSize',9,'FontColor',[.4 .4 .4]);
       h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnResetSel=uibutton(cGL,'Text','Reset this tile','ButtonPushedFcn',@cb_resetSel);
       btnResetSel.Layout.Row=r; btnResetSel.Layout.Column=1;
       btnResetAll=uibutton(cGL,'Text','Reset ALL nudges','ButtonPushedFcn',@cb_resetAll);
       btnResetAll.Layout.Row=r; btnResetAll.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','-- Display --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; cbSeams=uicheckbox(cGL,'Text','Show borders','Value',true,'ValueChangedFcn',@(~,~) drawOverlays());
       cbSeams.Layout.Row=r; cbSeams.Layout.Column=1;
       btnAuto=uibutton(cGL,'Text','Auto B/C','ButtonPushedFcn',@(~,~) autoBC());
       btnAuto.Layout.Row=r; btnAuto.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','-- Save --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnSave=uibutton(cGL,'Text','Save ch1 + coords','ButtonPushedFcn',@cb_save, ...
           'BackgroundColor',[.18 .55 .30],'FontColor','white');
       btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];
r=r+1; lblStat=uilabel(cGL,'Text','Ready. Click "Load dataset".','WordWrap','on','FontSize',10,'FontColor',[.3 .3 .3]);
       lblStat.Layout.Row=r; lblStat.Layout.Column=[1 2];
for ir=(r+1):nR, h2=uilabel(cGL,'Text',''); h2.Layout.Row=ir; h2.Layout.Column=[1 2]; end

ax=uiaxes(rootGL); ax.Layout.Column=2; colormap(ax,gray(256));
ax.XTick=[]; ax.YTick=[]; ax.DataAspectRatio=[1 1 1]; title(ax,'Manual stitch','Interpreter','none');
hImg=[]; hBord=[]; hSel=[];

%% ========================= CALLBACKS =========================
    function cb_load(~,~)
        d=uigetdir(char(S.datasetRoot),'Select dataset root (contains ch1/avgz)');
        if isequal(d,0), return; end
        S.datasetRoot=string(d); loadDataset();
    end

    function loadDataset()
        setStat('Loading ch1 tiles...'); drawnow;
        avgzDir=fullfile(char(S.datasetRoot),sprintf('ch%d',S.refChID),'avgz');
        if ~isfolder(avgzDir)
            uialert(fig,sprintf('Not found: %s',avgzDir),'Load error'); return;
        end
        files=sort_by_colrow(list_tifs(avgzDir));
        if isempty(files), uialert(fig,'No TIFFs in ch1/avgz.','Load error'); return; end
        S.files=files; S.nTiles=numel(files);

        sx_=nan(S.nTiles,1); sy_=nan(S.nTiles,1);
        ci_=nan(S.nTiles,1); ri_=nan(S.nTiles,1);
        for t=1:S.nTiles, [sx_(t),sy_(t),ci_(t),ri_(t)]=parse_name(files(t)); end

        haveCR = ~any(isnan(ci_)) && ~any(isnan(ri_));   % explicit colXX/rowYY
        haveXY = ~any(isnan(sx_)) && ~any(isnan(sy_));   % stage coords (microns)
        if haveCR
            uC=unique(ci_); uR=unique(ri_);
            S.colIdx=arrayfun(@(v) find(uC==v,1),ci_);
            S.rowIdx=arrayfun(@(v) find(uR==v,1),ri_);
            S.nCols=numel(uC); S.nRows=numel(uR);
        elseif haveXY
            % Derive grid from unique stage coords: ascending X -> col, Y -> row.
            [uX,~,S.colIdx]=unique(sx_); [uY,~,S.rowIdx]=unique(sy_);
            S.colIdx=S.colIdx(:); S.rowIdx=S.rowIdx(:);
            S.nCols=numel(uX); S.nRows=numel(uY);
            % Seed step from the REAL micron spacing (user can still tune).
            if numel(uX)>1, S.stepXpx=round(median(diff(uX))/umPerPx); end
            if numel(uY)>1, S.stepYpx=round(median(diff(uY))/umPerPx); end
        else
            uialert(fig,sprintf(['Filenames lack both colXX/rowYY indices and ' ...
                'tileNN_<X>_<Y> stage coordinates.\nExample: %s'], ...
                files(1)),'Parse error'); return;
        end

        info=imfinfo(files(1)); S.H0=info(1).Height; S.W0=info(1).Width;
        I0=imread(files(1),1); S.inClass=class(I0);

        S.manX=zeros(S.nTiles,1); S.manY=zeros(S.nTiles,1);
        S.coordDir=fullfile(char(S.datasetRoot),'matlab_stitch');
        if ~isfolder(S.coordDir), mkdir(S.coordDir); end
        lblPath.Text=char(S.datasetRoot);

        % Reload previous session?
        prevMat=fullfile(S.coordDir,'stitch_manual_coords.mat');
        if isfile(prevMat)
            sel=uiconfirm(fig,'Reload previous step/nudges from stitch_manual_coords.mat?', ...
                'Reload','Options',{'Reload','Start fresh'},'DefaultOption',1);
            if strcmp(sel,'Reload')
                p=load(prevMat);
                if isfield(p,'stepXpx'), S.stepXpx=p.stepXpx; end
                if isfield(p,'stepYpx'), S.stepYpx=p.stepYpx; end
                if isfield(p,'flipX'), S.flipX=p.flipX; cbFlipX.Value=S.flipX; end
                if isfield(p,'flipY'), S.flipY=p.flipY; cbFlipY.Value=S.flipY; end
                if isfield(p,'manX')&&numel(p.manX)==S.nTiles, S.manX=p.manX(:); S.manY=p.manY(:); end
            end
        end

        buildProjections();
        S.loaded=true;
        recomputePositions(); blendDisplay(); showMosaic(); autoBC();
        setStat(sprintf('Loaded %d tiles (%dx%d grid).',S.nTiles,S.nRows,S.nCols));
    end

    function buildProjections()
        ds=S.dsPreview;
        S.tileImg=cell(S.nTiles,1); S.tileImgDS=cell(S.nTiles,1);
        wb=uiprogressdlg(fig,'Title','Building avg projections','Message','...');
        for t=1:S.nTiles
            wb.Value=t/S.nTiles;
            tf=Tiff(S.files(t),'r');
            nZ=1; while ~tf.lastDirectory(), tf.nextDirectory(); nZ=nZ+1; end
            zu=S.zUseProj; if isempty(zu), zu=1:nZ; end
            zu=zu(zu>=1 & zu<=nZ);
            acc=zeros(S.H0,S.W0,'double');
            for k=zu(:)', tf.setDirectory(k); acc=acc+double(tf.read()); end
            tf.close();
            img=acc/numel(zu);
            S.tileImg{t}=img;
            S.tileImgDS{t}=single(imresize(img,1/ds,'bilinear'));
        end
        close(wb);
        S.Hd=size(S.tileImgDS{1},1); S.Wd=size(S.tileImgDS{1},2);
    end

    function recomputePositions()
        bx=(S.colIdx-1)*S.stepXpx;
        by=(S.rowIdx-1)*S.stepYpx;
        if S.flipX, bx=max(bx)-bx; end
        if S.flipY, by=max(by)-by; end
        S.x=bx+S.manX; S.y=by+S.manY;
        S.x=S.x-min(S.x); S.y=S.y-min(S.y);
        S.outW=ceil(max(S.x)+S.W0); S.outH=ceil(max(S.y)+S.H0);
    end

    function blendDisplay()
        % Fast downsampled preview blend (integer placement, single precision).
        ds=S.dsPreview; Hd=S.Hd; Wd=S.Wd;
        wTd=single(linear_blend_weights(Hd,Wd));
        outHd=ceil(S.outH/ds)+2; outWd=ceil(S.outW/ds)+2;
        acc=zeros(outHd,outWd,'single'); wac=zeros(outHd,outWd,'single');
        S.dsR0=zeros(S.nTiles,1); S.dsC0=zeros(S.nTiles,1);
        for t=1:S.nTiles
            r0=round(S.y(t)/ds); c0=round(S.x(t)/ds);
            r0=max(0,min(r0,outHd-Hd)); c0=max(0,min(c0,outWd-Wd));
            rr=r0+1:r0+Hd; cc=c0+1:c0+Wd;
            acc(rr,cc)=acc(rr,cc)+S.tileImgDS{t}.*wTd;
            wac(rr,cc)=wac(rr,cc)+wTd;
            S.dsR0(t)=r0; S.dsC0(t)=c0;
        end
        S.mosaic=acc./max(wac,eps('single'));
    end

    function cb_flip(~,~)
        if ~S.loaded, return; end
        S.flipX=cbFlipX.Value; S.flipY=cbFlipY.Value;
        recomputePositions(); blendDisplay(); showMosaic();
    end

    function assignMove(v), S.moveMode=v; end
    function assignMoveStep(v), S.moveStep=max(1,round(v)); end

    function cb_key(~,evt)
        if ~S.loaded || S.selectedTile==0, return; end
        step=S.moveStep; if any(strcmp(evt.Modifier,'shift')), step=S.moveStep*10; end
        dx=0; dy=0;
        switch evt.Key
            case 'leftarrow',  dx=-step; case 'rightarrow', dx=step;
            case 'uparrow',    dy=-step; case 'downarrow',  dy=step;
            otherwise, return;
        end
        g=selectionGroup();
        S.manX(g)=S.manX(g)+dx; S.manY(g)=S.manY(g)+dy;
        recomputePositions(); blendDisplay(); showMosaic(); updateSel();
    end

    function g=selectionGroup()
        t=S.selectedTile;
        switch S.moveMode
            case 'Row', g=find(S.rowIdx==S.rowIdx(t));
            case 'Col', g=find(S.colIdx==S.colIdx(t));
            otherwise,  g=t;
        end
    end

    function cb_resetSel(~,~)
        if ~S.loaded||S.selectedTile==0, return; end
        g=selectionGroup(); S.manX(g)=0; S.manY(g)=0;
        recomputePositions(); blendDisplay(); showMosaic(); updateSel();
    end
    function cb_resetAll(~,~)
        if ~S.loaded, return; end
        S.manX(:)=0; S.manY(:)=0;
        recomputePositions(); blendDisplay(); showMosaic(); updateSel();
        setStat('All nudges reset.');
    end

    function cb_imgClick(~,evt)
        if ~S.loaded, return; end
        cp=evt.IntersectionPoint; cx=cp(1); cy=cp(2);   % DS coords
        best=0; bd=Inf;
        for t=1:S.nTiles
            r0=S.dsR0(t); c0=S.dsC0(t);
            if cy>=r0+1 && cy<=r0+S.Hd && cx>=c0+1 && cx<=c0+S.Wd
                d=(cx-(c0+S.Wd/2))^2+(cy-(r0+S.Hd/2))^2;
                if d<bd, bd=d; best=t; end
            end
        end
        S.selectedTile=best; updateSel(); drawOverlays();
    end

    function autoBC()
        if ~S.loaded, return; end
        v=S.mosaic(S.mosaic>0); if isempty(v), v=S.mosaic(:); end
        lo=prctile(v,1); hi=prctile(v,99.5); if hi<=lo, hi=lo+1; end
        ax.CLim=double([lo hi]);
    end

    function cb_save(~,~)
        if ~S.loaded, return; end
        % --- 1) full-res AVG-PROJ mosaic (2D) ---
        setStat('Rendering full-res avg-proj mosaic...'); drawnow;
        M=blendFull();
        outProj=fullfile(S.coordDir,sprintf('stitched_ch%d_avgproj.tif',S.refChID));
        write_single_tiff(outProj,cast(round(M),S.inClass),S.inClass);
        % --- 2) coords ---
        x=S.x; y=S.y; manX=S.manX; manY=S.manY;
        stepXpx=S.stepXpx; stepYpx=S.stepYpx; flipX=S.flipX; flipY=S.flipY;
        nRows=S.nRows; nCols=S.nCols; colIdx=S.colIdx; rowIdx=S.rowIdx;
        H0=S.H0; W0=S.W0; outH=S.outH; outW=S.outW; files=S.files;
        save(fullfile(S.coordDir,'stitch_manual_coords.mat'), ...
            'x','y','manX','manY','stepXpx','stepYpx','flipX','flipY', ...
            'nRows','nCols','colIdx','rowIdx','H0','W0','outH','outW','files');
        % --- 3) full-res VOLUME mosaic (Z preserved, per-slice blend) ---
        outVol=fullfile(S.coordDir,sprintf('stitched_ch%d_volume.tif',S.refChID));
        blendVolumeWrite(outVol);
        setStat(sprintf('Saved avg-proj + volume + coords to %s',S.coordDir));
    end

    function blendVolumeWrite(outFile)
        % Per-Z-slice blend of the (multi-page) avgz tiles at the saved x,y.
        % Same placement & linear-blend weights as the avg-proj, but Z is kept.
        nZ=count_pages(S.files(1));
        wT=linear_blend_weights(S.H0,S.W0);
        % precompute per-tile windows
        rr=cell(S.nTiles,1); cc=cell(S.nTiles,1);
        for t=1:S.nTiles
            r0=round(S.y(t)); c0=round(S.x(t));
            r0=max(0,min(r0,S.outH-S.H0)); c0=max(0,min(c0,S.outW-S.W0));
            rr{t}=r0+1:r0+S.H0; cc{t}=c0+1:c0+S.W0;
        end
        % open all tiles
        T=cell(S.nTiles,1);
        for t=1:S.nTiles, T{t}=Tiff(char(S.files(t)),'r'); end
        cleanT=onCleanup(@() cellfun(@safeClose,T)); %#ok<NASGU>
        % output tiff (BigTIFF if large)
        [bps,sf]=class_to_tiff_format(S.inClass);
        estBytes=double(S.outH)*double(S.outW)*double(nZ)*double(bps/8);
        if isfile(outFile), delete(outFile); end
        if estBytes>3.5e9, tout=Tiff(outFile,'w8'); else, tout=Tiff(outFile,'w'); end
        cleanO=onCleanup(@() safeClose(tout)); %#ok<NASGU>
        tag.ImageLength=S.outH; tag.ImageWidth=S.outW;
        tag.Photometric=Tiff.Photometric.MinIsBlack; tag.SamplesPerPixel=1;
        tag.BitsPerSample=bps; tag.SampleFormat=sf;
        tag.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
        tag.Compression=Tiff.Compression.None; tag.RowsPerStrip=64;
        tag.Software=sprintf('stitch_manual_gui volume ch%d',S.refChID);
        wb=uiprogressdlg(fig,'Title','Rendering volume','Message','...');
        cleanWb=onCleanup(@() closeIfValid(wb)); %#ok<NASGU>
        for iz=1:nZ
            wb.Value=iz/nZ; wb.Message=sprintf('Z slice %d/%d',iz,nZ);
            acc=zeros(S.outH,S.outW,'double'); wac=zeros(S.outH,S.outW,'double');
            for t=1:S.nTiles
                T{t}.setDirectory(iz);
                img=double(T{t}.read());
                acc(rr{t},cc{t})=acc(rr{t},cc{t})+img.*wT;
                wac(rr{t},cc{t})=wac(rr{t},cc{t})+wT;
            end
            sliceImg=cast(round(acc./max(wac,eps)),S.inClass);
            tout.setTag(tag); tout.write(sliceImg);
            if iz<nZ, tout.writeDirectory(); end
        end
    end

    function M=blendFull()
        % Full-resolution blend (integer placement). Used only on save.
        wT=linear_blend_weights(S.H0,S.W0);
        acc=zeros(S.outH,S.outW,'double'); wac=zeros(S.outH,S.outW,'double');
        for t=1:S.nTiles
            r0=round(S.y(t)); c0=round(S.x(t));
            r0=max(0,min(r0,S.outH-S.H0)); c0=max(0,min(c0,S.outW-S.W0));
            rr=r0+1:r0+S.H0; cc=c0+1:c0+S.W0;
            acc(rr,cc)=acc(rr,cc)+S.tileImg{t}.*wT;
            wac(rr,cc)=wac(rr,cc)+wT;
        end
        M=acc./max(wac,eps);
    end

%% ========================= DISPLAY =========================
    function showMosaic()
        if ~S.loaded, return; end
        if isempty(hImg)||~isvalid(hImg)||any(size(hImg.CData)~=size(S.mosaic))
            cla(ax); hImg=imagesc(ax,S.mosaic); colormap(ax,gray(256));
            ax.DataAspectRatio=[1 1 1]; ax.XTick=[]; ax.YTick=[];
            hImg.ButtonDownFcn=@cb_imgClick; ax.ButtonDownFcn=@cb_imgClick;
            hBord=[]; hSel=[];
        else
            hImg.CData=S.mosaic;
        end
        title(ax,sprintf('ch%d | %d tiles | step(%d,%d)px | %dx%d (preview /%d)', ...
            S.refChID,S.nTiles,S.stepXpx,S.stepYpx,S.outH,S.outW,S.dsPreview),'Interpreter','none');
        drawOverlays();
    end

    function drawOverlays()
        if ~S.loaded, return; end
        % All tile borders as ONE line object (NaN-separated) for speed.
        if ~isempty(hBord)&&isvalid(hBord), delete(hBord); end; hBord=[];
        if cbSeams.Value
            X=[]; Y=[];
            for t=1:S.nTiles
                r1=S.dsR0(t)+.5; r2=S.dsR0(t)+S.Hd+.5; c1=S.dsC0(t)+.5; c2=S.dsC0(t)+S.Wd+.5;
                X=[X c1 c2 c2 c1 c1 NaN]; Y=[Y r1 r1 r2 r2 r1 NaN]; %#ok<AGROW>
            end
            hold(ax,'on');
            hBord=line(ax,X,Y,'Color',[.5 .5 .5],'LineWidth',.4,'HitTest','off');
            hold(ax,'off');
        end
        if ~isempty(hSel)&&isvalid(hSel), delete(hSel); end; hSel=[];
        if S.selectedTile>0
            t=S.selectedTile;
            r1=S.dsR0(t)+.5; r2=S.dsR0(t)+S.Hd+.5; c1=S.dsC0(t)+.5; c2=S.dsC0(t)+S.Wd+.5;
            hold(ax,'on');
            hSel=line(ax,[c1 c2 c2 c1 c1],[r1 r1 r2 r2 r1],'Color','c','LineWidth',2.5,'HitTest','off');
            hold(ax,'off');
        end
    end

    function updateSel()
        if S.selectedTile>0
            t=S.selectedTile;
            lblSel.Text=sprintf('Tile %d  grid(r%d,c%d)  nudge=(%+d,%+d)px  [%s]', ...
                t,S.rowIdx(t),S.colIdx(t),round(S.manX(t)),round(S.manY(t)),S.moveMode);
        else
            lblSel.Text='Selected: none (click a tile)';
        end
    end

    function setStat(m), lblStat.Text=m; end
end

%% ========================= LOCAL HELPERS =========================
function files=list_tifs(folder)
c=[dir(fullfile(folder,'*.tif')); dir(fullfile(folder,'*.tiff'))];
files=strings(numel(c),1);
for k=1:numel(c), files(k)=string(fullfile(c(k).folder,c(k).name)); end
end

function files=sort_by_colrow(files)
n=numel(files); key=nan(n,2);
for k=1:n
    [sx,sy,ci,ri]=parse_name(files(k));
    if ~isnan(ci)&&~isnan(ri), key(k,:)=[ri ci];
    elseif ~isnan(sy)&&~isnan(sx), key(k,:)=[sy sx]; end
end
if all(~isnan(key(:))), [~,ord]=sortrows(key); files=files(ord); else, files=sort(files); end
end

function [sx,sy,col,row]=parse_name(fullpath)
[~,bn]=fileparts(fullpath); bn=char(bn);
sx=num1(regexp(bn,'_x(-?\d+)','tokens','once'));
sy=num1(regexp(bn,'_y(-?\d+)','tokens','once'));
% Fallback: ScanImage "tileNN_<X>_<Y>_..." stage coords (microns).
if isnan(sx)||isnan(sy)
    tok=regexp(bn,'^tile\d+_(-?\d+)_(-?\d+)','tokens','once');
    if ~isempty(tok), sx=str2double(tok{1}); sy=str2double(tok{2}); end
end
col=num1(regexp(bn,'col(\d+)','tokens','once'));
row=num1(regexp(bn,'row(\d+)','tokens','once'));
end
function v=num1(tok), if isempty(tok), v=NaN; else, v=str2double(tok{1}); end, end

function w=linear_blend_weights(H,W)
[xg,yg]=meshgrid(1:W,1:H);
w=double(min(min(xg-1,W-xg),min(yg-1,H-yg))+1);
end

function write_single_tiff(outFile,img,cls)
[bps,sf]=class_to_tiff_format(cls);
if isfile(outFile), delete(outFile); end
t=Tiff(outFile,'w'); cln=onCleanup(@() safeClose(t)); %#ok<NASGU>
tag.ImageLength=size(img,1); tag.ImageWidth=size(img,2);
tag.Photometric=Tiff.Photometric.MinIsBlack; tag.SamplesPerPixel=1;
tag.BitsPerSample=bps; tag.SampleFormat=sf;
tag.PlanarConfiguration=Tiff.PlanarConfiguration.Chunky;
tag.Compression=Tiff.Compression.None; tag.RowsPerStrip=64;
tag.Software='stitch_manual_gui avgproj ch1';
t.setTag(tag); t.write(img);
end

function [bps,sf]=class_to_tiff_format(cls)
switch char(cls)
    case 'uint8',  bps=8;  sf=Tiff.SampleFormat.UInt;
    case 'uint16', bps=16; sf=Tiff.SampleFormat.UInt;
    case 'int16',  bps=16; sf=Tiff.SampleFormat.Int;
    case 'single', bps=32; sf=Tiff.SampleFormat.IEEEFP;
    case 'double', bps=64; sf=Tiff.SampleFormat.IEEEFP;
    otherwise, error('Unsupported class: %s',cls);
end
end
function safeClose(t), try t.close(); catch, end, end

function n=count_pages(file)
tf=Tiff(char(file),'r'); cln=onCleanup(@() safeClose(tf)); %#ok<NASGU>
n=1; while ~tf.lastDirectory(), tf.nextDirectory(); n=n+1; end
end

function closeIfValid(d), try if isvalid(d), close(d); end, catch, end, end
