function repr_video_gui_RZ()
% repr_video_gui_RZ  Synced GUI: Basler breathing-camera video (chest motion)
%   + 3 calcium ROI patches on the left; scrolling breath + dF/F traces on
%   the right. Play / scrub / export a composite clip.
%
% Left column (top→bottom):
%   [1] Breathing camera video  (Basler chest-motion capture;
%                                use a transcoded/labeled MP4 — FFV1 AVIs
%                                are not directly readable by VideoReader)
%   [2] Ca ROI 1 patch          (yellow perimeter, linear contrast)
%   [3] Ca ROI 2 patch
%   [4] Ca ROI 3 patch
%
% Right column:
%   Scrolling trace panel — Breath (cyan) + dF/F × 3 (green / blue / magenta)
%
% Sync: t0Offset_s sets breathing-cam lead relative to calcium (s).
%       beh_frame = round((ca_frame-1) × fps_beh/fps_ca + t0Offset_s×fps_beh)+1
%
% Export clip: renders the full composite to MP4.

%% ============================== USER PARAMS ==============================
breathCamDir = "D:\260615_thermalNbasler\basler_breath\videos";
%   ^ folder with Basler breathing-camera video (transcoded MP4 or labeled MP4)
%     The raw Basler AVI is FFV1 (unreadable by VideoReader). Transcode first:
%       ffmpeg -i input.avi -c:v libx264 -crf 18 output_breath.mp4
%     or use the DLC-labeled MP4 if DLC tracking was run on this video.
fovPath     = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\roi5_1400-1200-20_20lp_930_x4.4_512x512_6000f_00001";
%   ^ calcium FOV folder (must contain *_cpSAM_output.mat and *_MC*.tif)
roiDefaults = [1 2 3];      % initial 3 ROI indices to display
patchHalf_um = 20;          % Ca patch half-size (um)
t0Offset_s   = 0;           % behavior lead relative to calcium (s)
%   beh_t = ca_t + t0Offset_s  (positive = beh video starts later)
T0_init      = 0;           % initial playhead time (s)
WINSEC       = 3;           % trace window shown around playhead (s)
EXPORT_FPS   = 30;          % export target fps (speed box scales this)
EXPORT_WIN   = 3;           % rolling trace window during export (s)
EXPORT_DIR   = "";          % "" = same folder as this script
%% =========================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

%% ---- state ----
vBeh = []; nBeh = 0; fps_beh = 30;
fps_ca = 30; n_ca = 0; cf = 1;
ca_patches = {[], [], []};   % {[pH x pW x T_ca] single} raw crop per ROI
ca_lo = [0 0 0]; ca_hi = [1 1 1];
ca_perim_up = {[], [], []};  % [dispH x dispW] logical, upscaled perimeter
dFF_traces = [];             % [T_ca x 3]
breath_trace = [];           % [T_ca x 1]
t_ca = [];                   % [T_ca x 1] seconds
loading = false;
tmr = [];
roiIdxs = roiDefaults(:)';
DISP_CA = 150;   % Ca patch display size (px, square)
DISP_BHEH_H = 220; DISP_BHEH_W = 280;  % behavior video display size
caClrs = {[0.15 0.70 0.20], [0.10 0.45 0.90], [0.85 0.20 0.70]};  % g/b/m

%% ---- figure ----
fig = figure('Color','w','Position',[50 30 1540 920], ...
    'Name','Repr Video GUI — Ca + Behavior', 'NumberTitle','off', ...
    'CloseRequestFcn', @onClose);

% --- left column axes ---
% Behavior video: top
axBeh = axes(fig,'Units','normalized','Position',[0.016 0.525 0.255 0.445]);
axis(axBeh,'off'); box(axBeh,'off');

% Three Ca patches stacked below
caY   = [0.350  0.185  0.020];   % bottom edge of each Ca axes
caH   = 0.150;
axCa  = gobjects(3,1);
for k = 1:3
    axCa(k) = axes(fig,'Units','normalized','Position',[0.016, caY(k), 0.255, caH]);
    axis(axCa(k),'off'); box(axCa(k),'off');
end

% --- right: trace panel ---
axT = axes(fig,'Units','normalized','Position',[0.310 0.130 0.672 0.835]);
hold(axT,'on'); grid(axT,'on');

%% ---- controls ----
% Row A: FOV folder edit + ROI boxes + Load
uicontrol(fig,'Style','text','Units','normalized','Position',[0.014 0.096 0.035 0.024],...
    'String','FOV:','BackgroundColor','w','HorizontalAlignment','left','FontSize',9);
eFOV = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.052 0.093 0.175 0.033],...
    'String',char(fovPath),'HorizontalAlignment','left','FontSize',8);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.234 0.096 0.030 0.024],...
    'String','ROIs:','BackgroundColor','w','HorizontalAlignment','left','FontSize',9);
eROI  = gobjects(3,1);
roiX  = [0.266 0.310 0.354];
for k = 1:3
    eROI(k) = uicontrol(fig,'Style','edit','Units','normalized','Position',[roiX(k) 0.093 0.040 0.033],...
        'String',num2str(roiIdxs(k)),'FontSize',9);
end
uicontrol(fig,'Style','text','Units','normalized','Position',[0.397 0.096 0.042 0.024],...
    'String','t-offs:','BackgroundColor','w','HorizontalAlignment','left','FontSize',9);
eOff = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.438 0.093 0.044 0.033],...
    'String',num2str(t0Offset_s),'FontSize',9);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.488 0.091 0.070 0.038],...
    'String','Load','FontWeight','bold','FontSize',10,'Callback',@onLoad);

% Row B: slider + play + speed
sld = uicontrol(fig,'Style','slider','Units','normalized','Position',[0.014 0.058 0.530 0.025],...
    'Min',1,'Max',2,'Value',1,'Callback',@onSlide);
btnPlay = uicontrol(fig,'Style','togglebutton','Units','normalized','Position',[0.555 0.050 0.063 0.037],...
    'String','Play','FontWeight','bold','FontSize',10,'Callback',@onPlay);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.627 0.053 0.040 0.024],...
    'String','speed','BackgroundColor','w','HorizontalAlignment','left','FontSize',9);
eSpd = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.668 0.050 0.042 0.037],'String','1','FontSize',9);

% Row C: export controls
uicontrol(fig,'Style','text','Units','normalized','Position',[0.014 0.017 0.056 0.024],...
    'String','start (s)','BackgroundColor','w','HorizontalAlignment','left','FontSize',9);
eStart = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.070 0.014 0.050 0.033],...
    'String',num2str(T0_init),'FontSize',9);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.127 0.017 0.040 0.024],...
    'String','dur (s)','BackgroundColor','w','HorizontalAlignment','left','FontSize',9);
eDur = uicontrol(fig,'Style','edit','Units','normalized','Position',[0.167 0.014 0.050 0.033],...
    'String','3','FontSize',9);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.224 0.012 0.105 0.040],...
    'String','Export clip','FontWeight','bold','FontSize',10,'Callback',@onExport);
txStat = uicontrol(fig,'Style','text','Units','normalized','Position',[0.340 0.015 0.640 0.028],...
    'String','Set FOV path and click Load.','BackgroundColor','w',...
    'HorizontalAlignment','left','FontName','Consolas','FontSize',8);

% image / line handles
him_beh = [];
him_ca  = {[], [], []};
hPlay   = [];
hTr     = gobjects(4,1);

%% ---- auto-load if default path exists ----
if isfolder(char(fovPath))
    onLoad([],[]);
end

%% ============================================================
%  LOAD
%% ============================================================
    function onLoad(~,~)
        if loading, return; end
        loading = true; cl = onCleanup(@() setLd(false)); %#ok<NASGU>
        stopplay(); set(btnPlay,'Value',0,'String','Play');

        fovDir = strtrim(get(eFOV,'String'));
        if ~isfolder(fovDir), stat('FOV folder not found.'); return; end
        for k=1:3
            v = str2double(get(eROI(k),'String'));
            if isnan(v), stat('Bad ROI index.'); return; end
            roiIdxs(k) = round(v);
        end
        t0Offset_s = str2double(get(eOff,'String'));
        if isnan(t0Offset_s), t0Offset_s = 0; end

        stat('Detecting FPS...');

        %% FPS + pixel size
        try
            [fps_ca, sc_meta] = detect_session_fps(fovDir);
        catch
            fps_ca = 30; sc_meta = struct('pixelSize_um', 1.7778);
        end
        if ~isfield(sc_meta,'pixelSize_um') || ~isfinite(sc_meta.pixelSize_um)
            sc_meta.pixelSize_um = 1.7778;
        end
        pixSz = sc_meta.pixelSize_um;
        patchHalf_px = round(patchHalf_um / pixSz);

        %% cpSAM
        stat('Loading cpSAM...');
        sam_h = dir(fullfile(fovDir,'*_cpSAM_output.mat'));
        if isempty(sam_h), stat('No cpSAM output found.'); return; end
        sd = load(fullfile(sam_h(1).folder,sam_h(1).name),'maskL','F');
        if ~isfield(sd,'maskL')||~isfield(sd,'F'), stat('cpSAM missing maskL/F.'); return; end
        maskL_loc = sd.maskL;
        F_all = double(sd.F);
        [T_F, N_roi] = size(F_all);

        if any(roiIdxs<1) || any(roiIdxs>N_roi)
            stat(sprintf('ROI out of range [1..%d].',N_roi)); return;
        end

        %% dF/F
        stat('Computing dF/F...');
        dout = helper.dFF_RZ(F_all, 'FPS', fps_ca);
        dFF_all = dout.dFF;     % [T_dff x N]
        T_dff   = size(dFF_all,1);
        dFF_traces = dFF_all(:, roiIdxs);  % [T_dff x 3]

        %% Breath trace (from breath_master.mat or zeros)
        breath_trace = zeros(T_dff,1);
        bm_h = dir(fullfile(fovDir,'breath_master.mat'));
        if ~isempty(bm_h)
            try
                bm = load(fullfile(bm_h(1).folder,bm_h(1).name),'breath');
                if isfield(bm,'breath')
                    br = bm.breath(:);
                    if numel(br) >= T_F
                        nDrop = T_F - T_dff;
                        br = br(nDrop+1:T_F);
                    end
                    if numel(br) ~= T_dff
                        br = interp1(linspace(0,1,numel(br)),br,linspace(0,1,T_dff),'pchip');
                    end
                    breath_trace = br(:);
                end
            catch
            end
        end

        n_ca = T_dff;
        t_ca = (0:n_ca-1)' / fps_ca;

        %% Load Ca patches from MC TIFF
        stat('Loading TIFF patches...');
        try
            tifP = local_discover_mc_tiff(fovDir);
        catch ME
            stat(['TIFF: ' ME.message]); return;
        end
        tifInfo  = imfinfo(tifP);
        imgH = tifInfo(1).Height;  imgW = tifInfo(1).Width;
        T_tif = numel(tifInfo);
        if T_tif ~= T_F
            stat(sprintf('TIFF/F mismatch: %d vs %d.',T_tif,T_F)); return;
        end
        nDrop = T_F - T_dff;

        props = regionprops(maskL_loc, 'Centroid');
        if numel(props) < max(roiIdxs)
            stat('maskL has fewer ROIs than requested.'); return;
        end

        % Pre-compute crop bounds for each ROI
        r1v = zeros(3,1); r2v = zeros(3,1); c1v = zeros(3,1); c2v = zeros(3,1);
        for k = 1:3
            ri = roiIdxs(k);
            cx = round(props(ri).Centroid(1));
            cy = round(props(ri).Centroid(2));
            r1v(k) = max(1, cy-patchHalf_px);
            r2v(k) = min(imgH, cy+patchHalf_px);
            c1v(k) = max(1, cx-patchHalf_px);
            c2v(k) = min(imgW, cx+patchHalf_px);
        end

        % Allocate patch buffers
        for k = 1:3
            ph = r2v(k)-r1v(k)+1; pw = c2v(k)-c1v(k)+1;
            ca_patches{k} = zeros(ph, pw, n_ca, 'single');
        end

        % Read frames (load only rows needed — one imread per frame)
        for ff = 1:n_ca
            tidx = nDrop + ff;
            frm  = single(imread(tifP, tidx, 'Info', tifInfo));
            for k = 1:3
                ca_patches{k}(:,:,ff) = frm(r1v(k):r2v(k), c1v(k):c2v(k));
            end
            if mod(ff,500)==0
                stat(sprintf('Loading frames %d/%d...', ff, n_ca));
                drawnow limitrate;
            end
        end

        % Contrast limits + upscaled perimeter per ROI
        for k = 1:3
            ri = roiIdxs(k);
            msk = (maskL_loc(r1v(k):r2v(k), c1v(k):c2v(k)) == ri);
            roi_px = ca_patches{k}(repmat(msk,[1 1 n_ca]));
            lo_k = double(prctile(roi_px, 2));
            hi_k = double(prctile(roi_px, 99.9));
            if hi_k <= lo_k, hi_k = lo_k + 1; end
            ca_lo(k) = lo_k; ca_hi(k) = hi_k;
            pm = bwperim(msk,4);
            ca_perim_up{k} = imresize(double(pm), [DISP_CA DISP_CA], 'nearest') > 0.5;
        end

        %% Behavior video
        vBeh = []; nBeh = 0; fps_beh = 30;
        bv_hits = [dir(fullfile(char(breathCamDir),'*_labeled*.mp4')); ...
                   dir(fullfile(char(breathCamDir),'*.mp4'))];
        if ~isempty(bv_hits)
            try
                vBeh = VideoReader(fullfile(bv_hits(1).folder, bv_hits(1).name)); %#ok<TNMLP>
                nBeh = vBeh.NumFrames;
                fps_beh = vBeh.FrameRate;
                stat(sprintf('Behavior video: %d frames @ %.0f fps', nBeh, fps_beh));
            catch
                vBeh = []; stat('Could not open breathing camera video (will show blank).');
            end
        else
            stat('No MP4 found in breathCamDir — showing blank. Transcode the Basler AVI first.');
        end

        %% Trace panel
        cla(axT); hold(axT,'on'); grid(axT,'on');
        trClrs = {[0.30 0.80 0.95], caClrs{1}, caClrs{2}, caClrs{3}};
        trNames = {'Breath','ROI 1','ROI 2','ROI 3'};
        trOffsets = [3 2 1 0];
        trData = [breath_trace, dFF_traces];  % [n_ca x 4]
        for kk = 1:4
            sig = trData(:,kk);
            rng = max(sig)-min(sig); if rng<eps, rng=1; end
            sig_n = (sig-min(sig))/rng + trOffsets(kk);
            hTr(kk) = plot(axT, t_ca, sig_n, '-', ...
                'Color', trClrs{kk}, 'LineWidth', 1.0, 'DisplayName', trNames{kk});
        end
        ylim(axT,[-0.25 4.3]);
        yticks(axT, fliplr(trOffsets + 0.5));
        yticklabels(axT, fliplr(trNames));
        xlabel(axT,'Time (s)'); ylabel(axT,'');
        legend(axT, trNames, 'Location','northeast','Orientation','vertical','FontSize',8);
        hPlay = line(axT, [0 0], [-0.25 4.3], ...
            'Color',[0.9 0.1 0.1],'LineWidth',1.8,'HandleVisibility','off');

        %% Seed image handles at initial frame
        cf = max(1, min(n_ca, round(T0_init * fps_ca) + 1));

        % Behavior image
        cla(axBeh); axis(axBeh,'off');
        img_beh = get_beh_frame(cf);
        him_beh = imshow(img_beh, 'Parent', axBeh);
        title(axBeh, 'Breathing cam', 'FontSize', 9, 'Color', 'k');

        % Ca images
        for k = 1:3
            cla(axCa(k)); axis(axCa(k),'off');
            rgb = ca2rgb(ca_patches{k}(:,:,cf), ca_lo(k), ca_hi(k), ca_perim_up{k});
            him_ca{k} = imshow(rgb, 'Parent', axCa(k));
            title(axCa(k), sprintf('ROI %d', roiIdxs(k)), ...
                'FontSize', 9, 'Color', caClrs{k});
        end

        set(sld,'Min',1,'Max',n_ca,'Value',cf,'SliderStep',[1/max(n_ca-1,1), 30/max(n_ca-1,1)]);
        overlay();
        stat(sprintf('Loaded: %d Ca frames @ %d Hz | %d ROIs total', n_ca, fps_ca, N_roi));
    end

    function setLd(v), loading = v; end

%% ============================================================
%  TRANSPORT
%% ============================================================
    function onSlide(~,~)
        if n_ca==0, return; end
        stopplay(); cf = max(1,min(n_ca, round(get(sld,'Value')))); redraw();
    end

    function onPlay(src,~)
        if n_ca==0, set(src,'Value',0); return; end
        if get(src,'Value')
            set(src,'String','Pause');
            PERIOD = round(1000/fps_ca)/1000;
            PERIOD = max(PERIOD, round(1000/60)/1000);
            tmr = timer('ExecutionMode','fixedRate','Period',PERIOD,...
                'BusyMode','drop','TimerFcn',@tick);
            start(tmr);
        else
            stopplay();
        end
    end

    function tick(~,~)
        if loading || n_ca==0, return; end
        s = str2double(get(eSpd,'String')); if isnan(s)||s<=0, s=1; end
        cf = cf + max(1, round(s));
        if cf >= n_ca
            cf = n_ca; stopplay(); set(btnPlay,'Value',0,'String','Play');
        end
        redraw_fast();
    end

    function stopplay()
        if ~isempty(tmr) && isvalid(tmr), stop(tmr); delete(tmr); end
        tmr = [];
    end

    function redraw()
        if n_ca==0, return; end
        cf = max(1,min(n_ca,cf));
        if ~isempty(him_beh) && isgraphics(him_beh)
            set(him_beh,'CData', get_beh_frame(cf));
        end
        for k=1:3
            if ~isempty(him_ca{k}) && isgraphics(him_ca{k})
                set(him_ca{k},'CData', ca2rgb(ca_patches{k}(:,:,cf), ca_lo(k), ca_hi(k), ca_perim_up{k}));
            end
        end
        overlay(); drawnow limitrate;
    end

    function redraw_fast()
        if n_ca==0, return; end
        cf = max(1,min(n_ca,cf));
        if ~isempty(him_beh) && isgraphics(him_beh)
            set(him_beh,'CData', get_beh_frame(cf));
        end
        for k=1:3
            if ~isempty(him_ca{k}) && isgraphics(him_ca{k})
                set(him_ca{k},'CData', ca2rgb(ca_patches{k}(:,:,cf), ca_lo(k), ca_hi(k), ca_perim_up{k}));
            end
        end
        overlay(); drawnow;
    end

    function overlay()
        if isempty(t_ca) || ~isgraphics(hPlay), return; end
        tc = t_ca(cf);
        set(hPlay,'XData',[tc tc]);
        xlim(axT, tc + WINSEC*[-0.5 0.5]);
        set(sld,'Value',cf);
    end

    function img = get_beh_frame(ca_fr)
        if isempty(vBeh)
            img = zeros(DISP_BHEH_H, DISP_BHEH_W, 3, 'uint8');
            return;
        end
        bfr = max(1, min(nBeh, round((ca_fr-1)*fps_beh/fps_ca + t0Offset_s*fps_beh) + 1));
        try
            raw = read(vBeh, bfr);
            if size(raw,3)==1, raw = repmat(raw,[1 1 3]); end
            img = imresize(raw, [DISP_BHEH_H DISP_BHEH_W]);
        catch
            img = zeros(DISP_BHEH_H, DISP_BHEH_W, 3, 'uint8');
        end
    end

%% ============================================================
%  EXPORT
%% ============================================================
    function onExport(~,~)
        if n_ca==0, stat('Load data first.'); return; end
        t0e  = str2double(get(eStart,'String'));
        dure = str2double(get(eDur,'String'));
        if isnan(t0e)||isnan(dure)||dure<=0, stat('Invalid start/dur.'); return; end
        spd  = str2double(get(eSpd,'String')); if isnan(spd)||spd<=0, spd=1; end

        f0 = max(1, round(t0e*fps_ca)+1);
        f1 = min(n_ca, f0 + round(dure*fps_ca) - 1);
        if f1<=f0, stat('Range outside recording.'); return; end

        stopplay(); set(btnPlay,'Value',0,'String','Play');

        out_fps = max(1, round(EXPORT_FPS/spd));
        if strlength(EXPORT_DIR)>0, od=char(EXPORT_DIR);
        else, od=fileparts(mfilename('fullpath')); end
        outfile = fullfile(od, sprintf('repr_t%05.1f-%05.1fs_%gx.mp4', t0e, t0e+dure, spd));

        % Export frame layout
        EXPCA = 200;   % Ca patch size
        EXPBH = 240; EXPBW = 320;  % behavior video size
        GAP = 6;
        left_H = EXPBH + 3*(EXPCA+GAP) + 2*GAP;
        left_W = max(EXPBW, EXPCA);
        trace_W = round(left_H * 2.0);
        trace_H = left_H;
        total_W = left_W + GAP + trace_W;
        total_H = left_H;
        total_W = total_W + mod(total_W,2);  % even dimensions for H.264
        total_H = total_H + mod(total_H,2);

        %% Pre-render trace strip (once)
        stat('Rendering trace strip...');
        ws_t0 = t0e - EXPORT_WIN/2;
        ws_t1 = t0e + dure + EXPORT_WIN/2;
        trClrs_e = {[0.30 0.80 0.95], caClrs{1}, caClrs{2}, caClrs{3}};
        trNames_e = {'Breath','ROI 1','ROI 2','ROI 3'};
        trOff = [3 2 1 0];
        trData = [breath_trace, dFF_traces];
        tSig = cell(4,1);
        for kk = 1:4
            sig = trData(:,kk);
            r = max(sig)-min(sig); if r<eps, r=1; end
            tSig{kk} = (sig-min(sig))/r + trOff(kk);
        end
        idx = t_ca>=ws_t0 & t_ca<=ws_t1;

        STRIP_W = max(800, round((ws_t1-ws_t0)*400));
        STRIP_H = trace_H;
        LMARG = 0.11; RMARG = 0.02; AXBOT = 0.10; AXTOP = 0.04;
        tf2 = figure('Color','w','Position',[100 100 STRIP_W STRIP_H],'Visible','off');
        ax2 = axes(tf2,'Units','normalized','Position',[LMARG AXBOT 1-LMARG-RMARG 1-AXBOT-AXTOP]);
        hold(ax2,'on'); grid(ax2,'on');
        for kk=1:4
            plot(ax2, t_ca(idx), tSig{kk}(idx), '-', 'Color',trClrs_e{kk}, 'LineWidth',1.5);
        end
        xlim(ax2,[ws_t0 ws_t1]); ylim(ax2,[-0.25 4.3]);
        yticks(ax2, fliplr(trOff+0.5)); yticklabels(ax2, fliplr(trNames_e));
        xlabel(ax2,'Time (s)');
        set(ax2,'Color','w','XColor','k','YColor','k','FontSize',10);
        drawnow; fr2 = getframe(tf2); strip = imresize(fr2.cdata,[STRIP_H STRIP_W]); delete(tf2);

        Hs = size(strip,1); Ws = size(strip,2);
        Lpx = round(LMARG*Ws);
        plotWs = round((1-LMARG-RMARG)*Ws);
        colOfT = @(tt) Lpx + round((tt-ws_t0)/(ws_t1-ws_t0)*plotWs);
        wcols = colOfT(t0e+EXPORT_WIN/2) - colOfT(t0e-EXPORT_WIN/2);
        wcols = max(wcols, 4);
        rTop = max(1, round(AXTOP*Hs));
        rBot = min(Hs, round((1-AXBOT)*Hs));

        %% Write video
        vw = VideoWriter(outfile,'MPEG-4'); vw.FrameRate=out_fps; vw.Quality=95; open(vw);
        dstep = max(1, round(fps_ca/EXPORT_FPS));
        wb = waitbar(0, sprintf('Exporting %.1f-%.1f s ...', t0e, t0e+dure));
        kk_fr = f0; ii = 0;
        while kk_fr <= f1
            if mod(ii,dstep)==0
                full_frame = zeros(total_H, total_W, 3, 'uint8');

                %% Left column
                % Behavior video
                bfr_e = max(1,min(nBeh, round((kk_fr-1)*fps_beh/fps_ca + t0Offset_s*fps_beh)+1));
                if ~isempty(vBeh)
                    try
                        bhimg = read(vBeh,bfr_e);
                        if size(bhimg,3)==1, bhimg=repmat(bhimg,[1 1 3]); end
                        bhimg = imresize(bhimg,[EXPBH EXPBW]);
                    catch
                        bhimg = zeros(EXPBH,EXPBW,3,'uint8');
                    end
                else
                    bhimg = zeros(EXPBH,EXPBW,3,'uint8');
                end
                bh_c0 = max(1,round((left_W-EXPBW)/2)+1);
                full_frame(GAP:GAP+EXPBH-1, bh_c0:bh_c0+EXPBW-1, :) = bhimg;

                % Ca patches
                for k = 1:3
                    rgb_e = ca2rgb_sz(ca_patches{k}(:,:,kk_fr), ca_lo(k), ca_hi(k), ca_perim_up{k}, EXPCA);
                    ro = EXPBH + GAP + (k-1)*(EXPCA+GAP) + GAP;
                    co = max(1, round((left_W-EXPCA)/2)+1);
                    if ro+EXPCA-1<=total_H && co+EXPCA-1<=left_W
                        full_frame(ro:ro+EXPCA-1, co:co+EXPCA-1, :) = rgb_e;
                    end
                end

                %% Right: trace strip window + playhead
                tc_now = t_ca(kk_fr);
                c_ctr = colOfT(tc_now);
                c0s = max(Lpx+1, c_ctr - round(wcols/2));
                c1s = min(Ws, c0s + wcols - 1);
                ws_slice = [strip(:,1:Lpx,:), strip(:,c0s:c1s,:)];
                cx_in_slice = Lpx + round((c1s-c0s)/2);
                cxlo = max(1,cx_in_slice-1); cxhi = min(size(ws_slice,2),cx_in_slice+1);
                ws_slice(rTop:rBot, cxlo:cxhi, 1) = 220;
                ws_slice(rTop:rBot, cxlo:cxhi, 2) = 30;
                ws_slice(rTop:rBot, cxlo:cxhi, 3) = 30;
                tr_panel = imresize(ws_slice, [trace_H trace_W]);

                full_frame(:, left_W+GAP+1:left_W+GAP+trace_W, :) = tr_panel(1:total_H, 1:trace_W, :);

                writeVideo(vw, full_frame);
                if isvalid(wb), waitbar((kk_fr-f0)/max(f1-f0,1),wb); end
            end
            kk_fr = kk_fr+1; ii = ii+1;
        end
        close(vw); if isvalid(wb), delete(wb); end
        stat(sprintf('Exported: %s', outfile));
        fprintf('Exported: %s\n', outfile);
    end

    function onClose(~,~), stopplay(); delete(fig); end
    function stat(msg), if isgraphics(txStat), set(txStat,'String',msg); drawnow limitrate; end; end
end  % repr_video_gui_RZ


%% ============================================================
%  LOCAL HELPERS
%% ============================================================

function rgb = ca2rgb(ptch, lo, hi, perim_up_150)
% Normalize + gamma + upscale to DISP_CA×DISP_CA + yellow perimeter
    DISP_CA = 150;
    ptch_n  = max(0, min(1, (double(ptch)-lo)/(hi-lo))) .^ 0.6;
    ptch_up = imresize(ptch_n, [DISP_CA DISP_CA], 'bilinear');
    u8  = uint8(255 * ptch_up);
    rgb = repmat(u8,[1 1 3]);
    if ~isempty(perim_up_150) && ~isempty(ptch_n)
        pm = perim_up_150;
        rgb(:,:,1) = rgb(:,:,1).*uint8(~pm) + uint8(pm)*255;
        rgb(:,:,2) = rgb(:,:,2).*uint8(~pm) + uint8(pm)*255;
        rgb(:,:,3) = rgb(:,:,3).*uint8(~pm);
    end
end

function rgb = ca2rgb_sz(ptch, lo, hi, perim_up_orig, sz)
% Same but output is sz×sz (for export)
    ptch_n  = max(0, min(1, (double(ptch)-lo)/(hi-lo))) .^ 0.6;
    ptch_up = imresize(ptch_n, [sz sz], 'bilinear');
    u8  = uint8(255 * ptch_up);
    rgb = repmat(u8,[1 1 3]);
    if ~isempty(perim_up_orig)
        pm = imresize(double(perim_up_orig), [sz sz], 'nearest') > 0.5;
        rgb(:,:,1) = rgb(:,:,1).*uint8(~pm) + uint8(pm)*255;
        rgb(:,:,2) = rgb(:,:,2).*uint8(~pm) + uint8(pm)*255;
        rgb(:,:,3) = rgb(:,:,3).*uint8(~pm);
    end
end

function tifPath = local_discover_mc_tiff(fovFolder)
    mc_hits = dir(fullfile(fovFolder,'*_MC*.tif'));
    assert(~isempty(mc_hits), 'No *_MC*.tif in %s', fovFolder);
    names = {mc_hits.name};
    keep  = true(size(names));
    for ii = 1:numel(names)
        if contains(names{ii},'_AVG','IgnoreCase',true) || ...
           contains(names{ii},'_QC', 'IgnoreCase',true) || ...
           contains(names{ii},'_shifts','IgnoreCase',true) || ...
           contains(names{ii},'_SELROI','IgnoreCase',true)
            keep(ii) = false;
        end
    end
    mc_hits = mc_hits(keep);
    assert(~isempty(mc_hits),'No valid MC TIFF after filtering in %s',fovFolder);
    nMC = zeros(numel(mc_hits),1);
    for ii=1:numel(mc_hits)
        nMC(ii) = numel(strfind(upper(mc_hits(ii).name),'_MC'));
    end
    [~,bst] = max(nMC);
    tifPath = fullfile(mc_hits(bst).folder, mc_hits(bst).name);
end
