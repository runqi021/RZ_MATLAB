function drawROI_N_lpsub_gui()
% drawROI_N_lpsub_gui  Interactive GUI version of drawROI_N_lpsub_batch.
%   LEFT : the two nostril avg-projection images, each with a DRAGGABLE/RESIZABLE
%          ROI (ellipse or circle). RIGHT: LEFT / RIGHT / AVERAGE detrended
%          breathing traces that update LIVE as you move or resize either ROI.
%   - Auto-loads a previously saved <stem>_nostrilROI.mat ROI if it exists.
%   - Shape toggle (ellipse <-> circle), stat (mean/median/max), invert, LP cut.
%   - Prev/Next walk the DLC-csv list; Save writes <stem>_nostrilROI.mat (ROI) +
%     <stem>_breath.mat (avg breathing), same format as the batch script.
% LP-subtraction detrend, zero-phase, no flips. Run in MATLAB (interactive).

% ============================ USER-EDITABLE ============================
VIDEOS_DIR = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
dataRoot   = "D:\260615_thermalNbasler";
FINE_BP    = [2 10];      % breathing band saved as breath_bp ([] = skip); also the auto-ROI power band
DEF_SHAPE  = 'ellipse';   % 'ellipse' | 'circle'
DEF_STAT   = 'mean';      % per-frame ROI readout
DEF_LPCUT  = 1;           % Hz, low-pass baseline subtracted off
DEF_INVERT = true;        % inhale cools nostril -> invert so inhale reads as rise
START_AT   = 1;
% ======================================================================

here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
addpath(fullfile(repoRoot,'thermal_breathing'), fullfile(repoRoot,'mod','bluewhitered'));

d = dir(fullfile(char(VIDEOS_DIR), '*DLC*.csv'));
files = sort(string({d.name}));
assert(~isempty(files), 'no DLC csvs in %s', VIDEOS_DIR);

% ------------------------- shared state -------------------------
idx   = max(1,START_AT);
S=[]; P=[]; fps=[]; lb=[]; la=[]; lbd=[]; lad=[];
FL=[]; FR=[]; tvec=[]; FLd=[]; FRd=[]; tvd=[]; decim=1; fps_d=[];
FL0=[]; FR0=[]; cL=[]; cR=[]; szc=[]; avgL=[]; avgR=[];
roiL=[]; roiR=[];
lnL=[]; lnR=[]; lnA=[];
imgL=[]; imgR=[]; PmapL=[]; PmapR=[]; showPower=false;
BAND = [2 10];   % breathing band (Hz) for the per-pixel power map / auto-ROI
shape = DEF_SHAPE; stat = DEF_STAT; lpcut = DEF_LPCUT; invert = logical(DEF_INVERT);
roiSize = 5;   % default ROI half-size (px): circle radius / ellipse semi-axis (crop is only 24 px)
anchor = 'self';   % crop stabilization: 'self'=own dot | 'contra'=other dot | 'midpoint'=both

% ------------------------- build UI -------------------------
fig = uifigure('Name','nostril ROI breathing GUI','Position',[60 60 1520 840],'Color','w');
g  = uigridlayout(fig,[2 1]); g.RowHeight = {42,'1x'}; g.ColumnWidth = {'1x'};

cb = uigridlayout(g,[1 21]); cb.Layout.Row = 1;
cb.ColumnWidth = {180,56,56,40,80,34,74,112,42,48,32,46,28,28,68,74,84,78,54,98,'1x'};
cb.Padding = [6 4 6 4]; cb.ColumnSpacing = 5;

lblFile = uilabel(cb,'Text','','FontWeight','bold');                       lblFile.Layout.Column=1;
btnPrev = uibutton(cb,'Text','◀ Prev','ButtonPushedFcn',@(s,e)onPrev()); btnPrev.Layout.Column=2;
btnNext = uibutton(cb,'Text','Next ▶','ButtonPushedFcn',@(s,e)onNext()); btnNext.Layout.Column=3;
lbS=uilabel(cb,'Text','Shape','HorizontalAlignment','right'); lbS.Layout.Column=4;
ddShape = uidropdown(cb,'Items',{'ellipse','circle'},'Value',shape,'ValueChangedFcn',@(s,e)onShape()); ddShape.Layout.Column=5;
lbT=uilabel(cb,'Text','Stat','HorizontalAlignment','right'); lbT.Layout.Column=6;
ddStat  = uidropdown(cb,'Items',{'mean','median','max','min'},'Value',stat,'ValueChangedFcn',@(s,e)onStat()); ddStat.Layout.Column=7;
cbInv   = uicheckbox(cb,'Text','Invert (inhale up)','Value',invert,'ValueChangedFcn',@(s,e)onInvert()); cbInv.Layout.Column=8;
lbL=uilabel(cb,'Text','LP Hz','HorizontalAlignment','right'); lbL.Layout.Column=9;
efLP    = uieditfield(cb,'numeric','Value',lpcut,'Limits',[0.05 10],'ValueChangedFcn',@(s,e)onLP()); efLP.Layout.Column=10;
lbZ=uilabel(cb,'Text','Size','HorizontalAlignment','right'); lbZ.Layout.Column=11;
efSize  = uieditfield(cb,'numeric','Value',roiSize,'Limits',[1 30],'ValueChangedFcn',@(s,e)onSize()); efSize.Layout.Column=12;
btnSm   = uibutton(cb,'Text','−','ButtonPushedFcn',@(s,e)onStep(-1)); btnSm.Layout.Column=13;
btnBg   = uibutton(cb,'Text','+','ButtonPushedFcn',@(s,e)onStep(+1)); btnBg.Layout.Column=14;
btnAuto = uibutton(cb,'Text','Auto ROI','ButtonPushedFcn',@(s,e)onAuto()); btnAuto.Layout.Column=15;
ddAuto  = uidropdown(cb,'Items',{'self','contra','midpoint'},'Value',anchor,'ValueChangedFcn',@(s,e)onAnchor()); ddAuto.Layout.Column=16;
cbPow   = uicheckbox(cb,'Text','SNR view','Value',showPower,'ValueChangedFcn',@(s,e)onPowView()); cbPow.Layout.Column=17;
btnRst  = uibutton(cb,'Text','Reset ROIs','ButtonPushedFcn',@(s,e)onReset());  btnRst.Layout.Column=18;
btnSave = uibutton(cb,'Text','Save','ButtonPushedFcn',@(s,e)onSave(false));    btnSave.Layout.Column=19;
btnSavN = uibutton(cb,'Text','Save & Next','ButtonPushedFcn',@(s,e)onSave(true)); btnSavN.Layout.Column=20;
lblStat = uilabel(cb,'Text','','FontColor',[0.1 0.45 0.1]);                 lblStat.Layout.Column=21;

bg = uigridlayout(g,[1 2]); bg.Layout.Row = 2; bg.ColumnWidth = {'1.15x','1x'}; bg.ColumnSpacing = 8;
lg = uigridlayout(bg,[2 1]); lg.Layout.Column = 1; lg.RowSpacing = 8;
axL = uiaxes(lg); axL.Layout.Row = 1;
axR = uiaxes(lg); axR.Layout.Row = 2;
rg = uigridlayout(bg,[3 1]); rg.Layout.Column = 2; rg.RowSpacing = 6;
tL = uiaxes(rg); tL.Layout.Row = 1;
tR = uiaxes(rg); tR.Layout.Row = 2;
tA = uiaxes(rg); tA.Layout.Row = 3;

loadVideo(idx);

% ===================== nested callbacks/logic =====================
    function onPrev(),  if idx>1,            loadVideo(idx-1); end, end
    function onNext(),  if idx<numel(files), loadVideo(idx+1); end, end
    function onStat(),  stat   = ddStat.Value;  updateTraces(); end
    function onInvert(),invert = cbInv.Value;   updateTraces(); end
    function onLP(),    lpcut  = efLP.Value;     setFilter(); updateTraces(); end
    function onReset(), makeROIs([]);            updateTraces(); end
    function onSize()
        roiSize = efSize.Value;
        setROISize(roiL, roiSize); setROISize(roiR, roiSize); updateTraces();
    end
    function onStep(s)
        roiSize = min(30, max(1, roiSize + s)); efSize.Value = roiSize;
        setROISize(roiL, roiSize); setROISize(roiR, roiSize); updateTraces();
    end
    function setROISize(roi, sz)
        if isempty(roi) || ~isvalid(roi), return; end
        if isa(roi,'images.roi.Circle'), roi.Radius = sz; else, roi.SemiAxes = [sz sz]; end
    end
    function onPowView(), showPower = cbPow.Value; refreshImg(); end
    function onAnchor(), anchor = ddAuto.Value; applyStabilize(); end
    function onAuto()
        fitPowerROI(roiL, PmapL, roiSize); fitPowerROI(roiR, PmapR, roiSize);
        updateTraces(); drawnow; lblStat.Text = sprintf('auto ROI (%s) from breath SNR', anchor);
    end
    function onShape()
        shape = ddShape.Value;
        prev = struct('L',roiInfo(roiL),'R',roiInfo(roiR));   % keep centers
        makeROIs(prev); updateTraces();
    end

    function setFilter()
        [lb,la] = butter(2, lpcut/(fps/2), 'low');                 % full-res (used for Save)
        if ~isempty(fps_d) && fps_d > 0
            c = min(lpcut, 0.45*fps_d);
            [lbd,lad] = butter(2, c/(fps_d/2), 'low');             % decimated (live display)
        else
            lbd = lb; lad = la;
        end
    end

    function [cl,cr] = loadCenters()
        % per-stack-frame nostril dot positions (native px), aligned to the saved stack
        T0 = size(FL0,1);
        if isfield(S,'L_center') && isfield(S,'R_center')
            cl = double(S.L_center); cr = double(S.R_center);
            q = max(1, round(size(cl,1)/T0)); cl = cl(1:q:end,:); cr = cr(1:q:end,:);
            n = min([T0, size(cl,1), size(cr,1)]); cl = cl(1:n,:); cr = cr(1:n,:);
            if n < T0, cl(end+1:T0,:) = repmat(cl(end,:),T0-n,1); cr(end+1:T0,:) = repmat(cr(end,:),T0-n,1); end
        else
            cl = []; cr = [];     % no centers saved -> contra/midpoint disabled
        end
    end

    function applyStabilize()
        % Build working stacks FL/FR. 'self' = raw own-dot crops. 'contra'/'midpoint'
        % re-stabilize each crop per frame to a steadier reference (the contralateral
        % dot, or the both-dot midpoint) so the breathing pixel stops jittering.
        if strcmp(anchor,'self') || isempty(cL) || isempty(cR)
            FL = FL0; FR = FR0;
            if ~strcmp(anchor,'self'), lblStat.Text = 'no L/R_center in .mat -> need re-extract for contra/midpoint'; end
        elseif strcmp(anchor,'midpoint')
            mid = 0.5*(cL+cR); FL = stabil(FL0, cL, mid, szc); FR = stabil(FR0, cR, mid, szc);
        else  % contra
            FL = stabil(FL0, cL, cR, szc); FR = stabil(FR0, cR, cL, szc);
        end
        avgL = reshape(mean(FL,1,'omitnan'), szc); avgR = reshape(mean(FR,1,'omitnan'), szc);
        FLd = FL(1:decim:end,:); FRd = FR(1:decim:end,:);
        PmapL = breathPowerMap(FL, fps, BAND, szc); PmapR = breathPowerMap(FR, fps, BAND, szc);
        refreshImg(); updateTraces();
    end

    function Fs = stabil(F, thisC, refC, sz)
        % shift each frame to remove this-dot jitter relative to refC (keeps constant offset)
        Tn = size(F,1); d = thisC - refC; sh = round(median(d,1,'omitnan') - d); sh(~isfinite(sh)) = 0;
        Fs = F;
        for k = 1:Tn
            s = sh(k,:);
            if any(s)
                fr = circshift(reshape(F(k,:),sz), [s(2), s(1)]); dy = s(2); dx = s(1);
                if dy>0, fr(1:dy,:)=NaN; elseif dy<0, fr(end+dy+1:end,:)=NaN; end   % null wrapped border
                if dx>0, fr(:,1:dx)=NaN; elseif dx<0, fr(:,end+dx+1:end)=NaN; end
                Fs(k,:) = fr(:).';
            end
        end
    end

    function loadVideo(i)
        idx = i;
        csv = fullfile(VIDEOS_DIR, files(i));
        P   = thermal_resolve_paths(csv, dataRoot);
        S   = load(P.nostrilC);
        fps = double(S.L_stack_fps);
        stkL = double(S.L_stack); stkR = double(S.R_stack);
        FL0 = reshape(stkL, size(stkL,1), []);    % RAW own-dot crops [T x npix]
        FR0 = reshape(stkR, size(stkR,1), []);
        szc = size(double(S.L_avg));
        tvec = (0:size(FL0,1)-1)/fps;
        decim = max(1, round(numel(tvec)/3000));            % decimate traces for fast display
        tvd = tvec(1:decim:end);
        fps_d = fps/decim;
        setFilter();
        [cL,cR] = loadCenters();                  % per-frame nostril centers (native px), stack-aligned
        cla(axL); imgL = imagesc(axL, double(S.L_avg)); axis(axL,'image'); colorbar(axL);
        cla(axR); imgR = imagesc(axR, double(S.R_avg)); axis(axR,'image'); colorbar(axR);
        lnL = newTrace(tL, 'LEFT detrended');
        lnR = newTrace(tR, 'RIGHT detrended');
        lnA = newTrace(tA, 'L+R AVERAGE');  xlabel(tA,'s');
        prev = [];
        if isfile(P.nostrilROI), prev = load(P.nostrilROI); end
        if ~isempty(prev) && isfield(prev,'L') && isfield(prev.L,'shape')
            shape = prev.L.shape; ddShape.Value = shape;
        end
        makeROIs(prev);
        if ~isempty(roiL) && isvalid(roiL)               % reflect loaded ROI size in the field
            if isa(roiL,'images.roi.Circle'), roiSize = roiL.Radius; else, roiSize = mean(roiL.SemiAxes); end
            efSize.Value = round(roiSize,1);
        end
        applyStabilize();          % builds FL/FR (+decimated) + SNR maps per anchor, refresh + plot
        lblFile.Text = sprintf('%d/%d   %s', i, numel(files), files(i));
        if isempty(prev), lblStat.Text = 'new ROI'; else, lblStat.Text = 'loaded saved ROI'; end
    end

    function makeROIs(prev)
        if ~isempty(roiL) && isvalid(roiL), delete(roiL); end
        if ~isempty(roiR) && isvalid(roiR), delete(roiR); end
        roiL = makeOne(axL, getSaved(prev,'L'), size(double(S.L_avg)));
        roiR = makeOne(axR, getSaved(prev,'R'), size(double(S.R_avg)));
    end

    function roi = makeOne(ax, sv, sz)
        cc = [(sz(2)+1)/2, (sz(1)+1)/2];                       % [x y] image center
        c = cc; if ~isempty(sv) && isfield(sv,'center') && ~isempty(sv.center), c = sv.center; end
        if strcmpi(shape,'circle')
            r = roiSize; if ~isempty(sv) && isfield(sv,'radius') && ~isempty(sv.radius), r = sv.radius; end
            roi = drawcircle(ax,'Center',c,'Radius',r,'Color','c','FaceAlpha',0.08,'LineWidth',1.2);
        else
            s = [roiSize roiSize]; if ~isempty(sv) && isfield(sv,'semiaxes') && ~isempty(sv.semiaxes), s = sv.semiaxes; end
            ang = 0; if ~isempty(sv) && isfield(sv,'angle') && ~isempty(sv.angle), ang = sv.angle; end
            roi = drawellipse(ax,'Center',c,'SemiAxes',s,'RotationAngle',ang,'Color','c','FaceAlpha',0.08,'LineWidth',1.2);
        end
        addlistener(roi,'ROIMoved', @(s,e)updateTraces());     % plot ONLY when you finish moving/resizing (no laggy live update)
    end

    function updateTraces()
        if isempty(roiL) || ~isvalid(roiL) || isempty(roiR) || ~isvalid(roiR), return; end
        [~,dL] = sideTrace(roiL, FLd, lbd, lad);          % decimated -> fast live update
        [~,dR] = sideTrace(roiR, FRd, lbd, lad);
        dA = (dL + dR)/2;
        set(lnL,'XData',tvd,'YData',dL);
        set(lnR,'XData',tvd,'YData',dR);
        set(lnA,'XData',tvd,'YData',dA);
        title(tA, sprintf('L+R AVERAGE  (breathing, peak %.2f Hz)', peakHz(dA)));
    end

    function [tr,detr,mask] = sideTrace(roi, F, lbf, laf)
        mask = createMask(roi); cols = mask(:);
        if ~any(cols)
            tr = zeros(size(F,1),1);                 % empty ROI -> flat
        else
            roipix = F(:, cols);
            switch stat
                case 'mean',   tr = mean(roipix, 2, 'omitnan');
                case 'median', tr = median(roipix, 2, 'omitnan');
                case 'max',    tr = max(roipix, [], 2, 'omitnan');
                case 'min',    tr = min(roipix, [], 2, 'omitnan');
            end
        end
        tr = fillmissing(tr, 'linear', 'EndValues', 'nearest');   % edge/out-of-frame NaNs
        tr(~isfinite(tr)) = 0;                                    % all-NaN fallback -> filtfilt-safe
        base = filtfilt(lbf, laf, tr); detr = tr - base;
        if invert, detr = -detr; end
    end

    function onSave(advance)
        [trL,dL,maskL] = sideTrace(roiL, FL, lb, la);   % full-res for saving
        [trR,dR,maskR] = sideTrace(roiR, FR, lb, la);
        dA = (dL + dR)/2;
        R = struct();
        R.L = packSide('L', roiL, trL, dL, maskL);
        R.R = packSide('R', roiR, trR, dR, maskR);
        R.src = char(fullfile(VIDEOS_DIR, files(idx))); R.nostrilC = P.nostrilC;
        save(P.nostrilROI, '-struct', 'R');
        if ~isempty(FINE_BP)
            [fbb,faa] = butter(2, FINE_BP/(fps/2), 'bandpass'); dA_bp = filtfilt(fbb, faa, dA);
        else
            dA_bp = [];
        end
        B = struct('breath',dA(:), 'breath_bp',dA_bp(:), 'fps',fps, 't',tvec(:), ...
            'method','lpsub', 'lp_cut',lpcut, 'fine_bp',FINE_BP, 'inverted',invert, ...
            'roi_stat',stat, 'roi_shape',shape, 'animal',P.animal, 'run',P.k, ...
            'src_csv',char(fullfile(VIDEOS_DIR, files(idx))), 'src_ats',P.ats);
        save(P.breath, '-struct', 'B');
        lblStat.Text = sprintf('saved %s n%d', P.animal, P.k);
        if advance && idx < numel(files), loadVideo(idx+1); end
    end

    function sd = packSide(name, roi, tr, detr, mask)
        sd = struct('side',char(S.(name+"_side")), 'mask',mask, 'trace',tr, 'detr',detr, ...
            'fps',double(S.(name+"_stack_fps")), 'npix',nnz(mask), 'shape',shape, 'stat',stat, ...
            'center',roi.Center);
        if isa(roi,'images.roi.Circle')
            sd.radius = roi.Radius;
        else
            sd.semiaxes = roi.SemiAxes; sd.angle = roi.RotationAngle;
        end
    end

    % --------- small helpers (nested for shared fps/lpcut) ---------
    function f = peakHz(x)
        x = x - mean(x,'omitnan'); nf = 2^nextpow2(numel(x));
        Pw = abs(fft(x,nf)).^2; fr = (0:nf-1)*(fps_d/nf);          % x is the decimated trace
        inb = fr>=lpcut & fr<=min(10, 0.95*fps_d/2); frb = fr(inb); Pb = Pw(inb);
        [~,ip] = max(Pb); if isempty(ip), f = NaN; else, f = frb(ip); end
    end
    function refreshImg()
        if isempty(imgL) || ~isvalid(imgL), return; end
        if showPower
            set(imgL,'CData',PmapL); set(imgR,'CData',PmapR);
            mxL = max(PmapL(:)); if ~isfinite(mxL) || mxL<=0, mxL = 1; end
            mxR = max(PmapR(:)); if ~isfinite(mxR) || mxR<=0, mxR = 1; end
            axL.CLim = [0 mxL]; axR.CLim = [0 mxR];
            colormap(axL, hot(256)); colormap(axR, hot(256));     % set AFTER CLim
            title(axL, sprintf('%s breath SNR (2-10/10-40 Hz)', char(S.L_side)));
            title(axR, sprintf('%s breath SNR (2-10/10-40 Hz)', char(S.R_side)));
        else
            set(imgL,'CData',avgL); set(imgR,'CData',avgR);
            axL.CLimMode = 'auto'; axR.CLimMode = 'auto';
            colormap(axL, bluewhitered(256)); colormap(axR, bluewhitered(256));
            title(axL, sprintf('%s nostril avg (\\circC) — move ROI, release to plot', char(S.L_side)));
            title(axR, sprintf('%s nostril avg (\\circC) — move ROI, release to plot', char(S.R_side)));
        end
        drawnow limitrate;
    end
    function ln = newTrace(ax, ttl)
        cla(ax); ln = plot(ax, NaN, NaN, '-', 'LineWidth',1.0);
        grid(ax,'on'); ylabel(ax,'\circC'); title(ax, ttl);
    end
end

% ===================== file-scope helpers =====================
function sv = getSaved(prev, name)
    if ~isempty(prev) && isfield(prev, name), sv = prev.(name); else, sv = []; end
end
function info = roiInfo(roi)
    if isempty(roi) || ~isvalid(roi), info = []; return; end
    info = struct('center', roi.Center);
    if isa(roi,'images.roi.Circle')
        info.radius = roi.Radius; info.shape = 'circle';
    else
        info.semiaxes = roi.SemiAxes; info.angle = roi.RotationAngle; info.shape = 'ellipse';
    end
end

function Smap = breathPowerMap(F, fps, band, sz)
% Per-pixel BREATH SNR map = in-band power / out-of-band (noise) power.
%   in-band  = power in [2 10] Hz (the breathing band)
%   noise    = power in [10 .. min(40,0.95*Nyq)] Hz
%   SNR      = in-band / noise, with pixels whose ABSOLUTE in-band power is
%              < 10% of the max zeroed (so dead/cold tiny-over-tiny pixels don't win)
% This selects pixels where breathing DOMINATES, not just pixels that are noisy.
    nyq = fps/2;
    [bb,ab] = butter(2, band/nyq, 'bandpass');
    lo = band(2); hi = min(40, 0.95*nyq); useNoise = hi > lo*1.05;
    deadCol = all(isnan(F),1);
    Fc = F;
    if any(isnan(F(:))), Fc = fillmissing(Fc,'linear',1,'EndValues','nearest'); end
    Fc(:,deadCol) = 0;
    Pin = var(filtfilt(bb,ab,Fc), 0, 1);                 % breath-band power
    if useNoise
        [bn,an] = butter(2, [lo hi]/nyq, 'bandpass');
        Pno = var(filtfilt(bn,an,Fc), 0, 1);             % out-of-band noise power
    else
        Pno = max(var(Fc,0,1) - Pin, 0);                 % fallback: total - in-band
    end
    SNR = Pin ./ (Pno + eps);
    SNR(Pin < 0.10*max(Pin)) = 0;                        % absolute-signal floor
    SNR(deadCol) = 0;
    Smap = reshape(SNR, sz);
end

function fitPowerROI(roi, Smap, roiSize)
% Snap a small ROI onto the breath-SNR hot-spot NEAR the tracked nostril dot.
% The crop is centered on the DLC point, and high-contrast nose EDGES produce
% big 2-10 Hz swings (motion artifacts) with huge SNR -- so we ONLY search
% within MAX_OFF px of the crop center, then take the brightest SNR pixel there
% (refined by a tiny 3x3 SNR-weighted centroid). Size fixed (user Size, <=MAX_SEMI).
    if isempty(roi) || ~isvalid(roi), return; end
    MAX_SEMI = 3;       % px: auto ROI radius / semi-axis hard cap
    MAX_OFF  = 4;       % px: auto ROI center may not stray further from the tracked dot
    NBR      = 1;       % px: local refine half-window around the in-region peak
    P = Smap; P(~isfinite(P)) = 0;
    sz = size(P); ctr = [(sz(2)+1)/2, (sz(1)+1)/2];
    [X,Y] = meshgrid(1:sz(2), 1:sz(1));
    P(hypot(X-ctr(1), Y-ctr(2)) > MAX_OFF) = 0;     % search ONLY near the tracked dot
    [mx,im] = max(P(:));
    if mx <= 0, return; end
    [py,px] = ind2sub(sz, im);                       % brightest SNR pixel within MAX_OFF
    r0 = max(1,py-NBR); r1 = min(sz(1),py+NBR);
    c0 = max(1,px-NBR); c1 = min(sz(2),px+NBR);
    sub = P(r0:r1, c0:c1); [yy,xx] = ndgrid(r0:r1, c0:c1); sw = sum(sub(:));
    if sw > 0, cx = sum(xx(:).*sub(:))/sw; cy = sum(yy(:).*sub(:))/sw; else, cx = px; cy = py; end
    rsz = min(max(roiSize,1), MAX_SEMI);
    roi.Center = [cx cy];
    if isa(roi,'images.roi.Circle')
        roi.Radius = rsz;
    else
        roi.SemiAxes = [rsz rsz]; roi.RotationAngle = 0;
    end
end
