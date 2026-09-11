% io_repr_video_RZ.m
% -----------------------------------------------------------------------
% Representative video for INFERIOR-OLIVE calcium data (no behavior cam).
% Left  = manually-cropped calcium movie (free rectangle on the MC avg-proj).
% Right = population mean dF/F trace, rolling window, with a moving playhead
%         and red stimulus markers.
%
% Timeline: dF/F drops the first nDrop frames (helper.dFF_RZ default 30), so
% the dF/F sample t=0 corresponds to RAW time nDrop/fps. All user-facing times
% (clip, stim, x-axis) are RAW experiment seconds.
%
% Standalone; manual RGB compositing (consistent frame size) + ffmpeg mp4.
% -----------------------------------------------------------------------

clear; clc; close all;

%% ============================== USER PARAMS ==============================
fovPath   = 'C:\Users\Admin\Desktop\PO_z80_beep_5khz_00001';
repoPath  = 'C:\Users\Admin\Desktop\RZ_MATLAB';   % for helper.dFF_RZ + detect_session_fps

% clip window (RAW experiment seconds)
t_start   = 5;
t_end     = 30;

% stimulus (RAW seconds): 100 ms pulses at 1 Hz from 15 s, 10 of them
STIM_T0      = 15;
STIM_DUR     = 0.10;
STIM_PERIOD  = 1.0;
N_STIM       = 10;

% manual crop of the calcium field (px). [] = draw interactively on avg-proj.
CA_CROP            = [];      % [x y w h]
CA_AVGPROJ_NFRAMES = 300;

% calcium display
CA_PRC_LO   = 30;     CA_PRC_HI = 100;    % display window (lower lo = more black crush; ImageJ default ~0.35%)
gamma_ca    = 1;
CA_SMOOTH_SIGMAS = [0 0];   % ONE OUTPUT VIDEO PER sigma (px gaussian on the Ca video; 0 = off)
CA_CMAP_INVERT  = false;      % true = dark cells on light bg
CA_RESIZE_METHOD = 'lanczos3';% Ca tile resampling: 'lanczos3'(sharp) | 'nearest'(crisp ImageJ-like pixels) | 'bilinear'

% trace
TRACE_WIN   = 6.0;            % rolling window (s)
TRACE_YLIM  = [];             % [] = auto from clip; or [lo hi] dF/F  (MEAN mode only)

% --- trace display mode (revertable) ---
TRACE_MODE     = 'stack';     % 'mean' = population mean (one line) | 'stack' = all ROIs stacked offset
STACK_ROIS     = [];          % [] = all ROIs; or a vector of ROI indices to show (stack mode)
STACK_SORT     = 'latency';   % 'latency' = order lanes by stim-evoked peak time | 'none' = ROI index order
STACK_GAIN     = 0.9;         % each ROI trace spans ~this fraction of its lane (shared dF/F scale)
STACK_SCALEBAR = 1.0;         % dF/F height of the scale bar (0 = off)
STACK_LINE_W   = 1;           % stacked trace line width (px)

% --- video-only mode (no trace panel) ---
VIDEO_ONLY  = true;           % true = output ONLY the calcium movie + a stimulus indicator
STIM_MARKER = 'none';         % 'none' = no visual indicator (audio only) | 'dot' = corner dot | 'bar' = top bar
STIM_BAR_H  = 22;             % px height of the top bar (STIM_MARKER='bar')
STIM_DOT_R  = 14;             % px radius of the corner dot (STIM_MARKER='dot')

% --- stimulus audio (the DAO stimulus is a 5 kHz tone) ---
ADD_STIM_AUDIO = true;        % mux a sine tone at the stim pulse times into the mp4
AUDIO_FREQ     = 5000;        % Hz (5 kHz)
AUDIO_SR       = 44100;       % audio sample rate (Hz)
AUDIO_GAIN     = 0.25;        % tone amplitude, 0..1

% layout / style
CA_TILE_H   = 300;            % display height of the Ca crop (px); width keeps aspect. [] = NATIVE resolution (no rescale)
TRACE_PLOT_W= 1100;
X_AXIS_H    = 90;
GAP         = 20;
FONT_NAME   = 'Arial';   FONT_SIZE = 26;
% --- colors (R G B, 0-255) ---
PLAYHEAD_COLOR = [255 255 125];     % vertical line marking the CURRENT FRAME
STIM_COLOR     = [255 231 186];      % stimulus pulse band color
TRACE_COLOR    = [255 255 255];    % population dF/F trace
AXIS_COLOR     = [255 255 255];    % axes / ticks

outFPS    = 30;
OUT_CRF    = 12;              % x264 quality: lower = better (0 = lossless, 12 = near-lossless, 18 = old default)
OUT_PRESET = 'slow';         % x264 preset (slower = better quality/compression)
outFile   = '';               % '' = auto in fovPath
fallback_fps = 30;
RAW_TOSS_FRAMES = 30;         % frames removed in PREPROCESSING (MC_MC tiff/F are already post-toss);
                             % used only to put the x-axis / stim times in RAW experiment seconds

FFMPEG = 'C:\Users\Admin\AppData\Local\Microsoft\WinGet\Packages\Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe\ffmpeg-8.1-full_build\bin\ffmpeg.exe';
if ~isfile(FFMPEG), FFMPEG = 'ffmpeg'; end

MAKE_PREVIEW_ONLY = false;     % true = show one frame and stop
PREVIEW_T_SEC     = 15.5;     % preview at this RAW second
%% =========================================================================

addpath(repoPath); addpath(fullfile(repoPath,'2p_breathing_coherence'));

%% 1 — fps + pixel size
try
    [fps, sm] = detect_session_fps(fovPath);
catch
    fps = fallback_fps; sm = struct('pixelSize_um',NaN);
end
fprintf('fps=%g Hz\n', fps);

%% 2 — cpSAM + dF/F
sam = dir(fullfile(fovPath,'*_cpSAM_output.mat'));
assert(~isempty(sam),'no *_cpSAM_output.mat in %s', fovPath);
SD  = load(fullfile(sam(1).folder, sam(1).name),'F');
F   = double(SD.F);
dout = helper.dFF_RZ(F, 'FPS', fps);
dFF  = dout.dFF;
[T_dff, nROI] = size(dFF);
nDrop  = size(F,1) - T_dff;        % dFF-vs-tiff alignment offset (0 here: both already post-toss)
tossSec= RAW_TOSS_FRAMES/fps;      % preprocessing toss -> RAW experiment-time axis offset
t_raw  = tossSec + (0:T_dff-1)'/fps;     % RAW time of each dF/F sample
fprintf('  dFF [%d x %d], nDrop=%d (toss %.2fs), raw span %.1f-%.1f s\n', ...
    T_dff, nROI, nDrop, tossSec, t_raw(1), t_raw(end));

popTrace = mean(dFF, 2);                  % population mean dF/F

%% 3 — clip frame range (in dF/F frame index)
f0 = max(1, round((t_start - tossSec)*fps) + 1);
f1 = min(T_dff, round((t_end   - tossSec)*fps) + 1);
assert(f1 > f0, 'empty clip');
nFrames = f1 - f0 + 1;
fprintf('clip raw %.1f-%.1f s  -> dFF frames %d-%d (%d)\n', t_start, t_end, f0, f1, nFrames);

%% 4 — MC tiff + manual crop on avg-proj
tifP    = discover_mc_tiff(fovPath);
tifInfo = imfinfo(tifP);
imgH = tifInfo(1).Height; imgW = tifInfo(1).Width;
fprintf('movie %s  (%dx%d, %d dirs)\n', tifP, imgW, imgH, numel(tifInfo));

avgProj = ca_avgproj(tifP, tifInfo, nDrop, CA_AVGPROJ_NFRAMES);
if isempty(CA_CROP), CA_CROP = selectCropRect(avgProj); end
cx = max(1,round(CA_CROP(1))); cy = max(1,round(CA_CROP(2)));
cw = round(CA_CROP(3)); ch = round(CA_CROP(4));
c2 = min(imgW, cx+cw-1); r2 = min(imgH, cy+ch-1);
fprintf('Ca crop: x=%d y=%d w=%d h=%d\n', cx, cy, c2-cx+1, r2-cy+1);

%% 5 — preload Ca crop frames for the clip
caBuf = zeros(r2-cy+1, c2-cx+1, nFrames, 'single');
tic;
for ff = 1:nFrames
    frm = single(imread(tifP, nDrop+f0+ff-1, 'Info', tifInfo));
    caBuf(:,:,ff) = frm(cy:r2, cx:c2);
    if mod(ff,200)==0, fprintf('  loaded %d/%d (%.1fs)\n', ff, nFrames, toc); end
end
caLo = double(prctile(caBuf(:), CA_PRC_LO));
caHi = double(prctile(caBuf(:), CA_PRC_HI));
if caHi<=caLo, caHi=caLo+1; end

%% 6 — layout geometry
if isempty(CA_TILE_H)                 % native crop resolution (no rescale)
    caTileH = r2-cy+1;  caTileW = c2-cx+1;
else
    caTileH = CA_TILE_H;
    caTileW = round(caTileH * (c2-cx+1)/(r2-cy+1));
end
traceW  = TRACE_PLOT_W + 60;                 % +right pad
tileStackH = max(caTileH, caTileH);          % single row
if VIDEO_ONLY
    total_H = caTileH;                       % just the calcium movie
    total_W = caTileW;
else
    total_H = tileStackH + X_AXIS_H;
    total_W = caTileW + GAP + traceW;
end
total_W = total_W + mod(total_W,2); total_H = total_H + mod(total_H,2);
fprintf('frame %dx%d px\n', total_W, total_H);

% trace y-limits
if isempty(TRACE_YLIM)
    yy = popTrace(f0:f1); pad = 0.10*max(range(yy),eps);
    yl = [min(yy)-pad, max(yy)+pad];
else
    yl = TRACE_YLIM;
end

% trace plot geometry (within the right panel)
axisX = round(0.10*TRACE_PLOT_W);
plotL = axisX; plotR = round(0.99*TRACE_PLOT_W);
yTop  = 14; yBot = caTileH - 14;
xAxisY = tileStackH + round(0.18*X_AXIS_H);

% ---- stack-mode trace prep (display order, shared scale, lane geometry) ----
ST = struct('on',false);
traceBase = [];                         % unused in VIDEO_ONLY mode
if ~VIDEO_ONLY
if strcmpi(TRACE_MODE,'stack')
    if isempty(STACK_ROIS), stackIdx = 1:nROI; else, stackIdx = STACK_ROIS(:)'; end
    stackIdx = stackIdx(stackIdx>=1 & stackIdx<=nROI);
    assert(~isempty(stackIdx), 'STACK_ROIS selects no valid ROIs');

    if strcmpi(STACK_SORT,'latency')                 % order by stim-evoked peak time
        wi = find(t_raw >= STIM_T0 & t_raw <= STIM_T0 + N_STIM*STIM_PERIOD);
        if numel(wi) >= 2
            pk = zeros(numel(stackIdx),1);
            for ii = 1:numel(stackIdx)
                [~,im] = max(dFF(wi, stackIdx(ii))); pk(ii) = im;
            end
            [~,ord] = sort(pk); stackIdx = stackIdx(ord);
        end
    end

    Nst      = numel(stackIdx);
    stackMat = dFF(:, stackIdx);                      % [T x Nst] display order
    base0    = median(stackMat(f0:f1,:), 1, 'omitnan');
    stackMat = stackMat - base0;                      % per-ROI baseline-subtract
    av       = abs(stackMat(f0:f1,:)); av = av(isfinite(av));
    amp      = prctile(av(:), 99);                    % shared amplitude scale
    if isempty(amp) || ~isfinite(amp) || amp<=0, amp = 1; end

    spacing    = (yBot - yTop)/Nst;
    laneBaseY  = yTop + ((1:Nst)-0.5)*spacing;        % lane 1 (earliest) near top
    pixPerUnit = STACK_GAIN*spacing/amp;

    ST = struct('on',true, 'stackMat',stackMat, 'stackIdx',stackIdx, ...
                'laneBaseY',laneBaseY, 'pixPerUnit',pixPerUnit, 'lineW',STACK_LINE_W);
end

% precompute static trace base (axis, labels)
traceBase = zeros(total_H, traceW, 3, 'uint8');
traceBase = drawLineRGB(traceBase, axisX, yTop, axisX, yBot, AXIS_COLOR, 1);
traceBase = drawLineRGB(traceBase, plotL, xAxisY, plotR, xAxisY, AXIS_COLOR, 1);
traceBase = drawTextRGBCenter(traceBase, round((plotL+plotR)/2), xAxisY+34, 'Time (s)', AXIS_COLOR, FONT_SIZE, FONT_NAME);
if ST.on
    % scale bar (top-left of plot) + per-lane ROI labels when there is room
    if STACK_SCALEBAR > 0
        sbPix = round(pixPerUnit*STACK_SCALEBAR);
        sbX = plotL+10; sbY1 = yTop+6; sbY2 = yTop+6+sbPix;
        traceBase = drawLineRGB(traceBase, sbX, sbY1, sbX, sbY2, AXIS_COLOR, 3);
        traceBase = drawTextRGB(traceBase, sbX+6, round((sbY1+sbY2)/2)-round(FONT_SIZE*0.4), ...
            sprintf('%.2g dF/F',STACK_SCALEBAR), AXIS_COLOR, round(FONT_SIZE*0.8), FONT_NAME);
    end
    if spacing >= FONT_SIZE+2
        for r = 1:Nst
            traceBase = drawTextRGBRight(traceBase, axisX-10, round(laneBaseY(r))-round(FONT_SIZE*0.4), ...
                sprintf('%d',stackIdx(r)), AXIS_COLOR, round(FONT_SIZE*0.8), FONT_NAME);
        end
    end
else
    traceBase = drawTextRGBRight(traceBase, axisX-10, yTop-2,  sprintf('%.2g',yl(2)), AXIS_COLOR, FONT_SIZE, FONT_NAME);
    traceBase = drawTextRGBRight(traceBase, axisX-10, yBot-18, sprintf('%.2g',yl(1)), AXIS_COLOR, FONT_SIZE, FONT_NAME);
    traceBase = drawTextRGBRight(traceBase, axisX-44, round((yTop+yBot)/2)-14, 'dF/F', AXIS_COLOR, FONT_SIZE, FONT_NAME);
end
end   % ~VIDEO_ONLY

stim_on = STIM_T0 + (0:N_STIM-1)*STIM_PERIOD;   % RAW onset seconds

%% 7 — preview
if MAKE_PREVIEW_ONLY
    ffp = max(1, min(nFrames, round((PREVIEW_T_SEC - tossSec)*fps)+1 - (f0-1)));
    fr = ioFrame(ffp, total_H,total_W, caTileW,caTileH, GAP, traceBase, traceW, ...
        caBuf, caLo, caHi, gamma_ca, CA_SMOOTH_SIGMAS(end), CA_CMAP_INVERT, ...
        t_raw, f0, popTrace, yl, TRACE_WIN, plotL, plotR, yTop, yBot, xAxisY, ...
        stim_on, STIM_DUR, STIM_COLOR, PLAYHEAD_COLOR, TRACE_COLOR, FONT_NAME, FONT_SIZE, ST, ...
        VIDEO_ONLY, STIM_BAR_H, STIM_MARKER, STIM_DOT_R, CA_RESIZE_METHOD);
    figure('Color','k'); imshow(fr); title(sprintf('preview @ %.2fs', PREVIEW_T_SEC));
    return;
end

%% 7b — stimulus audio (5 kHz tone gated to the stim-active frames)
audPath = '';
if ADD_STIM_AUDIO
    vidDur = nFrames/outFPS;                       % seconds of output video
    nAud   = max(1, round(vidDur*AUDIO_SR));
    tt     = (0:nAud-1)'/AUDIO_SR;
    tone   = sin(2*pi*AUDIO_FREQ*tt);
    gate   = false(nAud,1);
    for ff = 1:nFrames                             % frame-accurate sync to the visual marker
        tc = t_raw(f0+ff-1);
        if any(tc >= stim_on & tc < stim_on + STIM_DUR)
            a = floor((ff-1)/outFPS*AUDIO_SR)+1;
            b = min(nAud, floor(ff/outFPS*AUDIO_SR));
            if b >= a, gate(a:b) = true; end
        end
    end
    g  = double(gate);                             % ~4 ms raised-cosine ramps to avoid clicks
    rN = round(0.004*AUDIO_SR);
    if rN >= 1, g = min(1, max(0, conv(g, ones(rN,1)/rN, 'same'))); end
    aud = AUDIO_GAIN * tone .* g;
    audPath = fullfile(fovPath, 'io_repr_stim_audio.wav');
    audiowrite(audPath, aud, AUDIO_SR);
    fprintf('stim audio: %d-pulse %g Hz tone, %.1fs -> %s\n', N_STIM, AUDIO_FREQ, vidDur, audPath);
end

%% 8 — write video (one per smoothing sigma)
for sigi = 1:numel(CA_SMOOTH_SIGMAS)
    caSig = CA_SMOOTH_SIGMAS(sigi);
    if caSig > 0, sufx = sprintf('_gauss%.2gpx', caSig); else, sufx = '_raw'; end

    modeTag = ''; if VIDEO_ONLY, modeTag = '_videoonly'; end
    if isempty(outFile)
        outFileV = fullfile(fovPath, sprintf('io_repr_t%d-%ds%s%s.mp4', round(t_start), round(t_end), modeTag, sufx));
    else
        [od,on] = fileparts(outFile); outFileV = fullfile(od, [on modeTag sufx '.mp4']);
    end

    % lossless intermediate: write PNG frames, then encode once with x264
    % (avoids the old Motion-JPEG AVI double-compression that softened the video)
    framesDir = strrep(outFileV,'.mp4','_frames');
    if isfolder(framesDir), rmdir(framesDir,'s'); end
    mkdir(framesDir);
    fprintf('Writing (%s): %s\n', sufx, outFileV);
    tic;
    for ff = 1:nFrames
        fr = ioFrame(ff, total_H,total_W, caTileW,caTileH, GAP, traceBase, traceW, ...
            caBuf, caLo, caHi, gamma_ca, caSig, CA_CMAP_INVERT, ...
            t_raw, f0, popTrace, yl, TRACE_WIN, plotL, plotR, yTop, yBot, xAxisY, ...
            stim_on, STIM_DUR, STIM_COLOR, PLAYHEAD_COLOR, TRACE_COLOR, FONT_NAME, FONT_SIZE, ST, ...
            VIDEO_ONLY, STIM_BAR_H, STIM_MARKER, STIM_DOT_R, CA_RESIZE_METHOD);
        imwrite(fr, fullfile(framesDir, sprintf('f%06d.png', ff)));
        if mod(ff,200)==0, fprintf('  %d/%d (%.1fs)\n', ff, nFrames, toc); end
    end
    fprintf('  encoding mp4 (lossless PNG -> x264 crf %d, preset %s)...\n', OUT_CRF, OUT_PRESET);
    inPat = fullfile(framesDir, 'f%06d.png');
    if ADD_STIM_AUDIO && ~isempty(audPath) && isfile(audPath)
        cmd = sprintf(['"%s" -y -framerate %d -i "%s" -i "%s" -c:v libx264 -crf %d -preset %s ' ...
                       '-pix_fmt yuv420p -c:a aac -b:a 192k -shortest "%s"'], ...
                       FFMPEG, outFPS, inPat, audPath, OUT_CRF, OUT_PRESET, outFileV);
    else
        cmd = sprintf(['"%s" -y -framerate %d -i "%s" -c:v libx264 -crf %d -preset %s ' ...
                       '-pix_fmt yuv420p "%s"'], FFMPEG, outFPS, inPat, OUT_CRF, OUT_PRESET, outFileV);
    end
    [st,out] = system(cmd);
    if st~=0, warning('ffmpeg failed:\n%s', out); else, rmdir(framesDir,'s'); end
    fprintf('  Done: %s\n', outFileV);
end
if ADD_STIM_AUDIO && ~isempty(audPath) && isfile(audPath), delete(audPath); end
fprintf('All %d version(s) written.\n', numel(CA_SMOOTH_SIGMAS));

%% ============================ LOCAL FUNCTIONS ============================
function fr = ioFrame(ff, total_H,total_W, caTileW,caTileH, GAP, traceBase, traceW, ...
    caBuf, caLo, caHi, gamma_ca, caSig, cmapInv, ...
    t_raw, f0, popTrace, yl, TRACE_WIN, plotL, plotR, yTop, yBot, xAxisY, ...
    stim_on, stimDur, stimCol, playCol, traceCol, fontName, fontSize, ST, ...
    videoOnly, stimBarH, stimMarker, stimDotR, caResize)

    fr = zeros(total_H, total_W, 3, 'uint8');

    % ---- Ca crop tile ----
    p = double(caBuf(:,:,ff));
    if caSig > 0, p = imgaussfilt(p, caSig); end
    pn = max(0,min(1,(p-caLo)/(caHi-caLo))).^gamma_ca;
    if cmapInv, pn = 1-pn; end
    if size(pn,1)==caTileH && size(pn,2)==caTileW       % native res: no resample
        u8 = uint8(255*pn);
    else
        u8 = uint8(255*imresize(pn,[caTileH caTileW], caResize));
    end
    fr(1:caTileH, 1:caTileW, :) = repmat(u8,[1 1 3]);

    % ---- VIDEO-ONLY: small stimulus indicator, no trace panel ----
    if videoOnly
        tc = t_raw(f0+ff-1);
        if any(tc >= stim_on & tc < stim_on + stimDur)     % stim active this frame
            switch lower(stimMarker)
                case 'bar'
                    h = max(1, min(stimBarH, caTileH));
                    for ch = 1:3, fr(1:h, 1:caTileW, ch) = stimCol(ch); end
                case 'dot'
                    r = stimDotR; cxd = caTileW-r-8; cyd = r+8;       % top-right corner
                    yy = max(1,cyd-r):min(caTileH,cyd+r);
                    xx = max(1,cxd-r):min(caTileW,cxd+r);
                    [Xg,Yg] = meshgrid(xx,yy); m = (Xg-cxd).^2 + (Yg-cyd).^2 <= r^2;
                    for ch = 1:3
                        sub = fr(yy,xx,ch); sub(m) = stimCol(ch); fr(yy,xx,ch) = sub;
                    end
                otherwise   % 'none'
            end
        end
        return;
    end

    % ---- trace panel ----
    panel = traceBase;
    tc = t_raw(f0+ff-1);
    x1 = tc - TRACE_WIN/2; x2 = tc + TRACE_WIN/2;
    xpix = @(tt) plotL + round((tt - x1)/max(eps,x2-x1)*(plotR-plotL));

    % stim bands (under trace)
    for s = stim_on(:)'
        a = max(s, x1); b = min(s+stimDur, x2);
        if b <= a, continue; end
        xa = max(plotL,min(plotR,xpix(a))); xb = max(plotL,min(plotR,xpix(b)));
        if xb<xa, [xa,xb]=deal(xb,xa); end
        xb = max(xb, xa+1);                              % ensure visible
        reg = panel(yTop:yBot, xa:xb, :);
        blend = 0.45*reshape(stimCol,1,1,3) + 0.55*double(reg);
        panel(yTop:yBot, xa:xb, :) = uint8(blend);
    end

    % trace line(s) in window
    idx = find(t_raw >= x1 & t_raw <= x2);
    if numel(idx) >= 2
        xp = plotL + round((t_raw(idx)-x1)/max(eps,x2-x1)*(plotR-plotL));
        if ST.on
            for r = 1:size(ST.stackMat,2)            % stacked offset traces
                yp = round(ST.laneBaseY(r) - ST.pixPerUnit*ST.stackMat(idx,r));
                panel = drawPolylineRGB(panel, xp, yp, round(traceCol), ST.lineW);
            end
        else
            yp = yBot - round((popTrace(idx)-yl(1))/max(eps,yl(2)-yl(1))*(yBot-yTop));
            panel = drawPolylineRGB(panel, xp, yp, round(traceCol), 2);
        end
    end

    % playhead
    cx = max(plotL,min(plotR,xpix(tc)));
    panel = drawLineRGB(panel, cx, yTop, cx, yBot, playCol, 2);

    % x ticks (integer raw seconds), centered labels
    for xt = ceil(x1):floor(x2)
        xp = max(plotL,min(plotR,xpix(xt)));
        panel = drawLineRGB(panel, xp, xAxisY, xp, xAxisY+6, [255 255 255], 1);
        panel = drawTextRGBCenter(panel, xp, xAxisY+8, sprintf('%d',xt), [255 255 255], fontSize, fontName);
    end

    fr(:, caTileW+GAP+1:caTileW+GAP+traceW, :) = panel;
end

function CA_CROP = selectCropRect(img)
    img = double(img); lo = prctile(img(:),1); hi = prctile(img(:),99.9);
    if hi<=lo, hi=lo+1; end
    hf = figure('Color','w'); imagesc(img,[lo hi]); axis image; colormap gray;
    title('Draw a rectangle crop for the Ca movie, then double-click / Enter');
    h = drawrectangle('Color','r'); wait(h);
    pos = round(h.Position);
    CA_CROP = [max(1,pos(1)), max(1,pos(2)), max(1,pos(3)), max(1,pos(4))];
    close(hf);
    fprintf('  Ca crop: x=%d y=%d w=%d h=%d\n', CA_CROP(1),CA_CROP(2),CA_CROP(3),CA_CROP(4));
end

function avgP = ca_avgproj(tifP, tifInfo, nDrop, nSample)
    nDir = numel(tifInfo); lo = nDrop+1; hi = nDir;
    nS = max(1, min(nSample, hi-lo+1)); idx = unique(round(linspace(lo,hi,nS)));
    avgP = zeros(tifInfo(1).Height, tifInfo(1).Width, 'single');
    for i=1:numel(idx), avgP = avgP + single(imread(tifP, idx(i), 'Info', tifInfo)); end
    avgP = avgP/numel(idx);
end

function tifPath = discover_mc_tiff(fovFolder)
    mc = dir(fullfile(fovFolder,'*_MC*.tif'));
    assert(~isempty(mc),'No *_MC*.tif in %s',fovFolder);
    nm = {mc.name}; keep = true(size(nm));
    for ii=1:numel(nm)
        if contains(nm{ii},'_AVG','IgnoreCase',true) || contains(nm{ii},'_QC','IgnoreCase',true) || ...
           contains(nm{ii},'_shifts','IgnoreCase',true) || contains(nm{ii},'_SELROI','IgnoreCase',true) || ...
           contains(nm{ii},'_ROI','IgnoreCase',true) || contains(nm{ii},'_masks','IgnoreCase',true)
            keep(ii) = false;
        end
    end
    mc = mc(keep); assert(~isempty(mc));
    nMC = arrayfun(@(x) numel(strfind(upper(x.name),'_MC')), mc);
    [~,b] = max(nMC); tifPath = fullfile(mc(b).folder, mc(b).name);
end

function img = drawPolylineRGB(img, x, y, col, lw)
    x = round(x(:)); y = round(y(:)); g = isfinite(x)&isfinite(y); x=x(g); y=y(g);
    for i=1:numel(x)-1, img = drawLineRGB(img, x(i),y(i),x(i+1),y(i+1), col, lw); end
end

function img = drawLineRGB(img, x1,y1,x2,y2, col, lw)
    [H,W,~]=size(img); x1=round(x1);y1=round(y1);x2=round(x2);y2=round(y2);
    n=max(abs(x2-x1),abs(y2-y1))+1; xs=round(linspace(x1,x2,n)); ys=round(linspace(y1,y2,n)); rad=floor(lw/2);
    for q=1:numel(xs)
        xx=xs(q); yy=ys(q); if xx<1||xx>W||yy<1||yy>H, continue; end
        xr=max(1,xx-rad):min(W,xx+rad); yr=max(1,yy-rad):min(H,yy+rad);
        img(yr,xr,1)=col(1); img(yr,xr,2)=col(2); img(yr,xr,3)=col(3);
    end
end

function img = drawTextRGB(img, x, y, str, col, fontSize, fontName)
    if nargin<7||isempty(fontName), fontName='Arial'; end
    if nargin<6||isempty(fontSize), fontSize=14; end
    try
        img = insertText(img,[x y],str,'TextColor',col,'BoxOpacity',0,'FontSize',fontSize,'Font',fontName);
    catch
        try, img = insertText(img,[x y],str,'TextColor',col,'BoxOpacity',0,'FontSize',fontSize); catch, end
    end
end

function img = drawTextRGBCenter(img, xCenter, yTop, str, col, fontSize, fontName)
    if nargin<7||isempty(fontName), fontName='Arial'; end
    if nargin<6||isempty(fontSize), fontSize=14; end
    persistent cache; if isempty(cache), cache=containers.Map('KeyType','char','ValueType','any'); end
    key = sprintf('%s|%d|%s', str, fontSize, fontName);
    if isKey(cache,key), wl=cache(key); w=wl(1); leftPad=wl(2);
    else
        tmp=zeros(90,500,3,'uint8'); ok=true;
        try, tmp=insertText(tmp,[2 2],str,'TextColor',[255 255 255],'BoxOpacity',0,'FontSize',fontSize,'Font',fontName);
        catch, try, tmp=insertText(tmp,[2 2],str,'TextColor',[255 255 255],'BoxOpacity',0,'FontSize',fontSize); catch, ok=false; end, end
        if ~ok, img=drawTextRGB(img,round(xCenter-0.28*fontSize*numel(str)),yTop,str,col,fontSize,fontName); return; end
        cc=find(any(any(tmp>0,3),1)); if isempty(cc), return; end
        w=max(cc)-min(cc)+1; leftPad=min(cc)-2; cache(key)=[w leftPad];
    end
    img = drawTextRGB(img, round(xCenter-w/2)-leftPad, yTop, str, col, fontSize, fontName);
end

function img = drawTextRGBRight(img, xRight, y, str, col, fontSize, fontName)
    if nargin<7||isempty(fontName), fontName='Arial'; end
    if nargin<6||isempty(fontSize), fontSize=14; end
    tmp=zeros(90,400,3,'uint8'); ok=true;
    try, tmp=insertText(tmp,[2 2],str,'TextColor',col,'BoxOpacity',0,'FontSize',fontSize,'Font',fontName);
    catch, try, tmp=insertText(tmp,[2 2],str,'TextColor',col,'BoxOpacity',0,'FontSize',fontSize); catch, ok=false; end, end
    if ~ok, return; end
    m=any(tmp>0,3); [rr,cc]=find(m); if isempty(rr), return; end
    r1=min(rr);r2=max(rr);c1=min(cc);c2=max(cc); crop=tmp(r1:r2,c1:c2,:); mk=m(r1:r2,c1:c2);
    [H,W,~]=size(img); [h,w,~]=size(crop); x=round(xRight-w+1); y=round(y);
    if x>W||y>H||x+w-1<1||y+h-1<1, return; end
    x1=max(1,x);y1=max(1,y);x2=min(W,x+w-1);y2=min(H,y+h-1);
    cx1=x1-x+1;cy1=y1-y+1;cx2=cx1+(x2-x1);cy2=cy1+(y2-y1);
    crop2=crop(cy1:cy2,cx1:cx2,:); mk2=mk(cy1:cy2,cx1:cx2);
    for ch=1:3, reg=img(y1:y2,x1:x2,ch); t=crop2(:,:,ch); reg(mk2)=t(mk2); img(y1:y2,x1:x2,ch)=reg; end
end
