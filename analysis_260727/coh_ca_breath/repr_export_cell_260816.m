% repr_export_cell_260816.m
% -----------------------------------------------------------------------
%  Representative composite MP4 -- calcium patches + breathing camera +
%  scrolling traces -- selected and saved the way the per-cell summary GUI does.
%
%  WHAT IS DIFFERENT from repr_export_RZ_260811.m (repo root), which this script
%  is otherwise a copy of:
%
%    SELECTION.  You type ids and nothing else. An id is read as a registry CELL
%    id first -- the same number as in the GUI's id box and in the batch filenames
%    (per-cell-summary_..._cell%03d.png) -- and THE RECORDING FOLLOWS FROM IT: the
%    registry already knows which recordings hold that cell, so fovPath is only a
%    preference. Failing that the id is read as a ROI index in fovPath, and if
%    that ROI is a registered cell you get the cell, labelled with its cell id.
%    So in Vgat/0730/deep, ids = 279 and ids = 29 are the same neuron. Plain
%    unregistered ROIs still work -- most ROIs are not registered cells.
%    Identity comes from
%      1. the registry (event_latency_data.mat + cell_merge_overrides_260814)
%      2. the session's cell_link.mat
%      3. plain ROI indices
%    in that order, exactly as temporal_phase_cell_gui_260812.m resolves it. The
%    registry is AUTHORITATIVE for any folder it covers: cell_link numbers cells
%    per session and the registry numbers them globally, so falling through to
%    cell_link on a miss would silently export a real but DIFFERENT cell.
%    idMode = 'cell' / 'roi' force one reading of an ambiguous number.
%
%    SAVING.  Everything lands in ONE collected outDir, named
%    <recName>_cell029-031_t35-65s.mp4, instead of being written beside whichever
%    recording happened to be open. Same reasoning as the GUI's gui_renders
%    folder: the filename already carries the recording, the ids and the window,
%    so a collected folder cannot collide and does not scatter renders through
%    the data archive.
%
%    Y-LIMIT OVERRIDES are keyed by the SELECTION ID (cell id in cell mode, ROI
%    index in roi mode), not by row position, so they stay attached to the right
%    neuron when the id list changes. The chosen limits are printed every run --
%    copy a line into YLIM_BY_ID and tweak it.
%
%  ONE RECORDING.  Unlike the summary figure, a video cannot pool: every tile
%  plays the same movie, so all the ids have to land in one recording. The script
%  picks one that holds them all -- preferring fovPath when it qualifies, and
%  saying so when it moves elsewhere. Only when NO single recording holds every id
%  does it stop, listing where each one lives.
%
%  TRAP: a passing preview does not mean the export passes. MAKE_LAYOUT_PREVIEW_ONLY
%  returns at the end of section 7c, so the naming, VideoWriter loop and ffmpeg
%  transcode in section 8 are never exercised by a preview run.
%
%  Requires in fovPath: *_cpSAM_output.mat, *_MC*.tif, a behavior AVI.
%
%  Runqi Zhang / 2026-08-16
% -----------------------------------------------------------------------

clear; clc; close all;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(fullfile(repoRoot, 'analysis_260806'));   % registry merges + pooled_obs

%% ============================== USER PARAMS ==============================
% WHICH CELLS. Just type the ids. Tile stack = 1 breath tile + numel(ids) Ca tiles.
%
% In the default 'auto' mode an id is whatever it can be, in this order:
%   1. a registry CELL id  -- the number on the summary figures (..._cell%03d.png)
%      and in the GUI's id box. The recording is then found FOR you.
%   2. a ROI index in fovPath that belongs to a registered cell -- you get that
%      cell, labelled with its cell id. (In Vgat/0730/deep, "29" is ROI 29, which
%      is cell 279; both spellings land on the same neuron.)
%   3. a plain ROI index in fovPath -- most ROIs are not registered cells, and
%      those still work, just labelled "ROI n".
% idMode = 'cell' or 'roi' forces one reading when a number could be both.
ids           = [53];
idMode        = 'cell';        % 'auto' | 'cell' | 'roi'

% Which recording. '' = let the ids choose it. When set it is a PREFERENCE, not a
% constraint: if the ids are not in this folder the script moves to the recording
% that holds them and says so, like the GUI does. A video cannot pool recordings
% -- every tile plays the same movie -- so all ids must end up in ONE folder.
fovPath       = '';

% Cell identity source. '' forces the old cell_link-then-ROI behaviour.
regFile       = 'D:\Ventral_surface_summary\event_latency_260811\event_latency_data.mat';

% Collected output folder for every clip this script renders. '' = write into
% fovPath instead.
outDir        = 'D:\Ventral_surface_summary\repr_renders_260816';

t_start       = 50;            % clip start (s)
t_end         = 80;            % clip end (s)

patchHalf_um  = 25;            % Ca ROI crop half-size (um)

BREATH_CROP   = [];            % [] = interactive select; or [x1 y1 ss]

% Behavioral video display contrast, applied AFTER crop
BEH_PRC_LO          = 5;       % lower percentile after crop
BEH_PRC_HI          = 99;    % upper percentile after crop
BEH_GAMMA           = 1;    % <1 brightens midtones; 1 = linear
BEH_CLIP_PER_FRAME  = false;   % false = stable sampled crop CLIM; true = reclip every frame
BEH_CLIM            = [];      % [] = auto estimate from cropped behavior region

% Breathing SVD options
BREATH_SVD_USE_CROP = true;    % true = compute SVD only from BREATH_CROP
BREATH_SVD_SIGN     = 1;       % use -1 to flip chest SVD trace sign

outFPS        = 30;            % output video fps
gamma_ca      = 1;             % Ca patch display gamma
CA_PRC_LO     = 25;            % Ca display clim lower percentile
CA_PRC_HI     = 99.5;          % Ca display clim upper percentile
TRACE_WIN     = 10.0;          % rolling trace window (s)
outFile       = '';            % '' = auto-named in outDir

%% --------------------------- LAYOUT PARAMETERS ---------------------------
% Main panel geometry
GAP            = 5;       % gap between left video panel and right trace panel
TILE_SZ        = 220;     % tile size, px
TILE_GAP       = 6;       % gap between video/trace rows, px
LEFT_LABEL_W   = 95;      % label gutter left of video tiles, px

% Right trace panel uniform scale + extra right background
TRACE_PANEL_SCALE   = 1.00;     % enlarge/shrink right trace plot uniformly
TRACE_PLOT_W_BASE   = 720;      % active plotting width before extra background
TRACE_PLOT_W        = round(TRACE_PLOT_W_BASE * TRACE_PANEL_SCALE);

TRACE_EXTRA_RIGHT_BG = 45;      % extra black space at right to prevent tick-label clipping
TRACE_PANEL_W        = TRACE_PLOT_W + TRACE_EXTRA_RIGHT_BG;

% Bottom x-axis reserved area
X_AXIS_H       = 72;      % bottom reserved space for x-axis tick labels + xlabel

% Left video label placement
LEFT_TILE_LABEL_X      = 6;
LEFT_TILE_LABEL_YSHIFT = -4;

% Right trace label / axis placement
% NOTE: these are computed relative to TRACE_PLOT_W, not the extra right background.
TRACE_AXIS_X_FRAC       = 0.170;   % y-axis / trace start x-position
TRACE_RIGHT_PAD_FRAC    = 0.020;   % right margin of active trace plot

TRACE_LABEL_RIGHT_DX    = 18;      % right edge of "Chest motion"/"dF/F" label is this px left of y-axis
TRACE_LABEL_YSHIFT      = -4;      % vertical shift for right trace labels

TRACE_YNUM_RIGHT_DX     = 8;       % right edge of y-number labels is this px left of y-axis
TRACE_YNUM_TOP_YSHIFT   = -2;
TRACE_YNUM_BOT_YSHIFT   = -18;

% Internal padding inside each trace row
TRACE_ROW_TOP_PAD       = 16;      % prevents top clipping
TRACE_ROW_BOT_PAD       = 12;      % prevents bottom clipping

% Trace line / axis line widths
TRACE_LINE_W            = 2;
TRACE_PLAYHEAD_W        = 2;
TRACE_AXIS_W            = 1;

% Bottom time axis placement
X_AXIS_Y_FRAC           = 0.20;    % location inside X_AXIS_H
X_TICK_LABEL_YOFF       = 8;
X_LABEL_YOFF            = 34;      % lower = closer to axis; higher = further down

% Font
FONT_NAME               = 'Arial';
FONT_SIZE_LABEL         = 18;      % left tile labels
FONT_SIZE_AX            = 18;      % y numeric ticks + x numeric ticks
FONT_SIZE_XLAB          = 18;      % Time (s)

% Right trace font scales with TRACE_PANEL_SCALE
FONT_SIZE_TRACE_LABEL   = round(FONT_SIZE_LABEL * TRACE_PANEL_SCALE);
FONT_SIZE_TRACE_AX      = round(FONT_SIZE_AX    * TRACE_PANEL_SCALE);
FONT_SIZE_TRACE_XLAB    = round(FONT_SIZE_XLAB  * TRACE_PANEL_SCALE);

% Y-limits.
% Auto by default. Overrides are keyed by SELECTION ID -- the cell id in cell
% mode, the ROI index in roi mode -- not by row position, so they stay attached to
% the right neuron when ids changes. The chosen limits are printed on every run.
YLIM_AUTO_MODE = 'robust';  % 'robust' = percentile floor + spike-proof peak
                            % 'minmax' = raw min/max (old behaviour)

YLIM_BREATH = [0 1.2];      % chest-motion row; [NaN NaN] = auto

YLIM_BY_ID = [              % [id  lo  hi];  ids not listed stay auto
    4   -0.2   0.4
];

YLIM_PAD_FRAC = 0.10;       % auto ylim padding, as a fraction of the range

% Preview workflow
MAKE_LAYOUT_PREVIEW_ONLY = false;  % true = show one frame and stop
PREVIEW_T_SEC            = 30;     % preview frame time, seconds
%% =========================================================================

%% 0 — Cell / ROI selection
% Resolve the ids -- and the recording they live in -- before anything is loaded,
% so a wrong id fails immediately rather than after the TIFF read.
fprintf('Resolving ids %s (%s mode)\n', mat2str(ids), lower(idMode));

[fovPath, roiIdxs, selIds, idLabels, idTag, idSrc, recName] = ...
        resolve_ids(fovPath, ids, idMode, regFile);

fprintf('  recording      : %s\n', fovPath);
fprintf('  identity source: %s\n', idSrc);
for k = 1:numel(roiIdxs)
    fprintf('  %-10s -> ROI %d\n', idLabels{k}, roiIdxs(k));
end

%% 1 — FPS + pixel size
fprintf('Detecting FPS...\n');
try
    [fps_ca, sc_meta] = detect_session_fps(fovPath);
catch
    fps_ca = 30; sc_meta = struct('pixelSize_um', 1.7778);
end

if ~isfield(sc_meta,'pixelSize_um') || ~isfinite(sc_meta.pixelSize_um)
    sc_meta.pixelSize_um = 1.7778;
end

pixSz        = sc_meta.pixelSize_um;
patchHalf_px = round(patchHalf_um / pixSz);

fprintf('  fps=%d Hz  %.4f um/px  patch_half=%d px\n', fps_ca, pixSz, patchHalf_px);

%% 2 — cpSAM + dF/F
fprintf('Loading cpSAM...\n');

sam_h = dir(fullfile(fovPath,'*_cpSAM_output.mat'));
assert(~isempty(sam_h), 'No *_cpSAM_output.mat in %s', fovPath);

sd    = load(fullfile(sam_h(1).folder, sam_h(1).name), 'maskL','F');
maskL = sd.maskL;
F_all = double(sd.F);

[T_F, N_roi] = size(F_all);
fprintf('  F: [%d x %d]\n', T_F, N_roi);

assert(max(roiIdxs) <= N_roi, 'ROI index %d > N_roi=%d', max(roiIdxs), N_roi);

fprintf('Computing dF/F...\n');

dout  = helper.dFF_RZ(F_all, 'FPS', fps_ca);
dFF   = dout.dFF;

T_dff = size(dFF,1);
nDrop = T_F - T_dff;
t_ca  = (0:T_dff-1)' / fps_ca;

fprintf('  T_dff=%d  nDrop=%d  duration=%.1fs\n', T_dff, nDrop, t_ca(end));

%% 3 — Clip frame range
f0      = max(1, round(t_start*fps_ca)+1);
f1      = min(T_dff, round(t_end*fps_ca));

assert(f1 > f0, 'Empty clip window.');

nFrames = f1 - f0 + 1;

fprintf('Clip: t=%.1f–%.1fs  frames %d–%d  (%d total)\n', ...
    t_start, t_end, f0, f1, nFrames);

%% 4 — Load Ca TIFF patches
fprintf('Loading TIFF patches...\n');

tifP    = discover_mc_tiff(fovPath);
tifInfo = imfinfo(tifP);

imgH = tifInfo(1).Height;
imgW = tifInfo(1).Width;

props = regionprops(maskL,'Centroid');
assert(numel(props) >= max(roiIdxs));

nCa    = numel(roiIdxs);   % number of Ca ROI tiles
nTiles = nCa + 1;          % + the breath camera tile on top

r1v=zeros(nCa,1); r2v=zeros(nCa,1);
c1v=zeros(nCa,1); c2v=zeros(nCa,1);

for k = 1:nCa
    ri = roiIdxs(k);
    cx = round(props(ri).Centroid(1));
    cy = round(props(ri).Centroid(2));

    r1v(k) = max(1, cy - patchHalf_px);
    r2v(k) = min(imgH, cy + patchHalf_px);
    c1v(k) = max(1, cx - patchHalf_px);
    c2v(k) = min(imgW, cx + patchHalf_px);
end

ca_buf = cell(nCa,1);
for k = 1:nCa
    ca_buf{k} = zeros(r2v(k)-r1v(k)+1, c2v(k)-c1v(k)+1, nFrames, 'single');
end

tic;
for ff = 1:nFrames
    frm = single(imread(tifP, nDrop+f0+ff-1, 'Info', tifInfo));

    for k = 1:nCa
        ca_buf{k}(:,:,ff) = frm(r1v(k):r2v(k), c1v(k):c2v(k));
    end

    if mod(ff,200)==0
        fprintf('  %d/%d (%.1fs)\n', ff, nFrames, toc);
    end
end

fprintf('  Loaded in %.1fs\n', toc);

ca_lo = zeros(nCa,1);
ca_hi = zeros(nCa,1);
ca_perim_up = cell(nCa,1);

for k = 1:nCa
    ri = roiIdxs(k);

    msk = (maskL(r1v(k):r2v(k), c1v(k):c2v(k)) == ri);
    pxv = ca_buf{k}(repmat(msk,[1 1 nFrames]));

    lo = double(prctile(pxv, CA_PRC_LO));
    hi = double(prctile(pxv, CA_PRC_HI));

    if hi <= lo
        hi = lo + 1;
    end

    ca_lo(k) = lo;
    ca_hi(k) = hi;

    ca_perim_up{k} = imresize(double(bwperim(msk,4)), [TILE_SZ TILE_SZ], 'nearest') > 0.5;
end

%% 5 — Breathing camera
fprintf('Setting up breathing camera...\n');

% Behavior AVI naming differs across rigs: older sessions are Basler_*.avi,
% the current Basler GUI writes cam1_*/cam2_*.avi. Try each in turn, then fall
% back to any AVI. Never pick up our own transcode or a DLC-labeled video.
AVI_PATTERNS = {'Basler_*.avi', 'cam*.avi', '*.avi'};

avi_h = [];
for ip = 1:numel(AVI_PATTERNS)
    hits = dir(fullfile(fovPath, AVI_PATTERNS{ip}));
    if ~isempty(hits)
        hits = hits(~contains({hits.name}, {'_mjpeg.avi','_labeled.avi'}));
    end
    if ~isempty(hits)
        avi_h = hits;
        break;
    end
end

assert(~isempty(avi_h), ...
    'No behavior AVI (Basler_*.avi / cam*.avi) in %s', fovPath);

if numel(avi_h) > 1
    [~, ibig] = max([avi_h.bytes]);     % largest = the actual recording
    fprintf('  %d AVIs found, using largest: %s\n', numel(avi_h), avi_h(ibig).name);
    avi_h = avi_h(ibig);
end

aviPath  = fullfile(avi_h(1).folder, avi_h(1).name);
fprintf('  Behavior AVI: %s\n', avi_h(1).name);

mjpgPath = [aviPath(1:end-4) '_mjpeg.avi'];

if ~isfile(mjpgPath)
    fprintf('  Transcoding FFV1 → MJPEG AVI (one-time)...\n');

    % Stamp the TRUE frame rate into the transcode. The camera is
    % 2P-frame-triggered, so its real rate is the calcium rate; the AVI
    % header's declared rate is unreliable (300 fps seen on 260806_sst).
    % -r goes BEFORE -i so it reinterprets the input timebase and frames pass
    % through 1:1. On the output side it would resample and drop frames.
    cmd = sprintf('"%s" -y -r %g -i "%s" -c:v mjpeg -q:v 3 "%s"', ...
                  resolve_ffmpeg(), fps_ca, aviPath, mjpgPath);
    [st, out] = system(cmd);

    if st ~= 0 || ~isfile(mjpgPath)
        error('ffmpeg failed.\nCmd: %s\nOut: %s', cmd, out);
    end

    fprintf('  Done: %s\n', mjpgPath);
end

vBeh    = VideoReader(mjpgPath);
fps_beh = vBeh.FrameRate;
nBeh    = round(vBeh.Duration * fps_beh);

toss_beh = nBeh - T_dff;

fprintf('  Breath cam: %d frames @ %.1f fps  (toss_beh=%d to match dFF)\n', ...
    nBeh, fps_beh, toss_beh);

if abs(fps_beh - fps_ca) / fps_ca > 0.1
    fprintf(['  NOTE: AVI header fps (%.1f) disagrees with calcium fps (%.1f) — ' ...
             'header is unreliable.\n        Frames are 1:1 (2P-triggered), so ' ...
             'seeking is done by frame index in the file''s own timebase.\n'], ...
            fps_beh, fps_ca);
end

if isempty(BREATH_CROP)
    BREATH_CROP = selectBreathCropInteractive(mjpgPath);
end

fprintf('  Breath crop: x=%d y=%d sz=%d\n', ...
    BREATH_CROP(1), BREATH_CROP(2), BREATH_CROP(3));

% Estimate behavior display contrast from the cropped region, not full frame
if isempty(BEH_CLIM) && ~BEH_CLIP_PER_FRAME
    BEH_CLIM = estimateBreathClipLimits(mjpgPath, BREATH_CROP, BEH_PRC_LO, BEH_PRC_HI);
    fprintf('  Behavior display CLIM from crop: [%.2f %.2f]\n', BEH_CLIM(1), BEH_CLIM(2));
end

%% 6 — SVD breathing trace
breath = [];

% Use a crop-specific cache file so old full-frame breath_pc1.mat is not reused.
if BREATH_SVD_USE_CROP
    pc1Name = sprintf('breath_pc1_crop_x%d_y%d_s%d.mat', ...
        round(BREATH_CROP(1)), round(BREATH_CROP(2)), round(BREATH_CROP(3)));
else
    pc1Name = 'breath_pc1_fullframe.mat';
end

pc1Path = fullfile(fovPath, pc1Name);

if isfile(pc1Path)
    try
        tmp = load(pc1Path, 'breathTrace');
        br  = tmp.breathTrace(:);

        if numel(br) ~= T_dff
            br = interp1(linspace(0,1,numel(br)), br, linspace(0,1,T_dff), 'pchip');
        end

        breath = br(:);
        fprintf('SVD trace loaded from %s\n', pc1Name);

    catch ME
        fprintf('%s failed: %s\n', pc1Name, ME.message);
    end
end

if isempty(breath)
    if BREATH_SVD_USE_CROP
        fprintf('Computing SVD breath trace from cropped region [x=%d y=%d s=%d]...\n', ...
            round(BREATH_CROP(1)), round(BREATH_CROP(2)), round(BREATH_CROP(3)));
    else
        fprintf('Computing SVD breath trace from full original FFV1 AVI...\n');
    end

    svdPy  = fullfile(fovPath, 'tmp_svd_breath.py');
    svdTxt = fullfile(fovPath, 'tmp_svd_breath.txt');

    pyLines = {
        'import cv2, numpy as np, sys'
        'avi_path = sys.argv[1]'
        'out_txt  = sys.argv[2]'
        'BIN      = int(sys.argv[3])'
        'USE_CROP = int(sys.argv[4])'
        'X1       = int(sys.argv[5])'
        'Y1       = int(sys.argv[6])'
        'SS       = int(sys.argv[7])'
        ''
        'cap = cv2.VideoCapture(avi_path)'
        'frames = []'
        'hb = wb = 1'
        ''
        'while True:'
        '    ret, fr = cap.read()'
        '    if not ret: break'
        ''
        '    if fr.ndim == 3:'
        '        fr = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY).astype(np.float32)'
        '    else:'
        '        fr = fr.astype(np.float32)'
        ''
        '    h, w = fr.shape'
        ''
        '    if USE_CROP:'
        '        # MATLAB crop is 1-based [x y side]; Python is 0-based'
        '        x0 = max(0, min(w-1, X1-1))'
        '        y0 = max(0, min(h-1, Y1-1))'
        '        ss = max(1, SS)'
        '        x2 = min(w, x0 + ss)'
        '        y2 = min(h, y0 + ss)'
        '        ss2 = min(x2-x0, y2-y0)'
        '        if ss2 > 1:'
        '            fr = fr[y0:y0+ss2, x0:x0+ss2]'
        ''
        '    h, w = fr.shape'
        '    hb = max(1, h//BIN)'
        '    wb = max(1, w//BIN)'
        '    fr_small = cv2.resize(fr, (wb, hb), interpolation=cv2.INTER_AREA)'
        '    frames.append(fr_small.ravel())'
        ''
        'cap.release()'
        ''
        'X = np.array(frames, dtype=np.float32).T'
        'X -= X.mean(axis=1, keepdims=True)'
        ''
        '# Power iteration for first temporal PC'
        'rng = np.random.default_rng(0)'
        'v = rng.standard_normal(X.shape[1]).astype(np.float32)'
        'for _ in range(30):'
        '    v = X.T @ (X @ v)'
        '    v /= (np.linalg.norm(v) + 1e-12)'
        ''
        '# Deterministic sign convention based on skewness'
        'mu = v.mean()'
        'sk = float(((v-mu)**3).mean() / (((v-mu)**2).mean()**1.5 + 1e-12))'
        'if sk < 0: v = -v'
        ''
        'np.savetxt(out_txt, v)'
        'print(f"SVD done: {len(frames)} frames, {hb}x{wb} binned, crop={USE_CROP}")'
    };

    fid = fopen(svdPy,'w');
    for li = 1:numel(pyLines)
        fprintf(fid,'%s\n', pyLines{li});
    end
    fclose(fid);

    pyExe = 'python';
    for c = {'C:\Program Files\Python314\python.exe', ...
             'C:\ProgramData\Miniconda3\python.exe', ...
             'C:\ProgramData\Miniconda3\envs\cellpose-gpu\python.exe', ...
             'C:\ProgramData\Miniconda3\envs\dlc310\python.exe'}
        if isfile(c{1})
            pyExe = c{1};
            break;
        end
    end

    if BREATH_SVD_USE_CROP
        useCropFlag = 1;
        cropX = round(BREATH_CROP(1));
        cropY = round(BREATH_CROP(2));
        cropS = round(BREATH_CROP(3));
    else
        useCropFlag = 0;
        cropX = 1;
        cropY = 1;
        cropS = 1;
    end

    cmd = sprintf('"%s" "%s" "%s" "%s" 8 %d %d %d %d', ...
        pyExe, svdPy, aviPath, svdTxt, useCropFlag, cropX, cropY, cropS);

    fprintf('  Running: %s\n', cmd);

    [st, out] = system(cmd);
    fprintf('  Python output: %s\n', strtrim(out));

    if st == 0 && isfile(svdTxt)
        br = readmatrix(svdTxt);
        br = br(:);

        if toss_beh > 0 && toss_beh < numel(br)
            br = br(toss_beh+1:end);
        end

        if numel(br) ~= T_dff
            br = interp1(linspace(0,1,numel(br)), br, linspace(0,1,T_dff), 'pchip');
        end

        breath = br(:);
        breathTrace = breath;

        save(pc1Path, 'breathTrace', 'BREATH_CROP', 'BREATH_SVD_USE_CROP');
        fprintf('  Saved SVD trace cache: %s\n', pc1Name);
        fprintf('  SVD breath trace: %d samples\n', numel(breath));
    else
        warning('Python SVD failed — breath trace will be zeros.\n%s', out);
        breath = zeros(T_dff,1);
    end

    if isfile(svdPy),  delete(svdPy);  end
    if isfile(svdTxt), delete(svdTxt); end
end

% Optional user-controlled sign flip
if BREATH_SVD_SIGN ~= 1
    breath = BREATH_SVD_SIGN * breath;
    fprintf('Applied BREATH_SVD_SIGN = %.1f\n', BREATH_SVD_SIGN);
end

%% 7 — Trace data prep + layout geometry

idx_clip = t_ca >= t_start & t_ca <= t_end;

% Tile labels carry the SELECTION id -- "Cell 29" in cell mode, "ROI 29" in roi
% mode -- so the video is labelled with the same number the summary figures use.
videoLabels = [{'Breath'}, idLabels(:).'];

traceLabels = [{'Chest motion'}, repmat({[char(916) 'F/F']}, 1, nCa)];

% Normalize breath to [0, 1] over the clip window
breathDisp = breath;
bclip = breathDisp(idx_clip);
bclip = bclip(isfinite(bclip));

if isempty(bclip)
    breathDisp(:) = 0;
else
    blo = prctile(bclip,1);
    bhi = prctile(bclip,99);

    if bhi <= blo
        bhi = blo + 1;
    end

    breathDisp = (breathDisp - blo) ./ (bhi - blo);
    breathDisp = max(0, min(1, breathDisp));
end

trData = [breathDisp, dFF(:,roiIdxs)];

% Auto Y limits first
trYLim = zeros(nTiles,2);
trYLim(1,:) = [0 1];

for kk = 2:nTiles
    yy = trData(idx_clip,kk);
    yy = yy(isfinite(yy));

    if isempty(yy)
        trYLim(kk,:) = [-1 1];
        continue;
    end

    switch lower(YLIM_AUTO_MODE)
        case 'robust'
            % Floor at the 1st percentile so a few deep noise dips can't drag
            % the axis down. Ceiling at the max of a 3-sample running median,
            % so one hot frame can't blow up the scale but a real transient
            % (many frames wide) still sets the top.
            lo = prctile(yy, 1);
            hi = max(movmedian(yy, 3));
            if hi <= lo
                lo = min(yy);  hi = max(yy);
            end
        otherwise
            lo = min(yy);
            hi = max(yy);
    end

    if hi <= lo
        ctr = lo;
        lo = ctr - 0.5;
        hi = ctr + 0.5;
    end

    r = hi - lo;
    trYLim(kk,:) = [lo - YLIM_PAD_FRAC*r, hi + YLIM_PAD_FRAC*r];
end

% Manual overrides: breath row, then per-cell keyed by SELECTION id.
if numel(YLIM_BREATH)==2 && all(isfinite(YLIM_BREATH)) && YLIM_BREATH(2)>YLIM_BREATH(1)
    trYLim(1,:) = YLIM_BREATH(:).';
end

isManual = false(nCa,1);
for k = 1:nCa
    if isempty(YLIM_BY_ID); break; end
    hit = find(YLIM_BY_ID(:,1) == selIds(k), 1);
    if ~isempty(hit) && YLIM_BY_ID(hit,3) > YLIM_BY_ID(hit,2)
        trYLim(k+1,:) = YLIM_BY_ID(hit,2:3);
        isManual(k)   = true;
    end
end

fprintf('Y-limits (auto mode: %s)\n', YLIM_AUTO_MODE);
fprintf('  Breath      [%8.3f %8.3f]\n', trYLim(1,1), trYLim(1,2));
for k = 1:nCa
    if isManual(k); src = 'manual'; else; src = 'auto'; end
    fprintf('  %-10s  [%8.3f %8.3f]  %s   (copy: %d %.3f %.3f)\n', ...
            idLabels{k}, trYLim(k+1,1), trYLim(k+1,2), src, ...
            selIds(k), trYLim(k+1,1), trYLim(k+1,2));
end

left_W = LEFT_LABEL_W + TILE_SZ;

tileStack_H = nTiles*TILE_SZ + (nTiles-1)*TILE_GAP;
trace_H     = tileStack_H + X_AXIS_H;
trace_W     = TRACE_PANEL_W;

total_W = left_W + GAP + trace_W;
total_H = trace_H;

total_W = total_W + mod(total_W,2);
total_H = total_H + mod(total_H,2);

fprintf('Frame size: %d x %d px\n', total_W, total_H);

% Exact row alignment shared by left video tiles and right trace bands
rowTop = zeros(nTiles,1);
rowTop(1) = 1;
for jj = 2:nTiles
    rowTop(jj) = rowTop(jj-1) + TILE_SZ + TILE_GAP;
end

tileX = LEFT_LABEL_W + 1;

% Pre-render static left-label background
leftBase = zeros(total_H, left_W, 3, 'uint8');

for jj = 1:nTiles
    yLabel = rowTop(jj) + round(TILE_SZ/2) - round(FONT_SIZE_LABEL/2) + LEFT_TILE_LABEL_YSHIFT;
    leftBase = drawTextRGB(leftBase, LEFT_TILE_LABEL_X, yLabel, videoLabels{jj}, ...
        [255 255 255], FONT_SIZE_LABEL, FONT_NAME);
end

% Pre-render static trace-panel background
traceBase = makeTraceBasePanel(trace_H, trace_W, TRACE_PLOT_W, ...
    rowTop, TILE_SZ, tileStack_H, X_AXIS_H, ...
    trYLim, traceLabels, FONT_NAME, FONT_SIZE_TRACE_LABEL, FONT_SIZE_TRACE_AX, FONT_SIZE_TRACE_XLAB, ...
    TRACE_AXIS_X_FRAC, TRACE_RIGHT_PAD_FRAC, ...
    TRACE_LABEL_RIGHT_DX, TRACE_LABEL_YSHIFT, ...
    TRACE_YNUM_RIGHT_DX, TRACE_YNUM_TOP_YSHIFT, TRACE_YNUM_BOT_YSHIFT, ...
    X_AXIS_Y_FRAC, X_LABEL_YOFF, ...
    TRACE_ROW_TOP_PAD, TRACE_ROW_BOT_PAD, TRACE_AXIS_W);

%% 7c — Layout preview frame

if MAKE_LAYOUT_PREVIEW_ONLY
    fprintf('Generating layout preview only...\n');

    ffPrev = round((PREVIEW_T_SEC - t_start) * fps_ca) + 1;
    ffPrev = max(1, min(nFrames, ffPrev));

    tc = t_ca(f0 + ffPrev - 1);

    % Seek by FRAME, in the video's OWN timebase. The behavior cam is
    % 2P-frame-triggered (1:1 with calcium), so the frame index is
    % f0+ffPrev-2+toss_beh — but CurrentTime is measured in the file's
    % declared fps, which is often bogus (300 here vs 30 for calcium).
    % Dividing by fps_ca seeks past the end of short-header files.
    vBeh.CurrentTime = (f0 + ffPrev - 2 + toss_beh) / fps_beh;

    previewFrame = makeOneFrame(ffPrev, tc, ...
        total_H, total_W, left_W, GAP, trace_W, TRACE_PLOT_W, ...
        leftBase, tileX, rowTop, TILE_SZ, ...
        vBeh, BREATH_CROP, BEH_CLIM, BEH_PRC_LO, BEH_PRC_HI, BEH_GAMMA, BEH_CLIP_PER_FRAME, ...
        ca_buf, ca_lo, ca_hi, ca_perim_up, gamma_ca, ...
        traceBase, t_ca, trData, trYLim, TRACE_WIN, ...
        tileStack_H, X_AXIS_H, FONT_NAME, FONT_SIZE_TRACE_AX, ...
        TRACE_AXIS_X_FRAC, TRACE_RIGHT_PAD_FRAC, X_AXIS_Y_FRAC, X_TICK_LABEL_YOFF, ...
        TRACE_ROW_TOP_PAD, TRACE_ROW_BOT_PAD, TRACE_LINE_W, TRACE_PLAYHEAD_W);

    figure('Color','k');
    imshow(previewFrame);
    set(gcf,'Name','Layout preview');
    return;
end

%% 8 — Write video

% Collected output folder, like the GUI's gui_renders_260815. The stem carries the
% recording, the ids and the window, so two clips never collide here.
if isempty(outFile)
    if isempty(outDir), outDirUse = fovPath; else, outDirUse = outDir; end
    if ~isfolder(outDirUse), mkdir(outDirUse); end

    stem = regexprep(sprintf('%s_%s_t%d-%ds', recName, idTag, ...
                             round(t_start), round(t_end)), '[\\/:*?"<>|]', '_');
    outFile = fullfile(outDirUse, [stem '.mp4']);
else
    outDirUse = fileparts(outFile);
    if ~isempty(outDirUse) && ~isfolder(outDirUse), mkdir(outDirUse); end
end

fprintf('Writing: %s\n', outFile);

% Seek by FRAME in the video's own timebase (see note in section 7c).
vBeh.CurrentTime = (f0 - 1 + toss_beh) / fps_beh;

tmpAvi = strrep(outFile, '.mp4', '_tmp.avi');

vw = VideoWriter(tmpAvi, 'Motion JPEG AVI');
vw.FrameRate = outFPS;
vw.Quality   = 98;
open(vw);

tic;

for ff = 1:nFrames
    tc = t_ca(f0+ff-1);

    frame = makeOneFrame(ff, tc, ...
        total_H, total_W, left_W, GAP, trace_W, TRACE_PLOT_W, ...
        leftBase, tileX, rowTop, TILE_SZ, ...
        vBeh, BREATH_CROP, BEH_CLIM, BEH_PRC_LO, BEH_PRC_HI, BEH_GAMMA, BEH_CLIP_PER_FRAME, ...
        ca_buf, ca_lo, ca_hi, ca_perim_up, gamma_ca, ...
        traceBase, t_ca, trData, trYLim, TRACE_WIN, ...
        tileStack_H, X_AXIS_H, FONT_NAME, FONT_SIZE_TRACE_AX, ...
        TRACE_AXIS_X_FRAC, TRACE_RIGHT_PAD_FRAC, X_AXIS_Y_FRAC, X_TICK_LABEL_YOFF, ...
        TRACE_ROW_TOP_PAD, TRACE_ROW_BOT_PAD, TRACE_LINE_W, TRACE_PLAYHEAD_W);

    writeVideo(vw, frame);

    if mod(ff,200)==0
        fprintf('  %d/%d (%.1fs)\n', ff, nFrames, toc);
    end
end

close(vw);

fprintf('  Wrote AVI in %.1fs, converting to mp4...\n', toc);

ffmpegExe = resolve_ffmpeg();
cmd = sprintf('"%s" -y -i "%s" -c:v libx264 -crf 18 -pix_fmt yuv420p "%s"', ...
              ffmpegExe, tmpAvi, outFile);
[st, fout] = system(cmd);

if st ~= 0
    warning('ffmpeg mp4 conversion failed:\n%s', fout);
else
    delete(tmpAvi);
end

fprintf('Done in %.1fs\nSaved: %s\n', toc, outFile);

%% ============================ LOCAL FUNCTIONS ============================

% =======================================================================
% ==================== CELL / ROI SELECTION (GUI RULES) =================
% =======================================================================

function [fovUse, roiIdxs, selIds, idLabels, idTag, idSrc, recName] = ...
        resolve_ids(fovPath, ids, idMode, regFile)
%RESOLVE_IDS  Typed ids -> one recording + the ROI indices to draw in it.
%
%  A number can mean two things and the script works out which: a registry CELL
%  id (global, the number on the summary figures) or a ROI index inside a
%  recording (local, 1..N in that FOV). 'auto' tries cell first, then ROI in the
%  preferred folder; 'cell' and 'roi' force one reading.
%
%  THE RECORDING FOLLOWS FROM THE IDS. Given a cell id the registry already knows
%  which recordings hold it, so fovPath is only a preference -- if the ids are not
%  there we move to a recording that has them and say so, the same thing the GUI
%  does when you load a cell that is not in the typed folder. A video cannot pool
%  recordings, so all ids must land in ONE folder; when they cannot, the error
%  lists where each one lives instead of guessing.
%
%  Identity priority is the GUI's: registry -> cell_link -> plain ROI. The
%  registry is authoritative because cell_link numbers cells per session while the
%  registry numbers them globally -- mixing the two returns a real but DIFFERENT
%  cell with no error.

    ids = ids(:).';
    assert(~isempty(ids), 'ids is empty -- nothing to export.');
    assert(all(isfinite(ids) & ids >= 1 & ids == round(ids)), ...
           'ids must be positive integers.');

    isAuto = strcmpi(idMode,'auto');
    assert(isAuto || strcmpi(idMode,'cell') || strcmpi(idMode,'roi'), ...
           'idMode must be ''auto'', ''cell'' or ''roi'' (got ''%s'')', idMode);

    hasPref = ~isempty(fovPath);
    if hasPref
        assert(isfolder(fovPath), 'Not a folder:\n  %s', fovPath);
        fovPath = regexprep(fovPath, '[\\/]+$', '');
    end

    % ---- roi mode: ids are ROI indices, no identity lookup at all ---------
    if strcmpi(idMode,'roi')
        assert(hasPref, 'idMode=''roi'' needs fovPath -- a ROI index means nothing on its own.');
        fovUse   = fovPath;
        roiIdxs  = ids;
        selIds   = ids;
        idSrc    = 'roi (as typed)';
        [idLabels, idTag] = label_ids(nan(size(ids)), roiIdxs);
        recName  = leaf_of(fovUse);
        return;
    end

    REG = get_registry(regFile);

    % ---- no registry: cell_link, then plain ROI ---------------------------
    if isempty(REG)
        assert(hasPref, ['no registry -- fovPath is required so the ids have a ' ...
                         'recording to be looked up in.']);
        [fovUse, roiIdxs, selIds, idLabels, idTag, idSrc] = ...
                resolve_via_cell_link(fovPath, ids, isAuto);
        recName = leaf_of(fovUse);
        return;
    end

    % ---- candidate readings of every id ----------------------------------
    % Rows of the preferred folder, if the registry covers it. Needed for the ROI
    % reading, which is only meaningful inside a known recording.
    kHere = [];
    if hasPref, kHere = reg_rows_for_folder(REG, fovPath); end

    cand = cell(1, numel(ids));
    for q = 1:numel(ids)
        c = struct('fp',{},'roi',{},'cellId',{},'how',{});

        % (a) the id as a CELL id -- every recording that holds it.
        %     Rows that ARE the preferred folder keep the path you typed: the
        %     registry stores archive paths, but the same recording is often
        %     opened from where it was acquired, and only that copy is guaranteed
        %     to have the behavior AVI next to it.
        for k = find([REG.cellNum] == ids(q))
            fpk = REG(k).fp;
            if any(kHere == k), fpk = fovPath; end
            c(end+1) = struct('fp',fpk, 'roi',REG(k).roi, ...
                              'cellId',ids(q), 'how','cell'); %#ok<AGROW>
        end

        % (b) the id as a ROI index in the preferred folder. Registered ROIs
        %     carry their cell id, so "29" and "279" label identically.
        if isAuto && ~isempty(kHere)
            for k = kHere([REG(kHere).roi] == ids(q))
                c(end+1) = struct('fp',fovPath, 'roi',ids(q), ...
                                  'cellId',REG(k).cellNum, 'how','roi'); %#ok<AGROW>
            end
        end

        % (c) an unregistered ROI in the preferred folder. Most ROIs are not
        %     registered cells -- only the ones that cleared the archive's gates
        %     are -- and those still have to be exportable.
        %
        %     This fires whenever the typed folder has no reading yet, NOT only
        %     when the id is unknown everywhere: a small number like 7 is a valid
        %     ROI in most FOVs and also a cell id somewhere in the archive, and
        %     silently leaving the recording you asked for is the worse answer.
        %     The guard is the folder's real ROI count from cpSAM -- 279 is not a
        %     ROI in a 40-ROI FOV, so a genuine cell id still moves you.
        if isAuto && hasPref && ~any(strcmpi({c.fp}, fovPath)) && ...
                ids(q) <= n_roi_in(fovPath)
            if ~isempty(c)
                fprintf(2, ['  NOTE %d is a ROI in the folder you typed AND cell %d ' ...
                            'in %s.\n       Staying in the folder you typed.\n'], ...
                        ids(q), ids(q), leaf_of(c(1).fp));
            end
            c(end+1) = struct('fp',fovPath, 'roi',ids(q), ...
                              'cellId',NaN, 'how','roi-unregistered'); %#ok<AGROW>
        end

        assert(~isempty(c), '%s', no_candidate_msg(REG, ids(q), fovPath, hasPref, isAuto));
        cand{q} = c;
    end

    % ---- one recording that satisfies every id ----------------------------
    common = unique({cand{1}.fp});
    for q = 2:numel(ids)
        common = intersect(common, unique({cand{q}.fp}));
    end
    if isempty(common)
        error('%s', no_common_folder_msg(ids, cand));
    end

    if hasPref && any(strcmpi(common, fovPath))
        fovUse = fovPath;
    else
        fovUse = common{1};
        if hasPref
            fprintf(2, ['  ids are not in the folder you typed -- using %s instead\n' ...
                        '  (typed: %s)\n'], leaf_of(fovUse), leaf_of(fovPath));
        end
        if numel(common) > 1
            fprintf(2, '  %d recordings hold all of these ids; using the first\n', ...
                    numel(common));
        end
    end

    % ---- pick one reading per id in that recording ------------------------
    roiIdxs = zeros(1,numel(ids));
    cellIds = nan(1,numel(ids));
    hows    = cell(1,numel(ids));
    for q = 1:numel(ids)
        c  = cand{q};
        c  = c(strcmpi({c.fp}, fovUse));
        % Cell reading wins when a number is both -- that is what idMode says it
        % is -- but an id that is genuinely ambiguous here gets said out loud.
        ic = find(strcmpi({c.how},'cell'), 1);
        if ~isempty(ic)
            if numel(c) > 1
                other = c(~strcmpi({c.how},'cell'));
                fprintf(2, ['  NOTE %d is BOTH cell %d and ROI %d (= cell %d) here; ' ...
                            'reading it as the cell.\n'], ids(q), c(ic).cellId, ...
                            other(1).roi, other(1).cellId);
            end
            pick = c(ic);
        else
            pick = c(1);
        end
        roiIdxs(q) = pick.roi;
        cellIds(q) = pick.cellId;
        hows{q}    = pick.how;
    end

    [idLabels, idTag] = label_ids(cellIds, roiIdxs);
    % Downstream keys on the CANONICAL id -- the cell id where there is one, the
    % ROI index otherwise -- so the number on the tile, in the filename and in the
    % YLIM_BY_ID copy-line are all the same number.
    selIds  = cellIds;
    selIds(~isfinite(selIds)) = roiIdxs(~isfinite(selIds));
    recName = leaf_of(fovUse);

    kUse = reg_rows_for_folder(REG, fovUse);
    idSrc = sprintf('registry (%d registered cells in this recording), read as: %s', ...
                    numel(unique([REG(kUse).cellNum])), strjoin(hows, ', '));
end

% -----------------------------------------------------------------------
function [fovUse, roiIdxs, selIds, idLabels, idTag, idSrc] = ...
        resolve_via_cell_link(fovPath, ids, isAuto)
%RESOLVE_VIA_CELL_LINK  Fallback identity for recordings the registry misses.
%  cell_link is per-session, so it can only ever resolve inside fovPath.
    fovUse = fovPath;
    lk = find_cell_link(fovPath);
    if isempty(lk)
        warning('repr:noCellIdentity', ...
            ['no registry and no cell_link for this recording -- treating the ' ...
             'ids as plain ROI indices.']);
        roiIdxs = ids;  selIds = ids;  idSrc = 'roi (automatic fallback)';
        [idLabels, idTag] = label_ids(nan(size(ids)), roiIdxs);
        return;
    end

    T    = lk.obsT;
    mm   = strcmpi(string(T.rec_path), fovPath);
    here = sort(T.cell_id(mm & ~isnan(T.cell_id))).';

    roiIdxs = zeros(1,numel(ids));
    cellIds = nan(1,numel(ids));
    for q = 1:numel(ids)
        hit = find(mm & T.cell_id == ids(q), 1);
        if ~isempty(hit)
            roiIdxs(q) = T.roi_index(hit);
            cellIds(q) = ids(q);
        elseif isAuto
            roiIdxs(q) = ids(q);          % read it as a ROI index instead
            hr = find(mm & T.roi_index == ids(q), 1);
            if ~isempty(hr), cellIds(q) = T.cell_id(hr); end
        else
            error(['cell %d is not in this recording (cell_link %s).\n' ...
                   '  cells here: %s'], ids(q), lk.srcFile, mat2str(here));
        end
    end

    [idLabels, idTag] = label_ids(cellIds, roiIdxs);
    selIds = cellIds;  selIds(~isfinite(selIds)) = roiIdxs(~isfinite(selIds));
    idSrc  = sprintf('cell_link (%s)', lk.srcFile);
end

% -----------------------------------------------------------------------
function [idLabels, idTag] = label_ids(cellIds, roiIdxs)
%LABEL_IDS  Tile labels and the filename tag, from the canonical id per tile.
%  A registered ROI is labelled by its CELL id, so "29" and "279" produce the same
%  tile and the same filename -- the clip is named after the neuron, not after how
%  it happened to be typed.
    n = numel(roiIdxs);
    idLabels = cell(1,n);
    tagParts = cell(1,n);
    for q = 1:n
        if isfinite(cellIds(q))
            idLabels{q} = sprintf('Cell %d', cellIds(q));
            tagParts{q} = sprintf('cell%03d', cellIds(q));
        else
            idLabels{q} = sprintf('ROI %d', roiIdxs(q));
            tagParts{q} = sprintf('ROI%02d', roiIdxs(q));
        end
    end
    idTag = strjoin(tagParts, '-');
end

% -----------------------------------------------------------------------
function n = n_roi_in(fovPath)
%N_ROI_IN  How many ROIs cpSAM found in a recording, for validating a ROI reading.
%  Read through matfile so the F matrix is never pulled into memory -- this is
%  asked once per id, before anything else is loaded. Inf when the answer cannot
%  be had, which lets the ROI reading through rather than blocking on a missing
%  file; the real bounds check against N_roi happens in section 2 either way.
    persistent CACHE
    if isempty(CACHE), CACHE = containers.Map('KeyType','char','ValueType','double'); end
    key = lower(regexprep(fovPath,'[\\/]+$',''));
    if CACHE.isKey(key), n = CACHE(key); return; end

    n = Inf;
    h = dir(fullfile(fovPath,'*_cpSAM_output.mat'));
    if ~isempty(h)
        f = fullfile(h(1).folder, h(1).name);
        try
            m  = matfile(f);
            sz = size(m,'F');
            if numel(sz) >= 2 && sz(2) > 0, n = sz(2); end
        catch
            try
                S = load(f,'F');  n = size(S.F,2);
            catch
            end
        end
    end
    CACHE(key) = n;
end

% -----------------------------------------------------------------------
function s = leaf_of(fp)
    p = strsplit(regexprep(fp,'[\\/]+$',''), {'\','/'});
    s = p{end};
end

% -----------------------------------------------------------------------
function msg = no_candidate_msg(REG, idVal, fovPath, hasPref, isAuto)
%NO_CANDIDATE_MSG  Say what the id could have been and what is actually there.
    L = sprintf('%d is not a registry cell id', idVal);
    if ~isAuto
        L = [L sprintf(['\n  idMode=''cell'' only accepts cell ids -- use ''auto'' ' ...
                        'to let it also mean a ROI index.'])];
    elseif ~hasPref
        L = [L sprintf(['\n  with fovPath = '''' there is no recording to read it ' ...
                        'as a ROI index in, either. Set fovPath, or use a cell id.'])];
    end
    if hasPref && ~isempty(REG)
        k = reg_rows_for_folder(REG, fovPath);
        if ~isempty(k)
            L = [L sprintf('\n  registered cells in %s:\n', leaf_of(fovPath))];
            for i = k(:)'
                L = [L sprintf('    cell %3d = ROI %d\n', REG(i).cellNum, REG(i).roi)]; %#ok<AGROW>
            end
        end
    end
    msg = L;
end

% -----------------------------------------------------------------------
function msg = no_common_folder_msg(ids, cand)
%NO_COMMON_FOLDER_MSG  No single recording holds all the ids -- show where each is.
    L = sprintf(['these ids are not all in one recording, and a video cannot pool ' ...
                 'recordings.\n']);
    for q = 1:numel(ids)
        L = [L sprintf('  %d is in:\n', ids(q))]; %#ok<AGROW>
        c = cand{q};
        [~, iu] = unique({c.fp});
        for i = iu(:)'
            L = [L sprintf('    %-55s ROI %d\n', leaf_of(c(i).fp), c(i).roi)]; %#ok<AGROW>
        end
    end
    msg = L;
end

% -----------------------------------------------------------------------
function REG = get_registry(regFile)
%GET_REGISTRY  Flatten event_latency_data.mat into one row per (cell, recording).
%  The SAME cell identity the batch renders with: the registry's cells with
%  cell_merge_overrides_260814 folded in, a merged cell taking the lowest of its
%  constituent numbers. So an id here is the same number as in
%  per-cell-summary_..._cell%03d.png and in the summary GUI's id box.
%  Returns [] when the registry is unavailable, which drops us to cell_link.
%
%  Ported from temporal_phase_cell_gui_260812.m -- keep in step with it.
    persistent CACHE CACHEKEY
    REG = [];
    if isempty(regFile) || ~isfile(regFile), return; end
    % Cache key includes SIZE and MODIFIED TIME, not just the path: keying on the
    % path alone meant a rebuilt registry was never picked up for the life of the
    % MATLAB session, with no error to say so.
    fi  = dir(regFile);
    key = sprintf('%s|%d|%.6f', regFile, fi.bytes, fi.datenum);
    if ~isempty(CACHE) && strcmp(CACHEKEY, key), REG = CACHE; return; end
    try
        D = load(regFile, 'CELL','OBS','REC');
        obsOf = pooled_obs_260814(D.CELL, D.OBS);
    catch ME
        warning('repr:noRegistry', ...
                'registry unavailable (%s) -- falling back to cell_link', ME.message);
        return;
    end
    rows = struct('cellNum',{},'fp',{},'recName',{},'roi',{},'group',{},'recDate',{});
    for c = 1:numel(obsOf)
        if isempty(obsOf{c}), continue; end       % merged away into a lower number
        for o = obsOf{c}(:)'
            p = regexp(D.OBS(o).label,'/','split');
            if numel(p) < 3, continue; end
            rows(end+1) = struct('cellNum',c, ...
                'fp',D.REC(D.OBS(o).rec).folder, ...
                'recName',strjoin(p(3:end-1),'/'), 'roi',str2double(p{end}), ...
                'group',p{1}, 'recDate',p{2}); %#ok<AGROW>
        end
    end
    REG = rows;  CACHE = REG;  CACHEKEY = key;
end

% -----------------------------------------------------------------------
function k = reg_rows_for_folder(REG, folderPath)
%REG_ROWS_FOR_FOLDER  Registry rows for a FOV folder, by path or by folder name.
%  The registry stores ARCHIVE paths (D:\Ventral_surface_summary\...), but the
%  same recording is often opened from where it was acquired (C:\260810_...\phys\).
%  Exact path first; if that misses, fall back to the leaf folder name, and only
%  when that name is unambiguous across the whole registry -- recording names do
%  repeat between sessions, and quietly matching the wrong session would put a
%  different animal's cell on screen.
    k = [];
    if isempty(REG), return; end
    k = find(strcmpi({REG.fp}, folderPath));
    if ~isempty(k), return; end

    leaf  = regexprep(folderPath, '[\\/]+$', '');
    parts = strsplit(leaf, {'\','/'});
    leaf  = parts{end};
    if isempty(leaf), return; end

    regLeaf = cell(1,numel(REG));
    for i = 1:numel(REG)
        q = strsplit(regexprep(REG(i).fp,'[\\/]+$',''), {'\','/'});
        regLeaf{i} = q{end};
    end
    hit = find(strcmpi(regLeaf, leaf));
    if isempty(hit), return; end
    if numel(unique({REG(hit).fp})) > 1
        warning('repr:ambiguousFolder', ...
            ['folder name "%s" appears in %d registry locations -- not matching ' ...
             'by name. Open the archive copy to use registry cell ids.'], ...
            leaf, numel(unique({REG(hit).fp})));
        return;
    end
    k = hit;
end

% -----------------------------------------------------------------------
function lk = find_cell_link(folderPath)
%FIND_CELL_LINK  Walk up from a FOV folder looking for the session cell_link.
%  Returns [] when there is none, which is the signal to fall back to ROI mode.
    lk = [];
    d = folderPath;
    for up = 1:4
        d = fileparts(d);
        if isempty(d); return; end
        cand = { fullfile(d,'analysis_260727','cell_pooled','cell_link.mat'), ...
                 fullfile(d,'cell_pooled','cell_link.mat'), ...
                 fullfile(d,'cell_link.mat') };
        for c = cand
            if isfile(c{1})
                S = load(c{1},'link');
                if isfield(S,'link') && isfield(S.link,'obsT')
                    lk = S.link; lk.srcFile = c{1};
                    return;
                end
            end
        end
    end
end

% =======================================================================
% ======================= VIDEO COMPOSITION =============================
% =======================================================================

function exe = resolve_ffmpeg()
% The ffmpeg on PATH on this machine is a broken conda build — it exits with
% 0xC0000139 (entrypoint not found) before doing anything. Prefer the working
% binary bundled with imageio_ffmpeg.

    persistent cached
    if ~isempty(cached)
        exe = cached;
        return;
    end

    exe = '';

    d = dir(fullfile(getenv('USERPROFILE'), '.conda', 'envs', 'dlc310', 'lib', ...
                     'site-packages', 'imageio_ffmpeg', 'binaries', 'ffmpeg*.exe'));
    if ~isempty(d)
        exe = fullfile(d(1).folder, d(1).name);
    else
        [s, p] = system(['python -c "import imageio_ffmpeg,sys; ' ...
                         'sys.stdout.write(imageio_ffmpeg.get_ffmpeg_exe())"']);
        if s == 0 && isfile(strtrim(p))
            exe = strtrim(p);
        end
    end

    if isempty(exe)
        warning('No working ffmpeg found; falling back to PATH.');
        exe = 'ffmpeg';
    end

    cached = exe;
end

function frame = makeOneFrame(ff, tc, ...
    total_H, total_W, left_W, GAP, trace_W, tracePlot_W, ...
    leftBase, tileX, rowTop, TILE_SZ, ...
    vBeh, BREATH_CROP, BEH_CLIM, BEH_PRC_LO, BEH_PRC_HI, BEH_GAMMA, BEH_CLIP_PER_FRAME, ...
    ca_buf, ca_lo, ca_hi, ca_perim_up, gamma_ca, ...
    traceBase, t_ca, trData, trYLim, TRACE_WIN, ...
    tileStack_H, X_AXIS_H, FONT_NAME, FONT_SIZE_AX, ...
    TRACE_AXIS_X_FRAC, TRACE_RIGHT_PAD_FRAC, X_AXIS_Y_FRAC, X_TICK_LABEL_YOFF, ...
    TRACE_ROW_TOP_PAD, TRACE_ROW_BOT_PAD, TRACE_LINE_W, TRACE_PLAYHEAD_W)

    frame = zeros(total_H, total_W, 3, 'uint8');

    leftPanel = leftBase;

    if hasFrame(vBeh)
        bimg = readFrame(vBeh);
    else
        bimg = zeros(TILE_SZ, TILE_SZ, 3, 'uint8');
    end

    breathTile = cropBreathTile(bimg, BREATH_CROP, TILE_SZ, ...
        BEH_CLIM, BEH_PRC_LO, BEH_PRC_HI, BEH_GAMMA, BEH_CLIP_PER_FRAME);
    leftPanel(rowTop(1):rowTop(1)+TILE_SZ-1, tileX:tileX+TILE_SZ-1, :) = breathTile;

    for k = 1:numel(ca_buf)
        ptch   = double(ca_buf{k}(:,:,ff));
        ptch_n = max(0,min(1,(ptch-ca_lo(k))/(ca_hi(k)-ca_lo(k)))).^gamma_ca;
        ptch_u = imresize(ptch_n,[TILE_SZ TILE_SZ],'bilinear');

        u8  = uint8(255*ptch_u);
        rgb = repmat(u8,[1 1 3]);

        pm = ca_perim_up{k};

        rgb(:,:,1) = rgb(:,:,1).*uint8(~pm) + uint8(pm)*255;
        rgb(:,:,2) = rgb(:,:,2).*uint8(~pm) + uint8(pm)*255;
        rgb(:,:,3) = rgb(:,:,3).*uint8(~pm);

        rr = rowTop(k+1):rowTop(k+1)+TILE_SZ-1;
        cc = tileX:tileX+TILE_SZ-1;

        leftPanel(rr, cc, :) = rgb;
    end

    frame(:, 1:left_W, :) = leftPanel;

    trPanel = makeTracePanelFast(traceBase, t_ca, trData, trYLim, tc, TRACE_WIN, ...
        rowTop, TILE_SZ, tileStack_H, X_AXIS_H, FONT_NAME, FONT_SIZE_AX, ...
        tracePlot_W, TRACE_AXIS_X_FRAC, TRACE_RIGHT_PAD_FRAC, X_AXIS_Y_FRAC, X_TICK_LABEL_YOFF, ...
        TRACE_ROW_TOP_PAD, TRACE_ROW_BOT_PAD, TRACE_LINE_W, TRACE_PLAYHEAD_W);

    frame(:, left_W+GAP+1:left_W+GAP+trace_W, :) = trPanel;
end

function tifPath = discover_mc_tiff(fovFolder)
    mc_hits = dir(fullfile(fovFolder,'*_MC*.tif'));
    assert(~isempty(mc_hits),'No *_MC*.tif in %s',fovFolder);

    names = {mc_hits.name};
    keep  = true(size(names));

    for ii = 1:numel(names)
        if contains(names{ii},'_AVG','IgnoreCase',true) || ...
           contains(names{ii},'_QC','IgnoreCase',true) || ...
           contains(names{ii},'_shifts','IgnoreCase',true) || ...
           contains(names{ii},'_SELROI','IgnoreCase',true)
            keep(ii) = false;
        end
    end

    mc_hits = mc_hits(keep);
    assert(~isempty(mc_hits));

    nMC = zeros(numel(mc_hits),1);
    for ii = 1:numel(mc_hits)
        nMC(ii) = numel(strfind(upper(mc_hits(ii).name),'_MC'));
    end

    [~,bst] = max(nMC);
    tifPath = fullfile(mc_hits(bst).folder, mc_hits(bst).name);
end

function BREATH_CROP = selectBreathCropInteractive(videoPath)
    v = VideoReader(videoPath);

    nSample = 200;
    sampleTimes = linspace(0, max(0, v.Duration - 1/v.FrameRate), nSample);

    avg = [];
    nUsed = 0;

    for i = 1:numel(sampleTimes)
        v.CurrentTime = sampleTimes(i);

        if ~hasFrame(v)
            continue;
        end

        fr = readFrame(v);

        if size(fr,3) == 3
            fr = rgb2gray(fr);
        end

        fr = single(fr);

        if isempty(avg)
            avg = zeros(size(fr),'single');
        end

        avg = avg + fr;
        nUsed = nUsed + 1;
    end

    avg = avg / max(1,nUsed);

    hf = figure('Color','w');
    imagesc(avg);
    axis image;
    colormap gray;
    title('Draw square crop for breathing tile, then double-click or press Enter');

    h = drawrectangle('AspectRatio',1,'Color','r');
    wait(h);

    pos = round(h.Position);
    ss  = round(min(pos(3), pos(4)));

    BREATH_CROP = [max(1,pos(1)), max(1,pos(2)), ss];

    close(hf);

    fprintf('  Selected crop: x=%d y=%d ss=%d\n', ...
        BREATH_CROP(1), BREATH_CROP(2), BREATH_CROP(3));
end

function BEH_CLIM = estimateBreathClipLimits(videoPath, BREATH_CROP, BEH_PRC_LO, BEH_PRC_HI)
    v = VideoReader(videoPath);

    nSample = 200;
    sampleTimes = linspace(0, max(0, v.Duration - 1/v.FrameRate), nSample);

    vals = [];

    for i = 1:numel(sampleTimes)
        v.CurrentTime = sampleTimes(i);

        if ~hasFrame(v)
            continue;
        end

        fr = readFrame(v);

        if size(fr,3) == 3
            fr = rgb2gray(fr);
        end

        fr = double(fr);

        [H,W] = size(fr);

        x1 = max(1, min(W, round(BREATH_CROP(1))));
        y1 = max(1, min(H, round(BREATH_CROP(2))));
        ss = round(BREATH_CROP(3));

        x2 = min(W, x1 + ss - 1);
        y2 = min(H, y1 + ss - 1);

        ss2 = min(x2-x1+1, y2-y1+1);

        crop = fr(y1:y1+ss2-1, x1:x1+ss2-1);

        vals = [vals; crop(:)]; %#ok<AGROW>
    end

    if isempty(vals)
        BEH_CLIM = [0 255];
        return;
    end

    lo = prctile(vals, BEH_PRC_LO);
    hi = prctile(vals, BEH_PRC_HI);

    if hi <= lo
        hi = lo + 1;
    end

    BEH_CLIM = [lo hi];
end

function tile = cropBreathTile(bimg, BREATH_CROP, TILE_SZ, ...
    BEH_CLIM, BEH_PRC_LO, BEH_PRC_HI, BEH_GAMMA, BEH_CLIP_PER_FRAME)

    if size(bimg,3) == 3
        bgray = rgb2gray(bimg);
    else
        bgray = bimg;
    end

    bgray = double(bgray);

    [H, W] = size(bgray);

    x1 = max(1, min(W, round(BREATH_CROP(1))));
    y1 = max(1, min(H, round(BREATH_CROP(2))));
    ss = round(BREATH_CROP(3));

    x2 = min(W, x1 + ss - 1);
    y2 = min(H, y1 + ss - 1);

    ss2 = min(x2-x1+1, y2-y1+1);

    crop = bgray(y1:y1+ss2-1, x1:x1+ss2-1);

    % Contrast clip AFTER crop.
    if BEH_CLIP_PER_FRAME || isempty(BEH_CLIM)
        lo = prctile(crop(:), BEH_PRC_LO);
        hi = prctile(crop(:), BEH_PRC_HI);
    else
        lo = BEH_CLIM(1);
        hi = BEH_CLIM(2);
    end

    if hi <= lo
        hi = lo + 1;
    end

    cropN = (crop - lo) ./ (hi - lo);
    cropN = max(0, min(1, cropN));

    if BEH_GAMMA ~= 1
        cropN = cropN .^ BEH_GAMMA;
    end

    tileGray = uint8(255 * cropN);
    tileGray = imresize(tileGray, [TILE_SZ TILE_SZ], 'bilinear');

    tile = repmat(tileGray, [1 1 3]);
end

function base = makeTraceBasePanel(trace_H, trace_W, tracePlot_W, ...
    rowTop, TILE_SZ, tileStack_H, X_AXIS_H, ...
    trYLim, traceLabels, fontName, fontSizeLabel, fontSizeAx, fontSizeXLab, ...
    traceAxisXFrac, traceRightPadFrac, ...
    traceLabelRightDX, traceLabelYShift, ...
    traceYNumRightDX, traceYNumTopYShift, traceYNumBotYShift, ...
    xAxisYFrac, xLabelYOff, ...
    traceRowTopPad, traceRowBotPad, traceAxisW)

    base = zeros(trace_H, trace_W, 3, 'uint8');

    axisX  = round(traceAxisXFrac * tracePlot_W);
    plotL  = axisX;
    plotR  = round((1-traceRightPadFrac) * tracePlot_W);

    labelRightX = axisX - traceLabelRightDX;
    yNumRightX  = axisX - traceYNumRightDX;

    axisCol = [255 255 255];
    textCol = [255 255 255];
    tickCol = [255 255 255];

    tickLen = 5;

    for kk = 1:numel(rowTop)
        rowYTop = rowTop(kk);
        rowYBot = rowTop(kk) + TILE_SZ - 1;

        yTop = rowYTop + traceRowTopPad;
        yBot = rowYBot - traceRowBotPad;
        yCtr = round((rowYTop + rowYBot)/2);

        ymin = trYLim(kk,1);
        ymax = trYLim(kk,2);

        base = drawLineRGB(base, axisX, yTop, axisX, yBot, axisCol, traceAxisW);
        base = drawLineRGB(base, axisX-tickLen, yTop, axisX+tickLen, yTop, axisCol, traceAxisW);
        base = drawLineRGB(base, axisX-tickLen, yBot, axisX+tickLen, yBot, axisCol, traceAxisW);

        base = drawTextRGBRight(base, labelRightX, ...
            yCtr - round(fontSizeLabel/2) + traceLabelYShift, ...
            traceLabels{kk}, textCol, fontSizeLabel, fontName);

        yTopText = sprintf('%.2g', ymax);
        yBotText = sprintf('%.2g', ymin);

        base = drawTextRGBRight(base, yNumRightX, yTop + traceYNumTopYShift, ...
            yTopText, tickCol, fontSizeAx, fontName);
        base = drawTextRGBRight(base, yNumRightX, yBot + traceYNumBotYShift, ...
            yBotText, tickCol, fontSizeAx, fontName);
    end

    xAxisY = tileStack_H + round(xAxisYFrac * X_AXIS_H);
    base = drawLineRGB(base, plotL, xAxisY, plotR, xAxisY, axisCol, traceAxisW);

    base = drawTextRGB(base, round((plotL+plotR)/2)-30, xAxisY+xLabelYOff, ...
        'Time (s)', textCol, fontSizeXLab, fontName);
end

function panel = makeTracePanelFast(traceBase, t_ca, trData, trYLim, tc, TRACE_WIN, ...
    rowTop, TILE_SZ, tileStack_H, X_AXIS_H, fontName, fontSizeAx, ...
    tracePlot_W, traceAxisXFrac, traceRightPadFrac, xAxisYFrac, xTickLabelYOff, ...
    traceRowTopPad, traceRowBotPad, traceLineW, tracePlayheadW)

    panel = traceBase;

    axisX = round(traceAxisXFrac * tracePlot_W);
    plotL = axisX;
    plotR = round((1-traceRightPadFrac) * tracePlot_W);

    halfWin = TRACE_WIN / 2;
    x1 = tc - halfWin;
    x2 = tc + halfWin;

    idx = find(t_ca >= x1 & t_ca <= x2);

    for kk = 1:numel(rowTop)
        rowYTop = rowTop(kk);
        rowYBot = rowTop(kk) + TILE_SZ - 1;

        yTop = rowYTop + traceRowTopPad;
        yBot = rowYBot - traceRowBotPad;

        ymin = trYLim(kk,1);
        ymax = trYLim(kk,2);

        if numel(idx) >= 2
            xp = plotL + round((t_ca(idx) - x1) ./ max(eps, x2-x1) * (plotR - plotL));
            yp = valueToYPix(trData(idx,kk), ymin, ymax, yTop, yBot);
            panel = drawPolylineRGB(panel, xp, yp, [255 255 255], traceLineW);
        end
    end

    cx = plotL + round((tc - x1) ./ max(eps, x2-x1) * (plotR - plotL));
    cx = max(plotL, min(plotR, cx));
    panel = drawLineRGB(panel, cx, rowTop(1)+traceRowTopPad, ...
        cx, rowTop(end)+TILE_SZ-1-traceRowBotPad, [255 35 35], tracePlayheadW);

    xAxisY = tileStack_H + round(xAxisYFrac * X_AXIS_H);
    tickLen = 6;
    xt = linspace(x1, x2, 6);

    for ii = 1:numel(xt)
        xp = plotL + round((xt(ii) - x1) ./ max(eps, x2-x1) * (plotR - plotL));
        xp = max(plotL, min(plotR, xp));

        panel = drawLineRGB(panel, xp, xAxisY, xp, xAxisY+tickLen, [255 255 255], 1);
        panel = drawTextRGB(panel, xp-14, xAxisY+xTickLabelYOff, sprintf('%.1f', xt(ii)), ...
            [255 255 255], fontSizeAx, fontName);
    end
end

function yp = valueToYPix(y, ymin, ymax, yTop, yBot)
    yp = yBot - round((y - ymin) ./ max(eps, ymax - ymin) * (yBot - yTop));
end

function img = drawPolylineRGB(img, x, y, col, lw)
    x = round(x(:));
    y = round(y(:));

    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    for i = 1:numel(x)-1
        img = drawLineRGB(img, x(i), y(i), x(i+1), y(i+1), col, lw);
    end
end

function img = drawLineRGB(img, x1, y1, x2, y2, col, lw)
    [H,W,~] = size(img);

    x1 = round(x1); y1 = round(y1);
    x2 = round(x2); y2 = round(y2);

    n = max(abs(x2-x1), abs(y2-y1)) + 1;
    xs = round(linspace(x1, x2, n));
    ys = round(linspace(y1, y2, n));

    rad = floor(lw/2);

    for q = 1:numel(xs)
        xx = xs(q);
        yy = ys(q);

        if xx < 1 || xx > W || yy < 1 || yy > H
            continue;
        end

        xr = max(1,xx-rad):min(W,xx+rad);
        yr = max(1,yy-rad):min(H,yy+rad);

        img(yr,xr,1) = col(1);
        img(yr,xr,2) = col(2);
        img(yr,xr,3) = col(3);
    end
end

function img = drawTextRGB(img, x, y, str, col, fontSize, fontName)
    if nargin < 7 || isempty(fontName)
        fontName = 'Arial';
    end

    if nargin < 6 || isempty(fontSize)
        fontSize = 14;
    end

    try
        img = insertText(img, [x y], str, ...
            'TextColor', col, ...
            'BoxOpacity', 0, ...
            'FontSize', fontSize, ...
            'Font', fontName);
    catch
        try
            img = insertText(img, [x y], str, ...
                'TextColor', col, ...
                'BoxOpacity', 0, ...
                'FontSize', fontSize);
        catch
        end
    end
end

function img = drawTextRGBRight(img, xRight, y, str, col, fontSize, fontName)
    if nargin < 7 || isempty(fontName)
        fontName = 'Arial';
    end

    if nargin < 6 || isempty(fontSize)
        fontSize = 14;
    end

    canvasW = 260;
    canvasH = 70;
    tmp = zeros(canvasH, canvasW, 3, 'uint8');

    try
        tmp = insertText(tmp, [2 2], str, ...
            'TextColor', col, ...
            'BoxOpacity', 0, ...
            'FontSize', fontSize, ...
            'Font', fontName);
    catch
        try
            tmp = insertText(tmp, [2 2], str, ...
                'TextColor', col, ...
                'BoxOpacity', 0, ...
                'FontSize', fontSize);
        catch
            return;
        end
    end

    mask = any(tmp > 0, 3);
    [rr, cc] = find(mask);

    if isempty(rr)
        return;
    end

    r1 = min(rr); r2 = max(rr);
    c1 = min(cc); c2 = max(cc);

    crop = tmp(r1:r2, c1:c2, :);
    maskCrop = mask(r1:r2, c1:c2);

    [H,W,~] = size(img);
    [h,w,~] = size(crop);

    x = round(xRight - w + 1);
    y = round(y);

    if x > W || y > H || x+w-1 < 1 || y+h-1 < 1
        return;
    end

    x1 = max(1,x);
    y1 = max(1,y);
    x2 = min(W,x+w-1);
    y2 = min(H,y+h-1);

    cx1 = x1 - x + 1;
    cy1 = y1 - y + 1;
    cx2 = cx1 + (x2-x1);
    cy2 = cy1 + (y2-y1);

    crop2 = crop(cy1:cy2, cx1:cx2, :);
    mask2 = maskCrop(cy1:cy2, cx1:cx2);

    for ch = 1:3
        region = img(y1:y2, x1:x2, ch);
        tmpCh = crop2(:,:,ch);
        region(mask2) = tmpCh(mask2);
        img(y1:y2, x1:x2, ch) = region;
    end
end
