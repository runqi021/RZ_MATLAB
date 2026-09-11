function falloff_projections_260805(dataDir, roiTag, powerLut, ellUm, zSurface, depthMode, mount, channel)
% FALLOFF_PROJECTIONS_260805  Depth-corrected max projections of a vessel stack.
%
% Companion to falloff_260804.m.  Builds ONE isotropic 1x1x1 um volume from the
% overlapping power rounds, corrects it for depth attenuation, and max-projects
% it three ways:
%
%   x-y   max projection in 100 um depth bins  (0-99, 100-199, 200-299, ...)
%   y-z   sagittal   -- max along x
%   z-x   coronal    -- max along y
%
% DEPTH CORRECTION
%   Each planar image is multiplied by A * exp(z / l), with A absorbing the
%   incident power (here 1/mW^2, so rounds at different powers land on one
%   scale) and l the attenuation length measured by falloff_260804.  Without
%   it the y-z and z-x views are unreadable: the top is saturated white and the
%   bottom is black, because the raw signal spans ~3 decades over the depth
%   range.  l is read from falloff_<roi>.mat unless passed in.
%
% GEOMETRY
%   Lateral pixel size = 1.7778 um at zoom 1 (rig constant, same convention as
%   detect_session_fps.m), scaled by the zoom and by pixelsPerLine.  Planes are
%   resampled to 1 um/px with BICUBIC interpolation.  The axial step is already
%   1 um so nothing is interpolated along z.  After that the volume is
%   isotropic, so reslicing to y-z and z-x is a pure permute -- no second
%   interpolation, which is what keeps the vessels from smearing.
%
% DEPTH AXIS comes from the filename z token plus the SET step size, exactly as
% in falloff_260804.m.  Motor zs is not used: its zero moves between stacks.
%
% Where rounds overlap, each 1 um plane is taken from the stack whose own depth
% range is centred nearest that depth -- i.e. the round that was exposed for it.
%
% USAGE
%   falloff_projections_260805('C:\fall-off\260804_shiver_dbh_vessel','roi1','post_260723')
%
% OUTPUT (into dataDir)
%   falloff_proj_<roiTag>.png / .pdf
%   falloff_proj_<roiTag>_vol.mat   cached isotropic volume; delete to re-read

%% --------------------------- USER PARAMETERS -------------------------------
DATA_DIR   = 'C:\fall-off\260804_shiver_dbh_vessel';
ROI_TAG    = 'roi1';
POWER_LUT  = 'post_260723';
P_EXPONENT = 2;

BIN_UM     = 100;          % depth bin for the x-y projections
DEPTH_MAX  = 600;          % FIXED depth axis for EVERY figure: 0-600 um at
                           % BIN_UM steps.  With 1x1x1 um voxels this puts x, y
                           % and z on one scale, so panel k of any two sessions
                           % shows the same slab at the same magnification.
                           % Shallower sessions get empty panels; deeper ones
                           % are cropped here.
CLIM_PCT   = [0.5 99.5];   % display window, percentiles of each panel
PIX_BASE   = 1.7778;       % um/px at zoom 1 for 512 px lines (rig constant)
PIX_TARGET = 1.0;          % target isotropic voxel, um
FLOOR_FRAC = 0.05;         % per-plane floor = mean of the dimmest this fraction;
                           % matches PCT_FRAC in falloff_260804
FIT_DEPTH  = [0 Inf];      % overwritten from falloff_<tag>.mat
Z_SURFACE  = 0;            % motor z of the pia (DEPTH_MODE 'motor' only)
DEPTH_MODE = 'filename';   % 'filename' | 'motor' -- see falloff_260804.m
MOUNT      = 'auto';       % 'auto' | 'old' | 'rot90' -- see below
CHANNEL    = [];           % SI channel carrying the vessel dye when more than one
                           % is saved (260806: channelSave [1;3], vessels on 3).
                           % Required there -- pages interleave, so the wrong
                           % channel projects the other label without an error.
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(dataDir),  DATA_DIR  = dataDir;  end
if nargin >= 2 && ~isempty(roiTag),   ROI_TAG   = roiTag;   end
if nargin >= 3 && ~isempty(powerLut),  POWER_LUT  = powerLut;  end
if nargin >= 5 && ~isempty(zSurface),  Z_SURFACE  = zSurface;  end
if nargin >= 6 && ~isempty(depthMode), DEPTH_MODE = depthMode; end
if nargin >= 7 && ~isempty(mount),     MOUNT      = mount;     end
if nargin >= 8 && ~isempty(channel),   CHANNEL    = channel;   end

% MOUNT decides what stage x and y MEAN anatomically, and therefore which of
% the two side views is coronal and which is sagittal.  Scan mirrors are fixed
% to the RIG, so an image column is ALWAYS stage x and a row ALWAYS stage y;
% the mount enters only at stage -> anatomy.  See [[project_rotated_mount_sessions]].
%   'old'   : stage x = medio-lateral, +stage y = rostral
%   'rot90' : -stage x = rostral,       stage y = medio-lateral   (>= 2026-07-21)
% The rule is TEMPORAL, never inferred from where the coordinates land.
if strcmpi(MOUNT, 'auto')
    [~, sess] = fileparts(DATA_DIR);
    t = regexp(sess, '^(\d{2})(\d{2})(\d{2})', 'tokens', 'once');
    assert(~isempty(t), ['Cannot read a YYMMDD date from folder ''%s'' -- pass ' ...
        'MOUNT explicitly as ''old'' or ''rot90'''], sess);
    acq = datetime(2000+str2double(t{1}), str2double(t{2}), str2double(t{3}));
    if acq >= datetime(2026,7,21), MOUNT = 'rot90'; else, MOUNT = 'old'; end
    fprintf('mount: %s (session %s, acquired %s)\n', MOUNT, sess, datestr(acq,'yyyy-mm-dd'));
else
    fprintf('mount: %s (supplied)\n', MOUNT);
end
rot = strcmpi(MOUNT, 'rot90');

[~, TAG] = match_stacks(DATA_DIR, ROI_TAG);   % filesystem-safe output name
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);              % laser_power_calibration lives here
addpath(fileparts(thisDir));   % repo root, unconditionally

%% --- attenuation length ------------------------------------------------------
if nargin >= 4 && ~isempty(ellUm)
    ell = ellUm;
    fprintf('attenuation length l = %.1f um (supplied)\n', ell);
else
    fp = fullfile(DATA_DIR, sprintf('falloff_%s.mat', TAG));
    assert(isfile(fp), ['%s not found -- run falloff_260804 first, or pass l ' ...
        'as the 4th argument.'], fp);
    S = load(fp, 'Lp', 'FIT_DEPTH');
    ell = S.Lp;
    if isfield(S, 'FIT_DEPTH'), FIT_DEPTH = S.FIT_DEPTH; end
    fprintf('attenuation length l = %.1f um (from %s), fit range %g-%g um\n', ...
        ell, fp, FIT_DEPTH(1), FIT_DEPTH(2));
end

%% --- isotropic volume --------------------------------------------------------
volFile = fullfile(DATA_DIR, sprintf('falloff_proj_%s_vol.mat', TAG));
useVol = false;
if isfile(volFile)
    S = load(volFile);
    % A volume built from another channel is a different label, not a variant.
    cCh = [];  if isfield(S, 'CHANNEL'), cCh = S.CHANNEL; end
    useVol = isequal(cCh, CHANNEL);
    if ~useVol
        fprintf('cached volume was channel %s, now asking for %s -- re-reading\n', ...
            mat2str(cCh), mat2str(CHANNEL));
    end
end
if useVol
    fprintf('using cached volume: %s  (delete to re-read)\n', volFile);
    V = S.V;  depth = S.depth;  pxUm = S.pxUm;  pctPlane = S.pctPlane;
else
    [V, depth, pxUm, pctPlane] = build_volume(DATA_DIR, ROI_TAG, ...
        PIX_BASE, PIX_TARGET, Z_SURFACE, DEPTH_MODE, CHANNEL);
    save(volFile, 'V', 'depth', 'pxUm', 'pctPlane', 'DATA_DIR', 'ROI_TAG', ...
        'CHANNEL', '-v7.3');
    fprintf('cached -> %s\n', volFile);
end
fprintf('volume %d x %d x %d  (y,x,z) at %g um isotropic, depth %g..%g um\n', ...
    size(V,1), size(V,2), size(V,3), PIX_TARGET, depth(1), depth(end));

%% --- depth correction --------------------------------------------------------
% A * exp(z/l) per plane.  The cached volume is RAW ADU and carries the Pockels
% setpoint per plane, so the power calibration is applied HERE -- a re-measured
% LUT then costs a replot, not a 25 GB re-read.
% (a) REMOVE THE FLOOR FIRST.  l was measured on (top - bot), i.e. with the
%     additive floor already subtracted, so exp(z/l) must not be applied to a
%     volume that still contains it.  The floor also GROWS with depth, so
%     amplifying it brightens the deep end for a reason unrelated to vessels --
%     that is what made the superficial layer look dark by comparison.
%     Done on the native planes, before any axial interpolation.
nFloor = max(1, round(FLOOR_FRAC * size(V,1) * size(V,2)));
for k = 1:size(V,3)
    pl = V(:,:,k);
    V(:,:,k) = pl - mean(mink(pl(:), nFloor));
end

% (b) CLAMP THE GAIN BEYOND THE FITTED RANGE.  The exponential was fitted over
%     FIT_DEPTH only.  Past that the signal stops decaying -- on these stacks it
%     even rises again -- so continuing to multiply by exp(z/l) explodes noise
%     and background.  Beyond the fit, hold the gain at its fit-edge value.
zGain = min(depth(:), FIT_DEPTH(2));
if any(depth(:) > FIT_DEPTH(2))
    fprintf('  gain clamped beyond %g um for %d of %d planes\n', ...
        FIT_DEPTH(2), nnz(depth(:) > FIT_DEPTH(2)), numel(depth));
end
mWplane = laser_power_calibration(pctPlane, POWER_LUT);
corr = reshape(single(exp(zGain / ell) ./ (mWplane(:) .^ P_EXPONENT)), 1, 1, []);
Vc   = V .* corr;
clear V

% --- axial interpolation to PIX_TARGET ---------------------------------------
% Lateral is already at PIX_TARGET (bicubic, done during the read).  If the
% axial step is coarser, interpolate along z with a shape-preserving cubic so
% the volume is genuinely isotropic before it is resliced.  pchip not spline:
% spline overshoots at vessel edges and can go negative.
dzNat = median(diff(depth));
if abs(dzNat - PIX_TARGET) > 1e-6
    dNew = (depth(1) : PIX_TARGET : depth(end)).';
    fprintf('  axial interpolation %g -> %g um: %d -> %d planes (pchip)\n', ...
        dzNat, PIX_TARGET, numel(depth), numel(dNew));
    sz = size(Vc);
    Vc = reshape(Vc, [], sz(3)).';                 % [nz, ny*nx]
    Vc = interp1(depth, Vc, dNew, 'pchip');   % keep single: double would ~4x RAM
    Vc = reshape(Vc.', [sz(1) sz(2) numel(dNew)]);
    depth = dNew;
end

%% --- projections -------------------------------------------------------------
edges = 0 : BIN_UM : DEPTH_MAX;
nBin  = numel(edges) - 1;

XY = cell(nBin,1);  lbl = cell(nBin,1);
for b = 1:nBin
    m = depth >= edges(b) & depth < edges(b+1);
    lbl{b} = sprintf('%g-%g \\mum', edges(b), edges(b+1));
    if any(m), XY{b} = max(Vc(:,:,m), [], 3); else, XY{b} = []; end
    fprintf('  bin %-14s %4d planes\n', lbl{b}, nnz(m));
end

% Image ROW = stage y, image COLUMN = stage x (scan mirrors are fixed to the
% rig).  Which of those is medio-lateral and which is antero-posterior depends
% on the MOUNT, so the two side views SWAP between old and rotated sessions:
%   coronal  = ML vs depth = max along the AP axis
%   sagittal = AP vs depth = max along the ML axis
inTissue = depth >= 0 & depth <= DEPTH_MAX;
if rot
    apDim = 2;  mlDim = 1;   % rot90: stage x (cols) = AP, stage y (rows) = ML
else
    apDim = 1;  mlDim = 2;   % old:   stage y (rows) = AP, stage x (cols) = ML
end
COR = squeeze(max(Vc(:,:,inTissue), [], apDim));   % (ML, z)  coronal
SAG = squeeze(max(Vc(:,:,inTissue), [], mlDim));   % (AP, z)  sagittal
dT  = depth(inTissue);

%% --- figure -------------------------------------------------------------------
f = figure('Color','w','Position',[30 30 max(1500, 300*nBin) 1000]);
tl = tiledlayout(f, 2, nBin, 'TileSpacing','compact','Padding','compact');
tl.Title.FontSize = 10;

axML = (0:size(COR,1)-1) * PIX_TARGET;   % in-plane axes in REAL um
axAP = (0:size(SAG,1)-1) * PIX_TARGET;
latX = axML;  latY = axAP;

for b = 1:nBin
    ax = nexttile(b);
    if isempty(XY{b})
        set(ax, 'XTick', [], 'YTick', [], 'Color', [0.93 0.93 0.93]);
        text(ax, 0.5, 0.5, 'no data', 'Units','normalized', ...
            'HorizontalAlignment','center', 'Color', [0.45 0.45 0.45]);
        box(ax,'on');
    else
        show(ax, XY{b}, CLIM_PCT, latX, latY);
        axis(ax,'image');
        if b == 1
            ylabel(ax, 'y (\mum)');
            add_scalebar(ax, latX(end), latY(end), 100);
        end
    end
    title(ax, lbl{b}, 'FontWeight','normal', 'FontSize', 10);
end

% Side views on a REAL um depth axis, fixed 0-DEPTH_MAX at BIN_UM ticks, so the
% z scale is identical to x and y and to every other session's figure.
% CORONAL first, then SAGITTAL.
nHalf = max(1, floor(nBin/2));
ax = nexttile(nBin+1, [1 nHalf]);
show(ax, COR.', CLIM_PCT, axML, dT);   % rows = depth, cols = ML
% DIRECTION of the ML axis, not just its name.  On a rot90 session the ML axis
% is the image ROW, and the image row runs along NEGATIVE stage y (see
% [[reference_stage_axis_convention]]), so increasing pixel index goes LATERAL
% -> MEDIAL -- the reverse of what 'medio-lateral' reads as.  The old mount puts
% ML on the image column and its direction has not been measured, so that label
% is left as it was rather than guessed at.
if rot
    mlLab = 'lateral \rightarrow medial (\mum)';
else
    mlLab = 'x, medio-lateral (\mum)';
end
side_axes(ax, mlLab, axML, DEPTH_MAX, BIN_UM);
title(ax, 'x-z  coronal view', 'FontWeight','normal');

ax = nexttile(nBin+1+nHalf, [1 nBin-nHalf]);
show(ax, SAG.', CLIM_PCT, axAP, dT);
side_axes(ax, 'y, antero-posterior (\mum)', axAP, DEPTH_MAX, BIN_UM);
title(ax, 'y-z  sagittal view', 'FontWeight','normal');

title(tl, {sprintf('%s  --  %s  --  MAX projections', DATA_DIR, ROI_TAG), ...
    sprintf('depth-corrected by exp(z/l), l = %.1f um   |   %g um isotropic, bicubic   |   display %g-%g%%   |   mount %s', ...
        ell, PIX_TARGET, CLIM_PCT(1), CLIM_PCT(2), MOUNT), ''}, ...
    'FontWeight','bold','Interpreter','none','FontSize',10);

if SAVE_FIG
    base = fullfile(DATA_DIR, sprintf('falloff_proj_%s', TAG));
    exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('\nsaved %s.png/.pdf\n', base);
end
end

% ============================== volume =======================================

function [V, depth, pxUm, pctPlane] = build_volume(DATA_DIR, ROI_TAG, ...
    PIX_BASE, PIX_TARGET, Z_SURFACE, DEPTH_MODE, CHANNEL)

files = match_stacks(DATA_DIR, ROI_TAG);
assert(~isempty(files), 'No stacks matching ''%s'' under %s', ROI_TAG, DATA_DIR);

% --- pass 1: metadata, depth ranges, power ---
S = struct('path',{},'d',{},'pct',{},'pctZ',{},'H',{},'W',{},'fpsl',{},'n',{}, ...
    'pxUm',{},'nCh',{},'chIdx',{});
for i = 1:numel(files)
    tp = fullfile(files(i).folder, files(i).name);
    info = imfinfo(tp);
    meta = si_meta(info);

    n    = si_num(meta, "SI.hStackManager.numSlices");
    fpsl = si_num(meta, "SI.hStackManager.framesPerSlice");
    dz   = abs(si_num(meta, "SI.hStackManager.stackZStepSize"));
    P0   = si_num(meta, "SI.hBeams.stackStartPower");
    zoom = si_num(meta, "SI.hRoiManager.scanZoomFactor");
    ppl  = si_num(meta, "SI.hRoiManager.pixelsPerLine");
    chSave = si_vec(meta, "SI.hChannels.channelSave");
    nCh    = max(numel(chSave), 1);
    chIdx  = pick_channel(chSave, CHANNEL, files(i).name);
    assert(numel(info) == n*fpsl*nCh, '%s: page count mismatch', files(i).name);

    z0set = si_num(meta, "SI.hStackManager.stackZStartPos");
    dzSgn = si_num(meta, "SI.hStackManager.stackZStepSize");
    if ~isfinite(dz) || dz == 0, dz = 1; end
    switch lower(DEPTH_MODE)
    case 'filename'
        t = regexp(files(i).name, '_z(-?\d+)-(\d+)_', 'tokens', 'once');
        assert(~isempty(t), ['No _z<start>-<end>_ token in %s -- if the ' ...
            'filename z is MOTOR z use DEPTH_MODE ''motor'''], files(i).name);
        dSt = str2double(t{1});  dStep = dz;
    case 'motor'
        if ~isfinite(dzSgn) || dzSgn == 0, dzSgn = -dz; end
        dSt   = Z_SURFACE - z0set;      % depth of the first slice
        dStep = -dzSgn;                 % motor going negative = depth increasing
    otherwise
        error('DEPTH_MODE must be ''filename'' or ''motor''');
    end

    H = info(1).Height;  W = info(1).Width;
    if ~isfinite(ppl) || ppl <= 0, ppl = W; end
    pxUm = PIX_BASE * (512/ppl) / zoom;

    k = numel(S) + 1;
    S(k).path = tp;  S(k).d = (dSt + (0:n-1)*dStep).';  S(k).n = n;
    % PER-SLICE setpoint.  A ramped stack (pzAdjust) changes power with depth;
    % storing only stackStartPower normalises every plane by the SURFACE power
    % and leaves the deep planes ~100x too bright.
    zs = si_vec(meta, "SI.hStackManager.zs");
    if numel(zs) ~= n, zs = (z0set + (0:n-1)*dzSgn).'; end
    P1   = si_num(meta, "SI.hBeams.stackEndPower");
    LC   = si_num(meta, "SI.hBeams.lengthConstants");
    pzOn = startsWith(lower(strtrim(si_line(meta, "SI.hBeams.pzAdjust"))), "true");
    if pzOn && ~(isfinite(P0) && isfinite(P1) && abs(P1-P0) < 1e-9)
        assert(isfinite(LC) && LC ~= 0, '%s: ramp on but lengthConstants unusable', files(i).name);
        zRef = si_num(meta, "SI.hStackManager.zPowerReference");
        if ~isfinite(zRef), zRef = zs(1); end
        S(k).pctZ = P0 .* exp((zs(:) - zRef) ./ LC);
    else
        S(k).pctZ = repmat(P0, n, 1);
    end
    S(k).pct = P0;
    S(k).H = H;  S(k).W = W;  S(k).fpsl = fpsl;  S(k).pxUm = pxUm;
    S(k).nCh = nCh;  S(k).chIdx = chIdx;
    fprintf('  %-36s %5g %%  %3d planes  %6.1f..%-6.1f um  %.4f um/px%s\n', ...
        files(i).name, S(k).pct, n, S(k).d(1), S(k).d(end), pxUm, ...
        string_if(nCh > 1, sprintf('  [ch %g]', chSave(chIdx))));
end
pxUm = S(1).pxUm;
assert(all(abs([S.pxUm] - pxUm) < 1e-9), 'Stacks have different pixel sizes');

% --- target grid ---
% Assemble on the NATIVE axial step, then interpolate to PIX_TARGET once.
% Sampling straight onto a 1 um grid from a 10 um stack would just repeat
% planes; and caching the interpolated volume would be ~10x the disk.
dzNat = min(arrayfun(@(s) median(abs(diff(s.d))), S));
depth = (0 : dzNat : max(arrayfun(@(s) max(s.d), S))).';
% Overlaps are split at their MIDPOINT: the shallower stack keeps the top half,
% the deeper stack takes the bottom half.  Same rule as falloff_260804.m, so the
% two figures describe the same volume.
[~, ord] = sort(arrayfun(@(s) min(s.d), S));
S = S(ord);
lo = -inf(numel(S),1);  hi = inf(numel(S),1);
for i = 1:numel(S)-1
    ovA = max(min(S(i).d), min(S(i+1).d));
    ovB = min(max(S(i).d), max(S(i+1).d));
    if ovB > ovA, mid = (ovA + ovB)/2;
    else,         mid = (max(S(i).d) + min(S(i+1).d))/2;
    end
    hi(i) = mid;  lo(i+1) = mid;
end
% Columns throughout: a row .* a column broadcasts into a matrix and find()
% then returns linear indices that overrun numel(S).
dMin = arrayfun(@(s) min(s.d), S).';
dMax = arrayfun(@(s) max(s.d), S).';
src = zeros(numel(depth),1);  idx = zeros(numel(depth),1);
for j = 1:numel(depth)
    ok = find(depth(j) >= dMin & depth(j) <= dMax & ...
              depth(j) >  lo   & depth(j) <= hi);
    if isempty(ok), continue; end
    src(j) = ok(1);
    [~, idx(j)] = min(abs(S(src(j)).d - depth(j)));
end
keep = src > 0;
depth = depth(keep);  src = src(keep);  idx = idx(keep);
pctPlane = arrayfun(@(j) S(src(j)).pctZ(idx(j)), (1:numel(src)).');
fprintf('  target grid: %d planes, %g..%g um\n', numel(depth), depth(1), depth(end));

sc = pxUm / PIX_TARGET;                       % >1 upsamples
probe = imresize(zeros(S(1).H, S(1).W, 'single'), sc, 'bicubic');
V = zeros(size(probe,1), size(probe,2), numel(depth), 'single');

for i = 1:numel(S)
    want = find(src == i);
    if isempty(want), continue; end
    fprintf('  reading %s: %d planes\n', S(i).path, numel(want));
    t = Tiff(S(i).path, 'r');  cl = onCleanup(@() t.close());
    for q = 1:numel(want)
        iz = idx(want(q));
        acc = zeros(S(i).H, S(i).W, 'single');
        for fr = 1:S(i).fpsl
            % Channels interleave page by page within each frame.
            setDirectory(t, ((iz-1)*S(i).fpsl + fr-1)*S(i).nCh + S(i).chIdx);
            acc = acc + single(t.read());
        end
        I = acc / S(i).fpsl;
        V(:,:,want(q)) = imresize(I, sc, 'bicubic');   % RAW; power applied later
    end
    clear cl
end
end

% --------------------------------------------------------------------------
function chIdx = pick_channel(chSave, CHANNEL, name)
% Position within the interleaved pages of the requested SI channel.  Every
% channel reads without error, so a multi-channel file must be told which.
nCh = max(numel(chSave), 1);
if nCh <= 1
    chIdx = 1;
    if ~isempty(CHANNEL) && ~isempty(chSave) && chSave(1) ~= CHANNEL
        warning('%s saves only channel %g; CHANNEL = %g ignored', name, chSave(1), CHANNEL);
    end
    return
end
assert(~isempty(CHANNEL), ['%s saves %d channels (%s) -- set CHANNEL to the ' ...
    'one carrying the vessel dye'], name, nCh, mat2str(chSave(:).'));
chIdx = find(chSave == CHANNEL, 1);
assert(~isempty(chIdx), '%s: channel %g not saved (saved: %s)', ...
    name, CHANNEL, mat2str(chSave(:).'));
end

function s = string_if(c, s)
if ~c, s = ''; end
end

% ============================== display ======================================

function show(ax, I, pct, xc, yc)
I = double(I);
v = I(isfinite(I));
lo = prctile(v, pct(1));  hi = prctile(v, pct(2));
if hi <= lo, hi = lo + 1; end
imagesc(ax, xc, yc, I, [lo hi]);
colormap(ax, gray);   % grayscale, NOT jet -- these are images, not measurements
set(ax, 'XTick', [], 'YTick', []);
end

function side_axes(ax, xlab, lat, depthMax, tick)
% Real-micron axes on a side view: depth fixed 0-depthMax at integer `tick`
% steps, lateral on the same micron scale.  Every session then plots the same
% box, which is the whole point of the summary figure.
ylim(ax, [0 depthMax]);
xlim(ax, [lat(1) lat(end)]);
set(ax, 'YTick', 0:tick:depthMax, 'YTickMode','manual', ...
        'XTick', 0:tick:lat(end), 'XTickMode','manual');
% The voxels are 1x1x1 um, so the DISPLAY must be isometric too -- 100 um of
% depth has to occupy the same distance on paper as 100 um of lateral extent,
% or the vessels are sheared and the two side views are not comparable.
daspect(ax, [1 1 1]);
xlabel(ax, xlab);  ylabel(ax, 'depth (\mum)');
box(ax,'on');
end

function add_scalebar(ax, W, H, umLen)
L = umLen;                       % axes are already in um
x0 = W*0.06;  y0 = H*0.94;
hold(ax,'on');
plot(ax, [x0 x0+L], [y0 y0], 'w-', 'LineWidth', 3);
text(ax, x0, y0-H*0.05, sprintf('%g \\mum', umLen), 'Color','w', ...
    'FontSize', 9, 'VerticalAlignment','bottom');
end

% ========================= ScanImage metadata helpers ========================

function meta = si_meta(info)
meta = "";
if isfield(info(1),'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1),'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
end
assert(strlength(meta) > 0, 'No ScanImage metadata');
meta = replace(meta, char(8629), newline);
end

function v = si_num(meta, key)
ln = si_line(meta, key);
if strlength(ln) == 0, v = NaN; return; end
tok = regexp(ln, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function vec = si_vec(meta, key)
vec = [];
ln = si_line(meta, key);
if strlength(ln) == 0, return; end
m = regexp(ln, '\[(.*)\]', 'tokens', 'once');
if isempty(m), src = ln; else, src = m{1}; end
nums = regexp(src, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(nums), return; end
vec = str2double(nums(:));
end

function ln = si_line(meta, key)
ln = "";
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
ln = extractBefore(s, newline);
if isempty(ln) || ismissing(ln), ln = s; end
end

function [files, tag] = match_stacks(DATA_DIR, ROI_TAG)
% ROI_TAG may be a plain prefix ('roi1', 'tile15') OR a glob when the part that
% distinguishes FOVs sits in the MIDDLE of the name ('*col00_row02*'), as on
% 251104 where every file starts with the same experiment string.
if contains(ROI_TAG, '*') || endsWith(lower(ROI_TAG), '.tif')
    pat = ROI_TAG;
    if ~endsWith(lower(pat), '.tif'), pat = [pat '.tif']; end
else
    pat = [ROI_TAG '*.tif'];
end
files = dir(fullfile(DATA_DIR, pat));
tag = regexprep(ROI_TAG, '[\*\.\s/\\]+', '_');
tag = regexprep(tag, '^_|_$', '');
end
