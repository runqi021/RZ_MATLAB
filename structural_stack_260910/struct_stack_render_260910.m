function struct_stack_render_260910(metaMat, varargin)
%STRUCT_STACK_RENDER_260910  Projection + 3-D figure for a structural stack.
%
%   struct_stack_render_260910(metaMat)
%   struct_stack_render_260910(metaMat, 'Name', value, ...)
%
% Stage 2 of the structural-map pipeline.  Takes the 1 um isotropic stacks
% written by struct_stack_prep_260910 and makes ONE FIGURE PER CHANNEL holding,
% in the layout the falloff project uses:
%
%   row 1   x-y MAX projection in 100 um depth bins
%   row 2   x-z coronal MAX  |  y-z sagittal MAX  |  true 3-D oblique render
%
% The 3-D panel is a REAL render, not two flat MIPs pasted on two faces of a
% cube: a screen basis is built from azimuth/elevation and the volume is sampled
% with interp3 ALONG the view direction in real microns, so structures at
% different depths shift relative to each other and the block has parallax.
% Orthographic, so distances stay measurable.  A projected wireframe of the
% bounding box is drawn over it, which is what makes the orientation readable.
%
% THE RAY ACCUMULATOR IS NOT A PLAIN MAX.  Any voxel pinned at the digitiser
% ceiling wins every ray that passes through it and smears into a bright band.
% The K brightest samples per ray are averaged instead (TopK); K = 1 restores a
% true MAX.  K is a COUNT, not a percentile -- there are only a few hundred
% steps per ray, so a 99.5th percentile rounds to 1 and silently gives back the
% MAX the accumulator exists to avoid.
%
% WHICH SIDE VIEW IS WHICH DEPENDS ON THE MOUNT.  Image row = stage y and image
% column = stage x always (the scan mirrors are bolted to the rig); the mount
% only sets what stage x and y MEAN anatomically.  Coronal = ML vs depth = max
% along the AP axis; sagittal = AP vs depth = max along ML.  So the two views
% SWAP between the old and the rotated mount.  'auto' decides from the YYMMDD
% prefix of the session folder against the 2026-07-21 boundary and prints its
% choice -- the rule is temporal, never inferred from where coordinates land.
%
% NOTHING IS DEPTH-CORRECTED HERE.  This stack was acquired with pzAdjust, i.e.
% the power ramp already compensates attenuation during acquisition, so applying
% exp(z/l) on top would double-count it.  'Norm' offers per-plane stretching if
% the deep half still reads dark, but it is OFF by default because it destroys
% the relative brightness that makes a structural map worth looking at.
%
% OPTIONS
%   'Channels'  []          SI channel numbers to render ([] = all in the meta)
%   'BinUm'     100         depth bin for the x-y panels
%   'ClimPct'   [0.5 99.5]  display window percentiles, per panel
%   'Norm'      'none'      'none' | 'perplane'  (see above)
%   'Az'        -32         3-D azimuth, deg
%   'El'        16          3-D elevation, deg
%   'TopK'      4           ray accumulator: mean of the K brightest samples
%   'RayStep'   2.0         um along the view direction
%   'OutPx'     2.0         um per pixel of the 3-D panel
%   'Mount'     'auto'      'auto' | 'old' | 'rot90'
%   'OutDir'    ''          default: beside the stacks
%
% Runqi Zhang / 2026-09-10

%% ------------------------------- options ----------------------------------
p = inputParser;
p.addParameter('Channels', []);
p.addParameter('BinUm',    100);
p.addParameter('ClimPct',  [0.5 99.5]);
p.addParameter('Norm',     'none');
p.addParameter('Az',       -32);
p.addParameter('El',        16);
p.addParameter('TopK',      4);
p.addParameter('RayStep',   2.0);
p.addParameter('OutPx',     2.0);
p.addParameter('Mount',    'auto');
p.addParameter('OutDir',   '');
p.parse(varargin{:});
opt = p.Results;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end     % NVIDIA hardware path drops ~1% of
                                        % exportgraphics PNGs as blank stubs
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

S = load(metaMat, 'meta'); meta = S.meta;
dataDir = fileparts(metaMat);
outDir  = opt.OutDir; if isempty(outDir), outDir = dataDir; end
if ~isfolder(outDir), mkdir(outDir); end

chList = opt.Channels; if isempty(chList), chList = meta.channels; end

fprintf('\n================== struct_stack_render_260910 ==================\n');
fprintf('meta   : %s\n', metaMat);
fprintf('volume : [%d %d %d] at %g um isotropic\n', meta.size(1), meta.size(2), meta.size(3), meta.pxOut);
fprintf('depth  : %+.0f -> %+.0f um  (depth = %g - motor z)\n', meta.depth(1), meta.depth(end), meta.zSurface);

%% -------------------------------- mount ------------------------------------
MOUNT = opt.Mount;
sess  = session_folder_name(meta.tifPath);
if strcmpi(MOUNT, 'auto')
    tok = regexp(sess, '^(\d{6})', 'tokens', 'once');
    assert(~isempty(tok), ['cannot read a YYMMDD prefix from session folder ''%s'' -- ' ...
        'pass Mount explicitly as ''old'' or ''rot90'''], sess);
    acq = datetime(tok{1}, 'InputFormat', 'yyMMdd');
    if acq >= datetime(2026,7,21), MOUNT = 'rot90'; else, MOUNT = 'old'; end
    fprintf('mount  : %s (session %s, acquired %s)\n', MOUNT, sess, datestr(acq,'yyyy-mm-dd')); %#ok<DATST>
else
    fprintf('mount  : %s (supplied)\n', MOUNT);
end
rot = strcmpi(MOUNT, 'rot90');
if rot
    apDim = 2; mlDim = 1;      % rot90: stage x (cols) = AP, stage y (rows) = ML
    mlLab = 'lateral \rightarrow medial (\mum)';
    apLab = 'x, antero-posterior (\mum)';
else
    apDim = 1; mlDim = 2;      % old:   stage y (rows) = AP, stage x (cols) = ML
    mlLab = 'x, medio-lateral (\mum)';
    apLab = 'y, antero-posterior (\mum)';
end

%% ------------------------------ per channel --------------------------------
for c = chList(:).'
    ci = find(meta.channels == c, 1);
    assert(~isempty(ci), 'channel %d is not in the meta file', c);
    tif = meta.outTif{ci};
    if ~isfile(tif)
        [~, nm, ex] = fileparts(tif);          % meta may carry the path it was
        tif = fullfile(dataDir, [nm ex]);      % written at, not where it is now
    end
    assert(isfile(tif), 'stack missing: %s', meta.outTif{ci});

    fprintf('\n--- channel %d ---\n  %s\n', c, tif);
    V = read_tiff_stack(tif);
    V = single(V);
    fprintf('  loaded [%d %d %d]\n', size(V,1), size(V,2), size(V,3));

    if strcmpi(opt.Norm, 'perplane')
        V = perplane_stretch(V, [5 99.8]);
        fprintf('  per-plane stretch applied (5-99.8%%)\n');
    end

    render_one_channel(V, meta, c, opt, rot, apDim, mlDim, mlLab, apLab, MOUNT, outDir, sess);
    clear V
end
fprintf('\ndone\n');
end

% =============================================================================
function render_one_channel(V, meta, ch, opt, rot, apDim, mlDim, mlLab, apLab, MOUNT, outDir, sess)

pt    = meta.pxOut;
depth = meta.depth(:);
BIN   = opt.BinUm;

%% --- depth bins --------------------------------------------------------------
% Bins are laid on the 0, 100, 200 ... grid so a panel means the same depth in
% every session.  Planes acquired ABOVE the surface are kept as their own
% leading bin rather than dropped -- on a structural map they are the only view
% of the surface itself.
hi = ceil(max(depth)/BIN)*BIN;
edges = 0 : BIN : hi;
if min(depth) < 0, edges = [floor(min(depth)/BIN)*BIN, edges]; end
nBin = numel(edges) - 1;

XY = cell(nBin,1); lbl = cell(nBin,1);
for b = 1:nBin
    m = depth >= edges(b) & depth < edges(b+1);
    if b == nBin, m = depth >= edges(b) & depth <= edges(b+1); end
    % Label the bin by what was ACTUALLY acquired in it, not by the nominal
    % edge.  Only the leading above-surface bin ever differs, and there a
    % nominal "-100 to 0" over 30 real planes reads as 70 um of missing data.
    e0 = edges(b);
    if b == 1 && e0 < 0 && any(m), e0 = floor(min(depth(m))); end
    lbl{b} = sprintf('%g to %g \\mum', e0, edges(b+1));
    if edges(b+1) <= 0, lbl{b} = [lbl{b} ' (above 0)']; end
    if any(m), XY{b} = max(V(:,:,m), [], 3); end
    fprintf('  bin %-22s %4d planes\n', sprintf('%g..%g', edges(b), edges(b+1)), nnz(m));
end

%% --- side views --------------------------------------------------------------
COR = squeeze(max(V, [], apDim));    % (ML, z)
SAG = squeeze(max(V, [], mlDim));    % (AP, z)
axML = (0:size(COR,1)-1) * pt;
axAP = (0:size(SAG,1)-1) * pt;
latX = (0:size(V,2)-1) * pt;
latY = (0:size(V,1)-1) * pt;

%% --- 3-D oblique render ------------------------------------------------------
fprintf('  ray-casting (az %g, el %g, K %d, step %g um, %g um/px) ...\n', ...
        opt.Az, opt.El, opt.TopK, opt.RayStep, opt.OutPx);
t3 = tic;
R = raycast_volume(V, pt, depth, opt.Az, opt.El, opt.TopK, opt.RayStep, opt.OutPx);
fprintf('    %.0f s -> [%d %d] px\n', toc(t3), size(R.img,1), size(R.img,2));

%% --- figure ------------------------------------------------------------------
% The canvas is capped at the screen, so the requested size is read BACK and the
% design scaled by one factor rather than positioning anything in raw pixels.
want = [max(1500, 260*nBin), 1000];
f = figure('Color','w','Units','pixels','Position',[30 30 want]);
drawnow; got = get(f,'Position'); FS = min(1, got(3)/want(1));
set(f,'DefaultAxesFontSize', max(6, 10*FS), 'DefaultTextFontSize', max(6, 10*FS));

% THREE tile rows, not two, with the side views spanning the lower two.  The
% x-y panels are square while a side view is as tall as the stack is deep, so
% equal rows leave a broad empty band under row 1 and squash the side views.
% A 1:2 split matches the geometry (455 um of field over 530 um of depth).
nCol = max(nBin, 6);
tl = tiledlayout(f, 3, nCol, 'TileSpacing','compact', 'Padding','compact');

for b = 1:nBin
    ax = nexttile(tl, b);
    if isempty(XY{b})
        set(ax, 'XTick', [], 'YTick', [], 'Color', [0.93 0.93 0.93]); box(ax,'on');
        text(ax, 0.5, 0.5, 'no data', 'Units','normalized', ...
             'HorizontalAlignment','center', 'Color', [0.45 0.45 0.45]);
    else
        show_img(ax, XY{b}, opt.ClimPct, latX, latY);
        axis(ax, 'image');
        if b == 1
            ylabel(ax, 'y (\mum)');
            add_scalebar(ax, 100);
        end
    end
    title(ax, lbl{b}, 'FontWeight','normal');
end

% Row 2 is split three ways: coronal | sagittal | 3-D.  The side views carry a
% REAL micron depth axis and daspect [1 1 1]; the voxels are isotropic, so the
% display has to be too or the structures come out sheared.
w = floor(nCol/3); w3 = nCol - 2*w;

ax = nexttile(tl, nCol+1, [2 w]);
show_img(ax, COR.', opt.ClimPct, axML, depth);
side_axes(ax, mlLab, axML, depth, BIN);
title(ax, 'coronal  (max along AP)', 'FontWeight','normal');

ax = nexttile(tl, nCol+1+w, [2 w]);
show_img(ax, SAG.', opt.ClimPct, axAP, depth);
side_axes(ax, apLab, axAP, depth, BIN);
title(ax, 'sagittal  (max along ML)', 'FontWeight','normal');

ax = nexttile(tl, nCol+1+2*w, [2 w3]);
show_img(ax, R.img, [1 99.7], R.uvec, R.vvec);
axis(ax, 'image'); set(ax, 'YDir', 'normal');
hold(ax, 'on');
plot_wireframe(ax, R);
hold(ax, 'off');
set(ax, 'XTick', [], 'YTick', []);
title(ax, sprintf('3-D oblique  (az %g\\circ, el %g\\circ)', opt.Az, opt.El), 'FontWeight','normal');
add_scalebar(ax, 100);

title(tl, { sprintf('%s   |   %s   |   channel %d   |   MAX projections', ...
                sess, meta_stem(meta), ch), ...
            sprintf(['%g um isotropic (%.4f um/px lateral, %g um step, %d frames/slice averaged)   |   ' ...
                     'display %g-%g%%   |   norm %s   |   mount %s   |   depth = %g - motor z'], ...
                meta.pxOut, meta.pxUm, abs(meta.dz), meta.framesPerSlice, ...
                opt.ClimPct(1), opt.ClimPct(2), opt.Norm, MOUNT, meta.zSurface), '' }, ...
      'FontWeight','bold', 'Interpreter','none', 'FontSize', max(7, 10*FS));

base = fullfile(outDir, sprintf('%s_ch%d_projections', meta_stem(meta), ch));
exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
d = dir([base '.png']);
fprintf('  saved %s.png (%.0f kB) / .pdf\n', base, d.bytes/1024);
if d.bytes < 20000
    warning('%s.png is only %d bytes -- likely an OpenGL export stub, re-run.', base, d.bytes);
end
end

% =============================================================================
function R = raycast_volume(V, pt, depth, AZ, EL, K, rayStep, outPx)
%RAYCAST_VOLUME  Orthographic oblique MIP-like render, top-K accumulator.
%
% World frame: X = image column (um), Y = image row (um), Z = -depth (um), so
% +Z is up and the shallow end of the stack sits at the top of the render.

[ny, nx, nz] = size(V);
xv = (0:nx-1) * pt;
yv = (0:ny-1) * pt;
zv = (-flip(depth(:))).';        % ascending ROW vector, as interp3 wants
Vz = flip(V, 3);                 % flipped to match

% Screen basis, MATLAB's view(az,el) convention.  w points from the volume
% toward the camera; rays march along -w.
az = AZ; el = EL;
wv = [ sind(az)*cosd(el), -cosd(az)*cosd(el), sind(el) ];
uv = [ cosd(az),           sind(az),          0        ];
vv = cross(wv, uv);
wv = wv/norm(wv); uv = uv/norm(uv); vv = vv/norm(vv);

% Bounding box of the volume, projected onto the basis, so the screen is exactly
% big enough and the ray range exactly spans the block.
[cx, cy, cz] = ndgrid([xv(1) xv(end)], [yv(1) yv(end)], [zv(1) zv(end)]);
P = [cx(:), cy(:), cz(:)];
C = mean([xv(1) xv(end); yv(1) yv(end); zv(1) zv(end)], 2).';
Q = P - C;
uu = Q*uv.'; vvp = Q*vv.'; ww = Q*wv.';

uvec = (min(uu) : outPx : max(uu));
vvec = (min(vvp): outPx : max(vvp));
svec = (min(ww) : rayStep : max(ww));
[UU, VV] = meshgrid(uvec, vvec);
npx = numel(UU);

TK = zeros(npx, K, 'single');
base = [UU(:), VV(:)] * [uv; vv] + C;      % [npx x 3] screen-plane points

for s = svec
    Pq = base + s*wv;
    q = interp3(xv, yv, zv, Vz, Pq(:,1), Pq(:,2), Pq(:,3), 'linear', 0);
    q = single(q);
    m = q > TK(:,K);
    if any(m)
        TK(m,K) = q(m);
        if K > 1, TK(m,:) = sort(TK(m,:), 2, 'descend'); end
    end
end

R.img  = reshape(mean(TK, 2), size(UU));
R.uvec = uvec;  R.vvec = vvec;
R.u = uv; R.v = vv; R.w = wv; R.C = C;
R.corners = P;
end

function plot_wireframe(ax, R)
%PLOT_WIREFRAME  The 12 edges of the bounding box, projected onto the screen.
% Without it an oblique render has no readable orientation.
P = R.corners;                                  % ndgrid order: x fastest
idx = [1 2; 3 4; 5 6; 7 8; 1 3; 2 4; 5 7; 6 8; 1 5; 2 6; 3 7; 4 8];
Q = P - R.C;
u = Q*R.u.'; v = Q*R.v.';
for e = 1:size(idx,1)
    plot(ax, u(idx(e,:)), v(idx(e,:)), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
end
end

% =============================================================================
function V = read_tiff_stack(fn)
t = Tiff(fn, 'r'); c = onCleanup(@() close(t));
A = t.read();
n = 1;
while ~t.lastDirectory()
    t.nextDirectory(); n = n + 1;
end
V = zeros(size(A,1), size(A,2), n, 'like', A);
t.setDirectory(1);
for k = 1:n
    V(:,:,k) = t.read();
    if k < n, t.nextDirectory(); end
end
end

function V = perplane_stretch(V, pct)
for k = 1:size(V,3)
    p = prctile(reshape(V(:,:,k), [], 1), pct);
    V(:,:,k) = (V(:,:,k) - p(1)) / max(eps, p(2) - p(1));
end
end

function show_img(ax, I, pct, xc, yc)
v = double(I(isfinite(I)));
lo = prctile(v, pct(1)); hi = prctile(v, pct(2));
if hi <= lo, hi = lo + 1; end
imagesc(ax, xc, yc, I, [lo hi]);
colormap(ax, gray);
set(ax, 'YDir', 'reverse', 'Box', 'on', 'Layer', 'top', 'TickDir', 'out');
end

function side_axes(ax, xlab, xv, depth, BIN)
% The voxels are isotropic, so 100 um of depth must occupy the same distance on
% paper as 100 um laterally -- otherwise the side view is sheared.
daspect(ax, [1 1 1]);
xlim(ax, [xv(1) xv(end)]);
ylim(ax, [depth(1) depth(end)]);
set(ax, 'YTick', ceil(depth(1)/BIN)*BIN : BIN : depth(end), 'YDir', 'reverse');
xlabel(ax, xlab); ylabel(ax, 'depth (\mum)');
if depth(1) < 0
    hold(ax,'on');
    plot(ax, [xv(1) xv(end)], [0 0], ':', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.75);
    hold(ax,'off');
end
end

function add_scalebar(ax, L)
%ADD_SCALEBAR  100 um bar in the corner of an image axes.
%
% The limits are read BACK off the axes rather than passed in: `axis image` has
% usually just changed them, and a bar placed from the pre-axis-image extent
% lands outside the panel.  The 3-D panel also runs YDir normal with a v axis
% that straddles zero, so the corner has to be found from the limits, not
% assumed to be at the largest coordinate.
xl = xlim(ax); yl = ylim(ax);
x0 = xl(2) - 0.06*diff(xl) - L;
if strcmp(get(ax, 'YDir'), 'reverse')
    y0 = yl(2) - 0.08*diff(yl);
else
    y0 = yl(1) + 0.08*diff(yl);
end
hold(ax, 'on');
plot(ax, [x0 x0+L], [y0 y0], '-w', 'LineWidth', 2.5);
text(ax, x0 + L/2, y0, sprintf('%g \\mum', L), 'Color','w', ...
     'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 8);
hold(ax, 'off');
end

function s = meta_stem(meta)
[~, s] = fileparts(meta.tifPath);
end

function s = session_folder_name(tifPath)
% <session>\<recording folder>\<file>.tif  -- the session is two levels up, but
% fall back to one level if the file was not filed into its own folder.
d = fileparts(tifPath);
[up, name] = fileparts(d);
if ~isempty(regexp(name, '^\d{6}', 'once'))
    s = name;
else
    [~, s] = fileparts(up);
end
end
