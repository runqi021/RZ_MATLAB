function map_subvolume_render_260530()

close all; clear all; clc;
%MAP_SUBVOLUME_RENDER_260530  Oblique MIP "volume render" of a sub-volume.
%
%   Extracts a BOX (X x Y um, full Z) from the ch?/avgz map tiles around a
%   stage coordinate, then renders it as a tilted maximum-intensity-
%   projection (the vasculature "volume render" look), with a red 3D
%   bounding box and a "X x Y x Z um^3" scale label.
%
%   The box is assembled directly from the per-tile avgz stacks using the
%   stage coordinates baked into the filenames (tileNN_<X>_<Y>_...), so no
%   prior stitch is required. Overlapping tiles are MAX-blended.
%
%   The render is an ORTHOGRAPHIC oblique MIP done by ray casting in real
%   microns (interp3): the view direction is set by azimuth/elevation, and
%   the wireframe box corners are projected through the SAME transform so
%   they always line up with the rendered tissue.
%
%   Standalone. MATLAB R2021a+ (Image Processing Toolbox for tiff/interp).
%   Runqi Zhang / 2026-05-30.

%% ========================= USER SETTINGS =========================

S.avgzDir   = 'C:\Users\Admin\Desktop\live_251104_definitely_not_sst_tdtomato_fitc\Map\ch1\avgz'; % multi-page avgz tiles
S.rawDir    = 'C:\Users\Admin\Desktop\live_251104_definitely_not_sst_tdtomato_fitc\Map';          % raw tiles (for metadata only)

%S.avgzDir   = 'C:\Users\Admin\Desktop\260114_homo_cy5\map\ch3\avgz'; % multi-page avgz tiles
%S.rawDir    = 'C:\Users\Admin\Desktop\260114_homo_cy5\map';          % raw tiles (for metadata only)

%C:\Users\Admin\Desktop\live_251104_definitely_not_sst_tdtomato_fitc\Map
%C:\Users\Admin\Desktop\live_251104_definitely_not_sst_tdtomato_fitc\Map\ch1\avgz

S.center_um = [1200 800];   % [X Y] stage coordinate (um) the box is centered on
S.box_um    = [400 100];    % [X Y] box footprint (um)
S.zUse      = [];           % Z pages to use ([] = all)

S.umPerPx   = 1.4648;       % XY um/px (calibrated; [] = auto from raw TIFF metadata)
S.dz_um     = 10;           % Z um/slice ([] = auto from stackZStepSize)
S.interp    = 'cubic';      % 'cubic'(bicubic/tricubic) | 'makima'(less overshoot) | 'spline' | 'linear'
S.xSign     = +1;           % within-tile: stage X vs image column
S.ySign     = -1;           % within-tile: stage Y vs image row (stage-Y usually opposite rows)

% --- volume preprocessing (applied to the 3D volume BEFORE projection) ---
S.gamma3d      = [];        % [] = off; e.g. 1.1 (>1 darkens midtones, pops bright vessels)
S.vesselness   = 'none';    % 'none' | 'fibermetric' (3D Frangi-style tube enhancement)
S.vesselThick_um = [4 24];  % vessel diameter range (um) for fibermetric scales
S.vesselSens   = [];        % fibermetric StructureSensitivity ([] = adaptive default)

% --- view / render ---
S.az        = -28;          % azimuth (deg) about the vertical axis (reveals depth/parallax)
S.el        =  14;          % elevation (deg): look down onto the top face
S.outPx_um  = 1.8;          % output image sampling (um/px)
S.rayStep_um= 2.0;          % MIP ray step along view direction (um)
S.clipPct   = [0 100];    % intensity percentile clip
S.gamma     = 0.60;         % display gamma (<1 brightens vessels)
S.invert    = false;        % true = dark vessels on white

% --- sides / orientation ---
S.bothSides  = true;        % render BOTH +X and -X (uses |center X|), saving each
S.medialLeft = true;        % orient so MEDIAL (midline, X->0) is on the LEFT, LATERAL on the RIGHT

% --- decoration / output ---
S.drawBox   = true;
S.boxColor  = [0.78 0.12 0.12];
S.label     = true;
S.savePNG   = true;
S.saveOrtho = true;         % also save flat XZ (Y-MIP) and YZ (X-MIP) projections
S.scalebar_um = 100;        % scale bar length on ortho planes (um; 0 = none)
S.outName   = '';           % '' = auto name in avgzDir's parent

%% ========================= GEOMETRY / METADATA =========================
if isempty(S.umPerPx) || isempty(S.dz_um)
    [umpp,dz] = read_si_meta(S.rawDir);
    if isempty(S.umPerPx), S.umPerPx = umpp; end
    if isempty(S.dz_um),   S.dz_um   = dz;   end
end
assert(~isempty(S.umPerPx) && ~isempty(S.dz_um), ...
    'Could not determine umPerPx / dz_um; set them in USER SETTINGS.');
fprintf('umPerPx = %.4f um   dz = %.3f um\n', S.umPerPx, S.dz_um);

%% ========================= LIST + LOCATE TILES =========================
c = [dir(fullfile(S.avgzDir,'*.tif')); dir(fullfile(S.avgzDir,'*.tiff'))];
assert(~isempty(c), 'No TIFFs in %s', S.avgzDir);
files = strings(numel(c),1); Xc = nan(numel(c),1); Yc = nan(numel(c),1);
for k = 1:numel(c)
    files(k) = string(fullfile(c(k).folder, c(k).name));
    [Xc(k), Yc(k)] = parse_xy(c(k).name);
end
keep = ~isnan(Xc) & ~isnan(Yc);
files = files(keep); Xc = Xc(keep); Yc = Yc(keep);
assert(~isempty(files), 'No tile filenames with tileNN_<X>_<Y> stage coords found.');

% one tile's page count + class
info  = imfinfo(char(files(1)));
nZall = numel(info); H0 = info(1).Height; W0 = info(1).Width;
zU = S.zUse; if isempty(zU), zU = 1:nZall; end
zU = zU(zU>=1 & zU<=nZall); nZ = numel(zU);

halfFovX = (W0/2)*S.umPerPx;  halfFovY = (H0/2)*S.umPerPx;

% ---- centers to render: BOTH +X and -X (mirror), or just the one set ----
if S.bothSides
    mag = abs(S.center_um(1));
    centers = [ mag, S.center_um(2);  -mag, S.center_um(2) ];   % right(+X) then left(-X)
else
    centers = S.center_um;
end

for ci = 1:size(centers,1)
cen = centers(ci,:);
bx0 = cen(1)-S.box_um(1)/2; bx1 = cen(1)+S.box_um(1)/2;
by0 = cen(2)-S.box_um(2)/2; by1 = cen(2)+S.box_um(2)/2;

% tiles whose FOV intersects the box
cover = (Xc+halfFovX>=bx0) & (Xc-halfFovX<=bx1) & ...
        (Yc+halfFovY>=by0) & (Yc-halfFovY<=by1);
tIdx = find(cover);
assert(~isempty(tIdx), 'No tiles cover the box around (%.0f, %.0f) um.', cen);
fprintf('Box X[%.0f %.0f] Y[%.0f %.0f] um, %d Z planes, %d contributing tiles.\n', ...
    bx0,bx1,by0,by1,nZ,numel(tIdx));

%% ========================= ASSEMBLE SUB-VOLUME =========================
Xg = bx0:S.umPerPx:bx1;  Yg = by0:S.umPerPx:by1;
W  = numel(Xg); Hh = numel(Yg);
[XX,YY] = meshgrid(Xg,Yg);                 % Hh x W  (stage um)
Vbox = -inf(Hh,W,nZ,'single');

cc = (W0+1)/2; rc = (H0+1)/2;
wb = waitbar(0,'Assembling sub-volume from tiles...');
for ii = 1:numel(tIdx)
    t = tIdx(ii);
    col = cc + S.xSign*(XX - Xc(t))/S.umPerPx;   % Hh x W
    row = rc + S.ySign*(YY - Yc(t))/S.umPerPx;
    inside = col>=1 & col<=W0 & row>=1 & row<=H0;
    if ~any(inside(:)), continue; end
    Tvol = read_pages(char(files(t)), zU, H0, W0);   % H0 x W0 x nZ single
    for z = 1:nZ
        samp = interp2(Tvol(:,:,z), col, row, S.interp, NaN);
        cur  = Vbox(:,:,z);
        m    = ~isnan(samp) & (samp>cur);
        cur(m) = samp(m);
        Vbox(:,:,z) = cur;
    end
    waitbar(ii/numel(tIdx), wb);
end
close(wb);
Vbox(~isfinite(Vbox)) = 0;

%% ============== VOLUME PREPROCESS (gamma3d / Frangi vesselness) ==============
depthVec = (0:nZ-1)*S.dz_um;               % um, page1 = surface (0)
depthTot = depthVec(end);                  % physical box depth (um)
Vol      = Vbox;
doVessel = strcmpi(S.vesselness,'fibermetric');

if doVessel
    % resample Z to (near-)isotropic so the 3D Hessian tube scales are physical
    ziso = linspace(0, depthTot, max(2,round(depthTot/S.umPerPx)+1));
    Vp   = permute(Vol,[3 1 2]);                    % nZ x Hh x W
    Vp   = interp1(depthVec(:), Vp, ziso(:), S.interp);
    Vol  = permute(Vp,[2 3 1]);                     % Hh x W x niso
    depthVec = ziso;
    fprintf('Vesselness: Z resampled %d -> %d planes (iso %.3f um)\n', nZ, numel(ziso), S.umPerPx);
end

% normalize to [0,1] before any nonlinear per-voxel op
if doVessel || ~isempty(S.gamma3d)
    vv = Vol(Vol>0); lo = prctile(vv,S.clipPct(1)); hi = prctile(vv,S.clipPct(2));
    if hi<=lo, hi = lo+1; end
    Vol = min(max((Vol-lo)/(hi-lo),0),1);
end

if doVessel
    th = max(1, round(S.vesselThick_um / S.umPerPx));   % vessel diameters (px)
    scales = unique(round(linspace(th(1), th(2), 5)));
    args = {'ObjectPolarity','bright'};
    if ~isempty(S.vesselSens), args = [args {'StructureSensitivity', S.vesselSens}]; end
    Vol = fibermetric(Vol, scales, args{:});
    fprintf('fibermetric: scales(px)=%s\n', mat2str(scales));
end

if ~isempty(S.gamma3d), Vol = Vol .^ S.gamma3d; end

%% ========================= OBLIQUE MIP RAY CAST =========================
C = [mean(Xg), mean(Yg), mean(depthVec)];

% screen basis: right(+x), up(toward surface), view(into screen, +Y)
right0 = [1 0 0]; up0 = [0 0 -1]; view0 = [0 1 0];
Raz   = rot3(deg2rad(S.az), [0 0 1]);      % azimuth about vertical (world Z)
right = Raz*right0(:);  view = Raz*view0(:);  up = up0(:);
Rel   = rot3(deg2rad(S.el), right);        % elevation about (rotated) right axis
up    = Rel*up;  view = Rel*view;
right = right/norm(right); up = up/norm(up); view = view/norm(view);

% orientation: keep MEDIAL (toward midline, X->0) on the LEFT, LATERAL on the RIGHT.
% +X hemisphere already satisfies this; -X hemisphere is mirrored across the YZ plane.
mirrorX = S.medialLeft && cen(1) < 0;
if mirrorX
    M = diag([-1 1 1]);  right = M*right; up = M*up; view = M*view;
end

% project the 8 box corners to set screen extents
[cx,cy,cz] = ndgrid([bx0 bx1],[by0 by1],[depthVec(1) depthVec(end)]);
corners = [cx(:) cy(:) cz(:)] - C;          % 8 x 3
us = corners*right; vs = corners*up; ws = corners*view;
m  = 6;                                      % um margin
uvec = (min(us)-m):S.outPx_um:(max(us)+m);
vvec = (min(vs)-m):S.outPx_um:(max(vs)+m);
wvec = (min(ws)-m):S.rayStep_um:(max(ws)+m);
[UU,VV] = meshgrid(uvec,vvec);               % nv x nu

img = zeros(size(UU),'single');
wb  = waitbar(0,'Ray casting oblique MIP...');
for iw = 1:numel(wvec)
    w  = wvec(iw);
    Xq = C(1) + UU*right(1) + VV*up(1) + w*view(1);
    Yq = C(2) + UU*right(2) + VV*up(2) + w*view(2);
    Zq = C(3) + UU*right(3) + VV*up(3) + w*view(3);
    samp = interp3(Xg, Yg, depthVec, Vol, Xq, Yq, Zq, S.interp, 0);
    img  = max(img, samp);
    if mod(iw,20)==0, waitbar(iw/numel(wvec), wb); end
end
close(wb);

%% ========================= DISPLAY + DECORATE =========================
v = img(img>0);
lo = prctile(v, S.clipPct(1)); hi = prctile(v, S.clipPct(2));
if hi<=lo, hi = lo+1; end
imn = (img-lo)/(hi-lo); imn = min(max(imn,0),1);
imn = imn.^S.gamma;
if S.invert, imn = 1-imn; bg=[1 1 1]; else, bg=[0 0 0]; end

fig = figure('Color',bg,'Name','Sub-volume oblique MIP','NumberTitle','off');
ax  = axes(fig); imagesc(ax, uvec, vvec, imn); colormap(ax, gray(256));
ax.YDir='normal'; axis(ax,'image'); ax.Color=bg; ax.XColor='none'; ax.YColor='none';
ax.XTick=[]; ax.YTick=[]; hold(ax,'on');

if S.drawBox
    drawWireBox(ax, [bx0 bx1],[by0 by1],[depthVec(1) depthVec(end)], C, right,up,view, S.boxColor);
end
if S.label
    txt = sprintf('%g \\times %g \\times %g \\mum^3', S.box_um(1), S.box_um(2), depthVec(end));
    tc  = 1-bg; % white on black
    text(ax, uvec(1)+0.04*range(uvec), vvec(1)+0.05*range(vvec), txt, ...
        'Color',tc,'FontSize',12,'FontWeight','bold','Interpreter','tex', ...
        'VerticalAlignment','bottom');
end
if cen(1) < 0, sideLbl = 'LEFT'; elseif cen(1) > 0, sideLbl = 'RIGHT'; else, sideLbl = 'MID'; end
title(ax, sprintf('%s  center (%.0f, %.0f) um  |  az %d^o el %d^o  |  medial<-  ->lateral', ...
    sideLbl,cen(1),cen(2),S.az,S.el),'Color',1-bg,'Interpreter','tex');
hold(ax,'off');

%% ========================= SAVE =========================
if S.savePNG
    if ~isempty(S.outName) && size(centers,1)==1
        outName = S.outName;
    else
        [par,~] = fileparts(S.avgzDir);          % .../chN
        [par,~] = fileparts(par);                % .../Map
        sfx = S.interp;
        if ~isempty(S.gamma3d), sfx = [sfx sprintf('_g3d%.2g',S.gamma3d)]; end
        if doVessel, sfx = [sfx '_frangi']; end
        outName = fullfile(par, sprintf('subvol_render_X%g_Y%g_%gx%gum_%s.png', ...
            cen(1),cen(2),S.box_um(1),S.box_um(2),sfx));
    end
    exportgraphics(fig, outName, 'BackgroundColor', bg, 'Resolution', 200);
    fprintf('Saved: %s\n', outName);

    % --- flat orthogonal MIP planes: XZ (collapse Y) and YZ (collapse X) ---
    if S.saveOrtho
        outBase = regexprep(outName,'\.png$','');
        % XZ: max over Y(dim1) -> [W(X) x nZ]; show as Z(down) x X.  X mirrored for -X side.
        Mxz = squeeze(max(Vol,[],1)).';
        save_plane_mip(Mxz, Xg, depthVec, sprintf('%s  XZ (Y-MIP)',sideLbl), ...
            [outBase '_XZ.png'], S, mirrorX, bg);
        % YZ: max over X(dim2) -> [Hh(Y) x nZ]; show as Z(down) x Y.  (no mirror; Y axis)
        Myz = squeeze(max(Vol,[],2)).';
        save_plane_mip(Myz, Yg, depthVec, sprintf('%s  YZ (X-MIP)',sideLbl), ...
            [outBase '_YZ.png'], S, false, bg);
    end
end

end   % ===== end of centers loop (both sides) =====
end

%% ========================= LOCAL HELPERS =========================
function save_plane_mip(M, hvec, vvec, ttl, outFile, S, mirrorH, bg)
% M: [numel(vvec) x numel(hvec)] image (rows = vertical/depth, cols = horizontal).
% hvec/vvec in um. Depth increases downward (surface on top). mirrorH flips X (medial-left).
v = M(M>0); if isempty(v), v = M(:); end
lo = prctile(v,S.clipPct(1)); hi = prctile(v,S.clipPct(2)); if hi<=lo, hi = lo+1; end
imn = min(max((M-lo)/(hi-lo),0),1); imn = imn.^S.gamma;
if S.invert, imn = 1-imn; end
fig = figure('Color',bg,'NumberTitle','off','Name',ttl);
ax  = axes(fig); imagesc(ax, hvec, vvec, imn); colormap(ax, gray(256));
axis(ax,'image'); ax.YDir = 'reverse';
if mirrorH, ax.XDir = 'reverse'; end
ax.Color = bg; ax.XColor = 'none'; ax.YColor = 'none'; ax.XTick = []; ax.YTick = [];
hold(ax,'on');
if S.drawBox
    rectangle(ax,'Position',[hvec(1) vvec(1) range(hvec) range(vvec)], ...
        'EdgeColor',S.boxColor,'LineWidth',1.4);
end
if S.scalebar_um > 0
    x1 = hvec(end) - 0.06*range(hvec); x0 = x1 - S.scalebar_um;
    yb = vvec(end) - 0.08*range(vvec);
    plot(ax,[x0 x1],[yb yb],'-','Color',1-bg,'LineWidth',3);
    text(ax,(x0+x1)/2, yb-0.02*range(vvec), sprintf('%g \\mum',S.scalebar_um), ...
        'Color',1-bg,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontWeight','bold','FontSize',10,'Interpreter','tex');
end
title(ax, ttl, 'Color',1-bg,'Interpreter','tex'); hold(ax,'off');
exportgraphics(fig, outFile, 'BackgroundColor', bg, 'Resolution', 200);
fprintf('Saved: %s\n', outFile);
end

function [X,Y] = parse_xy(name)
[~,bn] = fileparts(name); bn = char(bn);
% (1) explicit stage coords:  ..._x-1500_y-499_...
tx = regexp(bn,'_x(-?\d+\.?\d*)','tokens','once');
ty = regexp(bn,'_y(-?\d+\.?\d*)','tokens','once');
if ~isempty(tx) && ~isempty(ty)
    X = str2double(tx{1}); Y = str2double(ty{1}); return;
end
% (2) ScanImage tileNN_<X>_<Y>_...
tok = regexp(bn,'^tile\d+_(-?\d+)_(-?\d+)','tokens','once');
if ~isempty(tok), X = str2double(tok{1}); Y = str2double(tok{2});
else,            X = NaN;               Y = NaN;            end
end

function V = read_pages(file, zU, H0, W0)
V = zeros(H0,W0,numel(zU),'single');
tf = Tiff(file,'r'); cln = onCleanup(@() tf.close()); %#ok<NASGU>
for k = 1:numel(zU)
    tf.setDirectory(zU(k));
    V(:,:,k) = single(tf.read());
end
end

function R = rot3(theta, ax)
ax = ax(:)/norm(ax); x=ax(1); y=ax(2); z=ax(3);
c = cos(theta); s = sin(theta); C = 1-c;
R = [c+x^2*C,   x*y*C-z*s, x*z*C+y*s; ...
     y*x*C+z*s, c+y^2*C,   y*z*C-x*s; ...
     z*x*C-y*s, z*y*C+x*s, c+z^2*C];
end

function drawWireBox(ax, xr, yr, zr, C, right, up, view, col) %#ok<INUSD>
[cx,cy,cz] = ndgrid(xr,yr,zr);
P = [cx(:) cy(:) cz(:)];                 % 8 corners
u = (P-C)*right; v = (P-C)*up;           % project (orthographic)
% corner index from ndgrid([x],[y],[z]) is bit order: i + 2*(j-1)+4*(k-1)
idx = @(i,j,k) i + 2*(j-1) + 4*(k-1);
E = [ idx(1,1,1) idx(2,1,1); idx(1,2,1) idx(2,2,1); idx(1,1,2) idx(2,1,2); idx(1,2,2) idx(2,2,2); ... % x-edges
      idx(1,1,1) idx(1,2,1); idx(2,1,1) idx(2,2,1); idx(1,1,2) idx(1,2,2); idx(2,1,2) idx(2,2,2); ... % y-edges
      idx(1,1,1) idx(1,1,2); idx(2,1,1) idx(2,1,2); idx(1,2,1) idx(1,2,2); idx(2,2,1) idx(2,2,2) ];   % z-edges
for e = 1:size(E,1)
    plot(ax, u(E(e,:)), v(E(e,:)), '-', 'Color', col, 'LineWidth', 1.4);
end
end

function [umpp,dz] = read_si_meta(rawDir)
umpp = []; dz = [];
d = dir(fullfile(rawDir,'tile*.tif'));
if isempty(d), d = dir(fullfile(rawDir,'*.tif')); end
if isempty(d), return; end
info = imfinfo(fullfile(d(1).folder,d(1).name));
blob = '';
for f = {'Software','ImageDescription','Artist'}
    if isfield(info,f{1}) && ischar(info(1).(f{1}))
        blob = [blob, newline, info(1).(f{1})]; %#ok<AGROW>
    end
end
% XY: imagingFovUm = [x0 y0; x1 y1; ...]  -> full X span / pixelsPerLine
tok = regexp(blob,'imagingFovUm\s*=\s*\[([^\]]*)\]','tokens','once');
ppl = regexp(blob,'pixelsPerLine\s*=\s*(\d+)','tokens','once');
if ~isempty(tok)
    M = str2num(tok{1}); %#ok<ST2NM>
    if ~isempty(M)
        fovX = max(M(:,1))-min(M(:,1));
        np = 512; if ~isempty(ppl), np = str2double(ppl{1}); end
        umpp = fovX/np;
    end
end
% Z step
tz = regexp(blob,'stackZStepSize\s*=\s*(-?[\d.]+)','tokens','once');
if ~isempty(tz), dz = abs(str2double(tz{1})); end
end
