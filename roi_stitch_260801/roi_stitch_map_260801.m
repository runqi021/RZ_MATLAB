function roi_stitch_map_260801()
%% roi_stitch_map_260801  ONE XY map of every FOV. Click a tile, arrow it into place.
% -----------------------------------------------------------------------
% All FOVs of an experiment in one XY projection. Click a tile, nudge it with the
% arrow keys, or let phase correlation place it.
%
% Each FOV is its own image object positioned in micrometres, so a nudge just
% rewrites that object's XData/YData. Nothing is recomposited and nothing is
% resampled, so tiles stay at native resolution and 1 um steps are visible when
% you zoom in.
%
% DEPTH IS COLOUR. Tiles are tinted on a neon ramp against black: BLUE = DEEP,
% through cyan / green / yellow to MAGENTA = SHALLOW. Alpha is each tile's own
% brightness, so dark background is transparent and structure is opaque -- where
% two tiles at the same depth agree, structure reinforces; where they disagree
% you see it doubled, which is the thing you are nudging out.
%
% SITES ARE PACKED. Laterally separate sites (roi1 and roi2 sit ~1200 um apart)
% are given a DISPLAY-ONLY offset so they sit side by side instead of separated
% by dead space. This offset never reaches disk: it is excluded from tx/ty, so
% fov_offsets.mat and all matching stay in the true frame.
%
% ANCHOR = EARLIEST ACQUIRED. Stage error accumulates through a session, so the
% first recording is the least-drifted reference. It never moves, and
% auto-registration prefers to register a tile against EARLIER tiles whenever any
% overlap it -- so corrections propagate outward from the best coordinate rather
% than chaining through late, drifted ones.
%
% CONTROLS
%   click a tile    select it (topmost tile under the cursor)
%   arrow keys      move selected tile by STEP_UM
%   shift + arrow   move by 10x STEP_UM
%   g               auto-register selected tile (phase correlation)
%   shift + g       auto-register every tile, earliest first
%   [ / ]           previous / next depth band      0  show all bands
%   n / space       select next tile
%   r               reset selected tile to its raw stage position
%   h               hide / show selected tile
%   a               toggle ROI markers        o  toggle tile outlines
%   s               save now                  q / close  save and quit
%
% Writes <outDir>\fov_offsets.mat in the format roi_stitch_match_260801 reads.

%% ===================== USER-EDITABLE =====================
STEP_UM    = 1.0;     % arrow-key step
BIG_MULT   = 10;      % shift+arrow multiplier
AUTOREG_DZ = 15;      % um; only tiles within this depth register against each other
MAXSHIFT   = 60;      % um; a larger implied correction is not believed
SHOW_ROI   = true;
SHOW_BOX   = true;
% Display-only packing: slide laterally separate sites together so dead space
% between them does not dominate the view. OFF by default -- true coordinates are
% what you almost always want, and packing silently moves tiles by hundreds of um.
% It is also REFUSED automatically when sites straddle the midline (see below).
PACK_SITES = false;
PACK_GAP   = 80;      % um of space between packed sites
BAND_UM    = 25;      % depth-band height for the [ / ] filter
DIM_ALPHA  = 0.10;    % alpha multiplier for tiles outside the active band
%% =========================================================

cfg = stitch_cfg_260801();
tf  = fullfile(cfg.outDir,'fov_table.mat');
assert(isfile(tf), 'Run roi_stitch_scan_260801 first (%s missing).', tf);
S = load(tf); F = S.F; nF = numel(F);
assert(isfield(F,'acqOrder'), 'fov_table.mat predates acquisition ordering -- re-run roi_stitch_scan_260801.');

% ---- resume ----------------------------------------------------------------
tx = zeros(nF,1); ty = zeros(nF,1); placed = false(nF,1);
if isfile(cfg.offsetFile)
    O = load(cfg.offsetFile);
    if isfield(O,'fovTable')
        nm = string({F.name});
        for i = 1:nF
            j = find(O.fovTable.fov == nm(i), 1);
            if ~isempty(j)
                tx(i)=O.fovTable.tx_um(j); ty(i)=O.fovTable.ty_um(j);
                placed(i)=O.fovTable.placed(j);
            end
        end
        fprintf('resumed %d previously moved tile(s)\n', nnz(hypot(tx,ty)>0));
    end
end

[~,anchor] = min([F.acqOrder]);
fprintf('=== roi_stitch_map_260801 ===\n%s\n%d FOVs | step %g um\n', cfg.datasetPath, nF, STEP_UM);
fprintf('anchor (earliest acquired, never moves): %s\n', F(anchor).name);

% ---- lateral sites + display packing ---------------------------------------
lat = false(nF);
for a = 1:nF
    for b = 1:nF
        if a==b, continue; end
        w = max(0, min(F(a).motor(1)+F(a).halfW_um, F(b).motor(1)+F(b).halfW_um) - ...
                   max(F(a).motor(1)-F(a).halfW_um, F(b).motor(1)-F(b).halfW_um));
        h = max(0, min(F(a).motor(2)+F(a).halfH_um, F(b).motor(2)+F(b).halfH_um) - ...
                   max(F(a).motor(2)-F(a).halfH_um, F(b).motor(2)-F(b).halfH_um));
        lat(a,b) = w*h > 0;
    end
end
site = conncomp(graph(lat | lat'));
px_pack = zeros(nF,1); py_pack = zeros(nF,1);

% MIDLINE GUARD. When sites sit on opposite sides of zero they are opposite
% HEMISPHERES, and the distance between them is anatomy, not dead space. Packing
% them would slide one across the midline -- on Vglut2/0728 it moved site 1 by
% -1093 um and sites 2/3 by +2085/+2759 um, collapsing both hemispheres into one
% strip. Packing is refused in that case regardless of the flag above.
siteX = arrayfun(@(s) median(arrayfun(@(i) F(i).motor(1), find(site==s))), 1:max(site));
bilateral = any(siteX > 0) && any(siteX < 0);
if PACK_SITES && bilateral
    fprintf(['packing REFUSED: sites lie on both sides of the midline ' ...
             '(X = %s um)\n  -- that separation is anatomy; drawing true coordinates\n'], ...
        strjoin(compose('%+.0f', siteX), ', '));
end
if PACK_SITES && ~bilateral && max(site) > 1
    % Which axis to pack along is CHOSEN, not assumed. On this dataset the two
    % sites share an X range and are separated in Y, so packing along X spread the
    % map from 894 to 1679 um instead of compacting it. Both axes are tried and
    % the one giving the more compact (squarer) bounding box wins.
    [ox, by] = pack_axis(F, site, 1, PACK_GAP);
    [oy, bx] = pack_axis(F, site, 2, PACK_GAP);
    if by <= bx
        px_pack = ox;  axname = 'X';  bb = by;
    else
        py_pack = oy;  axname = 'Y';  bb = bx;
    end
    fprintf('packed %d sites along %s (longest side %.0f um; offsets are NOT saved)\n', ...
        max(site), axname, bb);
end

% ---- depth bands ------------------------------------------------------------
z = arrayfun(@(f) f.motor(3), F)';
band = floor((z - min(z))/BAND_UM) + 1;
[ub,~,band] = unique(band); nBand = numel(ub);
fprintf('%d depth bands of %g um\n\n', nBand, BAND_UM);

% ---- images, carried through the same transform as the coordinates ----------
G = cell(nF,1);
for i = 1:nF
    A = zeros(F(i).H, F(i).W);
    if ~isempty(F(i).avgPath) && isfile(F(i).avgPath)
        A = double(imread(F(i).avgPath));
        A = A - prctile(A(:),1); A = A / max(prctile(A(:),99.5),eps);
    end
    G{i} = roi_stitch_orient_260801(min(max(A,0),1), cfg);
end

%% ---- figure ------------------------------------------------------------------
fig = figure('Color','k','Name','roi_stitch_map_260801','NumberTitle','off', ...
             'WindowState','maximized','InvertHardcopy','off');
ax  = axes('Parent',fig,'Color','k'); hold(ax,'on');
set(ax,'YDir','normal','XColor','w','YColor','w','GridColor',[.3 .3 .3]);
axis(ax,'equal'); grid(ax,'on');
if cfg.rotatedMount
    xlabel(ax,'X (um)  lateral  [= +stage y]','Color','w');
    ylabel(ax,'Y (um)  rostral  [= -stage x]','Color','w');
else
    xlabel(ax,'X (um)  [= stage x]','Color','w'); ylabel(ax,'Y (um)  [= stage y]','Color','w');
end

st = struct('F',{F},'G',{G},'tx',tx,'ty',ty,'placed',placed,'nF',nF, ...
    'sel',anchor,'anchor',anchor,'step',STEP_UM,'big',BIG_MULT,'cfg',cfg, ...
    'showROI',SHOW_ROI,'showBox',SHOW_BOX,'hidden',false(nF,1),'ax',ax, ...
    'aregDz',AUTOREG_DZ,'maxShift',MAXSHIFT,'pxPack',px_pack,'pyPack',py_pack, ...
    'site',site,'band',band,'nBand',nBand,'curBand',0,'dimAlpha',DIM_ALPHA, ...
    'zrange',[min(z) max(z)]);

for i = 1:nF
    [xd,yd] = extent(st, i);
    col = depthcolor(F(i).motor(3), st.zrange);
    st.him(i) = image('Parent',ax,'XData',xd,'YData',yd, ...
        'CData',tint(G{i},col), 'AlphaData',G{i}, ...
        'ButtonDownFcn',@(~,~) pick(fig,i));
    st.hbox(i) = rectangle('Parent',ax,'Position',boxpos(st,i), ...
        'EdgeColor',col,'LineWidth',.5,'HitTest','off');
    st.hroi(i) = plot(ax, nan, nan, 'o', ...
        'MarkerSize',5,'MarkerEdgeColor',[1 1 1],'LineWidth',.6,'HitTest','off');
end
axis(ax,'tight');

% depth colourbar, so "blue = deep" is legible without reading the code
colormap(ax, neonmap(256)); caxis(ax, st.zrange);
cb = colorbar(ax,'Color','w'); cb.Label.String = 'depth z (um)   blue = deep';
cb.Label.Color = 'w';

guidata(fig, st);
for i = 1:nF, place_tile(fig,i); end
refresh_all(fig);
set(fig,'KeyPressFcn',@(~,e) onkey(fig,e), 'CloseRequestFcn',@(~,~) onquit(fig));
uiwait(fig);
end

%% ---------------------------------------------------------------------------
function [off, longside] = pack_axis(F, site, axis_, gap)
% Slide each site so they abut along axis_ (1 = X, 2 = Y), and report the longest
% side of the resulting bounding box so the caller can pick the better axis.
nF = numel(F); off = zeros(nF,1);
half = @(i,a) tern(a==1, F(i).halfW_um, F(i).halfH_um);
cursor = 0;
for s = 1:max(site)
    k  = find(site==s);
    lo = min(arrayfun(@(i) F(i).motor(axis_) - half(i,axis_), k));
    hi = max(arrayfun(@(i) F(i).motor(axis_) + half(i,axis_), k));
    off(k) = cursor - lo;
    cursor = cursor + (hi-lo) + gap;
end
other = 3 - axis_;
olo = min(arrayfun(@(i) F(i).motor(other) - half(i,other), 1:nF));
ohi = max(arrayfun(@(i) F(i).motor(other) + half(i,other), 1:nF));
longside = max(cursor - gap, ohi - olo);
end

function [xd,yd] = extent(st, i)
f = st.F(i);
cx = f.motor(1) + st.tx(i) + st.pxPack(i);
cy = f.motor(2) + st.ty(i) + st.pyPack(i);
xd = [cx-f.halfW_um, cx+f.halfW_um];
yd = [cy-f.halfH_um, cy+f.halfH_um];
end

function p = boxpos(st, i)
[xd,yd] = extent(st,i);
p = [xd(1) yd(1) diff(xd) diff(yd)];
end

function C = neonmap(n)
% Neon depth ramp on black: blue (deep) -> cyan -> green -> yellow -> magenta.
% Saturated hues so that overlapping tiles at different depths stay separable,
% which a luminance ramp cannot do once alpha-blended.
key = [0.20 0.30 1.00;    % deep   electric blue
       0.00 0.85 1.00;    %        cyan
       0.10 1.00 0.40;    %        neon green
       1.00 0.95 0.20;    %        neon yellow
       1.00 0.30 0.80];   % shallow neon magenta
C = interp1(linspace(0,1,size(key,1)), key, linspace(0,1,n));
end

function c = depthcolor(z, zr)
C = neonmap(256);
t = (z - zr(1)) / max(diff(zr), eps);
c = C(max(1,min(256, round(1 + 255*t))), :);
end

function C = tint(g, col)
C = cat(3, g*col(1), g*col(2), g*col(3));
end

function pick(fig, i)
st = guidata(fig); st.sel = i; guidata(fig, st); refresh_all(fig);
end

function move_tile(fig, dx, dy)
st = guidata(fig); i = st.sel;
if i == st.anchor
    fprintf('  %s is the ANCHOR (earliest acquired) -- it defines the frame and does not move.\n', ...
        st.F(i).name);
    return;
end
st.tx(i) = st.tx(i) + dx;  st.ty(i) = st.ty(i) + dy;  st.placed(i) = true;
guidata(fig, st); place_tile(fig, i); set_title(fig);
end

function place_tile(fig, i)
st = guidata(fig);
[xd,yd] = extent(st,i);
set(st.him(i), 'XData', xd, 'YData', yd);
set(st.hbox(i), 'Position', boxpos(st,i));
set(st.hroi(i), 'XData', st.F(i).roi_x_um + st.tx(i) + st.pxPack(i), ...
                'YData', st.F(i).roi_y_um + st.ty(i) + st.pyPack(i));
end

function refresh_all(fig)
st = guidata(fig);
for i = 1:st.nF
    isSel  = (i == st.sel);
    inBand = (st.curBand == 0) || (st.band(i) == st.curBand);
    vis    = ~st.hidden(i);
    a = st.G{i};
    if ~inBand, a = a * st.dimAlpha; end
    set(st.him(i), 'AlphaData', a, 'Visible', onoff(vis));
    if isSel, ec = [1 1 1]; lw = 2;
    elseif i == st.anchor, ec = [1 1 1]; lw = 1.2;
    else, ec = depthcolor(st.F(i).motor(3), st.zrange); lw = .5;
    end
    set(st.hbox(i), 'Visible', onoff(vis && st.showBox), 'EdgeColor', ec, 'LineWidth', lw);
    set(st.hroi(i), 'Visible', onoff(vis && st.showROI && inBand), ...
        'MarkerEdgeColor', pickc(isSel, [1 1 1], [.8 .8 .8]), ...
        'MarkerSize', pickn(isSel, 7, 4));
    if isSel
        uistack(st.him(i),'top'); uistack(st.hbox(i),'top'); uistack(st.hroi(i),'top');
    end
end
set_title(fig);
end

function set_title(fig)
st = guidata(fig); i = st.sel;
if st.curBand == 0
    bs = sprintf('ALL %d bands', st.nBand);
else
    k = st.band == st.curBand; zz = arrayfun(@(f) f.motor(3), st.F(k));
    bs = sprintf('band %d/%d  z %.0f..%.0f um  (%d tiles)', ...
        st.curBand, st.nBand, min(zz), max(zz), nnz(k));
end
if i == st.anchor, tag = '  [ANCHOR - earliest acquired, fixed]'; else, tag = ''; end
title(st.ax, sprintf([ ...
  'SELECTED: %s     z %.0f um     dx %+.1f  dy %+.1f um%s\n' ...
  'tile %d/%d (acq #%d)  |  %d moved  |  %s  |  step %.1f um (shift x%g)\n' ...
  'colour = depth (blue deep -> magenta shallow)   |   click select, arrows move, ' ...
  'g auto-reg, shift+g all, [ ] band, 0 all, n next, r reset, h hide, a ROIs, o boxes, s save, q quit'], ...
  st.F(i).name, st.F(i).motor(3), st.tx(i), st.ty(i), tag, ...
  i, st.nF, st.F(i).acqOrder, nnz(hypot(st.tx,st.ty)>0), bs, st.step, st.big), ...
  'Color','w','Interpreter','none','FontSize',9);
end

function onkey(fig, e)
st = guidata(fig);
big = any(strcmp(e.Modifier,'shift'));
s   = st.step * pickn(big, st.big, 1);
switch e.Key
    case 'leftarrow',  move_tile(fig, -s, 0);
    case 'rightarrow', move_tile(fig, +s, 0);
    case 'uparrow',    move_tile(fig, 0, +s);   % Y is up
    case 'downarrow',  move_tile(fig, 0, -s);
    case 'g'
        if big, autoreg_all(fig); else, autoreg(fig, st.sel, true); refresh_all(fig); end
    case 'leftbracket'
        st.curBand = max(1, tern(st.curBand==0, st.nBand, st.curBand-1));
        guidata(fig,st); refresh_all(fig);
    case 'rightbracket'
        st.curBand = tern(st.curBand>=st.nBand, 1, st.curBand+1);
        guidata(fig,st); refresh_all(fig);
    case '0'
        st.curBand = 0; guidata(fig,st); refresh_all(fig);
    case {'n','space'}
        st.sel = mod(st.sel, st.nF) + 1; guidata(fig,st); refresh_all(fig);
    case 'r'
        if st.sel ~= st.anchor
            st.tx(st.sel)=0; st.ty(st.sel)=0; st.placed(st.sel)=false;
            guidata(fig,st); place_tile(fig,st.sel); set_title(fig);
        end
    case 'h'
        st.hidden(st.sel) = ~st.hidden(st.sel); guidata(fig,st); refresh_all(fig);
    case 'a', st.showROI = ~st.showROI; guidata(fig,st); refresh_all(fig);
    case 'o', st.showBox = ~st.showBox; guidata(fig,st); refresh_all(fig);
    case 's', save_offsets(fig);
    case 'q', onquit(fig);
end
end

%% ---------------------------------------------------------------------------
function autoreg(fig, i, verbose)
% Phase-correlate tile i against the composite of the tiles it should trust.
%
% Reference selection is the substance here:
%   - same depth only (|dz| <= aregDz). Registering against tissue hundreds of um
%     deeper is meaningless, and phase correlation returns a confident wrong peak
%     rather than failing. Without this gate, corrections of 40-82 um appeared on
%     tiles whose true correction is a few um.
%   - EARLIER-ACQUIRED tiles preferred. Stage error accumulates, so earlier tiles
%     carry less of it; anchoring to them propagates corrections outward from the
%     best coordinate instead of chaining through drifted ones. Falls back to all
%     same-depth tiles only if nothing earlier overlaps.
%
% Phase correlation whitens magnitude and keeps only phase, so it locks onto
% structure rather than brightness (tiles differ in laser power and zoom). It
% assumes pure translation, so the result is clamped and always undoable with r.
st = guidata(fig);
if i == st.anchor, return; end

% PAIRWISE, not against a composite.
%
% Compositing every reference into one image and correlating once fails when
% overlap is partial: the canvas is built around tile i, so a partner overlapping
% only half of it is clipped and the zero-padded remainder drags the peak. On
% Sert/0721 that produced a 227 um "correction" for a pair whose true residual is
% 26 um. It never showed on Vgat/0730 because those tiles are nearly co-located.
%
% Each reference is instead correlated separately on a canvas covering the UNION
% of the two FOVs (the arrangement roi_stitch_check_axes_260801 validates), and
% the MEDIAN shift is taken. The median also means one bad reference cannot carry
% the result, which a composite average cannot promise.
[refs, usedEarlier] = pick_refs(st, i);
n = numel(refs);
if n == 0
    if verbose
        fprintf('  %-44s no tile within %g um of depth overlaps it\n', ...
            st.F(i).name(1:min(44,end)), st.aregDz);
    end
    return;
end
D = nan(n,2); PK = nan(n,1);
for k = 1:n
    [D(k,1), D(k,2), PK(k)] = pair_shift(st, i, refs(k));
end
ok = all(isfinite(D),2);
if ~any(ok)
    if verbose, fprintf('  %-44s no usable correlation\n', st.F(i).name(1:min(44,end))); end
    return;
end
dx = median(D(ok,1)); dy = median(D(ok,2)); pk = median(PK(ok)); n = nnz(ok);

% Disagreement among references is the honest warning that one of them is wrong.
if n > 1
    spread = max(hypot(D(ok,1)-dx, D(ok,2)-dy));
    if spread > st.maxShift/2 && verbose
        fprintf('  %-44s references disagree by %.0f um -- check this one by eye\n', ...
            st.F(i).name(1:min(44,end)), spread);
    end
end

if hypot(dx,dy) > st.maxShift
    if verbose
        fprintf('  %-44s REJECTED, implied shift %.1f um > %g um\n', ...
            st.F(i).name(1:min(44,end)), hypot(dx,dy), st.maxShift);
    end
    return;
end
st.tx(i) = st.tx(i)+dx; st.ty(i) = st.ty(i)+dy; st.placed(i) = true;
guidata(fig, st); place_tile(fig, i);
if verbose
    fprintf('  %-44s dx %+6.1f dy %+6.1f um  (%d %s tile(s), peak %.3f)\n', ...
        st.F(i).name(1:min(44,end)), dx, dy, n, ...
        tern(usedEarlier,'earlier','same-depth'), pk);
end
end

function [refs, earlierOnly] = pick_refs(st, i)
% Tiles that tile i should trust: same depth, and EARLIER-acquired when any such
% tile overlaps it (stage error accumulates, so earlier means less drifted).
% Falls back to all same-depth overlapping tiles only if nothing earlier does.
earlierOnly = true;
refs = scan_refs(st, i, true);
if isempty(refs)
    earlierOnly = false;
    refs = scan_refs(st, i, false);
end
end

function refs = scan_refs(st, i, earlierOnly)
refs = [];
for j = 1:st.nF
    if j==i || st.hidden(j), continue; end
    if abs(st.F(j).motor(3)-st.F(i).motor(3)) > st.aregDz, continue; end
    if earlierOnly && st.F(j).acqOrder > st.F(i).acqOrder, continue; end
    if abs(st.F(j).motor(1)+st.tx(j)-st.F(i).motor(1)-st.tx(i)) > ...
       st.F(i).halfW_um+st.F(j).halfW_um, continue; end
    if abs(st.F(j).motor(2)+st.ty(j)-st.F(i).motor(2)-st.ty(i)) > ...
       st.F(i).halfH_um+st.F(j).halfH_um, continue; end
    refs(end+1) = j; %#ok<AGROW>
end
end

function [dx, dy, pk] = pair_shift(st, i, j)
% Phase-correlate tiles i and j on a canvas covering the UNION of both, at a
% common pixel size. Returns the shift to apply to i, in um.
dx = NaN; dy = NaN; pk = NaN;
try
    px = max(st.F(i).px_um, st.F(j).px_um);
    A = imresize(st.G{i}, st.F(i).px_um/px);
    B = imresize(st.G{j}, st.F(j).px_um/px);
    aix = (st.F(i).motor(1)+st.tx(i))/px; aiy = (st.F(i).motor(2)+st.ty(i))/px;
    ajx = (st.F(j).motor(1)+st.tx(j))/px; ajy = (st.F(j).motor(2)+st.ty(j))/px;
    ha = [size(A,2) size(A,1)]/2;  hb = [size(B,2) size(B,1)]/2;
    x0 = min(aix-ha(1), ajx-hb(1)); x1 = max(aix+ha(1), ajx+hb(1));
    y0 = min(aiy-ha(2), ajy-hb(2)); y1 = max(aiy+ha(2), ajy+hb(2));
    W = round(x1-x0); H = round(y1-y0);
    if W<32||H<32||W>4000||H>4000, return; end
    CA = pasteMax(zeros(H,W), A, round(aiy-ha(2)-y0)+1, round(aix-ha(1)-x0)+1);
    CB = pasteMax(zeros(H,W), B, round(ajy-hb(2)-y0)+1, round(ajx-hb(1)-x0)+1);
    if nnz(CA)==0 || nnz(CB)==0, return; end
    CA = CA-mean(CA(:)); CB = CB-mean(CB(:));
    w = hann(H)*hann(W)';
    Rr = fft2(CA.*w) .* conj(fft2(CB.*w)); Rr = Rr./max(abs(Rr),eps);
    c  = fftshift(real(ifft2(Rr)));
    [pk,ix] = max(c(:)); [rr,cc] = ind2sub(size(c),ix);
    dx = -(cc-floor(W/2)-1)*px;  dy = -(rr-floor(H/2)-1)*px;
catch
end
end

function autoreg_all(fig)
st = guidata(fig);
% Earliest first, so each tile registers against neighbours that are themselves
% already anchored back toward the least-drifted coordinate.
[~,ord] = sort([st.F.acqOrder]);
fprintf('\nauto-registering %d tiles, earliest first (anchor fixed)...\n', st.nF);
for i = ord(:)'
    if st.hidden(i) || i == st.anchor, continue; end
    autoreg(fig, i, true);
end
refresh_all(fig);
st = guidata(fig);
d = hypot(st.tx,st.ty); d = d(d>0);
if isempty(d), fprintf('no corrections applied.\n');
else, fprintf('done: %d tiles moved, median %.1f um, max %.1f um. r resets, s saves.\n', ...
        numel(d), median(d), max(d));
end
end

function C = pasteMax(C, A, r0, c0)
[h,w] = size(A);
dr = max(r0,1):min(r0+h-1,size(C,1)); dc = max(c0,1):min(c0+w-1,size(C,2));
if isempty(dr)||isempty(dc), return; end
C(dr,dc) = max(C(dr,dc), A(dr-r0+1, dc-c0+1));
end

%% ---------------------------------------------------------------------------
function save_offsets(fig)
% tx/ty only. The display packing offset is deliberately NOT included: it exists
% to make the map readable and must never reach the coordinates used for matching.
st = guidata(fig); cfg = st.cfg;
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
fovTable = table(string({st.F.name})', st.tx, st.ty, st.placed, ...
                 arrayfun(@(f) f.motor(3), st.F)', ...
                 'VariableNames', {'fov','tx_um','ty_um','placed','z_um'});
save(cfg.offsetFile, 'fovTable', 'cfg');
d = hypot(st.tx,st.ty); d = d(d>0);
fprintf('saved %d moved tile(s)%s -> %s\n', numel(d), ...
    tern(isempty(d),'',sprintf(' (median %.1f, max %.1f um)', median(d), max(d))), ...
    cfg.offsetFile);
end

function onquit(fig)
save_offsets(fig);
fprintf('\nNext: roi_stitch_match_260801\n');
delete(fig);
end

function s = onoff(b), if b, s='on'; else, s='off'; end, end
function c = pickc(b,a,d), if b, c=a; else, c=d; end, end
function n = pickn(b,a,d), if b, n=a; else, n=d; end, end
function v = tern(c,a,b), if c, v=a; else, v=b; end, end
