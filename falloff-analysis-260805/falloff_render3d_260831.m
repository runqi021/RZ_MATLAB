function falloff_render3d_260831()
%FALLOFF_RENDER3D_260831  Oblique ray-cast render of the falloff sub-volumes.
%
% A REAL 3D VIEW, unlike falloff_cube_pair_260831 which pastes two flat MIPs on
% two faces. Here a screen basis is built from azimuth/elevation and the volume
% is sampled with interp3 ALONG the view direction, in real microns, so vessels
% at different depths shift relative to each other and the block has parallax.
% Orthographic, not perspective: distances stay measurable.
%
% Same ray-cast pattern as map_stitch_GUI\map_subvolume_render_260530.m, which
% was built for this look. Frangi vesselness was tried there and rejected as
% looking fake, so nothing is enhanced here either.
%
% THE ACCUMULATOR IS NOT A PLAIN MAX. Both stacks have planes with 5-10% of
% pixels pinned at the 32767 digitizer ceiling. Under a MIP those voxels win
% every ray they lie on and smear into bright bands. Instead the K brightest
% samples along each ray are kept and averaged (TOP_K), so a few saturated
% voxels cannot dominate while a real vessel -- bright over many consecutive
% samples -- is unchanged. TOP_K = 1 gives a true MAX.
%
% ABOVE THE PIA. falloff_projections_260805.m:364 builds its grid as
% (0 : dzNat : max) and so discards every plane above the surface; its cached
% volume starts at depth 0. VOL_FROM selects the EXTENDED cache written by
% falloff_extend_above_260831, which prepends those planes back by reading only
% the shallowest stack's first slices (~300-1500 pages, seconds) rather than
% re-reading all four stacks. How far up is set by acquisition, not by choice:
% vglut2 roi1 was scanned from -10 um and shiverer roi2 from -50, so -10 is the
% deepest common start and +20 was never acquired for the WT at all.
%
% Everything else matches the rest of the falloff pipeline and reuses its cache:
% 1x1x1 um isotropic bicubic, exp(min(z,fit)/ell)/mW^2 with the gain clamped
% past the fitted range, ell and FIT_DEPTH read from falloff_<tag>.mat.
%
% Runqi Zhang / 2026-08-31

%% --------------------------- USER PARAMETERS -------------------------------
% {dataDir, roiTag, label, ML_RANGE [lo hi], AP_RANGE [lo hi], powerLut}
DS = { ...
  'C:\fall-off\260728_vglut2_vessel',   'roi1', 'Vglut2 (WT)', [100 300], [50  250], 'post_260723'
  'C:\fall-off\260804_shiver_dbh_vessel','roi2','shiverer',    [0   200], [100 300], 'post_260723' };

AZ         = -32;          % azimuth (deg) about the vertical -- gives the parallax
EL         =  16;          % elevation (deg): look slightly down onto the top face
TOP_K      = 4;            % ray accumulator: average the K brightest samples on
                           % each ray. K = 1 is a true MAX. Specified as a COUNT,
                           % not a percentile: there are only ~175-210 steps per
                           % ray, so 99.5 pct rounds to K = 1 and silently gives
                           % back the plain MAX the accumulator exists to avoid.
VOL_FROM   = -10;          % um; 0 = the pia-start cache, -10 = the extended one
                           % built by falloff_extend_above_260831. Only depths
                           % the stack actually acquired are available: -10 for
                           % vglut2 roi1, -50 for shiverer roi2.
RAY_STEP   = 2.0;          % um along the view direction
OUT_PX     = 1.5;          % um per output pixel
NORM_MODE  = 'perplane';   % how brightness is made comparable:
                           %  'depthcorr' exp(z/l)/mW^2, the pipeline default.
                           %              Amplifies the WT by 446x at 200 um vs
                           %              26x for the shiverer, so its deep half
                           %              becomes bright haze and the two panels
                           %              cannot share a window.
                           %  'perplane'  no depth correction; each z plane
                           %              stretched to its OWN percentiles. Kills
                           %              the brightness ramp but NOT the SNR
                           %              collapse -- a noise plane normalises to
                           %              grain, a vessel plane to vessels -- so
                           %              penetration still reads, uncoupled from
                           %              exposure. Both datasets then share one
                           %              window legitimately.
                           %  'snr'       per plane, (I - median)/(1.4826*MAD),
                           %              i.e. units of background noise SD. Same
                           %              idea as perplane but quantitative: the
                           %              display window is then in sigma.
PP_PCT     = [5 99.8];     % per-plane stretch percentiles ('perplane' only)
CLIM_PCT   = [1 99.7];     % display window; per dataset for 'depthcorr',
                           % SHARED across datasets for 'perplane' / 'snr'
                           % (they are already in common units by then)
GAMMA      = 0.85;         % <1 lifts mid-tones; the depth-corrected data is
                           % heavily skewed toward the surface
P_EXPONENT = 2;
PIX_TARGET = 1.0;
OUT        = 'C:\fall-off\falloff_render3d_260831';
FS         = 0.62;
%% ---------------------------------------------------------------------------

here = fileparts(mfilename('fullpath')); addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

n = size(DS,1);  R = cell(n,1);
for k = 1:n
    R{k} = render_one(DS{k,1}, DS{k,2}, DS{k,4}, DS{k,5}, DS{k,6}, ...
                      P_EXPONENT, PIX_TARGET, AZ, EL, TOP_K, RAY_STEP, OUT_PX, VOL_FROM, NORM_MODE, PP_PCT);
end

%% --- figure ------------------------------------------------------------------
f = figure('Color','w','Units','pixels','Position',[30 30 round(1500*FS) round(1000*FS)]);
set(f,'DefaultAxesFontSize',max(5,10*FS),'DefaultTextFontSize',max(5,10*FS));
% ONE COMMON SCREEN EXTENT for both panels. Each render is already 1.5 um/px,
% but the blocks differ in depth (276 vs 511 um), so `axis image` inside
% equally-sized axes magnifies the shorter one and the 100 um bars come out
% different lengths on the page. Forcing identical limits makes one micron the
% same distance in both, which is the whole point of putting them side by side.
uLim = [min(cellfun(@(B) B.uvec(1), R)), max(cellfun(@(B) B.uvec(end), R))];
vLim = [min(cellfun(@(B) B.vvec(1), R)), max(cellfun(@(B) B.vvec(end), R))];
% For 'perplane'/'snr' the datasets are already in common units, so ONE window
% is used for both -- which is what makes them visually comparable. For
% 'depthcorr' they are not (the WT's 99.7 pct is 658 vs 142), so each keeps its
% own and a shared one would compare exposures rather than tissue.
shared = ~strcmpi(NORM_MODE,'depthcorr');
if shared
    allv = cell2mat(cellfun(@(B) double(B.img(B.img>0)), R(:).', 'uni', 0).');
    gLo = prctile(allv, CLIM_PCT(1)); gHi = prctile(allv, CLIM_PCT(2));
end
for k = 1:n
    ax = axes(f,'Position',[0.05 + (k-1)*0.48, 0.08, 0.43, 0.74]);
    B = R{k};
    I = B.img;
    if shared
        lo = gLo; hi = gHi;
    else
        v = double(I(I > 0));
        lo = prctile(v, CLIM_PCT(1)); hi = prctile(v, CLIM_PCT(2));
    end
    if hi <= lo, hi = lo + 1; end
    In = (double(I) - lo) / (hi - lo);
    In = min(max(In, 0), 1) .^ GAMMA;
    imagesc(ax, B.uvec, B.vvec, In, [0 1]);
    colormap(ax, gray(256)); axis(ax,'image'); set(ax,'YDir','normal');
    xlim(ax, uLim); ylim(ax, vLim);      % identical -> identical magnification
    set(ax,'XTick',[],'YTick',[],'Box','off','XColor','none','YColor','none');
    hold(ax,'on');
    % box wireframe, projected through the SAME basis so it always fits the data
    plot_box(ax, B);
    % scale bar, in real microns on the screen plane (orthographic, so one
    % length is valid everywhere in the image)
    x0 = B.uvec(1) + 0.06*range(B.uvec);  y0 = B.vvec(1) + 0.07*range(B.vvec);
    plot(ax,[x0 x0+100],[y0 y0],'-','Color','w','LineWidth',3);
    text(ax,x0+50,y0,'100 \mum','Color','w','HorizontalAlignment','center', ...
         'VerticalAlignment','bottom','FontSize',max(6,9*FS));
    title(ax, sprintf('%s  --  %s\nML %g-%g, AP %g-%g \\mum, depth %g to %g \\mum, l = %.1f \\mum', ...
        DS{k,3}, DS{k,2}, DS{k,4}(1), DS{k,4}(2), DS{k,5}(1), DS{k,5}(2), B.Z0, B.Zt, B.ell), ...
        'FontWeight','normal','FontSize',max(6,11*FS));
end
if TOP_K <= 1, pStr = 'MAX'; else, pStr = sprintf('mean of top %d samples', TOP_K); end
annotation(f,'textbox',[0.02 0.92 0.96 0.06],'String', ...
    sprintf(['oblique ray-cast render (orthographic MIP, interp3 in microns)   |   az %g, el %g   |   ' ...
             'ray accumulator: %s   |   step %g um, %g um/px   |   gamma %g   |   ' ...
             'norm: %s%s   |   depth from %g um (negative = above the pia)'], ...
             AZ, EL, pStr, RAY_STEP, OUT_PX, GAMMA, NORM_MODE, ...
             string_if_shared(shared), VOL_FROM), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',max(6,10*FS), ...
    'Interpreter','none','VerticalAlignment','middle');

exportgraphics(f,[OUT '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(f,[OUT '.pdf'],'ContentType','vector','BackgroundColor','white');
d = dir([OUT '.png']);
fprintf('\nsaved %s.png (%.0f KB) + .pdf\n', OUT, d.bytes/1024);
if d.bytes < 20000, fprintf(2,'  PNG suspiciously small -- possible blank stub\n'); end
end

% =========================================================================
function B = render_one(dataDir, tag, mlR, apR, lut, pExp, pixT, az, el, topK, rayStep, outPx, volFrom, normMode, ppPct)
S = load(pick_vol(dataDir, tag, volFrom));
V = S.V; depth = S.depth; pctPlane = S.pctPlane;
Ff = load(fullfile(dataDir, sprintf('falloff_%s.mat', tag)));
ell = Ff.Lp; fitD = Ff.FIT_DEPTH;

zGain   = min(depth(:), fitD(2));
mWplane = laser_power_calibration(pctPlane, lut);
corr    = reshape(single(exp(zGain/ell) ./ (mWplane(:).^pExp)), 1, 1, []);

[nR, nC, ~] = size(V);
rI = sort(round(mlR/pixT)+1); rI = max(1,min(nR,rI));
cI = sort(round(apR/pixT)+1); cI = max(1,min(nC,cI));

switch lower(normMode)
case 'depthcorr'
    Vol = V(rI(1):rI(2), cI(1):cI(2), :) .* corr;
otherwise
    % Crop BEFORE normalising, so each plane's statistics describe the box that
    % is actually rendered and not the whole field.
    Vol = V(rI(1):rI(2), cI(1):cI(2), :);
    for z = 1:size(Vol,3)
        pl = double(Vol(:,:,z));
        switch lower(normMode)
        case 'perplane'
            q  = prctile(pl(:), ppPct);
            if q(2) <= q(1), q(2) = q(1) + 1; end
            pl = (pl - q(1)) / (q(2) - q(1));
        case 'snr'
            md = median(pl(:));
            sg = 1.4826 * median(abs(pl(:) - md));
            if sg <= 0, sg = 1; end
            pl = (pl - md) / sg;
        otherwise
            error('NORM_MODE must be depthcorr | perplane | snr');
        end
        Vol(:,:,z) = single(pl);
    end
end
clear V

% Grids in REAL microns. Row = ML, col = AP, page = depth.
mlV = (0:size(Vol,1)-1)*pixT;
apV = (0:size(Vol,2)-1)*pixT;
dpV = depth(:).';
B.Z0 = dpV(1);  B.Zt = dpV(end);  B.ell = ell;

% screen basis: right, up (toward the surface), view (into the screen)
right0 = [1 0 0]; up0 = [0 0 -1]; view0 = [0 1 0];
Raz  = rot3(deg2rad(az), [0 0 1]);
rgt  = Raz*right0(:);  vw = Raz*view0(:);  up = up0(:);
Rel  = rot3(deg2rad(el), rgt);
up   = Rel*up;  vw = Rel*vw;
rgt = rgt/norm(rgt); up = up/norm(up); vw = vw/norm(vw);

C = [mean(mlV), mean(apV), mean(dpV)];
[cx,cy,cz] = ndgrid(mlV([1 end]), apV([1 end]), dpV([1 end]));
cor = [cx(:) cy(:) cz(:)] - C;
us = cor*rgt; vs = cor*up; ws = cor*vw;
m = 6;
uvec = (min(us)-m):outPx:(max(us)+m);
vvec = (min(vs)-m):outPx:(max(vs)+m);
wvec = (min(ws)-m):rayStep:(max(ws)+m);
[UU,VV] = meshgrid(uvec,vvec);

% Top-K accumulator. K is how many of the nSteps samples the percentile keeps;
% averaging them is what stops a lone saturated voxel owning the ray.
nStep = numel(wvec);
K = max(1, min(topK, nStep));
Buf = -inf(K, numel(UU), 'single');
fprintf('%-12s %s  box %dx%dx%d um, screen %dx%d, %d ray steps, K = %d\n', ...
    tag, dataDir, numel(mlV), numel(apV), numel(dpV), numel(uvec), numel(vvec), nStep, K);
for iw = 1:nStep
    w  = wvec(iw);
    Xq = C(1) + UU*rgt(1) + VV*up(1) + w*vw(1);
    Yq = C(2) + UU*rgt(2) + VV*up(2) + w*vw(2);
    Zq = C(3) + UU*rgt(3) + VV*up(3) + w*vw(3);
    s  = interp3(mlV, apV, dpV, permute(Vol,[2 1 3]), Xq, Yq, Zq, 'linear', 0);
    Buf = sort([Buf; single(s(:).')], 1, 'descend');
    Buf = Buf(1:K, :);
end
img = reshape(mean(Buf,1), size(UU));
img(~isfinite(img)) = 0;

B.img = img; B.uvec = uvec; B.vvec = vvec;
B.C = C; B.rgt = rgt; B.up = up;
B.mlV = mlV; B.apV = apV; B.dpV = dpV;
end

% =========================================================================
function plot_box(ax, B)
%PLOT_BOX  The 12 edges of the sub-volume, through the same screen basis, so the
% wireframe cannot drift from the rendered data.
x = B.mlV([1 end]); y = B.apV([1 end]); z = B.dpV([1 end]);
[cx,cy,cz] = ndgrid(x,y,z);
P = [cx(:) cy(:) cz(:)];
E = [1 2;1 3;1 5;2 4;2 6;3 4;3 7;4 8;5 6;5 7;6 8;7 8];
for e = 1:size(E,1)
    q = (P(E(e,:),:) - B.C);
    plot(ax, q*B.rgt, q*B.up, '-', 'Color',[0.75 0.15 0.15 0.8], 'LineWidth',0.9);
end
end

function Rm = rot3(th, ax)
ax = ax(:)/norm(ax); c = cos(th); s = sin(th); t = 1-c;
x=ax(1); y=ax(2); z=ax(3);
Rm = [t*x*x+c, t*x*y-s*z, t*x*z+s*y;
      t*x*y+s*z, t*y*y+c, t*y*z-s*x;
      t*x*z-s*y, t*y*z+s*x, t*z*z+c];
end

% =========================================================================
function fp = pick_vol(dataDir, tag, volFrom)
%PICK_VOL  The extended cache when one is asked for, the pia-start cache at 0.
% Errors rather than silently falling back: a figure captioned "from -10" that
% quietly started at the pia would be wrong in a way nothing on it would show.
if volFrom >= 0
    fp = fullfile(dataDir, sprintf('falloff_proj_%s_vol.mat', tag));
else
    fp = fullfile(dataDir, sprintf('falloff_proj_%s_vol_from%d.mat', tag, volFrom));
    assert(isfile(fp), ['no cache starting at %g um for %s -- run ' ...
        'falloff_extend_above_260831(''%s'',''%s'',%g) first.'], ...
        volFrom, tag, dataDir, tag, volFrom);
end
end

function s = string_if_shared(tf)
if tf, s = ' (shared window across datasets)'; else, s = ' (per-dataset window)'; end
end
