function falloff_cube_pair_260831()
%FALLOFF_CUBE_PAIR_260831  Two cube views side by side: WT next to shiverer.
%
% Each panel is the same object as falloff_cubeview_260831 -- the x-z and y-z
% MAX projections of a sub-volume, pasted onto two orthogonal planes that meet
% at a shared vertical edge, viewed from outside that corner. NOT a ray-cast
% volume render; nothing is accumulated with opacity.
%
% RANGES ARE IN THE COORDINATES OF THE EXISTING falloff_proj_<tag> FIGURE, i.e.
% the 'lateral -> medial' and 'antero-posterior' axes printed there, so a box can
% be read straight off that figure. Each range is [CORNER, FAR]: the first value
% is the edge where the two faces meet, the second is where the block ends. So
% [300 100] and [100 300] are the same slab viewed from opposite corners.
%
% ONE SHARED DEPTH AXIS across both panels, because the comparison IS depth.
% The WT stack ends at 265 um and the shiverer at 500, and drawing them on
% separate axes would hide exactly the difference the figure is about. Each
% block is only as tall as its own data.
%
% DISPLAY WINDOWS ARE PER-DATASET, not shared. Different animals, different
% laser powers, different attenuation lengths (24.6 vs 61.4 um) -- a common
% window would say more about the exposure than the tissue. Within one panel the
% two faces DO share a window, since they are sides of one block.
%
% SATURATION. Both stacks have planes with 5-10% of pixels pinned at the 32767
% digitizer ceiling (a clot, or dye pooling in a surface vessel). A MAX
% projection lets those voxels win every ray they lie on, so a few bad voxels
% smear into a bright band across the whole face. PROJ_PCT replaces the max with
% a high percentile along the projected axis: a handful of saturated voxels can
% no longer dominate, while a real vessel -- bright over many consecutive voxels
% -- is unchanged. Set PROJ_PCT = 100 for a true MAX, which is what
% falloff_projections_260805.m does.
%
% Everything else matches falloff_projections_260805.m and reuses its cached
% volume: 1x1x1 um isotropic bicubic, depth correction exp(min(z,fit))/mW^2 with
% the gain clamped past the fitted range, ell and FIT_DEPTH from falloff_<tag>.
%
% Runqi Zhang / 2026-08-31

%% --------------------------- USER PARAMETERS -------------------------------
% {dataDir, roiTag, label, ML_RANGE [corner far], AP_RANGE [corner far], powerLut}
% ML_RANGE is on the 'lateral -> medial' axis; AP_RANGE on 'antero-posterior'.
DS = { ...
  'C:\fall-off\260728_vglut2_vessel', 'roi1', 'Vglut2 (WT)',   [100 300], [50  250], 'post_260723'
  'C:\fall-off\260804_shiver_dbh_vessel','roi2','shiverer',    [0   200], [100 300], 'post_260723' };

PROJ_PCT   = 99.5;         % percentile along the projected axis; 100 = true MAX
CLIM_PCT   = [0.5 99.5];   % display window, per dataset
P_EXPONENT = 2;
PIX_TARGET = 1.0;
MOUNT      = 'rot90';
VIEW_AZEL  = [-45 22];
ZOOM       = 0.68;
FS         = 0.62;         % overall size factor
OUT        = 'C:\fall-off\falloff_cube_pair_260831';
%% ---------------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));  addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
assert(strcmpi(MOUNT,'rot90'), ['Only rot90 has a measured ML direction -- see ' ...
    'falloff_projections_260805.m:251']);

n = size(DS,1);
P = cell(n,1);  zMax = 0;
for k = 1:n
    P{k} = one_block(DS{k,1}, DS{k,2}, DS{k,4}, DS{k,5}, DS{k,6}, ...
                     P_EXPONENT, PIX_TARGET, CLIM_PCT, PROJ_PCT);
    zMax = max(zMax, P{k}.Zt);
end

%% --- figure ------------------------------------------------------------------
f = figure('Color','w','Units','pixels', ...
           'Position',[30 30 round(1500*FS) round(1150*FS)]);
set(f,'DefaultAxesFontSize',max(5,10*FS),'DefaultTextFontSize',max(5,10*FS));

for k = 1:n
    ax = axes(f,'Position',[0.07 + (k-1)*0.47, 0.10, 0.40, 0.72]); hold(ax,'on');
    B = P{k};
    U = B.U;  Aa = B.Aa;  Zt = B.Zt;

    % coronal face on the plane AP = 0 (ML x depth)
    surface(ax, [0 U; 0 U], [0 0; 0 0], [0 0; Zt Zt], ...
            'CData', B.corRGB, 'FaceColor','texturemap','EdgeColor','none');
    % sagittal face on the plane ML = 0 (AP x depth)
    surface(ax, [0 0; 0 0], [0 Aa; 0 Aa], [0 0; Zt Zt], ...
            'CData', B.sagRGB, 'FaceColor','texturemap','EdgeColor','none');

    e = [0 0 0; U 0 0; U Aa 0; 0 Aa 0];
    lw = 1.1*FS + 0.4;
    plot3(ax,[e(:,1);e(1,1)],[e(:,2);e(1,2)],[e(:,3);e(1,3)],'-','Color',[.75 .15 .15],'LineWidth',lw);
    for q = 1:4
        plot3(ax,[e(q,1) e(q,1)],[e(q,2) e(q,2)],[0 Zt],'-','Color',[.75 .15 .15],'LineWidth',lw);
    end
    eb = e; eb(:,3) = Zt;
    plot3(ax,[eb(:,1);eb(1,1)],[eb(:,2);eb(1,2)],[eb(:,3);eb(1,3)],'-','Color',[.75 .15 .15],'LineWidth',lw);

    set(ax,'ZDir','reverse'); axis(ax,'vis3d'); daspect(ax,[1 1 1]);
    view(ax, VIEW_AZEL(1), VIEW_AZEL(2));
    % SHARED depth limit, so the two blocks are on one scale and the WT block is
    % visibly the shorter one.
    xlim(ax,[-5 U+5]); ylim(ax,[-5 Aa+5]); zlim(ax,[-5 zMax+5]);
    camzoom(ax, ZOOM);
    xlabel(ax,'lateral \rightarrow medial (\mum)');
    ylabel(ax,'antero-posterior (\mum)');
    if k == 1, zlabel(ax,'depth (\mum)'); end
    set(ax,'ZTick',0:100:zMax,'XTick',0:100:U,'YTick',0:100:Aa, ...
           'Box','off','FontSize',max(5,10*FS));
    title(ax, sprintf('%s  --  %s\nML %g-%g, AP %g-%g \\mum, l = %.1f \\mum', ...
        DS{k,3}, DS{k,2}, DS{k,4}(1), DS{k,4}(2), DS{k,5}(1), DS{k,5}(2), B.ell), ...
        'FontWeight','normal','FontSize',max(6,11*FS));
end

if PROJ_PCT >= 100, pStr = 'MAX'; else, pStr = sprintf('%g pct', PROJ_PCT); end
annotation(f,'textbox',[0.02 0.93 0.96 0.06],'String', ...
    sprintf(['cube views: two MIP faces meeting at the near edge   |   %s projection ' ...
             'along each axis   |   depth-corrected exp(z/l), gain clamped past the fit   ' ...
             '|   1 um isotropic   |   display %g-%g pct per dataset   |   shared depth axis'], ...
             pStr, CLIM_PCT(1), CLIM_PCT(2)), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',max(6,10*FS), ...
    'Interpreter','none','VerticalAlignment','middle');

exportgraphics(f,[OUT '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(f,[OUT '.pdf'],'ContentType','vector','BackgroundColor','white');
d = dir([OUT '.png']);
fprintf('\nsaved %s.png (%.0f KB) + .pdf\n', OUT, d.bytes/1024);
if d.bytes < 20000, fprintf(2,'  PNG suspiciously small -- possible OpenGL blank stub\n'); end
end

% =========================================================================
function B = one_block(dataDir, tag, mlR, apR, lut, pExp, pixT, climPct, projPct)
volFile = fullfile(dataDir, sprintf('falloff_proj_%s_vol.mat', tag));
assert(isfile(volFile), 'no cached volume: %s', volFile);
S = load(volFile);  V = S.V; depth = S.depth; pctPlane = S.pctPlane;

Ff = load(fullfile(dataDir, sprintf('falloff_%s.mat', tag)));
ell = Ff.Lp;  fitD = Ff.FIT_DEPTH;

zGain   = min(depth(:), fitD(2));
mWplane = laser_power_calibration(pctPlane, lut);
corr    = reshape(single(exp(zGain/ell) ./ (mWplane(:).^pExp)), 1, 1, []);
Vc      = V .* corr;  clear V

dzNat = median(diff(depth));
if abs(dzNat - pixT) > 1e-6
    dNew = (depth(1):pixT:depth(end)).';
    sz = size(Vc); Vc = reshape(Vc,[],sz(3)).';
    Vc = interp1(depth, Vc, dNew, 'pchip');
    Vc = reshape(Vc.',[sz(1) sz(2) numel(dNew)]);  depth = dNew;
end

% Ranges are on the printed axes, which start at 0 for pixel index 1.
[nR, nC, ~] = size(Vc);
rI = sort(round(mlR/pixT) + 1);  rI = max(1,min(nR, rI));   % rows = ML
cI = sort(round(apR/pixT) + 1);  cI = max(1,min(nC, cI));   % cols = AP
Vs = Vc(rI(1):rI(2), cI(1):cI(2), :);
dT = depth;

COR = proj_pct(Vs, 2, projPct);    % (ML, z)  collapse AP -> coronal
SAG = proj_pct(Vs, 1, projPct);    % (AP, z)  collapse ML -> sagittal

% Index 0 of each face must be the CORNER the two meet at. mlR/apR are given
% [corner far], so a range whose corner is the HIGH coordinate must be flipped.
if mlR(1) > mlR(2), COR = flipud(COR); end
if apR(1) > apR(2), SAG = flipud(SAG); end

allv = double([COR(:); SAG(:)]);  allv = allv(isfinite(allv));
lo = prctile(allv, climPct(1));  hi = prctile(allv, climPct(2));
if hi <= lo, hi = lo + 1; end

B.corRGB = gray_rgb(COR.', lo, hi);
B.sagRGB = gray_rgb(SAG.', lo, hi);
B.U  = (size(COR,1)-1)*pixT;
B.Aa = (size(SAG,1)-1)*pixT;
B.Zt = dT(end);
B.ell = ell;
fprintf('%-34s %s  rows %d-%d, cols %d-%d, depth 0-%g um, window [%.3g %.3g]\n', ...
        tag, dataDir, rI(1), rI(2), cI(1), cI(2), B.Zt, lo, hi);
end

% =========================================================================
function M = proj_pct(Vs, dim, p)
%PROJ_PCT  Collapse one axis by a high percentile instead of the max, so a few
% saturated voxels cannot win every ray. p = 100 gives the true MAX.
if p >= 100
    M = squeeze(max(Vs, [], dim));
else
    M = squeeze(prctile(double(Vs), p, dim));
end
end

function RGB = gray_rgb(I, lo, hi)
I = double(I);  I = (I - lo) / (hi - lo);
I = min(max(I,0),1);  I(~isfinite(I)) = 0;
RGB = repmat(I,1,1,3);
end
