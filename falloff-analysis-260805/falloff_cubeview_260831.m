function falloff_cubeview_260831(dataDir, roiTag, powerLut)
% FALLOFF_CUBEVIEW_260831  Two orthogonal MIP faces of a sub-volume, drawn as
% the two visible sides of a cube.
%
% ONE PANEL. The x-z (coronal) and y-z (sagittal) max projections are pasted
% onto two orthogonal planes that share the vertical edge at the FOV CENTRE,
% and the camera looks in at that corner -- so both faces are seen at once, the
% way you would look at the corner of a block.
%
% THIS IS NOT A VOLUME RENDER. Nothing is ray-cast and no opacity is
% accumulated. It is exactly the two projections falloff_projections_260805.m
% already makes, placed in 3D. Read it that way: a bright feature on one face
% is somewhere along the whole thickness behind it, not at the corner.
%
% THE SUB-VOLUME is the quadrant reaching CUT_UM from the centre toward MEDIAL
% and CUT_UM toward ANTERIOR, over the full depth. Directions, for MOUNT rot90
% (>= 2026-07-21), from reference_stage_axis_convention + the mount rule:
%    image ROW = stage y = medio-lateral, increasing row  -> LATERAL to MEDIAL
%    image COL = stage x,  -stage x = rostral, so DECREASING col -> ANTERIOR
% Both directions are asserted against MOUNT below rather than assumed, because
% they swap on an 'old'-mount session and a silent swap would relabel the faces.
%
% THE IMAGE PROCESSING IS THE SAME as falloff_projections_260805.m, and reuses
% that script's cached volume so it cannot diverge:
%    * falloff_proj_<tag>_vol.mat   1x1x1 um isotropic, bicubic in plane
%    * depth correction  exp(min(z,FIT_DEPTH(2))/ell) ./ mW(z).^P_EXPONENT
%      with ell and FIT_DEPTH read from falloff_<tag>.mat  (roi2: 61.4, [10 500])
%    * gain CLAMPED beyond the fitted depth, or noise explodes
%    * display window = CLIM_PCT percentiles
%
% ONE SHARED DISPLAY WINDOW ACROSS BOTH FACES, unlike the multi-panel figure
% which scales every panel on its own. These two faces are adjacent sides of one
% block, so independent scaling would make the same vessel two brightnesses
% depending on which face it fell on.
%
% OUTPUT  falloff_cube_<roiTag>.png / .pdf   -- a NEW name; the multi-panel
%         falloff_proj_<roiTag>.* are not touched.
%
% Runqi Zhang / 2026-08-31

%% --------------------------- USER PARAMETERS -------------------------------
DATA_DIR   = 'C:\fall-off\260804_shiver_dbh_vessel';
ROI_TAG    = 'roi2';
POWER_LUT  = 'post_260723';
P_EXPONENT = 2;

CUT_UM     = 200;          % extent of the box from the FOV CENTRE, each way
DIR_ML     = 'lateral';    % 'lateral' | 'medial'      -- which way the ML face runs
DIR_AP     = 'posterior';  % 'posterior' | 'anterior'  -- which way the AP face runs
DEPTH_MAX  = 500;          % roi2's volume ends at 500; 600 would add empty space
CLIM_PCT   = [0.5 99.5];
PIX_TARGET = 1.0;
MOUNT      = 'rot90';      % 'old' | 'rot90'; sets which axis is ML and which AP
VIEW_AZEL  = [-45 22];     % camera. -45 puts it in the (-ML,-AP) octant, i.e.
                           % outside the corner the two faces meet at.
ZOOM       = 0.62;         % <1 pulls the camera back; 1 = block fills the axes
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(dataDir),  DATA_DIR  = dataDir;  end
if nargin >= 2 && ~isempty(roiTag),   ROI_TAG   = roiTag;   end
if nargin >= 3 && ~isempty(powerLut), POWER_LUT = powerLut; end

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

rot = strcmpi(MOUNT,'rot90');
assert(rot, ['This script only has a MEASURED anterior/medial direction for ' ...
    'MOUNT rot90. On an old-mount session the ML axis is the image COLUMN and ' ...
    'its direction was never measured -- see falloff_projections_260805.m:251.']);

%% --- cached isotropic volume -------------------------------------------------
volFile = fullfile(DATA_DIR, sprintf('falloff_proj_%s_vol.mat', ROI_TAG));
assert(isfile(volFile), ['no cached volume:\n  %s\nRun ' ...
    'falloff_projections_260805(''%s'',''%s'') once to build it.'], ...
    volFile, DATA_DIR, ROI_TAG);
S = load(volFile);
V = S.V; depth = S.depth; pctPlane = S.pctPlane;
fprintf('volume %s, depth %g-%g um, step %g\n', mat2str(size(V)), ...
        depth(1), depth(end), median(diff(depth)));

%% --- depth correction, identical to falloff_projections_260805 ---------------
fp = fullfile(DATA_DIR, sprintf('falloff_%s.mat', ROI_TAG));
assert(isfile(fp), 'no falloff_%s.mat -- needed for ell and FIT_DEPTH', ROI_TAG);
Ff = load(fp);
ell = Ff.Lp;  FIT_DEPTH = Ff.FIT_DEPTH;
fprintf('ell = %.1f um, fit %g-%g um\n', ell, FIT_DEPTH(1), FIT_DEPTH(2));

zGain   = min(depth(:), FIT_DEPTH(2));      % clamp the gain past the fit
mWplane = laser_power_calibration(pctPlane, POWER_LUT);
corr    = reshape(single(exp(zGain/ell) ./ (mWplane(:).^P_EXPONENT)), 1, 1, []);
Vc      = V .* corr;  clear V

dzNat = median(diff(depth));
if abs(dzNat - PIX_TARGET) > 1e-6
    dNew = (depth(1):PIX_TARGET:depth(end)).';
    sz = size(Vc);
    Vc = reshape(Vc, [], sz(3)).';
    Vc = interp1(depth, Vc, dNew, 'pchip');
    Vc = reshape(Vc.', [sz(1) sz(2) numel(dNew)]);
    depth = dNew;
end

%% --- the sub-volume ----------------------------------------------------------
[nR, nC, ~] = size(Vc);
nCut = round(CUT_UM / PIX_TARGET);
rC = round(nR/2);  cC = round(nC/2);
% Row index increases LATERAL -> MEDIAL (reference_stage_axis_convention: the
% image row runs along NEGATIVE stage y), so LATERAL is the LOW-index side.
if strcmpi(DIR_ML,'medial')
    rSel = rC : min(nR, rC + nCut - 1);
else
    rSel = max(1, rC - nCut + 1) : rC;
end
% Column index increases with +stage x, and on a rot90 mount -stage x = rostral,
% so increasing column is CAUDAL/POSTERIOR and decreasing is ANTERIOR.
if strcmpi(DIR_AP,'posterior')
    cSel = cC : min(nC, cC + nCut - 1);
else
    cSel = max(1, cC - nCut + 1) : cC;
end
inT  = depth >= 0 & depth <= DEPTH_MAX;
Vs   = Vc(rSel, cSel, inT);
dT   = depth(inT);
fprintf('sub-volume rows %d-%d (%s), cols %d-%d (%s) -> %s\n', ...
        rSel(1), rSel(end), DIR_ML, cSel(1), cSel(end), DIR_AP, mat2str(size(Vs)));

% rot90: cols = AP, rows = ML
COR = squeeze(max(Vs, [], 2));    % (ML, z)  max along AP  -> coronal
SAG = squeeze(max(Vs, [], 1));    % (AP, z)  max along ML  -> sagittal

% Index 0 must be the FOV CENTRE on BOTH faces, because that is the corner where
% they meet. A selection that runs TOWARD the centre rather than away from it
% must be flipped, or the shared edge would show the far side of the block.
axML = (0:size(COR,1)-1) * PIX_TARGET;
axAP = (0:size(SAG,1)-1) * PIX_TARGET;
if strcmpi(DIR_ML,'lateral'),  COR = flipud(COR); end   % rSel ran lateral -> centre
if strcmpi(DIR_AP,'anterior'), SAG = flipud(SAG); end   % cSel ran anterior -> centre

%% --- one shared display window ----------------------------------------------
allv = [COR(:); SAG(:)];  allv = double(allv(isfinite(allv)));
lo = prctile(allv, CLIM_PCT(1));  hi = prctile(allv, CLIM_PCT(2));
if hi <= lo, hi = lo + 1; end
fprintf('shared display window [%.3g %.3g] (%.4g-%.4g pct)\n', lo, hi, CLIM_PCT);
corRGB = gray_rgb(COR.', lo, hi);    % rows = depth, cols = ML
sagRGB = gray_rgb(SAG.', lo, hi);    % rows = depth, cols = AP

%% --- figure ------------------------------------------------------------------
% Half the original 1000x900. Font sizes scale with it, or the text would
% swallow the panel at the smaller canvas.
FS = 0.5;
% Portrait: the block is 200 x 200 um wide but 500 um deep, so a landscape
% canvas wastes width and crops the depth.
f = figure('Color','w','Units','pixels','Position',[40 40 round(1040*FS) round(1400*FS)]);
set(f,'DefaultAxesFontSize', max(5, 10*FS), 'DefaultTextFontSize', max(5, 10*FS));
ax = axes(f,'Position',[0.16 0.09 0.74 0.76]); hold(ax,'on');

U = axML(end);   % medial extent
Aa = axAP(end);  % anterior extent
Zt = dT(end);

% CORONAL face: the plane AP = 0, spanning ML x depth.
surface(ax, [0 U; 0 U], [0 0; 0 0], [dT(1) dT(1); Zt Zt], ...
        'CData', corRGB, 'FaceColor','texturemap','EdgeColor','none');
% SAGITTAL face: the plane ML = 0, spanning AP x depth.
surface(ax, [0 0; 0 0], [0 Aa; 0 Aa], [dT(1) dT(1); Zt Zt], ...
        'CData', sagRGB, 'FaceColor','texturemap','EdgeColor','none');

% wireframe of the block, so the two faces read as sides of one solid
e = [0 0 0; U 0 0; U Aa 0; 0 Aa 0];
plot3(ax, [e(:,1); e(1,1)], [e(:,2); e(1,2)], [e(:,3); e(1,3)], '-', ...
      'Color',[0.75 0.15 0.15], 'LineWidth',1.1*FS+0.4);
for k = 1:4
    plot3(ax, [e(k,1) e(k,1)], [e(k,2) e(k,2)], [0 Zt], '-', ...
          'Color',[0.75 0.15 0.15], 'LineWidth',1.1*FS+0.4);
end
eb = e; eb(:,3) = Zt;
plot3(ax, [eb(:,1); eb(1,1)], [eb(:,2); eb(1,2)], [eb(:,3); eb(1,3)], '-', ...
      'Color',[0.75 0.15 0.15], 'LineWidth',1.1*FS+0.4);

set(ax,'ZDir','reverse','YDir','normal');
axis(ax,'vis3d'); daspect(ax,[1 1 1]);
view(ax, VIEW_AZEL(1), VIEW_AZEL(2));
xlim(ax,[-5 U+5]); ylim(ax,[-5 Aa+5]); zlim(ax,[-5 Zt+5]);
% vis3d + equal aspect makes the block FILL the axes, which pushes the title and
% both lateral axis labels off the canvas. Pull the camera back so the whole
% block and its decorations sit inside.
camzoom(ax, ZOOM);
xlabel(ax, sprintf('centre \\rightarrow %s (\\mum)', DIR_ML));
ylabel(ax, sprintf('centre \\rightarrow %s (\\mum)', DIR_AP));
zlabel(ax,'depth (\mum)');
set(ax,'ZTick',0:100:Zt, 'XTick',0:100:U, 'YTick',0:100:Aa, ...
       'BoxStyle','full','Box','off','FontSize',max(5,10*FS));
grid(ax,'off');

% Interpreter 'none'. A Windows path is full of backslashes and the TeX parser
% reads them as control sequences, so a tex title fails to render at all
% ("Error updating Text") and the panel comes out untitled. Plain text means
% spelling out um and x rather than the TeX escapes.
title(ax, {sprintf('%s  --  %s  --  cube view, two MIP faces', DATA_DIR, ROI_TAG), ...
     sprintf(['%g x %g um from centre (%s, %s) x %g um deep   |   ' ...
              'depth-corrected exp(z/l), l = %.1f um   |   1 um isotropic   |   ' ...
              'shared %g-%g pct window   |   mount %s'], ...
              CUT_UM, CUT_UM, DIR_ML, DIR_AP, Zt, ell, CLIM_PCT(1), CLIM_PCT(2), MOUNT)}, ...
     'FontWeight','normal','FontSize',max(5,10*FS),'Interpreter','none');

if SAVE_FIG
    stem = fullfile(DATA_DIR, sprintf('falloff_cube_%s', ROI_TAG));
    exportgraphics(f, [stem '.png'], 'Resolution',300, 'BackgroundColor','white');
    exportgraphics(f, [stem '.pdf'], 'ContentType','vector','BackgroundColor','white');
    d = dir([stem '.png']);
    fprintf('saved %s.png (%.0f KB) + .pdf\n', stem, d.bytes/1024);
    if d.bytes < 20000
        fprintf(2,'  PNG suspiciously small -- check for an OpenGL blank stub\n');
    end
end
end

% =========================================================================
function RGB = gray_rgb(I, lo, hi)
%GRAY_RGB  Apply the display window explicitly and return truecolor, so both
% faces can share one window inside a single axes -- an axes has only one clim,
% and letting MATLAB scale each surface would give the same vessel two
% brightnesses depending on which face it landed on.
I = double(I);
I = (I - lo) / (hi - lo);
I = min(max(I, 0), 1);
I(~isfinite(I)) = 0;
RGB = repmat(I, 1, 1, 3);
end
