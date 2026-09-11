function falloff_proj_box_260901()
%FALLOFF_PROJ_BOX_260901  Old-style side projections, boxed and power-normalised.
%
% Four panels: x-z (coronal) and y-z (sagittal) for the WT and the shiverer.
% These are the same two side views falloff_projections_260805.m draws at the
% bottom of its multi-panel figure, with two differences:
%
%  1. RESTRICTED TO THE RENDERED BOX (200 x 200 um), not the whole 456 um field,
%     so this figure and the 3D renders describe the same block of tissue.
%
%  2. POWER-NORMALISED ONLY -- divided by mW^P_EXPONENT, with NO exp(z/l) depth
%     correction. That is the deliberate difference from the old figure, which
%     applies both. The consequence is that these panels go genuinely dark with
%     depth instead of being flattened to even brightness: what is shown is how
%     much signal per unit power actually came back, which is the quantity the
%     falloff plot fits. Set NORM = 'depthcorr' for the old look.
%
% ONE SHARED DISPLAY WINDOW across all four panels. After dividing by mW^2 both
% datasets are in the same units (ADU/mW^2) with the same dye, so a common window
% compares tissue rather than exposure. That is only true because the depth
% correction is OFF -- with it on, the WT is multiplied by 446x at 200 um against
% the shiverer's 26x and no shared window is meaningful.
%
% SATURATION. Both stacks have planes with 5-10% of pixels pinned at the 32767
% ceiling; under a true MAX those voxels win every ray and smear into bright
% bands. PROJ_PCT collapses each axis by a high percentile instead. Set
% PROJ_PCT = 100 for the exact MAX the old figure uses.
%
% Runqi Zhang / 2026-09-01

%% --------------------------- USER PARAMETERS -------------------------------
% {dataDir, roiTag, label, ML_RANGE [lo hi], AP_RANGE [lo hi], powerLut}
% Same boxes as falloff_render3d_260831 / falloff_metrics_box_260901.
DS = { ...
  'C:\fall-off\260728_vglut2_vessel',   'roi1', 'Vglut2 (WT)', [100 300], [50  250], 'post_260723'
  'C:\fall-off\260804_shiver_dbh_vessel','roi2','shiverer',    [0   200], [100 300], 'post_260723' };

NORM       = 'depthcorr';  % 'depthcorr' = exp(z/l)/mW^2, the old-style correction
                           %               with the gain clamped past FIT_DEPTH.
                           %               Window is PER DATASET -- the WT is
                           %               multiplied by 446x at 200 um against
                           %               the shiverer's 26x, so a shared window
                           %               would compare corrections, not tissue.
                           % 'power'     = /mW^2 only, no exp. Then the panels
                           %               ARE in common units and share one
                           %               window, but the 450x/3400x range over
                           %               depth leaves everything below ~80 um
                           %               black.
PROJ_PCT   = 99.5;         % percentile along the collapsed axis; 100 = true MAX
CLIM_PCT   = [0.5 99.5];   % display window, SHARED across panels when NORM='power'
DEPTH_LIM  = [0 500];      % shared depth axis, matching falloff_metrics_box_260901
P_EXPONENT = 2;
PIX_TARGET = 1.0;
VOL_FROM   = -10;          % cache to read; display is still clipped to DEPTH_LIM
OUT        = 'C:\fall-off\falloff_proj_box_260901';
FS         = 0.78;
%% ---------------------------------------------------------------------------

here = fileparts(mfilename('fullpath')); addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

n = size(DS,1);  P = cell(n,1);
for k = 1:n
    P{k} = box_proj(DS{k,1}, DS{k,2}, DS{k,4}, DS{k,5}, DS{k,6}, ...
                    VOL_FROM, NORM, PROJ_PCT, P_EXPONENT, PIX_TARGET, DEPTH_LIM);
end

% Shared window: legitimate only when the panels are in common units.
shared = strcmpi(NORM,'power');
if shared
    allv = [];
    for k = 1:n, allv = [allv; P{k}.COR(:); P{k}.SAG(:)]; end %#ok<AGROW>
    allv = double(allv(isfinite(allv)));
    gLo = prctile(allv, CLIM_PCT(1));  gHi = prctile(allv, CLIM_PCT(2));
    fprintf('shared display window [%.4g %.4g]\n', gLo, gHi);
end

%% --- figure: views down, datasets across -------------------------------------
f = figure('Color','w','Units','pixels','Position',[30 30 round(1250*FS) round(1000*FS)]);
set(f,'DefaultAxesFontSize',max(5,10*FS),'DefaultTextFontSize',max(5,10*FS));
tl = tiledlayout(f, 2, n, 'TileSpacing','compact','Padding','compact');

VIEWS = { 'COR', 'x-z  coronal',  'lateral \rightarrow medial (\mum)', 'axML'
          'SAG', 'y-z  sagittal', 'antero-posterior (\mum)',           'axAP' };

for v = 1:2
    for k = 1:n
        ax = nexttile(tl, (v-1)*n + k);
        B = P{k};  I = B.(VIEWS{v,1});  xax = B.(VIEWS{v,4});
        if shared, lo = gLo; hi = gHi;
        else
            % Per DATASET, not per panel: the coronal and sagittal views are two
            % sides of one block, and scaling them independently would give the
            % same vessel two brightnesses depending on which view it fell in.
            u = double([B.COR(:); B.SAG(:)]);  u = u(isfinite(u));
            lo = prctile(u, CLIM_PCT(1)); hi = prctile(u, CLIM_PCT(2));
        end
        if hi <= lo, hi = lo + 1; end
        % rows = depth, cols = the in-plane axis -- same orientation as the old
        % figure's side views
        imagesc(ax, xax, B.depth, I.', [lo hi]);
        colormap(ax, gray(256)); set(ax,'YDir','reverse');
        axis(ax,'image');              % z on the same scale as x and y
        ylim(ax, DEPTH_LIM); xlim(ax, [xax(1) xax(end)]);
        set(ax,'YTick', DEPTH_LIM(1):100:DEPTH_LIM(2), 'Box','on');
        xlabel(ax, VIEWS{v,3});
        if k == 1, ylabel(ax,'depth (\mum)'); end
        if v == 1
            title(ax, sprintf('%s  --  %s\nML %g-%g, AP %g-%g \\mum   (%s)', ...
                DS{k,3}, DS{k,2}, DS{k,4}(1), DS{k,4}(2), DS{k,5}(1), DS{k,5}(2), ...
                VIEWS{v,2}), 'FontWeight','normal');
        else
            title(ax, VIEWS{v,2}, 'FontWeight','normal','FontSize',max(6,10*FS));
        end
    end
end

if PROJ_PCT >= 100, pStr = 'MAX'; else, pStr = sprintf('%g pct', PROJ_PCT); end
if shared
    wStr = 'ONE window shared across all four panels';
    nStr = 'no exp(z/l) -- these darken with depth because the signal does';
else
    wStr = 'window per DATASET (both views of a block share one)';
    nStr = 'depth-corrected by exp(z/l), gain clamped past the fit range';
end
title(tl, {sprintf(['side projections inside the RENDERED 200 x 200 \\mum box   |   ' ...
        'normalisation: %s   |   %s projection along the collapsed axis'], NORM, pStr), ...
        sprintf('%s   |   display %g-%g%%, %s', nStr, CLIM_PCT(1), CLIM_PCT(2), wStr)}, ...
        'FontWeight','normal','FontSize',max(6,10*FS));

exportgraphics(f,[OUT '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(f,[OUT '.pdf'],'ContentType','vector','BackgroundColor','white');
d = dir([OUT '.png']);
fprintf('\nsaved %s.png (%.0f KB) + .pdf\n', OUT, d.bytes/1024);
if d.bytes < 20000, fprintf(2,'  PNG suspiciously small -- possible blank stub\n'); end
end

% =========================================================================
function B = box_proj(dataDir, tag, mlR, apR, lut, volFrom, norm, projPct, pExp, pixT, dLim)
if volFrom >= 0
    vf = fullfile(dataDir, sprintf('falloff_proj_%s_vol.mat', tag));
else
    vf = fullfile(dataDir, sprintf('falloff_proj_%s_vol_from%d.mat', tag, volFrom));
end
assert(isfile(vf), 'no cache: %s', vf);
S = load(vf);
Ff = load(fullfile(dataDir, sprintf('falloff_%s.mat', tag)));

mW = laser_power_calibration(S.pctPlane, lut);
switch lower(norm)
case 'power'
    corr = reshape(single(1 ./ (mW(:).^pExp)), 1, 1, []);
case 'depthcorr'
    zG   = min(S.depth(:), Ff.FIT_DEPTH(2));
    corr = reshape(single(exp(zG/Ff.Lp) ./ (mW(:).^pExp)), 1, 1, []);
otherwise
    error('NORM must be ''power'' or ''depthcorr''');
end

[nR, nC, ~] = size(S.V);
rI = sort(round(mlR/pixT)+1); rI = max(1,min(nR,rI));
cI = sort(round(apR/pixT)+1); cI = max(1,min(nC,cI));
keep = S.depth >= dLim(1) & S.depth <= dLim(2);
Vol = S.V(rI(1):rI(2), cI(1):cI(2), keep) .* corr(1,1,keep);

% rot90 mount: image ROW = stage y = ML, image COL = stage x = AP.
%   coronal  = ML vs depth = collapse the AP axis
%   sagittal = AP vs depth = collapse the ML axis
B.COR = proj_pct(Vol, 2, projPct);
B.SAG = proj_pct(Vol, 1, projPct);
B.depth = S.depth(keep);
B.axML  = (0:size(B.COR,1)-1)*pixT;
B.axAP  = (0:size(B.SAG,1)-1)*pixT;
fprintf('%-12s box %dx%d px, %d planes %g..%g um, norm %s\n', tag, ...
        diff(rI)+1, diff(cI)+1, nnz(keep), B.depth(1), B.depth(end), norm);
end

% =========================================================================
function M = proj_pct(Vs, dim, p)
if p >= 100
    M = squeeze(max(Vs, [], dim));
else
    M = squeeze(prctile(double(Vs), p, dim));
end
end
