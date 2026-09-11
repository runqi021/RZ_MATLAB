function T = roi_stitch_scan_260801()
%% roi_stitch_scan_260801  Build the FOV table this pipeline works from.
% -----------------------------------------------------------------------
% Reads the dataset ONCE and caches everything the aligner and the matcher need:
% per FOV -- motorPosition, pixel size, average projection, maskL, ROI centroids
% in stitched micrometres.
%
% STANDALONE: reads only *_cpSAM_output.mat / *_meta.mat / *_AVG_for_CP.tif from
% the dataset, and writes only into its own outDir. It does not open, and cannot
% modify, anything belonging to the appearance matcher.
%
% Output: <outDir>\fov_table.mat

cfg = stitch_cfg_260801();
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end

sam = dir(fullfile(cfg.datasetPath,'**','*_cpSAM_output.mat'));
sam = sam(~contains({sam.folder}, cfg.outDirName));
sam = sam(~contains({sam.folder}, 'roi_match_out'));   % never read the other pipeline
assert(~isempty(sam), 'No *_cpSAM_output.mat under %s', cfg.datasetPath);

fprintf('=== roi_stitch_scan_260801 ===\ndataset: %s\n%d FOVs\n\n', ...
    cfg.datasetPath, numel(sam));

F = struct('name',{},'folder',{},'site',{},'motor',{},'px_um',{}, ...
           'H',{},'W',{},'avgPath',{},'maskPath',{}, ...
           'roi_lab',{},'roi_cx_px',{},'roi_cy_px',{},'roi_area',{});
for i = 1:numel(sam)
    fp = sam(i).folder;
    q  = strsplit(fp, filesep);
    nm = q{end};  st = q{end-1};

    m = dir(fullfile(fp,'*_meta.mat'));
    if isempty(m), fprintf(2,'  skip (no _meta.mat): %s\n', nm); continue; end
    S = load(fullfile(m(1).folder, m(1).name));
    if ~isfield(S,'motorPosition') || numel(S.motorPosition) < 3
        fprintf(2,'  skip (no motorPosition): %s\n', nm); continue;
    end
    px = NaN; if isfield(S,'pixelSize_um'), px = S.pixelSize_um; end
    if ~isfinite(px), fprintf(2,'  skip (no pixelSize_um): %s\n', nm); continue; end

    L = load(fullfile(fp, sam(i).name), 'maskL');
    if ~isfield(L,'maskL') || isempty(L.maskL), fprintf(2,'  skip (no maskL): %s\n', nm); continue; end
    maskL = L.maskL;  [H,W] = size(maskL);
    lab = setdiff(unique(maskL(:)), 0);
    pr  = regionprops(maskL, 'Centroid', 'Area');
    pr  = pr(lab);

    a = dir(fullfile(fp,'*_AVG_for_CP.tif'));
    avgP = ''; if ~isempty(a), avgP = fullfile(a(1).folder, a(1).name); end

    n = numel(F) + 1;
    F(n).name    = nm;            F(n).folder = fp;      F(n).site = st;
    F(n).motor   = double(S.motorPosition(1:3));
    F(n).px_um   = px;            F(n).H = H;            F(n).W = W;
    F(n).avgPath = avgP;          F(n).maskPath = fullfile(fp, sam(i).name);
    F(n).roi_lab   = lab(:);
    F(n).roi_cx_px = arrayfun(@(p) p.Centroid(1), pr(:));
    F(n).roi_cy_px = arrayfun(@(p) p.Centroid(2), pr(:));
    F(n).roi_area  = arrayfun(@(p) p.Area, pr(:));

    % Acquisition ORDER. Stage error accumulates through a session, so the
    % earliest recording is the least-drifted reference and makes the best anchor.
    % _meta.mat has no timestamp field and the raw TIFFs named in source_tif are
    % not kept in the archive, so this is the meta file's mtime -- a proxy for
    % order, not acquisition time itself. Ordering is all that is used.
    F(n).acqTime = m(1).datenum;

end
assert(~isempty(F), 'No usable FOVs.');

% ---- FOV + ROI positions in the STITCH FRAME --------------------------------
% TWO transforms, applied in order, for two different physical reasons.
%
% 1. RIG: image -> stage.  X = stage x, Y = cfg.stageYSign * stage y.
%    The y sign is MEASURED (roi_stitch_check_axes_260801). With the wrong sign,
%    overlapping FOVs land tens to hundreds of um apart and no amount of manual
%    nudging closes the map. This one DOES change matching, because it sets where
%    each ROI sits inside its own FOV.
%
% 2. PREP: stage -> anatomy.  cfg.AXMAP_ROT90 for sessions mounted on the new
%    coordinate axis (260721 onward): lateral = +stage y, rostral = -stage x.
%    A rigid rotation, so it CANNOT change which ROIs match; it is applied so the
%    map reads in anatomical axes, matching the ventral cartoon.
%
%    The FOV FOOTPRINT TRANSPOSES with it: a 320x512 px frame spans its long side
%    laterally before rotation and rostrally after. Half-extents must swap, or
%    every overlap test and every drawn tile is the wrong shape.
%
% motorRaw keeps the untouched ScanImage value, so the applied transform stays
% auditable and nothing downstream has to guess what was done.
sy = cfg.stageYSign;
if cfg.rotatedMount, axmap = cfg.AXMAP_ROT90; else, axmap = cfg.AXMAP_ID; end
isTranspose = abs(axmap(1,1)) < abs(axmap(1,2));      % off-diagonal => axes swap

for i = 1:numel(F)
    F(i).motorRaw = F(i).motor;

    % (1) rig: FOV centre and ROI offsets into stage coordinates.
    %
    % sy applies to the ROI OFFSET ONLY, never to the centre. The measured fact
    % is relative -- within a FOV, increasing row moves toward DECREASING stage y,
    % while the centre stays where the stage says it is. Flipping both is a global
    % reflection: relative geometry is unchanged, so it silently does nothing.
    % (That exact mistake reverted matching to the wrong-axis numbers,
    % 787 cells / 78 seen >1x / max 3x, and the only thing that caught it was
    % re-running the matcher.)
    cx_st = F(i).motor(1);
    cy_st = F(i).motor(2);
    dx_st = (F(i).roi_cx_px - F(i).W/2) * F(i).px_um;
    dy_st = sy * (F(i).roi_cy_px - F(i).H/2) * F(i).px_um;
    hx_st = (F(i).W/2) * F(i).px_um;
    hy_st = (F(i).H/2) * F(i).px_um;

    % (2) prep: rotate centre, ROI offsets and half-extents together
    c = axmap * [cx_st; cy_st];
    d = axmap * [dx_st.'; dy_st.'];
    h = abs(axmap) * [hx_st; hy_st];                  % |.| : extents are unsigned

    F(i).motor(1:2) = c.';
    F(i).roi_x_um   = c(1) + d(1,:).';
    F(i).roi_y_um   = c(2) + d(2,:).';
    F(i).halfW_um   = h(1);
    F(i).halfH_um   = h(2);
    F(i).frameX_um  = 2*h(1);
    F(i).frameY_um  = 2*h(2);
end

fprintf('\nstitch frame\n');
fprintf('  rig  : X = stage x, Y = %+d * stage y   (measured)\n', sy);
if cfg.rotatedMount
    fprintf('  prep : ROTATED mount -- lateral = +stage y, rostral = -stage x\n');
else
    fprintf('  prep : standard mount -- no rotation\n');
end
if isTranspose, fprintf('  -> footprint transposes (half-extents swapped)\n'); end

fprintf('\n%-44s %-7s %9s %9s %8s %13s %8s %5s\n', ...
    'recording','site','X(um)','Y(um)','z(um)','frame(um)','um/px','nROI');
for i = 1:numel(F)
    fprintf('  %-42s %-7s %9.1f %9.1f %8.1f %6.0f x %-4.0f %8.4f %5d\n', ...
        F(i).name(1:min(42,end)), F(i).site, F(i).motor(1), F(i).motor(2), ...
        F(i).motor(3), F(i).frameX_um, F(i).frameY_um, F(i).px_um, numel(F(i).roi_lab));
end

[~, ord] = sort([F.acqTime]);
for k = 1:numel(F), F(ord(k)).acqOrder = k; end
fprintf('\nacquisition order (meta mtime): %s .. %s  -> anchor = %s\n', ...
    datestr(min([F.acqTime]),'mm-dd HH:MM'), datestr(max([F.acqTime]),'mm-dd HH:MM'), ...
    F(ord(1)).name);

nROI = sum(arrayfun(@(f) numel(f.roi_lab), F));
fprintf('\n%d FOVs, %d ROIs total\n', numel(F), nROI);
fprintf('X range %.0f .. %.0f um   Y range %.0f .. %.0f um   Z range %.0f .. %.0f um\n', ...
    min(arrayfun(@(f) f.motor(1),F)), max(arrayfun(@(f) f.motor(1),F)), ...
    min(arrayfun(@(f) f.motor(2),F)), max(arrayfun(@(f) f.motor(2),F)), ...
    min(arrayfun(@(f) f.motor(3),F)), max(arrayfun(@(f) f.motor(3),F)));
fprintf('frame sizes: %.0f x %.0f um (smallest) .. %.0f x %.0f um (largest)\n', ...
    min(arrayfun(@(f) f.frameX_um,F)), min(arrayfun(@(f) f.frameY_um,F)), ...
    max(arrayfun(@(f) f.frameX_um,F)), max(arrayfun(@(f) f.frameY_um,F)));
save(fullfile(cfg.outDir,'fov_table.mat'), 'F', 'cfg');
fprintf('Saved %s\n', fullfile(cfg.outDir,'fov_table.mat'));

if nargout, T = F; end
end
