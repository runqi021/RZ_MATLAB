% Ventral_surface_ROI_xy_map_260528.m
% -----------------------------------------------------------------------
%  Plot every ROI from every FOV in shared stage XY coordinates.
%  Each ROI = one dot at its centroid converted to (motor + offset*pxSize)
%  in microns. Color by group (matches coherence polar). Sig vs non-sig
%  distinguished by marker fill (filled = sig, hollow = non-sig).
%
%  Inputs:
%   - cpSAM_output.mat (per FOV)  for maskL -> ROI centroids in px
%   - *_ch1_meta.mat   (per FOV)  for motorPosition, zoomFactor, pixelSize_um
%   - coherence_polar_data.mat    for PP, labels, confC (sig flag)
%
%  Output:
%   - ventral_surface_ROI_xy_map.pdf / .png  in coherence_polar_260528/
%   - ROI_xy_table.csv                      per-ROI (group, sig, x_um, y_um, z_um)
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));

%% ===================== USER-EDITABLE =================================
rootPath = 'D:\Ventral_surface_summary';
dataMat  = fullfile(rootPath, 'coherence_polar_260528', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'coherence_polar_260528');

xSign         = +1;        % stage sign convention
ySign         = +1;
PixelSizeBase = 1.7778;    % um/px at zoom 1 (fallback)
fallback_fps  = 30;

% group label -> color (matches coherence polar)
groups       = {'IO','ChAT','Vglut2','Vgat','Sst'};
group_colors = [0    0    0   ;
                0.85 0.10 0.10;
                0.10 0.65 0.20;
                0.10 0.30 0.85;
                0.55 0.20 0.75];

scan_dirs = {'ChAT','Vglut2','Vgat','Sst'};

% marker style
ms_sig    = 25;            % size for sig dots (filled)
ms_nonsig = 12;            % size for non-sig dots (small hollow)
edgeW     = 0.4;

doSave = true;
% =====================================================================

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% ---- load coherence results: sig flag per ROI label ----
S = load(dataMat, 'PP','labels','confC');
PP = S.PP; labels = S.labels; confC = S.confC;
sig_lookup = containers.Map(labels, num2cell(PP.r >= confC));

%% ---- walk FOVs, compute ROI stage coords ----
R = struct('group',{},'date',{},'fov',{},'rid',{},'x',{},'y',{},'z',{},'gi',{},'is_sig',{});

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), continue; end
    sam_hits = dir(fullfile(gdir,'**','*cpSAM_output.mat'));
    fprintf('\n=== [%s] %d FOVs ===\n', sname, numel(sam_hits));

    for hh = 1:numel(sam_hits)
        fp = sam_hits(hh).folder;
        rel = strrep(fp, gdir, ''); rel = regexprep(rel,'^[\\/]+','');
        parts = regexp(rel,'[\\/]','split');
        if numel(parts) < 3, continue; end
        dateStr = parts{1};
        subcell = parts{2};
        [~,n,e] = fileparts(fp); fov = [n e];

        % assign group: ChAT split into IO vs ChAT by subcell name
        gname = sname;
        if strcmpi(sname,'ChAT') && ~isempty(regexpi(subcell,'IO','once'))
            gname = 'IO';
        end
        gi = find(strcmp(groups, gname),1);
        if isempty(gi), continue; end

        % --- meta ---
        mh = dir(fullfile(fp,'*_ch1_meta.mat'));
        if isempty(mh), mh = dir(fullfile(fp,'*_meta.mat')); end
        if isempty(mh), warning('no meta for %s', fov); continue; end
        M = load(fullfile(mh(1).folder, mh(1).name));
        if ~isfield(M,'motorPosition') || numel(M.motorPosition) < 3
            warning('no motorPosition for %s', fov); continue;
        end
        motor = M.motorPosition(:).';
        px_um = NaN;
        if isfield(M,'pixelSize_um') && isfinite(M.pixelSize_um) && M.pixelSize_um>0
            px_um = M.pixelSize_um;
        elseif isfield(M,'zoomFactor') && isfinite(M.zoomFactor) && M.zoomFactor>0
            px_um = PixelSizeBase / M.zoomFactor;
        end
        if ~isfinite(px_um), warning('no pixel size for %s', fov); continue; end

        % --- cpSAM mask -> per-ROI centroids ---
        S2 = load(fullfile(sam_hits(hh).folder, sam_hits(hh).name), 'maskL');
        if ~isfield(S2,'maskL') || isempty(S2.maskL), continue; end
        maskL = S2.maskL; [H, W] = size(maskL);
        lbls  = setdiff(unique(maskL(:)), 0);
        props = regionprops(maskL, 'Centroid');
        props = props(lbls);

        for k = 1:numel(props)
            c   = props(k).Centroid;
            rid = lbls(k);
            % stage XY in microns: centroid offset from image center, scaled by px size
            x_um = motor(1) + xSign*(c(1) - W/2) * px_um;
            y_um = motor(2) + ySign*(c(2) - H/2) * px_um;
            z_um = motor(3);

            % sig lookup from coherence labels (group/date/fov/rid)
            lab = sprintf('%s/%s/%s/%d', gname, dateStr, fov, rid);
            if isKey(sig_lookup, lab)
                is_sig = sig_lookup(lab);
            else
                is_sig = false;     % not in coherence dataset (filtered by minSpikes etc.)
            end

            R(end+1) = struct('group',gname,'date',dateStr,'fov',fov, ...
                              'rid',rid,'x',abs(x_um),'y',y_um,'z',z_um, ...
                              'gi',gi,'is_sig',logical(is_sig)); %#ok<SAGROW>
        end
        fprintf('  [%-6s/%s/%-7s] %-40s  %d ROIs   motor=[%.0f %.0f %.0f]\n', ...
                gname, dateStr, subcell, fov, numel(props), motor(1), motor(2), motor(3));
    end
end

assert(~isempty(R), 'No ROIs found.');
nROI = numel(R);
fprintf('\nTotal: %d ROIs across %d FOVs\n', nROI, numel(unique({R.fov})));

%% ---- plot ----
fig = figure('Color','w','Name','ROI XY map','Units','centimeters','Position',[2 2 24 22]);
ax = axes(fig); hold(ax,'on');

% non-sig first (under sig)
for gi = 1:numel(groups)
    sel = [R.gi]==gi & ~[R.is_sig];
    if any(sel)
        col = group_colors(gi,:);
        plot(ax, [R(sel).x], [R(sel).y], 'o', ...
             'MarkerFaceColor','none', 'MarkerEdgeColor', col, ...
             'MarkerSize', sqrt(ms_nonsig), 'LineWidth', edgeW);
    end
end
% sig on top, filled
for gi = 1:numel(groups)
    sel = [R.gi]==gi & [R.is_sig];
    if any(sel)
        col = group_colors(gi,:);
        plot(ax, [R(sel).x], [R(sel).y], 'o', ...
             'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
             'MarkerSize', sqrt(ms_sig)*1.5, 'LineWidth', edgeW);
    end
end

axis(ax,'equal'); grid(ax,'on'); box(ax,'on');
xMax = max([R.x]);                    % all x folded to abs -> single side
xlim(ax, [0 xMax * 1.05]);
xlabel(ax,'|stage X| (\mum)  -- L/R folded'); ylabel(ax,'stage Y (\mum)');
title(ax, sprintf('Ventral-surface ROI map  (N=%d ROIs, %d sig, L/R folded)', nROI, sum([R.is_sig])));

% legend by group, with sig/non-sig markers
hLeg = gobjects(2*numel(groups),1);
lblLeg = cell(2*numel(groups),1);
n_per_group = zeros(1,numel(groups));
n_sig_per_group = zeros(1,numel(groups));
for gi = 1:numel(groups)
    n_per_group(gi)     = sum([R.gi]==gi);
    n_sig_per_group(gi) = sum([R.gi]==gi & [R.is_sig]);
end
for gi = 1:numel(groups)
    col = group_colors(gi,:);
    hLeg(2*gi-1) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor',col, 'MarkerEdgeColor','k', ...
                        'MarkerSize', sqrt(ms_sig)*1.5);
    hLeg(2*gi)   = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor','none','MarkerEdgeColor',col, ...
                        'MarkerSize', sqrt(ms_nonsig));
    lblLeg{2*gi-1} = sprintf('%s sig (%d)', groups{gi}, n_sig_per_group(gi));
    lblLeg{2*gi}   = sprintf('%s n.s. (%d)', groups{gi}, n_per_group(gi)-n_sig_per_group(gi));
end
legend(ax, hLeg, lblLeg, 'Location','eastoutside', 'FontSize',8);

%% ---- 3D version (XY + Z) ----
fig3 = figure('Color','w','Name','ROI XYZ map','Units','centimeters','Position',[2 2 26 22]);
ax3 = axes(fig3); hold(ax3,'on');

% non-sig first
for gi = 1:numel(groups)
    sel = [R.gi]==gi & ~[R.is_sig];
    if any(sel)
        col = group_colors(gi,:);
        plot3(ax3, [R(sel).x], [R(sel).y], [R(sel).z], 'o', ...
              'MarkerFaceColor','none', 'MarkerEdgeColor', col, ...
              'MarkerSize', sqrt(ms_nonsig), 'LineWidth', edgeW);
    end
end
% sig on top, filled
for gi = 1:numel(groups)
    sel = [R.gi]==gi & [R.is_sig];
    if any(sel)
        col = group_colors(gi,:);
        plot3(ax3, [R(sel).x], [R(sel).y], [R(sel).z], 'o', ...
              'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
              'MarkerSize', sqrt(ms_sig)*1.5, 'LineWidth', edgeW);
    end
end
grid(ax3,'on'); box(ax3,'on'); axis(ax3,'vis3d');
xlim(ax3, [0 xMax * 1.05]);    % L/R folded -> single side
xlabel(ax3,'|stage X| (\mum)'); ylabel(ax3,'stage Y (\mum)'); zlabel(ax3,'stage Z (\mum)');
title(ax3, sprintf('Ventral-surface ROI 3D map  (N=%d ROIs, %d sig, L/R folded)', nROI, sum([R.is_sig])));
view(ax3, 3);

% legend by group (same as 2D)
hL = gobjects(2*numel(groups),1); lL = cell(2*numel(groups),1);
for gi = 1:numel(groups)
    col = group_colors(gi,:);
    hL(2*gi-1) = plot3(ax3, NaN,NaN,NaN, 'o', 'MarkerFaceColor',col, 'MarkerEdgeColor','k', 'MarkerSize', sqrt(ms_sig)*1.5);
    hL(2*gi)   = plot3(ax3, NaN,NaN,NaN, 'o', 'MarkerFaceColor','none','MarkerEdgeColor',col, 'MarkerSize', sqrt(ms_nonsig));
    lL{2*gi-1} = sprintf('%s sig (%d)', groups{gi}, n_sig_per_group(gi));
    lL{2*gi}   = sprintf('%s n.s. (%d)', groups{gi}, n_per_group(gi)-n_sig_per_group(gi));
end
legend(ax3, hL, lL, 'Location','eastoutside', 'FontSize',8);

%% ---- save ----
if doSave
    exportgraphics(fig,  fullfile(outDir,'ROI_xy_map.png'),  'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig,  fullfile(outDir,'ROI_xy_map.pdf'),  'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig3, fullfile(outDir,'ROI_xyz_map.png'), 'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig3, fullfile(outDir,'ROI_xyz_map.pdf'), 'ContentType','vector','BackgroundColor','white');

    T = struct2table(R);
    writetable(T, fullfile(outDir,'ROI_xy_table.csv'));
    fprintf('Saved ROI_xy_map + ROI_xyz_map + ROI_xy_table.csv to %s\n', outDir);
end
