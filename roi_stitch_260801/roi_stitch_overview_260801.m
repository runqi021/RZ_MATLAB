function roi_stitch_overview_260801()
%% roi_stitch_overview_260801  Overview map of every FOV, raw and stitched.
% -----------------------------------------------------------------------
% Two figures, written to <outDir>:
%
%   fov_overview_map      layout of all FOVs -- raw stage vs stitched, the
%                         correction applied to each, depth profile, and the
%                         ROI map coloured by how many times each cell was seen
%   fov_stitched_montage  the actual average projections composited, ONE PANEL
%                         PER DEPTH SLAB
%
% WHY PER SLAB. Z spans 330 um over 21 distinct planes here. Compositing all of
% that into one image would overlay tissue hundreds of microns apart and show a
% smear that means nothing. FOVs are grouped into slabs by single-linkage on z
% with a 20 um gap, so each panel is a set of planes that genuinely image the
% same cells.
%
% Runs before or after alignment: with no fov_offsets.mat it shows the raw layout
% and says so. Read-only apart from its own two figures.

%% ===================== USER-EDITABLE =====================
% Figure 2 panels = one lateral SITE x one depth BAND.
% Sites come from footprint overlap, so spatially separate sites (roi1 vs roi2,
% ~1200 um apart here) never share a panel. Depth is then cut into fixed-width
% bands rather than by gap-chaining: a finely sampled depth series has every
% consecutive gap below any sensible threshold, so single-linkage on z walks the
% whole series into one slab (it merged 80 um and 12 FOVs before this).
SLAB_SPAN_UM  = 25;      % max z span within one panel
LABEL_FOVS    = true;    % draw FOV names on the layout panels (off = less clutter)
MONTAGE_MAXPX = 1600;    % long side of each slab composite, px
SAVE_PDF      = true;    % vector PDF alongside the PNG
ROI_MARKER    = 4;       % ROI marker size on the montage
FIG1_INCHES   = [17 10];
%% =========================================================

cfg = stitch_cfg_260801();
tf  = fullfile(cfg.outDir,'fov_table.mat');
assert(isfile(tf), 'Run roi_stitch_scan_260801 first (%s missing).', tf);
S = load(tf); F = S.F; nF = numel(F);

% ---- offsets ----------------------------------------------------------------
tx = zeros(nF,1); ty = zeros(nF,1); placed = false(nF,1); haveOff = false;
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
        haveOff = true;
    end
end

% ---- ROI positions in the current frame -------------------------------------
X=[];Y=[];
for i=1:nF, X=[X;F(i).roi_x_um+tx(i)]; Y=[Y;F(i).roi_y_um+ty(i)]; end %#ok<AGROW>

% ---- grouping, if the matcher has run and is STILL VALID ---------------------
% The grouping comes from THIS pipeline's matcher (roi_stitch_results.mat, fields
% grp and cell_size) -- never from the appearance matcher's roi_match_out_* or
% cell_link.mat.
%
% A count check alone is not enough to trust it. Change the axis convention or the
% mount rotation, re-run the scan but not the matcher, and the ROI count is
% unchanged -- so stale grouping would be plotted on new coordinates while the log
% cheerfully says "loaded". That happened during development. The matcher saves
% the coordinates it actually used, so compare against them and refuse if they
% have moved by more than the offsets can explain.
grpSize = []; grpState = 'not run yet';
if isfile(cfg.resultsFile)
    R = load(cfg.resultsFile);
    if numel(R.roi_fov) ~= numel(X)
        grpState = sprintf('STALE - %d ROIs in results vs %d now; re-run roi_stitch_match_260801', ...
            numel(R.roi_fov), numel(X));
    else
        drift = max(hypot(R.roi_x_um - X, R.roi_y_um - Y));
        if drift > 1e-6
            grpState = sprintf(['STALE - ROI coordinates moved up to %.1f um since the ' ...
                'matcher ran; re-run roi_stitch_match_260801'], drift);
        else
            grpSize = R.cell_size(R.grp);
            grpState = 'loaded';
        end
    end
end

z  = arrayfun(@(f) f.motor(3), F)';
zr = [min(z) max(z)]; if diff(zr)==0, zr = zr + [-1 1]; end
cm = turbo(256);
zcol = @(v) cm(max(1,min(256, round(1 + 255*(v-zr(1))/diff(zr)))), :);

fprintf('=== roi_stitch_overview_260801 ===\n');
fprintf('%d FOVs | offsets: %s | grouping: %s\n', nF, ...
    tern(haveOff, sprintf('%d placed', nnz(placed)), 'NONE (raw stage layout)'), ...
    grpState);
if isempty(grpSize) && ~strcmp(grpState,'not run yet')
    fprintf('!! %s\n', grpState);
end

%% ============================ FIGURE 1 =====================================
f1 = figure('Color','w','Units','inches','Position',[0.5 0.5 FIG1_INCHES], ...
            'Name','fov_overview_map','NumberTitle','off');
t = tiledlayout(f1, 2, 3, 'Padding','compact','TileSpacing','compact');

% -- (1) raw stage layout ------------------------------------------------------
ax = nexttile(t,1); hold(ax,'on');
draw_footprints(ax, F, zeros(nF,1), zeros(nF,1), zcol, LABEL_FOVS);
title(ax, sprintf('RAW stage layout (%d FOVs)', nF)); fmt(ax);

% -- (2) stitched layout -------------------------------------------------------
ax = nexttile(t,2); hold(ax,'on');
draw_footprints(ax, F, tx, ty, zcol, LABEL_FOVS);
if haveOff
    title(ax, sprintf('STITCHED (%d/%d placed, max corr %.1f um)', ...
        nnz(placed), nF, max(hypot(tx,ty))));
else
    title(ax, 'STITCHED -- no fov\_offsets.mat yet (identical to raw)');
end
fmt(ax);
cb = colorbar(ax); cb.Label.String = 'z (um)';
caxis(ax, zr); colormap(ax, turbo);

% -- (3) correction applied ----------------------------------------------------
ax = nexttile(t,3); hold(ax,'on');
if haveOff && any(hypot(tx,ty) > 0)
    [~,o] = sort(hypot(tx,ty),'descend');
    barh(ax, hypot(tx(o),ty(o)), 'FaceColor',[.3 .5 .8], 'EdgeColor','none');
    set(ax,'YTick',1:nF,'YTickLabel',shortname({F(o).name}),'FontSize',6, ...
           'YDir','reverse');
    xlabel(ax,'stitch correction |d| (um)');
    title(ax, sprintf('correction per FOV (median %.1f um)', median(hypot(tx,ty))));
else
    text(ax,.5,.5,{'no stitch applied yet','run roi\_stitch\_align\_260801'}, ...
        'HorizontalAlignment','center'); axis(ax,'off');
    title(ax,'correction per FOV');
end
box(ax,'on');

% -- (4) depth profile ---------------------------------------------------------
ax = nexttile(t,4); hold(ax,'on');
[zs,o] = sort(z);
for i = 1:nF
    plot(ax, zs(i), i, 'o', 'MarkerSize',7, 'MarkerFaceColor',zcol(zs(i)), ...
        'MarkerEdgeColor','k','LineWidth',.3);
end
set(ax,'YTick',1:nF,'YTickLabel',shortname({F(o).name}),'FontSize',6);
xlabel(ax,'motor z (um)'); grid(ax,'on'); box(ax,'on');
title(ax, sprintf('depth: %d FOVs over %.0f um', nF, range(z)));

% -- (5) ROI map ---------------------------------------------------------------
ax = nexttile(t,5); hold(ax,'on');
draw_footprints(ax, F, tx, ty, zcol, false);
if isempty(grpSize)
    plot(ax, X, Y, '.', 'Color',[.85 .2 .2], 'MarkerSize',4);
    title(ax, sprintf('%d ROIs -- %s', numel(X), grpState));
else
    s1 = grpSize == 1;
    plot(ax, X(s1),  Y(s1),  '.', 'Color',[.7 .7 .7], 'MarkerSize',4);
    sc = scatter(ax, X(~s1), Y(~s1), 14, grpSize(~s1), 'filled');
    sc.MarkerEdgeColor = 'k'; sc.LineWidth = .2;
    colormap(ax, flipud(hot(8))); caxis(ax,[1 max(grpSize)]);
    c2 = colorbar(ax); c2.Label.String = 'times cell was seen';
    % Count CELLS, not ROIs: nnz(~s1) is the number of ROIs belonging to a
    % multi-observation cell, which is a different (larger) number and read as a
    % contradiction against the summary panel.
    title(ax, sprintf('%d ROIs in %d cells | %d cells seen >1x (%d ROIs, grey = seen once)', ...
        numel(X), numel(unique(R.grp)), nnz(R.cell_size>1), nnz(~s1)));
end
fmt(ax);

% -- (6) summary ---------------------------------------------------------------
ax = nexttile(t,6); axis(ax,'off');
L = { sprintf('dataset: %s', cfg.datasetPath), '', ...
      sprintf('FOVs                 %d', nF), ...
      sprintf('ROIs                 %d', numel(X)), ...
      sprintf('z range              %.0f .. %.0f um (span %.0f)', min(z), max(z), range(z)), ...
      sprintf('pixel size           %.3f .. %.3f um/px', ...
              min(arrayfun(@(f) f.px_um,F)), max(arrayfun(@(f) f.px_um,F))), '' };
if haveOff
    L{end+1} = sprintf('stitched             %d / %d FOVs', nnz(placed), nF);
    L{end+1} = sprintf('correction  median   %.1f um', median(hypot(tx,ty)));
    L{end+1} = sprintf('            max      %.1f um', max(hypot(tx,ty)));
else
    L{end+1} = 'stitched             NOT YET -- run roi_stitch_align_260801';
end
if ~isempty(grpSize)
    L{end+1} = ''; L{end+1} = sprintf('cells                %d', numel(unique(R.grp)));
    L{end+1} = sprintf('seen >1x             %d', nnz(R.cell_size>1));
    L{end+1} = sprintf('largest group        %dx', max(R.cell_size));
end
text(ax, 0, 1, L, 'VerticalAlignment','top','FontName','Consolas','FontSize',9, ...
     'Interpreter','none');
title(ax,'summary');

title(t, sprintf('FOV overview -- %s', cfg.datasetPath), 'Interpreter','none', ...
      'FontWeight','bold');
exportgraphics(f1, fullfile(cfg.outDir,'fov_overview_map.png'), 'Resolution',200);
if SAVE_PDF
    exportgraphics(f1, fullfile(cfg.outDir,'fov_overview_map.pdf'), 'ContentType','vector');
end
fprintf('Saved fov_overview_map.png%s\n', tern(SAVE_PDF,'/.pdf',''));

%% ============================ FIGURE 2 =====================================
% Panels = lateral site x depth band (see header).
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

slab = zeros(nF,1); k = 0;
for s = unique(site)
    idx = find(site == s);
    band = floor((z(idx) - min(z(idx))) / SLAB_SPAN_UM);   % fixed-width, no chaining
    for b = unique(band(:))'
        k = k + 1; slab(idx(band == b)) = k;
    end
end
nS = k;
fprintf('%d lateral sites -> %d panels (site x %g um depth band)\n', ...
    max(site), nS, SLAB_SPAN_UM);

nc = ceil(sqrt(nS)); nr = ceil(nS/nc);
f2 = figure('Color','w','Units','inches','Position',[0.5 0.5 min(20,5*nc) min(12,4.5*nr)], ...
            'Name','fov_stitched_montage','NumberTitle','off');
t2 = tiledlayout(f2, nr, nc, 'Padding','compact','TileSpacing','compact');

for s = 1:nS
    ax = nexttile(t2, s); hold(ax,'on');
    idx = find(slab == s);
    [C, x0, y0, px] = composite(F, idx, tx, ty, MONTAGE_MAXPX);
    if isempty(C), axis(ax,'off'); continue; end
    imagesc(ax, [x0 x0+size(C,2)*px], [y0 y0+size(C,1)*px], C);
    colormap(ax, gray); axis(ax,'image'); set(ax,'YDir','normal');   % Y up
    for i = idx(:)'
        rectangle(ax,'Position',[F(i).motor(1)+tx(i)-F(i).halfW_um, ...
                                 F(i).motor(2)+ty(i)-F(i).halfH_um, ...
                                 2*F(i).halfW_um, 2*F(i).halfH_um], ...
                  'EdgeColor',zcol(F(i).motor(3)),'LineWidth',1);
        plot(ax, F(i).roi_x_um+tx(i), F(i).roi_y_um+ty(i), 'o', ...
             'MarkerSize',ROI_MARKER,'MarkerEdgeColor',[1 .85 .1],'LineWidth',.7);
    end
    nR = sum(arrayfun(@(i) numel(F(i).roi_lab), idx));
    title(ax, sprintf('site %d  |  z %.0f .. %.0f um   (%d FOV, %d ROI)', ...
        site(idx(1)), min(z(idx)), max(z(idx)), numel(idx), nR), 'FontSize',9);
    [xl,yl] = framelabels(); xlabel(ax,xl); ylabel(ax,yl); set(ax,'FontSize',7);
end
title(t2, sprintf('Stitched average projections -- lateral site x depth band %s', ...
    tern(haveOff,'(STITCHED)','(RAW STAGE -- not yet stitched)')), 'FontWeight','bold');
exportgraphics(f2, fullfile(cfg.outDir,'fov_stitched_montage.png'), 'Resolution',200);
if SAVE_PDF
    exportgraphics(f2, fullfile(cfg.outDir,'fov_stitched_montage.pdf'), 'ContentType','vector');
end
fprintf('Saved fov_stitched_montage.png%s\n\nAll written to %s\n', ...
    tern(SAVE_PDF,'/.pdf',''), cfg.outDir);
end

%% ---------------------------------------------------------------------------
function draw_footprints(ax, F, tx, ty, zcol, lab)
for i = 1:numel(F)
    x = F(i).motor(1)+tx(i); y = F(i).motor(2)+ty(i);
    rectangle(ax,'Position',[x-F(i).halfW_um, y-F(i).halfH_um, ...
                             2*F(i).halfW_um, 2*F(i).halfH_um], ...
        'EdgeColor',zcol(F(i).motor(3)),'LineWidth',1.1);
    if lab
        sn = shortname({F(i).name});     % cannot index a call result directly
        text(ax, x, y, sprintf('%s\nz%.0f', sn{1}, F(i).motor(3)), ...
            'FontSize',5,'HorizontalAlignment','center','Color',[.15 .15 .15], ...
            'Interpreter','none');
    end
end
end

function fmt(ax)
% Y UP everywhere. Y is rostral in this frame, so up-is-rostral matches the
% ventral cartoon (which also uses YDir='normal'). It is also the setting that
% pairs correctly with roi_stitch_orient_260801, whose output has row index
% increasing WITH Y.
axis(ax,'equal'); set(ax,'YDir','normal'); grid(ax,'on'); box(ax,'on');
[xl,yl] = framelabels();
xlabel(ax,xl); ylabel(ax,yl); set(ax,'FontSize',8);
end

function [xl,yl] = framelabels()
% Label the axes for whatever frame the scan actually produced, rather than
% hard-coding "stage x/y" -- on a rotated mount those are the wrong names.
cfg = stitch_cfg_260801();
if cfg.rotatedMount
    xl = 'X (um)  lateral  [= +stage y]';
    yl = 'Y (um)  rostral  [= -stage x]';
else
    xl = 'X (um)  [= stage x]';
    yl = 'Y (um)  [= stage y]';
end
end

function s = shortname(c)
s = cellfun(@(n) regexprep(regexprep(n,'_\d{5}$',''),'_\d+f_[\d.]+lp',''), ...
            c, 'UniformOutput', false);
end

function [C, x0, y0, px] = composite(F, idx, tx, ty, maxpx)
% Max-blend the average projections of one slab onto a common grid.
C=[]; x0=0; y0=0; px=1;
if isempty(idx), return; end
x0 = min(arrayfun(@(i) F(i).motor(1)+tx(i)-F(i).halfW_um, idx));
x1 = max(arrayfun(@(i) F(i).motor(1)+tx(i)+F(i).halfW_um, idx));
y0 = min(arrayfun(@(i) F(i).motor(2)+ty(i)-F(i).halfH_um, idx));
y1 = max(arrayfun(@(i) F(i).motor(2)+ty(i)+F(i).halfH_um, idx));
px = max(min(arrayfun(@(i) F(i).px_um, idx)), max(x1-x0,y1-y0)/maxpx);
W = round((x1-x0)/px); H = round((y1-y0)/px);
if W<4 || H<4 || W>6000 || H>6000, return; end
C = zeros(H,W);
for i = idx(:)'
    if isempty(F(i).avgPath) || ~isfile(F(i).avgPath), continue; end
    A = double(imread(F(i).avgPath));
    A = A - prctile(A(:),1); A = A / max(prctile(A(:),99.5),eps); A = min(max(A,0),1);
    A = roi_stitch_orient_260801(A);      % same transform as the coordinates
    As = imresize(A, F(i).px_um/px);
    [h,w] = size(As);
    r0 = round((F(i).motor(2)+ty(i)-F(i).halfH_um - y0)/px)+1;
    c0 = round((F(i).motor(1)+tx(i)-F(i).halfW_um - x0)/px)+1;
    dr = max(r0,1):min(r0+h-1,H); dc = max(c0,1):min(c0+w-1,W);
    if isempty(dr)||isempty(dc), continue; end
    C(dr,dc) = max(C(dr,dc), As(dr-r0+1, dc-c0+1));
end
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
