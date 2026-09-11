% Ventral_surface_ROI_on_vessel_260528.m
% -----------------------------------------------------------------------
%  Overlay ALL ROIs from Ventral_surface_summary onto the stitched 260522
%  vessel mosaic, with a user-clicked branching point as the (0,0) origin.
%  2D only.
%
%  Workflow:
%   1. Show stitched ch3 avg-projection.
%   2. ginput(1)  -> click branching point.
%   3. Walk every ca_spike_data.mat folder under <sumRoot>\{ChAT,Vglut2,Vgat,Sst}\
%      compute every ROI centroid in stage um (fold by abs(x) so all
%      hemispheres collapse onto the side shown by the vessel image).
%   4. Map stage um -> canvas px using a fit derived from the 260522
%      stitch's per-tile motor positions (parsed from filenames) vs
%      canvas placements (x, y) in stitch_manual_coords.mat.
%   5. Subtract the click point, convert to um, plot.
%
%  Caveat: only 0522 FOVs are GEOMETRICALLY accurate on this vessel image.
%  Other mice are placed assuming their stage origin matches 0522's --
%  this is approximate cross-animal registration, useful for a population
%  overview, not for cell-level co-registration.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =================================
mouseRoot = 'C:\Users\Admin\Desktop\260522_ChAT_soma_G8s';
stitchTif = fullfile(mouseRoot,'matlab_stitch','stitched_ch3_avgproj.tif');
stitchMat = fullfile(mouseRoot,'matlab_stitch','stitch_manual_coords.mat');

sumRoot   = 'D:\Ventral_surface_summary';
dataMat   = fullfile(sumRoot,'coherence_polar_260528','coherence_polar_data.mat');

scan_dirs   = {'ChAT','Vglut2','Vgat','Sst'};
groups      = {'IO','ChAT','Vglut2','Vgat','Sst'};
group_colors = [0    0    0;
                0.85 0.10 0.10;
                0.10 0.65 0.20;
                0.10 0.30 0.85;
                0.55 0.20 0.75];

PixelSizeBase    = 1.7778;
um_per_px_stitch = 0.8889;

% --- vessel image preprocessing (matches earlier segmentation pipeline) ---
clip_pct   = [0.5 99.5];
gamma_val  = 0.6;
clahe_tiles = [32 32];
clahe_clip  = 0.01;
add_vessel_tint = false;   % true -> red-tint the vesselness map; false -> plain enhanced grayscale
invert_vessel  = true;     % true -> dark vessels on light background (anatomy-figure style)
add_fov_patches = true;    % overlay each FOV's AVG_for_CP image as a soft-circle patch
% per-group patch colors (RGB).  IO uses ChAT's red since it's a ChAT subset.
patch_tints = struct( ...
    'IO',     [0.85 0.10 0.10], ...   % red
    'ChAT',   [0.85 0.10 0.10], ...   % red
    'Vglut2', [1.00 0.85 0.10], ...   % yellow
    'Vgat',   [0.10 0.80 0.30], ...   % green
    'Sst',    [1.00 0.55 0.10]);      % orange
patch_alpha    = 0.85;     % patch opacity (0 = invisible, 1 = fully overrides bg)
patch_clip_pct = [1 99.5]; % per-patch percentile clamp
patch_gamma    = 0.6;

ms_sig    = 36;
ms_nonsig = 12;
% =====================================================================

%% ---- show stitched + click branching point ----
imS = double(imread(stitchTif));
lo = prctile(imS(:), clip_pct(1)); hi = prctile(imS(:), clip_pct(2));
imN_raw = max(0, min(1, (imS - lo)/max(hi-lo,eps))).^gamma_val;       % clip + gamma
imN = adapthisteq(imN_raw, 'NumTiles', clahe_tiles, 'ClipLimit', clahe_clip);  % local contrast
if add_vessel_tint
    V = fibermetric(imN, 2:2:12, 'StructureSensitivity', 0.05, 'ObjectPolarity','bright');
    V = V / max(V(:) + eps);
    % blend a red tint of the vesselness onto the grayscale background
    rgb = repmat(imN, [1 1 3]);
    rgb(:,:,1) = min(1, rgb(:,:,1) + 0.45*V);
    rgb(:,:,2) = rgb(:,:,2) .* (1 - 0.25*V);
    rgb(:,:,3) = rgb(:,:,3) .* (1 - 0.25*V);
    imN = rgb;       % imN is now HxWx3; downstream code handles either
end
if invert_vessel
    imN = 1 - imN;
end

figClick = figure('Color','w','Name','Click branch then midline-rostral point', ...
                  'Units','normalized','Position',[0.05 0.05 0.9 0.88]);
if size(imN,3) == 3
    imshow(imN); hold on;
else
    imshow(imN, 'InitialMagnification','fit'); hold on;
end
title({'Click TWO points:', ...
       '  1) branch point (bottom) = (0,0)', ...
       '  2) somewhere along midline toward the rostral side (defines +Y)'}, 'FontSize',11);
[cx_clk, cy_clk] = ginput(2);
bx = cx_clk(1); by = cy_clk(1);
mx_clk = cx_clk(2); my_clk = cy_clk(2);
plot(bx, by, 'r+', 'MarkerSize', 18, 'LineWidth', 2);
plot(mx_clk, my_clk, 'yo', 'MarkerSize', 12, 'LineWidth', 2);
plot([bx mx_clk], [by my_clk], 'y-', 'LineWidth', 1.5);

% rotation angle that brings the (bx,by)->(mx_clk,my_clk) vector to point UP
% (image y-down: up = -y). CCW-in-image rotation by theta uses
% R = [cos sin; -sin cos], and theta = atan2(vx, -vy) does the job.
vx = mx_clk - bx;  vy = my_clk - by;
theta_rad = atan2(vx, -vy);
theta_deg = rad2deg(theta_rad);
fprintf('Branch: canvas px = (%.1f, %.1f)\n', bx, by);
fprintf('Midline rotation: %.2f deg CCW (in image-y-down convention)\n', theta_deg);

%% ---- stage_um -> canvas_px (from 0522 stitch) ----
p = load(stitchMat);
motorX = nan(numel(p.files),1); motorY = motorX;
for t = 1:numel(p.files)
    nm = char(p.files(t));
    tok = regexp(nm, '_x(-?\d+)_y(-?\d+)_', 'tokens', 'once');
    if ~isempty(tok), motorX(t) = str2double(tok{1}); motorY(t) = str2double(tok{2}); end
end
canvasCenterX = p.x(:) + p.W0/2;
canvasCenterY = p.y(:) + p.H0/2;
ax = polyfit(motorX, canvasCenterX, 1);
ay = polyfit(motorY, canvasCenterY, 1);
fprintf('Stage-um -> canvas-px:  px_x = %.4f * mx + %.2f   px_y = %.4f * my + %.2f\n', ...
        ax(1), ax(2), ay(1), ay(2));

stage_to_canvas = @(mx, my) deal(ax(1)*mx + ax(2), ay(1)*my + ay(2));

% -x = x folding: cross-mouse plotting on one hemisphere.
% All FOVs are simply mapped to abs(motor_x) so left/right hemispheres
% collapse onto the same side, regardless of each mouse's stage zero.

%% ---- sig lookup ----
if isfile(dataMat)
    Scoh = load(dataMat,'PP','labels','confC');
    sig_map = containers.Map(Scoh.labels, num2cell(Scoh.PP.r >= Scoh.confC));
else
    warning('No coherence_polar_data.mat -- everyone marked non-sig.');
    sig_map = containers.Map('KeyType','char','ValueType','logical');
end

%% ---- walk every FOV ----
R = struct('group',{},'date',{},'fov',{},'rid',{},'mx',{},'my',{},'mz',{}, ...
           'cx',{},'cy',{},'gi',{},'is_sig',{});
F = struct('fov',{},'gname',{},'avgPath',{},'maskPath',{},'cx',{},'cy',{}, ...
           'H_fov',{},'W_fov',{},'px_um_fov',{});

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(sumRoot, sname);
    if ~isfolder(gdir), continue; end
    sam_hits = dir(fullfile(gdir,'**','*cpSAM_output.mat'));
    for hh = 1:numel(sam_hits)
        fp  = sam_hits(hh).folder;
        rel = strrep(fp, gdir, ''); rel = regexprep(rel,'^[\\/]+','');
        parts = regexp(rel,'[\\/]','split');
        if numel(parts) < 3, continue; end
        dateStr = parts{1}; subcell = parts{2};
        [~,n,e] = fileparts(fp); fov = [n e];

        gname = sname;
        if strcmpi(sname,'ChAT') && ~isempty(regexpi(subcell,'IO','once'))
            gname = 'IO';
        end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), continue; end

        mh = dir(fullfile(fp,'*_ch1_meta.mat')); if isempty(mh), mh = dir(fullfile(fp,'*_meta.mat')); end
        if isempty(mh), continue; end
        M = load(fullfile(mh(1).folder, mh(1).name));
        if ~isfield(M,'motorPosition') || numel(M.motorPosition) < 3, continue; end
        motor = M.motorPosition(:).';
        if isfield(M,'pixelSize_um') && M.pixelSize_um>0
            px_um_fov = M.pixelSize_um;
        elseif isfield(M,'zoomFactor') && M.zoomFactor>0
            px_um_fov = PixelSizeBase / M.zoomFactor;
        else, continue;
        end

        S2 = load(fullfile(sam_hits(hh).folder, sam_hits(hh).name), 'maskL');
        if ~isfield(S2,'maskL') || isempty(S2.maskL), continue; end
        maskL = S2.maskL; [H, W] = size(maskL);
        lbls  = setdiff(unique(maskL(:)), 0);
        props = regionprops(maskL, 'Centroid'); props = props(lbls);

        % per-FOV record: center on canvas, -x = x folded
        [cx_F, cy_F]  = stage_to_canvas(abs(motor(1)), motor(2));
        ah = dir(fullfile(fp,'*_AVG_for_CP.tif'));
        if ~isempty(ah), avgPath = fullfile(ah(1).folder, ah(1).name); else, avgPath = ''; end
        maskPath = fullfile(sam_hits(hh).folder, sam_hits(hh).name);
        F(end+1) = struct('fov',fov,'gname',gname,'avgPath',avgPath,'maskPath',maskPath, ...
                          'cx',cx_F,'cy',cy_F, ...
                          'H_fov',H,'W_fov',W,'px_um_fov',px_um_fov); %#ok<SAGROW>

        for j = 1:numel(props)
            c   = props(j).Centroid;
            rid = lbls(j);
            mx_raw = motor(1) + (c(1) - W/2) * px_um_fov;
            my     = motor(2) + (c(2) - H/2) * px_um_fov;
            % -x = x: fold around stage x=0 so L/R hemispheres collapse
            mx = abs(mx_raw);
            [cx, cy] = stage_to_canvas(mx, my);
            lab = sprintf('%s/%s/%s/%d', gname, dateStr, fov, rid);
            is_sig = isKey(sig_map, lab) && sig_map(lab);
            R(end+1) = struct('group',gname,'date',dateStr,'fov',fov,'rid',rid, ...
                              'mx',mx_raw,'my',my,'mz',motor(3), ...
                              'cx',cx,'cy',cy,'gi',gi,'is_sig',is_sig); %#ok<SAGROW>
        end
    end
end

assert(~isempty(R), 'No ROIs.');
fprintf('Total ROIs collected: %d\n', numel(R));

%% ---- composite FOV patches onto the (inverted) vessel canvas ----
if add_fov_patches && ~isempty(F)
    if ~ismatrix(imN), canvas = imN;             % already RGB
    else,              canvas = repmat(imN, [1 1 3]); end   % gray -> RGB
    [Hc, Wc, ~] = size(canvas);

    for u = 1:numel(F)
        if isempty(F(u).avgPath) || ~isfile(F(u).avgPath), continue; end

        % --- per-group tint ---
        gname_u = F(u).gname;
        if isfield(patch_tints, gname_u), tintU = patch_tints.(gname_u);
        else, tintU = [0.10 0.80 0.30]; end

        % --- patch image ---
        A = double(imread(F(u).avgPath));
        a_lo = prctile(A(:), patch_clip_pct(1));
        a_hi = prctile(A(:), patch_clip_pct(2));
        aN   = max(0, min(1, (A - a_lo)/max(a_hi-a_lo,eps))).^patch_gamma;

        hh_canvas = max(8, round(F(u).H_fov * F(u).px_um_fov / um_per_px_stitch));
        ww_canvas = max(8, round(F(u).W_fov * F(u).px_um_fov / um_per_px_stitch));
        aR  = imresize(aN, [hh_canvas ww_canvas]);
        [hp, wp] = size(aR);

        [Yp, Xp]  = ndgrid(linspace(-1, 1, hp), linspace(-1, 1, wp));
        rho       = sqrt(Xp.^2 + Yp.^2);
        circ_mask = max(0, min(1, (1 - rho)*4));

        r0 = round(F(u).cy - hp/2);  c0 = round(F(u).cx - wp/2);
        rr = max(1, r0):min(Hc, r0 + hp - 1);
        cc = max(1, c0):min(Wc, c0 + wp - 1);
        if isempty(rr) || isempty(cc), continue; end
        pr = (rr - r0 + 1);  pc = (cc - c0 + 1);

        alpha = patch_alpha * circ_mask(pr, pc) .* aR(pr, pc);

        for ch = 1:3
            canvas(rr, cc, ch) = canvas(rr, cc, ch) .* (1 - alpha) + ...
                                 tintU(ch) * alpha;
        end

        % --- white outline around SIG ROIs in this FOV (from cpSAM mask) ---
        sig_rids = [R(strcmp({R.fov}, F(u).fov) & [R.is_sig]).rid];
        if ~isempty(sig_rids) && ~isempty(F(u).maskPath) && isfile(F(u).maskPath)
            Sm = load(F(u).maskPath, 'maskL');
            maskL = Sm.maskL;
            sig_mask_native = ismember(maskL, sig_rids);
            % outline = mask minus its erosion
            outline_native = sig_mask_native & ~imerode(sig_mask_native, strel('disk',1));
            outline_resized = imresize(outline_native, [hp wp], 'nearest');
            outline_resized = imdilate(outline_resized, strel('disk',1));   % thicken so it survives rotation/AA
            % only show outline inside the circle window
            outline_resized = outline_resized & (circ_mask > 0.5);
            ol_sub = outline_resized(pr, pc);
            if any(ol_sub(:))
                for ch = 1:3
                    tmp = canvas(rr, cc, ch);
                    tmp(ol_sub) = 1;            % white
                    canvas(rr, cc, ch) = tmp;
                end
            end
        end
    end
    imN = canvas;
    fprintf('Composited %d FOV patches onto vessel canvas.\n', numel(F));
end

%% ---- apply midline-calibration rotation, then convert to um ----
% rotation matrix (image-y-down, CCW-in-image)
Rmat = [cos(theta_rad)  sin(theta_rad);
       -sin(theta_rad)  cos(theta_rad)];

dxy_px   = [[R.cx] - bx; [R.cy] - by];        % 2 x N raw offsets in canvas px
dxy_rot  = Rmat * dxy_px;                      % rotated offsets in canvas px
rel_um_x = dxy_rot(1,:) * um_per_px_stitch;
rel_um_y = dxy_rot(2,:) * um_per_px_stitch;

%% ---- rotate vessel image around the branch point ----
% Pad so branch becomes the geometric center, then imrotate around center,
% then convert displayed coordinates to um relative to branch.
H_im = size(imN,1); W_im = size(imN,2);
padL = max(0, round(W_im - 2*bx));   padR = max(0, round(2*bx - W_im));
padT = max(0, round(H_im - 2*by));   padB = max(0, round(2*by - H_im));
imP = padarray(imN, [padT padL 0], 0, 'pre');
imP = padarray(imP, [padB padR 0], 0, 'post');
imR = imrotate(imP, theta_deg, 'bilinear', 'loose');
% imrotate(angle>0) rotates the image VISUALLY counterclockwise by `angle`
% degrees. Our theta_rad = atan2(vx, -vy) is exactly the visual-CCW angle
% needed to bring the click vector (vx, vy) to point straight up in the
% display, so we pass it directly (no sign flip).
Hr = size(imR,1); Wr = size(imR,2);
x_um_axis = ((1:Wr) - Wr/2) * um_per_px_stitch;
y_um_axis = ((1:Hr) - Hr/2) * um_per_px_stitch;

%% ---- final figure: rotated vessel + rotated ROIs in um ----
figOut = figure('Color','w','Name','ROIs on vessel (rotated, origin = branch)', ...
                'Units','normalized','Position',[0.04 0.04 0.9 0.9]);
if size(imR,3) == 3
    image(x_um_axis, y_um_axis, imR); hold on; axis image;
else
    imagesc(x_um_axis, y_um_axis, imR); hold on; axis image;
    colormap(gray(256));
end
set(gca,'YDir','reverse');
plot(0, 0, 'r+', 'MarkerSize',18, 'LineWidth',2);
text(30, 0, '(0,0) branch', 'Color','r','FontSize',10,'FontWeight','bold');
% draw the midline reference (vertical after rotation)
ylimcur = ylim;
plot([0 0], ylimcur, 'y--', 'LineWidth', 0.7);

% Sig neurons get a WHITE soma outline burned into the FOV patch above;
% no extra marker is plotted here.  Build a legend purely for the group
% palette so readers know which color = which cell type.
hL = []; lL = {};
gn = fieldnames(patch_tints);
for gi = 1:numel(gn)
    nAll = sum(strcmp({R.group}, gn{gi}));
    if nAll == 0, continue; end
    col = patch_tints.(gn{gi});
    hL(end+1) = patch(NaN,NaN, col, 'EdgeColor','k'); %#ok<SAGROW>
    lL{end+1} = sprintf('%s (%d ROIs)', gn{gi}, nAll); %#ok<SAGROW>
end
if ~isempty(hL), legend(hL, lL, 'Location','eastoutside','FontSize',8,'TextColor','k'); end

xlabel('X (\mum) relative to branch'); ylabel('Y (\mum) relative to branch');
title(sprintf('All %d ROIs on 260522 vessel mosaic  (cross-mouse, folded to one hemisphere)', numel(R)));
grid on; box on;

%% ---- save ----
outDir = fullfile(mouseRoot,'matlab_stitch');
exportgraphics(figOut, fullfile(outDir,'ROI_on_vessel_origin_at_branch.png'), 'Resolution',200,'BackgroundColor','white');
exportgraphics(figOut, fullfile(outDir,'ROI_on_vessel_origin_at_branch.pdf'), 'ContentType','vector','BackgroundColor','white');
save(fullfile(outDir,'ROI_on_vessel_data.mat'), 'R','bx','by','mx_clk','my_clk', ...
     'theta_rad','theta_deg','um_per_px_stitch','ax','ay');
fprintf('Saved overlay PNG/PDF + .mat to %s\n', outDir);
