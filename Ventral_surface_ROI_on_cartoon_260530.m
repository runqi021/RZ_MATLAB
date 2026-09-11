% Ventral_surface_ROI_on_cartoon_260530.m
% -----------------------------------------------------------------------
%  FOV patches (per-FOV AVG image, circular soft-mask, group tint, white
%  outline on sig ROIs) overlaid on the brainstem cartoon background.
%  Coordinate convention = same as the working dot version:
%      data x = +lateral  (μm)
%      data y = +rostral  (μm)         y-up convention (YDir='normal')
%  ROIs from every ca_spike_data.mat under <sumRoot>\{groups} are computed
%  in μm relative to the branch via the vessel-derived calibration.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =====================
cartoonMat = 'C:\Users\Admin\Desktop\brainstem_map_cartoon_coordframe.mat';
mouseRoot  = 'C:\Users\Admin\Desktop\260522_ChAT_soma_G8s';
vesselMat  = fullfile(mouseRoot,'matlab_stitch','ROI_on_vessel_data.mat');

sumRoot    = 'D:\Ventral_surface_summary';
% coherence/significance source = SVD-motion pipeline output (was the older
% coherence_polar_260528). Same schema (PP.r, PP.th, labels, confC) and same
% label format 'group/date/recName/rid', so this is a drop-in source swap.
dataMat    = fullfile(sumRoot,'coherence_polar_svd_260606','coherence_polar_data.mat');

scan_dirs   = {'ChAT','Vglut2','Vgat','Sst'};
groups      = {'IO','ChAT','Vglut2','Vgat','Sst'};

PixelSizeBase = 1.7778;

% per-group tints (FOV patch colors)
patch_tints = struct( ...
    'IO',     [0.5    0.5    0.5   ], ...   % black
    'ChAT',   [0.85 0.10 0.10], ...   % red
    'Vglut2', [0.10 0.65 0.20], ...   % green
    'Vgat',   [0.10 0.30 0.85], ...   % blue
    'Sst',    [0.55 0.20 0.75]);      % purple
patch_alpha    = 0.85;
patch_clip_pct = [1 99.5];
patch_gamma    = 0.6;
% =========================================================

%% ---- load cartoon coordframe ----
C = load(cartoonMat);
assert(isfile(C.imgPath), 'Cartoon image missing: %s', C.imgPath);
Icart = im2double(imread(C.imgPath));
if size(Icart,3) == 1, Icart = repmat(Icart,1,1,3); end
% Light gaussian blur to soften pixelation (the cartoon was a raster)
bg_blur_sigma = 1.2;            % px ; bump to 2-3 if still pixelated
Icart = imgaussfilt(Icart, bg_blur_sigma);
[Hc, Wc, ~] = size(Icart);
% cartoon data extent (μm). pix_per_y carries sign of user's click.
xd = ([1 Wc] - C.origin_pix(1)) / C.pix_per_x;
yd = ([1 Hc] - C.origin_pix(2)) / C.pix_per_y;

%% ---- load vessel-derived calibration ----
V = load(vesselMat, 'bx','by','theta_rad','ax','ay','um_per_px_stitch');
bx = V.bx; by = V.by; theta_rad = V.theta_rad;
um_per_px_stitch = V.um_per_px_stitch;
ax = V.ax;  ay = V.ay;
stage_to_canvas = @(mx, my) deal(ax(1)*mx + ax(2), ay(1)*my + ay(2));
Rmat = [cos(theta_rad)  sin(theta_rad);
       -sin(theta_rad)  cos(theta_rad)];

%% ---- significance + coherence lookup ----
if isfile(dataMat)
    Scoh = load(dataMat,'PP','labels','confC');
    sig_map = containers.Map(Scoh.labels, num2cell(Scoh.PP.r >= Scoh.confC));
    r_map   = containers.Map(Scoh.labels, num2cell(Scoh.PP.r));
    th_map  = containers.Map(Scoh.labels, num2cell(Scoh.PP.th));
    confC_val = Scoh.confC;
else
    warning('No coherence_polar_data.mat -- everyone non-sig.');
    sig_map = containers.Map('KeyType','char','ValueType','logical');
    r_map   = containers.Map('KeyType','char','ValueType','double');
    th_map  = containers.Map('KeyType','char','ValueType','double');
    confC_val = 0;
end

%% ---- walk every FOV ----
R = struct('group',{},'date',{},'fov',{},'rid',{},'is_sig',{}, ...
           'x_um',{},'y_um',{},'r_coh',{},'th_coh',{});
F = struct('fov',{},'gname',{},'avgPath',{},'maskPath',{}, ...
           'x_um',{},'y_um',{}, ...
           'H_fov',{},'W_fov',{},'px_um_fov',{},'sig_rids',{});

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
        if ~isempty(regexpi(subcell,'IO','once')), gname = 'IO'; end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), continue; end

        mh = dir(fullfile(fp,'*_ch1_meta.mat'));
        if isempty(mh), mh = dir(fullfile(fp,'*_meta.mat')); end
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

        maskPath = fullfile(sam_hits(hh).folder, sam_hits(hh).name);
        S2 = load(maskPath, 'maskL');
        if ~isfield(S2,'maskL') || isempty(S2.maskL), continue; end
        maskL = S2.maskL; [Hf, Wf] = size(maskL);

        % FOV centre: stage->canvas->rotated->μm; then data convention y-up
        [cx_F, cy_F] = stage_to_canvas(-abs(motor(1)), motor(2));   % fold to LEFT
        dxy = Rmat * [cx_F - bx; cy_F - by];
        x_um = dxy(1) * um_per_px_stitch;
        y_um = -dxy(2) * um_per_px_stitch;            % negate so +y = rostral

        ah = dir(fullfile(fp,'*_AVG_for_CP.tif'));
        if ~isempty(ah), avgPath = fullfile(ah(1).folder, ah(1).name); else, avgPath = ''; end

        lbls = setdiff(unique(maskL(:)), 0);
        props = regionprops(maskL, 'Centroid'); props = props(lbls);
        sig_rids_here = [];
        for j = 1:numel(lbls)
            rid = lbls(j);
            lab = sprintf('%s/%s/%s/%d', gname, dateStr, fov, rid);
            is_sig = isKey(sig_map, lab) && sig_map(lab);
            if is_sig, sig_rids_here(end+1) = rid; end %#ok<AGROW>
            if isKey(r_map,  lab), r_coh  = r_map(lab);  else, r_coh  = NaN; end
            if isKey(th_map, lab), th_coh = th_map(lab); else, th_coh = NaN; end

            % per-ROI position: centroid -> stage um -> canvas -> rotated -> y-up um
            c = props(j).Centroid;
            mx_r = motor(1) + (c(1) - Wf/2) * px_um_fov;
            my_r = motor(2) + (c(2) - Hf/2) * px_um_fov;
            [cxR, cyR] = stage_to_canvas(-abs(mx_r), my_r);   % fold LEFT
            dxyR = Rmat * [cxR - bx; cyR - by];
            x_umR =  dxyR(1) * um_per_px_stitch;
            y_umR = -dxyR(2) * um_per_px_stitch;

            R(end+1) = struct('group',gname,'date',dateStr,'fov',fov, ...
                              'rid',rid,'is_sig',is_sig, ...
                              'x_um',x_umR,'y_um',y_umR, ...
                              'r_coh',r_coh,'th_coh',th_coh); %#ok<SAGROW>
        end

        F(end+1) = struct('fov',fov,'gname',gname,'avgPath',avgPath, ...
                          'maskPath',maskPath, ...
                          'x_um',x_um,'y_um',y_um, ...
                          'H_fov',Hf,'W_fov',Wf,'px_um_fov',px_um_fov, ...
                          'sig_rids',sig_rids_here); %#ok<SAGROW>
    end
end
assert(~isempty(F), 'No FOVs.');
fprintf('FOVs: %d   ROIs: %d\n', numel(F), numel(R));

%% ---- figure: cartoon background, y-up ----
figOut = figure('Color','w','Name','ROIs on cartoon', ...
                'Units','normalized','Position',[0.04 0.04 0.9 0.9]);
ax_p = axes(figOut);
image(ax_p, xd, yd, Icart); hold(ax_p,'on');
set(ax_p,'YDir','normal');
axis(ax_p,'image'); box(ax_p,'on'); grid(ax_p,'on');
xlabel(ax_p,'X (\mum)  (lateral)');
ylabel(ax_p,'Y (\mum)  (rostral +)');

% --- overlay each FOV patch as its own image() in data coords ---
for u = 1:numel(F)
    if isempty(F(u).avgPath) || ~isfile(F(u).avgPath), continue; end
    gname_u = F(u).gname;
    if isfield(patch_tints, gname_u), tintU = patch_tints.(gname_u);
    else, tintU = [0.10 0.80 0.30]; end

    A = double(imread(F(u).avgPath));
    a_lo = prctile(A(:), patch_clip_pct(1));
    a_hi = prctile(A(:), patch_clip_pct(2));
    aN   = max(0, min(1, (A - a_lo)/max(a_hi-a_lo,eps))).^patch_gamma;

    % patch tinted RGB in source pixel size
    [hp, wp] = size(aN);
    [Yp, Xp]  = ndgrid(linspace(-1,1,hp), linspace(-1,1,wp));
    rho       = sqrt(Xp.^2 + Yp.^2);
    % Strong low-frequency perturbation -> clearly non-circular blob
    rho_noise = imgaussfilt(randn(hp, wp), max(hp,wp)/8) * 0.7;
    circ_mask = max(0, min(1, (1 - (rho + rho_noise))*3));
    alphaPatch = patch_alpha .* circ_mask .* aN;

    rgbPatch = zeros(hp, wp, 3);
    for ch = 1:3, rgbPatch(:,:,ch) = tintU(ch); end

    % world extent in μm  (centre at FOV's x_um, y_um; y-up so flip vertically)
    halfW_um = (wp/2) * F(u).px_um_fov;
    halfH_um = (hp/2) * F(u).px_um_fov;
    xext = [F(u).x_um - halfW_um, F(u).x_um + halfW_um];
    yext = [F(u).y_um + halfH_um, F(u).y_um - halfH_um];  % image row 1 at top (high y)

    image(ax_p, xext, yext, rgbPatch, 'AlphaData', alphaPatch);

    % --- white sig-ROI outlines, same extent ---
    if ~isempty(F(u).sig_rids) && isfile(F(u).maskPath)
        Sm = load(F(u).maskPath,'maskL');
        sig_mask = ismember(Sm.maskL, F(u).sig_rids);
        outline = sig_mask & ~imerode(sig_mask, strel('disk',1));
        outline = imdilate(outline, strel('disk',1));
        outline = outline & (circ_mask > 0.5);
        if any(outline(:))
            rgbOut = ones(hp, wp, 3);          % white
            image(ax_p, xext, yext, rgbOut, 'AlphaData', double(outline));
        end
    end
end

% --- right hemisphere: per-ROI dots colored by preferred phase (HSV) ---
% Phase convention from coherency: foot=0, peak=pi.  HSV wraps cleanly.
phase_rot = 0;        % rotate HSV so phase 0 lands at desired hue (radians)
rx = -[R.x_um];                     % mirror LEFT->RIGHT
ry =  [R.y_um];
th = mod([R.th_coh] - phase_rot, 2*pi);
sigv = [R.is_sig];

cmap = hsv(256);
ci = round(th / (2*pi) * 255) + 1;
ci(isnan(ci)) = 1;
cols = cmap(min(max(ci,1),256), :);
cols(isnan(th),:) = repmat([0.6 0.6 0.6], sum(isnan(th)),1);

scatter(ax_p, rx( sigv), ry( sigv), 36, cols( sigv,:), 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.3);
scatter(ax_p, rx(~sigv), ry(~sigv), 10, cols(~sigv,:), ...
        'MarkerEdgeColor','none');

% Linear phase colorbar with breath-landmark ticks
colormap(ax_p, hsv(256));
caxis(ax_p, [0, 2*pi]);
cb = colorbar(ax_p,'Location','eastoutside');
cb.Ticks = [0 pi/2 pi 3*pi/2 2*pi];
cb.TickLabels = {'0 (foot)','\pi/2','\pi (peak)','3\pi/2','2\pi'};
cb.Label.String = 'preferred phase (rad)';

% origin + title
plot(ax_p, 0, 0, 'r+', 'MarkerSize',16, 'LineWidth',2);
title(ax_p, sprintf('%d ROIs / %d FOVs   (left: anatomy, right: coherence; confC=%.2f)', ...
      numel(R), numel(F), confC_val));

% legend
hL = []; lL = {};
gn = fieldnames(patch_tints);
for gi = 1:numel(gn)
    nAll = sum(strcmp({R.group}, gn{gi}));
    if nAll == 0, continue; end
    col = patch_tints.(gn{gi});
    hL(end+1) = patch(NaN,NaN, col, 'EdgeColor','k'); %#ok<SAGROW>
    lL{end+1} = sprintf('%s (%d ROIs)', gn{gi}, nAll); %#ok<SAGROW>
end
if ~isempty(hL), legend(ax_p, hL, lL, 'Location','eastoutside','FontSize',8); end

%% ---- save ----
outDir = fileparts(C.imgPath);
exportgraphics(figOut, fullfile(outDir,'ROI_on_cartoon.png'), 'Resolution',200,'BackgroundColor','white');
exportgraphics(figOut, fullfile(outDir,'ROI_on_cartoon.pdf'), 'ContentType','vector','BackgroundColor','white');
save(fullfile(outDir,'ROI_on_cartoon_data.mat'), 'R','F','xd','yd', ...
     'bx','by','theta_rad','um_per_px_stitch','ax','ay');
fprintf('Saved overlay PNG/PDF + .mat to %s\n', outDir);
