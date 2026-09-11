% Ventral_surface_ROI_phasehist_on_cartoon_260607.m
% -----------------------------------------------------------------------
%  Like Ventral_surface_ROI_on_cartoon, BUT each ROI is drawn as a little
%  ROSE glyph of its spike-phase HISTOGRAM (not a single mean-phase dot):
%      - 24 phase wedges, angular position = breath phase (onset=0 right,
%        peak=pi), radius proportional to spikes/cycle % in that bin,
%        hue = phase (cyclic HSV). So a sharply tuned cell shows a single
%        bright spike of color; a broad cell shows a full multi-color rose.
%  Left hemisphere = FOV anatomy patches; right hemisphere = phase-hist roses.
%
%  Histogram is computed with the EXACT svd-pipeline convention (occupancy-
%  compensated spk/cycle %, nDrop=30, ca_lag=0.1 s, Vglut2/1124 +1 shift,
%  piecewise onset=0/peak=pi) copied from
%  Ventral_surface_phase_polar_allROI_svd_260606, so glyphs match the polar
%  analyses exactly.
%
%  Source: coherence_polar_svd_260606/coherence_polar_data.mat
%  Dependencies: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); addpath(scriptDir);
addpath(fullfile(scriptDir,'2p_breathing_coherence'));

%% ===================== USER-EDITABLE =====================
cartoonMat = 'C:\Users\Admin\Desktop\brainstem_map_cartoon_coordframe.mat';
mouseRoot  = 'C:\Users\Admin\Desktop\260522_ChAT_soma_G8s';
vesselMat  = fullfile(mouseRoot,'matlab_stitch','ROI_on_vessel_data.mat');

sumRoot    = 'D:\Ventral_surface_summary';
dataMat    = fullfile(sumRoot,'coherence_polar_svd_260606','coherence_polar_data.mat');

scan_dirs   = {'ChAT','Vglut2','Vgat','Sst'};
groups      = {'IO','ChAT','Vglut2','Vgat','Sst'};

PixelSizeBase = 1.7778;

% --- phase histogram (svd-pipeline convention) ---
nPhaseBins  = 24;
ca_lag_sec  = 0.1;        % spikes lead-shifted 0.1 s (3 fr @30 Hz)
nDrop       = 30;
fallback_fps= 30;
minSpikes   = 1;          % ROIs with fewer spikes -> faint dot, no rose

% --- glyph appearance ---
glyphR_um   = 50;         % rose max radius (um, data coords)
normMode    = 'global';      % 'roi' = each rose scaled to its own max (shape);
                          % 'global' = common scale (compare magnitudes)
sigAlpha    = 0.92;       % wedge opacity for significant ROIs
nsigAlpha   = 0.45;       % wedge opacity for non-significant ROIs
showFOVpatch= true;       % draw FOV anatomy patches behind glyphs

% per-group tints (FOV patch colors)
patch_tints = struct( ...
    'IO',     [0.5    0.5    0.5   ], ...
    'ChAT',   [0.85 0.10 0.10], ...
    'Vglut2', [0.10 0.65 0.20], ...
    'Vgat',   [0.10 0.30 0.85], ...
    'Sst',    [0.55 0.20 0.75]);
patch_alpha    = 0.85;
patch_clip_pct = [1 99.5];
patch_gamma    = 0.6;
% =========================================================

%% ---- load cartoon coordframe ----
C = load(cartoonMat);
assert(isfile(C.imgPath), 'Cartoon image missing: %s', C.imgPath);
Icart = im2double(imread(C.imgPath));
if size(Icart,3) == 1, Icart = repmat(Icart,1,1,3); end
Icart = imgaussfilt(Icart, 1.2);
[Hc, Wc, ~] = size(Icart);
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

%% ---- significance lookup ----
assert(isfile(dataMat), 'No coherence data: %s', dataMat);
Scoh = load(dataMat,'PP','labels','confC');
sig_map = containers.Map(Scoh.labels, num2cell(Scoh.PP.r >= Scoh.confC));
confC_val = Scoh.confC;

edgesP = linspace(0, 2*pi, nPhaseBins+1);
ctrsP  = (edgesP(1:end-1)+edgesP(2:end))/2;

%% ---- walk every FOV: position + per-ROI phase histogram ----
R = struct('group',{},'date',{},'fov',{},'rid',{},'is_sig',{}, ...
           'x_um',{},'y_um',{},'phist',{},'nspk',{});
F = struct('fov',{},'gname',{},'avgPath',{},'maskPath',{}, ...
           'x_um',{},'y_um',{},'H_fov',{},'W_fov',{},'px_um_fov',{},'sig_rids',{});
Fcache = containers.Map();

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
        if isempty(find(strcmp(groups, gname), 1)), continue; end

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

        [cx_F, cy_F] = stage_to_canvas(-abs(motor(1)), motor(2));
        dxy = Rmat * [cx_F - bx; cy_F - by];
        x_um = dxy(1) * um_per_px_stitch;
        y_um = -dxy(2) * um_per_px_stitch;

        ah = dir(fullfile(fp,'*_AVG_for_CP.tif'));
        if ~isempty(ah), avgPath = fullfile(ah(1).folder, ah(1).name); else, avgPath = ''; end

        % per-FOV phase data (breath phase + spikes), cached. May be absent.
        Fdat = [];
        try
            if ~isKey(Fcache, fp), Fcache(fp) = load_folder(fp, nDrop, fallback_fps); end
            Fdat = Fcache(fp);
        catch
            Fdat = [];
        end

        lbls = setdiff(unique(maskL(:)), 0);
        props = regionprops(maskL, 'Centroid'); props = props(lbls);
        sig_rids_here = [];
        for j = 1:numel(lbls)
            rid = lbls(j);
            lab = sprintf('%s/%s/%s/%d', gname, dateStr, fov, rid);
            is_sig = isKey(sig_map, lab) && sig_map(lab);
            if is_sig, sig_rids_here(end+1) = rid; end %#ok<AGROW>

            phist = nan(1, nPhaseBins); nspk = 0;
            if ~isempty(Fdat)
                try
                    [phist, nspk] = roi_phase_hist(Fdat, rid, ca_lag_sec, edgesP);
                catch
                    phist = nan(1, nPhaseBins); nspk = 0;
                end
            end

            c = props(j).Centroid;
            mx_r = motor(1) + (c(1) - Wf/2) * px_um_fov;
            my_r = motor(2) + (c(2) - Hf/2) * px_um_fov;
            [cxR, cyR] = stage_to_canvas(-abs(mx_r), my_r);
            dxyR = Rmat * [cxR - bx; cyR - by];
            x_umR =  dxyR(1) * um_per_px_stitch;
            y_umR = -dxyR(2) * um_per_px_stitch;

            R(end+1) = struct('group',gname,'date',dateStr,'fov',fov, ...
                              'rid',rid,'is_sig',is_sig, ...
                              'x_um',x_umR,'y_um',y_umR, ...
                              'phist',phist,'nspk',nspk); %#ok<SAGROW>
        end

        F(end+1) = struct('fov',fov,'gname',gname,'avgPath',avgPath, ...
                          'maskPath',maskPath,'x_um',x_um,'y_um',y_um, ...
                          'H_fov',Hf,'W_fov',Wf,'px_um_fov',px_um_fov, ...
                          'sig_rids',sig_rids_here); %#ok<SAGROW>
    end
end
assert(~isempty(F), 'No FOVs.');
fprintf('FOVs: %d   ROIs: %d   (with phase hist: %d)\n', ...
        numel(F), numel(R), nnz(arrayfun(@(s) s.nspk>=minSpikes, R)));

%% ---- figure: cartoon background, y-up ----
figOut = figure('Color','w','Name','ROI phase-hist roses on cartoon', ...
                'Units','normalized','Position',[0.04 0.04 0.9 0.9], 'Renderer','opengl');
ax_p = axes(figOut);
image(ax_p, xd, yd, Icart); hold(ax_p,'on');
set(ax_p,'YDir','normal'); axis(ax_p,'image'); box(ax_p,'on'); grid(ax_p,'on');
xlabel(ax_p,'X (\mum)  (lateral)'); ylabel(ax_p,'Y (\mum)  (rostral +)');

% --- FOV anatomy patches (left hemisphere context) ---
if showFOVpatch
    for u = 1:numel(F)
        if isempty(F(u).avgPath) || ~isfile(F(u).avgPath), continue; end
        if isfield(patch_tints, F(u).gname), tintU = patch_tints.(F(u).gname);
        else, tintU = [0.10 0.80 0.30]; end
        A = double(imread(F(u).avgPath));
        a_lo = prctile(A(:), patch_clip_pct(1)); a_hi = prctile(A(:), patch_clip_pct(2));
        aN   = max(0, min(1, (A - a_lo)/max(a_hi-a_lo,eps))).^patch_gamma;
        [hp, wp] = size(aN);
        [Yp, Xp]  = ndgrid(linspace(-1,1,hp), linspace(-1,1,wp));
        rho       = sqrt(Xp.^2 + Yp.^2);
        rho_noise = imgaussfilt(randn(hp, wp), max(hp,wp)/8) * 0.7;
        circ_mask = max(0, min(1, (1 - (rho + rho_noise))*3));
        alphaPatch = patch_alpha .* circ_mask .* aN;
        rgbPatch = zeros(hp, wp, 3);
        for ch = 1:3, rgbPatch(:,:,ch) = tintU(ch); end
        halfW_um = (wp/2) * F(u).px_um_fov; halfH_um = (hp/2) * F(u).px_um_fov;
        xext = [F(u).x_um - halfW_um, F(u).x_um + halfW_um];
        yext = [F(u).y_um + halfH_um, F(u).y_um - halfH_um];
        image(ax_p, xext, yext, rgbPatch, 'AlphaData', alphaPatch);
    end
end

% --- right hemisphere: per-ROI phase-histogram ROSES ---
cmapHSV = hsv(256);
wedgeHalf = pi / nPhaseBins;          % half angular width of a wedge
% color per phase-bin center
binRGB = cmapHSV(min(max(round(ctrsP/(2*pi)*255)+1,1),256), :);
% global normalization value if requested
allP = cat(1, R(arrayfun(@(s) any(isfinite(s.phist)) && s.nspk>=minSpikes, R)).phist);
gmax = max(allP(:));  if isempty(gmax) || ~(gmax>0), gmax = 1; end

nRose = 0;
for k = 1:numel(R)
    p = R(k).phist;
    if R(k).nspk < minSpikes || ~any(isfinite(p)) || ~(max(p)>0)
        continue;        % weak/no-phase ROIs handled below as faint dots
    end
    rx = -R(k).x_um;  ry = R(k).y_um;        % mirror LEFT->RIGHT hemisphere
    if strcmp(normMode,'roi'), nv = max(p); else, nv = gmax; end
    aA = R(k).is_sig*sigAlpha + ~R(k).is_sig*nsigAlpha;
    for b = 1:nPhaseBins
        rr = glyphR_um * p(b) / nv;
        if ~(rr > 0), continue; end
        aa = linspace(ctrsP(b)-wedgeHalf, ctrsP(b)+wedgeHalf, 4);
        xs = [rx, rx + rr*cos(aa)];
        ys = [ry, ry + rr*sin(aa)];
        patch(ax_p, xs, ys, binRGB(b,:), 'EdgeColor','none', 'FaceAlpha', aA);
    end
    if R(k).is_sig                       % thin ring marks significant ROIs
        tt = linspace(0,2*pi,30);
        plot(ax_p, rx + glyphR_um*cos(tt), ry + glyphR_um*sin(tt), '-', ...
            'Color',[0 0 0 0.5], 'LineWidth',0.4);
    end
    nRose = nRose + 1;
end

% --- faint gray dots for ROIs without a usable phase histogram ---
noP = arrayfun(@(s) s.nspk < minSpikes || ~any(isfinite(s.phist)) || ~(max(s.phist)>0), R);
scatter(ax_p, -[R(noP).x_um], [R(noP).y_um], 6, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha',0.4);

% phase colorbar (cyclic)
colormap(ax_p, hsv(256)); caxis(ax_p, [0 2*pi]);
cb = colorbar(ax_p,'Location','eastoutside');
cb.Ticks = [0 pi/2 pi 3*pi/2 2*pi];
cb.TickLabels = {'0 (foot)','\pi/2','\pi (peak)','3\pi/2','2\pi'};
cb.Label.String = 'breath phase (rad)';

plot(ax_p, 0, 0, 'r+', 'MarkerSize',16, 'LineWidth',2);
title(ax_p, sprintf(['%d ROIs / %d FOVs   |   right: per-ROI spike-phase rose ' ...
    '(radius = spk/cyc %%, hue = phase; norm=%s; ring = sig, confC=%.2f)'], ...
    numel(R), numel(F), normMode, confC_val), 'Interpreter','tex');

%% ---- save ----
outDir = fileparts(C.imgPath);
exportgraphics(figOut, fullfile(outDir,'ROI_phasehist_on_cartoon.png'), ...
    'Resolution',220,'BackgroundColor','white');
exportgraphics(figOut, fullfile(outDir,'ROI_phasehist_on_cartoon.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
save(fullfile(outDir,'ROI_phasehist_on_cartoon_data.mat'), 'R','F','ctrsP','confC_val');
fprintf('Saved %d roses -> ROI_phasehist_on_cartoon.png/.pdf + .mat to %s\n', nRose, outDir);

%% ===================== LOCAL FUNCTIONS (copied from phase_polar_allROI_svd) =====
function F = load_folder(folderPath, nDrop, fallback_fps)
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
assert(~isempty(bp),'no breath_peak_pc1'); assert(~isempty(ip),'no breath_insp_start_pc1');
assert(isfile(fullfile(folderPath,'ca_spike_data.mat')),'no ca_spike_data');
F.fps = detect_session_fps(folderPath, fallback_fps);
BP = load(fullfile(bp(1).folder, bp(1).name));
IP = load(fullfile(ip(1).folder, ip(1).name));
CA = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes');
nB = numel(BP.breath);
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
ev(1:min(nDrop,numel(ev))) = [];
ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;
ef(1:min(nDrop,numel(ef))) = [];
if contains(folderPath, fullfile('Vglut2','1124'))
    ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
end
F.ev = ev; F.ef = ef; F.roi_spikes = CA.roi_spikes;
F.nDff = Inf;
if ~isempty(df), D = load(fullfile(df(1).folder, df(1).name),'dFF'); F.nDff = size(D.dFF,1); end
end

function [pctP, nspk] = roi_phase_hist(F, roi, ca_lag_sec, edgesP)
spk = double(F.roi_spikes(roi).spike_train(:));
T = min([F.nDff, numel(F.ev), numel(F.ef), numel(spk)]);
ev = F.ev(1:T); ef = F.ef(1:T);
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);
phi = piecewise_phase_local(find(ev>0), find(ef>0), T);
lag = round(ca_lag_sec*F.fps);
sp = find(spk>0) - lag; sp = sp(sp>=1 & sp<=T);
spk_phi = mod(phi(sp), 2*pi); spk_phi = spk_phi(~isnan(spk_phi));
phv = phi(~isnan(phi)); nCyc = (max(phv)-min(phv))/(2*pi);
occ = histcounts(mod(phv,2*pi), edgesP); Ob = mean(occ(occ>0));
wOcc = ones(1,numel(occ)); wOcc(occ>0) = Ob ./ occ(occ>0);
pctP = (100 * histcounts(spk_phi, edgesP) / max(nCyc,eps)) .* wOcc;
nspk = numel(spk_phi);
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];
[events, ord] = sort(events); types = types(ord);
keep = true(size(events));
for i = 2:numel(events), if types(i) == types(i-1), keep(i) = false; end, end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases = nan(size(events)); phi_cur = types(1) * pi;
for i = 1:numel(events), phases(i) = phi_cur; phi_cur = phi_cur + pi; end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
