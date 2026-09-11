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
mouseRoot  = 'C:\Users\Admin\Desktop\ChAT_MAP';
% REPOINTED. The old path (Desktop\260522_ChAT_soma_G8s) no longer exists -- that
% data moved to D:\260522_ChAT_soma_G8s, whose matlab_stitch holds only
% stitched_ch3_avgz.tif. ChAT_MAP\matlab_stitch is the only place with BOTH
% stitch_manual_coords.mat and stitched_ch3_avgproj.tif, which is what
% Ventral_surface_ROI_on_vessel_260528.m needs to regenerate the calibration.
%
% BLOCKER: ROI_on_vessel_data.mat does not exist anywhere on this machine, and it is
% what converts stage coordinates to cartoon-canvas coordinates. It is produced by
% Ventral_surface_ROI_on_vessel_260528.m, which is INTERACTIVE (you click the vessel
% branch to set the reference). So that has to be run once, by hand, before this
% script can draw anything -- including the Sert cells added below.
%
% Why the exclusion in coh_cfg_260727 matters here: stage_to_canvas is ONE affine
% applied to every animal, so the whole map assumes every session zeroed its stage
% on the same anatomical landmark. The first 260721_Sert FOV did not, which is
% exactly why it must stay excluded from anything spatial.
vesselMat  = fullfile(mouseRoot,'matlab_stitch','ROI_on_vessel_data.mat');

sumRoot    = 'D:\Ventral_surface_summary';
% coherence/significance source = SVD-motion pipeline output (was the older
% coherence_polar_260528). Same schema (PP.r, PP.th, labels, confC) and same
% label format 'group/date/recName/rid', so this is a drop-in source swap.
dataMat    = fullfile(sumRoot,'polar_coh_vs_rayleigh_260729','polar_coh_vs_rayleigh_data.mat');
% REPOINTED from the coherence run to the RAYLEIGH-on-linear-phase run.
%   significance : logZ >= logZ_crit, i.e. alpha = 0.05   (was coherence r >= confC
%                  at alpha = 0.001 -- only 14 of 368 ROIs cleared it)
%   phase        : the occupancy-weighted circular MEDIAN of the cell's event
%                  phases (RY.thMed), on the same convention the colorbar already
%                  uses -- 0 = insp onset, pi = insp peak.
% Same linear peak<->trough phase, same 100 ms lead correction, so the map and the
% polar plot now say the same thing about the same cells.
%
% NOTE the Rayleigh result is per CELL while the map draws per ROI observation, so
% obsOfCell is used to hand every observation of a cell that cell's phase and
% verdict. A cell imaged in three recordings therefore paints three FOVs the same
% colour, which is the point.

% ---- cell identity: report CELLS, not raw Cellpose masks --------------------
% A mask is not a cell. The same neuron re-imaged at a second Z or a second zoom
% produces one mask per recording, so summing masks over recordings overcounts.
% The cross-FOV matcher resolves that, and its cell_link.mat carries the mapping
%     (recording name, maskL label) -> cell_id
% which is exactly the portable key the curation is built on. Registered per
% session below.
%
% Two things the matcher output encodes that must be honoured:
%   status 'tossed'  (cell_id NaN)  -> rejected in curation, NOT a cell. Dropped.
%   status 'ungrouped'              -> a real cell, just seen in one recording.
%                                      Counts as one cell.
% Sessions with no matcher run fall back to one cell per ROI, which is what
% cell_link_260727 itself does without curation -- so the count degrades
% gracefully rather than silently reverting to mask counts.
% Read from the ARCHIVE's own copies. These pointed at D:, the sessions then
% moved to E:, and E: is removable -- so the lookup silently failed and every
% group fell back to "1 cell per mask", overcounting Sert 402 instead of 296.
% The failure is a warning, not an error, which is exactly why it has to be a
% path that cannot go missing. Copied in beside each session 2026-07-31.
cell_link_sources = { ...
    'Sert',   '0721', fullfile(sumRoot, 'Sert',   '0721', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0728', fullfile(sumRoot, 'Vglut2', '0728', 'cell_pooled', 'cell_link.mat')
    'Vgat',   '0730', fullfile(sumRoot, 'Vgat',   '0730', 'cell_pooled', 'cell_link.mat')
    };

scan_dirs   = {'ChAT','Vglut2','Vgat','Sst','Sert'};   % + Sert (260721 baseline)
groups      = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};

PixelSizeBase = 1.7778;

% ---- stage-axis -> anatomy map, PER SESSION --------------------------------
% The map's whole geometry assumes the animal was mounted the same way every
% session, so that stage x is the medio-lateral axis and stage y the
% rostro-caudal one. Two sessions were mounted ROTATED 90 deg, so their stage
% axes carry the opposite anatomical meaning:
%
%     -stage x = anterior (rostral)  ;  stage y = lateral (+y = left)
%
% Encoded as a 2x2 taking stage coordinates to anatomical ones:
%
%     [lateral ; rostral] = axmap * [stage x ; stage y]
%
% The lateral SIGN is not consequential -- everything is folded into the left
% hemisphere below with -abs() -- so the real content of the rotated row is
% "rostral comes from -stage x", and the FOV footprint transposes with it.
%
% Keyed on SESSION (scan dir + date), not genotype: how the animal was mounted
% is a property of that day's prep. Vglut2 0224 and 1124 are NOT rotated, only
% 0728 is -- so a genotype-level flag would have moved the older pFN data too.
%
% Independent corroboration for both, from the experimenter's own filenames:
% the rotated sessions write the coordinate pair y-first (`y1900x1320`), and
% their motorPosition satisfies motorX = -nameY, motorY = nameX.
AXMAP_DEFAULT = [ 1  0 ;  0  1 ];
AXMAP_ROT90   = [ 0  1 ; -1  0 ];
% Sessions acquired on the NEW mount / coordinate axis. This is a fact about the
% rig on the day, NOT something to infer from where the stage coordinates land --
% Vgat/0730's footprint happens to overlap the old-mount cloud, so guessing from
% coordinates got it wrong. Everything from 260721 onward is on the new axis.
rotated_sessions = { 'Sert','0721' ; 'Vglut2','0728' ; 'Vgat','0730' };

% ---- one patch per imaging SITE --------------------------------------------
% Several sites were re-imaged at a second Z plane, or (ro7/roi7) recorded twice
% under near-identical names. Each such recording used to stamp its own tinted
% patch at essentially the same place, which is what made a site look drawn
% "repetitively" -- 0.85-alpha patches piling up read as one dark smear.
% So the PATCH is drawn once per site, from the recording with the most ROIs.
% Every ROI of every recording is still kept in R: two Z planes are two
% different sets of cells and only the rendering was duplicated, not the data.
mergeSiteUm = 100;      % same site if FOV centres are closer than this (um)

% per-group tints (FOV patch colors)
patch_tints = struct( ...
    'IO',     [0.5    0.5    0.5   ], ...   % black
    'ChAT',   [0.85 0.10 0.10], ...   % red
    'Vglut2', [0.10 0.65 0.20], ...   % green
    'Vgat',   [0.10 0.30 0.85], ...   % blue
    'Sst',    [0.55 0.20 0.75], ...   % purple
    'Sert',   [0.90 0.45 0.10]);      % orange, matching cfg.genotype_color
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
sig_map = containers.Map('KeyType','char','ValueType','logical');
r_map   = containers.Map('KeyType','char','ValueType','double');
th_map  = containers.Map('KeyType','char','ValueType','double');
crit_val = 0;  critName = '';
% Significance level for THIS FIGURE, set here rather than inherited from the
% polar run. The polar draws two circles (alpha 0.05 inner, 0.001 outer) and
% saves the 0.05 one as logZ_crit; the map wants the stricter threshold, so it
% derives its own. Z = n*Rbar^2 with p ~ exp(-Z), hence logZ_crit = log(-log(alpha))
% -- n-independent, so one cut is valid for every cell regardless of event count.
sig_alpha = 0.001;      % 0.05 -> logZ 1.10 | 0.001 -> logZ 1.93
if isfile(dataMat)
    Sray = load(dataMat,'RY','labels','obsOfCell','logZ_crit','alpha_ray');
    critZ   = log(-log(sig_alpha));
    sigCell = Sray.RY.logZ >= critZ;                    % per CELL
    for i = 1:numel(Sray.labels)
        c = Sray.obsOfCell(i);                          % this observation's cell
        lab = Sray.labels{i};
        sig_map(lab) = sigCell(c);
        r_map(lab)   = Sray.RY.logZ(c);
        th_map(lab)  = mod(Sray.RY.thMed(c), 2*pi);     % circular median, 0..2pi
    end
    crit_val = critZ;
    critName = sprintf('logZ\\geq%.2f, \\alpha=%.3g', critZ, sig_alpha);
    fprintf('  phase/significance from %s\n', dataMat);
    fprintf('    %d observations -> %d cells, %d significant at alpha=%.3g (logZ>=%.2f)\n', ...
            numel(Sray.labels), numel(sigCell), nnz(sigCell), sig_alpha, critZ);
    fprintf('    (polar run itself used alpha=%.2g, logZ>=%.2f -> %d cells)\n', ...
            Sray.alpha_ray, Sray.logZ_crit, nnz(Sray.RY.logZ >= Sray.logZ_crit));
else
    warning('No polar_coh_vs_rayleigh_data.mat -- everyone non-sig.');
end

%% ---- cell-identity lookup ----
% keys are '<scanDir>/<date>/<recName>/<maskL label>'; values are a cell key that
% is namespaced by session, since cell ids restart at 1 in every matcher run.
cell_map   = containers.Map('KeyType','char','ValueType','char');
tossed_set = containers.Map('KeyType','char','ValueType','logical');
for s = 1:size(cell_link_sources,1)
    pre = sprintf('%s/%s', cell_link_sources{s,1}, cell_link_sources{s,2});
    f   = cell_link_sources{s,3};
    if ~isfile(f)
        warning('cell_link missing for %s (%s) -- that session will count masks as cells.', pre, f);
        continue;
    end
    Lk = load(f,'link');  T = Lk.link.obsT;
    for i = 1:height(T)
        kk = sprintf('%s/%s/%d', pre, T.rec_name(i), T.maskL_label(i));
        if isnan(T.cell_id(i)), tossed_set(kk) = true;
        else,                   cell_map(kk) = sprintf('%s#c%d', pre, T.cell_id(i));
        end
    end
    fprintf('  cell identity %-14s %4d masks -> %4d cells (%d tossed in curation)\n', ...
            pre, height(T), numel(unique(T.cell_id(~isnan(T.cell_id)))), nnz(isnan(T.cell_id)));
end

%% ---- walk every FOV ----
R = struct('group',{},'date',{},'fov',{},'rid',{},'is_sig',{}, ...
           'x_um',{},'y_um',{},'r_coh',{},'th_coh',{},'cell_key',{});
F = struct('fov',{},'gname',{},'avgPath',{},'maskPath',{}, ...
           'x_um',{},'y_um',{}, ...
           'H_fov',{},'W_fov',{},'px_um_fov',{},'sig_rids',{}, ...
           'sess',{},'site',{},'nroi',{},'rot90',{},'drawPatch',{});

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

        % which way was this animal mounted?
        isRot = any(strcmp(rotated_sessions(:,1), sname) & ...
                    strcmp(rotated_sessions(:,2), dateStr));
        if isRot, axmap = AXMAP_ROT90; else, axmap = AXMAP_DEFAULT; end

        % FOV centre: stage->anatomical axes->canvas->rotated->μm; then y-up
        [cx_F, cy_F] = stage_to_canvas_ax(stage_to_canvas, axmap, motor(1), motor(2));
        dxy = Rmat * [cx_F - bx; cy_F - by];
        x_um = dxy(1) * um_per_px_stitch;
        y_um = -dxy(2) * um_per_px_stitch;            % negate so +y = rostral

        ah = dir(fullfile(fp,'*_AVG_for_CP.tif'));
        if ~isempty(ah), avgPath = fullfile(ah(1).folder, ah(1).name); else, avgPath = ''; end

        lbls = setdiff(unique(maskL(:)), 0);
        props = regionprops(maskL, 'Centroid'); props = props(lbls);
        sig_rids_here = [];
        nTossHere = 0;
        for j = 1:numel(lbls)
            rid = lbls(j);
            lab = sprintf('%s/%s/%s/%d', gname, dateStr, fov, rid);
            % cell-identity key uses the SCAN DIR (sname), not gname: the matcher
            % ran per session, before the IO folder override was applied.
            ckin = sprintf('%s/%s/%s/%d', sname, dateStr, fov, rid);
            if isKey(tossed_set, ckin)
                nTossHere = nTossHere + 1;  continue;   % rejected in curation
            end
            if isKey(cell_map, ckin), cell_key = cell_map(ckin);
            else,                     cell_key = ['mask:' lab];   % no matcher run
            end
            is_sig = isKey(sig_map, lab) && sig_map(lab);
            if is_sig, sig_rids_here(end+1) = rid; end %#ok<AGROW>
            if isKey(r_map,  lab), r_coh  = r_map(lab);  else, r_coh  = NaN; end
            if isKey(th_map, lab), th_coh = th_map(lab); else, th_coh = NaN; end

            % per-ROI position: centroid -> stage um -> anatomical axes -> canvas
            % -> rotated -> y-up um.  The scan mirrors are fixed to the RIG, so a
            % FOV column is always stage x and a row always stage y no matter how
            % the animal was mounted; the mount only enters at the stage->anatomy
            % step, which is where axmap is applied.
            c = props(j).Centroid;                       % [column row]
            mx_r = motor(1) + (c(1) - Wf/2) * px_um_fov;
            my_r = motor(2) + (c(2) - Hf/2) * px_um_fov;
            [cxR, cyR] = stage_to_canvas_ax(stage_to_canvas, axmap, mx_r, my_r);
            dxyR = Rmat * [cxR - bx; cyR - by];
            x_umR =  dxyR(1) * um_per_px_stitch;
            y_umR = -dxyR(2) * um_per_px_stitch;

            R(end+1) = struct('group',gname,'date',dateStr,'fov',fov, ...
                              'rid',rid,'is_sig',is_sig, ...
                              'x_um',x_umR,'y_um',y_umR, ...
                              'r_coh',r_coh,'th_coh',th_coh, ...
                              'cell_key',cell_key); %#ok<SAGROW>
        end
        if nTossHere > 0
            fprintf('  %-42s dropped %d mask(s) tossed in curation\n', fov(1:min(42,end)), nTossHere);
        end

        F(end+1) = struct('fov',fov,'gname',gname,'avgPath',avgPath, ...
                          'maskPath',maskPath, ...
                          'x_um',x_um,'y_um',y_um, ...
                          'H_fov',Hf,'W_fov',Wf,'px_um_fov',px_um_fov, ...
                          'sig_rids',sig_rids_here, ...
                          'sess',[sname '/' dateStr],'site',subcell, ...
                          'nroi',numel(lbls),'rot90',isRot,'drawPatch',true); %#ok<SAGROW>
    end
end
assert(~isempty(F), 'No FOVs.');

%% ---- optional: tie-point registered coordinates -----------------------------
% Ventral_surface_tiepoint_register_260731 solves one translation per FOV from
% cells the matcher saw in more than one FOV. Applied HERE, before the site
% merge, deliberately: the merge groups FOVs by centre distance, so correcting
% centres afterwards regroups nothing and can pull a group apart (measured: one
% 8-member group went 80 -> 174 um when the correction was applied after).
% FOVs with no tie-point get zero shift and are left on stage coordinates, so a
% map made this way MIXES corrected and uncorrected placements -- F.reg_linked
% records which is which rather than hiding it.
USE_REGISTERED = false;      % <-- set true to draw from corrected coords
% Which offset file. Both carry the same fovTable schema:
%   ROI_on_cartoon_registered.mat        automatic, tie-points from matched cells
%   ROI_on_cartoon_manual_offsets.mat    hand-nudged in Ventral_surface_manual_align
REG_SOURCE = 'manual';       % 'auto' | 'manual'
if strcmpi(REG_SOURCE,'manual')
    regMat = fullfile(fileparts(cartoonMat), 'ROI_on_cartoon_manual_offsets.mat');
else
    regMat = fullfile(fileparts(cartoonMat), 'ROI_on_cartoon_registered.mat');
end
if USE_REGISTERED
    assert(isfile(regMat), 'USE_REGISTERED but %s is missing.', regMat);
    R2 = load(regMat, 'fovTable');  T2 = R2.fovTable;
    kF = string({F.sess}) + "/" + string({F.fov});
    nMoved = 0;
    for i = 1:numel(F)
        j = find(T2.fov == kF(i), 1);
        F(i).reg_linked = false;
        if isempty(j), continue; end
        F(i).x_um = F(i).x_um + T2.tx_um(j);
        F(i).y_um = F(i).y_um + T2.ty_um(j);
        F(i).reg_linked = T2.linked(j);
        if T2.linked(j), nMoved = nMoved + 1; end
    end
    kR = string({R.group}) + "/" + string({R.date}) + "/" + string({R.fov});
    for i = 1:numel(R)
        j = find(T2.fov == kR(i), 1);
        if isempty(j), continue; end
        R(i).x_um = R(i).x_um + T2.tx_um(j);
        R(i).y_um = R(i).y_um + T2.ty_um(j);
    end
    fprintf('  REGISTERED coordinates applied: %d of %d FOVs moved (rest have no tie-point)\n', ...
        nMoved, numel(F));
else
    for i = 1:numel(F), F(i).reg_linked = false; end
end

%% ---- collapse re-imaged sites to ONE patch each ----
% MERGE_SITES = false -> no collapsing at all: every FOV draws its own patch and
% overlapping fields simply overlay. That shows the true coverage instead of one
% representative field per site. The reason merging existed is that 0.85-alpha
% patches stack into a dark smear, so alpha is scaled down here to compensate --
% with 9 fields on one site the accumulated opacity is what needs controlling,
% not the number of fields.
MERGE_SITES = true;
if ~MERGE_SITES
    nOv = zeros(numel(F),1);
    for u = 1:numel(F)
        for v = 1:numel(F)
            if u==v || ~strcmp(F(u).sess,F(v).sess), continue; end
            if hypot(F(u).x_um-F(v).x_um, F(u).y_um-F(v).y_um) <= mergeSiteUm
                nOv(u) = nOv(u) + 1;
            end
        end
    end
    patch_alpha = patch_alpha / max(1, sqrt(median(nOv(nOv>0)+1)));
    fprintf('  MERGE_SITES=false: all %d patches drawn, alpha scaled to %.2f\n', ...
        numel(F), patch_alpha);
end
if MERGE_SITES
% Same session + same group + FOV centres within mergeSiteUm  ->  one site. The
% survivor is the recording with the most ROIs, so the patch shown is the one
% carrying the most cells. R is untouched: this suppresses duplicate DRAWING
% only, it does not drop any cell.
key = strcat({F.sess}, '|', {F.gname});
for k = unique(key)
    idx = find(strcmp(key, k{1}));
    for a = 1:numel(idx)
        ia = idx(a);
        if ~F(ia).drawPatch, continue; end
        for b = a+1:numel(idx)
            ib = idx(b);
            if ~F(ib).drawPatch, continue; end
            if hypot(F(ia).x_um-F(ib).x_um, F(ia).y_um-F(ib).y_um) > mergeSiteUm, continue; end
            if F(ib).nroi > F(ia).nroi, drop = ia; keepi = ib; else, drop = ib; keepi = ia; end
            F(drop).drawPatch = false;
            fprintf('  merged patch: %-42s -> %-42s (%s, %.0f um apart)\n', ...
                    F(drop).fov(1:min(42,end)), F(keepi).fov(1:min(42,end)), ...
                    F(drop).sess, hypot(F(ia).x_um-F(ib).x_um, F(ia).y_um-F(ib).y_um));
            if drop == ia, break; end
        end
    end
end
end   % if MERGE_SITES
nRot = nnz([F.rot90]);
fprintf('FOVs: %d (patches drawn: %d, %d on rotated-mount sessions)   ROIs: %d\n', ...
        numel(F), nnz([F.drawPatch]), nRot, numel(R));

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

    % Rotated mount -> the FOV's own axes are 90 deg off the anatomy, so the
    % footprint transposes: a 320x512 px FOV that spanned 569 um of stage x now
    % spans 569 um ROSTRALLY and only 356 um laterally. rot90 carries both the
    % texture and the extent across. (The residual handedness of the texture is
    % not resolved -- these patches are decorative soft blobs, and the existing
    % non-rotated ones already draw row 1 at high y. Position and footprint are
    % what carry information here, and both are correct.)
    if F(u).rot90, aN = rot90(aN, -1); end

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

    % Only the TINTED BLOB is suppressed at a re-imaged site -- that is what piled
    % up and smeared. The sig-ROI outlines below still draw for every recording,
    % each in its OWN extent (they carry which cells were significant, and a thin
    % white outline does not smear). Masks from different recordings cannot be
    % overlaid on one patch anyway: different zoom, pixel size and Z.
    if F(u).drawPatch
        image(ax_p, xext, yext, rgbPatch, 'AlphaData', alphaPatch);
    end

    % --- white sig-ROI outlines, same extent ---
    if ~isempty(F(u).sig_rids) && isfile(F(u).maskPath)
        Sm = load(F(u).maskPath,'maskL');
        sig_mask = ismember(Sm.maskL, F(u).sig_rids);
        if F(u).rot90, sig_mask = rot90(sig_mask, -1); end   % same frame as aN
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
cb.TickLabels = {'0 (onset)','\pi/2','\pi (peak)','3\pi/2','2\pi'};
cb.Label.String = 'breath phase of Ca events (circular median)';

% origin + title
plot(ax_p, 0, 0, 'r+', 'MarkerSize',16, 'LineWidth',2);
nCellAll  = numel(unique({R.cell_key}));
nSigCell  = numel(unique({R(logical([R.is_sig])).cell_key}));
title(ax_p, sprintf(['%d cells (%d masks) / %d FOVs, %d sites   ' ...
      '(left: anatomy;  right: Rayleigh phase, %d cells significant at %s)'], ...
      nCellAll, numel(R), numel(F), nnz([F.drawPatch]), nSigCell, critName));

% legend
hL = []; lL = {};
gn = fieldnames(patch_tints);
fprintf('\n--- masks vs CELLS per group ---\n');
fprintf('%-8s %8s %8s  %s\n','group','masks','cells','note');
for gi = 1:numel(gn)
    m    = strcmp({R.group}, gn{gi});
    nAll = sum(m);
    if nAll == 0, continue; end
    nCell = numel(unique({R(m).cell_key}));
    matched = nnz(~startsWith({R(m).cell_key}, 'mask:'));
    fprintf('%-8s %8d %8d  %s\n', gn{gi}, nAll, nCell, ...
        ternary(matched==0, 'no matcher run -- 1 cell per mask', ...
                sprintf('%d masks cross-FOV matched', matched)));
    col = patch_tints.(gn{gi});
    hL(end+1) = patch(NaN,NaN, col, 'EdgeColor','k'); %#ok<SAGROW>
    lL{end+1} = sprintf('%s (%d cells)', gn{gi}, nCell); %#ok<SAGROW>
end
% inside the axes, not 'eastoutside': outside it lands on top of the colorbar and
% hides its 'preferred phase (rad)' label. The caudal-right quadrant is empty.
if ~isempty(hL)
    legend(ax_p, hL, lL, 'Location','southeast','FontSize',8, 'Color','w');
end

%% ---- save ----
outDir = fileparts(C.imgPath);
% Registered runs write to their own names so the stage-coordinate map is never
% overwritten -- the two have to remain comparable.
if USE_REGISTERED, sfx = '_reg'; else, sfx = ''; end
exportgraphics(figOut, fullfile(outDir,['ROI_on_cartoon' sfx '.png']), 'Resolution',200,'BackgroundColor','white');
exportgraphics(figOut, fullfile(outDir,['ROI_on_cartoon' sfx '.pdf']), 'ContentType','vector','BackgroundColor','white');
save(fullfile(outDir,['ROI_on_cartoon_data' sfx '.mat']), 'R','F','xd','yd', ...
     'bx','by','theta_rad','um_per_px_stitch','ax','ay', ...
     'rotated_sessions','AXMAP_DEFAULT','AXMAP_ROT90','mergeSiteUm');
fprintf('Saved overlay PNG/PDF + .mat to %s\n', outDir);

%% ---- where every rotated-mount FOV landed (so the fix is checkable) ----
ir = find([F.rot90]);
if ~isempty(ir)
    fprintf('\n--- rotated-mount sessions: placement ---\n');
    fprintf('%-8s %-7s %-40s %9s %9s %6s %s\n', ...
            'sess','site','recording','lat(um)','ros(um)','patch','footprint');
    for i = ir
        fprintf('%-8s %-7s %-40s %9.0f %9.0f %6s %4.0f x %4.0f um\n', ...
            F(i).sess, F(i).site, F(i).fov(1:min(40,end)), F(i).x_um, F(i).y_um, ...
            string(F(i).drawPatch), F(i).H_fov*F(i).px_um_fov, F(i).W_fov*F(i).px_um_fov);
    end
end

%% ========================= local helpers =========================
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function [cx, cy] = stage_to_canvas_ax(s2c, axmap, mx, my)
% Stage coordinates -> stitched-canvas pixels, via the session's stage->anatomy
% axis map. axmap sends [stage x; stage y] to [lateral; rostral]; the lateral
% component is then folded into the LEFT hemisphere, exactly as before, so a
% non-rotated session reduces to the original stage_to_canvas(-abs(mx), my).
v = axmap * [mx; my];
[cx, cy] = s2c(-abs(v(1)), v(2));
end
